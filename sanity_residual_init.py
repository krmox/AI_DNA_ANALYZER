"""GATE: an untrained residual model must reproduce the binomial baseline.

Run before any training. Loads a freshly constructed
:class:`ResidualBinomialMamba`, runs it over the real 15x test loci with no
gradient updates, and compares its decision score against the binomial v1 LLR
computed by ``binomial_baseline``.

The claim being tested is not "similar" but "identical": the residual head is
zeroed, so the only permitted difference is float32 storage of the prior
(~1e-4 absolute on LLRs of magnitude up to 170). Correlations must be 1.0 and
every call-level count must match exactly.

Exit status is 0 only if every check passes. A non-zero exit means the
architecture does not reproduce the prior and training must not start.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import torch

from binomial_baseline import BinomialVariantCaller
from config import LABEL_SNP
from llr_features import base_features_from_counts, raw_llr_from_counts
from model_residual_binomial import ResidualBinomialMamba
from residual_metrics import prf, ranking_metrics
from residual_prior import NUM_CLASSES, prior_logits_from_llr, snp_decision_score

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

SEQ_LEN = 64

#: Tolerance on the logit difference. float32 holds ~7 significant digits, so
#: an LLR of 170 round-trips to within ~1e-5; 1e-3 is loose enough never to
#: fire on storage precision and tight enough that any real architectural
#: leak (a non-zero bias, a stray activation) would blow straight past it.
LOGIT_TOLERANCE = 1e-3


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    """Pearson correlation, guarding the degenerate zero-variance case."""
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman rank correlation via Pearson on ranks."""
    rank_a = np.empty(a.size); rank_a[np.argsort(a, kind="stable")] = np.arange(a.size)
    rank_b = np.empty(b.size); rank_b[np.argsort(b, kind="stable")] = np.arange(b.size)
    return pearson(rank_a, rank_b)


@torch.no_grad()
def score_untrained_model(counts: np.ndarray, llr: np.ndarray, device: torch.device,
                          seed: int) -> np.ndarray:
    """Decision score of a freshly initialized residual model over all loci.

    Args:
        counts: ``[N, 10]`` cached count matrix.
        llr: ``[N]`` binomial LLRs supplying the prior.
        device: Inference device.
        seed: Seed for the backbone's random initialization; varying it is
            part of the test, since a correct zero-init makes the result
            independent of the backbone's weights.

    Returns:
        ``[N]`` SNP-versus-Normal logit gap.
    """
    torch.manual_seed(seed)
    model = ResidualBinomialMamba().to(device).eval()

    features = torch.tensor(base_features_from_counts(counts), dtype=torch.float)
    features = features.reshape(-1, SEQ_LEN, features.shape[-1])
    prior = torch.tensor(prior_logits_from_llr(llr), dtype=torch.float)
    prior = prior.reshape(-1, SEQ_LEN, NUM_CLASSES)

    out = []
    for begin in range(0, features.shape[0], 256):
        logits = model(pileup_features=features[begin:begin + 256].to(device),
                       prior_logits=prior[begin:begin + 256].to(device))
        out.append(snp_decision_score(logits).reshape(-1).cpu())
    return torch.cat(out).numpy().astype(np.float64)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-counts", required=True)
    parser.add_argument("--binom-v1-threshold", type=float, required=True,
                        help="Frozen validation-derived binomial LLR threshold.")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 20260811])
    parser.add_argument("--out", default="results/residual_init_sanity.json")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    blob = np.load(args.test_counts, allow_pickle=True)
    counts, labels = blob["counts"], blob["labels"]
    snp = labels == LABEL_SNP
    primary = (labels == 0) | snp

    llr = raw_llr_from_counts(counts, BinomialVariantCaller())
    binomial_calls = llr >= args.binom_v1_threshold

    report: dict = {
        "test_counts": args.test_counts,
        "binom_v1_threshold": args.binom_v1_threshold,
        "logit_tolerance": LOGIT_TOLERANCE,
        "loci": int(labels.size), "snp": int(snp.sum()),
        "binomial_reference": {
            **prf(binomial_calls[primary], snp[primary]),
            **ranking_metrics(llr[primary], snp[primary]),
        },
        "seeds": {},
    }
    logger.info("Binomial v1 reference @thr=%.2f: %s", args.binom_v1_threshold,
                json.dumps(report["binomial_reference"]))

    passed = True
    for seed in args.seeds:
        score = score_untrained_model(counts, llr, device, seed)
        calls = score >= args.binom_v1_threshold
        max_abs = float(np.abs(score - llr).max())
        stats = {**prf(calls[primary], snp[primary]),
                 **ranking_metrics(score[primary], snp[primary])}

        entry = {
            "pearson_logit_vs_llr": pearson(score[primary], llr[primary]),
            "spearman_logit_vs_llr": spearman(score[primary], llr[primary]),
            "max_abs_logit_difference": max_abs,
            "identical_calls": bool(np.array_equal(calls, binomial_calls)),
            "n_call_disagreements": int(np.sum(calls != binomial_calls)),
            "at_frozen_binomial_threshold": stats,
            "delta_vs_binomial": {
                key: stats[key] - report["binomial_reference"][key]
                for key in ("precision", "recall", "f1", "roc_auc", "pr_auc")
            },
        }
        checks = {
            "max_abs_difference_within_tolerance": max_abs <= LOGIT_TOLERANCE,
            "pearson_is_one": abs(entry["pearson_logit_vs_llr"] - 1.0) < 1e-9,
            "spearman_is_one": abs(entry["spearman_logit_vs_llr"] - 1.0) < 1e-9,
            "calls_identical": entry["identical_calls"],
            "tp_fp_fn_identical": all(
                stats[key] == report["binomial_reference"][key] for key in ("tp", "fp", "fn")),
            "auc_identical": abs(entry["delta_vs_binomial"]["roc_auc"]) < 1e-9,
            "ap_identical": abs(entry["delta_vs_binomial"]["pr_auc"]) < 1e-9,
        }
        entry["checks"] = checks
        entry["passed"] = all(checks.values())
        passed &= entry["passed"]
        report["seeds"][str(seed)] = entry

        logger.info(
            "seed %d | pearson=%.12f spearman=%.12f max|dlogit|=%.2e | "
            "AUC=%.6f AP=%.6f F1=%.6f TP=%d FP=%d FN=%d | %s",
            seed, entry["pearson_logit_vs_llr"], entry["spearman_logit_vs_llr"], max_abs,
            stats["roc_auc"], stats["pr_auc"], stats["f1"], stats["tp"], stats["fp"],
            stats["fn"], "PASS" if entry["passed"] else "FAIL")
        for name, ok in checks.items():
            if not ok:
                logger.error("  failed check: %s", name)

    report["passed"] = bool(passed)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    logger.info("Wrote %s", args.out)

    if not passed:
        logger.error("INITIALIZATION GATE FAILED - do not train. Fix the architecture.")
        return 1
    logger.info("INITIALIZATION GATE PASSED - untrained residual model == binomial v1.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
