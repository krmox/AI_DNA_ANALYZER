"""The complete initialization gate for the 30x PB-residual experiment.

The residual model's entire claim to being a *residual* is that before any
gradient step it reproduces the Poisson-binomial caller exactly. If that is
only approximately true, then every later difference is contaminated by an
initialization offset and the experiment measures nothing. This gate is
therefore a hard precondition: it exits non-zero and training does not run.

``sanity_pb_residual.py`` already checks the maximum absolute deviation on a
sample. This script adds the four properties the 30x protocol additionally
requires, over the **whole validation split** rather than a sample:

* Pearson correlation of the two score vectors == 1
* Spearman correlation == 1 (the ranking is what threshold-free metrics read)
* zero call disagreements at the frozen operating point
* identical TP / FP / FN

plus the maximum absolute deviation, which is expected to be non-zero at the
scale of float32 rounding and is checked against that scale rather than
against zero. The model runs in float32; the PB LLR is computed in float64 and
cast, so agreement to ~1e-6 on scores of order 10 is the correct expectation
and exact bitwise equality is not.

Labels are used here only to compute TP/FP/FN on the validation split. No test
data is opened by this script.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import torch

from config import LABEL_SNP
from evaluate_binomial_baseline import pick_threshold, prf
from evaluate_30x_pb_residual import validation_slice
from evaluate_quality_error import LLR_GRID
from locus_evidence import FEATURE_DIM
from model_pb_residual import PoissonBinomialResidualMamba
from residual_prior import snp_decision_score
from train_pb_residual import build_tensors

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

#: Deviation consistent with float32 representation of scores of order 10.
FLOAT32_TOLERANCE = 1e-4

#: Spearman tolerance. NOTE: a fixed absolute tolerance on this statistic is
#: not well posed, because its deviation grows with the number of loci -- more
#: values means more pairs closer together than one float32 ULP, hence more
#: rounding-induced rank swaps, none of which is a property of the model. On
#: 190k loci the deviation is 3.4e-08; on 4k loci it is 2.0e-10. The binding
#: criterion is therefore ``rank_swaps_below_float32_resolution``, which is
#: exact and sample-size independent; this bound is kept only so a gross
#: ranking change still trips a second wire.
#:
#: Exact 1.0 is not a float32-achievable criterion: the
#: model computes in float32, so LLR values closer together than one float32
#: ULP are rounded to the same or to inverted representations and their ranks
#: may swap. That perturbs the rank correlation at the 1e-10 level while
#: changing no decision. The *meaningful* form of the requirement is
#: :func:`max_rank_swap_gap` below -- no swap may involve values that float32
#: could actually have distinguished -- which is checked exactly. This
#: tolerance only keeps the reported correlation honest alongside it.
SPEARMAN_TOLERANCE = 1e-6


def max_rank_swap_gap(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Largest true gap between two values whose ranks the model swapped.

    If this is below one float32 ULP at the working scale, then every
    disagreement in ordering is between values float32 cannot tell apart, and
    the ranking is preserved exactly as far as the representation allows.

    Args:
        reference: Exact (float64) scores.
        candidate: Model scores.

    Returns:
        The largest reference-value gap involved in a rank swap, or 0.0 if the
        two orderings agree everywhere.
    """
    reference_order = np.argsort(reference, kind="stable")
    candidate_order = np.argsort(candidate, kind="stable")
    differing = np.nonzero(reference_order != candidate_order)[0]
    if differing.size == 0:
        return 0.0
    return float(np.max(np.abs(reference[reference_order[differing]]
                               - reference[candidate_order[differing]])))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman correlation via Pearson on ranks, ties averaged.

    Args:
        x: First variable.
        y: Second variable.

    Returns:
        The rank correlation.
    """
    def rank(values: np.ndarray) -> np.ndarray:
        order = np.argsort(values, kind="stable")
        ranks = np.empty(values.size, dtype=np.float64)
        ranks[order] = np.arange(values.size, dtype=np.float64)
        sorted_values = values[order]
        start = 0
        for stop in range(1, values.size + 1):
            if stop == values.size or sorted_values[stop] != sorted_values[start]:
                if stop - start > 1:
                    ranks[order[start:stop]] = ranks[order[start:stop]].mean()
                start = stop
        return ranks

    return float(np.corrcoef(rank(np.asarray(x, dtype=np.float64)),
                             rank(np.asarray(y, dtype=np.float64)))[0, 1])


@torch.no_grad()
def untrained_scores(features: np.ndarray, llr: np.ndarray, seed: int,
                     device: torch.device, use_gate: bool = False,
                     batch_size: int = 256) -> np.ndarray:
    """Decision scores from a freshly constructed, untrained model.

    Args:
        features: ``[N, FEATURE_DIM]`` locus evidence, N a multiple of SEQ_LEN.
        llr: ``[N]`` Poisson-binomial LLR.
        seed: Seed used to build the model, so the gate covers the same random
            backbone the training run will start from.
        device: Inference device.
        use_gate: Whether to build the gated parameterization.
        batch_size: Windows per forward pass.

    Returns:
        ``[N]`` SNP-minus-Normal logit gap.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = PoissonBinomialResidualMamba(feature_dim=FEATURE_DIM,
                                         use_gate=use_gate).to(device).eval()

    # Windowing and prior construction come from the training script itself, so
    # the gate cannot pass on a tensor layout the training run does not use.
    feature_tensor, prior_tensor = build_tensors(features, llr)

    out = []
    for begin in range(0, feature_tensor.shape[0], batch_size):
        logits = model(pileup_features=feature_tensor[begin:begin + batch_size].to(device),
                       prior_logits=prior_tensor[begin:begin + batch_size].to(device))
        out.append(snp_decision_score(logits).reshape(-1).cpu())
    return torch.cat(out).numpy().astype(np.float64)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-evidence", default="cache/train_30x_evidence.npz")
    parser.add_argument("--seed", type=int, action="append", default=[],
                        help="Seeds to gate; repeatable. Defaults to both experiment seeds.")
    parser.add_argument("--use-gate", action="store_true")
    parser.add_argument("--whole-file", action="store_true",
                        help="Gate every locus in the file instead of its last 10%%. "
                             "Used when the file is already a validation split.")
    parser.add_argument("--out", default="results/30x_pb_residual/init_gate.json")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    seeds = args.seed or [20260811, 424242]
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    blob = np.load(args.train_evidence, allow_pickle=True)
    val = (slice(None) if args.whole_file else validation_slice(blob["labels"].size))
    features = blob["features"][val]
    llr = blob["pb_llr"][val].astype(np.float64)
    labels = blob["labels"][val]
    snp = labels == LABEL_SNP
    frame = (labels == 0) | snp

    threshold, pb_stats = pick_threshold(llr[frame], snp[frame], LLR_GRID)
    logger.info("PB on validation: thr=%.1f F1=%.4f TP=%d FP=%d FN=%d",
                threshold, pb_stats["f1"], pb_stats["tp"], pb_stats["fp"], pb_stats["fn"])

    report = {"validation_loci": int(labels.size), "pb_threshold": threshold,
              "pb_validation": pb_stats, "seeds": {}, "tolerance": FLOAT32_TOLERANCE}
    passed = True

    for seed in seeds:
        score = untrained_scores(features, llr, seed, device, use_gate=args.use_gate)
        deviation = float(np.max(np.abs(score - llr)))
        pearson = float(np.corrcoef(score, llr)[0, 1])
        rank_correlation = spearman(score, llr)
        model_stats = prf(score[frame] >= threshold, snp[frame])
        disagreements = int(np.sum((score[frame] >= threshold) != (llr[frame] >= threshold)))
        swap_gap = max_rank_swap_gap(llr, score)
        float32_ulp = float(np.spacing(np.float32(np.abs(llr).max())))

        checks = {
            "max_abs_deviation_within_float32": deviation <= FLOAT32_TOLERANCE,
            "pearson_is_one": abs(pearson - 1.0) <= 1e-12,
            "spearman_is_one": abs(rank_correlation - 1.0) <= SPEARMAN_TOLERANCE,
            # The exact form of the ranking requirement: no ordering
            # disagreement may involve values float32 could have distinguished.
            "rank_swaps_below_float32_resolution": swap_gap <= float32_ulp,
            "zero_call_disagreements": disagreements == 0,
            "identical_confusion": (model_stats["tp"] == pb_stats["tp"]
                                    and model_stats["fp"] == pb_stats["fp"]
                                    and model_stats["fn"] == pb_stats["fn"]),
        }
        seed_passed = all(checks.values())
        passed = passed and seed_passed
        report["seeds"][str(seed)] = {
            "max_abs_deviation": deviation, "pearson": pearson,
            "spearman": rank_correlation, "call_disagreements": disagreements,
            "max_rank_swap_gap": swap_gap, "float32_ulp": float32_ulp,
            "model_validation": model_stats, "checks": checks, "passed": seed_passed}

        logger.info("seed %d | max|d|=%.3e pearson=%.15f spearman=%.15f "
                    "disagreements=%d | TP=%d FP=%d FN=%d | %s",
                    seed, deviation, pearson, rank_correlation, disagreements,
                    model_stats["tp"], model_stats["fp"], model_stats["fn"],
                    "PASS" if seed_passed else "FAIL")

    report["passed"] = passed
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2, default=float))
    logger.info("Wrote %s | GATE %s", args.out, "PASSED" if passed else "FAILED")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
