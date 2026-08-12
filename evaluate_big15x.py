"""Powered 15x SNP benchmark: statistical baselines, error analysis, neural arms.

This is the evaluation for the regime the audit identified as the only one in
this project with a substantial error budget: 15x SNP calling, where the
Poisson-binomial caller errs at ~32.5 loci per megabase. The 0.47 Mb test
region used previously yielded about ten errors, which is why every earlier
comparison was underpowered; here the test split is 4 Mb.

The script runs in two modes, and the first is a gate on the second:

* **baseline mode** (no ``--checkpoint``) reports the statistical callers, the
  split's variant counts, the depth distribution and the error distribution.
  If the baseline turns out to be saturated here too, the experiment stops and
  no GPU time is spent.
* **full mode** adds trained checkpoints, the paired bootstrap, error overlap,
  residual diagnostics and the constant-shift null.

Thresholds for the statistical arms are selected on the validation blocks and
frozen before the test blocks are scored. Neural arms carry the threshold their
own training run selected on those same validation blocks.
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import torch

from analysis_30x import (
    constant_shift_null,
    depth_bins_30x,
    error_overlap,
    residual_diagnostics,
)
from binomial_baseline import BinomialVariantCaller
from config import LABEL_SNP
from evaluate_binomial_baseline import pick_threshold, prf
from evaluate_quality_error import LLR_GRID, VAF_GRID, paired_bootstrap
from genomic_split import assign_window_roles, split_summary
from residual_metrics import bootstrap_f1, vaf_from_counts
from train_raw_pileup import SEQ_LEN

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

SHIFT_GRID = np.arange(-20.0, 20.01, 0.5)


def load_positions(path: str) -> np.ndarray:
    """Load the cached per-locus genomic positions.

    Args:
        path: Path to the positions ``.npy``.

    Returns:
        ``[N]`` positions.
    """
    return np.load(path)


def depth_distribution(depth: np.ndarray) -> dict:
    """Summary statistics and strata counts for a depth vector.

    Args:
        depth: ``[N]`` per-locus depth.

    Returns:
        Mean, median, extrema and per-stratum loci counts.
    """
    return {
        "mean": float(depth.mean()), "median": float(np.median(depth)),
        "min": int(depth.min()), "max": int(depth.max()),
        "strata": {label: int(mask.sum()) for label, mask in depth_bins_30x(depth)},
    }


def error_distribution(score: np.ndarray, threshold: float, truth: np.ndarray,
                       counts: np.ndarray, depth: np.ndarray) -> dict:
    """Characterize an arm's mistakes by the evidence available at those loci.

    Args:
        score: ``[N]`` decision score.
        threshold: Frozen operating point.
        truth: ``[N]`` boolean truth.
        counts: ``[N, 10]`` count matrix restricted to the same loci.
        depth: ``[N]`` depth.

    Returns:
        False-negative and false-positive profiles.
    """
    called = score >= threshold
    vaf = vaf_from_counts(counts)
    alt = np.round(vaf * depth).astype(int)

    def profile(mask: np.ndarray) -> dict:
        if not mask.any():
            return {"n": 0}
        return {
            "n": int(mask.sum()),
            "depth_median": float(np.median(depth[mask])),
            "depth_min": int(depth[mask].min()), "depth_max": int(depth[mask].max()),
            "vaf_median": float(np.median(vaf[mask])),
            "alt_reads_median": float(np.median(alt[mask])),
            "alt_reads_histogram": {str(k): int((alt[mask] == k).sum())
                                    for k in range(0, 6)},
            "alt_reads_ge6": int((alt[mask] >= 6).sum()),
        }

    return {"false_negative": profile(~called & truth),
            "false_positive": profile(called & ~truth)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", default="cache/big15x_evidence.npz")
    parser.add_argument("--counts", default="cache/big15x_counts.npz")
    parser.add_argument("--positions", default="cache/big15x_positions.npy")
    parser.add_argument("--region-start", type=int, default=32000000)
    parser.add_argument("--checkpoint", action="append", default=[],
                        metavar="NAME=PATH")
    parser.add_argument("--out", default="results/big15x/baseline.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    blob = np.load(args.evidence, allow_pickle=True)
    features, labels, llr, depth = (blob["features"], blob["labels"],
                                    blob["pb_llr"], blob["depth"])
    counts = np.load(args.counts)["counts"]
    positions = load_positions(args.positions)
    assert positions.size == labels.size, "position/locus count mismatch"

    masks = assign_window_roles(positions, args.region_start, seq_len=SEQ_LEN)
    snp = labels == LABEL_SNP
    frame = (labels == 0) | snp

    report: dict = {
        "split": split_summary(positions, labels, args.region_start, seq_len=SEQ_LEN),
        "depth": {}, "thresholds": {}, "validation": {}, "primary": {},
        "error_distribution": {}, "by_depth": {}, "diagnostics": {}, "runtime": {},
    }
    for role, mask in masks.items():
        report["depth"][role] = depth_distribution(depth[mask])
    logger.info("split: %s", {r: (v["loci"], v["snp"]) for r, v in report["split"].items()})

    validation = masks["validation"] & frame
    test = masks["test"] & frame

    # ---- statistical arms, thresholds chosen on validation blocks only -----
    binom = BinomialVariantCaller(error_rate=0.01)
    start = time.perf_counter()
    binomial_llr = binom.score_counts(counts[:, 0:4], counts[:, 9].astype(int)).llr
    binomial_seconds = time.perf_counter() - start
    scores: dict[str, np.ndarray] = {
        "frequency": vaf_from_counts(counts),
        "binomial_v1": binomial_llr,
        "poisson_binomial": llr,
    }
    grids = {"frequency": VAF_GRID}
    report["runtime"]["binomial_v1"] = {
        "scoring_seconds": binomial_seconds,
        "loci_per_second": float(labels.size / max(binomial_seconds, 1e-9))}

    thresholds: dict[str, float] = {}
    for name, score in scores.items():
        threshold, stats = pick_threshold(score[validation], snp[validation],
                                          grids.get(name, LLR_GRID))
        thresholds[name] = threshold
        report["thresholds"][name] = {"threshold": threshold,
                                      "source": "validation blocks, F1-maximizing"}
        report["validation"][name] = stats
        logger.info("validation %-18s thr=%8.3f F1=%.4f TP=%d FP=%d FN=%d",
                    name, threshold, stats["f1"], stats["tp"], stats["fp"], stats["fn"])

    # ---- neural arms -------------------------------------------------------
    residuals: dict[str, np.ndarray] = {}
    if args.checkpoint:
        from evaluate_pb_residual import model_scores
        for entry in args.checkpoint:
            name, _, path = entry.partition("=")
            score, diagnostics, seconds = model_scores(Path(path), features, llr, device)
            scores[name] = score
            thresholds[name] = diagnostics["validation_threshold"]
            residuals[name] = score - llr
            report["diagnostics"][name] = diagnostics
            report["thresholds"][name] = {"threshold": diagnostics["validation_threshold"],
                                          "source": f"validation selection inside {path}"}
            report["runtime"][name] = {
                "inference_seconds": seconds,
                "loci_per_second": float(labels.size / max(seconds, 1e-9))}

    # ---- constant-shift null ----------------------------------------------
    shift, shift_stats = constant_shift_null(
        llr[validation] - thresholds["poisson_binomial"], snp[validation], SHIFT_GRID)
    scores["pb_constant_shift"] = llr
    thresholds["pb_constant_shift"] = thresholds["poisson_binomial"] - shift
    report["thresholds"]["pb_constant_shift"] = {
        "threshold": thresholds["pb_constant_shift"], "shift": shift,
        "source": "best additive constant on validation; algebraically a moved boundary"}
    report["validation"]["pb_constant_shift"] = shift_stats

    # ---- frozen test results ----------------------------------------------
    for name, score in scores.items():
        stats = prf(score[test] >= thresholds[name], snp[test])
        stats["f1_ci95"] = bootstrap_f1(score[test] >= thresholds[name], snp[test])
        report["primary"][name] = stats

    report["oracle_diagnostic"] = {
        "_note": "test-optimal thresholds; a ceiling diagnostic, NOT a result"}
    for name, score in scores.items():
        threshold, stats = pick_threshold(score[test], snp[test],
                                          grids.get(name, LLR_GRID))
        report["oracle_diagnostic"][name] = {"threshold": threshold, **stats}

    for name in ("poisson_binomial", "binomial_v1"):
        report["error_distribution"][name] = error_distribution(
            scores[name][test], thresholds[name], snp[test], counts[test], depth[test])

    for label, bin_mask in depth_bins_30x(depth):
        combined = bin_mask & test
        entry = {"loci": int(combined.sum()), "snp": int(snp[combined].sum())}
        for name, score in scores.items():
            entry[name] = prf(score[combined] >= thresholds[name], snp[combined])
        report["by_depth"][label] = entry

    baseline_calls = scores["poisson_binomial"][test] >= thresholds["poisson_binomial"]
    report["bootstrap_vs_poisson_binomial"] = {
        name: paired_bootstrap(scores[name][test] >= thresholds[name],
                               baseline_calls, snp[test])
        for name in scores if name != "poisson_binomial"}
    report["error_overlap"] = {
        name: error_overlap(baseline_calls, scores[name][test] >= thresholds[name],
                            snp[test], "pb", name)
        for name in scores if name != "poisson_binomial"}

    feature_names = list(blob["feature_names"])
    extra_index = {key: feature_names.index(key)
                   for key in ("alt_mean_phred", "alt_strand_bias", "alt_mean_mapq")
                   if key in feature_names}
    report["residual"] = {}
    for name, residual in residuals.items():
        report["residual"][name] = residual_diagnostics(
            residual[test], snp[test], depth[test], vaf_from_counts(counts)[test],
            llr[test], extra={k: features[test][:, i] for k, i in extra_index.items()})
        own = float(residual[test].mean())
        report["residual"][name]["pb_at_own_mean_residual"] = {
            "shift": own, "threshold": thresholds["poisson_binomial"] - own,
            **prf(llr[test] >= thresholds["poisson_binomial"] - own, snp[test])}

    report["dataset"] = {"loci": int(labels.size), "test_frame_loci": int(test.sum()),
                         "test_snp": int(snp[test].sum()),
                         "validation_snp": int(snp[validation].sum())}

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2, default=float))
    logger.info("Wrote %s", args.out)

    print("\n=== SPLIT ===")
    for role, entry in report["split"].items():
        print(f"  {role:<11} {entry['loci']:>9,} loci ({entry['megabases']:.2f} Mb)  "
              f"SNP={entry['snp']:>6}  INS={entry['insertion']:>4}  DEL={entry['deletion']:>4}  "
              f"blocks={entry['blocks']}")
    print("\n=== PRIMARY (15x powered): SNP vs Normal, validation-frozen thresholds ===")
    for name in scores:
        s = report["primary"][name]
        low, high = s["f1_ci95"]
        print(f"  {name:<22} thr={thresholds[name]:>8.3f} P={s['precision']:.4f} "
              f"R={s['recall']:.4f} F1={s['f1']:.4f} [{low:.4f},{high:.4f}] "
              f"TP={s['tp']} FP={s['fp']} FN={s['fn']}")
    if len(scores) > 4:
        print("\n=== Paired bootstrap vs Poisson-binomial (dF1, 95% CI) ===")
        for name, s in report["bootstrap_vs_poisson_binomial"].items():
            low, high = s["delta_f1_ci95"]
            print(f"  {name:<22} dF1={s['delta_f1']:+.4f} [{low:+.4f}, {high:+.4f}] "
                  f"dTP={s['delta_tp']:+d} dFP={s['delta_fp']:+d} dFN={s['delta_fn']:+d}")


if __name__ == "__main__":
    main()
