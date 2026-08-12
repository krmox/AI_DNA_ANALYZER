"""Frozen 30x test evaluation: Poisson-binomial versus the PB-residual Mamba.

This script is written and committed *before* any 30x test number is read, and
it is run exactly once per frozen configuration. It reports; it does not
select. Every threshold it applies was fixed on validation data beforehand:

* the statistical arms' thresholds are re-selected here on the **30x
  validation split** -- the last 10% of the training region's window batches,
  the same split the neural model trains against -- because a threshold frozen
  at 15x is not a validated threshold at 30x. Test labels are never consulted
  for this.
* each neural arm carries the threshold its own training run selected on that
  same validation split, stored inside its checkpoint.

Beyond the 15x report this adds the four things the 30x question actually
turns on (see :mod:`analysis_30x`): 30x-resolving depth strata, error-set
overlap, the constant-shift null, and residual diagnostics.

The residual is recovered exactly, not re-derived: the model computes
``final_logit = PB_logit + residual`` and the decision score is the SNP-Normal
logit gap, so ``residual_gap == model_score - pb_llr`` identically. No separate
forward pass and no extra assumption is involved.

Oracle thresholds are computed and reported under ``oracle_diagnostic`` as a
ceiling. They are never substituted for the frozen operating point.
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
from evaluate_pb_residual import model_scores
from evaluate_quality_error import LLR_GRID, VAF_GRID, paired_bootstrap
from residual_metrics import vaf_from_counts
from train_raw_pileup import BATCH_SIZE, SEQ_LEN, TRAIN_FRACTION

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

#: Constants searched for the constant-shift null, on the LLR axis.
SHIFT_GRID = np.arange(-20.0, 20.01, 0.5)


def validation_slice(n_loci: int) -> slice:
    """The validation loci of the training region.

    Byte-identical to the split in ``train_pb_residual`` and
    ``train_residual_binomial``: whole window-batches, sequential, no shuffle,
    so a locus is either fit or validation in every arm of the experiment.

    Args:
        n_loci: Loci in the training region.

    Returns:
        Slice selecting the validation loci.
    """
    n_windows = n_loci // SEQ_LEN
    n_batches = (n_windows + BATCH_SIZE - 1) // BATCH_SIZE
    split_locus = int(n_batches * TRAIN_FRACTION) * BATCH_SIZE * SEQ_LEN
    return slice(split_locus, n_loci)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-evidence", default="cache/train_30x_evidence.npz")
    parser.add_argument("--train-counts", default="cache/train_30x_counts.npz")
    parser.add_argument("--test-evidence", default="cache/test_30x_evidence.npz")
    parser.add_argument("--test-counts", default="cache/test_30x_counts.npz")
    parser.add_argument("--checkpoint", action="append", default=[],
                        metavar="NAME=PATH", help="Named checkpoint to evaluate.")
    parser.add_argument("--primary", default="pb_residual_seedA",
                        help="Arm compared against PB in the headline and overlap.")
    parser.add_argument("--out", default="results/30x_pb_residual/test_results.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # ---------- validation split: thresholds are chosen here, not on test ----
    train_blob = np.load(args.train_evidence, allow_pickle=True)
    train_counts = np.load(args.train_counts)["counts"]
    val = validation_slice(train_blob["labels"].size)
    val_labels = train_blob["labels"][val]
    val_llr = train_blob["pb_llr"][val]
    val_counts = train_counts[val]
    val_snp = val_labels == LABEL_SNP
    val_frame = (val_labels == 0) | val_snp

    val_binom_v1 = BinomialVariantCaller(error_rate=0.01).score_counts(
        val_counts[:, 0:4], val_counts[:, 9].astype(int)).llr
    val_frequency = vaf_from_counts(val_counts)

    validation_scores = {
        "frequency": (val_frequency, VAF_GRID),
        "binomial_v1": (val_binom_v1, LLR_GRID),
        "poisson_binomial": (val_llr, LLR_GRID),
    }

    report: dict = {
        "split": {
            "validation_loci": int(val_labels.size),
            "validation_frame_loci": int(val_frame.sum()),
            "validation_snp": int(val_snp[val_frame].sum()),
            "rule": "last 10% of training-region window batches; test labels unused",
        },
        "thresholds": {}, "validation": {}, "primary": {}, "by_depth": {},
        "diagnostics": {}, "runtime": {},
    }

    thresholds: dict[str, float] = {}
    for name, (score, grid) in validation_scores.items():
        threshold, stats = pick_threshold(score[val_frame], val_snp[val_frame], grid)
        thresholds[name] = threshold
        report["thresholds"][name] = {"threshold": threshold,
                                      "source": "30x validation split, F1-maximizing"}
        report["validation"][name] = stats
        logger.info("validation %-20s thr=%8.3f F1=%.4f (TP=%d FP=%d FN=%d)",
                    name, threshold, stats["f1"], stats["tp"], stats["fp"], stats["fn"])

    # ---------- test arrays ------------------------------------------------
    blob = np.load(args.test_evidence, allow_pickle=True)
    features, labels, llr, depth = (blob["features"], blob["labels"],
                                    blob["pb_llr"], blob["depth"])
    counts = np.load(args.test_counts)["counts"]
    snp = labels == LABEL_SNP
    frame = (labels == 0) | snp

    start = time.perf_counter()
    binom_v1 = BinomialVariantCaller(error_rate=0.01).score_counts(
        counts[:, 0:4], counts[:, 9].astype(int)).llr
    binomial_seconds = time.perf_counter() - start

    scores: dict[str, np.ndarray] = {
        "frequency": vaf_from_counts(counts),
        "binomial_v1": binom_v1,
        "poisson_binomial": llr,
    }
    report["runtime"]["binomial_v1"] = {
        "scoring_seconds": binomial_seconds,
        "loci_per_second": float(labels.size / max(binomial_seconds, 1e-9))}

    # ---------- neural arms, each with its own validation-frozen threshold --
    residuals: dict[str, np.ndarray] = {}
    for entry in args.checkpoint:
        name, _, path = entry.partition("=")
        score, diagnostics, seconds = model_scores(Path(path), features, llr, device)
        scores[name] = score
        thresholds[name] = diagnostics["validation_threshold"]
        # Exact by construction: final_logit = PB_logit + residual.
        residuals[name] = score - llr
        report["diagnostics"][name] = diagnostics
        report["thresholds"][name] = {"threshold": diagnostics["validation_threshold"],
                                      "source": f"validation selection inside {path}"}
        report["runtime"][name] = {"inference_seconds": seconds,
                                   "loci_per_second": float(labels.size / max(seconds, 1e-9))}

    # ---------- constant-shift null ----------------------------------------
    # Fitted on validation, applied to test, held to the same standard as the
    # neural arms. Reported per neural arm using that arm's own mean residual,
    # plus a free-parameter best shift.
    shift, shift_stats = constant_shift_null(
        val_llr[val_frame] - thresholds["poisson_binomial"], val_snp[val_frame], SHIFT_GRID)
    scores["pb_constant_shift"] = llr
    thresholds["pb_constant_shift"] = thresholds["poisson_binomial"] - shift
    report["thresholds"]["pb_constant_shift"] = {
        "threshold": thresholds["pb_constant_shift"],
        "source": "best additive constant on validation; algebraically a moved boundary",
        "shift": shift}
    report["validation"]["pb_constant_shift"] = shift_stats

    # ---------- primary, frozen-threshold results --------------------------
    for name, score in scores.items():
        report["primary"][name] = prf(score[frame] >= thresholds[name], snp[frame])

    report["oracle_diagnostic"] = {
        "_note": "test-optimal thresholds; a ceiling diagnostic, NOT a result"}
    for name, score in scores.items():
        grid = VAF_GRID if name == "frequency" else LLR_GRID
        threshold, stats = pick_threshold(score[frame], snp[frame], grid)
        report["oracle_diagnostic"][name] = {"threshold": threshold, **stats}

    # ---------- depth strata -----------------------------------------------
    for label, bin_mask in depth_bins_30x(depth):
        combined = bin_mask & frame
        entry = {"loci": int(combined.sum()), "snp": int(snp[combined].sum())}
        for name, score in scores.items():
            entry[name] = prf(score[combined] >= thresholds[name], snp[combined])
        report["by_depth"][label] = entry

    # ---------- paired bootstrap vs the statistical baseline ---------------
    baseline_calls = scores["poisson_binomial"][frame] >= thresholds["poisson_binomial"]
    report["bootstrap_vs_poisson_binomial"] = {
        name: paired_bootstrap(scores[name][frame] >= thresholds[name],
                               baseline_calls, snp[frame])
        for name in scores if name != "poisson_binomial"}

    # ---------- error overlap and residual diagnostics ---------------------
    report["error_overlap"] = {}
    for name in scores:
        if name == "poisson_binomial":
            continue
        report["error_overlap"][name] = error_overlap(
            baseline_calls, scores[name][frame] >= thresholds[name], snp[frame],
            "pb", name)

    feature_names = list(blob["feature_names"])
    extra_index = {key: feature_names.index(key)
                   for key in ("alt_mean_phred", "alt_strand_bias", "alt_mean_mapq")
                   if key in feature_names}
    report["residual"] = {}
    for name, residual in residuals.items():
        report["residual"][name] = residual_diagnostics(
            residual[frame], snp[frame], depth[frame],
            vaf_from_counts(counts)[frame], llr[frame],
            extra={key: features[frame][:, index] for key, index in extra_index.items()})
        # The arm's own mean residual, as a boundary move on PB.
        own = float(residual[frame].mean())
        moved = prf(llr[frame] >= thresholds["poisson_binomial"] - own, snp[frame])
        report["residual"][name]["pb_at_own_mean_residual"] = {
            "shift": own, "threshold": thresholds["poisson_binomial"] - own, **moved}

    report["dataset"] = {
        "loci": int(labels.size), "frame_loci": int(frame.sum()),
        "snp": int(snp[frame].sum()), "mean_depth": float(depth.mean()),
        "median_depth": float(np.median(depth))}

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2, default=float))
    logger.info("Wrote %s", args.out)

    print("\n=== PRIMARY (30x): SNP vs Normal, validation-frozen thresholds ===")
    for name in scores:
        stats = report["primary"][name]
        print(f"  {name:<28} thr={thresholds[name]:>8.3f} P={stats['precision']:.4f} "
              f"R={stats['recall']:.4f} F1={stats['f1']:.4f} "
              f"TP={stats['tp']} FP={stats['fp']} FN={stats['fn']}")
    print("\n=== Paired bootstrap vs Poisson-binomial (dF1, 95% CI) ===")
    for name, stats in report["bootstrap_vs_poisson_binomial"].items():
        low, high = stats["delta_f1_ci95"]
        print(f"  {name:<28} dF1={stats['delta_f1']:+.4f} [{low:+.4f}, {high:+.4f}] "
              f"dTP={stats['delta_tp']:+d} dFP={stats['delta_fp']:+d} "
              f"dFN={stats['delta_fn']:+d}")
    print("\n=== Error overlap vs Poisson-binomial ===")
    for name, stats in report["error_overlap"].items():
        print(f"  {name:<28} shared_FN={stats['shared_fn']} "
              f"PB_FN_recovered={stats['pb_fn_recovered_by_' + name]} "
              f"new_FP={stats[name + '_fp_not_made_by_pb']} "
              f"disagreements={stats['disagreements']}")


if __name__ == "__main__":
    main()
