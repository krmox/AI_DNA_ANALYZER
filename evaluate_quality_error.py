"""Experiment: does per-locus quality-derived epsilon beat fixed epsilon = 0.01?

Protocol
--------
Every arm is scored on exactly the same loci from the same cached tensors.
Thresholds -- and the choice of which new arm is the headline -- are fixed on
the validation split (the last 10% of the training region's windows, the same
split ``evaluate_binomial_baseline.py`` and the neural training used) before
the 15x test region is scored once.

The pre-registered primary new arm is **Estimator C (reference-only mean)**.
It is named in advance rather than chosen by validation F1, because it is the
one estimator with an a-priori argument behind it: it cannot be driven by the
ALT reads, so a gain from it is attributable to better uncertainty estimation
and not to a re-tuned decision boundary. A validation-selected arm is reported
alongside it, and both are frozen before the test region is touched.

Test labels are used for exactly one thing: computing the final numbers after
everything above is frozen. Oracle (test-optimal) thresholds are computed as a
labelled diagnostic and never reported as a result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import subprocess
import time
from pathlib import Path

import numpy as np

from binomial_baseline import BinomialVariantCaller
from config import LABEL_SNP
from evaluate_binomial_baseline import depth_bins, pick_threshold, prf
from quality_error_model import (
    EPSILON_CEILING,
    EPSILON_FLOOR,
    FIXED_EPSILON,
    binomial_llr,
    candidate_alt,
    estimate_epsilon,
    extract_quality_evidence,
    poisson_binomial_llr,
)
from read_level_pileup import reconstruct_counts

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

TRAIN_FRACTION = 0.9   # matches run_train.py / evaluate_binomial_baseline.py
SEQ_LEN = 64
BOOTSTRAP_RESAMPLES = 10000
BOOTSTRAP_SEED = 20260812

#: Frozen before the test region is scored. See the module docstring.
PRIMARY_NEW_ARM = "binomial_ref_mean"

#: Arms driven by a per-locus epsilon plug-in, in report order.
QUALITY_ARMS: tuple[tuple[str, str], ...] = (
    ("binomial_mean_q", "mean"),
    ("binomial_median_q", "median"),
    ("binomial_ref_mean", "ref_mean"),
    ("binomial_ref_trimmed", "ref_trimmed"),
)

LLR_GRID = np.arange(-20.0, 200.0, 0.5)
VAF_GRID = np.arange(0.01, 1.00, 0.005)


def file_digest(path: str | Path) -> str:
    """SHA-256 of a file, for the reproducibility record."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def array_digest(array: np.ndarray) -> str:
    """SHA-256 of an array's bytes, so every arm can prove it saw the same input."""
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def load_region(counts_npz: str, reads_npz: str, locus_slice: slice | None = None):
    """Load a region's cached counts, read tensor and labels, and cross-check them.

    The two caches were written by independent extraction passes. Rebuilding
    the aggregate counts from the read tensor and demanding exact equality is
    the control that makes every comparison in this script valid: if it holds,
    the quality-derived arms and the fixed-epsilon baseline are provably
    describing the same reads.

    Args:
        counts_npz: Path to the cached ``[N, 10]`` count matrix.
        reads_npz: Path to the cached ``[N, R, 8]`` read tensor.
        locus_slice: Optional slice restricting the loci returned.

    Returns:
        ``(counts, reads, labels)``.

    Raises:
        AssertionError: If the two representations disagree anywhere.
    """
    count_data = np.load(counts_npz)
    read_data = np.load(reads_npz)
    counts, count_labels = count_data["counts"], count_data["labels"]
    reads, read_labels = read_data["reads"], read_data["labels"]

    assert counts.shape[0] == reads.shape[0], "count/read locus-count mismatch"
    assert np.array_equal(count_labels, read_labels), "count/read label mismatch"

    if locus_slice is not None:
        counts, reads, count_labels = (counts[locus_slice], reads[locus_slice],
                                       count_labels[locus_slice])

    rebuilt = reconstruct_counts(reads, counts[:, 9])
    for column in (0, 1, 2, 3, 4, 6, 7):
        assert np.array_equal(rebuilt[:, column], counts[:, column]), \
            f"read tensor does not reproduce count column {column}"
    return counts, reads, count_labels


def pooled_error_rate(evidence, floor: float = EPSILON_FLOOR,
                      ceiling: float = EPSILON_CEILING) -> float:
    """One global error probability: the pooled mean over every counted read.

    This is the calibration control, and it is the most important number in the
    experiment. Binomial v1 assumes 0.01 (~Q20); if the reads are actually
    ~Q30, then *any* smaller constant would improve the caller and a
    per-locus estimator would win for reasons that have nothing to do with
    locus-specific uncertainty. Freezing this constant from the validation
    region separates "0.01 was the wrong number" from "the error rate varies
    per locus and knowing that helps".

    Label-free: it pools base qualities only.

    Args:
        evidence: Validation-region quality evidence.
        floor: Lower clamp, as for the per-locus estimators.
        ceiling: Upper clamp.

    Returns:
        The pooled mean error probability.
    """
    selected = evidence.counted
    pooled = float(np.sum(np.where(selected, evidence.error_probability, 0.0))
                   / max(int(selected.sum()), 1))
    return float(np.clip(pooled, floor, ceiling))


def score_all_arms(counts: np.ndarray, reads: np.ndarray,
                   global_epsilon: float | None = None,
                   shuffle_seed: int = BOOTSTRAP_SEED) -> tuple[dict, dict]:
    """Compute every arm's per-locus score from counts and the read tensor.

    Args:
        counts: ``[N, 10]`` count matrix (``pileup_counts.COUNT_COLUMNS``).
        reads: ``[N, R, 8]`` read tensor.
        global_epsilon: Calibrated global constant for the calibration-control
            arm. When None the arm is skipped (validation pass, where the
            constant is still being estimated).
        shuffle_seed: Seed for the locus-shuffled control arm.

    Returns:
        ``(scores, diagnostics)``; ``scores`` maps arm name to a ``[N]`` array.
    """
    base = counts[:, 0:4]
    reference_index = counts[:, 9].astype(int)

    start = time.perf_counter()
    evidence = extract_quality_evidence(reads, reference_index)
    evidence_seconds = time.perf_counter() - start

    alt_index, k, n = candidate_alt(base, reference_index)

    # Control: the ALT candidate must never be the reference base.
    has_reference = (reference_index >= 1) & (reference_index <= 4)
    assert not np.any(alt_index[has_reference] == reference_index[has_reference] - 1), \
        "candidate ALT collided with the reference base"

    scores: dict[str, np.ndarray] = {}
    diagnostics: dict[str, dict] = {"evidence_seconds": evidence_seconds}

    with np.errstate(invalid="ignore", divide="ignore"):
        scores["frequency"] = np.where(n > 0, k / np.maximum(n, 1), 0.0)

    start = time.perf_counter()
    epsilon_fixed = np.full(n.shape, FIXED_EPSILON)
    scores["binomial_v1"] = binomial_llr(k, n, epsilon_fixed)["llr"]
    fixed_seconds = time.perf_counter() - start

    epsilons = {"binomial_v1": epsilon_fixed}
    for arm, method in QUALITY_ARMS:
        start = time.perf_counter()
        epsilon, estimator_diagnostics = estimate_epsilon(evidence, method)
        scores[arm] = binomial_llr(k, n, epsilon)["llr"]
        estimator_diagnostics["seconds"] = time.perf_counter() - start
        diagnostics[arm] = estimator_diagnostics
        epsilons[arm] = epsilon

    # Control 1 -- calibration. A single global constant, estimated from the
    # validation region's base qualities. Isolates "0.01 was miscalibrated"
    # from "epsilon genuinely varies per locus".
    if global_epsilon is not None:
        epsilon_global = np.full(n.shape, float(global_epsilon))
        scores["binomial_fixed_calibrated"] = binomial_llr(k, n, epsilon_global)["llr"]
        epsilons["binomial_fixed_calibrated"] = epsilon_global
        diagnostics["binomial_fixed_calibrated"] = {
            "estimator": "global pooled mean, frozen from validation",
            "epsilon": float(global_epsilon), "seconds": 0.0}

    # Control 2 -- locus-specificity. The same epsilon values, randomly
    # reassigned to loci. The marginal distribution of epsilon is preserved
    # exactly; only the pairing with the locus is destroyed. If this scores
    # like the real arm, the gain is a scale effect, not information.
    shuffled = np.random.default_rng(shuffle_seed).permutation(epsilons["binomial_mean_q"])
    scores["binomial_mean_q_shuffled"] = binomial_llr(k, n, shuffled)["llr"]
    epsilons["binomial_mean_q_shuffled"] = shuffled
    diagnostics["binomial_mean_q_shuffled"] = {
        "estimator": "mean-Q epsilon permuted across loci", "seed": shuffle_seed,
        "seconds": 0.0}

    start = time.perf_counter()
    scores["poisson_binomial"] = poisson_binomial_llr(
        evidence, k.astype(np.int64), alt_index=alt_index)["llr"]
    diagnostics["poisson_binomial"] = {"seconds": time.perf_counter() - start}

    diagnostics["binomial_v1"] = {"seconds": fixed_seconds, "estimator": "fixed"}

    for name, score in scores.items():
        assert np.all(np.isfinite(score)), f"non-finite score in arm {name}"

    return scores, {"per_arm": diagnostics, "epsilon": epsilons,
                    "k": k, "n": n, "alt_index": alt_index, "evidence": evidence}


def grid_for(arm: str) -> np.ndarray:
    """Threshold search grid appropriate to an arm's score scale."""
    return VAF_GRID if arm == "frequency" else LLR_GRID


def epsilon_summary(epsilon: np.ndarray, mask: np.ndarray) -> dict:
    """Distribution summary of an epsilon array over a locus subset."""
    values = epsilon[mask]
    if values.size == 0:
        return {"loci": 0}
    percentiles = np.percentile(values, [1, 5, 25, 50, 75, 95, 99])
    return {
        "loci": int(values.size), "mean": float(values.mean()),
        "min": float(values.min()), "max": float(values.max()),
        "p1": float(percentiles[0]), "p5": float(percentiles[1]),
        "p25": float(percentiles[2]), "median": float(percentiles[3]),
        "p75": float(percentiles[4]), "p95": float(percentiles[5]),
        "p99": float(percentiles[6]),
        "implied_mean_phred": float(-10.0 * np.log10(values.mean())),
    }


def paired_bootstrap(predicted_a: np.ndarray, predicted_b: np.ndarray, truth: np.ndarray,
                     resamples: int = BOOTSTRAP_RESAMPLES, seed: int = BOOTSTRAP_SEED) -> dict:
    """Paired bootstrap CI for A - B on F1, TP, FP and FN.

    Loci are resampled with replacement; both arms are scored on the *same*
    resample, so the interval is for the paired difference and correctly
    absorbs the fact that the two arms agree on almost every locus.

    Only loci where some arm predicts positive or the truth is positive can
    contribute to TP/FP/FN. The rest are pooled into a single multinomial cell,
    which makes 10,000 resamples over 455k loci exact and cheap rather than
    approximate and slow.

    Args:
        predicted_a: Boolean calls from the arm under test.
        predicted_b: Boolean calls from the reference arm.
        truth: Boolean truth.
        resamples: Bootstrap resamples.
        seed: RNG seed.

    Returns:
        Point differences and 95% percentile intervals.
    """
    interesting = predicted_a | predicted_b | truth
    index = np.nonzero(interesting)[0]
    total = truth.size
    inert = total - index.size

    def counts(predicted: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return (predicted[index] & truth[index], predicted[index] & ~truth[index],
                ~predicted[index] & truth[index])

    tp_a, fp_a, fn_a = counts(predicted_a)
    tp_b, fp_b, fn_b = counts(predicted_b)

    probabilities = np.full(index.size + 1, 1.0 / total)
    probabilities[-1] = inert / total

    rng = np.random.default_rng(seed)
    draws = rng.multinomial(total, probabilities, size=resamples)[:, :-1].astype(np.float64)

    def f1_from(tp_mask, fp_mask, fn_mask) -> np.ndarray:
        tp = draws @ tp_mask
        fp = draws @ fp_mask
        fn = draws @ fn_mask
        return 2 * tp / np.maximum(2 * tp + fp + fn, 1e-12), tp, fp, fn

    f1_a, tp_a_boot, fp_a_boot, fn_a_boot = f1_from(tp_a, fp_a, fn_a)
    f1_b, tp_b_boot, fp_b_boot, fn_b_boot = f1_from(tp_b, fp_b, fn_b)

    def interval(values: np.ndarray) -> list[float]:
        low, high = np.percentile(values, [2.5, 97.5])
        return [float(low), float(high)]

    point = prf(predicted_a, truth), prf(predicted_b, truth)
    delta_f1 = f1_a - f1_b
    return {
        "resamples": resamples, "seed": seed,
        "delta_f1": point[0]["f1"] - point[1]["f1"],
        "delta_f1_ci95": interval(delta_f1),
        "delta_f1_p_two_sided": float(2 * min((delta_f1 <= 0).mean(), (delta_f1 >= 0).mean())),
        "delta_tp": point[0]["tp"] - point[1]["tp"],
        "delta_tp_ci95": interval(tp_a_boot - tp_b_boot),
        "delta_fp": point[0]["fp"] - point[1]["fp"],
        "delta_fp_ci95": interval(fp_a_boot - fp_b_boot),
        "delta_fn": point[0]["fn"] - point[1]["fn"],
        "delta_fn_ci95": interval(fn_a_boot - fn_b_boot),
    }


def disagreement_analysis(scores: dict, thresholds: dict, extras: dict, counts: np.ndarray,
                          truth: np.ndarray, mask: np.ndarray, arm: str,
                          limit: int = 25) -> dict:
    """Where the new arm and Binomial v1 disagree, and what those loci look like.

    Args:
        scores: Arm name -> per-locus score.
        thresholds: Arm name -> frozen threshold.
        extras: Output of :func:`score_all_arms`.
        counts: ``[N, 10]`` count matrix.
        truth: Boolean SNP truth.
        mask: Evaluation-frame mask (SNP vs Normal).
        arm: The arm to contrast against ``binomial_v1``.
        limit: Loci listed per category.

    Returns:
        Category counts plus per-locus detail for a bounded sample.
    """
    new = scores[arm] >= thresholds[arm]
    old = scores["binomial_v1"] >= thresholds["binomial_v1"]
    depth = counts[:, 6]

    categories = {
        "rescued_fn": mask & truth & new & ~old,       # v1 missed it, new finds it
        "lost_tp": mask & truth & ~new & old,          # v1 found it, new misses it
        "new_fp": mask & ~truth & new & ~old,          # new invents a false positive
        "removed_fp": mask & ~truth & ~new & old,      # new removes a false positive
    }

    report: dict = {}
    for name, selection in categories.items():
        indices = np.nonzero(selection)[0]
        report[name] = {
            "count": int(indices.size),
            "loci": [
                {
                    "index": int(i), "depth": float(depth[i]),
                    "n": float(extras["n"][i]), "k": float(extras["k"][i]),
                    "vaf": float(scores["frequency"][i]),
                    "mean_base_quality": float(counts[i, 7] / max(counts[i, 0:4].sum(), 1)),
                    "epsilon_v1": float(extras["epsilon"]["binomial_v1"][i]),
                    "epsilon_new": float(extras["epsilon"][arm][i]) if arm in
                    extras["epsilon"] else None,
                    "llr_v1": float(scores["binomial_v1"][i]),
                    "llr_new": float(scores[arm][i]),
                }
                for i in indices[:limit]
            ],
        }
        if indices.size:
            report[name]["depth_median"] = float(np.median(depth[indices]))
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-counts", default="cache/train_15x_counts.npz")
    parser.add_argument("--train-reads",
                        default="data/giab_hg002_15x_readlevel/train_15x_reads.npz")
    parser.add_argument("--test-counts", default="cache/test_15x_counts.npz")
    parser.add_argument("--test-reads",
                        default="data/giab_hg002_15x_readlevel/test_15x_reads.npz")
    parser.add_argument("--out", default="results/quality_error_15x_results.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report: dict = {
        "protocol": {
            "primary_new_arm": PRIMARY_NEW_ARM,
            "primary_new_arm_source": "pre-registered before test evaluation; chosen for "
                                      "its anti-circularity argument, not by validation F1",
            "threshold_source": "validation split = last 10% of training-region windows",
            "test_label_use": "final evaluation only",
            "epsilon_floor": EPSILON_FLOOR, "epsilon_ceiling": EPSILON_CEILING,
            "fixed_epsilon": FIXED_EPSILON,
            "bootstrap_resamples": BOOTSTRAP_RESAMPLES, "bootstrap_seed": BOOTSTRAP_SEED,
        },
        "inputs": {}, "validation": {}, "thresholds": {}, "primary": {},
        "by_depth": {}, "runtime": {}, "diagnostics": {},
    }

    try:
        report["reproducibility"] = {
            "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                  text=True).strip(),
            "git_status": subprocess.check_output(["git", "status", "--porcelain"],
                                                  text=True).strip().splitlines(),
        }
    except (subprocess.CalledProcessError, FileNotFoundError):
        report["reproducibility"] = {"git_commit": None}

    for name, path in (("train_counts", args.train_counts), ("train_reads", args.train_reads),
                       ("test_counts", args.test_counts), ("test_reads", args.test_reads)):
        report["inputs"][name] = {"path": path, "sha256": file_digest(path)}

    # ---------- validation: last 10% of the training region's windows ----------
    total_loci = int(np.load(args.train_counts)["labels"].shape[0])
    total_windows = total_loci // SEQ_LEN
    split_window = int(total_windows * TRAIN_FRACTION)
    val_slice = slice(split_window * SEQ_LEN, total_windows * SEQ_LEN)
    logger.info("Training region: %d loci / %d windows; validation = windows %d..%d",
                total_loci, total_windows, split_window, total_windows - 1)

    val_counts, val_reads, val_labels = load_region(args.train_counts, args.train_reads,
                                                    val_slice)
    val_scores, val_extras = score_all_arms(val_counts, val_reads)

    # The calibration constant is estimated here, on validation reads only, and
    # then frozen. It is label-free, so estimating it costs no test isolation.
    global_epsilon = pooled_error_rate(val_extras["evidence"])
    logger.info("Calibrated global epsilon (validation reads, label-free): %.3e (~Q%.1f)",
                global_epsilon, -10 * np.log10(global_epsilon))
    val_scores["binomial_fixed_calibrated"] = binomial_llr(
        val_extras["k"], val_extras["n"], np.full(val_extras["n"].shape, global_epsilon))["llr"]
    val_extras["epsilon"]["binomial_fixed_calibrated"] = np.full(val_extras["n"].shape,
                                                                 global_epsilon)
    report["protocol"]["calibrated_global_epsilon"] = global_epsilon

    val_snp = val_labels == LABEL_SNP
    val_frame = (val_labels == 0) | val_snp

    report["validation"] = {
        "loci": int(val_labels.size), "frame_loci": int(val_frame.sum()),
        "snp": int(val_snp[val_frame].sum()),
        "window_range": [split_window, total_windows - 1],
        "input_digest": array_digest(val_reads),
    }

    for arm in val_scores:
        threshold, stats = pick_threshold(val_scores[arm][val_frame], val_snp[val_frame],
                                          grid_for(arm))
        report["thresholds"][arm] = {"threshold": threshold, "validation": stats,
                                     "source": "validation split, test labels unused"}
        logger.info("  %-22s thr=%8.3f  val F1=%.4f (P=%.4f R=%.4f)", arm, threshold,
                    stats["f1"], stats["precision"], stats["recall"])

    thresholds = {arm: report["thresholds"][arm]["threshold"] for arm in val_scores}

    # Control 3 -- the strongest form of the calibration objection. Give the
    # fixed-epsilon caller a joint (epsilon, threshold) search on validation,
    # i.e. its best possible global constant, and carry that to test. If a
    # per-locus estimator still wins, no choice of constant explains the gain.
    best_f1, best_global_epsilon, best_global_threshold = -1.0, None, None
    sweep = []
    for candidate in 10 ** np.arange(-5.0, -1.4, 0.1):
        score = binomial_llr(val_extras["k"], val_extras["n"],
                             np.full(val_extras["n"].shape, candidate))["llr"]
        threshold, stats = pick_threshold(score[val_frame], val_snp[val_frame], LLR_GRID)
        sweep.append({"epsilon": float(candidate), "threshold": threshold,
                      "validation_f1": stats["f1"]})
        if stats["f1"] > best_f1:
            best_f1, best_global_epsilon, best_global_threshold = (
                stats["f1"], float(candidate), threshold)
    report["diagnostics"]["global_epsilon_sweep"] = {
        "note": "fixed-epsilon caller given a joint (epsilon, threshold) search on "
                "validation; the ceiling of the whole fixed-epsilon family",
        "best_epsilon": best_global_epsilon, "best_threshold": best_global_threshold,
        "best_validation_f1": best_f1, "sweep": sweep,
    }
    logger.info("Best GLOBAL epsilon on validation: %.3e (~Q%.1f) thr=%.1f valF1=%.4f",
                best_global_epsilon, -10 * np.log10(best_global_epsilon),
                best_global_threshold, best_f1)

    # Diagnostic arm: Binomial v1 at the threshold a previous run found better on
    # test. Reported explicitly as a test-informed diagnostic, never as a result.
    report["thresholds"]["binomial_v1_thr9.5"] = {
        "threshold": 9.5,
        "source": "DIAGNOSTIC ONLY -- known from a previous test-region observation; "
                  "not validation-derived and not a primary result",
    }

    # Validation-selected new arm, frozen here, before any test label is read.
    quality_arm_names = [arm for arm, _ in QUALITY_ARMS] + ["poisson_binomial"]
    validation_selected = max(quality_arm_names,
                              key=lambda arm: report["thresholds"][arm]["validation"]["f1"])
    report["protocol"]["validation_selected_arm"] = validation_selected
    logger.info("Pre-registered primary new arm: %s | validation-selected arm: %s",
                PRIMARY_NEW_ARM, validation_selected)

    report["diagnostics"]["validation_epsilon"] = {
        arm: {
            "all": epsilon_summary(val_extras["epsilon"][arm], val_frame),
            "normal": epsilon_summary(val_extras["epsilon"][arm], val_frame & ~val_snp),
            "snp": epsilon_summary(val_extras["epsilon"][arm], val_frame & val_snp),
        }
        for arm in val_extras["epsilon"]
    }
    report["diagnostics"]["estimators"] = {
        arm: val_extras["per_arm"][arm] for arm in val_extras["per_arm"]
        if isinstance(val_extras["per_arm"][arm], dict)
    }

    del val_reads, val_extras, val_scores

    # ---------- held-out test region, touched once ----------
    logger.info("Scoring the TEST region (thresholds now frozen)...")
    load_start = time.perf_counter()
    counts, reads, labels = load_region(args.test_counts, args.test_reads)
    load_seconds = time.perf_counter() - load_start

    score_start = time.perf_counter()
    scores, extras = score_all_arms(counts, reads, global_epsilon=global_epsilon)
    score_seconds = time.perf_counter() - score_start

    scores["binomial_best_global_eps"] = binomial_llr(
        extras["k"], extras["n"],
        np.full(extras["n"].shape, best_global_epsilon))["llr"]
    thresholds["binomial_best_global_eps"] = best_global_threshold
    report["thresholds"]["binomial_best_global_eps"] = {
        "threshold": best_global_threshold, "epsilon": best_global_epsilon,
        "source": "joint (epsilon, threshold) search on validation",
    }

    scores["binomial_v1_thr9.5"] = scores["binomial_v1"]
    thresholds["binomial_v1_thr9.5"] = 9.5

    snp = labels == LABEL_SNP
    frame = (labels == 0) | snp
    depth = counts[:, 6]

    report["dataset"] = {
        "test_loci": int(labels.size), "frame_loci": int(frame.sum()),
        "snp": int(snp[frame].sum()), "normal": int((labels == 0).sum()),
        "insertion": int((labels == 2).sum()), "deletion": int((labels == 3).sum()),
        "depth_mean": float(depth.mean()), "depth_median": float(np.median(depth)),
        "depth_max": float(depth.max()),
        "input_digest": array_digest(reads),
    }

    for arm, score in scores.items():
        report["primary"][arm] = prf(score[frame] >= thresholds[arm], snp[frame])

    report["oracle_diagnostic"] = {
        "_note": "test-optimal thresholds; a diagnostic ceiling, NOT a result",
    }
    for arm, score in scores.items():
        threshold, stats = pick_threshold(score[frame], snp[frame], grid_for(arm))
        report["oracle_diagnostic"][arm] = {"threshold": threshold, **stats}

    for label, bin_mask in depth_bins(depth):
        combined = bin_mask & frame
        entry = {"loci": int(combined.sum()), "snp": int(snp[combined].sum())}
        for arm, score in scores.items():
            entry[arm] = prf(score[combined] >= thresholds[arm], snp[combined])
        report["by_depth"][label] = entry

    # ---------- paired uncertainty against the fixed-epsilon baseline ----------
    baseline_calls = scores["binomial_v1"][frame] >= thresholds["binomial_v1"]
    report["bootstrap_vs_binomial_v1"] = {
        arm: paired_bootstrap(scores[arm][frame] >= thresholds[arm], baseline_calls,
                              snp[frame])
        for arm in scores if arm != "binomial_v1"
    }

    # ---------- error analysis ----------
    report["error_analysis"] = {
        arm: disagreement_analysis(scores, thresholds, extras, counts, snp, frame, arm)
        for arm in dict.fromkeys((PRIMARY_NEW_ARM, validation_selected, "poisson_binomial",
                                  "binomial_mean_q", "binomial_fixed_calibrated"))
    }

    report["diagnostics"]["test_epsilon"] = {
        arm: {
            "all": epsilon_summary(extras["epsilon"][arm], frame),
            "normal": epsilon_summary(extras["epsilon"][arm], frame & ~snp),
            "snp": epsilon_summary(extras["epsilon"][arm], frame & snp),
        }
        for arm in extras["epsilon"]
    }
    report["diagnostics"]["test_estimators"] = {
        arm: value for arm, value in extras["per_arm"].items() if isinstance(value, dict)
    }

    # Per-depth epsilon and quality context for the primary new arm.
    report["diagnostics"]["by_depth_context"] = {}
    for label, bin_mask in depth_bins(depth):
        combined = bin_mask & frame
        if not combined.any():
            continue
        report["diagnostics"]["by_depth_context"][label] = {
            "epsilon_ref_mean": epsilon_summary(extras["epsilon"][PRIMARY_NEW_ARM], combined),
            "mean_vaf_snp": float(scores["frequency"][combined & snp].mean())
            if (combined & snp).any() else None,
            "mean_k_snp": float(extras["k"][combined & snp].mean())
            if (combined & snp).any() else None,
            "mean_llr_v1_snp": float(scores["binomial_v1"][combined & snp].mean())
            if (combined & snp).any() else None,
        }

    report["runtime"] = {
        "cache_load_and_crosscheck_seconds": load_seconds,
        "all_arms_scoring_seconds": score_seconds,
        "quality_evidence_seconds": extras["per_arm"]["evidence_seconds"],
        "per_arm_seconds": {arm: value["seconds"] for arm, value in extras["per_arm"].items()
                            if isinstance(value, dict) and "seconds" in value},
        "loci": int(labels.size),
    }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2, default=float))
    logger.info("Wrote %s", args.out)

    print("\n=== PRIMARY: SNP vs Normal, validation-frozen thresholds ===")
    for arm in scores:
        stats = report["primary"][arm]
        print(f"  {arm:<24} thr={thresholds[arm]:>7.3f}  P={stats['precision']:.4f} "
              f"R={stats['recall']:.4f} F1={stats['f1']:.4f} "
              f"TP={stats['tp']} FP={stats['fp']} FN={stats['fn']}")

    print("\n=== Paired bootstrap vs binomial_v1 (dF1, 95% CI) ===")
    for arm, stats in report["bootstrap_vs_binomial_v1"].items():
        low, high = stats["delta_f1_ci95"]
        print(f"  {arm:<24} dF1={stats['delta_f1']:+.4f}  [{low:+.4f}, {high:+.4f}]  "
              f"dTP={stats['delta_tp']:+d} dFP={stats['delta_fp']:+d} "
              f"dFN={stats['delta_fn']:+d}")


if __name__ == "__main__":
    main()
