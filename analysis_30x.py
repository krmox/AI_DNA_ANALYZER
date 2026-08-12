"""Analysis primitives that the 30x regime needs and the 15x scripts did not.

Everything here is *new questions asked of old numbers*. No scoring rule, no
feature, no prior and no threshold-selection rule is redefined: the callers
still come from ``quality_error_model`` and ``binomial_baseline``, and the
confusion arithmetic still comes from ``evaluate_binomial_baseline.prf``. This
module only adds the four diagnostics that the 15x experiment did not need:

1. :func:`depth_bins_30x` -- strata that resolve the 30x depth distribution.
   The 15x ``depth_bins`` tops out at ``30+``, which at 30x mean depth is the
   *majority* of loci and therefore carries no information.

2. :func:`error_overlap` -- which loci each arm gets wrong, not how many. The
   decisive question at 30x is whether the neural arm recovers the statistical
   arm's false negatives *without* inventing false positives, and two equal F1
   scores can hide completely disjoint error sets.

3. :func:`constant_shift_null` -- the null hypothesis that a "learned" residual
   is nothing but a moved decision boundary. The 15x experiment found a
   near-constant residual twice, so this is fitted on validation and applied to
   test as a first-class competing arm rather than as an afterthought.

4. :func:`residual_diagnostics` -- whether the residual is constant, a
   threshold shift, or genuinely locus-specific, measured as correlations
   against the evidence it would have to be reading.

Leakage boundary
----------------
:func:`constant_shift_null` is the only function here that fits anything, and
it fits on the arrays it is handed. Callers pass validation data; the test
labels reach this module only through :func:`error_overlap` and the scoring
helpers, which are pure reporting and run after the experiment is frozen.
"""

from __future__ import annotations

import numpy as np

from evaluate_binomial_baseline import prf

#: Depth strata for the 30x regime, as specified for this experiment.
DEPTH_EDGES_30X: tuple[tuple[int, int], ...] = (
    (1, 4), (5, 9), (10, 14), (15, 19), (20, 29), (30, 39), (40, 49), (50, 10 ** 9),
)


def depth_bins_30x(depth: np.ndarray) -> list[tuple[str, np.ndarray]]:
    """Depth strata resolving the 30x distribution.

    Args:
        depth: ``[N]`` per-locus read depth.

    Returns:
        ``(label, mask)`` pairs, one per stratum.
    """
    out = []
    for low, high in DEPTH_EDGES_30X:
        label = f"{low}-{high}" if high < 10 ** 9 else f"{low}+"
        out.append((label, (depth >= low) & (depth <= high)))
    return out


def error_overlap(predicted_a: np.ndarray, predicted_b: np.ndarray,
                  truth: np.ndarray, name_a: str = "a", name_b: str = "b") -> dict:
    """Set-level comparison of two arms' mistakes.

    F1 counts errors; this counts *which* errors, which is the only way to tell
    a genuine recovery from a coincidence of totals. An arm that fixes three of
    the other's false negatives while adding three false positives of its own
    is scored identically by F1 and is a different scientific object.

    Args:
        predicted_a: Boolean calls from arm A.
        predicted_b: Boolean calls from arm B.
        truth: Boolean truth.
        name_a: Label for arm A in the returned keys.
        name_b: Label for arm B in the returned keys.

    Returns:
        False-negative and false-positive set sizes, their intersection, and
        the directional recoveries in each direction.
    """
    fn_a = ~predicted_a & truth
    fn_b = ~predicted_b & truth
    fp_a = predicted_a & ~truth
    fp_b = predicted_b & ~truth
    return {
        f"{name_a}_fn": int(fn_a.sum()),
        f"{name_b}_fn": int(fn_b.sum()),
        "shared_fn": int((fn_a & fn_b).sum()),
        f"{name_a}_fn_recovered_by_{name_b}": int((fn_a & ~fn_b).sum()),
        f"{name_b}_fn_recovered_by_{name_a}": int((fn_b & ~fn_a).sum()),
        f"{name_a}_fp": int(fp_a.sum()),
        f"{name_b}_fp": int(fp_b.sum()),
        "shared_fp": int((fp_a & fp_b).sum()),
        f"{name_b}_fp_not_made_by_{name_a}": int((fp_b & ~fp_a).sum()),
        f"{name_a}_fp_not_made_by_{name_b}": int((fp_a & ~fp_b).sum()),
        "disagreements": int((predicted_a != predicted_b).sum()),
    }


def constant_shift_null(validation_score: np.ndarray, validation_truth: np.ndarray,
                        shift_grid: np.ndarray) -> tuple[float, dict]:
    """Best single additive constant on the statistical score, fitted on validation.

    A residual that is constant in the loci is algebraically a threshold move:
    ``score + c >= t`` is ``score >= t - c``. If a neural arm's advantage is
    reproduced by this one free parameter, the network has demonstrated a
    better operating point and *not* additional per-locus information. Fitting
    the shift on validation and applying it to test holds it to exactly the
    same evidential standard as the neural arm.

    Args:
        validation_score: Statistical score on the validation split.
        validation_truth: Boolean truth on the validation split.
        shift_grid: Candidate constants to add.

    Returns:
        ``(shift, validation_stats)`` for the F1-maximizing constant.
    """
    best_shift, best = float(shift_grid[0]), {"f1": -1.0}
    for candidate in shift_grid:
        stats = prf(validation_score + candidate >= 0.0, validation_truth)
        if stats["f1"] > best["f1"]:
            best_shift, best = float(candidate), stats
    return best_shift, best


def _safe_correlation(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation that returns 0.0 when either side is constant.

    Args:
        x: First variable.
        y: Second variable.

    Returns:
        The correlation, or 0.0 if undefined.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size < 2 or x.std() < 1e-12 or y.std() < 1e-12:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def roc_auc(score: np.ndarray, truth: np.ndarray) -> float:
    """ROC-AUC by the rank identity, ties averaged.

    Args:
        score: Per-locus score.
        truth: Boolean truth.

    Returns:
        The AUC, or 0.5 when one class is absent.
    """
    truth = np.asarray(truth).astype(bool)
    n_pos = int(truth.sum())
    n_neg = int(truth.size - n_pos)
    if n_pos == 0 or n_neg == 0:
        return 0.5
    order = np.argsort(np.asarray(score, dtype=np.float64), kind="stable")
    ranks = np.empty(truth.size, dtype=np.float64)
    ranks[order] = np.arange(1, truth.size + 1, dtype=np.float64)
    # Average ranks within tie groups so a constant score scores exactly 0.5.
    sorted_score = np.asarray(score, dtype=np.float64)[order]
    start = 0
    for stop in range(1, sorted_score.size + 1):
        if stop == sorted_score.size or sorted_score[stop] != sorted_score[start]:
            if stop - start > 1:
                ranks[order[start:stop]] = ranks[order[start:stop]].mean()
            start = stop
    return float((ranks[truth].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def residual_diagnostics(residual: np.ndarray, truth: np.ndarray, depth: np.ndarray,
                         vaf: np.ndarray, prior_llr: np.ndarray,
                         extra: dict[str, np.ndarray] | None = None) -> dict:
    """Characterize a learned residual as constant, shift-like or locus-specific.

    The classification rule this feeds is stated before any test number is
    read: a residual whose across-loci standard deviation is negligible beside
    its mean is a constant; one that separates the classes but correlates with
    nothing is a threshold shift; one that tracks the read-level evidence is
    locus-specific.

    Args:
        residual: ``[N]`` residual contribution to the SNP decision score.
        truth: ``[N]`` boolean SNP truth.
        depth: ``[N]`` read depth.
        vaf: ``[N]`` alternate allele fraction.
        prior_llr: ``[N]`` Poisson-binomial LLR.
        extra: Optional further variables to correlate against.

    Returns:
        Means by class, dispersion, standalone ranking power and correlations.
    """
    truth = np.asarray(truth).astype(bool)
    residual = np.asarray(residual, dtype=np.float64)
    normal = ~truth

    correlations = {
        "depth": _safe_correlation(residual, depth),
        "vaf": _safe_correlation(residual, vaf),
        "pb_llr": _safe_correlation(residual, prior_llr),
    }
    for name, values in (extra or {}).items():
        correlations[name] = _safe_correlation(residual, values)

    by_depth = {}
    for label, mask in depth_bins_30x(np.asarray(depth)):
        if mask.any():
            by_depth[label] = {"loci": int(mask.sum()),
                               "mean": float(residual[mask].mean()),
                               "std": float(residual[mask].std())}

    mean_all = float(residual.mean())
    std_all = float(residual.std())
    return {
        "mean_all": mean_all,
        "std_all": std_all,
        "mean_snp": float(residual[truth].mean()) if truth.any() else 0.0,
        "mean_normal": float(residual[normal].mean()) if normal.any() else 0.0,
        "std_snp": float(residual[truth].std()) if truth.any() else 0.0,
        "std_normal": float(residual[normal].std()) if normal.any() else 0.0,
        "percentiles": {str(p): float(np.percentile(residual, p))
                        for p in (1, 25, 50, 75, 99)},
        "roc_auc": roc_auc(residual, truth),
        "correlations": correlations,
        "by_depth": by_depth,
        # A constant shift is exactly "dispersion negligible beside location".
        "dispersion_over_location": float(std_all / max(abs(mean_all), 1e-12)),
    }


__all__ = ["DEPTH_EDGES_30X", "depth_bins_30x", "error_overlap", "constant_shift_null",
           "residual_diagnostics", "roc_auc"]
