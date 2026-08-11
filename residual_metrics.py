"""Shared metric helpers for the residual-binomial experiment.

Factored out so the initialization sanity check, the training loop and the
final evaluation all score with the *same* code. A metric that is computed two
slightly different ways in two scripts is a comparison bug waiting to happen,
and this experiment turns on small differences.

Ranking metrics (AUC, AP) are delegated to ``run_train`` -- the same
implementations every previous experiment in this project reported.
"""

from __future__ import annotations

import numpy as np
import torch

from run_train import average_precision, binary_auc

#: Depth strata used throughout the 15x experiments.
DEPTH_BINS: tuple[tuple[int, int], ...] = (
    (1, 4), (5, 9), (10, 14), (15, 19), (20, 29), (30, 10 ** 9),
)


def prf(predicted: np.ndarray, truth: np.ndarray) -> dict:
    """Precision/recall/F1 plus the confusion counts.

    Args:
        predicted: Boolean call per locus.
        truth: Boolean truth per locus.

    Returns:
        Dict with precision, recall, f1, tp, fp, fn.
    """
    predicted = np.asarray(predicted, dtype=bool)
    truth = np.asarray(truth, dtype=bool)
    tp = int(np.sum(predicted & truth))
    fp = int(np.sum(predicted & ~truth))
    fn = int(np.sum(~predicted & truth))
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)
    return {"precision": precision, "recall": recall, "f1": f1, "tp": tp, "fp": fp, "fn": fn}


def ranking_metrics(score: np.ndarray, truth: np.ndarray) -> dict:
    """ROC-AUC and PR-AUC (average precision) for a continuous score.

    Args:
        score: Continuous decision score per locus; higher means more SNP-like.
        truth: Boolean truth per locus.

    Returns:
        Dict with ``roc_auc`` and ``pr_auc``.
    """
    scores = torch.tensor(np.asarray(score, dtype=np.float64), dtype=torch.float64)
    positives = torch.tensor(np.asarray(truth, dtype=bool))
    return {"roc_auc": float(binary_auc(scores, positives)),
            "pr_auc": float(average_precision(scores, positives))}


def bootstrap_f1(predicted: np.ndarray, truth: np.ndarray, draws: int = 2000,
                 seed: int = 20260811) -> list[float]:
    """Percentile bootstrap 95% CI for F1, resampling loci with replacement.

    The seed is fixed so that every method in a comparison is resampled on the
    *same* bootstrap draws, which makes the intervals comparable rather than
    independently noisy.

    Args:
        predicted: Boolean calls.
        truth: Boolean truth.
        draws: Number of bootstrap resamples.
        seed: Fixed RNG seed.

    Returns:
        ``[low, high]`` 2.5th and 97.5th percentiles of the F1 distribution.
    """
    rng = np.random.default_rng(seed)
    n = truth.size
    scores = np.empty(draws)
    for i in range(draws):
        idx = rng.integers(0, n, n)
        scores[i] = prf(predicted[idx], truth[idx])["f1"]
    return [float(np.percentile(scores, 2.5)), float(np.percentile(scores, 97.5))]


def paired_bootstrap_f1_delta(predicted_a: np.ndarray, predicted_b: np.ndarray,
                              truth: np.ndarray, draws: int = 2000,
                              seed: int = 20260811) -> dict:
    """Paired bootstrap of ``F1(a) - F1(b)`` on identical resamples.

    Pairing matters: the two methods share most of their errors, so comparing
    two independently-bootstrapped intervals badly overstates the uncertainty
    in their *difference*. This is the statistic the decision criteria use.

    Args:
        predicted_a: Boolean calls from the method under test.
        predicted_b: Boolean calls from the reference method.
        truth: Boolean truth.
        draws: Number of bootstrap resamples.
        seed: Fixed RNG seed.

    Returns:
        Dict with the point delta, its 95% CI, and the fraction of draws in
        which ``a`` beat ``b``.
    """
    rng = np.random.default_rng(seed)
    n = truth.size
    deltas = np.empty(draws)
    for i in range(draws):
        idx = rng.integers(0, n, n)
        deltas[i] = prf(predicted_a[idx], truth[idx])["f1"] - prf(predicted_b[idx], truth[idx])["f1"]
    return {
        "delta_f1": prf(predicted_a, truth)["f1"] - prf(predicted_b, truth)["f1"],
        "ci95": [float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))],
        "fraction_of_draws_favouring_a": float(np.mean(deltas > 0)),
    }


def depth_bin_masks(depth: np.ndarray):
    """Yield ``(label, mask)`` for each standard depth stratum."""
    for low, high in DEPTH_BINS:
        label = f"{low}-{high}" if high < 10 ** 9 else f"{low}+"
        yield label, (depth >= low) & (depth <= high)


def vaf_from_counts(counts: np.ndarray) -> np.ndarray:
    """Alternate-allele fraction, alt chosen as the top non-reference base.

    Mirrors the frequency baseline and ``BinomialVariantCaller``'s candidate
    choice: counts only, no truth.

    Args:
        counts: ``[N, 10]`` count matrix.

    Returns:
        ``[N]`` VAF in [0, 1].
    """
    bases = counts[:, 0:4].astype(np.float64)
    total = bases.sum(axis=1)
    reference_column = counts[:, 9].astype(int) - 1
    masked = np.where(np.arange(4)[None, :] == reference_column[:, None], -1.0, bases)
    k = masked.max(axis=1)
    k = np.where(k < 0, 0.0, k)
    return np.where(total > 0, k / np.maximum(total, 1.0), 0.0)


__all__ = [
    "DEPTH_BINS", "prf", "ranking_metrics", "bootstrap_f1", "paired_bootstrap_f1_delta",
    "depth_bin_masks", "vaf_from_counts", "binary_auc", "average_precision",
]
