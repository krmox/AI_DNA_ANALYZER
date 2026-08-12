"""Tests for the 30x analysis primitives.

The properties asserted here are the ones the experiment's conclusions rest
on: that the depth strata partition the loci exactly, that error overlap is a
true set operation rather than a difference of totals, that the constant-shift
null is algebraically a threshold move, and that a constant residual is
reported as carrying no ranking information.
"""

from __future__ import annotations

import numpy as np
import pytest

from analysis_30x import (
    DEPTH_EDGES_30X,
    constant_shift_null,
    depth_bins_30x,
    error_overlap,
    residual_diagnostics,
    roc_auc,
)


# --------------------------------------------------------------------------
# depth strata
# --------------------------------------------------------------------------

def test_depth_bins_partition_every_locus_exactly_once():
    depth = np.arange(0, 200)
    masks = [mask for _, mask in depth_bins_30x(depth)]
    covered = np.sum(masks, axis=0)
    # Depth 0 belongs to no stratum; every depth >= 1 belongs to exactly one.
    assert covered[0] == 0
    assert np.all(covered[1:] == 1)


def test_depth_bins_labels_match_specification():
    labels = [label for label, _ in depth_bins_30x(np.arange(60))]
    assert labels == ["1-4", "5-9", "10-14", "15-19", "20-29", "30-39", "40-49", "50+"]


def test_depth_bins_top_stratum_is_open_ended():
    _, mask = depth_bins_30x(np.array([49, 50, 5000]))[-1]
    assert mask.tolist() == [False, True, True]


def test_depth_edges_are_contiguous():
    for (_, high), (low, _) in zip(DEPTH_EDGES_30X, DEPTH_EDGES_30X[1:]):
        assert low == high + 1


# --------------------------------------------------------------------------
# error overlap
# --------------------------------------------------------------------------

def test_error_overlap_counts_are_set_operations_not_totals():
    truth = np.array([True, True, True, False, False])
    # A and B each miss exactly one SNP, but different ones.
    a = np.array([False, True, True, False, False])
    b = np.array([True, False, True, False, False])
    out = error_overlap(a, b, truth, "a", "b")
    assert out["a_fn"] == 1 and out["b_fn"] == 1
    assert out["shared_fn"] == 0
    assert out["a_fn_recovered_by_b"] == 1
    assert out["b_fn_recovered_by_a"] == 1
    assert out["disagreements"] == 2


def test_error_overlap_identical_arms_have_no_disagreement():
    rng = np.random.default_rng(0)
    truth = rng.random(500) < 0.1
    calls = rng.random(500) < 0.1
    out = error_overlap(calls, calls, truth, "a", "b")
    assert out["disagreements"] == 0
    assert out["a_fn"] == out["b_fn"] == out["shared_fn"]
    assert out["a_fp"] == out["b_fp"] == out["shared_fp"]
    assert out["a_fn_recovered_by_b"] == 0


def test_error_overlap_detects_recall_bought_with_false_positives():
    """The failure mode the experiment must not mistake for an improvement."""
    truth = np.array([True] + [False] * 20)
    pb = np.array([False] + [False] * 20)          # 0 TP, 0 FP, 1 FN
    mamba = np.array([True] + [True] * 12 + [False] * 8)  # 1 TP, 12 FP, 0 FN
    out = error_overlap(pb, mamba, truth, "pb", "mamba")
    assert out["pb_fn_recovered_by_mamba"] == 1
    assert out["mamba_fp_not_made_by_pb"] == 12


# --------------------------------------------------------------------------
# constant-shift null
# --------------------------------------------------------------------------

def test_constant_shift_is_algebraically_a_threshold_move():
    rng = np.random.default_rng(1)
    score = rng.normal(size=2000)
    truth = score + rng.normal(scale=0.5, size=2000) > 1.0
    grid = np.arange(-3.0, 3.0, 0.1)
    shift, stats = constant_shift_null(score, truth, grid)
    # Adding c and thresholding at 0 must equal thresholding the raw score at -c.
    from evaluate_binomial_baseline import prf
    assert stats == prf(score >= -shift, truth)


def test_constant_shift_recovers_a_known_offset():
    score = np.array([-5.0, -1.0, -0.5, 0.5, 1.0, 5.0])
    truth = np.array([False, False, False, True, True, True])
    shift, stats = constant_shift_null(score, truth, np.arange(-2.0, 2.01, 0.25))
    assert stats["f1"] == pytest.approx(1.0)
    # Any shift placing the boundary inside (-0.5, 0.5] separates the classes.
    assert -0.5 <= shift <= 0.5


# --------------------------------------------------------------------------
# residual diagnostics
# --------------------------------------------------------------------------

def _diag(residual, truth=None, n=1000, seed=0):
    rng = np.random.default_rng(seed)
    truth = rng.random(n) < 0.05 if truth is None else truth
    return residual_diagnostics(
        residual, truth,
        depth=rng.integers(1, 60, size=n),
        vaf=rng.random(n),
        prior_llr=rng.normal(size=n))


def test_constant_residual_carries_no_ranking_information():
    out = _diag(np.full(1000, 2.5))
    assert out["std_all"] == pytest.approx(0.0)
    assert out["roc_auc"] == pytest.approx(0.5)
    assert out["mean_snp"] == pytest.approx(out["mean_normal"])
    assert out["dispersion_over_location"] == pytest.approx(0.0)


def test_constant_residual_flagged_by_dispersion_ratio():
    """A near-constant shift has dispersion tiny beside its location."""
    rng = np.random.default_rng(3)
    residual = 3.0 + rng.normal(scale=1e-3, size=1000)
    assert _diag(residual)["dispersion_over_location"] < 0.01


def test_locus_specific_residual_has_dispersion_comparable_to_location():
    rng = np.random.default_rng(4)
    residual = rng.normal(loc=0.1, scale=2.0, size=1000)
    assert _diag(residual)["dispersion_over_location"] > 1.0


def test_residual_correlation_with_depth_is_detected():
    rng = np.random.default_rng(5)
    depth = rng.integers(1, 60, size=2000)
    residual = 0.1 * depth + rng.normal(scale=0.01, size=2000)
    out = residual_diagnostics(residual, rng.random(2000) < 0.05, depth,
                               vaf=rng.random(2000), prior_llr=rng.normal(size=2000))
    assert out["correlations"]["depth"] > 0.99
    assert abs(out["correlations"]["vaf"]) < 0.2


def test_residual_diagnostics_reports_every_populated_depth_stratum():
    rng = np.random.default_rng(6)
    depth = np.array([2, 7, 12, 17, 25, 35, 45, 55] * 10)
    out = residual_diagnostics(rng.normal(size=80), rng.random(80) < 0.2, depth,
                               vaf=rng.random(80), prior_llr=rng.normal(size=80))
    assert set(out["by_depth"]) == {"1-4", "5-9", "10-14", "15-19",
                                    "20-29", "30-39", "40-49", "50+"}


def test_extra_correlations_are_passed_through():
    rng = np.random.default_rng(7)
    residual = rng.normal(size=500)
    out = residual_diagnostics(residual, rng.random(500) < 0.1,
                               depth=rng.integers(1, 60, 500), vaf=rng.random(500),
                               prior_llr=rng.normal(size=500),
                               extra={"strand_bias": residual})
    assert out["correlations"]["strand_bias"] == pytest.approx(1.0)


# --------------------------------------------------------------------------
# roc_auc
# --------------------------------------------------------------------------

def test_roc_auc_perfect_and_inverted():
    score = np.arange(100.0)
    truth = np.arange(100) >= 50
    assert roc_auc(score, truth) == pytest.approx(1.0)
    assert roc_auc(-score, truth) == pytest.approx(0.0)


def test_roc_auc_single_class_is_defined():
    assert roc_auc(np.arange(10.0), np.zeros(10, bool)) == pytest.approx(0.5)
    assert roc_auc(np.arange(10.0), np.ones(10, bool)) == pytest.approx(0.5)


def test_roc_auc_matches_sklearn_on_random_data():
    sklearn_metrics = pytest.importorskip("sklearn.metrics")
    rng = np.random.default_rng(8)
    score = rng.normal(size=3000)
    truth = rng.random(3000) < 0.2
    assert roc_auc(score, truth) == pytest.approx(
        sklearn_metrics.roc_auc_score(truth, score), abs=1e-9)


def test_roc_auc_handles_ties_by_averaging_ranks():
    score = np.array([1.0, 1.0, 1.0, 1.0])
    truth = np.array([True, True, False, False])
    assert roc_auc(score, truth) == pytest.approx(0.5)
