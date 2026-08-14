"""Unit tests for the cheap pre-router.

The property that matters most is negative: the router's inference path must
not touch the Poisson-binomial LLR, because avoiding that computation is the
entire point. A pre-router that reads PB is not a pre-router, and that mistake
would be invisible in the accuracy numbers -- it would simply look like a very
good router.
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import cheap_router
from cheap_router import (
    CHEAP_FEATURE_NAMES,
    LogisticRouter,
    alt_count_uncertainty,
    binomial_llr,
    binomial_margin_uncertainty,
    cheap_features,
    depth_uncertainty,
    fit_logistic_router,
    oracle_target,
)
from hybrid_router import cutoff_for_coverage, route


def make_counts(n_loci: int, seed: int = 20260813) -> np.ndarray:
    """A synthetic count matrix with realistic 15x structure."""
    rng = np.random.default_rng(seed)
    counts = np.zeros((n_loci, 10), dtype=np.float64)
    depth = rng.integers(1, 30, size=n_loci)
    reference_index = rng.integers(1, 5, size=n_loci)
    alt_reads = rng.binomial(depth, rng.choice([0.0, 0.02, 0.5], size=n_loci))
    for locus in range(n_loci):
        reference_slot = reference_index[locus] - 1
        alt_slot = (reference_slot + 1 + rng.integers(0, 3)) % 4
        counts[locus, reference_slot] = depth[locus] - alt_reads[locus]
        counts[locus, alt_slot] = alt_reads[locus]
    counts[:, 6] = depth
    counts[:, 7] = depth * rng.uniform(25, 38, size=n_loci)
    counts[:, 8] = depth * rng.uniform(40, 60, size=n_loci)
    counts[:, 9] = reference_index
    return counts


@pytest.fixture
def counts():
    return make_counts(4096)


# ------------------------------------------------------ the pre-PB property ----

def test_no_router_entry_point_accepts_a_pb_score():
    """The whole thesis: routing must be computable without the PB LLR."""
    entry_points = [
        cheap_router.cheap_features, cheap_router.binomial_margin_uncertainty,
        cheap_router.alt_count_uncertainty, cheap_router.depth_uncertainty,
        LogisticRouter.score, LogisticRouter.design_matrix,
    ]
    for function in entry_points:
        parameters = set(inspect.signature(function).parameters)
        assert not parameters & {"pb_llr", "pb_score", "pb_threshold"}, function.__name__
        assert not parameters & {"truth", "labels", "snp"}, function.__name__


def test_scoring_never_calls_the_poisson_binomial_prior(counts, monkeypatch):
    """Fail loudly if any inference path reaches for the expensive convolution."""
    def forbidden(*args, **kwargs):
        raise AssertionError("the pre-router computed the Poisson-binomial prior")

    monkeypatch.setattr(cheap_router, "candidate_alt", cheap_router.candidate_alt)
    import quality_error_model
    monkeypatch.setattr(quality_error_model, "poisson_binomial_llr", forbidden)

    cheap_features(counts)
    binomial_margin_uncertainty(counts, 7.0)
    alt_count_uncertainty(counts)
    router, _ = fit_logistic_router(counts, make_target(counts), use_binomial=False)
    router.score(counts)


def make_target(counts: np.ndarray) -> np.ndarray:
    """A stand-in distillation target for tests that need a fitted router."""
    rng = np.random.default_rng(5)
    signal = counts[:, 6] + rng.normal(scale=3.0, size=counts.shape[0])
    return signal <= np.quantile(signal, 0.05)


# --------------------------------------------------------------- features ----

def test_feature_matrix_shape_and_finiteness(counts):
    features = cheap_features(counts)
    assert features.shape == (counts.shape[0], len(CHEAP_FEATURE_NAMES))
    assert np.all(np.isfinite(features))
    assert features.dtype == np.float32


def test_features_are_deterministic(counts):
    assert np.array_equal(cheap_features(counts), cheap_features(counts))


def test_zero_depth_locus_is_handled():
    counts = np.zeros((1, 10))
    counts[0, 9] = 1
    features = cheap_features(counts)
    assert np.all(np.isfinite(features))
    assert features[0, list(CHEAP_FEATURE_NAMES).index("has_alt")] == 0.0


def test_reference_count_matches_the_reference_column():
    counts = np.zeros((1, 10))
    counts[0, 0] = 8      # A
    counts[0, 1] = 2      # C
    counts[0, 6] = 10
    counts[0, 9] = 1      # reference is A
    features = cheap_features(counts)
    names = list(CHEAP_FEATURE_NAMES)
    assert features[0, names.index("vaf")] == pytest.approx(0.2)
    assert features[0, names.index("ref_fraction")] == pytest.approx(0.8)
    assert features[0, names.index("k")] == pytest.approx(0.2)


def test_pooled_quality_columns_are_read_from_the_count_matrix():
    counts = np.zeros((1, 10))
    counts[0, 0] = 10
    counts[0, 6] = 10
    counts[0, 7] = 300.0   # mean Phred 30
    counts[0, 8] = 600.0   # mean MAPQ 60
    counts[0, 9] = 1
    features = cheap_features(counts)
    names = list(CHEAP_FEATURE_NAMES)
    assert features[0, names.index("mean_base_quality")] == pytest.approx(30 / 40)
    assert features[0, names.index("mean_mapping_quality")] == pytest.approx(1.0)


# ------------------------------------------------------------- analytical ----

def test_binomial_margin_peaks_at_the_threshold(counts):
    llr = binomial_llr(counts)
    uncertainty = binomial_margin_uncertainty(counts, 7.0, llr)
    assert uncertainty.max() <= 0.0
    closest = np.abs(llr - 7.0).argmin()
    assert uncertainty.argmax() == closest


def test_alt_count_rule_ranks_ambiguous_loci_above_clean_ones():
    counts = np.zeros((3, 10))
    counts[:, 9] = 1
    counts[:, 6] = 10
    counts[0, 0], counts[0, 1] = 10, 0     # no alternate: decided
    counts[1, 0], counts[1, 1] = 5, 5      # 50/50: heterozygous, clean
    counts[2, 0], counts[2, 1] = 9, 1      # one alternate read: ambiguous
    uncertainty = alt_count_uncertainty(counts)
    assert uncertainty[0] == -1.0
    assert uncertainty[2] < uncertainty[1]


def test_depth_rule_is_monotone(counts):
    uncertainty = depth_uncertainty(counts)
    order = np.argsort(uncertainty, kind="stable")
    # Depth must be non-increasing along the routing order (ties may permute).
    assert np.all(np.diff(counts[order, 6]) <= 0)


# ---------------------------------------------------------------- learned ----

def test_router_is_small(counts):
    router, diagnostics = fit_logistic_router(counts, make_target(counts),
                                              use_binomial=False)
    assert router.weights.size == len(CHEAP_FEATURE_NAMES)
    assert diagnostics["parameter_count"] < 100


def test_router_with_binomial_has_two_more_features(counts):
    router, _ = fit_logistic_router(counts, make_target(counts), use_binomial=True,
                                    binomial_threshold=7.0)
    assert router.weights.size == len(CHEAP_FEATURE_NAMES) + 2
    assert router.feature_names[-2:] == ("binomial_llr", "binomial_margin")


def test_router_learns_its_target_better_than_chance(counts):
    target = make_target(counts)
    router, _ = fit_logistic_router(counts, target, use_binomial=False)
    scores = router.score(counts)
    routed = route(scores, cutoff_for_coverage(scores, 0.10))
    assert routed[target].mean() > 3 * routed.mean()


def test_router_scoring_is_deterministic_and_label_free(counts):
    target = make_target(counts)
    router, _ = fit_logistic_router(counts, target, use_binomial=False)
    first = router.score(counts)
    assert np.array_equal(first, router.score(counts))
    # Scoring cannot depend on the target: permuting it changes nothing, since
    # the target is not an argument of the inference path at all.
    assert "target" not in inspect.signature(LogisticRouter.score).parameters


def test_standardization_uses_stored_training_statistics(counts):
    router, _ = fit_logistic_router(counts[:2048], make_target(counts)[:2048],
                                    use_binomial=False)
    held_out = router.design_matrix(counts[2048:])
    assert held_out.shape[1] == len(CHEAP_FEATURE_NAMES)
    # Statistics come from the fit, so the held-out block is not re-centred.
    assert abs(held_out.mean()) > 1e-9


# ----------------------------------------------------------------- target ----

def test_oracle_target_selects_the_requested_fraction():
    rng = np.random.default_rng(11)
    pb_llr = rng.normal(10.5, 12.0, size=10000)
    target = oracle_target(pb_llr, 10.5, fraction=0.05)
    assert target.mean() == pytest.approx(0.05, abs=0.005)


def test_oracle_target_selects_loci_near_the_threshold():
    pb_llr = np.array([10.5, 11.0, 40.0, -30.0])
    target = oracle_target(pb_llr, 10.5, fraction=0.5)
    assert target.tolist() == [True, True, False, False]
