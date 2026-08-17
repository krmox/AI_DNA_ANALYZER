"""Regression tests for the devlog 15 cheap safety layer.

Three things must stay true no matter what else changes:

1. the safety layer never alters the frozen router's behaviour;
2. its LLR arithmetic is the *project's* binomial LLR, not a lookalike;
3. no truth label, position or PB score can reach a routing decision.

The third is enforced structurally by signature introspection, the same pattern
``test_cheap_router.py`` and ``test_cascade.py`` use, rather than by review.
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import cascade
import safety_layer
from binomial_baseline import BinomialVariantCaller
from safety_layer import (FROZEN_BINOMIAL_THRESHOLD, FROZEN_SAFETY_RULE,
                          FROZEN_SAFETY_RULE_LEAN, binomial_llr_at_epsilon,
                          binomial_llr_per_locus_epsilon, epsilon_flip,
                          guarded_route, quality_epsilon, quality_flip,
                          safety_escalate)

RNG = np.random.default_rng(20260817)


@pytest.fixture
def counts():
    """A random but reproducible (k, n) population spanning the interesting range."""
    n = RNG.integers(0, 120, size=4000).astype(np.float64)
    k = np.floor(n * RNG.random(4000) ** 3).astype(np.float64)
    return k, n


# -- the frozen constants are read, never redefined ---------------------------


def test_frozen_constants_come_from_cascade():
    assert safety_layer.FROZEN_BINOMIAL_THRESHOLD is cascade.FROZEN_BINOMIAL_THRESHOLD
    assert safety_layer.FROZEN_ROUTER_CUTOFF is cascade.FROZEN_ROUTER_CUTOFF
    assert cascade.FROZEN_BINOMIAL_THRESHOLD == 7.0
    assert cascade.FROZEN_ROUTER_CUTOFF == 5.411872376933351


def test_frozen_rule_matches_the_selection_record():
    """The production rule must be the one the §7.4 selection actually froze."""
    import json
    from pathlib import Path
    record = Path("results/bench_v15/rule_selection.json")
    if not record.exists():
        pytest.skip("rule_selection.json not present in this checkout")
    report = json.loads(record.read_text())
    for tag, rule in (("frozen_full", FROZEN_SAFETY_RULE),
                      ("frozen_lean", FROZEN_SAFETY_RULE_LEAN)):
        for key, value in rule.items():
            assert report[tag][key] == value, f"{tag}.{key} drifted from the frozen record"


# -- the arithmetic is the project's binomial LLR -----------------------------


def test_llr_matches_the_binomial_caller(counts):
    """``binomial_llr_at_epsilon`` must reproduce ``BinomialVariantCaller``."""
    k, n = counts
    for epsilon in (1e-4, 1e-3, 0.01, 0.05, 0.25):
        caller = BinomialVariantCaller(error_rate=epsilon, error_floor=min(epsilon, 1e-9))
        # reference base 'C' (index 2), so column 0 (A) is the candidate alternate
        matrix = np.stack([k, n - k, np.zeros_like(k), np.zeros_like(k)], axis=1)
        expected = caller.score_counts(matrix, np.full(k.size, 2)).llr
        np.testing.assert_allclose(binomial_llr_at_epsilon(k, n, epsilon), expected, atol=1e-9)


def test_quality_epsilon_matches_the_v2_caller(counts):
    """``quality_epsilon`` must reproduce ``BinomialVariantCaller._epsilon`` (v2)."""
    k, n = counts
    quality_sum = n * RNG.uniform(10, 45, size=n.size)
    for floor in (1e-4, 1e-3, 3e-3):
        caller = BinomialVariantCaller(quality_derived_epsilon=True, error_floor=floor)
        expected = caller._epsilon(n, quality_sum)
        np.testing.assert_allclose(quality_epsilon(n, quality_sum, floor), expected, atol=1e-12)


def test_per_locus_llr_agrees_with_scalar_llr(counts):
    k, n = counts
    scalar = binomial_llr_at_epsilon(k, n, 0.01)
    per_locus = binomial_llr_per_locus_epsilon(k, n, np.full(k.size, 0.01))
    np.testing.assert_allclose(scalar, per_locus, atol=1e-12)


def test_zero_depth_scores_zero():
    zeros = np.zeros(5)
    np.testing.assert_array_equal(binomial_llr_at_epsilon(zeros, zeros, 0.01), zeros)
    np.testing.assert_array_equal(
        binomial_llr_per_locus_epsilon(zeros, zeros, np.full(5, 0.01)), zeros)


def test_llr_is_always_finite(counts):
    k, n = counts
    for epsilon in (1e-6, 1e-4, 0.01, 0.4, 0.99):
        assert np.all(np.isfinite(binomial_llr_at_epsilon(k, n, epsilon)))


def test_llr_rejects_bad_epsilon(counts):
    k, n = counts
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            binomial_llr_at_epsilon(k, n, bad)


def test_llr_rejects_shape_mismatch():
    with pytest.raises(ValueError):
        binomial_llr_at_epsilon(np.zeros(3), np.zeros(4), 0.01)


def test_epsilon_flip_rejects_inverted_bracket(counts):
    k, n = counts
    with pytest.raises(ValueError):
        epsilon_flip(k, n, 0.05, 0.001)
    with pytest.raises(ValueError):
        epsilon_flip(k, n, 0.01, 0.01)


# -- the layer never touches the router ---------------------------------------


def test_escalation_is_disjoint_from_routing(counts):
    k, n = counts
    binomial_llr = binomial_llr_at_epsilon(k, n, 0.01)
    routed = cascade.route_stage1(binomial_llr)
    escalated = safety_escalate(k, n, routed)
    assert not np.any(escalated & routed)


def test_guarded_route_preserves_the_frozen_router(counts):
    """Every locus the router routes is still routed; none is ever un-routed."""
    k, n = counts
    binomial_llr = binomial_llr_at_epsilon(k, n, 0.01)
    routed = cascade.route_stage1(binomial_llr)
    to_pb, escalated = guarded_route(binomial_llr, k, n)
    assert np.all(to_pb[routed]), "the safety layer must not un-route a routed locus"
    np.testing.assert_array_equal(to_pb, routed | escalated)


def test_router_output_is_unchanged_by_the_layer(counts):
    """Calling the safety layer must not mutate its inputs or the router."""
    k, n = counts
    binomial_llr = binomial_llr_at_epsilon(k, n, 0.01)
    before = cascade.route_stage1(binomial_llr).copy()
    k_before, n_before = k.copy(), n.copy()
    guarded_route(binomial_llr, k, n)
    np.testing.assert_array_equal(cascade.route_stage1(binomial_llr), before)
    np.testing.assert_array_equal(k, k_before)
    np.testing.assert_array_equal(n, n_before)


def test_k_min_only_shrinks_the_escalated_set(counts):
    k, n = counts
    routed = cascade.route_stage1(binomial_llr_at_epsilon(k, n, 0.01))
    previous = None
    for k_min in (1, 2, 3, 4, 5):
        rule = {**FROZEN_SAFETY_RULE, "k_min": k_min}
        mask = safety_escalate(k, n, routed, rule)
        if previous is not None:
            assert np.all(mask <= previous), "raising k_min must not add escalations"
        previous = mask


def test_wider_bracket_only_grows_the_escalated_set(counts):
    """A superset bracket must escalate a superset of loci at fixed k_min."""
    k, n = counts
    routed = cascade.route_stage1(binomial_llr_at_epsilon(k, n, 0.01))
    narrow = safety_escalate(k, n, routed, {"epsilon_lo": 1e-3, "epsilon_hi": 1e-2, "k_min": 1})
    wide = safety_escalate(k, n, routed, {"epsilon_lo": 1e-4, "epsilon_hi": 5e-2, "k_min": 1})
    # A flip inside [1e-3, 1e-2] is a flip inside the enclosing bracket, because
    # the binomial call is monotone in epsilon at fixed (k, n).
    assert np.all(narrow <= wide)


# -- the quality (lean) rule --------------------------------------------------


def test_quality_rule_needs_quality_sum(counts):
    k, n = counts
    routed = np.zeros(k.size, dtype=bool)
    with pytest.raises(ValueError):
        safety_escalate(k, n, routed, FROZEN_SAFETY_RULE_LEAN)


def test_rule_without_parameterisation_is_rejected(counts):
    k, n = counts
    with pytest.raises(ValueError):
        safety_escalate(k, n, np.zeros(k.size, dtype=bool), {"k_min": 2})


def test_quality_flip_is_empty_when_quality_implies_the_fixed_rate(counts):
    """At mean Q20 the derived epsilon is 0.01, so nothing may flip."""
    k, n = counts
    quality_sum = n * 20.0
    assert not np.any(quality_flip(k, n, quality_sum, floor=1e-4))


# -- leakage boundary ---------------------------------------------------------

FORBIDDEN = ("label", "truth", "snp", "vcf", "position", "pb_llr", "pb", "answer")


@pytest.mark.parametrize("function", [
    safety_layer.binomial_llr_at_epsilon,
    safety_layer.binomial_llr_per_locus_epsilon,
    safety_layer.quality_epsilon,
    safety_layer.epsilon_flip,
    safety_layer.quality_flip,
    safety_layer.safety_escalate,
    safety_layer.guarded_route,
])
def test_no_routing_function_can_see_the_answer(function):
    """No parameter may name a truth label, a coordinate or a PB score.

    ``routed_to_pb`` is the one permitted 'pb' substring: it is the frozen
    router's own decision mask, not a PB score.
    """
    for name in inspect.signature(function).parameters:
        if name == "routed_to_pb":
            continue
        assert not any(token in name.lower() for token in FORBIDDEN), \
            f"{function.__name__} takes a forbidden parameter {name!r}"


def test_safety_layer_does_not_import_truth_or_pb():
    """The module must not reach a VCF, a label set or the PB caller."""
    source = inspect.getsource(safety_layer)
    for banned in ("quality_error_model", "read_level_pileup", "pysam",
                   "LABEL_SNP", "poisson_binomial"):
        assert banned not in source, f"safety_layer must not reference {banned}"


# -- the escalated set is a compute decision, and a small one -----------------


def test_escalation_is_a_small_minority(counts):
    """Sanity floor: the layer must not degenerate into 'escalate everything'."""
    k, n = counts
    routed = cascade.route_stage1(binomial_llr_at_epsilon(k, n, 0.01))
    escalated = safety_escalate(k, n, routed)
    assert escalated.mean() < 0.5


# -- sparse evaluation must be bit-identical to the dense computation ---------


def _dense_escalate(k, n, routed, rule, quality_sum=None):
    """The naive dense form the sparse implementation must reproduce exactly."""
    if "epsilon_lo" in rule:
        flip = epsilon_flip(k, n, rule["epsilon_lo"], rule["epsilon_hi"])
    else:
        flip = quality_flip(k, n, quality_sum, rule["quality_floor"])
    return flip & (np.asarray(k) >= float(rule["k_min"])) & ~routed


@pytest.mark.parametrize("rule", [
    FROZEN_SAFETY_RULE,
    {"epsilon_lo": 1e-3, "epsilon_hi": 1e-2, "k_min": 1},
    {"epsilon_lo": 1e-4, "epsilon_hi": 0.1, "k_min": 2},
    {"epsilon_lo": 3e-3, "epsilon_hi": 3e-2, "k_min": 4},
])
def test_sparse_escalation_matches_dense(rule, counts):
    k, n = counts
    routed = cascade.route_stage1(binomial_llr_at_epsilon(k, n, 0.01))
    np.testing.assert_array_equal(safety_escalate(k, n, routed, rule),
                                  _dense_escalate(k, n, routed, rule))


def test_sparse_quality_escalation_matches_dense(counts):
    k, n = counts
    quality_sum = n * RNG.uniform(8, 45, size=n.size)
    routed = cascade.route_stage1(binomial_llr_at_epsilon(k, n, 0.01))
    rule = FROZEN_SAFETY_RULE_LEAN
    np.testing.assert_array_equal(
        safety_escalate(k, n, routed, rule, quality_sum),
        _dense_escalate(k, n, routed, rule, quality_sum))


def test_zero_support_loci_can_never_escalate(counts):
    """The gate's correctness rests on this: k == 0 cannot flip at any epsilon."""
    k, n = counts
    n = np.maximum(n, 1.0)
    zero = np.zeros_like(k)
    for epsilon in (1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5):
        assert np.all(binomial_llr_at_epsilon(zero, n, epsilon) < FROZEN_BINOMIAL_THRESHOLD)


def test_all_zero_support_returns_an_empty_mask():
    n = np.full(50, 30.0)
    k = np.zeros(50)
    routed = np.zeros(50, dtype=bool)
    assert not safety_escalate(k, n, routed).any()
