"""Unit tests for the three-stage cascade router.

Follows the style of ``test_cheap_router.py``/``test_hybrid_router.py``: the
leakage property (no routing function accepts a truth label) matters more than
any single numeric assertion, and is checked by signature introspection.
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import cascade
from cascade import (
    FROZEN_BINOMIAL_THRESHOLD,
    FROZEN_PB_THRESHOLD,
    FROZEN_ROUTER_CUTOFF,
    STAGE_BINOMIAL,
    STAGE_MAMBA,
    STAGE_PB,
    route_stage1,
    route_stage2,
    run_cascade,
)


# ------------------------------------------------------------------ leakage ----

def test_no_routing_function_accepts_a_truth_label():
    functions = [cascade.route_stage1, cascade.route_stage2, cascade.run_cascade]
    for function in functions:
        parameters = set(inspect.signature(function).parameters)
        assert not parameters & {"truth", "labels", "snp", "y"}, function.__name__


def test_frozen_constants_match_devlog_13_14():
    assert FROZEN_BINOMIAL_THRESHOLD == 7.0
    assert FROZEN_PB_THRESHOLD == 10.5
    assert FROZEN_ROUTER_CUTOFF == pytest.approx(5.411872376933351)


# ---------------------------------------------------------------- stage 1 ----

def test_route_stage1_matches_frozen_rule():
    binomial = np.array([7.0, 12.411, 12.413, 1.589, 1.587])
    routed = route_stage1(binomial)
    assert routed.tolist() == [True, True, False, True, False]


def test_route_stage1_boundary_equality_is_inclusive():
    binomial = np.array([FROZEN_BINOMIAL_THRESHOLD + FROZEN_ROUTER_CUTOFF])
    assert route_stage1(binomial)[0]


# ---------------------------------------------------------------- stage 2 ----

def test_route_stage2_matches_margin_rule():
    pb = np.array([10.5, 12.0, 8.6, 20.0])
    routed = route_stage2(pb, pb_to_mamba_threshold=2.0)
    assert routed.tolist() == [True, True, True, False]


def test_route_stage2_boundary_equality_is_exclusive():
    """margin < threshold, so margin == threshold does NOT route (matches the
    task spec's ``>=`` for the confident/PB side)."""
    pb = np.array([FROZEN_PB_THRESHOLD + 3.0])
    routed = route_stage2(pb, pb_to_mamba_threshold=3.0)
    assert not routed[0]


def test_route_stage2_zero_threshold_routes_nothing():
    pb = np.array([10.5, 10.5, 100.0])
    assert not route_stage2(pb, pb_to_mamba_threshold=0.0).any()


def test_route_stage2_huge_threshold_routes_everything():
    pb = np.array([10.5, -100.0, 1e6])
    assert route_stage2(pb, pb_to_mamba_threshold=1e9).all()


# -------------------------------------------------------------- empty input ----

def test_run_cascade_empty_input():
    empty = np.array([])
    calls, stage = run_cascade(empty, empty, empty, pb_to_mamba_threshold=2.0)
    assert calls.size == 0
    assert stage.size == 0


# -------------------------------------------------------- endpoint behaviour ----

def test_all_confident_stage1_never_touches_pb_or_mamba():
    """Every locus far from the binomial threshold -> 100% binomial stage."""
    rng = np.random.default_rng(0)
    n = 512
    binomial = FROZEN_BINOMIAL_THRESHOLD + FROZEN_ROUTER_CUTOFF + 1.0 + rng.uniform(size=n)
    pb = rng.normal(size=n)  # should be irrelevant
    mamba = rng.normal(size=n)  # should be irrelevant
    calls, stage = run_cascade(binomial, pb, mamba, pb_to_mamba_threshold=100.0)
    assert (stage == STAGE_BINOMIAL).all()
    assert np.array_equal(calls, binomial >= FROZEN_BINOMIAL_THRESHOLD)


def test_all_uncertain_stage1_all_confident_stage2_is_100pct_pb():
    rng = np.random.default_rng(1)
    n = 512
    binomial = np.full(n, FROZEN_BINOMIAL_THRESHOLD)  # margin 0: always routed
    pb = FROZEN_PB_THRESHOLD + 100.0 + rng.uniform(size=n)  # far margin: stays PB
    mamba = rng.normal(size=n)  # should be irrelevant
    calls, stage = run_cascade(binomial, pb, mamba, pb_to_mamba_threshold=1.0)
    assert (stage == STAGE_PB).all()
    assert np.array_equal(calls, pb >= FROZEN_PB_THRESHOLD)


def test_all_uncertain_stage2_is_100pct_mamba():
    rng = np.random.default_rng(2)
    n = 512
    binomial = np.full(n, FROZEN_BINOMIAL_THRESHOLD)  # always routed to PB
    pb = np.full(n, FROZEN_PB_THRESHOLD)  # margin 0: always routed to Mamba
    mamba = rng.normal(loc=FROZEN_PB_THRESHOLD, scale=5.0, size=n)
    calls, stage = run_cascade(binomial, pb, mamba, pb_to_mamba_threshold=1e9)
    assert (stage == STAGE_MAMBA).all()
    assert np.array_equal(calls, mamba >= FROZEN_PB_THRESHOLD)


# --------------------------------------------------------------- boundaries ----

def test_boundary_equality_at_stage1_threshold_routes_to_pb():
    binomial = np.array([FROZEN_BINOMIAL_THRESHOLD])
    pb = np.array([0.0])
    mamba = np.array([0.0])
    _, stage = run_cascade(binomial, pb, mamba, pb_to_mamba_threshold=0.0)
    assert stage[0] == STAGE_PB  # margin 0 <= cutoff -> routed; stage2 margin(0)>=0 not <0


def test_boundary_equality_at_stage2_threshold():
    binomial = np.array([FROZEN_BINOMIAL_THRESHOLD])  # always PB-routed
    pb = np.array([FROZEN_PB_THRESHOLD + 2.0])  # margin exactly 2.0
    mamba = np.array([0.0])
    _, stage = run_cascade(binomial, pb, mamba, pb_to_mamba_threshold=2.0)
    assert stage[0] == STAGE_PB  # margin == threshold: not < threshold, stays PB


# --------------------------------------------------------------- properties ----

def make_random(n: int = 4096, seed: int = 20260814):
    rng = np.random.default_rng(seed)
    binomial = rng.normal(loc=FROZEN_BINOMIAL_THRESHOLD, scale=8.0, size=n)
    pb = rng.normal(loc=FROZEN_PB_THRESHOLD, scale=10.0, size=n)
    mamba = pb + rng.normal(scale=1.0, size=n)
    return binomial, pb, mamba


@pytest.mark.parametrize("threshold", [0.0, 0.5, 2.0, 5.0, 50.0])
def test_determinism(threshold):
    binomial, pb, mamba = make_random()
    calls1, stage1 = run_cascade(binomial, pb, mamba, threshold)
    calls2, stage2 = run_cascade(binomial, pb, mamba, threshold)
    assert np.array_equal(calls1, calls2)
    assert np.array_equal(stage1, stage2)


def test_stage_counts_sum_to_n_and_no_double_processing():
    binomial, pb, mamba = make_random()
    _, stage = run_cascade(binomial, pb, mamba, pb_to_mamba_threshold=3.0)
    n = binomial.size
    is_binomial = stage == STAGE_BINOMIAL
    is_pb = stage == STAGE_PB
    is_mamba = stage == STAGE_MAMBA
    assert (is_binomial.sum() + is_pb.sum() + is_mamba.sum()) == n
    # mutually exclusive: no locus counted in more than one stage
    assert not np.any(is_binomial & is_pb)
    assert not np.any(is_binomial & is_mamba)
    assert not np.any(is_pb & is_mamba)


def test_larger_threshold_never_decreases_mamba_coverage():
    binomial, pb, mamba = make_random()
    _, stage_small = run_cascade(binomial, pb, mamba, pb_to_mamba_threshold=1.0)
    _, stage_large = run_cascade(binomial, pb, mamba, pb_to_mamba_threshold=10.0)
    assert (stage_large == STAGE_MAMBA).sum() >= (stage_small == STAGE_MAMBA).sum()


def test_mamba_subset_is_always_within_stage1_routed_set():
    binomial, pb, mamba = make_random()
    for threshold in (0.0, 1.0, 5.0, 100.0):
        _, stage = run_cascade(binomial, pb, mamba, threshold)
        routed = route_stage1(binomial)
        assert np.all(routed[stage == STAGE_MAMBA])
        assert np.all(routed[stage == STAGE_PB])
        assert not np.any(routed[stage == STAGE_BINOMIAL])


def test_test_label_independence_by_signature():
    for function in (route_stage1, route_stage2, run_cascade):
        parameters = inspect.signature(function).parameters
        assert "truth" not in parameters
        assert "label" not in parameters
        assert "labels" not in parameters
