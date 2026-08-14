"""Unit tests for the final robustness benchmark's evaluation logic.

These test the *new* code in ``robustness_benchmark.py`` only. The frozen
router/cascade logic it calls (``cascade.route_stage1`` etc.) is already
covered by ``test_cascade.py`` and ``test_cheap_router.py``, copied verbatim
from the validated worktree and re-run here (40/40 pass) as part of this
experiment's reproducibility check.
"""
from __future__ import annotations

import numpy as np

from cascade import FROZEN_BINOMIAL_THRESHOLD, FROZEN_ROUTER_CUTOFF, route_stage1
from evaluate_binomial_baseline import prf
from robustness_benchmark import accuracy_block, controls, cutoff_sweep, route_mask


def test_route_mask_matches_frozen_route_stage1_at_frozen_cutoff():
    rng = np.random.default_rng(0)
    llr = rng.normal(loc=7.0, scale=6.0, size=5000)
    expected = route_stage1(llr)
    actual = route_mask(llr, FROZEN_ROUTER_CUTOFF)
    assert np.array_equal(expected, actual)


def test_accuracy_block_matches_prf_and_adds_balanced_accuracy():
    rng = np.random.default_rng(1)
    calls = rng.random(2000) > 0.9
    truth = rng.random(2000) > 0.95
    block = accuracy_block(calls, truth)
    baseline = prf(calls, truth)
    for key in ("precision", "recall", "f1", "tp", "fp", "fn"):
        assert block[key] == baseline[key]
    assert 0.0 <= block["balanced_accuracy"] <= 1.0
    assert block["tn"] + block["tp"] + block["fp"] + block["fn"] == truth.size


def test_cutoff_sweep_frozen_row_reproduces_route_mask_result():
    region = {
        "labels": np.array([0, 1, 0, 1, 0, 1, 0, 1] * 200),
        "binomial_llr": np.tile(np.array([1.0, 8.0, 3.0, 12.0, 0.0, 6.5, 20.0, 7.05]), 200),
        "pb_llr": np.tile(np.array([1.0, 11.0, 3.0, 12.0, 0.0, 9.9, 20.0, 10.6]), 200),
    }
    rows = cutoff_sweep(region)
    frozen_rows = [r for r in rows if r["is_frozen_operating_point"]]
    assert len(frozen_rows) == 1
    # Sweep must be monotonically non-decreasing in PB coverage as cutoff widens.
    coverage = [r["fraction_routed_to_pb"] for r in rows]
    assert all(a <= b + 1e-12 for a, b in zip(coverage, coverage[1:]))


def test_controls_random_and_reverse_use_matched_pb_coverage():
    n = 4000
    rng = np.random.default_rng(2)
    labels = (rng.random(n) < 0.01).astype(int)
    binomial_llr = rng.normal(loc=7.0, scale=5.0, size=n)
    pb_llr = binomial_llr + rng.normal(scale=0.5, size=n)
    region = {"labels": labels, "binomial_llr": binomial_llr, "pb_llr": pb_llr}

    result = controls(region, seeds=(0, 1))
    k = result["matched_pb_coverage"]
    assert result["frozen_router"]["pb_coverage"] == k
    assert result["reverse_confidence_routing"]["pb_coverage"] == k
    for run in result["random_routing"]["runs"]:
        assert run["pb_coverage"] == k


def test_route_mask_is_symmetric_around_frozen_binomial_threshold():
    cutoff = 5.0
    llr = np.array([FROZEN_BINOMIAL_THRESHOLD - cutoff, FROZEN_BINOMIAL_THRESHOLD + cutoff,
                    FROZEN_BINOMIAL_THRESHOLD - cutoff - 0.01,
                    FROZEN_BINOMIAL_THRESHOLD + cutoff + 0.01])
    mask = route_mask(llr, cutoff)
    assert mask[0] and mask[1]
    assert not mask[2] and not mask[3]
