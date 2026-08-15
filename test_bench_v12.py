"""Unit tests for the bench_v12 benchmark infrastructure.

These test the *new* code only. The frozen router/cascade/PB implementations
have their own suites (``test_cascade.py``, ``test_cheap_router.py``) and are
not re-tested here.
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import bench_v12_stage1 as s1
import bench_v12_stage23 as s23
import extract_bench_v12 as ex
from cascade import FROZEN_PB_THRESHOLD


def _region(n: int = 640, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    counts = np.zeros((n, 10), dtype=np.float32)
    counts[:, 0:4] = rng.integers(0, 12, size=(n, 4))
    counts[:, 6] = counts[:, 0:4].sum(axis=1)
    counts[:, 7] = counts[:, 6] * rng.uniform(20, 38, size=n)
    counts[:, 8] = counts[:, 6] * 60
    counts[:, 9] = rng.integers(1, 5, size=n)
    depth = counts[:, 6].astype(np.float64)
    return {
        "name": "synthetic",
        "counts": counts,
        "binomial_llr": rng.normal(7.0, 6.0, size=n),
        "pb_llr": rng.normal(10.5, 9.0, size=n),
        "labels": rng.integers(0, 2, size=n),
        "depth": depth,
        "mean_base_quality": np.where(depth > 0, counts[:, 7] / np.maximum(depth, 1), 0.0),
        "mean_mapq": np.where(depth > 0, counts[:, 8] / np.maximum(depth, 1), 0.0),
        "region": [0, n],
    }


def test_strata_masks_are_label_free():
    """No stratum may be a function of the truth labels."""
    body = inspect.getsource(s1.strata_masks).split('"""')[-1]
    assert "labels" not in body
    region = _region()
    masks_a = s1.strata_masks(region)
    region_b = dict(region, labels=1 - region["labels"])
    masks_b = s1.strata_masks(region_b)
    for key in masks_a:
        assert np.array_equal(masks_a[key], masks_b[key])


def test_strata_masks_definitions():
    region = _region()
    masks = s1.strata_masks(region)
    assert np.array_equal(masks["hard_low_depth"], region["depth"] <= 9)
    assert np.array_equal(masks["hard_pb_uncertain"],
                          np.abs(region["pb_llr"] - FROZEN_PB_THRESHOLD) <= 13.5)
    assert masks["hard_low_quality"].dtype == bool


def test_arms_on_mask_matches_manual_computation():
    region = _region()
    mask = region["depth"] <= 9
    entry = s1.arms_on_mask(region, mask)
    labels = region["labels"]
    frame = ((labels == 0) | (labels == 1)) & mask
    pb_calls = region["pb_llr"][frame] >= FROZEN_PB_THRESHOLD
    truth = labels[frame] == 1
    assert entry["pb_only"]["tp"] == int(np.sum(pb_calls & truth))
    assert entry["loci"] == int(frame.sum())


def test_extractor_pool_rule_is_the_pre_registered_margin():
    assert ex.FEATURE_POOL_MARGIN == s23.CANDIDATE_MARGIN == 13.5


def test_extractor_frozen_constants_match_cascade():
    from cascade import FROZEN_BINOMIAL_THRESHOLD, FROZEN_ROUTER_CUTOFF
    assert ex.FROZEN_BINOMIAL_THRESHOLD == FROZEN_BINOMIAL_THRESHOLD
    assert ex.FROZEN_PB_THRESHOLD == FROZEN_PB_THRESHOLD
    assert ex.FROZEN_ROUTER_CUTOFF == FROZEN_ROUTER_CUTOFF


def test_breakdown_partitions_the_routed_subset():
    truth = np.array([1, 1, 0, 0, 1, 0], dtype=bool)
    pb = np.array([1, 0, 0, 1, 1, 0], dtype=bool)
    final = np.array([1, 1, 1, 1, 0, 0], dtype=bool)
    routed = np.array([1, 1, 1, 1, 1, 0], dtype=bool)
    out = s23.breakdown(pb, final, truth, routed)
    assert out["routed"] == 5
    assert (out["both_correct"] + out["mamba_fixes"] + out["mamba_introduces"]
            + out["both_wrong"]) == out["routed"]
    assert out["mamba_fixes"] == 1        # locus 1: pb wrong, final right
    assert out["mamba_introduces"] == 2   # loci 2 and 4
    assert out["decisions_changed"] == 3


def test_sweep_is_nested():
    """Larger t must route a superset of the loci a smaller t routes."""
    margin = np.abs(np.linspace(-20, 30, 500) - FROZEN_PB_THRESHOLD)
    previous = np.zeros(margin.size, dtype=bool)
    for t in s23.THRESHOLD_SWEEP:
        routed = margin < t
        assert np.all(previous <= routed)
        previous = routed


def test_threshold_sweep_matches_devlog16_grid():
    assert s23.THRESHOLD_SWEEP == [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 13.5]


def test_no_routing_function_sees_labels():
    for fn in (s23.stage2, s23.stage3):
        assert "truth" not in inspect.signature(fn).parameters
