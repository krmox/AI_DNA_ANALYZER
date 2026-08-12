"""Tests for the 30x initialization gate's numerical helpers.

The gate decides whether the experiment may run at all, so its two helpers are
tested for the properties the decision rests on: that the rank correlation is
computed correctly including ties, and that ``max_rank_swap_gap`` really does
distinguish "reordered values float32 cannot tell apart" from "reordered values
it could have".
"""

from __future__ import annotations

import numpy as np
import pytest

from sanity_30x_init_gate import max_rank_swap_gap, spearman


# --------------------------------------------------------------------------
# spearman
# --------------------------------------------------------------------------

def test_spearman_is_one_for_any_monotone_transform():
    x = np.linspace(0.1, 10.0, 500)
    assert spearman(x, np.log(x)) == pytest.approx(1.0)
    assert spearman(x, x ** 3) == pytest.approx(1.0)


def test_spearman_is_minus_one_when_reversed():
    x = np.linspace(0.0, 1.0, 200)
    assert spearman(x, -x) == pytest.approx(-1.0)


def test_spearman_matches_scipy_with_ties():
    scipy_stats = pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(0)
    x = rng.integers(0, 5, size=400).astype(float)   # heavy ties
    y = rng.integers(0, 5, size=400).astype(float)
    assert spearman(x, y) == pytest.approx(scipy_stats.spearmanr(x, y).statistic, abs=1e-12)


# --------------------------------------------------------------------------
# max_rank_swap_gap
# --------------------------------------------------------------------------

def test_no_swap_gap_when_orderings_agree():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert max_rank_swap_gap(x, x * 2.0) == 0.0


def test_swap_gap_reports_the_true_gap_of_a_real_reordering():
    reference = np.array([1.0, 2.0, 10.0])
    # A candidate that ranks the 10.0 locus below the 1.0 locus.
    candidate = np.array([5.0, 6.0, 0.0])
    assert max_rank_swap_gap(reference, candidate) == pytest.approx(9.0)


def test_swap_gap_is_tiny_when_only_float32_indistinguishable_values_swap():
    """The situation the gate must accept: rounding reorders near-ties."""
    rng = np.random.default_rng(1)
    reference = np.sort(rng.uniform(0.0, 100.0, size=2000))
    # Collide pairs so they differ far below one float32 ULP at this scale.
    reference[1::2] = reference[::2] + 1e-9
    candidate = reference.astype(np.float32).astype(np.float64)
    gap = max_rank_swap_gap(reference, candidate)
    assert gap <= float(np.spacing(np.float32(np.abs(reference).max())))


def test_swap_gap_catches_a_genuine_ranking_change():
    """The situation the gate must reject: a real, representable reordering."""
    reference = np.linspace(0.0, 100.0, 1000)
    candidate = reference.copy()
    candidate[10], candidate[900] = candidate[900], candidate[10]
    gap = max_rank_swap_gap(reference, candidate)
    assert gap > float(np.spacing(np.float32(100.0)))
    assert gap == pytest.approx(89.0, abs=0.5)
