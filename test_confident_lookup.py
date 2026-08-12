"""Equivalence tests for the high-confidence window lookup.

``GiabAlignmentProvider._is_confident`` was changed from a linear scan over
every BED interval to a binary search. The predicate must be *identical* --
this is core code that decides which loci enter every dataset in the project,
so a behaviour change would silently alter every experiment. These tests pin
the equivalence against a brute-force reference implementation, including the
boundary cases where an off-by-one would hide.
"""

from __future__ import annotations

import numpy as np
import pytest

from providers import GiabAlignmentProvider, load_bed_regions


class _Lookup(GiabAlignmentProvider):
    """Bare provider exposing only the confidence predicate."""

    def __init__(self, intervals, contig: str = "chr21") -> None:  # noqa: D107
        self.high_confidence_bed = "unused"
        self.contig = contig
        self._confident_regions = intervals
        self._confident_starts = [s for s, _ in intervals]
        self._confident_ends = [e for _, e in intervals]


def _brute_force(intervals, start: int, end: int) -> bool:
    """The original predicate, kept as the reference."""
    return any(a <= start and end <= b for a, b in intervals)


@pytest.fixture
def intervals():
    return [(100, 200), (300, 400), (1000, 5000), (10_000, 10_001)]


def test_matches_brute_force_across_a_dense_sweep(intervals):
    lookup = _Lookup(intervals)
    for start in range(0, 11_000, 7):
        for width in (1, 64, 200):
            assert lookup._is_confident(start, start + width) == \
                _brute_force(intervals, start, start + width)


def test_window_exactly_filling_an_interval_is_confident(intervals):
    assert _Lookup(intervals)._is_confident(100, 200)


def test_window_one_past_an_interval_end_is_not(intervals):
    assert not _Lookup(intervals)._is_confident(100, 201)


def test_window_starting_one_before_an_interval_is_not(intervals):
    assert not _Lookup(intervals)._is_confident(99, 200)


def test_window_spanning_a_gap_between_intervals_is_not(intervals):
    assert not _Lookup(intervals)._is_confident(150, 350)


def test_window_before_every_interval_is_not(intervals):
    assert not _Lookup(intervals)._is_confident(0, 64)


def test_window_after_every_interval_is_not(intervals):
    assert not _Lookup(intervals)._is_confident(50_000, 50_064)


def test_empty_interval_list_is_never_confident():
    assert not _Lookup([])._is_confident(0, 64)


def test_matches_brute_force_on_randomised_interval_sets():
    rng = np.random.default_rng(0)
    for _ in range(20):
        starts = np.sort(rng.integers(0, 50_000, size=40))
        ivs = []
        for s in starts:
            e = int(s) + int(rng.integers(1, 800))
            if ivs and s <= ivs[-1][1]:          # keep them merged/disjoint
                ivs[-1] = (ivs[-1][0], max(ivs[-1][1], e))
            else:
                ivs.append((int(s), e))
        lookup = _Lookup(ivs)
        for _ in range(300):
            start = int(rng.integers(0, 51_000))
            end = start + int(rng.integers(1, 300))
            assert lookup._is_confident(start, end) == _brute_force(ivs, start, end)


def test_matches_brute_force_on_the_real_giab_bed():
    """The case that motivated the change: thousands of real intervals."""
    path = "data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed"
    try:
        real = load_bed_regions(path, "chr21")
    except FileNotFoundError:  # pragma: no cover - dataset not present
        pytest.skip("12 Mb BED not available")
    lookup = _Lookup(real)
    rng = np.random.default_rng(1)
    for start in rng.integers(32_000_000, 44_000_000, size=4000):
        assert lookup._is_confident(int(start), int(start) + 64) == \
            _brute_force(real, int(start), int(start) + 64)
