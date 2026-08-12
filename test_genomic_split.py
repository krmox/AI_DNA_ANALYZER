"""Tests for the genomic block split.

The experiment's validity rests on the split being a true partition with no
locus in two roles, so that is what is asserted here, along with the
interleaving property that keeps regional difficulty from confounding role.
"""

from __future__ import annotations

import numpy as np
import pytest

from genomic_split import (
    ROLE_CYCLE,
    assign_roles,
    assign_window_roles,
    block_index,
    split_summary,
)


def _positions(start: int = 32_000_000, span: int = 8_000_000, step: int = 1000):
    return np.arange(start, start + span, step)


def test_roles_partition_every_locus_exactly_once():
    positions = _positions()
    masks = assign_roles(positions, 32_000_000)
    covered = np.sum(list(masks.values()), axis=0)
    assert np.all(covered == 1)


def test_no_locus_appears_in_two_roles():
    positions = _positions()
    masks = assign_roles(positions, 32_000_000)
    roles = list(masks)
    for i, a in enumerate(roles):
        for b in roles[i + 1:]:
            assert not np.any(masks[a] & masks[b])


def test_expected_block_allocation_for_eight_blocks():
    positions = _positions()
    index = block_index(positions, 32_000_000)
    masks = assign_roles(positions, 32_000_000)
    assert sorted(set(index[masks["train"]].tolist())) == [0, 4]
    assert sorted(set(index[masks["validation"]].tolist())) == [1, 5]
    assert sorted(set(index[masks["test"]].tolist())) == [2, 3, 6, 7]


def test_train_is_two_megabases_matching_previous_experiments():
    positions = _positions(step=1)
    masks = assign_roles(positions, 32_000_000)
    assert masks["train"].sum() == 2_000_000
    assert masks["validation"].sum() == 2_000_000
    assert masks["test"].sum() == 4_000_000


def test_roles_interleave_rather_than_occupying_contiguous_halves():
    """Each role must draw from both the first and second half of the span."""
    positions = _positions(step=1)
    midpoint = 32_000_000 + 4_000_000
    masks = assign_roles(positions, 32_000_000)
    for role, mask in masks.items():
        chosen = positions[mask]
        assert np.any(chosen < midpoint), role
        assert np.any(chosen >= midpoint), role


def test_block_index_is_zero_based_from_region_start():
    positions = np.array([32_000_000, 32_999_999, 33_000_000, 39_999_999])
    assert block_index(positions, 32_000_000).tolist() == [0, 0, 1, 7]


def test_split_summary_counts_variants_per_role():
    positions = _positions(step=1000)
    labels = np.zeros(positions.size, dtype=int)
    labels[block_index(positions, 32_000_000) == 2] = 1        # SNPs, all in test
    summary = split_summary(positions, labels, 32_000_000)
    assert summary["test"]["snp"] == int((labels == 1).sum())
    assert summary["train"]["snp"] == 0
    assert summary["validation"]["snp"] == 0


def test_summary_megabases_match_loci():
    positions = _positions(step=1)
    summary = split_summary(positions, np.zeros(positions.size, int), 32_000_000)
    assert summary["train"]["megabases"] == pytest.approx(2.0)
    assert summary["test"]["megabases"] == pytest.approx(4.0)


def test_cycle_is_the_documented_one():
    assert ROLE_CYCLE == ("train", "validation", "test", "test")


# --------------------------------------------------------------------------
# window-level assignment
# --------------------------------------------------------------------------

def test_window_roles_are_constant_within_each_window():
    positions = np.arange(32_000_000, 32_000_000 + 64 * 500)
    masks = assign_window_roles(positions, 32_000_000, seq_len=64)
    for mask in masks.values():
        assert np.all(mask.reshape(-1, 64).std(axis=1) == 0)


def test_window_roles_partition_and_stay_window_aligned():
    positions = np.arange(32_000_000, 32_000_000 + 64 * 500)
    masks = assign_window_roles(positions, 32_000_000, seq_len=64)
    assert np.all(np.sum(list(masks.values()), axis=0) == 1)
    for mask in masks.values():
        assert mask.sum() % 64 == 0


def test_window_assignment_rejects_a_partial_window():
    with pytest.raises(ValueError):
        assign_window_roles(np.arange(32_000_000, 32_000_100), 32_000_000, seq_len=64)


def test_no_training_window_contains_a_test_locus():
    """The property that makes the split leak-free at window granularity."""
    positions = np.arange(32_000_000, 32_000_000 + 64 * 2000)
    masks = assign_window_roles(positions, 32_000_000, seq_len=64)
    train_windows = masks["train"].reshape(-1, 64).any(axis=1)
    test_windows = masks["test"].reshape(-1, 64).any(axis=1)
    assert not np.any(train_windows & test_windows)
