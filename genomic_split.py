"""Deterministic genomic block split for the powered 15x SNP benchmark.

Why this exists
---------------
Every previous experiment split by *position in the locus array*: the last 10%
of the training region's window batches became validation. That was adequate
while train and test were separate genomic regions, but it left validation with
251 SNPs and -- at 30x -- zero errors, which made checkpoint selection an
unmeasurable operation. The binding constraint on the whole project is now
evaluation power, not model capacity.

The split
---------
The region is cut into 1 Mb blocks and each block is assigned a role by its
index modulo 4:

    block i % 4 == 0  -> train
    block i % 4 == 1  -> validation
    block i % 4 in (2, 3) -> test

For chr21:32-40 Mb that is 8 blocks giving 2 Mb train, 2 Mb validation, 4 Mb
test. Two properties are deliberate:

* **Train stays at 2 Mb**, the exact size every previous experiment trained on.
  Evaluation power is therefore the only variable that moves, and a change in
  outcome cannot be attributed to more training data.
* **Roles interleave** rather than occupying contiguous halves. Difficulty is
  not uniform along a chromosome -- repeat density, GC content and coverage all
  drift -- so contiguous blocks would confound "role" with "region". Round-robin
  assignment gives each role a sample spread across the whole span.

Leakage
-------
The blocks are disjoint in genomic coordinate, so no locus appears in two roles.
The only contact is at block boundaries, where a read up to ~250 bp long can
overlap the neighbouring block. That affects at most 250 bp of each 1,000,000 bp
block (0.025%), no locus is shared, and every locus is labelled from the truth
VCF rather than from its neighbours. Windows are built inside a single block
because window construction never spans a block edge by more than the 64 bp
window width, which is likewise negligible and is reported rather than hidden.
"""

from __future__ import annotations

import numpy as np

#: Block width in base pairs.
BLOCK_BP = 1_000_000

#: Role for each block index modulo the cycle length.
ROLE_CYCLE: tuple[str, ...] = ("train", "validation", "test", "test")


def block_index(positions: np.ndarray, region_start: int,
                block_bp: int = BLOCK_BP) -> np.ndarray:
    """Block index of each locus, counted from the region start.

    Args:
        positions: ``[N]`` 0-based genomic positions.
        region_start: First position of the region.
        block_bp: Block width.

    Returns:
        ``[N]`` integer block indices.
    """
    return (np.asarray(positions) - region_start) // block_bp


def assign_roles(positions: np.ndarray, region_start: int,
                 block_bp: int = BLOCK_BP,
                 cycle: tuple[str, ...] = ROLE_CYCLE) -> dict[str, np.ndarray]:
    """Boolean masks selecting the loci of each role.

    Args:
        positions: ``[N]`` 0-based genomic positions.
        region_start: First position of the region.
        block_bp: Block width.
        cycle: Role for each block index modulo ``len(cycle)``.

    Returns:
        Dict mapping role name to a ``[N]`` boolean mask. The masks partition
        the loci exactly.
    """
    index = block_index(positions, region_start, block_bp)
    role_of_block = index % len(cycle)
    masks = {}
    for slot, role in enumerate(cycle):
        mask = role_of_block == slot
        masks[role] = mask if role not in masks else (masks[role] | mask)
    return masks


def assign_window_roles(positions: np.ndarray, region_start: int, seq_len: int = 64,
                        block_bp: int = BLOCK_BP,
                        cycle: tuple[str, ...] = ROLE_CYCLE) -> dict[str, np.ndarray]:
    """Role masks assigned per *window* rather than per locus.

    The model consumes windows of ``seq_len`` consecutive loci, so a locus-level
    mask could split a window across two roles and would break the reshape.
    Assigning by the window's first position keeps windows intact, keeps the
    locus count a multiple of ``seq_len``, and guarantees that no training
    window contains a validation or test locus.

    Args:
        positions: ``[N]`` genomic positions, N a multiple of ``seq_len``.
        region_start: First position of the region.
        seq_len: Window width in loci.
        block_bp: Block width in base pairs.
        cycle: Role assignment cycle.

    Returns:
        Dict mapping role to a ``[N]`` boolean locus mask, constant within each
        window.

    Raises:
        ValueError: If ``positions`` is not a whole number of windows.
    """
    positions = np.asarray(positions)
    if positions.size % seq_len:
        raise ValueError(f"{positions.size} loci is not a multiple of {seq_len}")
    window_start = positions.reshape(-1, seq_len)[:, 0]
    role_of_window = block_index(window_start, region_start, block_bp) % len(cycle)

    masks: dict[str, np.ndarray] = {}
    for slot, role in enumerate(cycle):
        window_mask = role_of_window == slot
        locus_mask = np.repeat(window_mask, seq_len)
        masks[role] = locus_mask if role not in masks else (masks[role] | locus_mask)
    return masks


def split_summary(positions: np.ndarray, labels: np.ndarray, region_start: int,
                  snp_label: int = 1, block_bp: int = BLOCK_BP,
                  cycle: tuple[str, ...] = ROLE_CYCLE,
                  seq_len: int | None = None) -> dict:
    """Loci and variant counts per role, for the power calculation.

    Args:
        positions: ``[N]`` genomic positions.
        labels: ``[N]`` class labels.
        region_start: First position of the region.
        snp_label: Label id counted as a SNP.
        block_bp: Block width.
        cycle: Role assignment cycle.

    Returns:
        Per-role loci, SNP, insertion and deletion counts plus the blocks used.
    """
    masks = (assign_roles(positions, region_start, block_bp, cycle) if seq_len is None
             else assign_window_roles(positions, region_start, seq_len, block_bp, cycle))
    index = block_index(positions, region_start, block_bp)
    out = {}
    for role, mask in masks.items():
        out[role] = {
            "loci": int(mask.sum()),
            "megabases": float(mask.sum() / 1e6),
            "snp": int((labels[mask] == snp_label).sum()),
            "insertion": int((labels[mask] == 2).sum()),
            "deletion": int((labels[mask] == 3).sum()),
            "blocks": sorted(set(index[mask].tolist())),
        }
    return out


__all__ = ["BLOCK_BP", "ROLE_CYCLE", "block_index", "assign_roles",
           "assign_window_roles", "split_summary"]
