"""Tensor-level controls that isolate *what kind* of information the model uses.

Each control transforms the cached ``uint8`` read tensor before features are
built, so the model, the optimizer and the evaluation are untouched and the
only changed variable is the information content of the input.

Two of the controls act on the association between a read's *base* and its
*attributes*, which is precisely the thing an aggregate caller cannot see:

* :func:`aggregate_only` (CONTROL 1) replaces every per-read attribute with its
  locus mean, so the tensor retains only what the 14-channel representation
  could already express. A read-level model that still improves under this
  transform is winning through architecture or optimization, not read-level
  information.
* :func:`decouple_base_from_attributes` (CONTROL 6) keeps the observed bases
  and the aggregate counts exactly, but reassigns each read's attribute bundle
  (quality, strand, position) to a different read within the same locus. The
  marginal distribution of every attribute is preserved and so is every count
  the binomial caller uses; only the answer to "*which* reads carry the
  alternate allele, and under what conditions" is destroyed.

:func:`permute_reads` (CONTROL 2) permutes whole rows. Given an
order-invariant aggregator this provably changes nothing, and
``test_read_level.py`` proves the invariance exactly; it is run anyway as a
numerical confirmation that the deployed pipeline has no hidden order
dependence.

Every transform is deterministic given its seed.
"""

from __future__ import annotations

import numpy as np

from read_level_pileup import FLAG_COUNTED, FLAG_NEAR_END, FLAG_VALID, read_mask

#: Columns holding per-read attributes, as opposed to the observed base and the
#: flags that determine the aggregate counts.
ATTRIBUTE_COLUMNS = (1, 2, 3, 4, 5)     # base_quality, mapping_quality, strand,
                                        # read_position, end_distance


def permute_reads(reads: np.ndarray, seed: int = 1234) -> np.ndarray:
    """CONTROL 2: permute whole read rows within each locus.

    Preserves everything -- number of reads, base counts, every marginal and
    every within-read association. Only row order changes.

    Args:
        reads: ``[N, R, C]`` uint8 tensor.
        seed: RNG seed.

    Returns:
        A new tensor with rows permuted per locus.
    """
    rng = np.random.default_rng(seed)
    out = reads.copy()
    mask = read_mask(reads)
    for locus in range(reads.shape[0]):
        count = int(mask[locus].sum())
        if count > 1:
            out[locus, :count] = reads[locus, :count][rng.permutation(count)]
    return out


def decouple_base_from_attributes(reads: np.ndarray, seed: int = 4321) -> np.ndarray:
    """CONTROL 6: keep bases and counts, randomize which read carries which.

    The attribute bundle ``(base_quality, mapping_quality, strand,
    read_position, end_distance)`` plus the near-end flag moves as a unit, so
    each bundle stays internally coherent -- a real read's quality profile is
    never fabricated, only reattached to a different observation in the same
    locus.

    Base codes and all count-determining flags stay on their original rows, and
    the permutation is **stratified by counted status**: attributes are shuffled
    among the reads the aggregate caller counts, and separately among the reads
    it discards. Stratifying is what makes the control airtight -- it preserves
    not only the allele counts and depth but ``quality_sum`` and ``mapping_sum``
    exactly, because a sum over a fixed set of rows is invariant to permuting
    values within that set. Every one of the ten aggregate count columns is
    therefore byte-identical, so *everything the binomial caller and the
    14-channel model can see is unchanged*, and the only destroyed information
    is which individual read carried which base under which conditions.

    Args:
        reads: ``[N, R, C]`` uint8 tensor.
        seed: RNG seed.

    Returns:
        A new tensor with attributes shuffled across reads within each locus.
    """
    rng = np.random.default_rng(seed)
    out = reads.copy()
    mask = read_mask(reads)
    counted = (reads[..., 6].astype(np.int64) & FLAG_COUNTED) > 0

    for locus in range(reads.shape[0]):
        valid = np.nonzero(mask[locus])[0]
        if valid.size < 2:
            continue
        for stratum in (valid[counted[locus, valid]], valid[~counted[locus, valid]]):
            if stratum.size < 2:
                continue
            order = stratum[rng.permutation(stratum.size)]
            for column in ATTRIBUTE_COLUMNS:
                out[locus, stratum, column] = reads[locus, order, column]
            # The near-end flag is a function of read position, so it travels
            # with the bundle; every other flag bit stays put.
            near_end = (reads[locus, order, 6] & FLAG_NEAR_END) > 0
            flags = out[locus, stratum, 6] & ~np.uint8(FLAG_NEAR_END)
            out[locus, stratum, 6] = flags | (near_end.astype(np.uint8) * FLAG_NEAR_END)
    return out


def aggregate_only(reads: np.ndarray) -> np.ndarray:
    """CONTROL 1: keep only information the 14-channel representation retains.

    Every valid read at a locus is given that locus' mean base quality and mean
    mapping quality, a constant strand, and a constant read position, wiping
    the near-end flag. Base codes, gaps, insertions and the counted flags are
    untouched, so allele counts, depth and the aggregate means survive while
    all per-read variation is destroyed.

    Mean base quality is taken over counted bases, matching channel 7 of the
    aggregate representation; integer rounding perturbs ``quality_sum`` by at
    most half a Phred unit per read, which is far below the resolution of any
    downstream feature.

    Args:
        reads: ``[N, R, C]`` uint8 tensor.

    Returns:
        A new tensor carrying only aggregate-derivable per-read information.
    """
    out = reads.copy()
    mask = read_mask(reads)
    counts = mask.sum(axis=1)

    quality = np.where(mask, reads[..., 1].astype(np.float64), 0.0)
    mapping = np.where(mask, reads[..., 2].astype(np.float64), 0.0)
    with np.errstate(invalid="ignore"):
        mean_quality = np.rint(quality.sum(axis=1) / np.maximum(counts, 1)).astype(np.uint8)
        mean_mapping = np.rint(mapping.sum(axis=1) / np.maximum(counts, 1)).astype(np.uint8)

    out[..., 1] = np.where(mask, mean_quality[:, None], 0)
    out[..., 2] = np.where(mask, mean_mapping[:, None], 0)
    out[..., 3] = np.where(mask, 0, 0)               # constant strand
    out[..., 4] = np.where(mask, 128, 0)             # constant read position
    out[..., 5] = np.where(mask, 255, 0)             # constant end distance
    out[..., 6] = out[..., 6] & ~np.uint8(FLAG_NEAR_END)
    return out


CONTROLS = {
    "permute_reads": permute_reads,
    "decouple": decouple_base_from_attributes,
    "aggregate_only": aggregate_only,
}


def apply_control(reads: np.ndarray, name: str, seed: int = 4321) -> np.ndarray:
    """Dispatch to a named control, or return the tensor unchanged for ''."""
    if not name:
        return reads
    if name == "aggregate_only":
        return aggregate_only(reads)
    return CONTROLS[name](reads, seed)


__all__ = ["permute_reads", "decouple_base_from_attributes", "aggregate_only",
           "apply_control", "CONTROLS", "ATTRIBUTE_COLUMNS", "FLAG_VALID"]
