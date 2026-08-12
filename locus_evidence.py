"""Per-locus evidence features that preserve what the 14-channel view destroys.

Why this representation exists
------------------------------
The quality-error experiment (``History/8_DEVLOG.md``) established two facts
that condemn the 14-channel ``raw_pileup`` representation as a *representation*
rather than as a model:

1. **The decisive variable is the base quality of the ALT-supporting reads.**
   The reference-only estimator, which cannot see ALT-read quality, was null
   (ΔF1 +0.0033, CI spanning zero) and added 9 false positives at loci whose
   ALT reads were poor. Every arm that *did* see ALT-read quality won
   (ΔF1 +0.015 to +0.018). ``raw_pileup`` supplies exactly one quality channel,
   ``mean_base_quality``, pooled over reference *and* alternate reads together.
   At 15x with k=2 of n=10 the ALT reads contribute a fifth of that mean, so
   the load-bearing signal arrives diluted 5:1 and confounded with the
   background. The 14-channel model is on the wrong side of the split the
   experiment proved matters.

2. **Integer ``n`` and ``k`` are not recoverable from the 14 channels.** The
   fractions use ``called + gap_reads`` as denominator while ``log_depth``
   counts reads admitted before the base-quality filter. Measured on the test
   region: the two disagree at 427,737 of 455,488 loci (93.9%). The model
   cannot reconstruct the binomial statistic, let alone the Poisson-binomial
   one, from its own inputs.

That is a sufficient explanation for the residual experiment's collapse to a
constant shift: the residual head had no locus-specific information that was
not already inside the prior, so a constant was the correct thing to learn.

What this module adds, and the hypothesis behind each group
-----------------------------------------------------------
Everything is derived from the **already cached** read-level tensor, so no new
extraction pass runs and every arm consumes byte-identical input.

* ``count`` -- integer semantics restored: ``n``, ``k``, reference count, the
  runner-up base count, gap and insertion counts, depth. *Hypothesis: the model
  should not have to rediscover sampling uncertainty from fractions.*
* ``alt_quality`` / ``ref_quality`` -- the quality **distribution** on each side
  of the ref/alt split: mean error probability, min and max Phred, the fraction
  of reads clearing Q30, and the alt-minus-ref contrast. *Hypothesis: this is
  the variable the statistical experiment proved decisive, and a single pooled
  mean cannot express it.*
* ``strand_position`` -- strand balance, read-position and end-distance
  summaries, computed **separately for ALT and reference reads**, plus mapping
  quality. *Hypothesis: this is genuine information the Poisson-binomial model
  cannot represent at all. PB treats reads as independent Bernoulli draws
  indexed only by quality; it has no way to say "these 3 high-quality ALT reads
  are all on one strand and all within 5 bp of a read end", which is the
  signature of a systematic artefact rather than a variant.*
* ``reference`` -- the reference base one-hot, unchanged from ``raw_pileup``.

A deliberate, documented change of philosophy
---------------------------------------------
``raw_pileup`` refused to reorder counts into reference/alternate slots, on the
grounds that "that reordering is itself the comparison the model must perform".
This module *does* reorder. That is a real change and it makes the model's task
easier by construction, so it is stated plainly rather than buried: the earlier
stance was appropriate while the open question was whether the network could
discover the comparison. That question is settled -- it could not, from
fractions that do not contain the comparison's operands. The remaining question
is whether anything exists *above* the statistical caller, and for that the
model must be handed the caller's operands, with the caller itself supplied as
a prior so that merely re-deriving it earns nothing.

The candidate ALT is selected by delegating to the established caller, so it is
bit-identical to every previous experiment.

Leakage boundary
----------------
No VCF, label, or model output is read here. The reference base enters only as
it does in the binomial caller: to mask the ALT candidate and to split reads
into reference-supporting and not. Labels travel in a separate array.
"""

from __future__ import annotations

import numpy as np

from quality_error_model import (
    EPSILON_CEILING,
    EPSILON_FLOOR,
    candidate_alt,
    extract_quality_evidence,
)
from read_level_pileup import FLAG_NEAR_END

#: Feature layout. Order is fixed and asserted by the sanity checks.
FEATURE_NAMES: tuple[str, ...] = (
    # -- count group: integer semantics, which the fractions destroyed --------
    "log_n", "log_k", "log_ref_count", "log_other_count",
    "log_depth", "log_gap", "log_insertion",
    "k_small", "vaf", "ref_fraction", "other_fraction", "gap_fraction",
    # -- alt-read quality distribution ---------------------------------------
    "alt_mean_phred", "alt_min_phred", "alt_max_phred",
    "alt_mean_error", "alt_fraction_q30", "alt_has_reads",
    # -- reference-read quality distribution ---------------------------------
    "ref_mean_phred", "ref_min_phred", "ref_max_phred",
    "ref_mean_error", "ref_fraction_q30",
    # -- the contrast the pooled mean cannot express -------------------------
    "alt_minus_ref_phred", "alt_minus_ref_error",
    # -- strand / position: information Poisson-binomial cannot represent ----
    "alt_forward_fraction", "alt_strand_bias", "alt_mean_read_position",
    "alt_mean_end_distance", "alt_near_end_fraction", "alt_mean_mapq",
    "ref_forward_fraction", "ref_strand_bias", "ref_mean_read_position",
    "ref_mean_end_distance", "ref_near_end_fraction", "ref_mean_mapq",
    "alt_minus_ref_mapq", "alt_minus_ref_near_end",
    # -- reference identity, unchanged from raw_pileup ------------------------
    "ref_is_N", "ref_is_A", "ref_is_C", "ref_is_G", "ref_is_T",
)
FEATURE_DIM: int = len(FEATURE_NAMES)

#: Semantic groups, for the grouped-encoder variant and for ablations.
FEATURE_GROUPS: dict[str, tuple[str, ...]] = {
    "count": FEATURE_NAMES[0:12],
    "alt_quality": FEATURE_NAMES[12:18],
    "ref_quality": FEATURE_NAMES[18:25],
    "strand_position": FEATURE_NAMES[25:39],
    "reference": FEATURE_NAMES[39:44],
}

#: Normalisation constants. Chosen so typical values land near 1 without
#: clipping, matching the spirit of ``raw_pileup``'s scales.
DEPTH_SCALE = 4.0
PHRED_SCALE = 40.0
MAPPING_QUALITY_SCALE = 60.0
K_SCALE = 10.0


def _masked_mean(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Row-wise mean over ``mask``; 0 where the row selects nothing."""
    count = mask.sum(axis=1)
    total = np.sum(np.where(mask, values, 0.0), axis=1)
    return np.where(count > 0, total / np.maximum(count, 1), 0.0)


def _masked_fraction(condition: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Row-wise fraction of masked entries satisfying ``condition``."""
    count = mask.sum(axis=1)
    return np.where(count > 0, np.sum(condition & mask, axis=1) / np.maximum(count, 1), 0.0)


def _masked_extreme(values: np.ndarray, mask: np.ndarray, largest: bool) -> np.ndarray:
    """Row-wise max (or min) over ``mask``; 0 where the row selects nothing."""
    filler = -np.inf if largest else np.inf
    filled = np.where(mask, values, filler)
    extreme = filled.max(axis=1) if largest else filled.min(axis=1)
    return np.where(mask.any(axis=1), extreme, 0.0)


def build_locus_features(reads: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """Build the ``[N, FEATURE_DIM]`` per-locus evidence matrix.

    Args:
        reads: ``[N, R, READ_DIM]`` uint8 read tensor from
            :mod:`read_level_pileup`.
        counts: ``[N, 10]`` count matrix following
            ``pileup_counts.COUNT_COLUMNS``, used for the reference index and
            for the depth/gap/insertion columns.

    Returns:
        ``[N, FEATURE_DIM]`` float32 features, in :data:`FEATURE_NAMES` order.
    """
    reads = np.asarray(reads)
    reference_index = counts[:, 9].astype(int)
    n_loci = reads.shape[0]

    evidence = extract_quality_evidence(reads, reference_index)
    alt_index, k, n = candidate_alt(counts[:, 0:4], reference_index)

    flags = reads[..., 6].astype(np.int64)
    base_code = reads[..., 0].astype(np.int64)
    phred = reads[..., 1].astype(np.float64)
    mapq = reads[..., 2].astype(np.float64)
    strand = reads[..., 3].astype(np.float64)
    read_position = reads[..., 4].astype(np.float64) / 255.0
    end_distance = reads[..., 5].astype(np.float64) / 255.0

    counted = evidence.counted
    supports_alt = counted & (base_code == alt_index[:, None])
    supports_reference = evidence.supports_reference
    # "Other": counted, but neither the reference nor the candidate ALT. A
    # locus with several competing non-reference bases looks like noise, and
    # neither the fractions nor the binomial statistic expose that directly.
    supports_other = counted & ~supports_alt & ~supports_reference

    depth = counts[:, 6]
    gap = counts[:, 4]
    insertion = counts[:, 5]
    reference_count = supports_reference.sum(axis=1).astype(np.float64)
    other_count = supports_other.sum(axis=1).astype(np.float64)
    safe_n = np.maximum(n, 1)

    error_probability = np.clip(evidence.error_probability, EPSILON_FLOOR, EPSILON_CEILING)
    is_q30 = phred >= 30.0

    def side(mask: np.ndarray) -> dict[str, np.ndarray]:
        """Distribution summaries for one side of the ref/alt split."""
        return {
            "mean_phred": _masked_mean(phred, mask) / PHRED_SCALE,
            "min_phred": _masked_extreme(phred, mask, largest=False) / PHRED_SCALE,
            "max_phred": _masked_extreme(phred, mask, largest=True) / PHRED_SCALE,
            # Scaled by the fixed-epsilon assumption the old caller made, so
            # "1.0" means "exactly as bad as Binomial v1 assumed".
            "mean_error": _masked_mean(error_probability, mask) / 0.01,
            "fraction_q30": _masked_fraction(is_q30, mask),
            "forward_fraction": 1.0 - _masked_mean(strand, mask),
            "mean_read_position": _masked_mean(read_position, mask),
            "mean_end_distance": _masked_mean(end_distance, mask),
            "near_end_fraction": _masked_fraction((flags & FLAG_NEAR_END) > 0, mask),
            "mean_mapq": _masked_mean(mapq, mask) / MAPPING_QUALITY_SCALE,
        }

    alt = side(supports_alt)
    ref = side(supports_reference)

    columns = {
        "log_n": np.log1p(n) / DEPTH_SCALE,
        "log_k": np.log1p(k) / DEPTH_SCALE,
        "log_ref_count": np.log1p(reference_count) / DEPTH_SCALE,
        "log_other_count": np.log1p(other_count) / DEPTH_SCALE,
        "log_depth": np.log1p(depth) / DEPTH_SCALE,
        "log_gap": np.log1p(gap) / DEPTH_SCALE,
        "log_insertion": np.log1p(insertion) / DEPTH_SCALE,
        # A linear copy of k as well as its log: at k = 1, 2, 3 -- the entire
        # regime where the 5-19x decisions live -- the difference between
        # consecutive integers is the whole signal, and log1p compresses it.
        "k_small": np.minimum(k, K_SCALE) / K_SCALE,
        "vaf": np.where(n > 0, k / safe_n, 0.0),
        "ref_fraction": np.where(n > 0, reference_count / safe_n, 0.0),
        "other_fraction": np.where(n > 0, other_count / safe_n, 0.0),
        "gap_fraction": np.where(depth > 0, gap / np.maximum(depth, 1), 0.0),

        "alt_mean_phred": alt["mean_phred"],
        "alt_min_phred": alt["min_phred"],
        "alt_max_phred": alt["max_phred"],
        "alt_mean_error": alt["mean_error"],
        "alt_fraction_q30": alt["fraction_q30"],
        "alt_has_reads": (k > 0).astype(np.float64),

        "ref_mean_phred": ref["mean_phred"],
        "ref_min_phred": ref["min_phred"],
        "ref_max_phred": ref["max_phred"],
        "ref_mean_error": ref["mean_error"],
        "ref_fraction_q30": ref["fraction_q30"],

        "alt_minus_ref_phred": alt["mean_phred"] - ref["mean_phred"],
        "alt_minus_ref_error": alt["mean_error"] - ref["mean_error"],

        "alt_forward_fraction": alt["forward_fraction"],
        "alt_strand_bias": np.abs(alt["forward_fraction"] - 0.5) * 2.0,
        "alt_mean_read_position": alt["mean_read_position"],
        "alt_mean_end_distance": alt["mean_end_distance"],
        "alt_near_end_fraction": alt["near_end_fraction"],
        "alt_mean_mapq": alt["mean_mapq"],
        "ref_forward_fraction": ref["forward_fraction"],
        "ref_strand_bias": np.abs(ref["forward_fraction"] - 0.5) * 2.0,
        "ref_mean_read_position": ref["mean_read_position"],
        "ref_mean_end_distance": ref["mean_end_distance"],
        "ref_near_end_fraction": ref["near_end_fraction"],
        "ref_mean_mapq": ref["mean_mapq"],
        "alt_minus_ref_mapq": alt["mean_mapq"] - ref["mean_mapq"],
        "alt_minus_ref_near_end": alt["near_end_fraction"] - ref["near_end_fraction"],
    }

    for slot, token in enumerate(("N", "A", "C", "G", "T")):
        columns[f"ref_is_{token}"] = (reference_index == slot).astype(np.float64)

    features = np.empty((n_loci, FEATURE_DIM), dtype=np.float32)
    for position, name in enumerate(FEATURE_NAMES):
        features[:, position] = columns[name]

    if not np.all(np.isfinite(features)):
        raise ValueError("locus features contain a non-finite value")
    return features


def resolve_group_indices(groups: tuple[str, ...]) -> tuple[int, ...]:
    """Map feature-group names to their column indices."""
    return tuple(sorted(FEATURE_NAMES.index(name)
                        for group in groups for name in FEATURE_GROUPS[group]))


__all__ = [
    "FEATURE_NAMES", "FEATURE_DIM", "FEATURE_GROUPS",
    "build_locus_features", "resolve_group_indices",
]
