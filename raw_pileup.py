"""Raw-pileup input representation: evidence in, no pre-called consensus.

Why this module exists
----------------------
``GiabAlignmentProvider`` hands the model ``input_ids``, a *consensus* base
produced by thresholding allele fraction at ``vaf_threshold`` (0.15). At any
realistic depth a sequencing error never clears that bar and a true
heterozygous site always does, so the consensus already *is* the variant
call. Measured on real HG002 chr21:30.0-30.47Mb: of 454,754 reference-matching
loci only 5 had ``input_ids != reference_ids``, and the rule "call SNP wherever
the two differ" scored recall 1.0000. The model was copying a decision made
for it upstream.

This provider removes that decision from the input path. Each locus becomes a
vector of **raw counts and quality statistics pooled from the reads**, plus the
reference identity. Nothing in it says "this position differs from the
reference"; separating a true variant from sequencing noise requires comparing
the observed base distribution against the reference channel, which is the
thing the model is supposed to learn.

What is deliberately NOT here
-----------------------------
* ``vaf_threshold`` is never consulted when building model input.
* ``consensus_base`` / ``consensus_insertion`` are never read.
* ``indel_support_threshold`` is never consulted when building model input.
* No VCF is opened while constructing features. Truth enters only through
  :meth:`GiabAlignmentProvider._fetch_labels`, which is called separately and
  cannot influence the tensor.
* Counts are kept in fixed ``A,C,G,T`` order and the reference is supplied as
  a separate one-hot. They are *not* reordered into "reference allele" and
  "alternate allele" slots, because that reordering is itself the comparison
  the model must perform.

One column per reference position
---------------------------------
The consensus path emitted extra gap-padded columns for insertions, which
required ``indel_support_threshold`` to decide an insertion had happened. Here
each window emits exactly ``seq_len`` columns, one per reference base, and
insertion evidence rides along as a per-locus feature. Labels already sit on
reference coordinates (``_classify_allele`` anchors insertions at the anchor
base and deletions on the deleted bases), so no label expansion is needed and
none is done.
"""

from __future__ import annotations

from typing import Any, Dict, List

import torch
from torch.utils.data import Dataset

from config import NUCLEOTIDES
from providers import GiabAlignmentProvider, GenomicWindow

#: Per-locus feature layout. Order is fixed and asserted by the sanity checks.
FEATURE_NAMES: tuple[str, ...] = (
    "frac_A", "frac_C", "frac_G", "frac_T",   # base composition of the read stack
    "frac_del",                                # reads asserting a deletion here
    "frac_ins",                                # reads carrying an insertion after here
    "log_depth",                               # log1p(depth) / 5
    "mean_base_quality",                       # mean Phred of counted bases / 40
    "mean_mapping_quality",                    # mean MAPQ of counted reads / 60
    "ref_is_N", "ref_is_A", "ref_is_C", "ref_is_G", "ref_is_T",
)
FEATURE_DIM: int = len(FEATURE_NAMES)

#: Normalisation constants. Chosen to put typical values near 1 without
#: clipping: depth spans ~20x (synthetic) to ~120x (real HG002), Illumina
#: Phred tops out near 40, novoalign MAPQ near 70.
DEPTH_SCALE = 5.0
BASE_QUALITY_SCALE = 40.0
MAPPING_QUALITY_SCALE = 60.0

_REFERENCE_TOKENS = ("N",) + NUCLEOTIDES  # N, A, C, G, T -> one-hot slots 9..13


class RawPileupProvider(GiabAlignmentProvider):
    """Serves windows as raw per-locus pileup evidence instead of a consensus.

    Inherits window tiling, BED confidence masking, contig-name reconciliation,
    reference fetching and label construction from
    :class:`~providers.GiabAlignmentProvider` unchanged; only the observed
    channel is replaced. The parent's ``vaf_threshold`` and
    ``indel_support_threshold`` attributes still exist but are never read by
    this class.
    """

    def _raw_locus_features(
        self, column: "Any", reference_base: str
    ) -> List[float]:
        """Pool one pileup column into raw evidence features.

        Applies the same read-admission filters as the parent
        (``min_mapping_quality``, ``min_base_quality``, duplicate/secondary/
        supplementary/QC-fail exclusion) so that read quality control is held
        constant against the previous experiment. Those filters decide *which
        reads count*, never *what the locus is*.

        Args:
            column: A pysam pileup column.
            reference_base: The reference base at this locus.

        Returns:
            A list of :data:`FEATURE_DIM` floats.
        """
        counts = {base: 0 for base in NUCLEOTIDES}
        gap_reads = 0
        insertion_reads = 0
        depth = 0
        quality_sum = 0.0
        quality_n = 0
        mapping_sum = 0.0

        for read in column.pileups:
            alignment = read.alignment

            if alignment.mapping_quality < self.min_mapping_quality:
                continue
            if alignment.is_duplicate or alignment.is_qcfail or alignment.is_unmapped:
                continue
            if alignment.is_secondary or alignment.is_supplementary:
                continue

            depth += 1
            mapping_sum += float(alignment.mapping_quality)

            if read.is_del or read.is_refskip:
                gap_reads += 1
                continue

            position = read.query_position
            if position is None:
                continue

            quality = alignment.query_qualities[position] if alignment.query_qualities else 0
            if quality < self.min_base_quality:
                continue

            base = alignment.query_sequence[position].upper()
            if base not in counts:
                continue

            counts[base] += 1
            quality_sum += float(quality)
            quality_n += 1

            if read.indel > 0:
                insertion_reads += 1

        called = sum(counts.values())
        observed = called + gap_reads

        if observed > 0:
            fractions = [counts[base] / observed for base in NUCLEOTIDES]
            frac_del = gap_reads / observed
            frac_ins = insertion_reads / observed
        else:
            fractions = [0.0, 0.0, 0.0, 0.0]
            frac_del = 0.0
            frac_ins = 0.0

        mean_quality = (quality_sum / quality_n / BASE_QUALITY_SCALE) if quality_n else 0.0
        mean_mapping = (mapping_sum / depth / MAPPING_QUALITY_SCALE) if depth else 0.0

        reference_onehot = [
            1.0 if reference_base.upper() == token else 0.0 for token in _REFERENCE_TOKENS
        ]
        # A reference base outside N/ACGT (never seen on this data) would give
        # an all-zero one-hot rather than silently aliasing onto another base.

        return (
            fractions
            + [frac_del, frac_ins]
            + [torch.log1p(torch.tensor(float(depth))).item() / DEPTH_SCALE]
            + [mean_quality, mean_mapping]
            + reference_onehot
        )

    def _window_features(self, window: GenomicWindow, reference: str) -> torch.Tensor:
        """Build the ``[seq_len, FEATURE_DIM]`` evidence tensor for a window.

        Loci with no coverage keep an all-zero evidence block and only their
        reference one-hot set, which reads as "no evidence" rather than as
        evidence of a reference match.

        Args:
            window: Interval to build.
            reference: Reference sequence for the window.

        Returns:
            Float tensor ``[len(window), FEATURE_DIM]``.
        """
        self._ensure_open()
        assert self._bam is not None and self._bam_contig is not None

        rows: List[List[float]] = []
        for offset in range(len(window)):
            base = reference[offset]
            reference_onehot = [
                1.0 if base.upper() == token else 0.0 for token in _REFERENCE_TOKENS
            ]
            rows.append([0.0] * 9 + reference_onehot)

        for column in self._bam.pileup(
            self._bam_contig,
            window.start,
            window.end,
            truncate=True,
            min_base_quality=0,
            stepper="samtools",
        ):
            offset = column.reference_pos - window.start
            if not 0 <= offset < len(window):
                continue
            rows[offset] = self._raw_locus_features(column, reference[offset])

        return torch.tensor(rows, dtype=torch.float)

    def __getitem__(self, index: int) -> Dict[str, torch.Tensor]:  # type: ignore[override]
        """Assemble one window as raw evidence plus labels.

        Args:
            index: Positional index into the tiled windows.

        Returns:
            ``{"pileup_features": [seq_len, FEATURE_DIM], "labels": [seq_len],
            "depth": [seq_len]}``. ``depth`` is carried for diagnostics only;
            the model's forward signature does not accept it.

        Raises:
            IndexError: If ``index`` is out of range.
        """
        if self._windows is None:
            self._windows = self._build_windows()
        if not 0 <= index < len(self._windows):
            raise IndexError(f"index {index} out of range for {len(self._windows)} windows")

        window = self._windows[index]
        reference = self._fetch_reference(window)

        features = self._window_features(window, reference)
        labels = torch.tensor(self._fetch_labels(window), dtype=torch.long)

        # Recover raw depth from the normalised channel for diagnostics.
        depth = torch.expm1(features[:, 6] * DEPTH_SCALE)

        return {"pileup_features": features, "labels": labels, "depth": depth}


class RawPileupDataset(Dataset):
    """Thin ``Dataset`` over a :class:`RawPileupProvider`."""

    def __init__(self, provider: RawPileupProvider) -> None:
        """Wrap a provider.

        Args:
            provider: The provider to serve windows from.
        """
        self.provider = provider

    def __len__(self) -> int:
        """Return the number of windows."""
        return len(self.provider)

    def __getitem__(self, index: int) -> Dict[str, torch.Tensor]:
        """Return the window at ``index``."""
        return self.provider[index]


__all__ = [
    "FEATURE_NAMES",
    "FEATURE_DIM",
    "RawPileupProvider",
    "RawPileupDataset",
]
