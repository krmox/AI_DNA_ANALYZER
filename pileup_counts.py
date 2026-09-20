"""Raw integer pileup counts, aligned 1:1 with ``raw_pileup.RawPileupProvider``.

Why a separate extractor
------------------------
``RawPileupProvider`` emits *fractions* (``frac_A`` ... ``frac_del``) with a
denominator of ``called + gap_reads``, plus a normalised ``log1p(depth)/5``.
Integer ``n`` and ``k`` cannot be recovered from those channels: ``depth``
counts reads that passed the MAPQ/duplicate filters including ones later
dropped by the base-quality filter, so ``depth`` is not the fraction
denominator and the product does not round-trip to counts. A binomial caller
needs exact ``n`` and ``k``, so it needs its own extractor.

This module subclasses :class:`~providers.GiabAlignmentProvider` and applies
the *identical* read-admission filters as ``RawPileupProvider``
(``min_mapping_quality``, ``min_base_quality``, duplicate / secondary /
supplementary / QC-fail exclusion) so the two views describe the same reads.
Window tiling, BED masking and label construction are inherited unchanged,
which guarantees locus-for-locus alignment with the neural model's inputs.

Inference-time information only: no VCF is opened while counting. Labels come
from the inherited ``_fetch_labels`` and are returned in a separate array that
the caller may use for evaluation only.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

from config import NUCLEOTIDES
from providers import GenomicWindow, GiabAlignmentProvider

#: Directory holding the optional ``pileup_native`` C/htslib extension (relative to this file).
_NATIVE_DIR = Path(__file__).resolve().parent / "native"
if _NATIVE_DIR.is_dir() and str(_NATIVE_DIR) not in sys.path:
    sys.path.insert(0, str(_NATIVE_DIR))

try:
    import pileup_native as _pileup_native  # type: ignore[import-not-found]
except ImportError:
    _pileup_native = None

#: ``AI_DNA_ANALYZER_DISABLE_NATIVE_PILEUP=1`` forces the pure-Python/pysam path.
NATIVE_AVAILABLE = _pileup_native is not None and not os.environ.get(
    "AI_DNA_ANALYZER_DISABLE_NATIVE_PILEUP"
)

#: Column order of the returned count matrix.
COUNT_COLUMNS: tuple[str, ...] = (
    "count_A", "count_C", "count_G", "count_T",
    "gap_reads",        # reads asserting a deletion at this locus
    "insertion_reads",  # reads carrying an insertion immediately after
    "depth",            # reads passing MAPQ/flag filters (base-quality agnostic)
    "quality_sum",      # summed Phred of the counted A/C/G/T observations
    "mapping_sum",      # summed MAPQ over the depth reads
    "reference_index",  # 0=N, 1=A, 2=C, 3=G, 4=T
)
COUNT_DIM = len(COUNT_COLUMNS)

_REFERENCE_TOKENS = ("N",) + NUCLEOTIDES


class PileupCountsProvider(GiabAlignmentProvider):
    """Serves windows as raw integer pileup counts."""

    def _count_locus(self, column: Any, reference_base: str) -> List[float]:
        """Pool one pileup column into raw integer counts.

        Args:
            column: A pysam pileup column.
            reference_base: The reference base at this locus.

        Returns:
            A list of :data:`COUNT_DIM` numbers following
            :data:`COUNT_COLUMNS`.
        """
        counts = {base: 0 for base in NUCLEOTIDES}
        gap_reads = 0
        insertion_reads = 0
        depth = 0
        quality_sum = 0.0
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

            if read.indel > 0:
                insertion_reads += 1

        reference_index = (
            _REFERENCE_TOKENS.index(reference_base.upper())
            if reference_base.upper() in _REFERENCE_TOKENS else 0
        )
        return [
            counts["A"], counts["C"], counts["G"], counts["T"],
            gap_reads, insertion_reads, depth, quality_sum, mapping_sum, reference_index,
        ]

    def window_counts(self, window: GenomicWindow, reference: str) -> np.ndarray:
        """Build the ``[seq_len, COUNT_DIM]`` count matrix for one window."""
        self._ensure_open()
        assert self._bam is not None and self._bam_contig is not None

        rows: List[List[float]] = []
        for offset in range(len(window)):
            base = reference[offset].upper()
            reference_index = (
                _REFERENCE_TOKENS.index(base) if base in _REFERENCE_TOKENS else 0
            )
            rows.append([0] * 9 + [reference_index])

        for column in self._bam.pileup(
            self._bam_contig, window.start, window.end,
            truncate=True, min_base_quality=0, stepper="samtools",
        ):
            offset = column.reference_pos - window.start
            if not 0 <= offset < len(window):
                continue
            rows[offset] = self._count_locus(column, reference[offset])

        return np.asarray(rows, dtype=np.float64)

    def __getitem__(self, index: int) -> Dict[str, np.ndarray]:  # type: ignore[override]
        """Return counts and labels for one window.

        Args:
            index: Positional index into the tiled windows.

        Returns:
            ``{"counts": [seq_len, COUNT_DIM], "labels": [seq_len]}``.

        Raises:
            IndexError: If ``index`` is out of range.
        """
        if self._windows is None:
            self._windows = self._build_windows()
        if not 0 <= index < len(self._windows):
            raise IndexError(f"index {index} out of range for {len(self._windows)} windows")

        window = self._windows[index]
        reference = self._fetch_reference(window)
        return {
            "counts": self.window_counts(window, reference),
            "labels": np.asarray(self._fetch_labels(window), dtype=np.int64),
            "positions": np.arange(window.start, window.end, dtype=np.int64),
        }


def load_counts(
    fasta: str, bam: str, vcf: str, bed: str | None, region: tuple[int, int],
    seq_len: int = 64, window_slice: slice | None = None, contig: str = "chr21",
    return_positions: bool = False,
) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Materialize a region's raw counts and labels.

    Args:
        fasta: Reference FASTA.
        bam: Alignment BAM.
        vcf: Truth VCF (labels only; never consulted while counting).
        bed: Optional high-confidence BED.
        region: ``(start, end)`` in the FASTA's coordinate system.
        seq_len: Window width.
        window_slice: Optional slice restricting which tiled windows are materialized.
        contig: Contig to tile.
        return_positions: Also return each row's true 0-based genomic coordinate.
            A row index is NOT a coordinate once ``_build_windows`` has dropped any
            window (BED gap / all-N slice), so ``chunk_start + row_index`` silently
            desyncs; this returns the coordinate from the surviving windows.

    Returns:
        ``(counts, labels)`` or, with ``return_positions``, ``(counts, labels, positions)``.
    """
    provider = PileupCountsProvider(
        fasta_path=fasta, bam_path=bam, vcf_path=vcf, contig=contig,
        region=region, seq_len=seq_len, high_confidence_bed=bed,
    )
    with provider:
        total = len(provider)
        indices = range(total)[window_slice] if window_slice is not None else range(total)
        windows = [provider._windows[i] for i in indices]  # type: ignore[index]

        if NATIVE_AVAILABLE and windows:
            counts, labels, positions = _native_window_batch(provider, windows)
        else:
            counts, labels, positions = [], [], []
            for window in windows:
                reference = provider._fetch_reference(window)
                counts.append(provider.window_counts(window, reference))
                labels.append(np.asarray(provider._fetch_labels(window), dtype=np.int64))
                if return_positions:
                    positions.append(np.arange(window.start, window.end, dtype=np.int64))

    counts_arr = np.concatenate(counts) if counts else np.zeros((0, COUNT_DIM), dtype=np.float64)
    labels_arr = np.concatenate(labels) if labels else np.zeros(0, dtype=np.int64)
    if return_positions:
        positions_arr = np.concatenate(positions) if positions else np.zeros(0, dtype=np.int64)
        return counts_arr, labels_arr, positions_arr
    return counts_arr, labels_arr


def _native_window_batch(
    provider: "PileupCountsProvider", windows: List[GenomicWindow],
) -> tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """Fill ``windows`` from ONE ``pileup_native.count_region`` call over the whole span
    (columns 0-8), attaching the FASTA-derived reference_index column (9)."""
    assert _pileup_native is not None
    span_start = min(w.start for w in windows)
    span_end = max(w.end for w in windows)
    buf = _pileup_native.count_region(
        str(provider.bam_path), provider._bam_contig, span_start, span_end,
        provider.min_mapping_quality, provider.min_base_quality,
    )
    block9 = np.frombuffer(buf, dtype=np.float64).reshape(span_end - span_start, 9)

    counts: List[np.ndarray] = []
    labels: List[np.ndarray] = []
    positions: List[np.ndarray] = []
    for window in windows:
        reference = provider._fetch_reference(window)
        ref_idx_col = np.fromiter(
            (
                _REFERENCE_TOKENS.index(base) if (base := reference[i].upper()) in _REFERENCE_TOKENS else 0
                for i in range(len(window))
            ),
            dtype=np.float64, count=len(window),
        )
        rows9 = block9[window.start - span_start: window.end - span_start]
        counts.append(np.concatenate([rows9, ref_idx_col[:, None]], axis=1))
        labels.append(np.asarray(provider._fetch_labels(window), dtype=np.int64))
        positions.append(np.arange(window.start, window.end, dtype=np.int64))
    return counts, labels, positions


__all__ = ["COUNT_COLUMNS", "COUNT_DIM", "PileupCountsProvider", "load_counts", "NATIVE_AVAILABLE"]
