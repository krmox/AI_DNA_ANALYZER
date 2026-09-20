"""Read-level pileup extraction: individual observations, not summary statistics.

The 14-channel aggregate representation compresses every read at a locus into
fractions and means. That compression is exactly what the binomial caller also
does -- it sees only ``(n, k, epsilon)`` -- so an aggregate neural model cannot
in principle express "10 alternate reads out of 30, but all 10 are low-quality,
clustered at read ends, and on one strand". This module keeps the individual
observations so that question can be asked.

Layout
------
Each locus becomes an ``[R, C]`` matrix of ``uint8`` columns (see
:data:`READ_COLUMNS`) plus a validity mask. ``R = 48`` is chosen from the data,
not by convention: maximum observed depth is 43 (train) and 38 (test) under the
project's MAPQ/flag filters, so **every locus in both regions fits without
truncation** and no observation is discarded. The deterministic sampling path
is implemented and unit-tested regardless, so the representation stays correct
on deeper data, but it never activates here.

Storage is ``uint8`` because a float32 ``[N, R, F]`` tensor for the 1.88 M-locus
training region would be ~7 GB. The compact columns are ~720 MB and expand to
float features per batch (:func:`build_features`).

Read admission
--------------
Identical to :class:`pileup_counts.PileupCountsProvider`: ``min_mapping_quality``,
and exclusion of duplicate / QC-fail / unmapped / secondary / supplementary
reads. Every admitted read becomes a row, *including* reads whose base fails the
base-quality filter -- those are dropped entirely by the binomial caller, so
retaining them (flagged) is one of the concrete pieces of evidence this
representation adds. The set of admitted reads is exactly the ``depth`` column
of the aggregate representation, which is what makes
:func:`reconstruct_counts` an exact check rather than an approximate one.

Read order
----------
Rows are ordered by a deterministic 64-bit hash of the read name. The hash is
stable across runs and machines (BLAKE2b, not Python's salted ``hash``), and is
statistically independent of every feature, so BAM traversal order cannot leak
into the model and the ordering itself carries no information. When depth
exceeds ``R``, the ``R`` smallest hashes are kept -- a deterministic,
feature-independent subsample rather than "whatever pysam yielded first".

Leakage boundary
----------------
No VCF is opened during extraction. Labels come from the inherited
``_fetch_labels`` and are returned in a separate array for evaluation only. No
feature encodes "differs from reference": the reference base is never consulted
while building read rows, and the candidate alternate allele, where one is
needed downstream, is derived from observed counts alone.
"""

from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

from config import NUCLEOTIDES
from providers import GenomicWindow, GiabAlignmentProvider

_NATIVE_DIR = Path(__file__).resolve().parent / "native"
if _NATIVE_DIR.is_dir() and str(_NATIVE_DIR) not in sys.path:
    sys.path.insert(0, str(_NATIVE_DIR))
try:
    import reads_native as _reads_native  # type: ignore[import-not-found]
except ImportError:
    _reads_native = None

#: ``AI_DNA_ANALYZER_DISABLE_NATIVE_READS=1`` (or the shared ``..._NATIVE_PILEUP=1``) forces the
#: original pure-Python ``window_reads`` path, the reference for the equivalence gate.
NATIVE_READS_AVAILABLE = _reads_native is not None and not (
    os.environ.get("AI_DNA_ANALYZER_DISABLE_NATIVE_READS")
    or os.environ.get("AI_DNA_ANALYZER_DISABLE_NATIVE_PILEUP"))

#: Reads retained per locus. Justified from the data: max depth is 43 (train)
#: and 38 (test), so R=48 truncates nothing in either region.
MAX_READS = 48

#: Compact per-read storage columns, all uint8.
READ_COLUMNS: tuple[str, ...] = (
    "base_code",        # 0=A 1=C 2=G 3=T 4=gap(del/refskip) 5=other/N 255=pad
    "base_quality",     # Phred, 0 when the read has no base here
    "mapping_quality",  # Phred, clipped to 255
    "strand",           # 0 = forward, 1 = reverse
    "read_position",    # round(255 * query_position / query_length)
    "end_distance",     # round(255 * min(pos, len-1-pos) / (len/2)), 255 = centre
    "flags",            # see the FLAG_* constants
    "mismatch_density",  # round(255 * 10 * NM / query_length), clipped
)
READ_DIM = len(READ_COLUMNS)

#: Bit positions in the ``flags`` column.
FLAG_VALID = 1 << 0        # row holds a real read (not padding)
FLAG_COUNTED = 1 << 1      # base is A/C/G/T and passed the base-quality filter
FLAG_GAP = 1 << 2          # read asserts a deletion / ref-skip at this locus
FLAG_INSERTION = 1 << 3    # read carries an insertion immediately after
FLAG_DELETION_NEXT = 1 << 4  # read carries a deletion immediately after
FLAG_NEAR_END = 1 << 5     # observation within NEAR_END_BASES of a read end

#: Distance from a read end below which an observation is flagged as terminal.
#: 5 bp is the conventional cutoff for end-repair / adapter artefacts.
NEAR_END_BASES = 5

_BASE_CODES = {base: index for index, base in enumerate(NUCLEOTIDES)}
_GAP_CODE, _OTHER_CODE, _PAD_CODE = 4, 5, 255


def read_sort_key(read_name: str) -> int:
    """Deterministic, feature-independent 64-bit ordering key for a read.

    Python's built-in ``hash`` is salted per process and would make extraction
    non-reproducible across runs, which is exactly the failure this function
    exists to prevent.

    Args:
        read_name: The read's query name.

    Returns:
        A stable integer in ``[0, 2**64)``.
    """
    digest = hashlib.blake2b(read_name.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big")


class ReadLevelPileupProvider(GiabAlignmentProvider):
    """Serves windows as ``[seq_len, MAX_READS, READ_DIM]`` read matrices."""

    def __init__(self, *args, max_reads: int = MAX_READS, **kwargs) -> None:
        """Configure the provider.

        Args:
            *args: Forwarded to :class:`GiabAlignmentProvider`.
            max_reads: Rows retained per locus.
            **kwargs: Forwarded to :class:`GiabAlignmentProvider`.
        """
        super().__init__(*args, **kwargs)
        self.max_reads = max_reads

    def _read_row(self, read: Any) -> List[int] | None:
        """Encode one pileup read into its compact uint8 row.

        Args:
            read: A pysam ``PileupRead``.

        Returns:
            The row, or None if the read fails admission.
        """
        alignment = read.alignment

        # Admission filters, identical to PileupCountsProvider._count_locus.
        if alignment.mapping_quality < self.min_mapping_quality:
            return None
        if alignment.is_duplicate or alignment.is_qcfail or alignment.is_unmapped:
            return None
        if alignment.is_secondary or alignment.is_supplementary:
            return None

        flags = FLAG_VALID
        if read.indel > 0:
            flags |= FLAG_INSERTION
        elif read.indel < 0:
            flags |= FLAG_DELETION_NEXT

        mapping_quality = min(int(alignment.mapping_quality), 255)
        strand = 1 if alignment.is_reverse else 0

        length = alignment.query_length or (
            len(alignment.query_sequence) if alignment.query_sequence else 0)
        try:
            mismatches = int(alignment.get_tag("NM"))
        except (KeyError, ValueError):
            mismatches = 0
        mismatch_density = min(255, int(round(2550.0 * mismatches / max(length, 1))))

        if read.is_del or read.is_refskip:
            # No base at this locus; position within the read is undefined.
            return [_GAP_CODE, 0, mapping_quality, strand, 0, 0,
                    flags | FLAG_GAP, mismatch_density]

        position = read.query_position
        if position is None:
            return [_OTHER_CODE, 0, mapping_quality, strand, 0, 0, flags, mismatch_density]

        quality = int(alignment.query_qualities[position]) if alignment.query_qualities else 0
        base = alignment.query_sequence[position].upper()
        code = _BASE_CODES.get(base, _OTHER_CODE)

        # FLAG_COUNTED reproduces the aggregate caller's admission decision:
        # an A/C/G/T base whose Phred clears the threshold. Bases failing it
        # are kept as rows -- invisible to the binomial caller, visible here.
        if code <= 3 and quality >= self.min_base_quality:
            flags |= FLAG_COUNTED

        end_distance = min(position, max(length - 1 - position, 0))
        if end_distance < NEAR_END_BASES:
            flags |= FLAG_NEAR_END

        return [
            code, min(quality, 255), mapping_quality, strand,
            int(round(255.0 * position / max(length, 1))),
            min(255, int(round(255.0 * end_distance / max(length / 2.0, 1.0)))),
            flags, mismatch_density,
        ]

    def _locus_reads(self, column: Any) -> np.ndarray:
        """Build one locus' ``[max_reads, READ_DIM]`` matrix.

        Args:
            column: A pysam pileup column.

        Returns:
            uint8 array; unused rows are padded with ``_PAD_CODE`` bases and
            cleared flags.
        """
        rows: List[List[int]] = []
        keys: List[int] = []
        for read in column.pileups:
            row = self._read_row(read)
            if row is not None:
                rows.append(row)
                keys.append(read_sort_key(read.alignment.query_name))

        matrix = self._empty_locus()
        if not rows:
            return matrix

        # Deterministic order; on overflow keep the R smallest hashes, which is
        # a reproducible feature-independent subsample.
        order = np.argsort(np.asarray(keys, dtype=np.uint64), kind="stable")
        selected = order[: self.max_reads]
        matrix[: selected.size] = np.asarray(rows, dtype=np.uint8)[selected]
        return matrix

    def _empty_locus(self) -> np.ndarray:
        """A fully padded locus matrix."""
        matrix = np.zeros((self.max_reads, READ_DIM), dtype=np.uint8)
        matrix[:, 0] = _PAD_CODE
        return matrix

    def window_reads(self, window: GenomicWindow) -> np.ndarray:
        """Build the ``[seq_len, max_reads, READ_DIM]`` matrix for one window."""
        self._ensure_open()
        assert self._bam is not None and self._bam_contig is not None

        matrix = np.stack([self._empty_locus() for _ in range(len(window))])
        for column in self._bam.pileup(
            self._bam_contig, window.start, window.end,
            truncate=True, min_base_quality=0, stepper="samtools",
        ):
            offset = column.reference_pos - window.start
            if 0 <= offset < len(window):
                matrix[offset] = self._locus_reads(column)
        return matrix

    def __getitem__(self, index: int) -> Dict[str, np.ndarray]:  # type: ignore[override]
        """Return read matrices and labels for one window.

        Args:
            index: Positional index into the tiled windows.

        Returns:
            ``{"reads": [seq_len, R, READ_DIM], "labels": [seq_len]}``.

        Raises:
            IndexError: If ``index`` is out of range.
        """
        if self._windows is None:
            self._windows = self._build_windows()
        if not 0 <= index < len(self._windows):
            raise IndexError(f"index {index} out of range for {len(self._windows)} windows")
        window = self._windows[index]
        return {"reads": self.window_reads(window),
                "labels": np.asarray(self._fetch_labels(window), dtype=np.int64)}


def reconstruct_counts(reads: np.ndarray, reference_index: np.ndarray) -> np.ndarray:
    """Rebuild the aggregate count matrix from the read-level tensor.

    This is the correctness check the whole experiment rests on: if the
    read-level tensor is a faithful superset of the aggregate representation,
    the aggregate must fall out of it exactly. Any mismatch means the two views
    describe different reads and no comparison between them would be valid.

    Args:
        reads: ``[N, R, READ_DIM]`` uint8 read matrices.
        reference_index: ``[N]`` reference base index (0=N, 1=A..4=T), carried
            through unchanged; it is not used to build any read feature.

    Returns:
        ``[N, 10]`` float64 matrix following ``pileup_counts.COUNT_COLUMNS``.
    """
    reads = np.asarray(reads)
    flags = reads[..., 6].astype(np.int64)
    valid = (flags & FLAG_VALID) > 0
    counted = (flags & FLAG_COUNTED) > 0
    gap = (flags & FLAG_GAP) > 0
    code = reads[..., 0].astype(np.int64)

    out = np.zeros((reads.shape[0], 10), dtype=np.float64)
    for base in range(4):
        out[:, base] = np.sum(counted & (code == base), axis=1)
    out[:, 4] = np.sum(gap, axis=1)
    out[:, 5] = np.sum(counted & ((flags & FLAG_INSERTION) > 0), axis=1)
    out[:, 6] = np.sum(valid, axis=1)
    out[:, 7] = np.sum(np.where(counted, reads[..., 1].astype(np.float64), 0.0), axis=1)
    out[:, 8] = np.sum(np.where(valid, reads[..., 2].astype(np.float64), 0.0), axis=1)
    out[:, 9] = reference_index
    return out


#: Names of the float features produced by :func:`build_features`, in order.
FEATURE_NAMES: tuple[str, ...] = (
    "is_A", "is_C", "is_G", "is_T", "is_gap",
    "base_quality", "mapping_quality", "strand",
    "read_position", "end_distance", "near_end",
    "insertion_next", "deletion_next", "low_quality_base",
)
FEATURE_DIM = len(FEATURE_NAMES)

#: Feature groups the ablations switch off. ``mismatch_density`` is extracted
#: and cached but deliberately excluded from the primary feature set: it is
#: alignment-derived rather than truth-derived, but it partially encodes
#: reference-discordance, which the protocol restricts, so it is kept available
#: for a dedicated follow-up rather than folded into the headline result.
FEATURE_GROUPS: dict[str, tuple[str, ...]] = {
    "quality": ("base_quality", "mapping_quality", "low_quality_base"),
    "strand": ("strand",),
    "position": ("read_position", "end_distance", "near_end"),
    "indel": ("insertion_next", "deletion_next"),
}

BASE_QUALITY_SCALE = 60.0
MAPPING_QUALITY_SCALE = 60.0


def build_features(reads: np.ndarray, drop_groups: tuple[str, ...] = ()) -> np.ndarray:
    """Expand compact uint8 rows into the float per-read feature tensor.

    Args:
        reads: ``[..., R, READ_DIM]`` uint8 array.
        drop_groups: Names from :data:`FEATURE_GROUPS` to zero out, used by the
            ablations. Zeroing rather than removing keeps the tensor shape and
            therefore the parameter count identical across ablations, so an
            ablation measures the information and not the capacity.

    Returns:
        ``[..., R, FEATURE_DIM]`` float32 features.
    """
    reads = np.asarray(reads)
    flags = reads[..., 6].astype(np.int64)
    code = reads[..., 0].astype(np.int64)
    valid = (flags & FLAG_VALID) > 0

    features = np.zeros(reads.shape[:-1] + (FEATURE_DIM,), dtype=np.float32)
    for base in range(4):
        features[..., base] = (code == base) & valid
    features[..., 4] = (flags & FLAG_GAP) > 0
    features[..., 5] = reads[..., 1].astype(np.float32) / BASE_QUALITY_SCALE
    features[..., 6] = reads[..., 2].astype(np.float32) / MAPPING_QUALITY_SCALE
    features[..., 7] = np.where(valid, reads[..., 3].astype(np.float32), 0.0)
    features[..., 8] = np.where(valid, reads[..., 4].astype(np.float32) / 255.0, 0.0)
    features[..., 9] = np.where(valid, reads[..., 5].astype(np.float32) / 255.0, 0.0)
    features[..., 10] = (flags & FLAG_NEAR_END) > 0
    features[..., 11] = (flags & FLAG_INSERTION) > 0
    features[..., 12] = (flags & FLAG_DELETION_NEXT) > 0
    # A base present but not counted: invisible to the aggregate caller.
    features[..., 13] = valid & (code <= 3) & ((flags & FLAG_COUNTED) == 0)

    for group in drop_groups:
        for name in FEATURE_GROUPS[group]:
            features[..., FEATURE_NAMES.index(name)] = 0.0

    features *= valid[..., None].astype(np.float32)
    return features


def read_mask(reads: np.ndarray) -> np.ndarray:
    """Boolean validity mask ``[..., R]`` for the read rows."""
    return (np.asarray(reads)[..., 6].astype(np.int64) & FLAG_VALID) > 0


def build_features_torch(reads: "torch.Tensor", drop_indices: tuple[int, ...] = ()):
    """GPU-side twin of :func:`build_features`, for per-batch expansion.

    The cached ``uint8`` tensor is ~720 MB for the training region; its float32
    expansion would be ~7 GB, so the conversion happens per batch on the
    device. This function must agree with :func:`build_features` exactly --
    ``test_read_level.py`` asserts bit equality between the two, because a
    silent divergence would mean the model trains on different data than the
    sanity checks validated.

    Args:
        reads: ``[..., R, READ_DIM]`` uint8 tensor.
        drop_indices: Feature indices to zero (ablations), pre-resolved from
            :data:`FEATURE_GROUPS` so the hot path does no name lookups.

    Returns:
        ``[..., R, FEATURE_DIM]`` float32 tensor on the same device.
    """
    import torch

    flags = reads[..., 6].to(torch.int64)
    code = reads[..., 0].to(torch.int64)
    valid = (flags & FLAG_VALID) > 0
    valid_f = valid.to(torch.float32)

    channels = [
        ((code == 0) & valid).to(torch.float32),
        ((code == 1) & valid).to(torch.float32),
        ((code == 2) & valid).to(torch.float32),
        ((code == 3) & valid).to(torch.float32),
        ((flags & FLAG_GAP) > 0).to(torch.float32),
        reads[..., 1].to(torch.float32) / BASE_QUALITY_SCALE,
        reads[..., 2].to(torch.float32) / MAPPING_QUALITY_SCALE,
        reads[..., 3].to(torch.float32) * valid_f,
        reads[..., 4].to(torch.float32) / 255.0 * valid_f,
        reads[..., 5].to(torch.float32) / 255.0 * valid_f,
        ((flags & FLAG_NEAR_END) > 0).to(torch.float32),
        ((flags & FLAG_INSERTION) > 0).to(torch.float32),
        ((flags & FLAG_DELETION_NEXT) > 0).to(torch.float32),
        (valid & (code <= 3) & ((flags & FLAG_COUNTED) == 0)).to(torch.float32),
    ]
    features = torch.stack(channels, dim=-1)
    for index in drop_indices:
        features[..., index] = 0.0
    return features * valid_f.unsqueeze(-1)


def resolve_drop_indices(drop_groups: tuple[str, ...]) -> tuple[int, ...]:
    """Map ablation group names to feature indices."""
    return tuple(sorted(FEATURE_NAMES.index(name)
                        for group in drop_groups for name in FEATURE_GROUPS[group]))


def load_reads(fasta: str, bam: str, vcf: str, bed: str | None,
               region: tuple[int, int], seq_len: int = 64,
               max_reads: int = MAX_READS,
               contig: str = "chr21") -> tuple[np.ndarray, np.ndarray]:
    """Materialize a region's read-level tensor and labels.

    Args:
        fasta: Reference FASTA.
        bam: Alignment BAM.
        vcf: Truth VCF (labels only; never consulted while extracting reads).
        bed: Optional high-confidence BED.
        region: ``(start, end)`` coordinates.
        seq_len: Window width.
        max_reads: Rows per locus.
        contig: Contig to tile. Defaults to ``"chr21"`` for backward
            compatibility; see :func:`pileup_counts.load_counts`.

    Returns:
        ``(reads [N, R, READ_DIM] uint8, labels [N])``.
    """
    provider = ReadLevelPileupProvider(
        fasta_path=fasta, bam_path=bam, vcf_path=vcf, contig=contig,
        region=region, seq_len=seq_len, high_confidence_bed=bed, max_reads=max_reads)
    with provider:
        if NATIVE_READS_AVAILABLE and len(provider):
            try:
                return _native_load_reads(provider)
            except RuntimeError:
                # The C path refuses inputs it cannot reproduce bit-for-bit
                # (l_qseq == 0, non-integer NM); fall back rather than diverge.
                pass
        reads, labels = [], []
        for index in range(len(provider)):
            window = provider[index]
            reads.append(window["reads"])
            labels.append(window["labels"])
    return np.concatenate(reads), np.concatenate(labels)


def _native_threads() -> int:
    """Worker threads inside the C call (output is bit-identical for any count)."""
    return max(1, int(os.environ.get("AI_DNA_ANALYZER_READ_THREADS", "1")))


def _native_load_reads(provider: ReadLevelPileupProvider) -> tuple[np.ndarray, np.ndarray]:
    """Fill every window from ``reads_native.read_windows``: the same per-window pileup passes
    ``window_reads`` makes, run in C over one open BAM handle."""
    assert _reads_native is not None and provider._windows is not None
    windows = provider._windows
    starts = np.fromiter((w.start for w in windows), dtype=np.int64, count=len(windows))
    buf = _reads_native.read_windows(
        str(provider.bam_path), provider._bam_contig, starts, provider.seq_len,
        provider.min_mapping_quality, provider.min_base_quality, provider.max_reads,
        _native_threads())
    reads = np.frombuffer(buf, dtype=np.uint8).reshape(
        len(windows) * provider.seq_len, provider.max_reads, READ_DIM)
    labels = np.concatenate([np.asarray(provider._fetch_labels(w), dtype=np.int64) for w in windows])
    return reads, labels


__all__ = [
    "MAX_READS", "READ_COLUMNS", "READ_DIM", "FEATURE_NAMES", "FEATURE_DIM",
    "FEATURE_GROUPS", "NEAR_END_BASES", "ReadLevelPileupProvider", "build_features",
    "build_features_torch", "resolve_drop_indices",
    "read_mask", "reconstruct_counts", "load_reads", "read_sort_key",
    "FLAG_VALID", "FLAG_COUNTED", "FLAG_GAP", "FLAG_INSERTION", "FLAG_DELETION_NEXT",
    "FLAG_NEAR_END",
]
