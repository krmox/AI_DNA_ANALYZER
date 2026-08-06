"""Data provider interfaces for synthetic and real sequencing input.

A single contract — :class:`VariantDataProvider` — is satisfied by both the
synthetic simulator and the real GIAB BAM/VCF/FASTA pipeline, so everything
downstream consumes windows without knowing which produced them.

Coordinate convention
---------------------
Internally **every** coordinate is 0-based half-open, ``[start, end)``. This
matters because the three file formats disagree:

* FASTA via ``pysam.FastaFile.fetch`` — 0-based half-open. Matches.
* BAM via ``pysam`` ``reference_start`` / pileup ``reference_pos`` — 0-based.
  Matches.
* VCF — **1-based inclusive** in the file. ``pysam.VariantRecord.start`` is
  already converted to 0-based by pysam, while ``.pos`` is the raw 1-based
  value. This module uses ``.start`` exclusively and never ``.pos``, which
  removes the usual off-by-one class of bug at the source.

Pileup consensus
----------------
The single-read representation made depth and VAF meaningless, because one
read cannot tell you how many reads agree. :meth:`GiabAlignmentProvider._fetch_observed`
now builds a consensus across all reads covering each locus, which is what
makes the depth channel carry signal: a true variant recurs across reads at
one locus, a sequencing error does not. That difference is exactly the
allele fraction.

What real data will break that synthetic data did not
-----------------------------------------------------
Stated up front so the first real run is not a surprise:

* Class density falls from 1:20 to roughly 1:1000. Loss weighting calibrated
  on synthetic density will not transfer; run
  ``summarise_label_distribution.py`` and recalibrate before training.
* Roughly half of chr21 is repetitive, and its acrocentric short arm is
  effectively uncallable. Windows there produce low-MAPQ pileups and are the
  main source of false positives.
* GIAB high-confidence BED coverage is not the whole chromosome. Regions
  outside it are not reliably labelled, so training on them injects false
  negatives. Pass ``high_confidence_bed`` unless you have a reason not to.
"""

from __future__ import annotations

import gzip
from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterator, List, Sequence, Tuple

import torch

from .config import (
    GAP_TOKEN,
    LABEL_DELETION,
    LABEL_INSERTION,
    LABEL_NORMAL,
    LABEL_SNP,
    NUCLEOTIDES,
)
from .dataset import DNATokenizer, SyntheticVariantDataset, VariantSample

if TYPE_CHECKING:  # pragma: no cover - typing only
    import pysam


@dataclass(frozen=True)
class GenomicWindow:
    """A half-open reference interval on one contig.

    Attributes:
        contig: Reference contig name, e.g. ``"chr21"``.
        start: 0-based inclusive start coordinate.
        end: 0-based exclusive end coordinate.
    """

    contig: str
    start: int
    end: int

    def __post_init__(self) -> None:
        """Validate interval invariants.

        Raises:
            ValueError: If the interval is empty, reversed or negative.
        """
        if self.start < 0:
            raise ValueError("start must be non-negative")
        if self.end <= self.start:
            raise ValueError("end must be strictly greater than start")

    def __len__(self) -> int:
        """Return the window width in reference bases."""
        return self.end - self.start


@dataclass
class LocusPileup:
    """Consensus summary of all reads covering one reference locus.

    Attributes:
        depth: Number of reads passing filters at this locus.
        base_counts: Read counts per observed base, including ``-`` for
            deletion gaps.
        consensus_base: The called base after applying the VAF threshold.
        reference_base: The reference base at this locus.
        vaf: Fraction of passing reads carrying a non-reference allele.
        quality: Mean Phred of the reads supporting ``consensus_base``.
        insertion_sequences: Inserted sequences observed immediately after
            this locus, one entry per supporting read.
    """

    depth: int
    base_counts: Dict[str, int]
    consensus_base: str
    reference_base: str
    vaf: float
    quality: float
    insertion_sequences: List[str]

    @property
    def insertion_support(self) -> float:
        """Fraction of covering reads showing an insertion after this locus.

        Returns:
            Support in ``[0, 1]``; 0 when depth is 0.
        """
        if self.depth <= 0:
            return 0.0
        return len(self.insertion_sequences) / self.depth

    def consensus_insertion(self) -> str:
        """Return the most frequently observed inserted sequence.

        Returns:
            The modal inserted sequence, or an empty string if none.
        """
        if not self.insertion_sequences:
            return ""
        return Counter(self.insertion_sequences).most_common(1)[0][0]


@dataclass
class WindowTensors:
    """Model-ready tensors for one window, plus provenance.

    The five tensor fields match :class:`~.dataset.VariantSample` exactly, so
    a provider's output drops into the existing training loop unchanged.

    Attributes:
        reference_ids: Long tensor ``[seq_len]`` of reference column tokens.
        input_ids: Long tensor ``[seq_len]`` of observed column tokens.
        labels: Long tensor ``[seq_len]`` of per-column class ids.
        base_quality: Float tensor ``[seq_len]`` of Phred scores.
        depth: Float tensor ``[seq_len]`` of read depth.
        vaf: Optional float tensor ``[seq_len]`` of variant allele fractions.
            Populated from real pileups; ``None`` for synthetic data.
        window: The genomic interval this came from, when known. Carried so
            predictions can be written back to VCF coordinates.
    """

    reference_ids: torch.Tensor
    input_ids: torch.Tensor
    labels: torch.Tensor
    base_quality: torch.Tensor
    depth: torch.Tensor
    vaf: torch.Tensor | None = None
    window: GenomicWindow | None = None

    def to_sample(self) -> VariantSample:
        """Reduce to the five-tensor form the model consumes.

        VAF is retained on this object rather than discarded at source, so
        adding a VAF input channel later is a local change to
        :class:`~.model.InputFusion` alone.

        Returns:
            A :class:`~.dataset.VariantSample`.
        """
        return VariantSample(
            reference_ids=self.reference_ids,
            input_ids=self.input_ids,
            labels=self.labels,
            base_quality=self.base_quality,
            depth=self.depth,
        )


class VariantDataProvider(ABC):
    """Contract for anything that yields fixed-width labelled windows.

    Implementations must guarantee that every emitted tensor has exactly
    ``seq_len`` elements and that labels use the package's class ids.
    """

    def __init__(self, seq_len: int) -> None:
        """Initialise the provider.

        Args:
            seq_len: Fixed window width in tokens.

        Raises:
            ValueError: If ``seq_len`` is not positive.
        """
        if seq_len <= 0:
            raise ValueError("seq_len must be positive")
        self.seq_len = seq_len
        self.tokenizer = DNATokenizer()

    @abstractmethod
    def __len__(self) -> int:
        """Return the number of windows available."""

    @abstractmethod
    def __getitem__(self, index: int) -> WindowTensors:
        """Return the window at ``index``."""

    def __iter__(self) -> Iterator[WindowTensors]:
        """Iterate over all windows in order.

        Yields:
            One :class:`WindowTensors` per window.
        """
        for index in range(len(self)):
            yield self[index]

    def _validate_widths(self, tensors: WindowTensors) -> WindowTensors:
        """Assert that a window matches the declared width.

        Args:
            tensors: The window to check.

        Returns:
            The same object, unchanged, for call chaining.

        Raises:
            ValueError: If any required tensor has the wrong length.
        """
        for name in ("reference_ids", "input_ids", "labels", "base_quality", "depth"):
            tensor = getattr(tensors, name)
            if tensor.shape != (self.seq_len,):
                raise ValueError(
                    f"{name} has shape {tuple(tensor.shape)}, expected ({self.seq_len},)"
                )
        return tensors


class SyntheticProvider(VariantDataProvider):
    """Adapter exposing the synthetic simulator through the provider contract.

    Exists so that swapping simulated for real data is a constructor change
    at the call site rather than a rewrite of the training loop.
    """

    def __init__(self, dataset: SyntheticVariantDataset) -> None:
        """Wrap an existing synthetic dataset.

        Args:
            dataset: The simulator to adapt.
        """
        super().__init__(dataset.seq_len)
        self.dataset = dataset

    def __len__(self) -> int:
        """Return the number of simulated windows."""
        return len(self.dataset)

    def __getitem__(self, index: int) -> WindowTensors:
        """Generate one simulated window.

        Args:
            index: Positional index.

        Returns:
            A :class:`WindowTensors` with ``vaf`` unset, since a single-read
            simulator has no allele fraction to report.
        """
        sample = self.dataset[index]
        return self._validate_widths(
            WindowTensors(
                reference_ids=sample["reference_ids"],
                input_ids=sample["input_ids"],
                labels=sample["labels"],
                base_quality=sample["base_quality"],
                depth=sample["depth"],
            )
        )


def load_bed_regions(bed_path: str | Path, contig: str) -> List[Tuple[int, int]]:
    """Read half-open intervals for one contig from a BED file.

    BED is already 0-based half-open, so no conversion is applied. Both
    plain and gzipped files are accepted.

    Args:
        bed_path: Path to the BED file.
        contig: Contig to filter on. Matching ignores any ``chr`` prefix, so
            a ``21`` BED works against a ``chr21`` request and vice versa.

    Returns:
        Sorted, merged list of ``(start, end)`` intervals.
    """
    bed_path = Path(bed_path)
    wanted = contig[3:] if contig.startswith("chr") else contig

    opener = gzip.open if bed_path.suffix == ".gz" else open
    intervals: List[Tuple[int, int]] = []

    with opener(bed_path, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            name = fields[0]
            bare = name[3:] if name.startswith("chr") else name
            if bare != wanted:
                continue
            intervals.append((int(fields[1]), int(fields[2])))

    if not intervals:
        return []

    intervals.sort()
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


class GiabAlignmentProvider(VariantDataProvider):
    """Reads GIAB benchmark data on a contig into training windows.

    Expects:

    * **FASTA** — reference genome with a ``.fai`` index.
    * **BAM** — coordinate-sorted, indexed aligned reads.
    * **VCF** — GIAB benchmark calls, ideally bgzipped and tabix-indexed.
    * **BED** (optional but recommended) — GIAB high-confidence regions.

    Handles are opened lazily so constructing the provider is cheap and safe
    inside a DataLoader worker.
    """

    def __init__(
        self,
        fasta_path: str | Path,
        bam_path: str | Path,
        vcf_path: str | Path,
        seq_len: int = 64,
        contig: str = "chr21",
        stride: int | None = None,
        min_base_quality: int = 13,
        min_mapping_quality: int = 20,
        min_depth: int = 5,
        vaf_threshold: float = 0.15,
        indel_support_threshold: float = 0.25,
        high_confidence_bed: str | Path | None = None,
        max_windows: int | None = None,
        region: Tuple[int, int] | None = None,
    ) -> None:
        """Configure the provider without opening any files.

        Args:
            fasta_path: Indexed reference FASTA.
            bam_path: Coordinate-sorted, indexed BAM.
            vcf_path: Benchmark VCF.
            seq_len: Window width in alignment columns.
            contig: Contig to tile, e.g. ``"chr21"``.
            stride: Step between window starts, defaulting to ``seq_len``
                (non-overlapping). A stride below ``seq_len`` overlaps
                windows so variants near an edge are seen with full context
                at least once.
            min_base_quality: Phred floor below which a base is not counted
                toward the consensus.
            min_mapping_quality: MAPQ floor below which a read is skipped.
                Low-MAPQ reads in repeats are the main false-positive source
                on real chr21.
            min_depth: Loci with fewer passing reads are emitted as ``N``
                rather than guessed at.
            vaf_threshold: Minimum non-reference allele fraction required
                before the consensus departs from the reference base. The
                default suits a germline diploid sample, where a true
                heterozygous site sits near 0.5 and errors near 0.01.
            indel_support_threshold: Minimum fraction of covering reads that
                must show an indel before it enters the consensus. Set above
                ``vaf_threshold`` because indel alignments are noisier than
                substitutions.
            high_confidence_bed: GIAB high-confidence BED. When given, only
                windows fully inside it are kept.
            max_windows: Optional cap for quick diagnostic runs.
            region: Optional ``(start, end)`` restricting tiling to a slice
                of the contig, for testing without scanning 46 Mb.

        Raises:
            ValueError: If a numeric argument is outside its valid range.
        """
        super().__init__(seq_len)

        self.fasta_path = Path(fasta_path)
        self.bam_path = Path(bam_path)
        self.vcf_path = Path(vcf_path)
        self.contig = contig
        self.stride = stride if stride is not None else seq_len
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality
        self.min_depth = min_depth
        self.vaf_threshold = vaf_threshold
        self.indel_support_threshold = indel_support_threshold
        self.high_confidence_bed = Path(high_confidence_bed) if high_confidence_bed else None
        self.max_windows = max_windows
        self.region = region

        if self.stride <= 0:
            raise ValueError("stride must be positive")
        if not 0.0 <= vaf_threshold <= 1.0:
            raise ValueError("vaf_threshold must lie in [0, 1]")
        if not 0.0 <= indel_support_threshold <= 1.0:
            raise ValueError("indel_support_threshold must lie in [0, 1]")
        if min_depth < 0:
            raise ValueError("min_depth must be non-negative")

        self._fasta: "pysam.FastaFile | None" = None
        self._bam: "pysam.AlignmentFile | None" = None
        self._vcf: "pysam.VariantFile | None" = None
        self._windows: List[GenomicWindow] | None = None
        self._confident_regions: List[Tuple[int, int]] | None = None

        # Per-handle contig spellings, resolved on open().
        self._fasta_contig: str | None = None
        self._bam_contig: str | None = None
        self._vcf_contig: str | None = None

    # -- lifecycle ---------------------------------------------------------

    def open(self) -> None:
        """Open the FASTA, BAM and VCF handles and resolve contig naming.

        Raises:
            ImportError: If ``pysam`` is not installed.
            ValueError: If the contig cannot be resolved in any file. This is
                deliberately loud: a silent naming mismatch yields an empty
                dataset that trains to a degenerate all-``Normal`` model,
                which is far harder to diagnose than an exception here.
        """
        import pysam

        if self._fasta is None:
            self._fasta = pysam.FastaFile(str(self.fasta_path))
        if self._bam is None:
            self._bam = pysam.AlignmentFile(str(self.bam_path), "rb")
        if self._vcf is None:
            self._vcf = pysam.VariantFile(str(self.vcf_path))

        self._fasta_contig = self._normalise_contig(self.contig, self._fasta.references)
        self._bam_contig = self._normalise_contig(self.contig, self._bam.references)
        self._vcf_contig = self._normalise_contig(
            self.contig, tuple(self._vcf.header.contigs.keys())
        )

    def close(self) -> None:
        """Close any open handles."""
        for attribute in ("_fasta", "_bam", "_vcf"):
            handle = getattr(self, attribute)
            if handle is not None:
                handle.close()
                setattr(self, attribute, None)

    def __enter__(self) -> "GiabAlignmentProvider":
        """Open handles on context entry.

        Returns:
            This provider.
        """
        self.open()
        return self

    def __exit__(self, *exc_info: object) -> None:
        """Close handles on context exit."""
        self.close()

    def _ensure_open(self) -> None:
        """Open handles if they are not already open."""
        if self._fasta is None or self._bam is None or self._vcf is None:
            self.open()

    # -- coordinate and naming handling ------------------------------------

    @staticmethod
    def _normalise_contig(name: str, available: Sequence[str]) -> str:
        """Reconcile ``chr21`` against ``21`` naming for one file.

        GRCh38 analysis sets and GIAB releases disagree on the ``chr``
        prefix, and mixing them silently produces zero windows.

        Args:
            name: Contig name as supplied by the caller.
            available: Contig names present in the file.

        Returns:
            The spelling valid for that file.

        Raises:
            ValueError: If neither spelling is present.
        """
        if name in available:
            return name

        alternative = name[3:] if name.startswith("chr") else f"chr{name}"
        if alternative in available:
            return alternative

        preview = ", ".join(list(available)[:5])
        raise ValueError(
            f"contig {name!r} (and {alternative!r}) not found; file has: {preview}..."
        )

    def _contig_length(self) -> int:
        """Return the contig length as reported by the FASTA index.

        Returns:
            Length in bases.
        """
        self._ensure_open()
        assert self._fasta is not None and self._fasta_contig is not None
        index = self._fasta.references.index(self._fasta_contig)
        return int(self._fasta.lengths[index])

    def _is_confident(self, start: int, end: int) -> bool:
        """Check whether a window lies entirely in high-confidence regions.

        Args:
            start: 0-based inclusive start.
            end: 0-based exclusive end.

        Returns:
            ``True`` if no BED was supplied, or the window is fully covered.
        """
        if self.high_confidence_bed is None:
            return True

        if self._confident_regions is None:
            self._confident_regions = load_bed_regions(self.high_confidence_bed, self.contig)

        return any(
            region_start <= start and end <= region_end
            for region_start, region_end in self._confident_regions
        )

    def _build_windows(self) -> List[GenomicWindow]:
        """Tile the contig into windows, skipping uncallable regions.

        Windows whose reference slice is mostly ``N`` are dropped: the chr21
        acrocentric short arm is a large such block and carries no signal.

        Returns:
            The ordered list of windows to serve.
        """
        self._ensure_open()

        contig_start, contig_end = self.region or (0, self._contig_length())
        windows: List[GenomicWindow] = []

        for start in range(contig_start, contig_end - self.seq_len + 1, self.stride):
            end = start + self.seq_len

            if not self._is_confident(start, end):
                continue

            reference = self._fetch_reference(GenomicWindow(self.contig, start, end))
            if reference.count("N") > len(reference) // 2:
                continue

            windows.append(GenomicWindow(self.contig, start, end))
            if self.max_windows is not None and len(windows) >= self.max_windows:
                break

        return windows

    # -- per-window extraction --------------------------------------------

    def _fetch_reference(self, window: GenomicWindow) -> str:
        """Read the reference slice for a window.

        Uppercased because soft-masked repeats are stored lowercase and
        would otherwise tokenise to ``N``.

        Args:
            window: Interval to read, 0-based half-open.

        Returns:
            An uppercase reference string of length ``len(window)``, right
            padded with ``N`` if the contig ends early.
        """
        self._ensure_open()
        assert self._fasta is not None and self._fasta_contig is not None

        sequence = self._fasta.fetch(self._fasta_contig, window.start, window.end).upper()
        if len(sequence) < len(window):
            sequence += "N" * (len(window) - len(sequence))
        return sequence

    def _pileup_locus(
        self, column: "pysam.PileupColumn", reference_base: str
    ) -> LocusPileup:
        """Summarise one pileup column into a consensus call.

        The consensus rule is deliberately conservative: the column stays at
        the reference base unless a non-reference allele clears
        ``vaf_threshold``. On a germline diploid sample a true heterozygous
        site sits near 0.5 and sequencing errors near 0.01, so a threshold of
        0.15 separates them with wide margin while still admitting sites with
        skewed allelic balance.

        Args:
            column: A pysam pileup column.
            reference_base: The reference base at this locus.

        Returns:
            The consensus summary for this locus.
        """
        base_counts: Counter[str] = Counter()
        quality_by_base: Dict[str, List[int]] = {}
        insertion_sequences: List[str] = []
        depth = 0

        for read in column.pileups:
            alignment = read.alignment

            if alignment.mapping_quality < self.min_mapping_quality:
                continue
            if alignment.is_duplicate or alignment.is_qcfail or alignment.is_unmapped:
                continue
            if alignment.is_secondary or alignment.is_supplementary:
                continue

            depth += 1

            if read.is_del or read.is_refskip:
                base_counts[GAP_TOKEN] += 1
                continue

            position = read.query_position
            if position is None:
                continue

            quality = alignment.query_qualities[position] if alignment.query_qualities else 0
            if quality < self.min_base_quality:
                continue

            base = alignment.query_sequence[position].upper()
            if base not in NUCLEOTIDES:
                continue

            base_counts[base] += 1
            quality_by_base.setdefault(base, []).append(int(quality))

            # A positive indel on a pileup read means the bases immediately
            # following this locus are inserted relative to the reference.
            if read.indel > 0:
                inserted = alignment.query_sequence[
                    position + 1 : position + 1 + read.indel
                ].upper()
                if inserted:
                    insertion_sequences.append(inserted)

        called_bases = sum(count for base, count in base_counts.items() if base in NUCLEOTIDES)
        non_reference = sum(
            count
            for base, count in base_counts.items()
            if base in NUCLEOTIDES and base != reference_base
        )
        vaf = non_reference / called_bases if called_bases else 0.0

        consensus_base = reference_base
        if depth < self.min_depth:
            # Not enough evidence to call anything; N marks "read but
            # untrustworthy", which is distinct from the '-' alignment gap.
            consensus_base = "N"
        elif vaf >= self.vaf_threshold:
            alternatives = [
                (count, base)
                for base, count in base_counts.items()
                if base in NUCLEOTIDES and base != reference_base
            ]
            if alternatives:
                consensus_base = max(alternatives)[1]

        supporting = quality_by_base.get(consensus_base, [])
        quality = sum(supporting) / len(supporting) if supporting else 0.0

        return LocusPileup(
            depth=depth,
            base_counts=dict(base_counts),
            consensus_base=consensus_base,
            reference_base=reference_base,
            vaf=vaf,
            quality=quality,
            insertion_sequences=insertion_sequences,
        )

    def _collect_pileups(self, window: GenomicWindow, reference: str) -> List[LocusPileup]:
        """Build a consensus summary for every locus in the window.

        Loci with no coverage at all still receive an entry, so the returned
        list is always exactly ``len(window)`` long and indexable by offset.

        Args:
            window: Interval to summarise.
            reference: Reference sequence for the window.

        Returns:
            One :class:`LocusPileup` per reference position.
        """
        self._ensure_open()
        assert self._bam is not None and self._bam_contig is not None

        empty = [
            LocusPileup(
                depth=0,
                base_counts={},
                consensus_base="N",
                reference_base=reference[offset],
                vaf=0.0,
                quality=0.0,
                insertion_sequences=[],
            )
            for offset in range(len(window))
        ]

        # truncate keeps pysam from returning columns outside the requested
        # span; min_base_quality is applied per read below rather than here
        # so that depth still counts low-quality reads.
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
            empty[offset] = self._pileup_locus(column, reference[offset])

        return empty

    def _fetch_labels(self, window: GenomicWindow) -> List[int]:
        """Derive per-locus ground-truth labels from the benchmark VCF.

        Records are mapped by comparing ref and alt lengths:

        * equal length 1 → :data:`~.config.LABEL_SNP`
        * ``len(alt) > len(ref)`` → :data:`~.config.LABEL_INSERTION`
        * ``len(alt) < len(ref)`` → :data:`~.config.LABEL_DELETION`, spanning
          the deleted reference bases
        * equal length above 1 → an MNP, decomposed into per-base SNPs

        Two decisions worth naming rather than defaulting: coordinates come
        from ``record.start``, which pysam has already converted to 0-based
        (``record.pos`` is the raw 1-based value and is never used here); and
        records whose genotype is homozygous reference are skipped, since a
        site present in the VCF but called ``0/0`` is not a variant in this
        sample.

        Args:
            window: Interval to label.

        Returns:
            A list of ``len(window)`` class ids.
        """
        self._ensure_open()
        assert self._vcf is not None and self._vcf_contig is not None

        labels = [LABEL_NORMAL] * len(window)

        for record in self._vcf.fetch(self._vcf_contig, window.start, window.end):
            if record.alts is None:
                continue
            if self._is_homozygous_reference(record):
                continue

            reference_allele = record.ref or ""
            for alternate in record.alts:
                if alternate is None or alternate.startswith("<"):
                    continue  # symbolic allele, e.g. <DEL>; not handled here

                label, span = self._classify_allele(reference_allele, alternate)
                if label is None:
                    continue

                for offset in range(span):
                    index = record.start + offset - window.start
                    if 0 <= index < len(window):
                        labels[index] = label

        return labels

    @staticmethod
    def _is_homozygous_reference(record: "pysam.VariantRecord") -> bool:
        """Check whether every sample in the record is called ``0/0``.

        Args:
            record: A VCF record.

        Returns:
            ``True`` if all genotypes are homozygous reference. Records with
            no genotype information return ``False``, so a sites-only VCF is
            treated as carrying real variants.
        """
        if not record.samples:
            return False

        for sample in record.samples.values():
            genotype = sample.get("GT")
            if genotype is None:
                return False
            if any(allele not in (0, None) for allele in genotype):
                return False
        return True

    @staticmethod
    def _classify_allele(reference: str, alternate: str) -> Tuple[int | None, int]:
        """Map one ref/alt pair to a class id and the span it covers.

        Args:
            reference: The reference allele string.
            alternate: The alternate allele string.

        Returns:
            Tuple of ``(label, span)``. ``label`` is ``None`` when the allele
            should be ignored. ``span`` is the number of reference positions
            the event covers, starting at the record's own position.
        """
        if not reference or not alternate:
            return None, 0

        if len(reference) == len(alternate):
            if len(reference) == 1:
                return LABEL_SNP, 1
            # MNP: decomposed into per-base substitutions across its span.
            return LABEL_SNP, len(reference)

        if len(alternate) > len(reference):
            # VCF anchors an insertion on the preceding reference base; the
            # inserted sequence itself has no reference coordinate, so the
            # event is marked on the anchor.
            return LABEL_INSERTION, 1

        # Deletion: the anchor base remains, the following bases are deleted.
        return LABEL_DELETION, max(1, len(reference) - len(alternate))

    def _fetch_observed(
        self, window: GenomicWindow, reference: str
    ) -> Tuple[str, str, List[int], List[float], List[float], List[float]]:
        """Build gap-padded observed and reference columns from the pileup.

        Insertions add columns that carry ``-`` in the reference stream, so
        the two streams stay in 1-to-1 correspondence and no downstream frame
        shift occurs. Deletions carry ``-`` in the observed stream instead.
        Labels are carried through the same expansion so they never drift out
        of alignment with the columns they describe.

        Args:
            window: Interval being built.
            reference: Reference sequence for the window.

        Returns:
            Tuple of ``(reference_columns, observed_columns, column_offsets,
            qualities, depths, vafs)``. ``column_offsets`` maps each emitted
            column back to its 0-based offset within the window, which is
            what lets labels be expanded alongside.
        """
        pileups = self._collect_pileups(window, reference)

        reference_columns: List[str] = []
        observed_columns: List[str] = []
        column_offsets: List[int] = []
        qualities: List[float] = []
        depths: List[float] = []
        vafs: List[float] = []

        for offset, locus in enumerate(pileups):
            reference_columns.append(locus.reference_base)
            depths.append(float(locus.depth))
            vafs.append(locus.vaf)
            column_offsets.append(offset)

            deletion_support = (
                locus.base_counts.get(GAP_TOKEN, 0) / locus.depth if locus.depth else 0.0
            )
            if deletion_support >= self.indel_support_threshold:
                observed_columns.append(GAP_TOKEN)
                qualities.append(0.0)
            else:
                observed_columns.append(locus.consensus_base)
                qualities.append(locus.quality)

            # Insertion columns follow the anchor locus and are gap-padded in
            # the reference stream.
            if locus.insertion_support >= self.indel_support_threshold:
                for base in locus.consensus_insertion():
                    reference_columns.append(GAP_TOKEN)
                    observed_columns.append(base if base in NUCLEOTIDES else "N")
                    column_offsets.append(offset)
                    qualities.append(locus.quality)
                    depths.append(float(locus.depth))
                    vafs.append(locus.insertion_support)

        return (
            "".join(reference_columns),
            "".join(observed_columns),
            column_offsets,
            qualities,
            depths,
            vafs,
        )

    # -- provider contract -------------------------------------------------

    def __len__(self) -> int:
        """Return the number of tiled windows.

        Returns:
            Window count.
        """
        if self._windows is None:
            self._windows = self._build_windows()
        return len(self._windows)

    def __getitem__(self, index: int) -> WindowTensors:
        """Assemble one window into model-ready tensors.

        Args:
            index: Positional index into the tiled windows.

        Returns:
            A fully populated :class:`WindowTensors`, including VAF.

        Raises:
            IndexError: If ``index`` is out of range.
        """
        if self._windows is None:
            self._windows = self._build_windows()
        if not 0 <= index < len(self._windows):
            raise IndexError(f"index {index} out of range for {len(self._windows)} windows")

        window = self._windows[index]
        reference = self._fetch_reference(window)
        locus_labels = self._fetch_labels(window)

        (
            reference_columns,
            observed_columns,
            column_offsets,
            qualities,
            depths,
            vafs,
        ) = self._fetch_observed(window, reference)

        # Expand locus labels across the emitted columns. An inserted column
        # inherits its anchor's label, which is Insertion whenever the VCF
        # agreed there was one.
        labels = [locus_labels[offset] for offset in column_offsets]

        reference_columns, observed_columns, labels, qualities, depths, vafs = self._fit_to_length(
            reference_columns, observed_columns, labels, qualities, depths, vafs
        )

        return self._validate_widths(
            WindowTensors(
                reference_ids=self.tokenizer.encode(reference_columns),
                input_ids=self.tokenizer.encode(observed_columns),
                labels=torch.tensor(labels, dtype=torch.long),
                base_quality=torch.tensor(qualities, dtype=torch.float),
                depth=torch.tensor(depths, dtype=torch.float),
                vaf=torch.tensor(vafs, dtype=torch.float),
                window=window,
            )
        )

    def _fit_to_length(
        self,
        reference_columns: str,
        observed_columns: str,
        labels: List[int],
        qualities: List[float],
        depths: List[float],
        vafs: List[float],
    ) -> Tuple[str, str, List[int], List[float], List[float], List[float]]:
        """Truncate or pad every channel together to ``seq_len`` columns.

        Insertions add columns, so a window can overrun; padding covers the
        contig-end case. Padded columns use ``N`` and ``Normal`` so they are
        inert rather than spurious variants.

        Args:
            reference_columns: Reference column string.
            observed_columns: Observed column string.
            labels: Per-column class ids.
            qualities: Per-column Phred scores.
            depths: Per-column read depths.
            vafs: Per-column variant allele fractions.

        Returns:
            The six channels, each of exactly ``seq_len`` elements.
        """
        length = len(reference_columns)

        if length > self.seq_len:
            end = self.seq_len
            return (
                reference_columns[:end],
                observed_columns[:end],
                labels[:end],
                qualities[:end],
                depths[:end],
                vafs[:end],
            )

        if length < self.seq_len:
            pad = self.seq_len - length
            return (
                reference_columns + "N" * pad,
                observed_columns + "N" * pad,
                labels + [LABEL_NORMAL] * pad,
                qualities + [0.0] * pad,
                depths + [0.0] * pad,
                vafs + [0.0] * pad,
            )

        return reference_columns, observed_columns, labels, qualities, depths, vafs


class ProviderDataset(torch.utils.data.Dataset[VariantSample]):
    """Adapts any :class:`VariantDataProvider` to a torch ``Dataset``.

    Lets a real-data provider feed the existing training loop with no change
    to the loop itself.
    """

    def __init__(self, provider: VariantDataProvider) -> None:
        """Wrap a provider.

        Args:
            provider: The provider to serve samples from.
        """
        self.provider = provider

    def __len__(self) -> int:
        """Return the number of windows."""
        return len(self.provider)

    def __getitem__(self, index: int) -> VariantSample:
        """Return the window at ``index`` in five-tensor form.

        Args:
            index: Positional index.

        Returns:
            A :class:`~.dataset.VariantSample`.
        """
        return self.provider[index].to_sample()


def summarise_label_distribution(labels: Sequence[int]) -> Dict[str, int]:
    """Count label occurrences by class name.

    A near-total absence of a class on real data usually means a coordinate
    or contig-naming bug rather than a genuinely quiet region.

    Args:
        labels: Flat sequence of class ids.

    Returns:
        Mapping of class name to count.
    """
    names = {
        LABEL_NORMAL: "Normal",
        LABEL_SNP: "SNP",
        LABEL_INSERTION: "Insertion",
        LABEL_DELETION: "Deletion",
    }
    counts = {name: 0 for name in names.values()}
    for label in labels:
        counts[names[int(label)]] += 1
    return counts


__all__ = [
    "GenomicWindow",
    "LocusPileup",
    "WindowTensors",
    "VariantDataProvider",
    "SyntheticProvider",
    "GiabAlignmentProvider",
    "ProviderDataset",
    "load_bed_regions",
    "summarise_label_distribution",
]
