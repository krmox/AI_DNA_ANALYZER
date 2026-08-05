"""Data provider interfaces for synthetic and real sequencing input.

The point of this module is a single contract — :class:`VariantDataProvider`
— that both the synthetic simulator and a real GIAB BAM/VCF/FASTA pipeline
satisfy. Anything downstream (model, loss, training loop) consumes windows
through that contract and does not care which side produced them.

Status: :class:`GiabAlignmentProvider` is a **scaffold**. Its structure,
coordinate handling and CIGAR semantics are specified in full, but the
``pysam`` calls that read actual files raise
:class:`NotImplementedError` where real I/O belongs. This is deliberate —
implementing them against no real data would produce code that type-checks
and is silently wrong. Every such site carries a docstring describing the
exact expected behaviour.

Dependency note: ``pysam`` is imported lazily inside the provider so the
synthetic path and all existing tests keep running on a bare torch install.

**The signal that is currently missing.** The synthetic generator cannot
distinguish a sequencing error from a true SNP, because both reduce to one
mismatched base (see :meth:`~.dataset.SyntheticVariantDataset._apply_noise`).
Real BAM input carries the two signals that break that tie, and both are
surfaced by this interface:

* ``base_quality`` — Phred scores per base. Sequencing errors concentrate
  at low Q; a Q40 mismatch is far likelier to be real than a Q15 one.
* ``depth`` and ``allele_fraction`` — a true variant recurs across the
  reads covering a locus, while an error appears in one read.

:class:`WindowTensors` therefore carries optional fields for both. They are
``None`` for synthetic data and populated from BAM. Feeding them to the
model requires additional embedding or projection channels, which is a
model change and is intentionally *not* made here.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Iterator, List, Sequence

import torch

from .config import LABEL_DELETION, LABEL_INSERTION, LABEL_NORMAL, LABEL_SNP
from .dataset import DNATokenizer, SyntheticVariantDataset, VariantSample

if TYPE_CHECKING:  # pragma: no cover - typing only
    import pysam


@dataclass(frozen=True)
class GenomicWindow:
    """A half-open reference interval on one contig.

    Attributes:
        contig: Reference contig name, e.g. ``"chr21"``. Note that GIAB
            releases differ in whether contigs carry the ``chr`` prefix; see
            :meth:`GiabAlignmentProvider._normalise_contig`.
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
class WindowTensors:
    """Model-ready tensors for one window, plus optional real-data signals.

    The three required fields match :class:`~.dataset.VariantSample` exactly,
    so a provider's output drops into the existing training loop unchanged.

    Attributes:
        reference_ids: Long tensor ``[seq_len]`` of reference token ids.
        input_ids: Long tensor ``[seq_len]`` of observed token ids.
        labels: Long tensor ``[seq_len]`` of per-token class ids.
        base_quality: Optional float tensor ``[seq_len]`` of Phred scores,
            normalised to roughly ``[0, 1]``. ``None`` for synthetic data.
        depth: Optional float tensor ``[seq_len]`` of read depth per
            position. ``None`` for synthetic data.
        window: The genomic interval this came from, when known. Carried so
            predictions can be written back to VCF coordinates.
    """

    reference_ids: torch.Tensor
    input_ids: torch.Tensor
    labels: torch.Tensor
    base_quality: torch.Tensor | None = None
    depth: torch.Tensor | None = None
    window: GenomicWindow | None = None

    def to_sample(self) -> VariantSample:
        """Reduce to the three-tensor form the current model consumes.

        Zeros stand in for absent quality or depth, which the model reads
        as "no confidence information" — the same value gap columns carry.

        Returns:
            A :class:`~.dataset.VariantSample`.
        """
        zeros = torch.zeros_like(self.labels, dtype=torch.float)
        return VariantSample(
            reference_ids=self.reference_ids,
            input_ids=self.input_ids,
            labels=self.labels,
            base_quality=self.base_quality if self.base_quality is not None else zeros,
            depth=self.depth if self.depth is not None else zeros,
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
        for name in ("reference_ids", "input_ids", "labels"):
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
            A fully populated :class:`WindowTensors`. The simulator now
            models Phred quality and coverage explicitly, so these are no
            longer ``None`` on the synthetic path.
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


class GiabAlignmentProvider(VariantDataProvider):
    """Scaffold for reading GIAB chr21 truth sets into training windows.

    Intended inputs:

    * **FASTA** — reference genome, for the reference channel.
    * **BAM** — aligned reads, for the observed channel plus base qualities
      and depth. Must be coordinate-sorted and indexed (``.bai``).
    * **VCF** — GIAB benchmark calls, for ground-truth labels. Should be
      intersected with the GIAB high-confidence BED before use; regions
      outside it are not reliably labelled and training on them injects
      false negatives.

    Every method that touches a file raises :class:`NotImplementedError`.
    The docstrings specify what each must do, including the CIGAR handling
    that is the genuinely delicate part.

    Warning:
        Do not treat metrics from the synthetic simulator as predictive of
        performance here. Real chr21 differs in ways that matter: reference
        base composition is not uniform, roughly half the sequence is
        repetitive, variants cluster rather than occurring independently,
        and chr21 carries a large acrocentric region that is effectively
        uncallable. Class balance also shifts by orders of magnitude — real
        variant density is on the order of 1 in 1000 bases, against the
        1 in 20 used in simulation.
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
        high_confidence_bed: str | Path | None = None,
    ) -> None:
        """Configure the provider without opening any files.

        File handles are opened lazily so that constructing a provider is
        cheap and safe in a DataLoader worker.

        Args:
            fasta_path: Indexed reference FASTA.
            bam_path: Coordinate-sorted, indexed BAM of aligned reads.
            vcf_path: GIAB benchmark VCF, ideally bgzipped and tabix-indexed.
            seq_len: Window width in tokens.
            contig: Contig to tile, e.g. ``"chr21"``.
            stride: Step between window starts. Defaults to ``seq_len``
                (non-overlapping). A stride below ``seq_len`` overlaps
                windows, which helps variants near a window edge get seen
                with full context at least once.
            min_base_quality: Phred floor below which a base is masked to
                ``N`` rather than trusted.
            min_mapping_quality: MAPQ floor below which a read is skipped
                entirely. Low-MAPQ reads in repetitive regions are the main
                source of false positives on real chr21.
            high_confidence_bed: Optional GIAB high-confidence BED. When
                given, windows not fully contained in it are skipped.

        Raises:
            ValueError: If ``seq_len`` or ``stride`` is not positive.
        """
        super().__init__(seq_len)

        self.fasta_path = Path(fasta_path)
        self.bam_path = Path(bam_path)
        self.vcf_path = Path(vcf_path)
        self.contig = contig
        self.stride = stride if stride is not None else seq_len
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality
        self.high_confidence_bed = Path(high_confidence_bed) if high_confidence_bed else None

        if self.stride <= 0:
            raise ValueError("stride must be positive")

        self._fasta: "pysam.FastaFile | None" = None
        self._bam: "pysam.AlignmentFile | None" = None
        self._vcf: "pysam.VariantFile | None" = None
        self._windows: List[GenomicWindow] | None = None

    # -- lifecycle ---------------------------------------------------------

    def open(self) -> None:
        """Open the FASTA, BAM and VCF handles.

        Implementation notes:
            * ``self._fasta = pysam.FastaFile(str(self.fasta_path))``
            * ``self._bam = pysam.AlignmentFile(str(self.bam_path), "rb")``
            * ``self._vcf = pysam.VariantFile(str(self.vcf_path))``
            * Verify ``self.contig`` resolves in all three, applying
              :meth:`_normalise_contig`, and fail loudly if it does not.
              A silent contig-naming mismatch yields an empty dataset that
              trains to a degenerate all-``Normal`` model, which is far
              harder to diagnose than an exception here.

        Raises:
            NotImplementedError: Always, in this scaffold.
        """
        raise NotImplementedError("Open pysam handles for FASTA, BAM and VCF.")

    def close(self) -> None:
        """Close any open handles.

        Raises:
            NotImplementedError: Always, in this scaffold.
        """
        raise NotImplementedError("Close any open pysam handles.")

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

    # -- coordinate handling ----------------------------------------------

    def _normalise_contig(self, name: str) -> str:
        """Reconcile ``chr21`` against ``21`` naming between files.

        GRCh38 analysis sets and GIAB releases disagree on the ``chr``
        prefix, and mixing them silently produces zero windows.

        Implementation notes:
            Inspect the contig list of each open handle and return whichever
            spelling that handle uses.

        Args:
            name: Contig name as supplied by the caller.

        Returns:
            The spelling valid for the currently open files.

        Raises:
            NotImplementedError: Always, in this scaffold.
        """
        raise NotImplementedError("Resolve contig naming against open file headers.")

    def _build_windows(self) -> List[GenomicWindow]:
        """Tile the contig into windows, skipping uncallable regions.

        Implementation notes:
            * Step from 0 to the contig length by ``self.stride``.
            * Skip windows whose reference slice is mostly ``N`` — the chr21
              acrocentric arm is a large such block and carries no signal.
            * When ``high_confidence_bed`` is set, keep only windows fully
              inside it.
            * Optionally skip zero-coverage windows by consulting the BAM
              index, which is much cheaper than discovering it per window.

        Returns:
            The ordered list of windows to serve.

        Raises:
            NotImplementedError: Always, in this scaffold.
        """
        raise NotImplementedError("Tile the contig into callable windows.")

    # -- per-window extraction --------------------------------------------

    def _fetch_reference(self, window: GenomicWindow) -> str:
        """Read the reference slice for a window.

        Implementation notes:
            ``self._fasta.fetch(window.contig, window.start, window.end)``,
            uppercased — soft-masked repeats come back lowercase and would
            otherwise tokenise to ``N``.

        Args:
            window: Interval to read.

        Returns:
            An uppercase reference string of length ``len(window)``.

        Raises:
            NotImplementedError: Always, in this scaffold.
        """
        raise NotImplementedError("Fetch the reference slice via pysam.FastaFile.")

    def _fetch_observed(self, window: GenomicWindow) -> tuple[str, List[float], List[float]]:
        """Build the observed channel from aligned reads.

        This is the part that most repays care, because it is where the
        synthetic simulator's shortcuts are repaid. Reads carry CIGAR
        strings, and translating them into a fixed-width per-position view
        is the actual alignment work:

        * ``M``/``=``/``X`` — consume both read and reference; the read base
          occupies the reference position.
        * ``I`` — consumes read only. The inserted bases have no reference
          coordinate. Choose a representation and hold it consistently with
          whatever the label builder does; anchoring the insertion on the
          preceding reference position is the usual convention and matches
          how VCF represents it.
        * ``D`` — consumes reference only. Emit ``N`` at the deleted
          positions, matching the synthetic gap convention.
        * ``S``/``H`` — soft and hard clips; skip, and do not let soft-clipped
          bases leak into the window.
        * ``N`` — reference skip, relevant for spliced RNA alignments.

        With multiple reads covering a position, decide explicitly between
        taking the highest-MAPQ read and taking a consensus across reads.
        Consensus is closer to what a real caller does and is what makes
        depth and allele fraction meaningful; single-read selection throws
        away the very signal that distinguishes an error from a variant.

        Apply ``min_base_quality`` by masking failing bases to ``N``, and
        skip reads below ``min_mapping_quality`` outright.

        Args:
            window: Interval to build.

        Returns:
            Tuple of ``(observed, base_qualities, depths)``, each of length
            ``len(window)``. Qualities are normalised Phred scores; depths
            are raw read counts.

        Raises:
            NotImplementedError: Always, in this scaffold.
        """
        raise NotImplementedError("Build the observed channel from BAM via CIGAR walk.")

    def _fetch_labels(self, window: GenomicWindow) -> List[int]:
        """Derive ground-truth labels from the benchmark VCF.

        Implementation notes:
            Fetch records overlapping the window and map each to this
            package's classes by comparing ``ref`` and ``alt`` lengths:

            * equal lengths, length 1 → :data:`~.config.LABEL_SNP`
            * ``len(alt) > len(ref)`` → :data:`~.config.LABEL_INSERTION`,
              spanning the inserted zone
            * ``len(alt) < len(ref)`` → :data:`~.config.LABEL_DELETION`,
              spanning the deleted positions
            * equal lengths above 1 → an MNP; the current four-class scheme
              has no slot for it. Decide explicitly whether to decompose it
              into per-base SNPs or drop it, and record the choice, because
              it shifts the SNP class balance either way.

            Two further cases need a deliberate decision rather than a
            default: VCF is 1-based while :class:`GenomicWindow` is 0-based,
            so subtract one when converting; and genotype matters, since a
            record present as ``0/0`` is not a variant in this sample.
            Everything not covered by a record is
            :data:`~.config.LABEL_NORMAL`.

        Args:
            window: Interval to label.

        Returns:
            A list of ``len(window)`` class ids.

        Raises:
            NotImplementedError: Always, in this scaffold.
        """
        raise NotImplementedError("Derive per-position labels from the benchmark VCF.")

    # -- provider contract -------------------------------------------------

    def __len__(self) -> int:
        """Return the number of tiled windows.

        Returns:
            Window count.

        Raises:
            NotImplementedError: Until :meth:`_build_windows` is implemented.
        """
        if self._windows is None:
            self._windows = self._build_windows()
        return len(self._windows)

    def __getitem__(self, index: int) -> WindowTensors:
        """Assemble one window into model-ready tensors.

        Args:
            index: Positional index into the tiled windows.

        Returns:
            A fully populated :class:`WindowTensors`, including quality and
            depth.

        Raises:
            IndexError: If ``index`` is out of range.
            NotImplementedError: Until the extraction methods are implemented.
        """
        if self._windows is None:
            self._windows = self._build_windows()
        if not 0 <= index < len(self._windows):
            raise IndexError(f"index {index} out of range for {len(self._windows)} windows")

        window = self._windows[index]
        reference = self._fetch_reference(window)
        observed, qualities, depths = self._fetch_observed(window)
        labels = self._fetch_labels(window)

        return self._validate_widths(
            WindowTensors(
                reference_ids=self.tokenizer.encode(reference),
                input_ids=self.tokenizer.encode(observed),
                labels=torch.tensor(labels, dtype=torch.long),
                base_quality=torch.tensor(qualities, dtype=torch.float),
                depth=torch.tensor(depths, dtype=torch.float),
                window=window,
            )
        )


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
        """Return the window at ``index`` in three-tensor form.

        Args:
            index: Positional index.

        Returns:
            A :class:`~.dataset.VariantSample`.
        """
        return self.provider[index].to_sample()


def summarise_label_distribution(labels: Sequence[int]) -> dict[str, int]:
    """Count label occurrences by class name.

    Useful as a first sanity check on real data, where a near-total absence
    of a class usually means a coordinate or contig-naming bug rather than a
    genuinely quiet region.

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
    "WindowTensors",
    "VariantDataProvider",
    "SyntheticProvider",
    "GiabAlignmentProvider",
    "ProviderDataset",
    "summarise_label_distribution",
]
