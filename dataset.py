"""Gap-padded read simulation with quality and depth channels.

Two upgrades over the frame-shifted generator, addressing the two ceilings
found in stress testing.

**CIGAR-style padded alignment.** Both streams are now emitted at the same
coordinates, with an explicit ``-`` gap token filling the stream that has no
base at a given column:

======================================================================
Kind       Reference column(s)   Observed column(s)   Labels
======================================================================
Normal     ``ref``               ``ref``              ``[0]``
SNP        ``ref``               ``base != ref``      ``[1]``
Insertion  ``ref`` + ``-`` * L   ``ref`` + L bases    ``[0]`` + ``[2]`` * L
Deletion   ``ref``               ``-``                ``[3]``
======================================================================

This eliminates the downstream frame shift entirely: ``reference_ids[i]``
and ``input_ids[i]`` always describe the same alignment column.

Be clear about what this does and does not buy. It does not solve
alignment — it assumes alignment, exactly as a pileup-based caller does
downstream of BWA or minimap2. A consequence is that indels become
near-trivially detectable, since ``reference_ids[i] == GAP_ID`` fully
determines an insertion column and ``input_ids[i] == GAP_ID`` a deletion
column. Indel F1 from this generator therefore measures whether the model
can read a marker, not whether it can find an indel. The genuinely hard
problem that remains is separating sequencing noise from true SNPs.

**Quality and depth channels.** Each column carries a Phred score and a
read depth. Errors are drawn from a low-Phred distribution and correct
bases from a high-Phred one, which is what makes noise and SNPs separable
at all — see :meth:`SyntheticVariantDataset._draw_quality`. That
correlation is a modelling assumption planted by this simulator; real
Q-score distributions overlap considerably more.
"""

from __future__ import annotations

import random
from typing import Dict, List, Sequence, Tuple, TypedDict

import torch
from torch.utils.data import DataLoader, Dataset

from config import (
    GAP_ID,
    GAP_TOKEN,
    INV_VOCAB,
    LABEL_DELETION,
    LABEL_INSERTION,
    LABEL_NORMAL,
    LABEL_SNP,
    MUTATION_KINDS,
    NUCLEOTIDES,
    VOCAB,
    DataConfig,
)


class VariantSample(TypedDict):
    """One aligned training example.

    Attributes:
        reference_ids: Long tensor ``[seq_len]`` of reference column tokens,
            gap-padded opposite insertions.
        input_ids: Long tensor ``[seq_len]`` of observed column tokens,
            gap-padded opposite deletions.
        labels: Long tensor ``[seq_len]`` of per-column class ids in ``0..3``.
            Sequencing errors are deliberately absent from this array.
        base_quality: Float tensor ``[seq_len]`` of Phred scores. Gap columns
            carry 0.
        depth: Float tensor ``[seq_len]`` of read depth. Gap columns carry 0.
        vaf: Float tensor ``[seq_len]`` of variant allele fractions. The
            synthetic simulator emits **all zeros** here: it models a single
            read, and one read cannot express an allele fraction. Deriving a
            VAF from the labels instead would hand the model a copy of its
            own target, so ``use_vaf`` is inert on synthetic data by design
            and only carries signal from a real pileup.
    """

    reference_ids: torch.Tensor
    input_ids: torch.Tensor
    labels: torch.Tensor
    base_quality: torch.Tensor
    depth: torch.Tensor
    vaf: torch.Tensor


class DNATokenizer:
    """Character-level tokenizer over ``{N, A, C, G, T, -}``.

    ``N`` (0) is padding and unknown base; ``-`` (5) is the alignment gap.
    The two are kept distinct because they mean different things: ``N`` is a
    base that was read but is untrustworthy, ``-`` is the absence of any base
    at that alignment column.
    """

    def __init__(self) -> None:
        """Initialise the tokenizer from the package-level vocabulary."""
        self.vocab: Dict[str, int] = dict(VOCAB)
        self.inv_vocab: Dict[int, str] = dict(INV_VOCAB)
        self.vocab_size: int = len(self.vocab)

    def encode(self, sequence: str) -> torch.Tensor:
        """Encode a nucleotide string into token ids.

        Characters outside the alphabet fall back to ``N`` rather than
        raising, so malformed upstream reads degrade instead of crashing.

        Args:
            sequence: String over the alphabet, e.g. ``"ACGT-N"``.

        Returns:
            Long tensor of shape ``[len(sequence)]``.
        """
        return torch.tensor(
            [self.vocab.get(char, self.vocab["N"]) for char in sequence],
            dtype=torch.long,
        )

    def decode(self, token_ids: torch.Tensor | Sequence[int]) -> str:
        """Decode token ids back into a string.

        Args:
            token_ids: 1-D tensor or sequence of integer token ids.

        Returns:
            The decoded string. Unknown ids render as ``N``.
        """
        ids = token_ids.tolist() if isinstance(token_ids, torch.Tensor) else list(token_ids)
        return "".join(self.inv_vocab.get(int(index), "N") for index in ids)


class SyntheticVariantDataset(Dataset[VariantSample]):
    """Simulates gap-padded alignment columns with quality and depth.

    Seeding is *per index*, not per stream. A shared ``random.Random(seed)``
    is reproducible across runs but advances on every access, so a second
    pass yields different data — which silently makes each epoch validate on
    a different set. With a seed supplied, each index derives its own
    generator, so ``dataset[i]`` is stable across passes and shuffling.

    Note:
        With ``seed=None`` the generator streams fresh samples on every
        access, useful as augmentation for training. Never use an unseeded
        dataset for validation.
    """

    #: Multiplier mixing the base seed with the sample index. A large odd
    #: constant keeps neighbouring indices from producing correlated streams.
    _INDEX_SEED_STRIDE: int = 1_000_003

    def __init__(
        self,
        num_samples: int,
        seq_len: int,
        mutation_rate: float,
        seed: int | None = None,
        noise_rate: float = 0.0,
        max_insertion_len: int = 6,
        high_quality_mean: float = 36.0,
        low_quality_mean: float = 12.0,
        quality_spread: float = 5.0,
        max_phred: int = 60,
        mean_depth: float = 30.0,
        min_depth: int = 1,
    ) -> None:
        """Initialise the generator.

        Args:
            num_samples: Nominal dataset length reported by ``__len__``.
            seq_len: Fixed number of alignment columns per sample.
            mutation_rate: Per-reference-position probability of a variant.
            seed: Seed for per-index generators, or ``None`` to stream.
            noise_rate: Per-base probability of a sequencing substitution.
                Never labelled.
            max_insertion_len: Inclusive upper bound of the insertion length
                draw.
            high_quality_mean: Mean Phred for a correctly read base.
            low_quality_mean: Mean Phred for a sequencing error.
            quality_spread: Standard deviation of the Phred draw.
            max_phred: Upper clamp on Phred scores.
            mean_depth: Poisson mean for per-column read depth.
            min_depth: Floor applied to the depth draw.

        Raises:
            ValueError: If any argument is outside its permitted range.
        """
        if num_samples <= 0:
            raise ValueError("num_samples must be positive")
        if seq_len <= 0:
            raise ValueError("seq_len must be positive")
        if not 0.0 <= mutation_rate <= 1.0:
            raise ValueError("mutation_rate must lie in [0, 1]")
        if not 0.0 <= noise_rate <= 1.0:
            raise ValueError("noise_rate must lie in [0, 1]")
        if max_insertion_len < 1:
            raise ValueError("max_insertion_len must be at least 1")
        if max_phred < 1:
            raise ValueError("max_phred must be at least 1")
        if quality_spread < 0.0:
            raise ValueError("quality_spread must be non-negative")
        if mean_depth <= 0.0:
            raise ValueError("mean_depth must be positive")
        if min_depth < 0:
            raise ValueError("min_depth must be non-negative")

        self.num_samples = num_samples
        self.seq_len = seq_len
        self.mutation_rate = mutation_rate
        self.noise_rate = noise_rate
        self.max_insertion_len = max_insertion_len
        self.high_quality_mean = high_quality_mean
        self.low_quality_mean = low_quality_mean
        self.quality_spread = quality_spread
        self.max_phred = max_phred
        self.mean_depth = mean_depth
        self.min_depth = min_depth
        self.seed = seed
        self.tokenizer = DNATokenizer()

        self._rng = random.Random(seed) if seed is None else random.Random()

    @classmethod
    def from_config(cls, config: DataConfig, split: str) -> "SyntheticVariantDataset":
        """Build the train or validation split described by a config.

        Args:
            config: Data configuration.
            split: Either ``"train"`` or ``"val"``.

        Returns:
            A dataset instance for the requested split.

        Raises:
            ValueError: If ``split`` is not ``"train"`` or ``"val"``.
        """
        common = {
            "seq_len": config.seq_len,
            "mutation_rate": config.mutation_rate,
            "noise_rate": config.noise_rate,
            "max_insertion_len": config.max_insertion_len,
            "high_quality_mean": config.high_quality_mean,
            "low_quality_mean": config.low_quality_mean,
            "quality_spread": config.quality_spread,
            "max_phred": config.max_phred,
            "mean_depth": config.mean_depth,
            "min_depth": config.min_depth,
        }
        if split == "train":
            return cls(config.train_samples, seed=config.train_seed, **common)
        if split == "val":
            return cls(config.val_samples, seed=config.val_seed, **common)
        raise ValueError("split must be 'train' or 'val'")

    def __len__(self) -> int:
        """Return the nominal number of samples."""
        return self.num_samples

    def _random_base(self, exclude: str | None = None) -> str:
        """Draw a uniformly random base, optionally excluding one.

        Args:
            exclude: Base to omit, or ``None`` to allow all four.

        Returns:
            A single character from :data:`~.config.NUCLEOTIDES`.
        """
        choices = [base for base in NUCLEOTIDES if base != exclude] if exclude else list(NUCLEOTIDES)
        return self._rng.choice(choices)

    def _generate_reference(self) -> str:
        """Draw an i.i.d. uniform reference window.

        The window is drawn longer than ``seq_len`` so that insertion
        columns can be added without the sample running short.

        Returns:
            A reference nucleotide string.
        """
        return "".join(self._rng.choice(NUCLEOTIDES) for _ in range(self.seq_len))

    def _draw_quality(self, is_error: bool) -> float:
        """Draw a Phred score conditioned on whether the base is miscalled.

        This is the crux of noise/SNP separation. A sequencing error and a
        true SNP are identical in the base channels — both are one observed
        base disagreeing with the reference. They differ only in confidence:
        a sequencer reports low Phred where it is unsure, and errors
        concentrate there. Without this channel the two classes are not
        separable by any model.

        The separation planted here (means 36 against 12) is cleaner than
        reality. On real Illumina data the distributions overlap
        substantially, so expect quality to help materially but not to make
        the problem trivial.

        Args:
            is_error: Whether this column carries a sequencing error.

        Returns:
            A Phred score clamped to ``[0, max_phred]``.
        """
        mean = self.low_quality_mean if is_error else self.high_quality_mean
        score = self._rng.gauss(mean, self.quality_spread)
        return float(min(max(score, 0.0), float(self.max_phred)))

    def _draw_depth(self) -> float:
        """Draw a per-column read depth.

        Implemented as a Poisson via a Knuth sampler over the stdlib RNG, so
        the whole generator remains seedable through one object.

        Note:
            Depth is drawn independently of whether a column carries a
            variant or an error, which is honest for a single-read
            representation but means depth alone carries little
            discriminative signal here. It becomes informative in a pileup
            representation, where a true variant recurs across reads at a
            locus and an error does not. The ablation in ``stress_test.py``
            reports its actual contribution rather than assuming one.

        Returns:
            A depth value floored at ``min_depth``.
        """
        target = pow(2.718281828459045, -self.mean_depth)
        count, product = 0, 1.0
        while True:
            product *= self._rng.random()
            if product <= target:
                break
            count += 1
            if count > 10 * self.mean_depth + 50:  # numerical guard
                break
        return float(max(count, self.min_depth))

    def _build_alignment(
        self, reference: str
    ) -> Tuple[List[str], List[str], List[int], List[bool]]:
        """Emit gap-padded reference and observed columns with labels.

        Both lists grow in lockstep, one entry per alignment column, so the
        1-to-1 coordinate guarantee holds by construction rather than by a
        later reconciliation step.

        Args:
            reference: The unmutated reference window to walk.

        Returns:
            Tuple of ``(reference_columns, observed_columns, labels,
            error_flags)``. ``error_flags`` marks which columns carry a
            sequencing error, and is consumed only by the quality draw — it
            never reaches the labels.
        """
        reference_columns: List[str] = []
        observed_columns: List[str] = []
        labels: List[int] = []

        for reference_base in reference:
            if self._rng.random() >= self.mutation_rate:
                reference_columns.append(reference_base)
                observed_columns.append(reference_base)
                labels.append(LABEL_NORMAL)
                continue

            kind = self._rng.choice(MUTATION_KINDS)

            if kind == "snp":
                reference_columns.append(reference_base)
                observed_columns.append(self._random_base(exclude=reference_base))
                labels.append(LABEL_SNP)

            elif kind == "insertion":
                # Anchor column carries the reference base in both streams
                # and stays Normal; the inserted bases occupy gap columns.
                reference_columns.append(reference_base)
                observed_columns.append(reference_base)
                labels.append(LABEL_NORMAL)

                for _ in range(self._rng.randint(1, self.max_insertion_len)):
                    reference_columns.append(GAP_TOKEN)
                    observed_columns.append(self._random_base())
                    labels.append(LABEL_INSERTION)

            else:  # deletion
                reference_columns.append(reference_base)
                observed_columns.append(GAP_TOKEN)
                labels.append(LABEL_DELETION)

        # Sequencing noise is overlaid after variant simulation, on real
        # base calls only: a gap column has no base for the instrument to
        # miscall. Labels are deliberately left untouched.
        error_flags = [False] * len(observed_columns)
        if self.noise_rate > 0.0:
            for column, base in enumerate(observed_columns):
                if base == GAP_TOKEN:
                    continue
                if self._rng.random() < self.noise_rate:
                    observed_columns[column] = self._random_base(exclude=base)
                    error_flags[column] = True

        return reference_columns, observed_columns, labels, error_flags

    def _fit_to_length(
        self,
        reference_columns: List[str],
        observed_columns: List[str],
        labels: List[int],
        error_flags: List[bool],
    ) -> Tuple[List[str], List[str], List[int], List[bool]]:
        """Truncate or pad every channel together to ``seq_len`` columns.

        Padding uses ``N`` in both streams and ``Normal`` in the labels, so
        padded columns are inert rather than spurious variants.

        Args:
            reference_columns: Reference column characters.
            observed_columns: Observed column characters.
            labels: Per-column class ids.
            error_flags: Per-column sequencing-error markers.

        Returns:
            The four channels, each of exactly ``seq_len`` elements.
        """
        length = len(reference_columns)

        if length > self.seq_len:
            end = self.seq_len
            return (
                reference_columns[:end],
                observed_columns[:end],
                labels[:end],
                error_flags[:end],
            )

        if length < self.seq_len:
            pad = self.seq_len - length
            return (
                reference_columns + ["N"] * pad,
                observed_columns + ["N"] * pad,
                labels + [LABEL_NORMAL] * pad,
                error_flags + [False] * pad,
            )

        return reference_columns, observed_columns, labels, error_flags

    def __getitem__(self, index: int) -> VariantSample:
        """Generate one gap-padded sample with all four channels.

        Args:
            index: Positional index, used to derive the per-sample seed.

        Returns:
            A :class:`VariantSample`.

        Raises:
            IndexError: If ``index`` is outside ``[0, len(self))``.
        """
        if not 0 <= index < self.num_samples:
            raise IndexError(f"index {index} out of range for {self.num_samples} samples")

        if self.seed is not None:
            self._rng = random.Random(self.seed * self._INDEX_SEED_STRIDE + index)

        reference = self._generate_reference()
        columns = self._build_alignment(reference)
        reference_columns, observed_columns, labels, error_flags = self._fit_to_length(*columns)

        # Gap and pad columns have no base call, so they carry neither a
        # meaningful Phred score nor coverage; zero is the natural sentinel
        # and the model learns to read it alongside the gap token itself.
        qualities: List[float] = []
        depths: List[float] = []
        for column, base in enumerate(observed_columns):
            if base in (GAP_TOKEN, "N"):
                qualities.append(0.0)
                depths.append(0.0)
            else:
                qualities.append(self._draw_quality(error_flags[column]))
                depths.append(self._draw_depth())

        return VariantSample(
            reference_ids=self.tokenizer.encode("".join(reference_columns)),
            input_ids=self.tokenizer.encode("".join(observed_columns)),
            labels=torch.tensor(labels, dtype=torch.long),
            base_quality=torch.tensor(qualities, dtype=torch.float),
            depth=torch.tensor(depths, dtype=torch.float),
            # Zeros, not label-derived values: see the VariantSample docstring.
            vaf=torch.zeros(self.seq_len, dtype=torch.float),
        )


def build_dataloaders(
    config: DataConfig,
) -> Tuple[
    DataLoader[VariantSample],
    DataLoader[VariantSample],
    SyntheticVariantDataset,
    SyntheticVariantDataset,
]:
    """Construct the training and validation datasets and their loaders.

    Args:
        config: Data configuration.

    Returns:
        Tuple of ``(train_loader, val_loader, train_dataset, val_dataset)``.
    """
    train_dataset = SyntheticVariantDataset.from_config(config, "train")
    val_dataset = SyntheticVariantDataset.from_config(config, "val")

    train_loader: DataLoader[VariantSample] = DataLoader(
        train_dataset,
        batch_size=config.batch_size,
        shuffle=True,
        num_workers=config.num_workers,
    )
    val_loader: DataLoader[VariantSample] = DataLoader(
        val_dataset,
        batch_size=config.batch_size,
        shuffle=False,
        num_workers=config.num_workers,
    )
    return train_loader, val_loader, train_dataset, val_dataset


__all__ = [
    "DNATokenizer",
    "SyntheticVariantDataset",
    "VariantSample",
    "build_dataloaders",
]
