"""Synthetic pairwise read generation and PyTorch dataset plumbing.

The generator emits *aligned triples* — a reference channel, an observed
channel and a per-token label — all of exactly ``seq_len`` tokens. Keeping
the two channels index-aligned is what makes per-position classification
solvable at all: a single-channel model cannot distinguish a substituted
base from an ordinary one, because both are just a uniformly random letter.

Event encoding (per reference position):

======================================================================
Kind       Observed tokens        Reference tokens     Labels
======================================================================
Normal     ``[ref]``              ``[ref]``            ``[0]``
SNP        ``[base != ref]``      ``[ref]``            ``[1]``
Insertion  ``[ref, random]``      ``[ref, "N"]``       ``[2, 2]``
Deletion   ``["N"]``              ``[ref]``            ``[3]``
======================================================================

Insertions emit two tokens — a tandem duplicate of the reference base
followed by a fresh random base whose reference slot is an ``N`` filler.
That makes an insertion a *bigram* pattern (match immediately followed by
filler) rather than a lone mismatch, which is what separates it from a SNP.
Because insertions lengthen the stream, all three channels are truncated
(or padded) together to ``seq_len`` as the final step.

Caveat worth carrying into any real-data port: this fixed two-token shape
is a considerable simplification. Real insertions vary in length and are
not reliably tandem duplications, so metrics measured here are an upper
bound, not a transferable estimate.
"""

from __future__ import annotations

import random
from typing import Dict, List, Sequence, Tuple, TypedDict

import torch
from torch.utils.data import DataLoader, Dataset

from .config import (
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
        reference_ids: Long tensor ``[seq_len]`` of reference token ids.
        input_ids: Long tensor ``[seq_len]`` of observed token ids.
        labels: Long tensor ``[seq_len]`` of per-token class ids in ``0..3``.
    """

    reference_ids: torch.Tensor
    input_ids: torch.Tensor
    labels: torch.Tensor


class DNATokenizer:
    """Character-level nucleotide tokenizer over ``{N, A, C, G, T}``.

    ``N`` maps to id 0 and serves three distinct roles: padding, unknown
    base, and alignment gap. The model learns to read it from context.
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
            sequence: Nucleotide string such as ``"ACGTN"``.

        Returns:
            Long tensor of shape ``[len(sequence)]``.
        """
        return torch.tensor(
            [self.vocab.get(char, self.vocab["N"]) for char in sequence],
            dtype=torch.long,
        )

    def decode(self, token_ids: torch.Tensor | Sequence[int]) -> str:
        """Decode token ids back into a nucleotide string.

        Args:
            token_ids: 1-D tensor or sequence of integer token ids.

        Returns:
            The decoded nucleotide string. Unknown ids render as ``N``.
        """
        ids = token_ids.tolist() if isinstance(token_ids, torch.Tensor) else list(token_ids)
        return "".join(self.inv_vocab.get(int(index), "N") for index in ids)


class SyntheticVariantDataset(Dataset[VariantSample]):
    """Generates aligned (reference, observed, labels) triples on the fly.

    Samples are produced lazily in :meth:`__getitem__` and never cached.

    Seeding semantics are deliberately *per index*, not per stream. A single
    shared ``random.Random(seed)`` would be reproducible across runs but its
    state advances on every access, so a second pass over the same dataset
    yields different sequences. That silently made every epoch validate on
    fresh data and broke ``val_loss`` as an early-stopping signal. With a
    seed supplied, each index derives its own generator, so ``dataset[i]`` is
    stable across passes, across epochs and under shuffling.

    Note:
        With ``seed=None`` the generator streams fresh samples on every
        access, which acts as free augmentation for training. Do not use an
        unseeded dataset for validation.
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
    ) -> None:
        """Initialise the generator.

        Args:
            num_samples: Nominal dataset length reported by ``__len__``.
            seq_len: Fixed token length of every emitted channel.
            mutation_rate: Per-position probability of a mutation event.
            seed: Seed for the private RNG, or ``None`` for a fresh stream.

        Raises:
            ValueError: If ``num_samples`` or ``seq_len`` is non-positive, or
                ``mutation_rate`` falls outside ``[0, 1]``.
        """
        if num_samples <= 0:
            raise ValueError("num_samples must be positive")
        if seq_len <= 0:
            raise ValueError("seq_len must be positive")
        if not 0.0 <= mutation_rate <= 1.0:
            raise ValueError("mutation_rate must lie in [0, 1]")

        self.num_samples = num_samples
        self.seq_len = seq_len
        self.mutation_rate = mutation_rate
        self.seed = seed
        self.tokenizer = DNATokenizer()

        # Unseeded: one shared streaming generator. Seeded: a per-index
        # generator is built on demand in __getitem__ instead.
        self._rng = random.Random(seed) if seed is None else random.Random()

    @classmethod
    def from_config(cls, config: DataConfig, split: str) -> "SyntheticVariantDataset":
        """Build the train or validation split described by a config.

        Args:
            config: Data configuration to read sizes and seeds from.
            split: Either ``"train"`` or ``"val"``.

        Returns:
            A dataset instance for the requested split.

        Raises:
            ValueError: If ``split`` is not ``"train"`` or ``"val"``.
        """
        if split == "train":
            return cls(config.train_samples, config.seq_len, config.mutation_rate, config.train_seed)
        if split == "val":
            return cls(config.val_samples, config.seq_len, config.mutation_rate, config.val_seed)
        raise ValueError("split must be 'train' or 'val'")

    def __len__(self) -> int:
        """Return the nominal number of samples."""
        return self.num_samples

    def _random_base(self, exclude: str | None = None) -> str:
        """Draw a uniformly random base, optionally excluding one.

        Args:
            exclude: Base to omit from the draw, or ``None`` to allow all.

        Returns:
            A single character from :data:`~.config.NUCLEOTIDES`.
        """
        choices = [base for base in NUCLEOTIDES if base != exclude] if exclude else list(NUCLEOTIDES)
        return self._rng.choice(choices)

    def _generate_reference(self) -> str:
        """Draw an i.i.d. uniform reference sequence of length ``seq_len``.

        Returns:
            A reference nucleotide string.
        """
        return "".join(self._rng.choice(NUCLEOTIDES) for _ in range(self.seq_len))

    def _apply_mutations(self, reference: str) -> Tuple[str, str, List[int]]:
        """Walk the reference and emit aligned observed/reference/label streams.

        Insertions append two tokens to every stream, so the intermediate
        length may exceed ``seq_len``; :meth:`_fit_to_length` reconciles that
        afterwards. All three streams are appended to in lockstep, which is
        what guarantees index alignment.

        Args:
            reference: The unmutated reference string to walk.

        Returns:
            Tuple of ``(observed, aligned_reference, labels)``, each of the
            same (not yet truncated) length.
        """
        observed: List[str] = []
        aligned_reference: List[str] = []
        labels: List[int] = []

        for reference_base in reference:
            if self._rng.random() >= self.mutation_rate:
                observed.append(reference_base)
                aligned_reference.append(reference_base)
                labels.append(LABEL_NORMAL)
                continue

            kind = self._rng.choice(MUTATION_KINDS)

            if kind == "snp":
                # Guaranteed mismatch: the model must key on ref != obs.
                observed.append(self._random_base(exclude=reference_base))
                aligned_reference.append(reference_base)
                labels.append(LABEL_SNP)

            elif kind == "insertion":
                # Token 1 — tandem duplicate; reference matches observed here,
                # so this token alone looks Normal. Only the pair is decisive.
                observed.append(reference_base)
                aligned_reference.append(reference_base)
                labels.append(LABEL_INSERTION)

                # Token 2 — novel base with no reference counterpart, marked
                # by an N filler in the reference channel.
                observed.append(self._random_base())
                aligned_reference.append("N")
                labels.append(LABEL_INSERTION)

            else:  # deletion
                # Gap in the observed channel; reference retains the base.
                observed.append("N")
                aligned_reference.append(reference_base)
                labels.append(LABEL_DELETION)

        return "".join(observed), "".join(aligned_reference), labels

    def _fit_to_length(
        self, observed: str, aligned_reference: str, labels: List[int]
    ) -> Tuple[str, str, List[int]]:
        """Truncate or pad all three channels together to ``seq_len``.

        Padding uses ``N`` in both sequence channels and ``Normal`` in the
        labels, so padded positions are inert rather than spurious variants.

        Args:
            observed: Observed nucleotide string.
            aligned_reference: Reference string aligned to ``observed``.
            labels: Per-token class ids aligned to both channels.

        Returns:
            The three channels, each of exactly ``seq_len`` elements.
        """
        length = len(observed)

        if length > self.seq_len:
            return observed[: self.seq_len], aligned_reference[: self.seq_len], labels[: self.seq_len]

        if length < self.seq_len:
            pad = self.seq_len - length
            return (
                observed + "N" * pad,
                aligned_reference + "N" * pad,
                labels + [LABEL_NORMAL] * pad,
            )

        return observed, aligned_reference, labels

    def __getitem__(self, index: int) -> VariantSample:
        """Generate one aligned sample.

        Args:
            index: Positional index. Ignored beyond bounds checking — samples
                are generated fresh rather than looked up.

        Returns:
            A :class:`VariantSample` with three ``[seq_len]`` long tensors.

        Raises:
            IndexError: If ``index`` is outside ``[0, len(self))``.
        """
        if not 0 <= index < self.num_samples:
            raise IndexError(f"index {index} out of range for {self.num_samples} samples")

        if self.seed is not None:
            # Rebind to an index-derived generator so this sample is identical
            # on every access, regardless of traversal order.
            self._rng = random.Random(self.seed * self._INDEX_SEED_STRIDE + index)

        reference = self._generate_reference()
        observed, aligned_reference, labels = self._apply_mutations(reference)
        observed, aligned_reference, labels = self._fit_to_length(observed, aligned_reference, labels)

        return VariantSample(
            reference_ids=self.tokenizer.encode(aligned_reference),
            input_ids=self.tokenizer.encode(observed),
            labels=torch.tensor(labels, dtype=torch.long),
        )


def build_dataloaders(
    config: DataConfig,
) -> Tuple[DataLoader[VariantSample], DataLoader[VariantSample], SyntheticVariantDataset, SyntheticVariantDataset]:
    """Construct the training and validation datasets and their loaders.

    The validation loader is unshuffled so that per-epoch metrics are
    computed over an identically ordered set.

    Args:
        config: Data configuration.

    Returns:
        Tuple of ``(train_loader, val_loader, train_dataset, val_dataset)``.
        The datasets are returned alongside the loaders because class-weight
        estimation needs direct indexed access.
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
