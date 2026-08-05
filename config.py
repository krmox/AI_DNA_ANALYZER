"""Configuration dataclasses for the BiMamba variant caller.

All tunable knobs live here so that experiments differ only by config,
never by edits scattered across modules. Every dataclass validates its
own invariants in ``__post_init__`` so an invalid experiment fails fast
at construction time rather than midway through training.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Final, List, Literal, Tuple

# ---------------------------------------------------------------------------
# Alphabet / label constants (single source of truth for the whole package)
# ---------------------------------------------------------------------------

#: Nucleotide vocabulary. ``N`` (id 0) is padding and unknown base. ``-``
#: (id 5) is the alignment gap introduced by CIGAR-style padding: it appears
#: in the reference stream opposite an insertion and in the observed stream
#: opposite a deletion. Keeping gap distinct from ``N`` matters — ``N`` means
#: "a base was read but is unreliable", ``-`` means "no base exists at this
#: coordinate in this stream". Collapsing them would make a low-quality call
#: look like a deletion.
VOCAB: Final[Dict[str, int]] = {"N": 0, "A": 1, "C": 2, "G": 3, "T": 4, "-": 5}

#: Inverse of :data:`VOCAB`, used for decoding tensors back to strings.
INV_VOCAB: Final[Dict[int, str]] = {index: token for token, index in VOCAB.items()}

#: The alignment gap character and its token id.
GAP_TOKEN: Final[str] = "-"
GAP_ID: Final[int] = VOCAB[GAP_TOKEN]

#: Real bases that the generator may sample. ``N`` and ``-`` are deliberately
#: excluded: they are structural markers, not base calls.
NUCLEOTIDES: Final[Tuple[str, ...]] = ("A", "C", "G", "T")

#: Per-token variant classes predicted by the model.
LABEL_NORMAL: Final[int] = 0
LABEL_SNP: Final[int] = 1
LABEL_INSERTION: Final[int] = 2
LABEL_DELETION: Final[int] = 3

#: Human-readable names, indexed by class id. Used for reporting only.
LABEL_NAMES: Final[Tuple[str, ...]] = ("Normal", "SNP", "Insertion", "Deletion")

#: Mutation kinds the synthetic generator can emit, in a fixed order so that
#: a seeded RNG reproduces the same dataset across runs and platforms.
MUTATION_KINDS: Final[Tuple[str, ...]] = ("snp", "insertion", "deletion")

#: Fusion strategies for combining the reference and observed embeddings.
FusionMode = Literal["concat", "sum"]


@dataclass(frozen=True)
class DataConfig:
    """Parameters controlling synthetic pairwise read generation.

    Attributes:
        train_samples: Number of training sequences drawn per epoch.
        val_samples: Number of validation sequences.
        seq_len: Fixed token length of every emitted sequence pair.
        mutation_rate: Per-reference-position probability of introducing a
            true variant. The three mutation kinds are then chosen uniformly.
        noise_rate: Per-base probability of a sequencing substitution error
            in the observed channel. Illumina platforms sit around 0.001 to
            0.01 depending on cycle and chemistry; 0.005 to 0.01 is the range
            requested here. These errors are deliberately *not* labelled, so
            they place a hard ceiling on attainable SNP precision — see
            :meth:`~.dataset.SyntheticVariantDataset._apply_noise`.
        max_insertion_len: Inclusive upper bound on insertion length, drawn
            as ``Uniform(1, max_insertion_len)``.
        high_quality_mean: Mean Phred score emitted for a correctly read
            base. Typical Illumina Q30-Q40.
        low_quality_mean: Mean Phred score emitted for a sequencing error.
            The gap between this and ``high_quality_mean`` is precisely the
            signal that lets the model separate noise from a true SNP. It is
            a modelling assumption: on real data the two distributions
            overlap far more, so quality helps but does not fully separate.
        quality_spread: Standard deviation of the Phred draw around either
            mean, before clamping to ``[0, max_phred]``.
        max_phred: Upper clamp on emitted Phred scores.
        mean_depth: Mean per-position read depth, drawn from a Poisson.
        min_depth: Floor applied after the Poisson draw.
        val_seed: Seed for the validation generator. Fixing it keeps the
            validation set identical across epochs and across runs, which is
            what makes ``val_loss`` a trustworthy early-stopping signal.
        train_seed: Seed for the training generator, or ``None`` for a fresh
            random stream (the default: fresh samples every epoch act as
            free data augmentation).
        batch_size: Mini-batch size for both loaders.
        num_workers: DataLoader worker processes. Keep at 0 unless the
            generator becomes a bottleneck; the generator is CPU-cheap.
    """

    train_samples: int = 1000
    val_samples: int = 200
    seq_len: int = 64
    mutation_rate: float = 0.05
    noise_rate: float = 0.0075
    max_insertion_len: int = 6
    high_quality_mean: float = 36.0
    low_quality_mean: float = 12.0
    quality_spread: float = 5.0
    max_phred: int = 60
    mean_depth: float = 30.0
    min_depth: int = 1
    val_seed: int = 123
    train_seed: int | None = None
    batch_size: int = 16
    num_workers: int = 0

    def __post_init__(self) -> None:
        """Validate data invariants.

        Raises:
            ValueError: If any field is outside its permitted range.
        """
        if self.train_samples <= 0 or self.val_samples <= 0:
            raise ValueError("train_samples and val_samples must be positive")
        if self.seq_len <= 0:
            raise ValueError("seq_len must be positive")
        if not 0.0 <= self.mutation_rate <= 1.0:
            raise ValueError("mutation_rate must lie in [0, 1]")
        if not 0.0 <= self.noise_rate <= 1.0:
            raise ValueError("noise_rate must lie in [0, 1]")
        if self.max_insertion_len < 1:
            raise ValueError("max_insertion_len must be at least 1")
        if self.max_phred < 1:
            raise ValueError("max_phred must be at least 1")
        if self.quality_spread < 0.0:
            raise ValueError("quality_spread must be non-negative")
        if not 0.0 <= self.low_quality_mean <= self.max_phred:
            raise ValueError("low_quality_mean must lie in [0, max_phred]")
        if not 0.0 <= self.high_quality_mean <= self.max_phred:
            raise ValueError("high_quality_mean must lie in [0, max_phred]")
        if self.mean_depth <= 0.0:
            raise ValueError("mean_depth must be positive")
        if self.min_depth < 0:
            raise ValueError("min_depth must be non-negative")
        if self.batch_size <= 0:
            raise ValueError("batch_size must be positive")
        if self.num_workers < 0:
            raise ValueError("num_workers must be non-negative")


@dataclass(frozen=True)
class ModelConfig:
    """Architecture hyperparameters for :class:`~.model.VariantCaller`.

    Attributes:
        vocab_size: Size of the nucleotide vocabulary (5: N, A, C, G, T).
        d_model: Hidden width carried through every BiMamba block.
        num_layers: Number of stacked BiMamba blocks.
        num_classes: Number of per-token output classes (4).
        conv_kernel_size: Depthwise convolution width inside each directional
            core. This is the block's entire receptive field per layer, so it
            bounds how far apart two tokens can be and still interact.
        dropout: Dropout probability applied at the embedding fusion and
            inside every block.
        use_base_quality: Whether to consume the Phred channel. Turning this
            off is the ablation that shows how much of the SNP/noise
            separation quality is actually responsible for.
        use_depth: Whether to consume the coverage channel.
        quality_bins: Number of buckets the Phred embedding discretises into.
            Bucket width is ``max_phred / quality_bins``.
        max_phred: Upper Phred bound, must match the data config.
        fusion: How reference and observed embeddings are combined.
            ``"concat"`` projects ``[obs; ref]`` through a linear layer and is
            the validated default; ``"sum"`` adds them elementwise (cheaper,
            strictly less expressive).
        padding_idx: Vocabulary id treated as padding by the embeddings.
    """

    vocab_size: int = len(VOCAB)
    d_model: int = 128
    num_layers: int = 6
    num_classes: int = len(LABEL_NAMES)
    conv_kernel_size: int = 3
    dropout: float = 0.15
    use_base_quality: bool = True
    use_depth: bool = True
    quality_bins: int = 16
    max_phred: int = 60
    fusion: FusionMode = "concat"
    padding_idx: int = 0

    def __post_init__(self) -> None:
        """Validate architecture invariants.

        Raises:
            ValueError: If any field is outside its permitted range.
        """
        if self.d_model <= 0 or self.num_layers <= 0:
            raise ValueError("d_model and num_layers must be positive")
        if self.conv_kernel_size <= 0:
            raise ValueError("conv_kernel_size must be positive")
        if not 0.0 <= self.dropout < 1.0:
            raise ValueError("dropout must lie in [0, 1)")
        if self.fusion not in ("concat", "sum"):
            raise ValueError("fusion must be 'concat' or 'sum'")
        if self.quality_bins < 1:
            raise ValueError("quality_bins must be at least 1")
        if self.max_phred < 1:
            raise ValueError("max_phred must be at least 1")
        if self.num_classes != len(LABEL_NAMES):
            raise ValueError(f"num_classes must be {len(LABEL_NAMES)}")


@dataclass(frozen=True)
class TrainingConfig:
    """Optimisation, loss-calibration and early-stopping settings.

    Attributes:
        num_epochs: Upper bound on epochs; early stopping usually fires first.
        learning_rate: Peak AdamW learning rate reached after warmup.
        weight_decay: AdamW decoupled weight decay.
        warmup_epochs: Epochs of linear warmup before cosine annealing.
        grad_clip_norm: Global gradient-norm clip, or ``None`` to disable.
        focal_gamma: Focal-loss focusing exponent. Higher values down-weight
            already-confident tokens more aggressively.
        label_smoothing: Target smoothing. Kept at 0.0: smoothing caps the
            achievable confidence on the rare mutation classes and measurably
            slowed convergence in this setup.
        max_weight_cap: Upper clamp on per-class alpha weights after sqrt
            smoothing and normalisation. Prevents the runaway ~16x mutation
            weighting that previously drove mass false positives.
        weight_sample_size: Sequences sampled when estimating class frequencies.
        early_stopping_patience: Epochs without ``val_loss`` improvement before
            stopping.
        early_stopping_min_delta: Smallest ``val_loss`` decrease counted as an
            improvement.
        seed: Global torch seed for reproducible initialisation.
        device: ``"cuda"``, ``"cpu"``, or ``"auto"`` to pick CUDA when present.
    """

    num_epochs: int = 15
    learning_rate: float = 1e-3
    weight_decay: float = 1e-2
    warmup_epochs: int = 3
    grad_clip_norm: float | None = 1.0
    focal_gamma: float = 2.0
    label_smoothing: float = 0.0
    max_weight_cap: float = 6.0
    weight_sample_size: int = 200
    early_stopping_patience: int = 5
    early_stopping_min_delta: float = 1e-4
    seed: int = 42
    device: str = "auto"

    def __post_init__(self) -> None:
        """Validate training invariants.

        Raises:
            ValueError: If any field is outside its permitted range.
        """
        if self.num_epochs <= 0:
            raise ValueError("num_epochs must be positive")
        if self.learning_rate <= 0.0:
            raise ValueError("learning_rate must be positive")
        if self.weight_decay < 0.0:
            raise ValueError("weight_decay must be non-negative")
        if not 0 <= self.warmup_epochs < self.num_epochs:
            raise ValueError("warmup_epochs must satisfy 0 <= warmup_epochs < num_epochs")
        if self.grad_clip_norm is not None and self.grad_clip_norm <= 0.0:
            raise ValueError("grad_clip_norm must be positive or None")
        if self.focal_gamma < 0.0:
            raise ValueError("focal_gamma must be non-negative")
        if not 0.0 <= self.label_smoothing < 1.0:
            raise ValueError("label_smoothing must lie in [0, 1)")
        if self.max_weight_cap <= 0.0:
            raise ValueError("max_weight_cap must be positive")
        if self.weight_sample_size <= 0:
            raise ValueError("weight_sample_size must be positive")
        if self.early_stopping_patience <= 0:
            raise ValueError("early_stopping_patience must be positive")
        if self.device not in ("auto", "cpu", "cuda"):
            raise ValueError("device must be 'auto', 'cpu' or 'cuda'")


@dataclass(frozen=True)
class ExperimentConfig:
    """Top-level bundle tying the three configs into one experiment.

    Attributes:
        data: Synthetic data generation settings.
        model: Architecture settings.
        training: Optimisation and early-stopping settings.
    """

    data: DataConfig = field(default_factory=DataConfig)
    model: ModelConfig = field(default_factory=ModelConfig)
    training: TrainingConfig = field(default_factory=TrainingConfig)


__all__: List[str] = [
    "VOCAB",
    "GAP_TOKEN",
    "GAP_ID",
    "INV_VOCAB",
    "NUCLEOTIDES",
    "LABEL_NORMAL",
    "LABEL_SNP",
    "LABEL_INSERTION",
    "LABEL_DELETION",
    "LABEL_NAMES",
    "MUTATION_KINDS",
    "FusionMode",
    "DataConfig",
    "ModelConfig",
    "TrainingConfig",
    "ExperimentConfig",
]
