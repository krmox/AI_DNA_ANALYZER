"""Reference-guided bidirectional Mamba variant caller (proof of concept).

A per-token classifier that labels aligned DNA positions as Normal, SNP,
Insertion or Deletion by comparing an observed read against a reference.

Scope note: the ``BiMambaBlock`` is a causal gated-convolutional stand-in
for a true selective state-space model, and the training data is synthetic
with a deliberately regular insertion pattern. Metrics from this package are
a ceiling for that synthetic setup, not a prediction of real-read accuracy.

Typical use::

    from bimamba_variant_caller import ExperimentConfig, Trainer

    trainer = Trainer(ExperimentConfig())
    trainer.fit()
    val_loss, metrics = trainer.evaluate()
"""

from __future__ import annotations

from .config import (
    LABEL_NAMES,
    DataConfig,
    ExperimentConfig,
    ModelConfig,
    TrainingConfig,
)
from .dataset import DNATokenizer, SyntheticVariantDataset, VariantSample, build_dataloaders
from .loss import FocalLoss, compute_calibrated_class_weights
from .model import BiMambaBlock, CausalGatedConvCore, VariantCaller
from .train import EarlyStopping, EpochReport, Trainer, compute_metrics, train_epoch, validate

__version__ = "1.0.0"

__all__ = [
    "__version__",
    "LABEL_NAMES",
    "DataConfig",
    "ModelConfig",
    "TrainingConfig",
    "ExperimentConfig",
    "DNATokenizer",
    "SyntheticVariantDataset",
    "VariantSample",
    "build_dataloaders",
    "CausalGatedConvCore",
    "BiMambaBlock",
    "VariantCaller",
    "FocalLoss",
    "compute_calibrated_class_weights",
    "EarlyStopping",
    "EpochReport",
    "Trainer",
    "compute_metrics",
    "train_epoch",
    "validate",
]
