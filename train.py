"""Training loop, metrics and early stopping.

Metric choice deserves a note, because it was the subject of a real bug.
Early stopping monitors ``val_loss``, not ``f1_mutations_macro``. The
mutation F1 starts at exactly zero and can sit there while the model is in
fact improving, which caused premature stops. ``macro_f1`` is no better as a
stopping signal: ``Normal`` alone reaches ~0.99 almost immediately and drags
the mean up regardless of mutation performance. ``val_loss`` moves from the
first epoch and reflects progress independently of whether the rare classes
have started separating yet.

All F1 computation is implemented directly on tensors — no scikit-learn, in
keeping with the torch-only dependency constraint.
"""

from __future__ import annotations

import copy
import time
from dataclasses import dataclass
from typing import Dict, List, Literal, Tuple

import torch
import torch.nn as nn
from torch.optim import Optimizer
from torch.optim.lr_scheduler import CosineAnnealingLR, LinearLR, LRScheduler, SequentialLR
from torch.utils.data import DataLoader

from .config import LABEL_NAMES, ExperimentConfig
from .dataset import VariantSample, build_dataloaders
from .loss import FocalLoss, compute_calibrated_class_weights
from .model import VariantCaller


@dataclass
class EpochReport:
    """Metrics captured for a single epoch.

    Attributes:
        epoch: 1-based epoch number.
        learning_rate: Learning rate in effect during the epoch.
        train_loss: Token-weighted mean training loss.
        val_loss: Token-weighted mean validation loss.
        metrics: Metric name to value, as returned by :func:`compute_metrics`.
        seconds: Wall-clock duration of the epoch.
    """

    epoch: int
    learning_rate: float
    train_loss: float
    val_loss: float
    metrics: Dict[str, float]
    seconds: float

    def format_line(self) -> str:
        """Render the report as two aligned console lines.

        Returns:
            A formatted multi-line summary string.
        """
        per_class = ", ".join(
            f"F1[{name}]={self.metrics[f'f1_class_{index}']:.3f}"
            for index, name in enumerate(LABEL_NAMES)
        )
        header = (
            f"Epoch {self.epoch:3d} | lr {self.learning_rate:.6f} | "
            f"train {self.train_loss:.4f} | val {self.val_loss:.4f} | "
            f"macro F1 {self.metrics['macro_f1']:.4f} | "
            f"mutations F1 {self.metrics['f1_mutations_macro']:.4f} | "
            f"{self.seconds:.1f}s"
        )
        return f"{header}\n            {per_class}"


def compute_metrics(
    predictions: torch.Tensor,
    targets: torch.Tensor,
    num_classes: int = len(LABEL_NAMES),
) -> Dict[str, float]:
    """Compute per-class precision, recall and F1 plus two macro averages.

    Args:
        predictions: Long tensor ``[num_tokens]`` of predicted class ids.
        targets: Long tensor ``[num_tokens]`` of true class ids.
        num_classes: Number of classes.

    Returns:
        Mapping containing ``precision_class_{i}``, ``recall_class_{i}`` and
        ``f1_class_{i}`` for each class, plus ``macro_f1`` over all classes
        and ``f1_mutations_macro`` over the mutation classes only. The latter
        is the number to report: ``macro_f1`` is inflated by the dominant
        ``Normal`` class.

    Raises:
        ValueError: If the two tensors differ in length.
    """
    if predictions.numel() != targets.numel():
        raise ValueError(
            f"predictions ({predictions.numel()}) and targets ({targets.numel()}) differ in length"
        )

    predictions = predictions.reshape(-1)
    targets = targets.reshape(-1)

    metrics: Dict[str, float] = {}
    f1_scores: List[float] = []

    for class_id in range(num_classes):
        predicted_positive = predictions == class_id
        actual_positive = targets == class_id

        true_positives = float((predicted_positive & actual_positive).sum().item())
        false_positives = float((predicted_positive & ~actual_positive).sum().item())
        false_negatives = float((~predicted_positive & actual_positive).sum().item())

        precision = (
            true_positives / (true_positives + false_positives)
            if true_positives + false_positives > 0
            else 0.0
        )
        recall = (
            true_positives / (true_positives + false_negatives)
            if true_positives + false_negatives > 0
            else 0.0
        )
        f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0

        metrics[f"precision_class_{class_id}"] = precision
        metrics[f"recall_class_{class_id}"] = recall
        metrics[f"f1_class_{class_id}"] = f1
        f1_scores.append(f1)

    metrics["macro_f1"] = sum(f1_scores) / len(f1_scores)
    mutation_f1 = f1_scores[1:]
    metrics["f1_mutations_macro"] = sum(mutation_f1) / len(mutation_f1) if mutation_f1 else 0.0
    return metrics


class EarlyStopping:
    """Stops training once a monitored metric stops improving.

    Keeps a deep copy of the best-scoring weights so the final model is the
    best one seen, not merely the last.
    """

    def __init__(
        self,
        patience: int = 5,
        min_delta: float = 1e-4,
        mode: Literal["min", "max"] = "min",
    ) -> None:
        """Initialise the tracker.

        Args:
            patience: Epochs without improvement tolerated before stopping.
            min_delta: Smallest change counted as an improvement.
            mode: ``"min"`` for losses, ``"max"`` for scores.

        Raises:
            ValueError: If ``patience`` is non-positive or ``mode`` is invalid.
        """
        if patience <= 0:
            raise ValueError("patience must be positive")
        if mode not in ("min", "max"):
            raise ValueError("mode must be 'min' or 'max'")

        self.patience = patience
        self.min_delta = min_delta
        self.mode = mode
        self.best_score: float = float("inf") if mode == "min" else float("-inf")
        self.best_epoch: int = 0
        self.best_state_dict: Dict[str, torch.Tensor] | None = None
        self.counter: int = 0
        self.should_stop: bool = False

    def step(self, score: float, model: nn.Module, epoch: int) -> bool:
        """Record a score and update stopping state.

        Args:
            score: The monitored value for this epoch.
            model: Model whose weights are snapshotted on improvement.
            epoch: 1-based epoch number, stored for reporting.

        Returns:
            ``True`` if this epoch improved on the best score so far.
        """
        improved = (
            score < self.best_score - self.min_delta
            if self.mode == "min"
            else score > self.best_score + self.min_delta
        )

        if improved:
            self.best_score = score
            self.best_epoch = epoch
            self.best_state_dict = copy.deepcopy(model.state_dict())
            self.counter = 0
        else:
            self.counter += 1
            self.should_stop = self.counter >= self.patience

        return improved

    def restore_best(self, model: nn.Module) -> None:
        """Load the best snapshot back into the model, if one exists.

        Args:
            model: Model to restore in place.
        """
        if self.best_state_dict is not None:
            model.load_state_dict(self.best_state_dict)


def build_scheduler(optimizer: Optimizer, warmup_epochs: int, total_epochs: int) -> LRScheduler:
    """Build a linear-warmup then cosine-annealing schedule.

    Warmup keeps the first, highest-variance updates small; cosine decay then
    anneals smoothly rather than stepping.

    Args:
        optimizer: Optimizer to schedule.
        warmup_epochs: Epochs spent ramping from 10% to 100% of the base rate.
        total_epochs: Total planned epochs.

    Returns:
        A scheduler to be stepped once per epoch.

    Raises:
        ValueError: If ``warmup_epochs`` is not strictly less than
            ``total_epochs``.
    """
    if not 0 <= warmup_epochs < total_epochs:
        raise ValueError("warmup_epochs must satisfy 0 <= warmup_epochs < total_epochs")

    if warmup_epochs == 0:
        return CosineAnnealingLR(optimizer, T_max=total_epochs)

    return SequentialLR(
        optimizer,
        schedulers=[
            LinearLR(optimizer, start_factor=0.1, end_factor=1.0, total_iters=warmup_epochs),
            CosineAnnealingLR(optimizer, T_max=max(1, total_epochs - warmup_epochs)),
        ],
        milestones=[warmup_epochs],
    )


def train_epoch(
    model: VariantCaller,
    dataloader: DataLoader[VariantSample],
    optimizer: Optimizer,
    criterion: nn.Module,
    device: torch.device,
    grad_clip_norm: float | None = 1.0,
) -> float:
    """Run one training pass over the loader.

    Args:
        model: Model to train, switched to train mode internally.
        dataloader: Training batches.
        optimizer: Optimizer to step.
        criterion: Loss taking flattened ``[tokens, classes]`` logits.
        device: Device to move batches onto.
        grad_clip_norm: Global gradient-norm clip, or ``None`` to skip.

    Returns:
        Token-weighted mean training loss. Weighting by token count keeps a
        short final batch from skewing the average.
    """
    model.train()
    total_loss = 0.0
    total_tokens = 0

    for batch in dataloader:
        input_ids = batch["input_ids"].to(device)
        reference_ids = batch["reference_ids"].to(device)
        labels = batch["labels"].to(device)
        base_quality = batch["base_quality"].to(device)
        depth = batch["depth"].to(device)

        optimizer.zero_grad(set_to_none=True)
        logits = model(input_ids, reference_ids, base_quality, depth)
        loss = criterion(logits.reshape(-1, logits.size(-1)), labels.reshape(-1))
        loss.backward()

        if grad_clip_norm is not None:
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=grad_clip_norm)
        optimizer.step()

        token_count = labels.numel()
        total_loss += loss.item() * token_count
        total_tokens += token_count

    return total_loss / max(total_tokens, 1)


@torch.no_grad()
def validate(
    model: VariantCaller,
    dataloader: DataLoader[VariantSample],
    criterion: nn.Module,
    device: torch.device,
    num_classes: int = len(LABEL_NAMES),
) -> Tuple[float, Dict[str, float]]:
    """Evaluate the model over the validation loader.

    Args:
        model: Model to evaluate, switched to eval mode internally.
        dataloader: Validation batches.
        criterion: Loss taking flattened ``[tokens, classes]`` logits.
        device: Device to move batches onto.
        num_classes: Number of classes for metric computation.

    Returns:
        Tuple of ``(val_loss, metrics)``.
    """
    model.eval()
    total_loss = 0.0
    total_tokens = 0
    all_predictions: List[torch.Tensor] = []
    all_targets: List[torch.Tensor] = []

    for batch in dataloader:
        input_ids = batch["input_ids"].to(device)
        reference_ids = batch["reference_ids"].to(device)
        labels = batch["labels"].to(device)
        base_quality = batch["base_quality"].to(device)
        depth = batch["depth"].to(device)

        logits = model(input_ids, reference_ids, base_quality, depth)
        loss = criterion(logits.reshape(-1, logits.size(-1)), labels.reshape(-1))

        token_count = labels.numel()
        total_loss += loss.item() * token_count
        total_tokens += token_count

        all_predictions.append(logits.argmax(dim=-1).reshape(-1).cpu())
        all_targets.append(labels.reshape(-1).cpu())

    metrics = compute_metrics(torch.cat(all_predictions), torch.cat(all_targets), num_classes)
    return total_loss / max(total_tokens, 1), metrics


def resolve_device(preference: str) -> torch.device:
    """Resolve a device preference string to a concrete device.

    Args:
        preference: ``"auto"``, ``"cpu"`` or ``"cuda"``.

    Returns:
        The resolved device; ``"auto"`` selects CUDA when available.
    """
    if preference == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(preference)


class Trainer:
    """Owns the full training lifecycle for one experiment config."""

    def __init__(self, config: ExperimentConfig, verbose: bool = True) -> None:
        """Assemble data, model, loss, optimizer and schedule.

        Args:
            config: Experiment configuration.
            verbose: Whether to print per-epoch progress.
        """
        self.config = config
        self.verbose = verbose
        self.device = resolve_device(config.training.device)

        torch.manual_seed(config.training.seed)

        (
            self.train_loader,
            self.val_loader,
            self.train_dataset,
            self.val_dataset,
        ) = build_dataloaders(config.data)

        self.model = VariantCaller(config.model).to(self.device)

        self.class_weights = compute_calibrated_class_weights(
            self.train_dataset,
            num_classes=config.model.num_classes,
            sample_size=config.training.weight_sample_size,
            max_weight_cap=config.training.max_weight_cap,
            alpha=config.training.class_weight_alpha,
        ).to(self.device)

        self.criterion = FocalLoss(
            alpha=self.class_weights,
            gamma=config.training.focal_gamma,
            label_smoothing=config.training.label_smoothing,
        ).to(self.device)

        self.optimizer: Optimizer = torch.optim.AdamW(
            self.model.parameters(),
            lr=config.training.learning_rate,
            weight_decay=config.training.weight_decay,
        )
        self.scheduler = build_scheduler(
            self.optimizer, config.training.warmup_epochs, config.training.num_epochs
        )
        self.early_stopping = EarlyStopping(
            patience=config.training.early_stopping_patience,
            min_delta=config.training.early_stopping_min_delta,
            mode="min",
        )
        self.history: List[EpochReport] = []

    def _log(self, message: str) -> None:
        """Print a message when verbose mode is on.

        Args:
            message: Text to print.
        """
        if self.verbose:
            print(message)

    def fit(self) -> List[EpochReport]:
        """Train until convergence or the epoch budget is exhausted.

        Restores the best-``val_loss`` weights before returning.

        Returns:
            The per-epoch history.
        """
        self._log(f"Device: {self.device}")
        self._log(f"Parameters: {self.model.count_parameters():,}")
        weights = ", ".join(
            f"{name}={weight:.3f}" for name, weight in zip(LABEL_NAMES, self.class_weights.tolist())
        )
        self._log(f"Class weights: {weights}\n")

        for epoch in range(1, self.config.training.num_epochs + 1):
            started = time.perf_counter()
            learning_rate = self.optimizer.param_groups[0]["lr"]

            train_loss = train_epoch(
                self.model,
                self.train_loader,
                self.optimizer,
                self.criterion,
                self.device,
                self.config.training.grad_clip_norm,
            )
            val_loss, metrics = validate(
                self.model,
                self.val_loader,
                self.criterion,
                self.device,
                self.config.model.num_classes,
            )
            self.scheduler.step()

            report = EpochReport(
                epoch=epoch,
                learning_rate=learning_rate,
                train_loss=train_loss,
                val_loss=val_loss,
                metrics=metrics,
                seconds=time.perf_counter() - started,
            )
            self.history.append(report)
            self._log(report.format_line())

            self.early_stopping.step(val_loss, self.model, epoch)
            if self.early_stopping.should_stop:
                self._log(
                    f"\nEarly stopping at epoch {epoch}: no val_loss improvement for "
                    f"{self.config.training.early_stopping_patience} epochs."
                )
                break

        self.early_stopping.restore_best(self.model)
        self._log(
            f"Restored best weights from epoch {self.early_stopping.best_epoch} "
            f"(val_loss={self.early_stopping.best_score:.4f})."
        )
        return self.history

    def evaluate(self) -> Tuple[float, Dict[str, float]]:
        """Evaluate the current weights on the validation split.

        Returns:
            Tuple of ``(val_loss, metrics)``.
        """
        return validate(
            self.model,
            self.val_loader,
            self.criterion,
            self.device,
            self.config.model.num_classes,
        )


__all__ = [
    "EpochReport",
    "EarlyStopping",
    "Trainer",
    "build_scheduler",
    "compute_metrics",
    "resolve_device",
    "train_epoch",
    "validate",
]
