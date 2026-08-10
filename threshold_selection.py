"""Post-hoc, per-checkpoint decision-threshold selection (Option 2).

Read-only analysis utility: given a trained checkpoint and a validation
set, sweep candidate probability thresholds and report the one that
maximizes ``mutation_macro_f1`` for that specific checkpoint.

This does not change how ``run_train.py`` picks a checkpoint during
training (that still keys off ``VARIANT_PROB_THRESHOLD`` = 0.5, fixed).
It answers a narrower question after the fact: at its own best operating
point, how good was a given checkpoint actually?

Motivation: AUC/AP measure ranking quality and are threshold-independent,
so they stay flat even as F1@0.5 collapses when the background score
distribution's tail drifts upward under continued training (see
History/3_DEVLOG.md, 2026-08-08 and 2026-08-10). A per-checkpoint optimal
threshold isolates "is the model's ranking still good" from "did we pick
the right cutoff for this particular epoch's calibration."
"""

from __future__ import annotations

from typing import NamedTuple

import torch
from torch.utils.data import DataLoader

from config import LABEL_NAMES
from run_train import (
    filter_batch_for_model,
    move_batch_to_device,
    per_class_precision_recall_f1,
    predict_with_threshold,
    torch_confusion_matrix,
)

NUM_CLASSES = len(LABEL_NAMES)

#: Default sweep grid. Finer than 0.01 near the low end, where the
#: tail-drift failure mode lives (see module docstring); coarser above
#: 0.9 where nothing interesting has been observed to happen.
DEFAULT_THRESHOLDS: tuple[float, ...] = tuple(
    [round(0.5 + 0.01 * i, 2) for i in range(0, 45)]  # 0.50 .. 0.94
    + [round(0.95 + 0.01 * i, 2) for i in range(0, 5)]  # 0.95 .. 0.99
)


class ThresholdSweepResult(NamedTuple):
    """One checkpoint's threshold sweep outcome.

    Attributes:
        f1_at_fixed: mutation_macro_f1 at the fixed reference threshold
            (0.5), for continuity with existing devlog entries.
        best_threshold: Threshold in the sweep grid that maximized
            mutation_macro_f1 for this checkpoint.
        f1_at_best: mutation_macro_f1 achieved at ``best_threshold``.
        all_scores: Mapping of every swept threshold to its
            mutation_macro_f1, for inspection/plotting.
    """

    f1_at_fixed: float
    best_threshold: float
    f1_at_best: float
    all_scores: dict[float, float]


@torch.no_grad()
def collect_logits(
    model: torch.nn.Module, loader: DataLoader, device: torch.device
) -> tuple[torch.Tensor, torch.Tensor]:
    """Run one full pass over ``loader`` and buffer flat logits + labels.

    Buffering once and re-scoring at many thresholds in memory is far
    cheaper than re-running the forward pass per threshold.

    Args:
        model: Trained model, will be set to eval mode.
        loader: Validation (or any) DataLoader yielding the same batch
            shape ``run_train.run_epoch`` consumes.
        device: Device to run inference on.

    Returns:
        ``(logits, labels)``: float tensor ``[N, num_classes]`` and long
        tensor ``[N]``, both on CPU.
    """
    model.eval()
    all_logits: list[torch.Tensor] = []
    all_labels: list[torch.Tensor] = []
    for batch in loader:
        batch = move_batch_to_device(batch, device)
        labels = batch["labels"].reshape(-1)
        model_inputs = filter_batch_for_model(model, batch)
        logits = model(**model_inputs)
        flat_logits = logits.reshape(-1, logits.size(-1))
        all_logits.append(flat_logits.detach().cpu())
        all_labels.append(labels.cpu())
    return torch.cat(all_logits), torch.cat(all_labels)


def mutation_macro_f1_at_threshold(
    logits: torch.Tensor, labels: torch.Tensor, threshold: float
) -> float:
    """mutation_macro_f1 (classes 1..N-1) at one decision threshold.

    Args:
        logits: Float tensor ``[N, num_classes]``.
        labels: Long tensor ``[N]`` of true class ids.
        threshold: Minimum softmax probability a non-background class
            needs to be predicted; see ``run_train.predict_with_threshold``.

    Returns:
        mutation_macro_f1, averaging F1 over classes 1..num_classes-1.
    """
    preds = predict_with_threshold(logits, threshold)
    confusion = torch_confusion_matrix(labels, preds, NUM_CLASSES)
    _, _, f1 = per_class_precision_recall_f1(confusion)
    return float(f1[1:].mean().item())


def sweep_thresholds(
    logits: torch.Tensor,
    labels: torch.Tensor,
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS,
    fixed_reference: float = 0.5,
) -> ThresholdSweepResult:
    """Find the threshold that maximizes mutation_macro_f1 for one checkpoint.

    Args:
        logits: Float tensor ``[N, num_classes]`` from one checkpoint's
            forward pass over the validation set (see :func:`collect_logits`).
        labels: Long tensor ``[N]`` of true class ids, same order as
            ``logits``.
        thresholds: Candidate thresholds to sweep.
        fixed_reference: Threshold to also report for continuity with
            existing fixed-threshold metrics (default matches
            ``run_train.VARIANT_PROB_THRESHOLD``).

    Returns:
        :class:`ThresholdSweepResult` for this checkpoint.
    """
    all_scores = {t: mutation_macro_f1_at_threshold(logits, labels, t) for t in thresholds}
    best_threshold = max(all_scores, key=all_scores.get)
    f1_at_fixed = mutation_macro_f1_at_threshold(logits, labels, fixed_reference)
    return ThresholdSweepResult(
        f1_at_fixed=f1_at_fixed,
        best_threshold=best_threshold,
        f1_at_best=all_scores[best_threshold],
        all_scores=all_scores,
    )
