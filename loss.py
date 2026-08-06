"""Focal loss and calibrated class-weight estimation.

Roughly 95% of tokens carry the ``Normal`` label, so an uncalibrated
objective either collapses to predicting ``Normal`` everywhere or, if the
weights over-correct, floods the output with false positives. Two
mechanisms share the work here:

* **Focal modulation** — ``(1 - p_t) ** gamma`` shrinks the contribution of
  tokens the model already classifies confidently, concentrating gradient
  on hard positions without inflating any class weight.
* **Calibrated alpha weights** — inverse frequency passed through a square
  root, renormalised to mean 1, then hard-clamped. The square root matters:
  raw inverse frequency produced roughly 16x weighting on mutation classes
  against 0.26x on ``Normal``, a ~60x asymmetry that made guessing
  mutations everywhere the cheaper strategy.

Label smoothing is supported but defaults to off. It caps attainable
confidence on exactly the rare classes that need it most.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset

from .config import LABEL_NAMES


class FocalLoss(nn.Module):
    """Multi-class focal loss with per-class alpha weights.

    Reduces to weighted cross-entropy when ``gamma`` is 0.
    """

    def __init__(
        self,
        alpha: torch.Tensor | None = None,
        gamma: float = 2.0,
        label_smoothing: float = 0.0,
        reduction: str = "mean",
    ) -> None:
        """Initialise the criterion.

        Args:
            alpha: Optional float tensor ``[num_classes]`` of per-class
                weights, registered as a buffer so it follows ``.to(device)``.
            gamma: Focusing exponent; must be non-negative.
            label_smoothing: Target smoothing in ``[0, 1)``.
            reduction: One of ``"mean"``, ``"sum"`` or ``"none"``.

        Raises:
            ValueError: If ``gamma``, ``label_smoothing`` or ``reduction`` is
                outside its permitted range.
        """
        super().__init__()
        if gamma < 0.0:
            raise ValueError("gamma must be non-negative")
        if not 0.0 <= label_smoothing < 1.0:
            raise ValueError("label_smoothing must lie in [0, 1)")
        if reduction not in ("mean", "sum", "none"):
            raise ValueError("reduction must be 'mean', 'sum' or 'none'")

        self.gamma = gamma
        self.label_smoothing = label_smoothing
        self.reduction = reduction

        if alpha is None:
            self.register_buffer("alpha", None)
        else:
            self.register_buffer("alpha", alpha.detach().clone().float())

    def forward(self, logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        """Compute the focal loss.

        Args:
            logits: Float tensor ``[num_tokens, num_classes]`` of raw logits.
            targets: Long tensor ``[num_tokens]`` of class ids.

        Returns:
            Scalar loss, or ``[num_tokens]`` when ``reduction`` is ``"none"``.

        Raises:
            ValueError: If ``logits`` is not 2-D or the batch dimensions
                disagree.
        """
        if logits.dim() != 2:
            raise ValueError(f"expected 2-D logits [num_tokens, num_classes], got {logits.dim()}-D")
        if logits.size(0) != targets.size(0):
            raise ValueError(
                f"logits ({logits.size(0)}) and targets ({targets.size(0)}) disagree on token count"
            )

        num_classes = logits.size(-1)
        log_probs = F.log_softmax(logits, dim=-1)

        if self.label_smoothing > 0.0:
            with torch.no_grad():
                smooth_targets = torch.full_like(
                    log_probs, self.label_smoothing / (num_classes - 1)
                )
                smooth_targets.scatter_(1, targets.unsqueeze(1), 1.0 - self.label_smoothing)
            cross_entropy = -(smooth_targets * log_probs).sum(dim=1)
        else:
            cross_entropy = F.nll_loss(log_probs, targets, reduction="none")

        # Focal factor keyed on the true (hard) class probability.
        true_class_prob = log_probs.exp().gather(1, targets.unsqueeze(1)).squeeze(1)
        focal_factor = (1.0 - true_class_prob.clamp(min=1e-6, max=1.0)) ** self.gamma

        loss = focal_factor * cross_entropy
        if self.alpha is not None:
            loss = self.alpha[targets] * loss

        if self.reduction == "mean":
            return loss.mean()
        if self.reduction == "sum":
            return loss.sum()
        return loss


def estimate_class_counts(
    dataset: Dataset,
    num_classes: int = len(LABEL_NAMES),
    sample_size: int = 200,
) -> torch.Tensor:
    """Tally label frequencies over a prefix of the dataset.

    Args:
        dataset: Any dataset whose items expose a ``"labels"`` tensor.
        num_classes: Number of classes to tally.
        sample_size: Maximum number of samples to inspect.

    Returns:
        Float tensor ``[num_classes]`` of raw counts, floored at 1 so that
        an unobserved class cannot produce a division by zero downstream.

    Raises:
        ValueError: If ``num_classes`` or ``sample_size`` is non-positive.
    """
    if num_classes <= 0:
        raise ValueError("num_classes must be positive")
    if sample_size <= 0:
        raise ValueError("sample_size must be positive")

    counts = torch.zeros(num_classes, dtype=torch.float)
    inspected = min(sample_size, len(dataset))  # type: ignore[arg-type]

    for index in range(inspected):
        labels = dataset[index]["labels"]
        counts += torch.bincount(labels.reshape(-1), minlength=num_classes)[:num_classes].float()

    return counts.clamp(min=1.0)


def weights_from_counts(
    counts: torch.Tensor,
    alpha: float = 0.5,
    max_weight_cap: float = 10.0,
    normalise: bool = True,
) -> torch.Tensor:
    """Derive alpha weights from class counts, given the counts directly.

    Computes ``w_c = min((N_total / N_c) ** alpha, max_weight_cap)``, with
    optional normalisation to mean 1 before the clamp.

    Separated from :func:`compute_calibrated_class_weights` so real-data
    counts produced by ``summarise_label_distribution.py`` can be fed in
    without re-scanning a dataset — scanning a real chr21 BAM is minutes of
    work, not milliseconds.

    On the exponent. At real genomic density (roughly 1:1000) raw inverse
    frequency asks for a weight near 1000x on the variant classes, which
    makes predicting variants everywhere the cheaper strategy and reproduces
    the false-positive flood seen early in this project. ``alpha=0.5``
    reduces that to about 32x, and the cap then brings it inside a range the
    optimiser can work with. The cap is doing real work at this density and
    is not a formality: without it, the sqrt alone is still too aggressive.

    Args:
        counts: Float tensor ``[num_classes]`` of per-class token counts.
        alpha: Smoothing exponent. ``0.5`` is sqrt weighting; ``1.0`` is raw
            inverse frequency; ``0.0`` disables weighting.
        max_weight_cap: Upper clamp on the resulting weights.
        normalise: Whether to rescale to mean 1 before clamping. Keeps the
            overall loss magnitude comparable across datasets.

    Returns:
        Float tensor ``[num_classes]`` of per-class weights.

    Raises:
        ValueError: If ``alpha`` is negative or ``max_weight_cap`` is not
            positive.
    """
    if alpha < 0.0:
        raise ValueError("alpha must be non-negative")
    if max_weight_cap <= 0.0:
        raise ValueError("max_weight_cap must be positive")

    counts = counts.float().clamp(min=1.0)
    raw_weights = (counts.sum() / counts) ** alpha

    if normalise:
        raw_weights = raw_weights / raw_weights.mean()

    return raw_weights.clamp(max=max_weight_cap)


def compute_calibrated_class_weights(
    dataset: Dataset,
    num_classes: int = len(LABEL_NAMES),
    sample_size: int = 200,
    max_weight_cap: float = 10.0,
    alpha: float = 0.5,
) -> torch.Tensor:
    """Estimate class frequencies from a dataset and derive alpha weights.

    Thin wrapper over :func:`weights_from_counts`; see that function for the
    reasoning behind the exponent and the cap.

    Args:
        dataset: Dataset to estimate frequencies from.
        num_classes: Number of classes.
        sample_size: Maximum number of samples to inspect.
        max_weight_cap: Upper clamp applied after normalisation.
        alpha: Smoothing exponent passed through.

    Returns:
        Float tensor ``[num_classes]`` of per-class weights.
    """
    counts = estimate_class_counts(dataset, num_classes, sample_size)
    return weights_from_counts(counts, alpha=alpha, max_weight_cap=max_weight_cap)


__all__ = [
    "FocalLoss",
    "estimate_class_counts",
    "weights_from_counts",
    "compute_calibrated_class_weights",
]
