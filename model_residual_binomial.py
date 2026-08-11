"""Residual Mamba variant caller initialized to reproduce the binomial caller.

    final_logits = binomial_prior_logits + residual(pileup_features)

The backbone is the *unmodified* :class:`~model_raw_pileup.RawPileupMamba`
stack over the same 14-channel raw-pileup representation, at the same width,
depth and dropout, so any difference against the earlier 14-channel experiment
is attributable to the residual parameterization and not to capacity.

Zero initialization
-------------------
The final classifier's weight *and* bias are set to zero, which makes
``residual(x) == 0`` identically for every input, at every position, for any
backbone state. Consequences:

* before a single gradient step the model's logits are exactly the binomial
  prior, so its ranking, AUC, AP and calls are *identical* -- not merely
  similar -- to the binomial baseline. ``sanity_residual_init.py`` verifies
  this numerically and refuses to proceed otherwise;
* the classifier is not dead: its gradient is ``delta @ final_norm(h)``, which
  is non-zero as soon as the prediction is wrong, so training escapes zero on
  the first step. Only the *output* is zeroed, never the backbone, whose
  weights keep their standard random initialization and keep producing a
  varied ``h`` for the classifier to read.

Zeroing the classifier is what makes this a *residual* model rather than an
independent classifier with a bias term. Randomly initializing that layer --
as the earlier experiments did -- puts an arbitrary O(1) perturbation on top
of a statistic whose useful range is ~[-20, 170], which is exactly the
"distorting the same statistic" failure this experiment is designed to avoid.
"""

from __future__ import annotations

import torch
import torch.nn as nn

from model_raw_pileup import RawPileupMamba
from raw_pileup import FEATURE_DIM
from residual_prior import NUM_CLASSES


class ResidualBinomialMamba(nn.Module):
    """Binomial prior plus a zero-initialized neural residual.

    Attributes:
        backbone: The unmodified 14-channel :class:`RawPileupMamba`, whose
            classifier layer is zeroed so it emits the residual.
    """

    def __init__(
        self,
        feature_dim: int = FEATURE_DIM,
        d_model: int = 128,
        num_layers: int = 6,
        num_classes: int = NUM_CLASSES,
        conv_kernel_size: int = 3,
        dropout: float = 0.15,
    ) -> None:
        """Build the model and zero the residual head.

        Args:
            feature_dim: Width of the per-locus evidence vector (14).
            d_model: Hidden width; kept at the previous experiments' 128.
            num_layers: Number of Mamba blocks; kept at 6.
            num_classes: Output classes; kept at 4 with unchanged meaning.
            conv_kernel_size: Passed through to each block.
            dropout: Passed through to each block and the input dropout.
        """
        super().__init__()
        self.backbone = RawPileupMamba(
            feature_dim=feature_dim, d_model=d_model, num_layers=num_layers,
            num_classes=num_classes, conv_kernel_size=conv_kernel_size, dropout=dropout,
        )
        self.reset_residual_head()

    def reset_residual_head(self) -> None:
        """Force ``residual(x) == 0`` for every input."""
        nn.init.zeros_(self.backbone.classifier.weight)
        nn.init.zeros_(self.backbone.classifier.bias)

    def residual(self, pileup_features: torch.Tensor) -> torch.Tensor:
        """The learned correction term.

        Args:
            pileup_features: Float tensor ``[B, L, feature_dim]``.

        Returns:
            Residual logits ``[B, L, num_classes]``.
        """
        return self.backbone(pileup_features)

    def forward(
        self,
        pileup_features: torch.Tensor,
        prior_logits: torch.Tensor,
        return_residual: bool = False,
    ) -> torch.Tensor | tuple[torch.Tensor, torch.Tensor]:
        """Add the neural residual to the binomial prior.

        Args:
            pileup_features: Float tensor ``[B, L, feature_dim]``.
            prior_logits: Float tensor ``[B, L, num_classes]`` from
                ``residual_prior.prior_logits_from_counts``.
            return_residual: Also return the raw residual, which the training
                loop needs for its ``lambda * mean(residual^2)`` penalty and
                the diagnostics need for their distribution analysis.

        Returns:
            ``final_logits`` or ``(final_logits, residual)``.
        """
        residual = self.residual(pileup_features)
        final_logits = prior_logits + residual
        if return_residual:
            return final_logits, residual
        return final_logits


__all__ = ["ResidualBinomialMamba"]
