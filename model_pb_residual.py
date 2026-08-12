"""Poisson-binomial prior plus a neural residual over count/quality evidence.

    final_logits = PB_prior_logits + [gate(x) *] residual(x)

Relationship to the previous residual experiment
------------------------------------------------
``model_residual_binomial.ResidualBinomialMamba`` used the *same* backbone,
width, depth and dropout over the 14-channel ``raw_pileup`` features with the
fixed-epsilon binomial LLR as prior, and learned essentially a constant shift.

Exactly two things change here, and they are separable:

1. **The representation** -- 44 channels from :mod:`locus_evidence`, which
   restore integer count semantics and split base quality, strand, read
   position and mapping quality by reference-versus-alternate support. This is
   Variant 1, and it is the only change when ``use_gate=False``.
2. **The residual parameterization** -- an optional scalar gate. This is
   Variant 2 and is off by default, so a Variant 1 run isolates representation
   as the single moved variable against the previous experiment.

The backbone is the unmodified :class:`~model_raw_pileup.RawPileupMamba` stack
at ``d_model=128``, ``num_layers=6``, so parameter count moves only by the
input projection's width (14 -> 44 inputs, +3,840 weights) and, when enabled,
the gate head. A gain therefore cannot be attributed to capacity.

Zero initialization, preserved
------------------------------
The residual head's weight and bias are zeroed, so ``residual(x) == 0``
identically for every input and every backbone state. Before a single gradient
step the model's SNP-minus-Normal logit gap is *exactly* the Poisson-binomial
LLR -- not approximately -- so the PB caller's frozen threshold transfers
without recalibration and any deviation from PB is something the model learned.
``sanity_pb_residual.py`` verifies this numerically and refuses to proceed
otherwise.

The head is not dead at zero: its gradient is ``delta @ final_norm(h)``, which
is non-zero as soon as a prediction is wrong.

The gate
--------
``gate(x) = sigmoid(g(h))``, with ``g`` zero-initialized so the gate starts at
exactly 0.5 everywhere -- open enough to pass gradient, and *uninformative*, so
it carries no prior opinion about where the statistical caller is reliable. The
gate multiplies a residual that is identically zero at initialization, so the
initialization guarantee above is unaffected by enabling it.

The scientific point of the gate is that a constant shift is no longer the
cheapest way to reduce loss: the model can instead modulate *where* its
correction applies, which is the behaviour the previous experiment failed to
produce. Whether it does so is measured, not assumed -- the training loop logs
the gate's variance across loci.
"""

from __future__ import annotations

import torch
import torch.nn as nn

from locus_evidence import FEATURE_DIM
from model_raw_pileup import RawPileupMamba
from residual_prior import NUM_CLASSES


class PoissonBinomialResidualMamba(nn.Module):
    """PB prior plus a zero-initialized, optionally gated neural residual.

    Attributes:
        backbone: The unmodified ``RawPileupMamba`` stack, its classifier zeroed
            so it emits the residual.
        gate_head: Optional scalar gate over the backbone's hidden state.
    """

    def __init__(
        self,
        feature_dim: int = FEATURE_DIM,
        d_model: int = 128,
        num_layers: int = 6,
        num_classes: int = NUM_CLASSES,
        conv_kernel_size: int = 3,
        dropout: float = 0.15,
        use_gate: bool = False,
    ) -> None:
        """Build the model and zero the residual head.

        Args:
            feature_dim: Width of the per-locus evidence vector (44).
            d_model: Hidden width; kept at the previous experiments' 128.
            num_layers: Number of Mamba blocks; kept at 6.
            num_classes: Output classes; kept at 4 with unchanged meaning.
            conv_kernel_size: Passed through to each block.
            dropout: Passed through to each block and the input dropout.
            use_gate: Enable the Variant 2 scalar gate.
        """
        super().__init__()
        self.backbone = RawPileupMamba(
            feature_dim=feature_dim, d_model=d_model, num_layers=num_layers,
            num_classes=num_classes, conv_kernel_size=conv_kernel_size, dropout=dropout,
        )
        self.use_gate = use_gate
        self.gate_head = nn.Linear(d_model, 1) if use_gate else None
        self.reset_residual_head()

    def reset_residual_head(self) -> None:
        """Force ``residual(x) == 0`` and, if gated, ``gate(x) == 0.5``."""
        nn.init.zeros_(self.backbone.classifier.weight)
        nn.init.zeros_(self.backbone.classifier.bias)
        if self.gate_head is not None:
            nn.init.zeros_(self.gate_head.weight)
            nn.init.zeros_(self.gate_head.bias)

    def _hidden(self, pileup_features: torch.Tensor) -> torch.Tensor:
        """Backbone hidden state after the final norm, ``[B, L, d_model]``."""
        backbone = self.backbone
        x = backbone.input_norm(pileup_features)
        x = backbone.input_proj(x)
        x = backbone.embed_dropout(x)
        for block in backbone.mamba_blocks:
            x = block(x)
        return backbone.final_norm(x)

    def forward(
        self,
        pileup_features: torch.Tensor,
        prior_logits: torch.Tensor,
        return_parts: bool = False,
    ):
        """Add the (optionally gated) neural residual to the PB prior.

        Args:
            pileup_features: Float tensor ``[B, L, feature_dim]``.
            prior_logits: Float tensor ``[B, L, num_classes]``, the PB prior
                from ``residual_prior.prior_logits_from_llr``.
            return_parts: Also return the raw residual and the gate, which the
                training loop needs for its penalty and its diagnostics.

        Returns:
            ``final_logits``, or ``(final_logits, residual, gate)``.
        """
        hidden = self._hidden(pileup_features)
        residual = self.backbone.classifier(hidden)

        if self.gate_head is not None:
            gate = torch.sigmoid(self.gate_head(hidden))
            applied = gate * residual
        else:
            gate = torch.ones_like(residual[..., :1])
            applied = residual

        final_logits = prior_logits + applied
        if return_parts:
            return final_logits, residual, gate
        return final_logits


__all__ = ["PoissonBinomialResidualMamba"]
