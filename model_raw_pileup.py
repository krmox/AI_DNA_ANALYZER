"""Mamba variant caller adapted to a raw-pileup input instead of token ids.

Only the input projection changes. The backbone is the *same*
:class:`~dna_mamba_baseline_v5.MambaBlock` stack at the same width and depth,
and the classifier keeps the same 4-class semantics (Normal/SNP/Insertion/
Deletion), so a difference in results against the token model is attributable
to the input representation rather than to capacity.

The two ``nn.Embedding`` tables (observed + reference) and their fusion
projection are replaced by a single ``nn.Linear(FEATURE_DIM -> d_model)``,
because the input is now continuous evidence rather than a categorical base.
That makes existing checkpoints structurally incompatible, which is expected:
they were trained on a different input space. Old checkpoints are not touched
and remain loadable by the old model class.
"""

from __future__ import annotations

import torch
import torch.nn as nn

from dna_mamba_baseline_v5 import MambaBlock
from raw_pileup import FEATURE_DIM


class RawPileupMamba(nn.Module):
    """Per-locus variant classifier over raw pileup evidence.

    Attributes:
        input_proj: Projects the per-locus feature vector into model width.
        mamba_blocks: Unmodified bidirectional gated-conv blocks.
        classifier: Per-position 4-class head.
    """

    def __init__(
        self,
        feature_dim: int = FEATURE_DIM,
        d_model: int = 128,
        num_layers: int = 6,
        num_classes: int = 4,
        conv_kernel_size: int = 3,
        dropout: float = 0.15,
    ) -> None:
        """Build the model.

        Args:
            feature_dim: Width of the per-locus evidence vector.
            d_model: Hidden width; kept at the token model's 128.
            num_layers: Number of Mamba blocks; kept at the token model's 6.
            num_classes: Output classes; kept at 4 with unchanged meaning.
            conv_kernel_size: Passed through to each block.
            dropout: Passed through to each block and the input dropout.
        """
        super().__init__()

        # LayerNorm before the projection so that channels on different natural
        # scales (fractions in [0,1], normalised depth, one-hots) enter the
        # backbone comparably without hand-tuned per-channel weights.
        self.input_norm = nn.LayerNorm(feature_dim)
        self.input_proj = nn.Linear(feature_dim, d_model)
        self.embed_dropout = nn.Dropout(dropout)

        self.mamba_blocks = nn.ModuleList(
            [
                MambaBlock(d_model=d_model, conv_kernel_size=conv_kernel_size, dropout=dropout)
                for _ in range(num_layers)
            ]
        )

        self.final_norm = nn.LayerNorm(d_model)
        self.classifier = nn.Linear(d_model, num_classes)

    def forward(self, pileup_features: torch.Tensor) -> torch.Tensor:
        """Classify every locus in a batch of windows.

        Args:
            pileup_features: Float tensor ``[B, L, feature_dim]``.

        Returns:
            Logits ``[B, L, num_classes]``.
        """
        x = self.input_norm(pileup_features)
        x = self.input_proj(x)
        x = self.embed_dropout(x)

        for block in self.mamba_blocks:
            x = block(x)

        x = self.final_norm(x)
        return self.classifier(x)


__all__ = ["RawPileupMamba"]
