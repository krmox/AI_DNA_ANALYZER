"""Read-level variant caller: per-read encoder -> masked pooling -> Mamba.

    [B, L, R, F] reads
        -> per-read MLP encoder            (shared across reads and loci)
        -> masked attention + mean pooling (permutation-invariant)
        -> [B, L, d_model] locus embedding
        -> the *unmodified* MambaBlock stack
        -> per-locus 4-class logits

Design constraints, all of them there to keep the comparison honest:

* **The backbone is unchanged.** Same :class:`~dna_mamba_baseline_v5.MambaBlock`
  stack, same width 128, same depth 6, same dropout. Only the way a locus
  embedding is formed differs, so a difference against the 14-channel model is
  attributable to the input representation rather than to capacity.

* **Parameter count stays comparable.** 420,193 against the aggregate model's
  412,064 (+1.97%): the aggregate model's ``Linear(14 -> 128)`` input projection
  is replaced by a two-layer read encoder plus a scoring vector. The read
  encoder is shared across all ``R`` reads and all ``L`` loci, so widening the
  input from 14 numbers to 48x14 numbers costs almost nothing.

* **The aggregation respects the mask and is permutation-invariant.** Padded
  rows are excluded from both the softmax and the mean, so "how many reads" can
  only enter through the explicit depth term, never as a side effect of summing
  over padding. Invariance is proved by unit test rather than assumed, which
  also disposes of the read-order control: with an order-invariant aggregator,
  permuting whole read rows provably cannot change the output.

* **Depth is supplied explicitly.** ``log1p(n_valid) / 5`` is projected and
  added to the locus embedding. Without it the read-level model would be
  strictly *missing* the single most important input the binomial caller uses,
  and any loss would be uninterpretable. This mirrors channel 6 of the
  aggregate representation exactly.

Flattening ``[R, F]`` into one vector was rejected: it would make the model
order-dependent (the thing we are explicitly controlling for), scale parameters
with ``R``, and force a fixed depth. The pooling formulation keeps the model
agnostic to how many reads a locus happens to have.
"""

from __future__ import annotations

import torch
import torch.nn as nn

from dna_mamba_baseline_v5 import MambaBlock
from read_level_pileup import FEATURE_DIM, MAX_READS

#: Matches ``log1p(depth)/5`` in the aggregate representation.
DEPTH_SCALE = 5.0


class ReadLevelMamba(nn.Module):
    """Per-locus variant classifier over individual read observations.

    Attributes:
        read_encoder: Shared MLP applied to every read independently.
        attention_score: Produces the pooling logit for each read.
        depth_proj: Projects normalized depth into the locus embedding.
        mamba_blocks: Unmodified bidirectional gated-conv blocks.
        classifier: Per-position 4-class head.
    """

    def __init__(
        self,
        feature_dim: int = FEATURE_DIM,
        read_hidden: int = 64,
        d_model: int = 128,
        num_layers: int = 6,
        num_classes: int = 4,
        conv_kernel_size: int = 3,
        dropout: float = 0.15,
    ) -> None:
        """Build the model.

        Args:
            feature_dim: Per-read feature width.
            read_hidden: Hidden width of the shared read encoder.
            d_model: Locus embedding width; kept at the previous models' 128.
            num_layers: Number of Mamba blocks; kept at 6.
            num_classes: Output classes; kept at 4 with unchanged meaning.
            conv_kernel_size: Passed through to each block.
            dropout: Passed through to each block and the embedding dropout.
        """
        super().__init__()
        self.read_norm = nn.LayerNorm(feature_dim)
        self.read_encoder = nn.Sequential(
            nn.Linear(feature_dim, read_hidden),
            nn.GELU(),
            nn.LayerNorm(read_hidden),
            nn.Linear(read_hidden, d_model),
        )
        self.attention_score = nn.Linear(d_model, 1)
        self.depth_proj = nn.Linear(1, d_model)
        # The reference base, one-hot over (N, A, C, G, T). Supplied at the
        # locus level exactly as channels 9-13 of the aggregate representation
        # supply it, and for the same reason: without it "every read says G" is
        # uninterpretable, since whether that is a variant depends entirely on
        # what the reference says. It is reference-genome information, not
        # truth-derived, and every baseline in the comparison already has it.
        # No *per-read* discordance feature is derived from it.
        self.reference_proj = nn.Linear(5, d_model)

        self.locus_norm = nn.LayerNorm(d_model)
        self.embed_dropout = nn.Dropout(dropout)
        self.mamba_blocks = nn.ModuleList(
            [MambaBlock(d_model=d_model, conv_kernel_size=conv_kernel_size, dropout=dropout)
             for _ in range(num_layers)]
        )
        self.final_norm = nn.LayerNorm(d_model)
        self.classifier = nn.Linear(d_model, num_classes)

    def encode_locus(self, read_features: torch.Tensor, read_valid: torch.Tensor,
                     reference_onehot: torch.Tensor) -> torch.Tensor:
        """Pool ``[B, L, R, F]`` reads into ``[B, L, d_model]`` locus embeddings.

        Args:
            read_features: Float tensor ``[B, L, R, F]``.
            read_valid: Bool tensor ``[B, L, R]``; True for real reads.
            reference_onehot: Float tensor ``[B, L, 5]`` over (N, A, C, G, T).

        Returns:
            ``[B, L, d_model]`` locus embeddings.
        """
        mask = read_valid.unsqueeze(-1).to(read_features.dtype)
        encoded = self.read_encoder(self.read_norm(read_features)) * mask

        # Masked attention pooling. Padded rows get -inf before the softmax so
        # they receive exactly zero weight; a locus with no reads at all would
        # produce a uniform-over-nothing softmax, so its pooled value is
        # forced to zero explicitly rather than left as NaN.
        scores = self.attention_score(encoded).masked_fill(~read_valid.unsqueeze(-1), -1e9)
        weights = torch.softmax(scores, dim=-2) * mask
        attention_pool = (encoded * weights).sum(dim=-2)

        count = read_valid.sum(dim=-1, keepdim=True).to(read_features.dtype)
        mean_pool = encoded.sum(dim=-2) / count.clamp(min=1.0)

        has_reads = (count > 0).to(read_features.dtype)
        depth_term = self.depth_proj(torch.log1p(count) / DEPTH_SCALE)

        locus = ((attention_pool + mean_pool) * has_reads + depth_term
                 + self.reference_proj(reference_onehot))
        return self.locus_norm(locus)

    def forward(self, read_features: torch.Tensor, read_valid: torch.Tensor,
                reference_onehot: torch.Tensor) -> torch.Tensor:
        """Classify every locus in a batch of windows.

        Args:
            read_features: Float tensor ``[B, L, R, F]``.
            read_valid: Bool tensor ``[B, L, R]``.
            reference_onehot: Float tensor ``[B, L, 5]``.

        Returns:
            Logits ``[B, L, num_classes]``.
        """
        x = self.embed_dropout(
            self.encode_locus(read_features, read_valid, reference_onehot))
        for block in self.mamba_blocks:
            x = block(x)
        return self.classifier(self.final_norm(x))


__all__ = ["ReadLevelMamba", "MAX_READS"]
