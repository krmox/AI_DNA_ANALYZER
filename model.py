"""Reference-guided bidirectional Mamba variant caller with quality/depth input.

Architecture naming, stated plainly: :class:`BiMambaBlock` is **not** a
selective state-space model. A true Mamba block uses input-dependent
discretisation of a continuous SSM (dynamic ``B``, ``C`` and step ``delta``)
together with a hardware-aware parallel scan supplied by the ``mamba-ssm``
and ``causal-conv1d`` CUDA extensions. What follows is a *causal gated
convolutional* approximation in stock PyTorch so the model runs anywhere,
CPU included. Its receptive field per layer is bounded by
``conv_kernel_size``, not unbounded state.

The stack is:

.. code-block:: text

    input_ids, reference_ids [B, L]      base_quality, depth [B, L]
      -> Embedding(obs) , Embedding(ref)   -> quality / depth features
      -> InputFusion: concat all active channels -> Linear -> [B, L, D]
      -> num_layers x BiMambaBlock                             [B, L, D]
      -> LayerNorm -> Linear(D -> num_classes)                 [B, L, C]

Why four channels rather than two. The base channels alone cannot separate
a sequencing error from a true SNP — both are a single observed base
disagreeing with the reference, so no function of ``(ref, obs)`` can tell
them apart. Phred quality is the signal that breaks the tie, because
miscalls concentrate at low confidence. Depth is included as specified,
though in a single-read representation it carries far less signal than
quality; ``stress_test.py`` ablates both rather than assuming.
"""

from __future__ import annotations

import torch
import torch.nn as nn

from config import ModelConfig


class CausalGatedConvCore(nn.Module):
    """Single-direction causal gated convolution over a sequence.

    Applies pre-norm, a depthwise causal convolution, post-convolution norm,
    SiLU, and an input-dependent gate mixing the convolved signal with the
    normalised input. Causality comes from left-padding by
    ``kernel_size - 1`` and discarding the overhang, so position ``t`` sees
    only positions ``<= t``.

    This core deliberately omits a residual connection; the surrounding
    :class:`BiMambaBlock` owns residual and normalisation for the merged
    forward and backward paths.
    """

    def __init__(self, d_model: int, conv_kernel_size: int, dropout: float) -> None:
        """Initialise the directional core.

        Args:
            d_model: Channel width, preserved end to end.
            conv_kernel_size: Depthwise convolution width.
            dropout: Dropout applied after the activation.
        """
        super().__init__()
        self.norm_pre = nn.LayerNorm(d_model)
        self.depthwise_conv = nn.Conv1d(
            in_channels=d_model,
            out_channels=d_model,
            kernel_size=conv_kernel_size,
            padding=conv_kernel_size - 1,
            groups=d_model,
        )
        self.norm_post_conv = nn.LayerNorm(d_model)
        self.activation = nn.SiLU()
        self.selective_gate = nn.Linear(d_model, d_model)
        self.dropout = nn.Dropout(dropout)

    def forward(self, hidden_states: torch.Tensor) -> torch.Tensor:
        """Run the causal gated convolution.

        Args:
            hidden_states: Float tensor ``[batch, seq_len, d_model]``.

        Returns:
            Float tensor ``[batch, seq_len, d_model]``.
        """
        seq_len = hidden_states.size(1)
        normed = self.norm_pre(hidden_states)

        # Conv1d expects [B, C, L]; trim the right overhang to stay causal.
        convolved = self.depthwise_conv(normed.transpose(1, 2))[:, :, :seq_len]
        convolved = convolved.transpose(1, 2)

        convolved = self.dropout(self.activation(self.norm_post_conv(convolved)))

        gate = torch.sigmoid(self.selective_gate(normed))
        return gate * convolved + (1.0 - gate) * normed


class BiMambaBlock(nn.Module):
    """Bidirectional causal gated convolutional block.

    Runs two independently parameterised causal cores: one over the sequence
    as given (left context only) and one over the reversed sequence whose
    output is flipped back (right context only). Concatenating and projecting
    gives every position access to both sides while each core stays causal.

    Weights are *not* shared between directions, so each can specialise.
    """

    def __init__(self, d_model: int, conv_kernel_size: int, dropout: float) -> None:
        """Initialise the bidirectional block.

        Args:
            d_model: Channel width, preserved end to end.
            conv_kernel_size: Depthwise convolution width for both cores.
            dropout: Dropout applied inside the cores and after fusion.
        """
        super().__init__()
        self.forward_core = CausalGatedConvCore(d_model, conv_kernel_size, dropout)
        self.backward_core = CausalGatedConvCore(d_model, conv_kernel_size, dropout)
        self.combine_proj = nn.Linear(2 * d_model, d_model)
        self.combine_dropout = nn.Dropout(dropout)
        self.out_norm = nn.LayerNorm(d_model)

    def forward(self, hidden_states: torch.Tensor) -> torch.Tensor:
        """Merge left-context and right-context representations.

        Args:
            hidden_states: Float tensor ``[batch, seq_len, d_model]``.

        Returns:
            Float tensor ``[batch, seq_len, d_model]``.
        """
        forward_states = self.forward_core(hidden_states)

        reversed_states = torch.flip(hidden_states, dims=[1])
        backward_states = torch.flip(self.backward_core(reversed_states), dims=[1])

        merged = self.combine_proj(torch.cat([forward_states, backward_states], dim=-1))
        return self.out_norm(hidden_states + self.combine_dropout(merged))


class InputFusion(nn.Module):
    """Embeds and fuses the reference, observed, quality and depth channels.

    Quality is discretised into learnable buckets rather than fed as a raw
    float. The relationship between Phred and reliability is non-linear —
    the gap between Q10 and Q20 matters far more than between Q40 and Q50 —
    and bucketing lets the model learn that shape instead of having a linear
    layer approximate it. Depth is passed as ``log1p`` because coverage is
    heavy-tailed and its informative variation is multiplicative.

    VAF is already a bounded fraction in ``[0, 1]``, so it needs no
    transform and is projected directly. It is the channel that most
    directly separates a true variant from a sequencing error: a germline
    heterozygous site sits near 0.5 and a homozygous one near 1.0, while an
    error sits near the per-base error rate. Deletion gaps count toward it,
    so all three variant classes land on a comparable scale.

    All three auxiliary channels are optional, controlled by
    :attr:`~.config.ModelConfig.use_base_quality`,
    :attr:`~.config.ModelConfig.use_depth` and
    :attr:`~.config.ModelConfig.use_vaf`, so an ablation is a config change
    rather than a code change.
    """

    def __init__(self, config: ModelConfig) -> None:
        """Build the fusion layer.

        Args:
            config: Validated :class:`~.config.ModelConfig`.
        """
        super().__init__()
        self.config = config

        self.observed_embedding = nn.Embedding(
            config.vocab_size, config.d_model, padding_idx=config.padding_idx
        )
        self.reference_embedding = nn.Embedding(
            config.vocab_size, config.d_model, padding_idx=config.padding_idx
        )

        # Width of one Phred bucket. Registered as a buffer so it moves with
        # the model and is captured in the state dict.
        self.register_buffer(
            "quality_bin_width",
            torch.tensor(config.max_phred / config.quality_bins, dtype=torch.float),
            persistent=False,
        )

        self.quality_embedding: nn.Embedding | None = None
        self.depth_proj: nn.Linear | None = None
        self.vaf_proj: nn.Linear | None = None
        channels = 2

        if config.use_base_quality:
            self.quality_embedding = nn.Embedding(config.quality_bins, config.d_model)
            channels += 1

        if config.use_depth:
            self.depth_proj = nn.Linear(1, config.d_model)
            channels += 1

        if config.use_vaf:
            self.vaf_proj = nn.Linear(1, config.d_model)
            channels += 1

        self.channels = channels
        self.fusion_proj: nn.Linear | None = (
            nn.Linear(channels * config.d_model, config.d_model)
            if config.fusion == "concat"
            else None
        )
        self.dropout = nn.Dropout(config.dropout)

    def _bucket_quality(self, base_quality: torch.Tensor) -> torch.Tensor:
        """Discretise Phred scores into embedding bucket indices.

        Args:
            base_quality: Float tensor ``[batch, seq_len]`` of Phred scores.

        Returns:
            Long tensor ``[batch, seq_len]`` of bucket indices.
        """
        buckets = (base_quality / self.quality_bin_width).long()
        return buckets.clamp(0, self.config.quality_bins - 1)

    def forward(
        self,
        input_ids: torch.Tensor,
        reference_ids: torch.Tensor,
        base_quality: torch.Tensor | None = None,
        depth: torch.Tensor | None = None,
        vaf: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Fuse all active channels into one hidden stream.

        Args:
            input_ids: Long tensor ``[batch, seq_len]`` of observed tokens.
            reference_ids: Long tensor ``[batch, seq_len]`` of reference tokens.
            base_quality: Float tensor ``[batch, seq_len]`` of Phred scores.
                Required when ``use_base_quality`` is set.
            depth: Float tensor ``[batch, seq_len]`` of read depth. Required
                when ``use_depth`` is set.
            vaf: Float tensor ``[batch, seq_len]`` of variant allele
                fractions in ``[0, 1]``. Required when ``use_vaf`` is set.

        Returns:
            Float tensor ``[batch, seq_len, d_model]``.

        Raises:
            ValueError: If a channel the config enables was not supplied.
        """
        features = [
            self.observed_embedding(input_ids),
            self.reference_embedding(reference_ids),
        ]

        if self.quality_embedding is not None:
            if base_quality is None:
                raise ValueError("base_quality is required when use_base_quality is enabled")
            features.append(self.quality_embedding(self._bucket_quality(base_quality)))

        if self.depth_proj is not None:
            if depth is None:
                raise ValueError("depth is required when use_depth is enabled")
            features.append(self.depth_proj(torch.log1p(depth).unsqueeze(-1)))

        if self.vaf_proj is not None:
            if vaf is None:
                raise ValueError("vaf is required when use_vaf is enabled")
            features.append(self.vaf_proj(vaf.unsqueeze(-1)))

        if self.fusion_proj is not None:
            fused = self.fusion_proj(torch.cat(features, dim=-1))
        else:
            fused = torch.stack(features, dim=0).sum(dim=0)

        return self.dropout(fused)


class VariantCaller(nn.Module):
    """Per-column variant classifier over gap-padded alignment channels.

    Consumes reference and observed token streams that are already
    coordinate-aligned (gap-padded CIGAR style), plus per-column Phred
    quality and read depth, and emits a class distribution per column.
    """

    def __init__(self, config: ModelConfig) -> None:
        """Build the model from an architecture config.

        Args:
            config: Validated :class:`~.config.ModelConfig`.
        """
        super().__init__()
        self.config = config

        self.input_fusion = InputFusion(config)

        self.blocks = nn.ModuleList(
            BiMambaBlock(config.d_model, config.conv_kernel_size, config.dropout)
            for _ in range(config.num_layers)
        )

        self.final_norm = nn.LayerNorm(config.d_model)
        self.classifier = nn.Linear(config.d_model, config.num_classes)

    def forward(
        self,
        input_ids: torch.Tensor,
        reference_ids: torch.Tensor,
        base_quality: torch.Tensor | None = None,
        depth: torch.Tensor | None = None,
        vaf: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Classify every alignment column.

        Args:
            input_ids: Long tensor ``[batch, seq_len]`` of observed tokens.
            reference_ids: Long tensor ``[batch, seq_len]`` of reference tokens.
            base_quality: Optional float tensor ``[batch, seq_len]`` of Phred
                scores.
            depth: Optional float tensor ``[batch, seq_len]`` of read depth.
            vaf: Optional float tensor ``[batch, seq_len]`` of variant allele
                fractions.

        Returns:
            Logit tensor ``[batch, seq_len, num_classes]``.

        Raises:
            ValueError: If the two token channels differ in shape, are not
                2-D, or an auxiliary channel's shape disagrees with them.
        """
        if input_ids.shape != reference_ids.shape:
            raise ValueError(
                f"input_ids {tuple(input_ids.shape)} and reference_ids "
                f"{tuple(reference_ids.shape)} must have identical shapes"
            )
        if input_ids.dim() != 2:
            raise ValueError(f"expected 2-D [batch, seq_len] inputs, got {input_ids.dim()}-D")

        for name, tensor in (
            ("base_quality", base_quality),
            ("depth", depth),
            ("vaf", vaf),
        ):
            if tensor is not None and tensor.shape != input_ids.shape:
                raise ValueError(
                    f"{name} {tuple(tensor.shape)} must match input_ids {tuple(input_ids.shape)}"
                )

        hidden_states = self.input_fusion(input_ids, reference_ids, base_quality, depth, vaf)

        for block in self.blocks:
            hidden_states = block(hidden_states)

        return self.classifier(self.final_norm(hidden_states))

    def count_parameters(self, trainable_only: bool = True) -> int:
        """Count model parameters.

        Args:
            trainable_only: Count only parameters requiring gradients.

        Returns:
            The parameter count.
        """
        return sum(
            parameter.numel()
            for parameter in self.parameters()
            if parameter.requires_grad or not trainable_only
        )


__all__ = ["CausalGatedConvCore", "BiMambaBlock", "InputFusion", "VariantCaller"]
