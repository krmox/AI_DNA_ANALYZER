"""Reference-guided dual-channel bidirectional Mamba variant caller.

Architecture naming, stated plainly: :class:`BiMambaBlock` is **not** a
selective state-space model. A true Mamba block uses input-dependent
discretisation of a continuous SSM (dynamic ``B``, ``C`` and step ``delta``)
together with a hardware-aware parallel scan supplied by the ``mamba-ssm``
and ``causal-conv1d`` CUDA extensions. What follows is a *causal gated
convolutional* approximation written in stock PyTorch so the proof of
concept runs anywhere, CPU included. It borrows Mamba's shape — depthwise
causal convolution plus an input-dependent gate — but its receptive field
per layer is bounded by ``conv_kernel_size``, not unbounded state.

The stack is:

.. code-block:: text

    input_ids [B, L], reference_ids [B, L]
      -> Embedding(observed) + Embedding(reference)      two separate tables
      -> fusion (concat -> Linear, or elementwise sum)   [B, L, D]
      -> num_layers x BiMambaBlock                       [B, L, D]
      -> LayerNorm -> Linear(D -> num_classes)           [B, L, C]

Two embedding tables rather than one shared table: the same base means
something different depending on which channel it appears in, and the
comparison the model must learn is easier when the two roles start from
independent parameters.
"""

from __future__ import annotations

import torch
import torch.nn as nn

from .config import ModelConfig


class CausalGatedConvCore(nn.Module):
    """Single-direction causal gated convolution over a sequence.

    Applies pre-norm, a depthwise causal convolution, post-convolution norm,
    SiLU, and an input-dependent gate that mixes the convolved signal with
    the normalised input. Causality comes from left-padding by
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

        # Input-dependent gate: the selective-update stand-in.
        gate = torch.sigmoid(self.selective_gate(normed))
        return gate * convolved + (1.0 - gate) * normed


class BiMambaBlock(nn.Module):
    """Bidirectional causal gated convolutional block.

    Runs two independently parameterised causal cores: one over the sequence
    as given (left context only) and one over the reversed sequence, whose
    output is flipped back (right context only). Concatenating and projecting
    the two gives every position genuine access to both sides while each core
    stays strictly causal.

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


class VariantCaller(nn.Module):
    """Per-token variant classifier over aligned reference/observed channels.

    Consumes two index-aligned token streams and emits a class distribution
    for every position. Supplying the reference is not optional decoration:
    without it the substitution class is information-theoretically
    indistinguishable from background, since both are uniform random bases.
    """

    def __init__(self, config: ModelConfig) -> None:
        """Build the model from an architecture config.

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

        self.fusion_proj: nn.Linear | None = (
            nn.Linear(2 * config.d_model, config.d_model) if config.fusion == "concat" else None
        )
        self.embed_dropout = nn.Dropout(config.dropout)

        self.blocks = nn.ModuleList(
            BiMambaBlock(config.d_model, config.conv_kernel_size, config.dropout)
            for _ in range(config.num_layers)
        )

        self.final_norm = nn.LayerNorm(config.d_model)
        self.classifier = nn.Linear(config.d_model, config.num_classes)

    def fuse_embeddings(self, input_ids: torch.Tensor, reference_ids: torch.Tensor) -> torch.Tensor:
        """Embed both channels and combine them into one hidden stream.

        Args:
            input_ids: Long tensor ``[batch, seq_len]`` of observed tokens.
            reference_ids: Long tensor ``[batch, seq_len]`` of reference tokens.

        Returns:
            Float tensor ``[batch, seq_len, d_model]``.
        """
        observed = self.observed_embedding(input_ids)
        reference = self.reference_embedding(reference_ids)

        if self.fusion_proj is not None:
            return self.fusion_proj(torch.cat([observed, reference], dim=-1))
        return observed + reference

    def forward(self, input_ids: torch.Tensor, reference_ids: torch.Tensor) -> torch.Tensor:
        """Classify every aligned position.

        Args:
            input_ids: Long tensor ``[batch, seq_len]`` of observed tokens.
            reference_ids: Long tensor ``[batch, seq_len]`` of reference tokens.

        Returns:
            Logit tensor ``[batch, seq_len, num_classes]``.

        Raises:
            ValueError: If the two channels differ in shape or are not 2-D.
        """
        if input_ids.shape != reference_ids.shape:
            raise ValueError(
                f"input_ids {tuple(input_ids.shape)} and reference_ids "
                f"{tuple(reference_ids.shape)} must have identical shapes"
            )
        if input_ids.dim() != 2:
            raise ValueError(f"expected 2-D [batch, seq_len] inputs, got {input_ids.dim()}-D")

        hidden_states = self.embed_dropout(self.fuse_embeddings(input_ids, reference_ids))

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


__all__ = ["CausalGatedConvCore", "BiMambaBlock", "VariantCaller"]
