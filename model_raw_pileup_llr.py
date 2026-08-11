"""RawPileupMamba variant that consumes the 15-channel (pileup + LLR) input.

Deliberately a thin subclass: the backbone, width, depth, dropout and
classifier are inherited unchanged from :class:`~model_raw_pileup.RawPileupMamba`,
and only ``feature_dim`` differs (15 instead of 14). Keeping it a separate
class rather than passing ``feature_dim=15`` to the old one means the old
checkpoints stay loadable by the old class with no ambiguity about which
representation a file belongs to.
"""

from __future__ import annotations

from model_raw_pileup import RawPileupMamba
from llr_features import FEATURE_DIM_LLR


class RawPileupMambaLLR(RawPileupMamba):
    """Same architecture as :class:`RawPileupMamba`, 15-dimensional input."""

    def __init__(self, feature_dim: int = FEATURE_DIM_LLR, **kwargs) -> None:
        """Build the model.

        Args:
            feature_dim: Input width; 15 for pileup + LLR.
            **kwargs: Forwarded unchanged (``d_model``, ``num_layers``,
                ``num_classes``, ``conv_kernel_size``, ``dropout``).
        """
        super().__init__(feature_dim=feature_dim, **kwargs)


__all__ = ["RawPileupMambaLLR"]
