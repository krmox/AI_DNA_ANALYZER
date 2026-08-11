"""Feature assembly for the "raw pileup + binomial LLR" experiment.

Builds the model input from raw integer pileup counts:

* channels 0-13 -- byte-identical to ``raw_pileup.RawPileupProvider``'s 14
  features, recomputed here from counts so that a single pileup pass feeds
  both the neural features and the binomial caller;
* channel 14 -- the binomial log-likelihood ratio, transformed.

Reconstructing rather than re-extracting is verified, not assumed:
:func:`verify_against_raw_pileup_provider` compares this module's channels
0-13 against the original provider on real data and is run as a sanity check
before training.

Leakage boundary
----------------
Everything here is a function of ``counts``, which comes from
``pileup_counts.PileupCountsProvider`` and is built from reads and the
reference only. No VCF is opened, no label is read, and the candidate
alternate allele is chosen by ``BinomialVariantCaller`` from observed counts
alone. The LLR enters as a *continuous* value; no thresholded or binarised
form of it is ever exposed to the model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from binomial_baseline import BinomialVariantCaller

#: Index of the LLR channel in the 15-dimensional feature vector.
LLR_CHANNEL = 14
FEATURE_DIM_LLR = 15

# Column indices into the count matrix (see pileup_counts.COUNT_COLUMNS).
C_A, C_C, C_G, C_T, C_GAP, C_INS, C_DEPTH, C_QSUM, C_MSUM, C_REF = range(10)

DEPTH_SCALE = 5.0
BASE_QUALITY_SCALE = 40.0
MAPPING_QUALITY_SCALE = 60.0


def base_features_from_counts(counts: np.ndarray) -> np.ndarray:
    """Rebuild ``raw_pileup``'s 14 channels from raw counts.

    Args:
        counts: ``[N, 10]`` matrix following ``pileup_counts.COUNT_COLUMNS``.

    Returns:
        ``[N, 14]`` float array matching ``RawPileupProvider`` channel order.
    """
    counts = np.asarray(counts, dtype=np.float64)
    n = counts.shape[0]
    called = counts[:, C_A:C_T + 1].sum(axis=1)
    observed = called + counts[:, C_GAP]
    safe = np.maximum(observed, 1.0)

    features = np.zeros((n, 14), dtype=np.float64)
    for offset, column in enumerate((C_A, C_C, C_G, C_T)):
        features[:, offset] = np.where(observed > 0, counts[:, column] / safe, 0.0)
    features[:, 4] = np.where(observed > 0, counts[:, C_GAP] / safe, 0.0)
    features[:, 5] = np.where(observed > 0, counts[:, C_INS] / safe, 0.0)
    features[:, 6] = np.log1p(counts[:, C_DEPTH]) / DEPTH_SCALE
    features[:, 7] = np.where(
        called > 0, counts[:, C_QSUM] / np.maximum(called, 1.0) / BASE_QUALITY_SCALE, 0.0)
    features[:, 8] = np.where(
        counts[:, C_DEPTH] > 0,
        counts[:, C_MSUM] / np.maximum(counts[:, C_DEPTH], 1.0) / MAPPING_QUALITY_SCALE, 0.0)

    reference_index = counts[:, C_REF].astype(int)
    for slot in range(5):                       # 9..13 -> N, A, C, G, T
        features[:, 9 + slot] = (reference_index == slot).astype(np.float64)
    return features


def raw_llr_from_counts(counts: np.ndarray, caller: BinomialVariantCaller) -> np.ndarray:
    """Untransformed binomial LLR per locus.

    Args:
        counts: ``[N, 10]`` count matrix.
        caller: Configured :class:`BinomialVariantCaller`.

    Returns:
        ``[N]`` LLR values.
    """
    counts = np.asarray(counts, dtype=np.float64)
    score = caller.score_counts(
        counts[:, C_A:C_T + 1], counts[:, C_REF].astype(int),
        quality_sum=counts[:, C_QSUM] if caller.quality_derived_epsilon else None,
    )
    return score.llr


@dataclass
class LLRTransform:
    """Deterministic squash-and-standardize for the LLR channel.

    Raw LLRs span roughly -35 to +300 on this data, which would dominate the
    other channels (all in [0, 2]) at the input LayerNorm. The transform is
    ``sign(x) * log1p(|x|)`` followed by standardization; the signed log keeps
    the ordering and the sign of the evidence while compressing the long
    positive tail, and it is monotone, so no information about the ranking of
    loci is lost.

    ``mean`` and ``scale`` are estimated on TRAINING data only and then frozen;
    validation and test are transformed with the stored constants.

    Attributes:
        mean: Mean of the signed-log values on the training split.
        scale: Standard deviation of the same, floored away from zero.
    """

    mean: float
    scale: float

    @staticmethod
    def squash(llr: np.ndarray) -> np.ndarray:
        """Monotone signed-log compression of raw LLR values."""
        llr = np.asarray(llr, dtype=np.float64)
        return np.sign(llr) * np.log1p(np.abs(llr))

    @classmethod
    def fit(cls, train_llr: np.ndarray) -> "LLRTransform":
        """Estimate the transform from training LLRs only."""
        squashed = cls.squash(train_llr)
        return cls(mean=float(squashed.mean()), scale=float(max(squashed.std(), 1e-6)))

    def apply(self, llr: np.ndarray) -> np.ndarray:
        """Apply the frozen transform."""
        return (self.squash(llr) - self.mean) / self.scale

    def to_dict(self) -> dict:
        """Serializable form, recorded in every checkpoint."""
        return {"kind": "sign(x)*log1p(|x|) then standardize",
                "mean": self.mean, "scale": self.scale}


def assemble(counts: np.ndarray, caller: BinomialVariantCaller,
             transform: LLRTransform) -> np.ndarray:
    """Build the full ``[N, 15]`` feature matrix.

    Args:
        counts: ``[N, 10]`` count matrix.
        caller: Binomial caller supplying the LLR channel.
        transform: Frozen, train-fitted LLR transform.

    Returns:
        ``[N, 15]`` float array.
    """
    base = base_features_from_counts(counts)
    llr = transform.apply(raw_llr_from_counts(counts, caller))
    return np.concatenate([base, llr[:, None]], axis=1)


def verify_against_raw_pileup_provider(
    counts: np.ndarray, provider_features: np.ndarray, tolerance: float = 1e-6
) -> tuple[bool, float]:
    """Check that channels 0-13 reproduce the original provider exactly.

    The tolerance is 1e-6, not 0: ``RawPileupProvider`` materializes its
    features as ``torch.float`` (float32) while this module works in float64,
    so a round-trip differs by float32 epsilon (~1e-7) on the fraction
    channels. Anything above 1e-6 would indicate a genuine formula mismatch
    rather than storage precision.

    Args:
        counts: ``[N, 10]`` counts for the same loci, same order.
        provider_features: ``[N, 14]`` features from ``RawPileupProvider``.
        tolerance: Maximum permitted absolute difference.

    Returns:
        ``(passed, max_absolute_difference)``.
    """
    rebuilt = base_features_from_counts(counts)
    difference = float(np.abs(rebuilt - np.asarray(provider_features, np.float64)).max())
    return difference <= tolerance, difference


__all__ = [
    "LLR_CHANNEL", "FEATURE_DIM_LLR", "LLRTransform", "assemble",
    "base_features_from_counts", "raw_llr_from_counts",
    "verify_against_raw_pileup_provider",
]
