"""Depth-aware binomial SNP caller.

The flat allele-fraction rule treats ``k/n = 0.4`` identically at ``n=5`` and
``n=50``, but those observations carry very different evidence: 2-of-5 is an
unremarkable draw under a 1% error model at low depth, while 20-of-50 is
overwhelming. This module scores each locus by a likelihood ratio in which
``n`` and ``k`` enter explicitly.

Hypotheses
----------
Let ``eps`` be the per-base sequencing error probability at the locus and
``k`` the count of the most-supported non-reference base out of ``n`` usable
A/C/G/T observations.

* **H0 (noise)** -- no alternate allele; the alternate base arises only from
  error. A miscall lands on one specific wrong base, so ``p0 = eps / 3``,
  floored to keep the log finite.
* **H1_het** -- heterozygous: ``p = 0.5(1 - eps) + 0.5(eps / 3)``. The second
  term is the chance that a reference read is miscalled *into* the alternate
  base; it is tiny but keeps the model consistent.
* **H1_hom** -- homozygous alternate: ``p = 1 - eps + eps / 3``, i.e. almost
  every read carries the alternate.

``score = max(logpmf(k; n, p_het), logpmf(k; n, p_hom)) - logpmf(k; n, p0)``

Taking the max over the two alternate hypotheses matters here: the GIAB test
region is 472/660 homozygous-alt, and scoring those under a het-only model
would penalise exactly the loci with the *strongest* evidence. Reporting a
het-only variant is supported via ``include_homozygous=False`` so the two can
be compared rather than conflated.

Error model
-----------
``epsilon`` is either a fixed configurable floor (v1) or derived from the
locus' mean base quality as ``10 ** (-Q / 10)`` clamped to
``[error_floor, 0.25]`` (v2). The quality-derived form uses only the Phred
scores htslib reports for the reads -- no truth labels, no VCF, no model
output. Both variants are kept so v2 can be evaluated against v1 rather than
silently replacing it.

Everything is vectorized over loci with ``scipy.stats.binom.logpmf``; no
factorials are formed and no Python loop runs over loci.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import binom

#: Smallest error probability allowed, so ``log p0`` stays finite even when a
#: locus reports perfect base quality.
DEFAULT_ERROR_FLOOR = 1e-3

#: Fixed error rate for the v1 caller: ~Q20, a conservative Illumina figure.
DEFAULT_ERROR_RATE = 0.01


@dataclass
class LocusScore:
    """Per-locus diagnostics from the caller.

    Attributes:
        n: Usable A/C/G/T observations.
        k: Observations supporting the candidate alternate base.
        alt_index: 0-3 index (A,C,G,T) of the candidate alternate base.
        vaf: ``k / n``, 0 where ``n == 0``.
        epsilon: Error probability used at this locus.
        loglik_h0: ``log P(k | n, H0)``.
        loglik_het: ``log P(k | n, H1_het)``.
        loglik_hom: ``log P(k | n, H1_hom)``.
        llr: ``max(loglik_het, loglik_hom) - loglik_h0``.
    """

    n: np.ndarray
    k: np.ndarray
    alt_index: np.ndarray
    vaf: np.ndarray
    epsilon: np.ndarray
    loglik_h0: np.ndarray
    loglik_het: np.ndarray
    loglik_hom: np.ndarray
    llr: np.ndarray


class BinomialVariantCaller:
    """Likelihood-ratio SNP caller over raw pileup counts."""

    def __init__(
        self,
        error_rate: float = DEFAULT_ERROR_RATE,
        include_homozygous: bool = True,
        quality_derived_epsilon: bool = False,
        error_floor: float = DEFAULT_ERROR_FLOOR,
    ) -> None:
        """Configure the caller.

        Args:
            error_rate: Fixed per-base error probability, used when
                ``quality_derived_epsilon`` is False.
            include_homozygous: Score ``max(H1_het, H1_hom)`` rather than
                ``H1_het`` alone.
            quality_derived_epsilon: Derive epsilon per locus from mean base
                quality instead of using ``error_rate``.
            error_floor: Lower clamp on epsilon.

        Raises:
            ValueError: If a probability argument is outside (0, 1).
        """
        if not 0.0 < error_rate < 1.0:
            raise ValueError("error_rate must lie in (0, 1)")
        if not 0.0 < error_floor < 1.0:
            raise ValueError("error_floor must lie in (0, 1)")
        self.error_rate = error_rate
        self.include_homozygous = include_homozygous
        self.quality_derived_epsilon = quality_derived_epsilon
        self.error_floor = error_floor

    # -- core statistics ---------------------------------------------------

    def _epsilon(self, n: np.ndarray, quality_sum: np.ndarray | None) -> np.ndarray:
        """Per-locus error probability.

        Args:
            n: Usable observation counts.
            quality_sum: Summed Phred scores of those observations, or None.

        Returns:
            Array of epsilon values, same shape as ``n``.
        """
        if not self.quality_derived_epsilon or quality_sum is None:
            return np.full(n.shape, self.error_rate, dtype=np.float64)

        with np.errstate(invalid="ignore", divide="ignore"):
            mean_quality = np.where(n > 0, quality_sum / np.maximum(n, 1), 0.0)
        epsilon = np.power(10.0, -mean_quality / 10.0)
        # Loci with no reads get the fixed rate rather than 10**0 == 1.
        epsilon = np.where(n > 0, epsilon, self.error_rate)
        return np.clip(epsilon, self.error_floor, 0.25)

    def score_counts(
        self,
        base_counts: np.ndarray,
        reference_index: np.ndarray,
        quality_sum: np.ndarray | None = None,
    ) -> LocusScore:
        """Score every locus from raw A/C/G/T counts.

        The candidate alternate allele is the most-supported non-reference
        nucleotide, chosen from the counts alone. The truth VCF plays no part.

        Args:
            base_counts: ``[N, 4]`` integer counts in A,C,G,T order.
            reference_index: ``[N]`` with 0=N, 1=A, 2=C, 3=G, 4=T.
            quality_sum: Optional ``[N]`` summed Phred of counted bases.

        Returns:
            A :class:`LocusScore` with per-locus diagnostics.
        """
        base_counts = np.asarray(base_counts, dtype=np.float64)
        n = base_counts.sum(axis=1)

        # Mask the reference column so argmax cannot select it. Loci whose
        # reference base is N keep all four candidates, which is the
        # conservative reading of "no reference base to exclude".
        masked = base_counts.copy()
        has_reference = (reference_index >= 1) & (reference_index <= 4)
        rows = np.nonzero(has_reference)[0]
        masked[rows, (reference_index[rows] - 1).astype(int)] = -1.0

        alt_index = masked.argmax(axis=1)
        k = base_counts[np.arange(base_counts.shape[0]), alt_index]
        k = np.where(masked.max(axis=1) < 0, 0.0, k)

        epsilon = self._epsilon(n, quality_sum)

        p0 = np.clip(epsilon / 3.0, 1e-12, 1 - 1e-12)
        p_het = np.clip(0.5 * (1.0 - epsilon) + 0.5 * (epsilon / 3.0), 1e-12, 1 - 1e-12)
        p_hom = np.clip(1.0 - epsilon + epsilon / 3.0, 1e-12, 1 - 1e-12)

        loglik_h0 = binom.logpmf(k, n, p0)
        loglik_het = binom.logpmf(k, n, p_het)
        loglik_hom = binom.logpmf(k, n, p_hom)

        alternative = np.maximum(loglik_het, loglik_hom) if self.include_homozygous \
            else loglik_het
        llr = alternative - loglik_h0

        # n == 0 gives logpmf(0, 0, p) == 0 for every hypothesis, so the LLR is
        # already 0 -- no evidence either way. Force it explicitly so a NaN can
        # never leak from an unusual scipy edge case.
        no_data = n <= 0
        llr = np.where(no_data, 0.0, llr)
        llr = np.nan_to_num(llr, nan=0.0, posinf=0.0, neginf=0.0)

        with np.errstate(invalid="ignore", divide="ignore"):
            vaf = np.where(n > 0, k / np.maximum(n, 1), 0.0)

        return LocusScore(
            n=n, k=k, alt_index=alt_index, vaf=vaf, epsilon=epsilon,
            loglik_h0=loglik_h0, loglik_het=loglik_het, loglik_hom=loglik_hom, llr=llr,
        )

    def call_counts(
        self,
        base_counts: np.ndarray,
        reference_index: np.ndarray,
        threshold: float,
        quality_sum: np.ndarray | None = None,
    ) -> np.ndarray:
        """Boolean SNP call per locus at a fixed LLR threshold."""
        return self.score_counts(base_counts, reference_index, quality_sum).llr >= threshold

    # -- single-locus convenience (used by the unit tests) -----------------

    def score_locus(self, counts: dict[str, int], reference_base: str,
                    quality_sum: float | None = None) -> LocusScore:
        """Score one locus from a ``{'A': .., 'C': .., ...}`` count dict.

        Args:
            counts: Per-base observation counts.
            reference_base: Reference nucleotide.
            quality_sum: Optional summed Phred of the counted bases.

        Returns:
            A :class:`LocusScore` of length-1 arrays.
        """
        order = ("A", "C", "G", "T")
        tokens = ("N",) + order
        matrix = np.array([[counts.get(base, 0) for base in order]], dtype=np.float64)
        base = reference_base.upper()
        index = np.array([tokens.index(base) if base in tokens else 0])
        return self.score_counts(
            matrix, index, None if quality_sum is None else np.array([quality_sum], float)
        )


def binomial_tail_pvalue(k: np.ndarray, n: np.ndarray, epsilon: np.ndarray) -> np.ndarray:
    """One-sided P(K >= k | n, p0) under the noise hypothesis.

    Provided as the multiple-testing-friendly companion to the LLR; it is a
    reported diagnostic, not the primary decision score.

    Args:
        k: Alternate-supporting counts.
        n: Usable observations.
        epsilon: Per-locus error probability.

    Returns:
        Upper-tail probabilities in [0, 1].
    """
    p0 = np.clip(epsilon / 3.0, 1e-12, 1 - 1e-12)
    return np.where(n > 0, binom.sf(k - 1, n, p0), 1.0)


def benjamini_hochberg(pvalues: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg rejection mask at level ``alpha``.

    Args:
        pvalues: Array of p-values.
        alpha: Target false-discovery rate.

    Returns:
        Boolean array, True where the null is rejected.
    """
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    m = ranked.size
    thresholds = alpha * np.arange(1, m + 1) / m
    passing = np.nonzero(ranked <= thresholds)[0]
    rejected = np.zeros(m, dtype=bool)
    if passing.size:
        rejected[order[: passing[-1] + 1]] = True
    return rejected


__all__ = [
    "BinomialVariantCaller",
    "LocusScore",
    "binomial_tail_pvalue",
    "benjamini_hochberg",
    "DEFAULT_ERROR_RATE",
    "DEFAULT_ERROR_FLOOR",
]
