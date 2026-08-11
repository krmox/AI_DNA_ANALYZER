"""Unit tests for the depth-aware binomial caller.

These test the statistics in isolation -- no BAM, no VCF, no model. Every
expectation is a property the likelihood ratio must satisfy for the
downstream benchmark to be interpretable.
"""

from __future__ import annotations

import numpy as np
import pytest

from binomial_baseline import (
    BinomialVariantCaller,
    benjamini_hochberg,
    binomial_tail_pvalue,
)


def llr_for(n: int, k: int, reference: str = "A", alt: str = "C", **kwargs) -> float:
    """LLR for a locus with ``k`` alternate observations out of ``n``."""
    caller = BinomialVariantCaller(**kwargs)
    counts = {reference: n - k, alt: k}
    return float(caller.score_locus(counts, reference).llr[0])


class TestHypothesisPreference:
    def test_zero_alternate_reads_favours_noise(self) -> None:
        """Test 1: n=20, k=0 must strongly favour H0."""
        assert llr_for(20, 0) < -5.0

    def test_heterozygous_like_favours_h1(self) -> None:
        """Test 2: n=20, k=10 must strongly favour H1_het."""
        assert llr_for(20, 10) > 20.0

    def test_high_depth_error_favours_noise(self) -> None:
        """Test 3: n=100, k=2 at a realistic error rate should favour H0."""
        assert llr_for(100, 2) < 0.0

    def test_homozygous_like_prefers_hom_over_het(self) -> None:
        """Test 5: at k==n the homozygous hypothesis must dominate."""
        caller = BinomialVariantCaller(include_homozygous=True)
        score = caller.score_locus({"A": 0, "C": 20}, "A")
        assert score.loglik_hom[0] > score.loglik_het[0]

    def test_het_only_variant_underscores_homozygous_sites(self) -> None:
        """The het-only caller must score a hom-alt site strictly lower."""
        with_hom = llr_for(20, 20, include_homozygous=True)
        het_only = llr_for(20, 20, include_homozygous=False)
        assert with_hom > het_only


class TestDepthAwareness:
    def test_low_depth_evidence_is_weaker_than_high_depth(self) -> None:
        """Test 4: n=5,k=2 must carry far less evidence than n=50,k=20."""
        weak = llr_for(5, 2)
        strong = llr_for(50, 20)
        assert strong > weak
        assert strong > 3 * max(weak, 1.0)

    def test_evidence_grows_with_depth_at_fixed_vaf(self) -> None:
        """At a fixed VAF of 0.5 the LLR must increase with depth."""
        values = [llr_for(n, n // 2) for n in (10, 20, 40, 80)]
        assert all(b > a for a, b in zip(values, values[1:]))


class TestMonotonicity:
    def test_llr_non_decreasing_in_k_at_fixed_n(self) -> None:
        """Test 6: at fixed n, more alternate support cannot reduce evidence."""
        values = [llr_for(30, k) for k in range(0, 16)]
        assert all(b >= a - 1e-9 for a, b in zip(values, values[1:]))

    def test_small_k_at_high_depth_is_not_a_variant(self) -> None:
        """A couple of stray reads out of many must not look like a SNP."""
        assert llr_for(200, 3) < llr_for(200, 100)
        assert llr_for(200, 3) < 0.0


class TestNumericalStability:
    def test_large_depth_stays_finite(self) -> None:
        """Test 7: extreme depths must not overflow or produce NaN."""
        for n in (1_000, 100_000, 1_000_000):
            for k in (0, n // 2, n):
                assert np.isfinite(llr_for(n, k))

    def test_zero_depth_is_neutral(self) -> None:
        """A locus with no reads carries no evidence either way."""
        assert llr_for(0, 0) == 0.0

    def test_vectorised_matches_scalar(self) -> None:
        """Batched scoring must agree with per-locus scoring."""
        caller = BinomialVariantCaller()
        counts = np.array([[10, 10, 0, 0], [20, 0, 0, 0], [0, 0, 5, 15]], dtype=float)
        reference = np.array([1, 1, 4])  # A, A, T
        batched = caller.score_counts(counts, reference).llr
        scalar = [
            float(caller.score_locus({"A": 10, "C": 10}, "A").llr[0]),
            float(caller.score_locus({"A": 20}, "A").llr[0]),
            float(caller.score_locus({"G": 5, "T": 15}, "T").llr[0]),
        ]
        assert np.allclose(batched, scalar)


class TestCandidateAllele:
    def test_reference_base_is_never_the_candidate_alt(self) -> None:
        """The most-supported base being the reference must give k from the runner-up."""
        caller = BinomialVariantCaller()
        score = caller.score_counts(np.array([[30.0, 4.0, 1.0, 0.0]]), np.array([1]))
        assert int(score.k[0]) == 4          # C, not the A majority
        assert int(score.alt_index[0]) == 1

    def test_epsilon_from_quality_is_used(self) -> None:
        """Quality-derived epsilon must change the score versus the fixed rate."""
        caller = BinomialVariantCaller(quality_derived_epsilon=True)
        counts = np.array([[15.0, 5.0, 0.0, 0.0]])
        high = caller.score_counts(counts, np.array([1]), quality_sum=np.array([20 * 40.0]))
        low = caller.score_counts(counts, np.array([1]), quality_sum=np.array([20 * 13.0]))
        assert high.epsilon[0] < low.epsilon[0]
        assert high.llr[0] > low.llr[0]


class TestMultipleTesting:
    def test_tail_pvalue_is_a_probability(self) -> None:
        p = binomial_tail_pvalue(np.array([0.0, 5.0, 20.0]), np.array([20.0, 20.0, 20.0]),
                                 np.array([0.01, 0.01, 0.01]))
        assert np.all((p >= 0) & (p <= 1))
        assert p[0] > p[1] > p[2]

    def test_benjamini_hochberg_rejects_the_obvious(self) -> None:
        pvalues = np.array([1e-12, 1e-10, 0.4, 0.9])
        rejected = benjamini_hochberg(pvalues, alpha=0.05)
        assert rejected[0] and rejected[1]
        assert not rejected[2] and not rejected[3]


def test_invalid_configuration_rejected() -> None:
    with pytest.raises(ValueError):
        BinomialVariantCaller(error_rate=0.0)
    with pytest.raises(ValueError):
        BinomialVariantCaller(error_floor=1.0)
