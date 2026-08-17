"""Unit tests for the per-locus quality-derived error model.

No BAM, no VCF, no cache: every locus here is a hand-built read matrix, so a
failure localises to the statistics and not to extraction. The properties
tested are the ones the experiment's interpretation depends on -- if the Phred
conversion, the reference-read selection, the ALT-candidate rule or the
Poisson-binomial reduction were wrong, the comparison against Binomial v1 would
be measuring something other than what it claims.
"""

from __future__ import annotations

import numpy as np
import pytest

from binomial_baseline import BinomialVariantCaller
from quality_error_model import (
    EPSILON_CEILING,
    EPSILON_FLOOR,
    ESTIMATORS,
    FIXED_EPSILON,
    binomial_llr,
    candidate_alt,
    estimate_epsilon,
    extract_quality_evidence,
    phred_to_error,
    poisson_binomial_llr,
)
from read_level_pileup import (
    FLAG_COUNTED,
    FLAG_GAP,
    FLAG_VALID,
    READ_DIM,
    reconstruct_counts,
)

BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
REFERENCE_INDEX = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 0}


def make_locus(observations: list[tuple[str, int]], gaps: int = 0,
               uncounted: list[tuple[str, int]] | None = None,
               max_reads: int = 48) -> np.ndarray:
    """Build one ``[1, max_reads, READ_DIM]`` read matrix.

    Args:
        observations: ``(base, phred)`` pairs that pass the base-quality filter.
        gaps: Number of deletion/ref-skip reads to add. These are valid reads
            but carry no base, so the caller must exclude them.
        uncounted: ``(base, phred)`` pairs present in the BAM but failing the
            base-quality filter, i.e. invisible to the binomial denominator.
        max_reads: Rows per locus.

    Returns:
        The uint8 read matrix.
    """
    matrix = np.zeros((1, max_reads, READ_DIM), dtype=np.uint8)
    matrix[0, :, 0] = 255  # pad code
    row = 0
    for base, quality in observations:
        matrix[0, row] = [BASE_CODE[base], quality, 60, 0, 128, 255,
                          FLAG_VALID | FLAG_COUNTED, 0]
        row += 1
    for base, quality in (uncounted or []):
        matrix[0, row] = [BASE_CODE[base], quality, 60, 0, 128, 255, FLAG_VALID, 0]
        row += 1
    for _ in range(gaps):
        matrix[0, row] = [4, 0, 60, 0, 0, 0, FLAG_VALID | FLAG_GAP, 0]
        row += 1
    return matrix


def evidence_for(observations, reference: str = "A", **kwargs):
    """Quality evidence for a single hand-built locus."""
    reads = make_locus(observations, **kwargs)
    return extract_quality_evidence(reads, np.array([REFERENCE_INDEX[reference]]))


class TestPhredConversion:
    @pytest.mark.parametrize("quality,expected", [(10, 0.1), (20, 0.01), (30, 0.001),
                                                  (0, 1.0), (40, 1e-4)])
    def test_known_values(self, quality: int, expected: float) -> None:
        """Q10 -> 0.1, Q20 -> 0.01, Q30 -> 0.001, exactly."""
        assert phred_to_error(quality) == pytest.approx(expected, rel=1e-12)

    def test_vectorized_matches_scalar(self) -> None:
        """The array path must agree elementwise with the scalar path."""
        qualities = np.arange(0, 61)
        vector = phred_to_error(qualities)
        scalar = np.array([phred_to_error(int(q)) for q in qualities])
        assert np.array_equal(vector, scalar)


class TestEvidenceExtraction:
    def test_gap_reads_excluded_from_denominator(self) -> None:
        """Deletion reads are valid reads but must not enter ``n``."""
        evidence = evidence_for([("A", 30), ("C", 30)], gaps=5)
        assert int(evidence.n_counted[0]) == 2

    def test_base_quality_filtered_reads_excluded(self) -> None:
        """Reads without FLAG_COUNTED are excluded, matching pileup_counts."""
        evidence = evidence_for([("A", 30)], uncounted=[("C", 2)])
        assert int(evidence.n_counted[0]) == 1

    def test_reference_reads_identified(self) -> None:
        """Only reads carrying the reference base count as reference support."""
        evidence = evidence_for([("A", 30), ("A", 30), ("C", 30)], reference="A")
        assert int(evidence.n_reference[0]) == 2

    def test_reference_n_has_no_supporting_reads(self) -> None:
        """A locus whose reference is N has no reference-supporting reads."""
        evidence = evidence_for([("A", 30), ("C", 30)], reference="N")
        assert int(evidence.n_reference[0]) == 0

    def test_matches_reconstructed_counts(self) -> None:
        """``n`` must equal the A/C/G/T sum the aggregate representation gives."""
        reads = make_locus([("A", 30)] * 7 + [("C", 25)] * 3, gaps=2,
                           uncounted=[("G", 2)])
        counts = reconstruct_counts(reads, np.array([REFERENCE_INDEX["A"]]))
        evidence = extract_quality_evidence(reads, np.array([REFERENCE_INDEX["A"]]))
        assert int(evidence.n_counted[0]) == int(counts[0, 0:4].sum()) == 10


class TestEstimators:
    def test_high_quality_reference_gives_low_epsilon(self) -> None:
        """Q40 reference reads must yield a smaller epsilon than Q20 ones."""
        high = estimate_epsilon(evidence_for([("A", 40)] * 10), "ref_mean")[0][0]
        low = estimate_epsilon(evidence_for([("A", 20)] * 10), "ref_mean")[0][0]
        assert high < low

    def test_low_quality_reference_gives_high_epsilon(self) -> None:
        """Q10 reference reads must yield epsilon near 0.1."""
        epsilon = estimate_epsilon(evidence_for([("A", 10)] * 10), "ref_mean")[0][0]
        assert epsilon == pytest.approx(0.1, rel=1e-9)

    def test_mean_is_arithmetic_not_geometric(self) -> None:
        """Estimator A averages probabilities, not Phred scores.

        This is the concrete difference from the existing ``binomial_v2``,
        which computes ``10 ** (-mean(Q) / 10)``. Half Q40 and half Q10 gives
        an arithmetic mean of ~0.05 but a geometric mean of ~0.0032.
        """
        evidence = evidence_for([("A", 40)] * 5 + [("A", 10)] * 5)
        epsilon = estimate_epsilon(evidence, "mean")[0][0]
        assert epsilon == pytest.approx((5 * 1e-4 + 5 * 0.1) / 10, rel=1e-9)
        assert epsilon != pytest.approx(10 ** (-25 / 10), rel=1e-3)

    def test_median_resists_one_terrible_read(self) -> None:
        """One Q2 read moves the mean far more than the median."""
        evidence = evidence_for([("A", 35)] * 9 + [("A", 2)])
        mean = estimate_epsilon(evidence, "mean")[0][0]
        median = estimate_epsilon(evidence, "median")[0][0]
        assert median < mean / 10

    def test_reference_only_ignores_alt_read_quality(self) -> None:
        """Estimator C must not change when the ALT reads' quality changes.

        This is the anti-circularity property: the error rate the likelihood
        assumes cannot be driven by the very reads whose status is in question.
        """
        good_alt = evidence_for([("A", 35)] * 8 + [("C", 40)] * 2)
        bad_alt = evidence_for([("A", 35)] * 8 + [("C", 5)] * 2)
        assert estimate_epsilon(good_alt, "ref_mean")[0][0] == \
            estimate_epsilon(bad_alt, "ref_mean")[0][0]

    def test_trimmed_drops_the_extremes(self) -> None:
        """Estimator D must discard the worst read that the mean retains."""
        evidence = evidence_for([("A", 35)] * 9 + [("A", 2)])
        mean = estimate_epsilon(evidence, "ref_mean")[0][0]
        trimmed = estimate_epsilon(evidence, "ref_trimmed")[0][0]
        assert trimmed < mean

    def test_trimmed_falls_back_below_min_reads(self) -> None:
        """With fewer than 5 reads, D equals C rather than trimming to nothing."""
        evidence = evidence_for([("A", 35), ("A", 20), ("A", 10)])
        assert estimate_epsilon(evidence, "ref_trimmed")[0][0] == \
            estimate_epsilon(evidence, "ref_mean")[0][0]

    def test_fallback_when_no_reference_reads(self) -> None:
        """A locus with only ALT reads falls back to the all-read mean."""
        evidence = evidence_for([("C", 20)] * 4, reference="A")
        epsilon, diagnostics = estimate_epsilon(evidence, "ref_mean")
        assert epsilon[0] == pytest.approx(0.01, rel=1e-9)
        assert diagnostics["fallback_all_read"] == 1

    def test_fallback_when_no_reads_at_all(self) -> None:
        """A zero-depth locus falls back to the fixed constant."""
        evidence = evidence_for([])
        epsilon, diagnostics = estimate_epsilon(evidence, "ref_mean")
        assert epsilon[0] == FIXED_EPSILON
        assert diagnostics["fallback_fixed"] == 1

    def test_clamped_only_at_the_declared_bounds(self) -> None:
        """Q60 clamps to the floor; Q1 clamps to the ceiling; Q20 does neither."""
        assert estimate_epsilon(evidence_for([("A", 60)] * 10), "ref_mean")[0][0] == \
            EPSILON_FLOOR
        assert estimate_epsilon(evidence_for([("A", 1)] * 10), "ref_mean")[0][0] == \
            EPSILON_CEILING
        assert estimate_epsilon(evidence_for([("A", 20)] * 10), "ref_mean")[0][0] == \
            pytest.approx(0.01, rel=1e-9)

    @pytest.mark.parametrize("method", ESTIMATORS)
    def test_epsilon_finite_and_in_range(self, method: str) -> None:
        """Every estimator stays finite and inside the declared range."""
        rng = np.random.default_rng(0)
        reads = np.zeros((200, 48, READ_DIM), dtype=np.uint8)
        reads[:, :, 0] = 255
        for locus in range(200):
            depth = int(rng.integers(0, 30))
            for row in range(depth):
                reads[locus, row] = [int(rng.integers(0, 4)), int(rng.integers(0, 61)),
                                     60, 0, 128, 255, FLAG_VALID | FLAG_COUNTED, 0]
        reference = rng.integers(0, 5, size=200)
        evidence = extract_quality_evidence(reads, reference)
        epsilon, _ = estimate_epsilon(evidence, method)
        assert np.all(np.isfinite(epsilon))
        assert np.all(epsilon >= EPSILON_FLOOR) and np.all(epsilon <= EPSILON_CEILING)

    @pytest.mark.parametrize("method", ESTIMATORS)
    def test_deterministic(self, method: str) -> None:
        """Identical input must give bit-identical output."""
        evidence = evidence_for([("A", 33)] * 6 + [("C", 21)] * 2)
        first = estimate_epsilon(evidence, method)[0]
        second = estimate_epsilon(evidence, method)[0]
        assert np.array_equal(first, second)

    def test_unknown_estimator_rejected(self) -> None:
        with pytest.raises(ValueError):
            estimate_epsilon(evidence_for([("A", 30)]), "not_an_estimator")


class TestCandidateAlt:
    def test_alt_selection_unchanged_from_baseline(self) -> None:
        """``candidate_alt`` must reproduce the established caller exactly."""
        rng = np.random.default_rng(1)
        counts = rng.integers(0, 40, size=(500, 4)).astype(np.float64)
        reference = rng.integers(0, 5, size=500)
        alt, k, n = candidate_alt(counts, reference)
        baseline = BinomialVariantCaller().score_counts(counts, reference)
        assert np.array_equal(alt, baseline.alt_index)
        assert np.array_equal(k, baseline.k)
        assert np.array_equal(n, baseline.n)

    def test_alt_never_equals_reference(self) -> None:
        """With a real reference base, the candidate ALT is never the reference."""
        rng = np.random.default_rng(2)
        counts = rng.integers(0, 40, size=(500, 4)).astype(np.float64)
        reference = rng.integers(1, 5, size=500)   # A..T only
        alt, _, _ = candidate_alt(counts, reference)
        assert np.all(alt != reference - 1)


class TestBinomialLikelihood:
    def test_matches_baseline_at_fixed_epsilon(self) -> None:
        """At epsilon = 0.01 the rebuilt likelihood must equal Binomial v1."""
        rng = np.random.default_rng(3)
        counts = rng.integers(0, 40, size=(2000, 4)).astype(np.float64)
        reference = rng.integers(0, 5, size=2000)
        baseline = BinomialVariantCaller(error_rate=FIXED_EPSILON).score_counts(
            counts, reference)
        rebuilt = binomial_llr(baseline.k, baseline.n,
                               np.full(2000, FIXED_EPSILON))["llr"]
        assert np.allclose(rebuilt, baseline.llr, atol=1e-12)

    def test_llr_monotonic_in_k(self) -> None:
        """At fixed n and epsilon, more ALT reads must never lower the LLR."""
        n = np.full(21, 20.0)
        k = np.arange(21, dtype=np.float64)
        llr = binomial_llr(k, n, np.full(21, 0.01))["llr"]
        assert np.all(np.diff(llr) > 0)

    def test_llr_monotonic_in_epsilon(self) -> None:
        """A higher assumed error rate must weaken the evidence at fixed k/n."""
        epsilon = np.array([0.001, 0.005, 0.01, 0.05, 0.1])
        llr = binomial_llr(np.full(5, 2.0), np.full(5, 10.0), epsilon)["llr"]
        assert np.all(np.diff(llr) < 0)

    def test_zero_depth_is_neutral(self) -> None:
        """A locus with no reads must score exactly 0, not NaN."""
        result = binomial_llr(np.zeros(1), np.zeros(1), np.array([0.01]))
        assert result["llr"][0] == 0.0

    def test_large_n_stays_finite(self) -> None:
        """Numerical stability at depths far beyond this dataset's maximum."""
        n = np.array([1e4, 1e5, 1e6])
        result = binomial_llr(n * 0.5, n, np.full(3, 0.01))
        assert np.all(np.isfinite(result["llr"]))

    def test_batch_matches_scalar(self) -> None:
        """Vectorized evaluation must equal per-locus evaluation."""
        rng = np.random.default_rng(4)
        n = rng.integers(1, 45, size=300).astype(np.float64)
        k = np.floor(n * rng.random(300))
        epsilon = 10 ** (-rng.uniform(10, 45, size=300) / 10)
        batch = binomial_llr(k, n, epsilon)["llr"]
        scalar = np.array([
            binomial_llr(k[i: i + 1], n[i: i + 1], epsilon[i: i + 1])["llr"][0]
            for i in range(300)])
        assert np.allclose(batch, scalar, atol=1e-12)


class TestPoissonBinomial:
    def test_reduces_to_binomial_when_qualities_equal(self) -> None:
        """With identical per-read quality the two formulations must coincide.

        This is the justification for treating Poisson-binomial as a strictly
        weaker assumption rather than a different model.
        """
        for depth, alt in ((10, 0), (10, 2), (10, 5), (10, 10), (25, 3), (7, 7)):
            evidence = evidence_for([("C", 20)] * alt + [("A", 20)] * (depth - alt))
            poisson = poisson_binomial_llr(evidence, np.array([alt]))["llr"][0]
            binomial = binomial_llr(np.array([float(alt)]), np.array([float(depth)]),
                                    np.array([0.01]))["llr"][0]
            assert poisson == pytest.approx(binomial, abs=1e-9), f"{depth=} {alt=}"

    def test_pmf_normalises(self) -> None:
        """The exact log-pmf must sum to 1 over all possible counts."""
        from quality_error_model import _log_poisson_binomial_pmf

        rng = np.random.default_rng(5)
        q = rng.uniform(0.001, 0.9, size=(20, 48))
        active = np.arange(48)[None, :] < rng.integers(1, 49, size=(20, 1))
        log_pmf = _log_poisson_binomial_pmf(q, active)
        assert np.allclose(np.exp(log_pmf).sum(axis=1), 1.0, atol=1e-10)

    def test_mixed_quality_differs_from_plug_in(self) -> None:
        """With heterogeneous qualities the two formulations must disagree.

        If they agreed, the Poisson-binomial arm could not be informative and
        there would be nothing to test.
        """
        evidence = evidence_for([("C", 35)] * 3 + [("A", 8)] * 7)
        poisson = poisson_binomial_llr(evidence, np.array([3]))["llr"][0]
        epsilon = estimate_epsilon(evidence, "mean")[0]
        plug_in = binomial_llr(np.array([3.0]), np.array([10.0]), epsilon)["llr"][0]
        assert abs(poisson - plug_in) > 1e-3

    def test_zero_depth_is_neutral(self) -> None:
        """No reads means no evidence, scored as exactly 0."""
        assert poisson_binomial_llr(evidence_for([]), np.array([0]))["llr"][0] == 0.0

    def test_finite_over_random_loci(self) -> None:
        """No NaN or infinity anywhere over a randomised sweep."""
        rng = np.random.default_rng(6)
        reads = np.zeros((500, 48, READ_DIM), dtype=np.uint8)
        reads[:, :, 0] = 255
        alt_counts = np.zeros(500, dtype=np.int64)
        for locus in range(500):
            depth = int(rng.integers(0, 44))
            for row in range(depth):
                reads[locus, row] = [int(rng.integers(0, 4)), int(rng.integers(0, 61)),
                                     60, 0, 128, 255, FLAG_VALID | FLAG_COUNTED, 0]
            alt_counts[locus] = int(rng.integers(0, depth + 1))
        evidence = extract_quality_evidence(reads, rng.integers(0, 5, size=500))
        result = poisson_binomial_llr(evidence, alt_counts)
        assert np.all(np.isfinite(result["llr"]))

    def test_chunking_does_not_change_results(self) -> None:
        """Block size is a memory knob, never a numerical one."""
        rng = np.random.default_rng(7)
        reads = np.zeros((300, 48, READ_DIM), dtype=np.uint8)
        reads[:, :, 0] = 255
        for locus in range(300):
            for row in range(int(rng.integers(1, 30))):
                reads[locus, row] = [int(rng.integers(0, 4)), int(rng.integers(5, 45)),
                                     60, 0, 128, 255, FLAG_VALID | FLAG_COUNTED, 0]
        evidence = extract_quality_evidence(reads, rng.integers(1, 5, size=300))
        k = np.minimum(evidence.n_counted, 3)
        assert np.array_equal(poisson_binomial_llr(evidence, k, chunk=17)["llr"],
                              poisson_binomial_llr(evidence, k, chunk=4096)["llr"])


class TestPoissonBinomialDepthClipping:
    """Regression tests for the >=48x depth-clipping defect (devlog 12 §15, fixed devlog 13).

    The defect: ``evidence`` covers at most ``MAX_READS = 48`` reads while ``k``
    comes from the uncapped count matrix, and the two were reconciled by
    clipping ``k`` to the tensor width. Where that pushed ``k`` past a locus'
    number of *counted* slots, every hypothesis returned ``-inf``, the LLR came
    out ``nan``, and ``nan_to_num`` mapped the strongest possible variant
    evidence to exactly ``0.0``.

    Each test below states the property that must hold rather than a magic
    number, except where the number is the point (``!= 0``).
    """

    def test_legacy_path_reproduces_the_defect(self) -> None:
        """Guard on the A/B control itself: without the fix, a strong ALT locus scores 0.

        If this ever stops failing the way it used to, the benchmark's
        before/after comparison would be measuring nothing.
        """
        # 40 counted ALT reads + 8 gap reads fills all 48 slots; the uncapped
        # pileup at ~70x reports k = 54.
        evidence = evidence_for([("C", 35)] * 40, gaps=8, reference="A")
        legacy = poisson_binomial_llr(evidence, np.array([54]), legacy_truncation=True)
        assert legacy["llr"][0] == 0.0
        assert legacy["n_non_finite"] == 1

    def test_high_depth_strong_evidence_is_no_longer_zeroed(self) -> None:
        """The previously failing case: k above the retained read count."""
        evidence = evidence_for([("C", 35)] * 40, gaps=8, reference="A")
        result = poisson_binomial_llr(evidence, np.array([54]))
        assert np.isfinite(result["llr"][0])
        assert result["llr"][0] > 0, "40/40 ALT reads is evidence *for* a variant"
        assert result["n_out_of_support"] == 1
        assert result["n_non_finite"] == 0

    def test_high_depth_strong_evidence_clears_the_frozen_threshold(self) -> None:
        """The defect's consequence was a confident false negative; it must be gone."""
        from cascade import FROZEN_PB_THRESHOLD

        evidence = evidence_for([("C", 35)] * 40, gaps=8, reference="A")
        assert poisson_binomial_llr(evidence, np.array([54]))["llr"][0] > FROZEN_PB_THRESHOLD

    def test_alt_index_recounts_over_the_retained_reads(self) -> None:
        """With ``alt_index``, ``k`` is taken from the evidence, not from the caller.

        The retained sample here is 30 ALT + 10 REF, so a caller-supplied
        ``k = 54`` (the uncapped count at ~70x) must be replaced by 30, and the
        answer must equal the answer for a locus that simply had k = 30.
        """
        evidence = evidence_for([("C", 35)] * 30 + [("A", 35)] * 10, reference="A")
        with_alt = poisson_binomial_llr(evidence, np.array([54]), alt_index=np.array([1]))
        assert int(with_alt["k_effective"][0]) == 30
        direct = poisson_binomial_llr(evidence, np.array([30]))
        assert with_alt["llr"][0] == direct["llr"][0]

    def test_clipping_no_longer_overstates_evidence(self) -> None:
        """The other half of the defect: clipping asserted 'all retained reads are ALT'.

        With 30 of 48 retained reads carrying ALT, the pre-fix path scored the
        locus as though all 48 did, i.e. strictly more confident than the truth.
        """
        evidence = evidence_for([("C", 35)] * 30 + [("A", 35)] * 18, reference="A")
        legacy = poisson_binomial_llr(evidence, np.array([54]), legacy_truncation=True)["llr"][0]
        fixed = poisson_binomial_llr(evidence, np.array([54]), alt_index=np.array([1]))["llr"][0]
        assert legacy > fixed, "the pre-fix path was over-confident here, not zeroed"

    @pytest.mark.parametrize("depth,alt", [(10, 0), (10, 5), (30, 15), (30, 30), (48, 24),
                                           (48, 48), (1, 1), (2, 0)])
    def test_at_or_below_tensor_width_the_fix_is_a_no_op(self, depth: int, alt: int) -> None:
        """Every result at <=48x must be reproduced bit-for-bit.

        This is what licenses carrying the <=30x results of devlogs 3-12
        forward unchanged: below the tensor width, ``k <= n_counted`` always,
        so all three code paths coincide exactly.
        """
        evidence = evidence_for([("C", 32)] * alt + [("A", 32)] * (depth - alt), reference="A")
        legacy = poisson_binomial_llr(evidence, np.array([alt]), legacy_truncation=True)["llr"][0]
        fixed = poisson_binomial_llr(evidence, np.array([alt]))["llr"][0]
        exact = poisson_binomial_llr(evidence, np.array([alt]), alt_index=np.array([1]))["llr"][0]
        assert fixed == legacy
        assert exact == legacy

    def test_random_sweep_below_tensor_width_is_bit_identical(self) -> None:
        """The same no-op claim, over 500 randomised <=44x loci at mixed qualities."""
        rng = np.random.default_rng(1301)
        reads = np.zeros((500, 48, READ_DIM), dtype=np.uint8)
        reads[:, :, 0] = 255
        alt_counts = np.zeros(500, dtype=np.int64)
        for locus in range(500):
            depth = int(rng.integers(0, 44))
            for row in range(depth):
                reads[locus, row] = [int(rng.integers(0, 4)), int(rng.integers(0, 61)),
                                     60, 0, 128, 255, FLAG_VALID | FLAG_COUNTED, 0]
            alt_counts[locus] = int(rng.integers(0, depth + 1))
        evidence = extract_quality_evidence(reads, rng.integers(0, 5, size=500))
        assert np.array_equal(poisson_binomial_llr(evidence, alt_counts)["llr"],
                              poisson_binomial_llr(evidence, alt_counts,
                                                   legacy_truncation=True)["llr"])

    def test_k_never_leaves_the_pmf_support(self) -> None:
        """Numerical edge: whatever the caller passes, the scored count is admissible."""
        evidence = evidence_for([("C", 30)] * 5 + [("A", 30)] * 5, gaps=10, reference="A")
        for k in (-3, 0, 10, 11, 48, 500):
            result = poisson_binomial_llr(evidence, np.array([k]))
            assert 0 <= int(result["k_effective"][0]) <= int(evidence.n_counted[0])
            assert np.isfinite(result["llr"][0])

    def test_zero_counted_reads_stays_neutral_at_high_k(self) -> None:
        """A locus with only gap reads scores 0 -- the one place 0 is correct."""
        evidence = evidence_for([], gaps=48, reference="A")
        assert poisson_binomial_llr(evidence, np.array([54]))["llr"][0] == 0.0

    def test_llr_is_monotone_in_retained_alt_count(self) -> None:
        """Sanity on the repaired axis: more ALT support cannot mean less evidence."""
        scores = [poisson_binomial_llr(
            evidence_for([("C", 33)] * alt + [("A", 33)] * (48 - alt), reference="A"),
            np.array([alt]))["llr"][0] for alt in range(0, 49, 4)]
        assert all(later >= earlier for earlier, later in zip(scores, scores[1:]))

    def test_out_of_support_diagnostic_counts_affected_loci(self) -> None:
        """The diagnostic that tells a benchmark how much data the defect touched."""
        reads = np.concatenate([make_locus([("C", 35)] * 40, gaps=8),
                                make_locus([("C", 35)] * 10)])
        evidence = extract_quality_evidence(reads, np.array([1, 1]))
        result = poisson_binomial_llr(evidence, np.array([54, 10]))
        assert result["n_out_of_support"] == 1


class TestLeakage:
    def test_module_never_mentions_truth_sources(self) -> None:
        """Static guard: the estimator module must not touch labels or VCFs.

        Cheap, but it is the check that would actually fire if a future edit
        reached for the truth data to 'improve' the estimate.
        """
        from pathlib import Path

        source = Path(__file__).with_name("quality_error_model.py").read_text()
        # Segments alternate code / docstring on triple quotes; keep the code.
        code = "".join(segment for index, segment in enumerate(source.split('"""'))
                       if index % 2 == 0)
        body = "\n".join(line for line in code.splitlines()
                         if not line.strip().startswith("#"))
        for forbidden in ("vcf", "VCF", "label", "truth", "LABEL_SNP"):
            assert forbidden not in body, f"{forbidden!r} appears in executable code"
