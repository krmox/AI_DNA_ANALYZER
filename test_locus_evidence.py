"""Unit tests for the 44-channel locus evidence and the PB-residual model.

The features are hand-built from synthetic read matrices, so a failure
localises to the feature definition rather than to extraction. Every test
corresponds to a property the experiment's interpretation depends on: if the
ref/alt split, the count semantics, or the model's zero-initialization were
wrong, a comparison against the Poisson-binomial baseline would be measuring
something other than what it claims.
"""

from __future__ import annotations

import numpy as np
import pytest
import torch

from locus_evidence import (
    FEATURE_DIM,
    FEATURE_GROUPS,
    FEATURE_NAMES,
    build_locus_features,
    resolve_group_indices,
)
from model_pb_residual import PoissonBinomialResidualMamba
from quality_error_model import candidate_alt
from read_level_pileup import (
    FLAG_COUNTED,
    FLAG_GAP,
    FLAG_NEAR_END,
    FLAG_VALID,
    READ_DIM,
    reconstruct_counts,
)
from residual_prior import prior_logits_from_llr, snp_decision_score
from test_quality_error import BASE_CODE, REFERENCE_INDEX, make_locus


def feature(name: str, features: np.ndarray) -> float:
    """One named feature from a single-locus feature matrix."""
    return float(features[0, FEATURE_NAMES.index(name)])


def locus(observations, reference: str = "A", gaps: int = 0, strands=None,
          near_end=None, mapq=None) -> np.ndarray:
    """Build a one-locus ``[1, 48, READ_DIM]`` read matrix with rich attributes.

    Args:
        observations: ``(base, phred)`` pairs passing the base-quality filter.
        reference: Reference base.
        gaps: Deletion reads to append.
        strands: Optional per-observation strand (0 forward, 1 reverse).
        near_end: Optional per-observation near-end flag.
        mapq: Optional per-observation mapping quality.

    Returns:
        ``(reads, counts)`` ready for :func:`build_locus_features`.
    """
    matrix = np.zeros((1, 48, READ_DIM), dtype=np.uint8)
    matrix[0, :, 0] = 255
    for row, (base, quality) in enumerate(observations):
        flags = FLAG_VALID | FLAG_COUNTED
        if near_end is not None and near_end[row]:
            flags |= FLAG_NEAR_END
        matrix[0, row] = [
            BASE_CODE[base], quality,
            60 if mapq is None else mapq[row],
            0 if strands is None else strands[row],
            128, 255, flags, 0,
        ]
    for offset in range(gaps):
        matrix[0, len(observations) + offset] = [4, 0, 60, 0, 0, 0,
                                                 FLAG_VALID | FLAG_GAP, 0]
    counts = reconstruct_counts(matrix, np.array([REFERENCE_INDEX[reference]]))
    return matrix, counts


class TestLayout:
    def test_groups_partition_the_features(self) -> None:
        """Every feature belongs to exactly one semantic group."""
        covered = [name for group in FEATURE_GROUPS.values() for name in group]
        assert sorted(covered) == sorted(FEATURE_NAMES)
        assert len(covered) == FEATURE_DIM

    def test_group_indices_resolve(self) -> None:
        indices = resolve_group_indices(("alt_quality",))
        assert [FEATURE_NAMES[i] for i in indices] == list(FEATURE_GROUPS["alt_quality"])


class TestCountSemantics:
    def test_k_and_n_recoverable(self) -> None:
        """The failure that motivated this module: integer n and k must survive."""
        reads, counts = locus([("A", 35)] * 8 + [("C", 35)] * 2)
        features = build_locus_features(reads, counts)
        assert np.expm1(feature("log_n", features) * 4.0) == pytest.approx(10.0, abs=1e-4)
        assert np.expm1(feature("log_k", features) * 4.0) == pytest.approx(2.0, abs=1e-4)
        assert feature("k_small", features) == pytest.approx(0.2)
        assert feature("vaf", features) == pytest.approx(0.2)

    def test_gap_reads_excluded_from_n(self) -> None:
        """Deletion reads must not enter the substitution denominator."""
        reads, counts = locus([("A", 35)] * 6 + [("C", 35)] * 2, gaps=4)
        features = build_locus_features(reads, counts)
        assert np.expm1(feature("log_n", features) * 4.0) == pytest.approx(8.0, abs=1e-4)
        assert feature("gap_fraction", features) > 0.0

    def test_other_base_count(self) -> None:
        """A third competing base is counted separately from ref and alt."""
        reads, counts = locus([("A", 35)] * 6 + [("C", 35)] * 3 + [("G", 35)] * 1)
        features = build_locus_features(reads, counts)
        assert np.expm1(feature("log_other_count", features) * 4.0) == \
            pytest.approx(1.0, abs=1e-4)

    def test_zero_depth_is_all_zero_evidence(self) -> None:
        """An uncovered locus reads as 'no evidence', not as reference support."""
        reads, counts = locus([])
        features = build_locus_features(reads, counts)
        assert feature("log_n", features) == 0.0
        assert feature("vaf", features) == 0.0
        assert feature("alt_has_reads", features) == 0.0
        assert feature("ref_is_A", features) == 1.0


class TestRefAltQualitySplit:
    def test_alt_quality_isolated_from_reference_quality(self) -> None:
        """The decisive property: ALT quality must not be diluted by REF reads.

        This is precisely what the 14-channel `mean_base_quality` could not do.
        """
        reads, counts = locus([("A", 20)] * 8 + [("C", 40)] * 2)
        features = build_locus_features(reads, counts)
        assert feature("alt_mean_phred", features) == pytest.approx(40 / 40.0)
        assert feature("ref_mean_phred", features) == pytest.approx(20 / 40.0)
        # The pooled mean the old representation supplied would be 24.
        pooled = (8 * 20 + 2 * 40) / 10
        assert feature("alt_mean_phred", features) * 40 != pytest.approx(pooled)

    def test_alt_quality_changes_while_pooled_mean_is_held_constant(self) -> None:
        """Two loci with identical (n, k) and identical pooled mean quality.

        The old representation cannot tell these apart at all. The new one must.
        """
        reads_a, counts_a = locus([("A", 20)] * 8 + [("C", 40)] * 2)
        reads_b, counts_b = locus([("A", 25)] * 8 + [("C", 20)] * 2)
        a = build_locus_features(reads_a, counts_a)
        b = build_locus_features(reads_b, counts_b)
        pooled_a = (8 * 20 + 2 * 40) / 10
        pooled_b = (8 * 25 + 2 * 20) / 10
        assert pooled_a == pooled_b == 24.0
        assert feature("alt_mean_phred", a) != pytest.approx(feature("alt_mean_phred", b))
        assert feature("alt_minus_ref_phred", a) > feature("alt_minus_ref_phred", b)

    def test_min_and_max_capture_distribution_not_just_mean(self) -> None:
        """Same mean ALT quality, different spread, must give different features."""
        tight, tight_counts = locus([("A", 35)] * 8 + [("C", 30), ("C", 30)])
        wide, wide_counts = locus([("A", 35)] * 8 + [("C", 10), ("C", 50)])
        t = build_locus_features(tight, tight_counts)
        w = build_locus_features(wide, wide_counts)
        assert feature("alt_mean_phred", t) == pytest.approx(feature("alt_mean_phred", w))
        assert feature("alt_min_phred", w) < feature("alt_min_phred", t)
        assert feature("alt_max_phred", w) > feature("alt_max_phred", t)

    def test_fraction_q30(self) -> None:
        reads, counts = locus([("A", 35)] * 4 + [("C", 35), ("C", 35), ("C", 10), ("C", 10)])
        features = build_locus_features(reads, counts)
        assert feature("alt_fraction_q30", features) == pytest.approx(0.5)
        assert feature("ref_fraction_q30", features) == pytest.approx(1.0)


class TestStrandAndPosition:
    def test_strand_bias_detected_on_alt_reads_only(self) -> None:
        """One-sided ALT reads against balanced REF reads: the artefact signature.

        Poisson-binomial cannot represent this at all -- it sees only qualities.
        """
        strands = [0, 1, 0, 1, 0, 1, 0, 1] + [0, 0]
        reads, counts = locus([("A", 35)] * 8 + [("C", 35)] * 2, strands=strands)
        features = build_locus_features(reads, counts)
        assert feature("alt_strand_bias", features) == pytest.approx(1.0)
        assert feature("ref_strand_bias", features) == pytest.approx(0.0)

    def test_near_end_clustering_of_alt_reads(self) -> None:
        near_end = [False] * 8 + [True, True]
        reads, counts = locus([("A", 35)] * 8 + [("C", 35)] * 2, near_end=near_end)
        features = build_locus_features(reads, counts)
        assert feature("alt_near_end_fraction", features) == pytest.approx(1.0)
        assert feature("ref_near_end_fraction", features) == pytest.approx(0.0)
        assert feature("alt_minus_ref_near_end", features) == pytest.approx(1.0)

    def test_mapping_quality_split(self) -> None:
        mapq = [60] * 8 + [10, 10]
        reads, counts = locus([("A", 35)] * 8 + [("C", 35)] * 2, mapq=mapq)
        features = build_locus_features(reads, counts)
        assert feature("alt_mean_mapq", features) < feature("ref_mean_mapq", features)
        assert feature("alt_minus_ref_mapq", features) < 0


class TestConsistencyWithTheEstablishedCaller:
    def test_alt_candidate_matches_binomial_caller(self) -> None:
        """Feature construction must use the baseline's ALT, unchanged."""
        rng = np.random.default_rng(11)
        reads = np.zeros((300, 48, READ_DIM), dtype=np.uint8)
        reads[:, :, 0] = 255
        for index in range(300):
            for row in range(int(rng.integers(0, 25))):
                reads[index, row] = [int(rng.integers(0, 4)), int(rng.integers(5, 45)),
                                     60, int(rng.integers(0, 2)), 128, 255,
                                     FLAG_VALID | FLAG_COUNTED, 0]
        reference = rng.integers(0, 5, size=300)
        counts = reconstruct_counts(reads, reference)
        features = build_locus_features(reads, counts)
        _, k, n = candidate_alt(counts[:, 0:4], reference)
        assert np.allclose(np.expm1(features[:, FEATURE_NAMES.index("log_k")] * 4.0),
                           k, atol=1e-3)
        assert np.allclose(np.expm1(features[:, FEATURE_NAMES.index("log_n")] * 4.0),
                           n, atol=1e-3)

    def test_features_finite_and_deterministic(self) -> None:
        rng = np.random.default_rng(12)
        reads = np.zeros((200, 48, READ_DIM), dtype=np.uint8)
        reads[:, :, 0] = 255
        for index in range(200):
            for row in range(int(rng.integers(0, 30))):
                reads[index, row] = [int(rng.integers(0, 6)), int(rng.integers(0, 61)),
                                     int(rng.integers(0, 61)), int(rng.integers(0, 2)),
                                     128, 255, FLAG_VALID | FLAG_COUNTED, 0]
        counts = reconstruct_counts(reads, rng.integers(0, 5, size=200))
        first = build_locus_features(reads, counts)
        second = build_locus_features(reads, counts)
        assert np.all(np.isfinite(first))
        assert np.array_equal(first, second)

    def test_chunking_is_exact(self) -> None:
        """Every feature is a function of one locus, so blocking cannot change it.

        The cache builder processes the training region in blocks to bound
        memory; this is what makes that an implementation detail rather than a
        numerical choice.
        """
        rng = np.random.default_rng(13)
        reads = np.zeros((250, 48, READ_DIM), dtype=np.uint8)
        reads[:, :, 0] = 255
        for index in range(250):
            for row in range(int(rng.integers(0, 20))):
                reads[index, row] = [int(rng.integers(0, 4)), int(rng.integers(5, 45)),
                                     60, int(rng.integers(0, 2)), 128, 255,
                                     FLAG_VALID | FLAG_COUNTED, 0]
        counts = reconstruct_counts(reads, rng.integers(0, 5, size=250))
        whole = build_locus_features(reads, counts)
        blocked = np.concatenate([build_locus_features(reads[begin:begin + 37],
                                                       counts[begin:begin + 37])
                                  for begin in range(0, 250, 37)])
        assert np.array_equal(whole, blocked)

    def test_no_feature_encodes_the_label(self) -> None:
        """Static guard: the feature module must not reach for truth data."""
        from pathlib import Path

        source = Path(__file__).with_name("locus_evidence.py").read_text()
        code = "".join(part for index, part in enumerate(source.split('"""'))
                       if index % 2 == 0)
        body = "\n".join(line for line in code.splitlines()
                         if not line.strip().startswith("#"))
        for forbidden in ("vcf", "VCF", "label", "truth", "LABEL_SNP"):
            assert forbidden not in body


class TestModelInitialization:
    def test_residual_is_exactly_zero_at_init(self) -> None:
        """The model must start life exactly equal to the PB caller."""
        torch.manual_seed(0)
        model = PoissonBinomialResidualMamba().eval()
        features = torch.randn(2, 64, FEATURE_DIM)
        prior = torch.tensor(prior_logits_from_llr(
            np.random.default_rng(0).normal(size=128) * 10), dtype=torch.float
        ).reshape(2, 64, 4)
        logits, residual, _ = model(features, prior, return_parts=True)
        assert torch.all(residual == 0)
        assert torch.equal(logits, prior)
        assert torch.allclose(snp_decision_score(logits), snp_decision_score(prior))

    def test_gated_model_also_starts_at_the_prior(self) -> None:
        """Enabling the gate must not disturb the initialization guarantee."""
        torch.manual_seed(0)
        model = PoissonBinomialResidualMamba(use_gate=True).eval()
        features = torch.randn(2, 64, FEATURE_DIM)
        prior = torch.zeros(2, 64, 4)
        prior[..., 1] = torch.randn(2, 64) * 10
        logits, residual, gate = model(features, prior, return_parts=True)
        assert torch.all(residual == 0)
        assert torch.equal(logits, prior)
        assert torch.allclose(gate, torch.full_like(gate, 0.5))

    def test_residual_head_receives_gradient(self) -> None:
        """Zeroed output must not mean a dead layer."""
        torch.manual_seed(0)
        model = PoissonBinomialResidualMamba()
        features = torch.randn(2, 64, FEATURE_DIM)
        prior = torch.zeros(2, 64, 4)
        logits = model(features, prior)
        logits.sum().backward()
        gradient = model.backbone.classifier.weight.grad
        assert gradient is not None and torch.any(gradient != 0)

    def test_capacity_matches_the_previous_residual_model(self) -> None:
        """Only the input projection widens, so a gain is not extra capacity."""
        from model_residual_binomial import ResidualBinomialMamba

        torch.manual_seed(0)
        new = sum(p.numel() for p in PoissonBinomialResidualMamba().parameters())
        old = sum(p.numel() for p in ResidualBinomialMamba().parameters())
        # 14 -> 44 input channels at d_model=128, plus the two LayerNorm params.
        assert new - old == (44 - 14) * 128 + 2 * (44 - 14)

    def test_forward_shape(self) -> None:
        model = PoissonBinomialResidualMamba().eval()
        logits = model(torch.randn(3, 64, FEATURE_DIM), torch.zeros(3, 64, 4))
        assert logits.shape == (3, 64, 4)
