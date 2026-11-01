"""Tests for the residual binomial variant caller.

Three properties have to hold before the experiment means anything:

1. the residual is *exactly* zero at initialization, so the untrained model is
   the binomial caller;
2. the prior is a faithful, monotone re-expression of the binomial LLR, so the
   binomial's threshold transfers;
3. the LLR this experiment computes is byte-identical to the one the existing
   ``binomial_baseline`` / ``evaluate_binomial_baseline`` pipeline computes,
   so the "baseline" being improved upon is the real published baseline and
   not a re-implementation that happens to be close.
"""

from __future__ import annotations

import numpy as np
import pytest
import torch

from binomial_baseline import BinomialVariantCaller
from config import LABEL_NORMAL, LABEL_SNP
from llr_features import raw_llr_from_counts
from model_residual_binomial import ResidualBinomialMamba
from residual_prior import (
    NUM_CLASSES,
    prior_logits_from_counts,
    prior_logits_from_llr,
    snp_decision_score,
)

SEQ_LEN = 8
FEATURE_DIM = 14


@pytest.fixture
def counts() -> np.ndarray:
    """A small synthetic count matrix spanning depth and VAF regimes."""
    rng = np.random.default_rng(0)
    n = 256
    matrix = np.zeros((n, 10), dtype=np.float64)
    depth = rng.integers(0, 40, n)
    for i in range(n):
        reference = rng.integers(1, 5)
        alt = rng.integers(1, 5)
        k = rng.integers(0, depth[i] + 1)
        matrix[i, reference - 1] += depth[i] - k
        matrix[i, alt - 1] += k
        matrix[i, 6] = depth[i]
        matrix[i, 7] = 35.0 * depth[i]
        matrix[i, 8] = 60.0 * depth[i]
        matrix[i, 9] = reference
    return matrix


@pytest.fixture
def model() -> ResidualBinomialMamba:
    torch.manual_seed(1234)
    return ResidualBinomialMamba().eval()


# --- 1. residual starts at exactly zero --------------------------------------


def test_residual_is_exactly_zero_at_initialization(model):
    """The whole experiment rests on this: no tolerance, exact equality."""
    features = torch.randn(4, SEQ_LEN, FEATURE_DIM)
    residual = model.residual(features)
    assert torch.equal(residual, torch.zeros_like(residual))


def test_residual_head_parameters_are_zero(model):
    """Both weight and bias, not just the weight."""
    assert torch.equal(model.backbone.classifier.weight,
                       torch.zeros_like(model.backbone.classifier.weight))
    assert torch.equal(model.backbone.classifier.bias,
                       torch.zeros_like(model.backbone.classifier.bias))


def test_residual_is_zero_for_extreme_and_degenerate_inputs(model):
    """Zero output must not depend on the input being well-scaled."""
    for scale in (0.0, 1.0, 1e3, 1e6):
        features = torch.full((2, SEQ_LEN, FEATURE_DIM), scale)
        residual = model.residual(features)
        assert torch.equal(residual, torch.zeros_like(residual))


def test_forward_returns_prior_untouched_at_initialization(model, counts):
    """final_logits == prior_logits, exactly, before any gradient step."""
    prior = torch.tensor(
        prior_logits_from_counts(counts[:4 * SEQ_LEN], BinomialVariantCaller()),
        dtype=torch.float32).reshape(4, SEQ_LEN, NUM_CLASSES)
    features = torch.randn(4, SEQ_LEN, FEATURE_DIM)
    final = model(pileup_features=features, prior_logits=prior)
    assert torch.equal(final, prior)


def test_decision_score_at_initialization_equals_llr(model, counts):
    """The model's SNP score is the binomial LLR itself, to float32 precision."""
    llr = raw_llr_from_counts(counts[:4 * SEQ_LEN], BinomialVariantCaller())
    prior = torch.tensor(prior_logits_from_llr(llr),
                         dtype=torch.float32).reshape(4, SEQ_LEN, NUM_CLASSES)
    features = torch.randn(4, SEQ_LEN, FEATURE_DIM)
    with torch.no_grad():
        score = snp_decision_score(model(pileup_features=features, prior_logits=prior))
    assert np.allclose(score.reshape(-1).numpy(), llr, atol=1e-4)


def test_residual_head_can_escape_zero(model):
    """Zero init must not be a dead gradient: one step has to move it."""
    model.train()
    features = torch.randn(2, SEQ_LEN, FEATURE_DIM)
    prior = torch.zeros(2, SEQ_LEN, NUM_CLASSES)
    labels = torch.randint(0, NUM_CLASSES, (2, SEQ_LEN))
    logits, residual = model(pileup_features=features, prior_logits=prior,
                             return_residual=True)
    loss = torch.nn.functional.cross_entropy(
        logits.reshape(-1, NUM_CLASSES), labels.reshape(-1))
    loss.backward()
    grad = model.backbone.classifier.weight.grad
    assert grad is not None and grad.abs().max() > 0


def test_backbone_is_not_zeroed(model):
    """Only the output head is zeroed; the backbone stays randomly initialized."""
    assert model.backbone.input_proj.weight.abs().max() > 0
    penultimate = model.backbone.final_norm(
        model.backbone.input_proj(model.backbone.input_norm(
            torch.randn(2, SEQ_LEN, FEATURE_DIM))))
    assert penultimate.abs().max() > 0


# --- 2. the prior faithfully re-expresses the LLR ----------------------------


def test_prior_places_llr_in_the_snp_channel_only(counts):
    llr = raw_llr_from_counts(counts, BinomialVariantCaller())
    prior = prior_logits_from_llr(llr)
    assert np.array_equal(prior[:, LABEL_SNP], llr)
    for cls in range(NUM_CLASSES):
        if cls != LABEL_SNP:
            assert np.all(prior[:, cls] == 0.0)


def test_prior_snp_minus_normal_gap_is_the_llr(counts):
    llr = raw_llr_from_counts(counts, BinomialVariantCaller())
    prior = prior_logits_from_llr(llr)
    assert np.array_equal(prior[:, LABEL_SNP] - prior[:, LABEL_NORMAL], llr)


def test_indel_channels_carry_no_smuggled_prior(counts):
    """Indel logits equal the Normal logit: the binomial has no indel opinion."""
    prior = prior_logits_from_counts(counts, BinomialVariantCaller())
    assert np.array_equal(prior[:, 2], prior[:, LABEL_NORMAL])
    assert np.array_equal(prior[:, 3], prior[:, LABEL_NORMAL])


def test_decision_score_is_monotone_in_llr(counts):
    """Ranking by decision score == ranking by LLR, so AUC/AP are preserved."""
    llr = raw_llr_from_counts(counts, BinomialVariantCaller())
    score = snp_decision_score(prior_logits_from_llr(llr))
    assert np.array_equal(np.argsort(score, kind="stable"),
                          np.argsort(llr, kind="stable"))


# --- 3. the binomial computation matches the existing implementation ---------


def test_prior_llr_matches_binomial_baseline_score_counts(counts):
    """Exact equality against ``BinomialVariantCaller.score_counts``."""
    caller = BinomialVariantCaller()
    expected = caller.score_counts(counts[:, 0:4], counts[:, 9].astype(int)).llr
    produced = prior_logits_from_counts(counts, caller)[:, LABEL_SNP]
    assert np.array_equal(produced, expected)


def test_prior_llr_matches_llr_features_helper(counts):
    """Exact equality against the helper the previous experiment trained on."""
    caller = BinomialVariantCaller()
    assert np.array_equal(prior_logits_from_counts(counts, caller)[:, LABEL_SNP],
                          raw_llr_from_counts(counts, caller))


def test_prior_uses_v1_parameters_by_default(counts):
    """v1 means fixed eps=0.01 and max(het, hom); v2 must differ, proving the
    default is not silently quality-derived."""
    v1 = BinomialVariantCaller()
    assert v1.error_rate == 0.01 and v1.include_homozygous
    assert not v1.quality_derived_epsilon
    v2 = BinomialVariantCaller(quality_derived_epsilon=True)
    a = prior_logits_from_counts(counts, v1)[:, LABEL_SNP]
    b = prior_logits_from_counts(counts, v2, )[:, LABEL_SNP]
    assert not np.allclose(a, b)


def test_calls_at_threshold_match_the_binomial_caller(counts):
    """The model's calls at a threshold equal the binomial caller's calls."""
    caller = BinomialVariantCaller()
    threshold = 10.5
    expected = caller.call_counts(counts[:, 0:4], counts[:, 9].astype(int), threshold)
    score = snp_decision_score(prior_logits_from_counts(counts, caller))
    assert np.array_equal(score >= threshold, expected)


def test_no_nan_or_inf_in_the_prior(counts):
    prior = prior_logits_from_counts(counts, BinomialVariantCaller())
    assert np.all(np.isfinite(prior))


def test_zero_depth_loci_get_a_neutral_prior():
    """No reads means no evidence: the prior must be flat, not a call."""
    empty = np.zeros((3, 10), dtype=np.float64)
    empty[:, 9] = 1
    prior = prior_logits_from_counts(empty, BinomialVariantCaller())
    assert np.all(prior == 0.0)
