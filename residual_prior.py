"""Binomial evidence expressed as a 4-class prior logit vector.

The residual experiment asks whether the Mamba backbone can *improve on* the
binomial caller rather than relearn it. That question only has a clean answer
if the neural model starts life exactly equal to the binomial caller, which
means the binomial statistic has to be injected as a logit, not as an input
channel (the previous ``llr_features`` experiment did the latter, and the
model promptly distorted the statistic).

Construction
------------
For each locus the binomial LLR

    LLR = max(log P(k|n,p_het), log P(k|n,p_hom)) - log P(k|n,p0)

is a log-likelihood ratio of "alternate allele present" against "sequencing
noise". Under equal priors that *is* the log-odds, so it can be dropped
straight into a softmax as the SNP-versus-Normal logit gap:

    prior[NORMAL]    = 0
    prior[SNP]       = LLR
    prior[INSERTION] = 0
    prior[DELETION]  = 0

Two properties follow, and both are load-bearing:

* ``prior[SNP] - prior[NORMAL] == LLR`` exactly, so the model's SNP decision
  score is the binomial LLR itself and the binomial's frozen threshold
  transfers to the neural model with no recalibration;
* the indel classes are set to the *same* value as Normal rather than a large
  negative number, because the binomial caller has no opinion about indels.
  Pinning them at -20 would have smuggled in a non-binomial prior.

Why the score is read in logit space
------------------------------------
``softmax(prior)[SNP] = sigmoid(LLR - log 3)`` saturates to exactly 1.0 in
float32 above LLR ~= 17, which would tie together most true SNPs and destroy
the ranking that the whole comparison rests on. Every score in this experiment
is therefore the *logit gap* ``logit[SNP] - logit[NORMAL]``, which is exact,
unsaturated, and monotone-equivalent to the probability. See
:func:`snp_decision_score`.

Leakage boundary
----------------
Everything here is a function of pileup counts and the reference base, via
``BinomialVariantCaller``. No VCF, no label, no test statistic.
"""

from __future__ import annotations

import numpy as np
import torch

from binomial_baseline import BinomialVariantCaller
from config import LABEL_NORMAL, LABEL_SNP
from llr_features import raw_llr_from_counts

#: Number of output classes, unchanged from every earlier experiment.
NUM_CLASSES = 4


def prior_logits_from_llr(llr: np.ndarray, num_classes: int = NUM_CLASSES) -> np.ndarray:
    """Expand per-locus LLR values into a ``[N, num_classes]`` prior.

    Args:
        llr: ``[N]`` binomial log-likelihood ratios.
        num_classes: Output width; the SNP column carries the LLR and every
            other column is zero.

    Returns:
        ``[N, num_classes]`` float64 prior logits.
    """
    llr = np.asarray(llr, dtype=np.float64)
    prior = np.zeros((llr.shape[0], num_classes), dtype=np.float64)
    prior[:, LABEL_SNP] = llr
    return prior


def prior_logits_from_counts(
    counts: np.ndarray, caller: BinomialVariantCaller, num_classes: int = NUM_CLASSES
) -> np.ndarray:
    """Compute the binomial prior straight from a ``[N, 10]`` count matrix.

    Args:
        counts: Count matrix following ``pileup_counts.COUNT_COLUMNS``.
        caller: Configured binomial caller (v1 by default in this experiment).
        num_classes: Output width.

    Returns:
        ``[N, num_classes]`` prior logits.
    """
    return prior_logits_from_llr(raw_llr_from_counts(counts, caller), num_classes)


def snp_decision_score(logits: torch.Tensor | np.ndarray) -> torch.Tensor | np.ndarray:
    """The SNP-versus-Normal logit gap used as the decision score everywhere.

    At initialization this equals the binomial LLR exactly, so a threshold
    selected for the binomial caller and a threshold selected for the residual
    model live on the same axis and are directly comparable.

    Args:
        logits: ``[..., num_classes]`` logits, torch or numpy.

    Returns:
        ``[...]`` score of the same type as the input.
    """
    return logits[..., LABEL_SNP] - logits[..., LABEL_NORMAL]


__all__ = [
    "NUM_CLASSES",
    "prior_logits_from_llr",
    "prior_logits_from_counts",
    "snp_decision_score",
]
