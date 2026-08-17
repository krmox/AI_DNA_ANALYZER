"""Per-locus sequencing-error estimation from observed base qualities.

Why
---
``binomial_baseline.BinomialVariantCaller`` scores every locus under a *single*
fixed error probability ``epsilon = 0.01`` (~Q20). The reads already carry a
per-observation Phred score, so the error rate the likelihood assumes need not
be a constant chosen once for the whole genome. This module estimates
``epsilon`` locus by locus from those Phred scores and hands it to the
*unchanged* binomial likelihood, so the only variable that moves between arms
is how ``epsilon`` is obtained.

It also implements the Poisson-binomial formulation, which drops the plug-in
step entirely: if the reads at a locus have error probabilities
``p_1 .. p_n`` then the ALT count is a sum of *non-identically* distributed
Bernoulli variables, and Binomial(n, epsilon) is an approximation to that sum.
The Poisson-binomial arm asks whether collapsing the qualities into one number
was costing information.

Input
-----
Everything is derived from the cached read-level ``uint8`` tensor of
:mod:`read_level_pileup` (``[N, R, READ_DIM]``). That tensor is already the
byte-identical input the read-level neural experiment consumed, and
``reconstruct_counts`` proves it reproduces the aggregate count matrix exactly,
so the two representations describe the same reads by construction. No new
extraction pass is added and no existing cache is modified.

Leakage boundary
----------------
No VCF, no label, and no model output is read anywhere in this module. The
reference base enters only through ``reference_index`` -- already used by the
existing caller to mask the ALT candidate -- and the candidate ALT is chosen by
delegating to :meth:`BinomialVariantCaller.score_counts`, so ALT selection is
bit-identical to the established baseline.

Circularity
-----------
Estimators C and D use *only reference-supporting* reads. This is deliberate.
Estimating ``epsilon`` from the reads that disagree with the reference would
make the statement "these reads look like errors" partly a consequence of
"these reads disagree with the reference", which is the very thing the
likelihood is supposed to test. Reference-supporting reads are also the far
larger sample at low depth, where 1-2 ALT reads carry no usable estimate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import binom

from binomial_baseline import BinomialVariantCaller
from read_level_pileup import FLAG_COUNTED

#: Fixed error probability of the established Binomial v1 arm (~Q20).
FIXED_EPSILON = 0.01

#: Clamp applied to every estimated epsilon.
#:
#: Lower bound: two independent scientific reasons, not numerical convenience.
#: (1) ``log(epsilon / 3)`` must stay finite, and a locus whose reads are all
#: Q40 would otherwise assert a 1-in-10,000 error rate as certain truth.
#: (2) Base-quality scores describe *base miscall* only; they say nothing about
#: misalignment, which does not vanish at high Q and puts a floor on the real
#: rate at which a non-reference base appears for non-variant reasons. 1e-4
#: (Q40) is the most permissive floor still defensible on that second ground.
#: Upper bound: above 0.25 an "error" is more likely than any specific base and
#: the substitution model stops meaning anything.
EPSILON_FLOOR = 1e-4
EPSILON_CEILING = 0.25

#: Minimum reference-supporting reads before the trimmed estimator trims.
#: Below this, trimming 10% of each tail removes nothing or removes the entire
#: sample, so estimator D falls back to the untrimmed mean of the same reads.
MIN_READS_FOR_TRIM = 5

#: Fraction removed from *each* tail by estimator D.
TRIM_FRACTION = 0.1

#: Estimator identifiers accepted by :func:`estimate_epsilon`.
ESTIMATORS: tuple[str, ...] = ("fixed", "mean", "median", "ref_mean", "ref_trimmed")


def phred_to_error(quality: np.ndarray | float) -> np.ndarray:
    """Convert a Phred score to an error probability.

    ``p = 10 ** (-Q / 10)``: Q10 -> 0.1, Q20 -> 0.01, Q30 -> 0.001.

    Args:
        quality: Phred score(s).

    Returns:
        Error probabilit(ies), same shape as the input.
    """
    return np.power(10.0, -np.asarray(quality, dtype=np.float64) / 10.0)


@dataclass
class QualityEvidence:
    """Per-locus read-quality evidence, all arrays ``[N, R]`` unless noted.

    Attributes:
        error_probability: ``10 ** (-Q / 10)`` per read slot.
        counted: True where the slot holds an observation the binomial caller
            counts -- an A/C/G/T base that passed the base-quality filter.
            Gap / deletion reads are False, exactly as they are excluded from
            the substitution denominator in ``pileup_counts``.
        supports_reference: ``counted`` and the base equals the reference base.
        base_code: 0=A 1=C 2=G 3=T per slot.
        n_counted: ``[N]`` number of counted observations (the binomial ``n``).
        n_reference: ``[N]`` number of reference-supporting observations.
    """

    error_probability: np.ndarray
    counted: np.ndarray
    supports_reference: np.ndarray
    base_code: np.ndarray
    n_counted: np.ndarray
    n_reference: np.ndarray


def extract_quality_evidence(reads: np.ndarray, reference_index: np.ndarray) -> QualityEvidence:
    """Decode the compact read tensor into per-read error probabilities.

    Args:
        reads: ``[N, R, READ_DIM]`` uint8 tensor from :mod:`read_level_pileup`.
        reference_index: ``[N]`` reference base, 0=N and 1..4 = A,C,G,T.

    Returns:
        A :class:`QualityEvidence`.
    """
    reads = np.asarray(reads)
    reference_index = np.asarray(reference_index).astype(np.int64)

    flags = reads[..., 6].astype(np.int64)
    counted = (flags & FLAG_COUNTED) > 0
    base_code = reads[..., 0].astype(np.int64)
    error_probability = phred_to_error(reads[..., 1].astype(np.float64))

    # reference_index is 1-based over A,C,G,T; base_code is 0-based. A locus
    # with reference N (index 0) has no reference base, so no read supports it.
    reference_code = reference_index - 1
    supports_reference = counted & (base_code == reference_code[:, None]) \
        & (reference_code >= 0)[:, None]

    return QualityEvidence(
        error_probability=error_probability,
        counted=counted,
        supports_reference=supports_reference,
        base_code=base_code,
        n_counted=counted.sum(axis=1),
        n_reference=supports_reference.sum(axis=1),
    )


def _masked_mean(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Row-wise mean of ``values`` over ``mask``; 0 where the row is empty."""
    total = np.sum(np.where(mask, values, 0.0), axis=1)
    count = mask.sum(axis=1)
    return np.where(count > 0, total / np.maximum(count, 1), 0.0)


def _masked_median(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Row-wise median of ``values`` over ``mask``; 0 where the row is empty.

    Masked-out slots are pushed past every real value by substituting ``inf``
    before the sort, so the first ``m`` sorted entries of a row are exactly its
    ``m`` selected values.
    """
    ordered = np.sort(np.where(mask, values, np.inf), axis=1)
    count = mask.sum(axis=1)
    safe = np.maximum(count, 1)
    lower = (safe - 1) // 2
    upper = safe // 2
    rows = np.arange(values.shape[0])
    median = 0.5 * (ordered[rows, lower] + ordered[rows, upper])
    return np.where(count > 0, median, 0.0)


def _masked_trimmed_mean(values: np.ndarray, mask: np.ndarray,
                         trim: float = TRIM_FRACTION,
                         min_reads: int = MIN_READS_FOR_TRIM) -> np.ndarray:
    """Row-wise symmetric trimmed mean over ``mask``; 0 where the row is empty.

    Rows with fewer than ``min_reads`` selected values are averaged untrimmed.

    Args:
        values: ``[N, R]`` values.
        mask: ``[N, R]`` selection mask.
        trim: Fraction discarded from each tail.
        min_reads: Below this many values, do not trim.

    Returns:
        ``[N]`` trimmed means.
    """
    ordered = np.sort(np.where(mask, values, np.inf), axis=1)
    count = mask.sum(axis=1)
    cut = np.where(count >= min_reads, np.floor(trim * count).astype(np.int64), 0)
    low = cut
    high = count - cut                      # exclusive
    index = np.arange(values.shape[1])[None, :]
    keep = (index >= low[:, None]) & (index < high[:, None])
    kept = keep.sum(axis=1)
    total = np.sum(np.where(keep, ordered, 0.0), axis=1)
    return np.where(kept > 0, total / np.maximum(kept, 1), 0.0)


def estimate_epsilon(
    evidence: QualityEvidence,
    method: str,
    floor: float = EPSILON_FLOOR,
    ceiling: float = EPSILON_CEILING,
    fixed: float = FIXED_EPSILON,
) -> tuple[np.ndarray, dict]:
    """Estimate the per-locus sequencing error probability.

    Estimators
    ----------
    ``fixed``
        The Binomial v1 constant. Present so every arm runs the same code path.
    ``mean``   (A)
        ``mean(p_i)`` over all counted observations. Note this is *not* the
        existing ``binomial_v2``, which computes ``10 ** (-mean(Q) / 10)`` --
        a geometric mean of the ``p_i``. The arithmetic mean is the estimator
        implied by treating the locus as one pooled Bernoulli error process.
    ``median`` (B)
        ``median(p_i)`` over all counted observations; robust to a single
        very-low-quality read dragging the mean up.
    ``ref_mean`` (C)
        ``mean(p_i)`` over reference-supporting observations only. Avoids the
        circularity of estimating the error rate from the reads whose status as
        errors is the question, and uses the larger sample at low depth.
    ``ref_trimmed`` (D)
        Symmetric 10%-trimmed mean over reference-supporting observations.

    Fallbacks (all label-free): a locus with no reads of the required kind
    falls back to the all-read mean, and failing that to ``fixed``. Both events
    are counted in the returned diagnostics.

    Args:
        evidence: Output of :func:`extract_quality_evidence`.
        method: One of :data:`ESTIMATORS`.
        floor: Lower clamp; see :data:`EPSILON_FLOOR`.
        ceiling: Upper clamp.
        fixed: Constant used by ``fixed`` and as the final fallback.

    Returns:
        ``(epsilon [N], diagnostics)``.

    Raises:
        ValueError: If ``method`` is not recognised.
    """
    if method not in ESTIMATORS:
        raise ValueError(f"unknown estimator {method!r}; expected one of {ESTIMATORS}")

    n_loci = evidence.n_counted.shape[0]
    if method == "fixed":
        return np.full(n_loci, float(fixed)), {
            "estimator": "fixed", "fallback_all_read": 0, "fallback_fixed": 0,
            "clamped_low": 0, "clamped_high": 0,
        }

    p = evidence.error_probability
    all_mask = evidence.counted
    reference_mask = evidence.supports_reference

    if method == "mean":
        raw, primary_mask = _masked_mean(p, all_mask), all_mask
    elif method == "median":
        raw, primary_mask = _masked_median(p, all_mask), all_mask
    elif method == "ref_mean":
        raw, primary_mask = _masked_mean(p, reference_mask), reference_mask
    else:  # ref_trimmed
        raw, primary_mask = _masked_trimmed_mean(p, reference_mask), reference_mask

    # Fallback chain: preferred reads -> any counted read -> the fixed constant.
    has_primary = primary_mask.sum(axis=1) > 0
    all_read_mean = _masked_mean(p, all_mask)
    has_any = all_mask.sum(axis=1) > 0

    epsilon = np.where(has_primary, raw, np.where(has_any, all_read_mean, fixed))
    diagnostics = {
        "estimator": method,
        "fallback_all_read": int(np.sum(~has_primary & has_any)),
        "fallback_fixed": int(np.sum(~has_any)),
    }

    clamped = np.clip(epsilon, floor, ceiling)
    diagnostics["clamped_low"] = int(np.sum(epsilon < floor))
    diagnostics["clamped_high"] = int(np.sum(epsilon > ceiling))
    diagnostics["floor"] = floor
    diagnostics["ceiling"] = ceiling

    if not np.all(np.isfinite(clamped)):
        raise ValueError(f"estimator {method!r} produced a non-finite epsilon")
    return clamped, diagnostics


# -- candidate ALT and the unchanged binomial likelihood --------------------


def candidate_alt(base_counts: np.ndarray, reference_index: np.ndarray) -> tuple[
        np.ndarray, np.ndarray, np.ndarray]:
    """Candidate ALT allele, ``k`` and ``n``, delegated to the existing caller.

    Delegation is the point: the candidate-ALT rule must be bit-identical to
    the established baseline, so this experiment cannot accidentally change two
    things at once.

    Args:
        base_counts: ``[N, 4]`` A,C,G,T counts.
        reference_index: ``[N]`` 0=N, 1..4 = A,C,G,T.

    Returns:
        ``(alt_index [N], k [N], n [N])``.
    """
    score = BinomialVariantCaller().score_counts(base_counts, reference_index)
    return score.alt_index, score.k, score.n


def binomial_llr(k: np.ndarray, n: np.ndarray, epsilon: np.ndarray,
                 include_homozygous: bool = True) -> dict[str, np.ndarray]:
    """Binomial log-likelihood ratio at a per-locus epsilon.

    Identical in form to Binomial v1; ``epsilon`` is the only thing that moves::

        p0    = eps / 3
        p_het = 0.5 (1 - eps) + 0.5 (eps / 3)
        p_hom = 1 - eps + eps / 3
        llr   = max(logpmf(k; n, p_het), logpmf(k; n, p_hom)) - logpmf(k; n, p0)

    Args:
        k: ALT-supporting counts.
        n: Usable A/C/G/T observations.
        epsilon: Per-locus error probability.
        include_homozygous: Take the max over het and hom, as v1 does.

    Returns:
        Dict with ``llr``, ``loglik_h0``, ``loglik_het``, ``loglik_hom``.
    """
    epsilon = np.asarray(epsilon, dtype=np.float64)
    p0 = np.clip(epsilon / 3.0, 1e-12, 1 - 1e-12)
    p_het = np.clip(0.5 * (1.0 - epsilon) + 0.5 * (epsilon / 3.0), 1e-12, 1 - 1e-12)
    p_hom = np.clip(1.0 - epsilon + epsilon / 3.0, 1e-12, 1 - 1e-12)

    loglik_h0 = binom.logpmf(k, n, p0)
    loglik_het = binom.logpmf(k, n, p_het)
    loglik_hom = binom.logpmf(k, n, p_hom)

    alternative = np.maximum(loglik_het, loglik_hom) if include_homozygous else loglik_het
    llr = alternative - loglik_h0
    llr = np.where(n <= 0, 0.0, llr)
    llr = np.nan_to_num(llr, nan=0.0, posinf=0.0, neginf=0.0)
    return {"llr": llr, "loglik_h0": loglik_h0,
            "loglik_het": loglik_het, "loglik_hom": loglik_hom}


# -- Poisson-binomial arm ---------------------------------------------------


def _log_poisson_binomial_pmf(q: np.ndarray, active: np.ndarray) -> np.ndarray:
    """Exact log-pmf of a sum of independent, non-identical Bernoulli draws.

    Direct convolution in log space: after processing ``j`` reads, row ``i`` of
    ``dp`` holds ``log P(sum of first j = k)``. The recursion is

        dp'[k] = logaddexp(dp[k] + log(1 - q_j), dp[k-1] + log(q_j))

    which is numerically stable because every term is a log-probability in
    ``[-inf, 0]`` and ``logaddexp`` subtracts the max internally. There is no
    subtraction of nearly-equal quantities anywhere, unlike the FFT/characteristic
    -function methods for the Poisson binomial, which lose precision in the
    tails -- and the tails are exactly where a variant call lives.

    Inactive slots are skipped by carrying ``dp`` forward unchanged.

    Args:
        q: ``[N, R]`` per-read success probabilities.
        active: ``[N, R]`` mask of real observations.

    Returns:
        ``[N, R + 1]`` log-pmf over counts 0..R.
    """
    n_loci, n_slots = q.shape
    dp = np.full((n_loci, n_slots + 1), -np.inf, dtype=np.float64)
    dp[:, 0] = 0.0

    q = np.clip(np.where(active, q, 0.0), 1e-12, 1 - 1e-12)
    log_q = np.log(q)
    log_not_q = np.log1p(-q)

    for slot in range(n_slots):
        stay = dp + log_not_q[:, slot: slot + 1]
        step = np.full_like(dp, -np.inf)
        step[:, 1:] = dp[:, :-1] + log_q[:, slot: slot + 1]
        updated = np.logaddexp(stay, step)
        dp = np.where(active[:, slot: slot + 1], updated, dp)
    return dp


def poisson_binomial_llr(
    evidence: QualityEvidence,
    k: np.ndarray,
    floor: float = EPSILON_FLOOR,
    ceiling: float = EPSILON_CEILING,
    include_homozygous: bool = True,
    chunk: int = 65536,
    alt_index: np.ndarray | None = None,
    legacy_truncation: bool = False,
) -> dict[str, np.ndarray]:
    """LLR under read-specific error probabilities, no plug-in epsilon.

    Each counted read ``i`` carries its own ``p_i = 10 ** (-Q_i / 10)``. The
    probability that read ``i`` shows the candidate ALT base is then

        H0     : p_i / 3
        H1_het : 0.5 (1 - p_i) + 0.5 (p_i / 3)
        H1_hom : 1 - p_i + p_i / 3

    and the ALT count is a Poisson-binomial sum, not a binomial one. When every
    ``p_i`` is equal this reduces exactly to :func:`binomial_llr` -- a property
    the unit tests assert, and the reason this is a strictly weaker assumption
    rather than a different model.

    Scope note: unlike estimators C and D, this arm necessarily uses the ALT
    reads' own qualities, because it models each observation individually.
    That is not circular -- no read's probability depends on whether it matches
    the reference -- but it does mean a genuine variant carried by low-quality
    reads is penalised, which the error analysis should check for.

    Which ``k`` this is a likelihood *of* (devlog 13 fix)
    -----------------------------------------------------
    ``evidence`` describes at most ``read_level_pileup.MAX_READS = 48`` reads
    per locus, because that is the width of the cached read tensor. ``k`` as
    supplied by :func:`candidate_alt` is computed from the *uncapped* pileup
    count matrix. Above ~48x depth the two are on different scales, and the
    pre-fix code reconciled them with ``np.clip(k, 0, tensor_width)``. That is
    not a valid reconciliation:

    * if ``k`` exceeded the locus' number of counted slots, the pmf was indexed
      outside its own support, every hypothesis returned ``-inf``, the LLR came
      out ``nan`` and the trailing ``nan_to_num`` turned the *strongest possible
      evidence for a variant* into exactly ``0.0`` (devlog 12 §15);
    * if it did not, clipping still asserted "48 of 48 retained reads are ALT"
      for a locus where only a fraction of the retained reads were.

    A likelihood must be evaluated on the observation the likelihood's own
    parameters describe, so the ALT count is now taken over the retained reads:

    * ``alt_index`` given (preferred): the exact count of retained, counted
      reads carrying the candidate ALT base;
    * otherwise: ``min(k, n_counted)``, which cannot leave the support and so
      cannot produce the ``nan -> 0`` collapse.

    Both are identities when every read fits the tensor (``depth <= 48``), so
    every result obtained at <=30x is reproduced bit-for-bit.

    Args:
        evidence: Output of :func:`extract_quality_evidence`.
        k: ``[N]`` observed ALT-supporting counts, from the full pileup.
        floor: Lower clamp on each ``p_i``, same rationale as
            :data:`EPSILON_FLOOR`.
        ceiling: Upper clamp on each ``p_i``.
        include_homozygous: Take the max over het and hom, as v1 does.
        chunk: Loci per block, to bound peak memory.
        alt_index: ``[N]`` candidate ALT base code (0=A..3=T) from
            :func:`candidate_alt`. When given, the ALT count is recounted over
            the retained reads rather than taken from ``k``.
        legacy_truncation: Reproduce the pre-fix ``clip(k, 0, tensor_width)``
            behaviour. **Benchmark use only** -- it exists so the fix can be
            A/B-compared on identical reads, and must never be enabled in a
            calling path.

    Returns:
        Dict with ``llr``, ``loglik_h0``, ``loglik_het``, ``loglik_hom``,
        ``k_effective`` (the count actually scored), and the diagnostic
        integers ``n_out_of_support`` (loci where the supplied ``k`` exceeded
        the retained read count) and ``n_non_finite`` (loci whose raw LLR was
        not finite before the ``nan_to_num`` guard).
    """
    p = np.clip(evidence.error_probability, floor, ceiling)
    active = evidence.counted
    k = np.asarray(k).astype(np.int64)
    n_loci = p.shape[0]
    n_counted = np.asarray(evidence.n_counted).astype(np.int64)
    out_of_support = int(np.count_nonzero(k > n_counted))

    if legacy_truncation:
        k_effective = np.clip(k, 0, p.shape[1])
    elif alt_index is None:
        k_effective = np.clip(k, 0, n_counted)
    else:
        alt_index = np.asarray(alt_index).astype(np.int64)
        k_effective = np.sum(active & (evidence.base_code == alt_index[:, None]),
                             axis=1).astype(np.int64)

    out = {name: np.zeros(n_loci, dtype=np.float64)
           for name in ("loglik_h0", "loglik_het", "loglik_hom")}

    for start in range(0, n_loci, chunk):
        stop = min(start + chunk, n_loci)
        block_p = p[start:stop]
        block_active = active[start:stop]
        rows = np.arange(stop - start)
        block_k = k_effective[start:stop]

        for name, q in (
            ("loglik_h0", block_p / 3.0),
            ("loglik_het", 0.5 * (1.0 - block_p) + 0.5 * (block_p / 3.0)),
            ("loglik_hom", 1.0 - block_p + block_p / 3.0),
        ):
            out[name][start:stop] = _log_poisson_binomial_pmf(q, block_active)[rows, block_k]

    alternative = np.maximum(out["loglik_het"], out["loglik_hom"]) if include_homozygous \
        else out["loglik_het"]
    with np.errstate(invalid="ignore"):
        llr = alternative - out["loglik_h0"]
    llr = np.where(evidence.n_counted <= 0, 0.0, llr)
    non_finite = int(np.count_nonzero(~np.isfinite(llr)))
    llr = np.nan_to_num(llr, nan=0.0, posinf=0.0, neginf=0.0)
    return {"llr": llr, "k_effective": k_effective,
            "n_out_of_support": out_of_support, "n_non_finite": non_finite, **out}


__all__ = [
    "FIXED_EPSILON", "EPSILON_FLOOR", "EPSILON_CEILING", "ESTIMATORS",
    "MIN_READS_FOR_TRIM", "TRIM_FRACTION",
    "QualityEvidence", "phred_to_error", "extract_quality_evidence",
    "estimate_epsilon", "candidate_alt", "binomial_llr", "poisson_binomial_llr",
]
