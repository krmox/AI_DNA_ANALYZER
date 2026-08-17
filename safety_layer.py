"""Cheap safety layer: escalate loci where the fixed-epsilon binomial is *confidently wrong*.

Why this module exists
----------------------
Devlog 14 §17.3 established the frozen cascade's one structural failure mode:

    The frozen cascade departs from PB exactly at loci where the fixed-epsilon
    binomial is confidently wrong -- and the router, which measures
    *uncertainty* rather than *error*, cannot see them by construction.

The frozen router routes on ``|binomial_llr - 7.0| <= 5.412``. A locus where the
binomial is confidently wrong has, by definition, a *large* margin, so widening
the cutoff cannot reach it (devlog 14 §13 measured that directly). Detecting the
regime needs a signal that is not a function of the binomial's own confidence.

The signal
----------
The binomial caller's confidence at these loci is an artefact of one fixed
constant: ``epsilon = 0.01``. Devlog 14's failure dump shows both failure
directions are epsilon misspecification --

* below ~Q25 the true error rate exceeds 0.01, so the caller over-calls
  (``router_fp``);
* above ~Q30 it is far below 0.01, so a handful of high-quality non-reference
  reads is stronger evidence than the caller credits, and it under-calls
  (``router_fn`` -- the five lost true SNP).

So the safety layer asks a question the router cannot: **is this call robust to
plausible misspecification of epsilon?** It re-runs the *same* binomial
likelihood ratio at two bracketing error rates and escalates when the call flips:

    escalate  <=>  (llr(eps_lo) >= 7.0) != (llr(eps_hi) >= 7.0)  AND  k >= k_min
                   AND the frozen router did not already route the locus

``k >= k_min`` requires a minimum of raw alternate-allele support, which removes
flips driven by a single read without changing which true variants are caught
(devlog 15 §9).

Cost
----
Two extra ``scipy.stats.binom.logpmf`` evaluations over the ``(k, n)`` the
binomial caller has already computed. No pileup is re-read, no read tensor is
built, no per-read quality is touched, and PB is never consulted to decide
whether PB should run. Measured at ~2x the binomial caller and ~1/175 of PB per
locus (devlog 15 §12), so running it on *every* locus costs well under 1% of the
cascade's wall clock.

Frozen boundary
---------------
:data:`FROZEN_BINOMIAL_THRESHOLD` and :data:`FROZEN_ROUTER_CUTOFF` are imported
from :mod:`cascade` and only ever read. This module adds a stage *after* the
router and never alters what the router does: :func:`safety_escalate` takes the
router's own mask as an input and returns a mask disjoint from it.

Leakage boundary
----------------
Every function here takes caller statistics and nothing else. No truth label, no
VCF, no genomic position and no PB score ever appears in a signature -- enforced
by ``test_safety_layer.py`` via signature introspection, the same pattern
``test_cheap_router.py`` and ``test_cascade.py`` use.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import binom

from cascade import FROZEN_BINOMIAL_THRESHOLD, FROZEN_ROUTER_CUTOFF, route_stage1

#: Candidate grid searched on validation data only (devlog 15 §7). The frozen
#: choice is :data:`FROZEN_SAFETY_RULE`; nothing else in the grid is production.
CANDIDATE_EPSILON_BRACKETS = (
    (1e-3, 1e-2), (1e-3, 3e-2), (1e-3, 5e-2),
    (1e-4, 1e-2), (1e-4, 3e-2), (1e-4, 5e-2), (1e-4, 1e-1),
    (3e-3, 3e-2), (1e-2, 5e-2),
)
CANDIDATE_K_MIN = (1, 2, 3, 4)

#: Frozen by the devlog 15 §7.4 selection rule on the validation cells, before
#: any test locus was extracted or scored. See
#: ``results/bench_v15/rule_selection.json`` for the full 36-candidate grid this
#: was chosen from and the criterion that chose it.
#:
#: FULL is the minimum-budget candidate that loses **zero** true SNP relative to
#: PB-only on validation; LEAN is the minimum-budget candidate allowed to lose
#: one. LEAN happens to fall in family Sq (the quality-derived-epsilon straddle)
#: rather than family SK, so both code paths below are production paths.
FROZEN_SAFETY_RULE: dict = {
    "family": "SK",
    "name": "SK[0.0001,0.05]&k>=3",
    "epsilon_lo": 1e-4,
    "epsilon_hi": 5e-2,
    "k_min": 3,
}

FROZEN_SAFETY_RULE_LEAN: dict = {
    "family": "Sq",
    "name": "Sq(floor=0.001)",
    "quality_floor": 1e-3,
    "k_min": 1,
}


def binomial_llr_at_epsilon(k: np.ndarray, n: np.ndarray, epsilon: float) -> np.ndarray:
    """The project's binomial LLR recomputed at an arbitrary fixed error rate.

    This is ``BinomialVariantCaller.score_counts``'s likelihood ratio with
    ``(k, n)`` supplied directly rather than re-derived from the count matrix,
    which is what makes the safety layer cheap: the caller has already chosen
    the candidate alternate allele and counted it, so the only remaining work is
    three ``logpmf`` evaluations.

    The three hypotheses, the ``max(het, hom)`` alternative and the ``n == 0``
    convention are identical to :class:`binomial_baseline.BinomialVariantCaller`;
    ``test_safety_layer.py`` asserts agreement with it to floating-point
    tolerance rather than trusting this comment.

    Args:
        k: ``[N]`` observations supporting the candidate alternate base.
        n: ``[N]`` usable A/C/G/T observations.
        epsilon: Per-base error probability to score under, in ``(0, 1)``.

    Returns:
        ``[N]`` float64 log-likelihood ratios.

    Raises:
        ValueError: If ``epsilon`` is outside ``(0, 1)`` or the shapes differ.
    """
    if not 0.0 < float(epsilon) < 1.0:
        raise ValueError("epsilon must lie in (0, 1)")
    k = np.asarray(k, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    if k.shape != n.shape:
        raise ValueError(f"k and n must have the same shape, got {k.shape} and {n.shape}")

    p0 = np.clip(epsilon / 3.0, 1e-12, 1 - 1e-12)
    p_het = np.clip(0.5 * (1.0 - epsilon) + 0.5 * (epsilon / 3.0), 1e-12, 1 - 1e-12)
    p_hom = np.clip(1.0 - epsilon + epsilon / 3.0, 1e-12, 1 - 1e-12)

    llr = np.maximum(binom.logpmf(k, n, p_het), binom.logpmf(k, n, p_hom)) \
        - binom.logpmf(k, n, p0)
    llr = np.where(n <= 0, 0.0, llr)
    return np.nan_to_num(llr, nan=0.0, posinf=0.0, neginf=0.0)


def epsilon_flip(k: np.ndarray, n: np.ndarray, epsilon_lo: float,
                 epsilon_hi: float) -> np.ndarray:
    """Does the binomial's SNP call change between two bracketing error rates?

    Args:
        k: ``[N]`` alternate-allele support.
        n: ``[N]`` usable observations.
        epsilon_lo: Optimistic error rate (the low end of the bracket).
        epsilon_hi: Pessimistic error rate (the high end).

    Returns:
        ``[N]`` boolean mask, True where the two calls differ.

    Raises:
        ValueError: If ``epsilon_lo >= epsilon_hi``.
    """
    if not float(epsilon_lo) < float(epsilon_hi):
        raise ValueError("epsilon_lo must be strictly less than epsilon_hi")
    call_lo = binomial_llr_at_epsilon(k, n, epsilon_lo) >= FROZEN_BINOMIAL_THRESHOLD
    call_hi = binomial_llr_at_epsilon(k, n, epsilon_hi) >= FROZEN_BINOMIAL_THRESHOLD
    return call_lo != call_hi


def quality_epsilon(n: np.ndarray, quality_sum: np.ndarray, floor: float) -> np.ndarray:
    """The v2 caller's per-locus error rate: ``10**(-meanQ/10)``, clamped.

    Identical arithmetic to ``BinomialVariantCaller._epsilon`` with
    ``quality_derived_epsilon=True`` — including its ``quality_sum / n``
    denominator and its ``[floor, 0.25]`` clamp — and asserted equal to it in
    ``test_safety_layer.py``. It is restated here only so the safety layer can
    take ``(n, quality_sum)`` directly instead of a count matrix.

    Args:
        n: ``[N]`` usable A/C/G/T observations.
        quality_sum: ``[N]`` summed Phred of those observations.
        floor: Lower clamp on the returned probability.

    Returns:
        ``[N]`` per-locus error probabilities in ``[floor, 0.25]``.
    """
    n = np.asarray(n, dtype=np.float64)
    quality_sum = np.asarray(quality_sum, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_quality = np.where(n > 0, quality_sum / np.maximum(n, 1), 0.0)
    epsilon = np.power(10.0, -mean_quality / 10.0)
    epsilon = np.where(n > 0, epsilon, 0.01)
    return np.clip(epsilon, floor, 0.25)


def binomial_llr_per_locus_epsilon(k: np.ndarray, n: np.ndarray,
                                   epsilon: np.ndarray) -> np.ndarray:
    """:func:`binomial_llr_at_epsilon` with a *per-locus* epsilon array.

    Args:
        k: ``[N]`` alternate-allele support.
        n: ``[N]`` usable observations.
        epsilon: ``[N]`` per-locus error probabilities.

    Returns:
        ``[N]`` float64 log-likelihood ratios.
    """
    k = np.asarray(k, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    epsilon = np.asarray(epsilon, dtype=np.float64)
    p0 = np.clip(epsilon / 3.0, 1e-12, 1 - 1e-12)
    p_het = np.clip(0.5 * (1.0 - epsilon) + 0.5 * (epsilon / 3.0), 1e-12, 1 - 1e-12)
    p_hom = np.clip(1.0 - epsilon + epsilon / 3.0, 1e-12, 1 - 1e-12)
    llr = np.maximum(binom.logpmf(k, n, p_het), binom.logpmf(k, n, p_hom)) \
        - binom.logpmf(k, n, p0)
    llr = np.where(n <= 0, 0.0, llr)
    return np.nan_to_num(llr, nan=0.0, posinf=0.0, neginf=0.0)


def quality_flip(k: np.ndarray, n: np.ndarray, quality_sum: np.ndarray,
                 floor: float) -> np.ndarray:
    """Family Sq: does the call differ between fixed ε = 0.01 and quality-derived ε?

    Args:
        k: ``[N]`` alternate-allele support.
        n: ``[N]`` usable observations.
        quality_sum: ``[N]`` summed Phred of the counted observations.
        floor: Lower clamp on the quality-derived epsilon.

    Returns:
        ``[N]`` boolean mask, True where the two calls differ.
    """
    fixed = binomial_llr_at_epsilon(k, n, 0.01) >= FROZEN_BINOMIAL_THRESHOLD
    derived = binomial_llr_per_locus_epsilon(
        k, n, quality_epsilon(n, quality_sum, floor)) >= FROZEN_BINOMIAL_THRESHOLD
    return fixed != derived


def safety_escalate(k: np.ndarray, n: np.ndarray, routed_to_pb: np.ndarray,
                    rule: dict | None = None,
                    quality_sum: np.ndarray | None = None) -> np.ndarray:
    """Which *unrouted* loci the safety layer sends to PB anyway.

    Args:
        k: ``[N]`` alternate-allele support, as counted by the binomial caller.
        n: ``[N]`` usable A/C/G/T observations.
        routed_to_pb: ``[N]`` the frozen router's own mask. Loci it already
            routed are never escalated again, so the returned mask is disjoint
            from it and the router's behaviour is untouched.
        rule: Rule parameters; defaults to :data:`FROZEN_SAFETY_RULE`. A rule
            carrying ``epsilon_lo``/``epsilon_hi`` uses the ε-bracket flip; one
            carrying ``quality_floor`` uses the quality-derived-ε flip.
        quality_sum: ``[N]`` summed Phred, required only by a
            ``quality_floor`` rule.

    Returns:
        ``[N]`` boolean mask, True meaning "escalate this locus to PB".

    Raises:
        ValueError: If the rule is a quality rule and ``quality_sum`` is absent,
            or if the rule carries neither parameterisation.
    """
    rule = FROZEN_SAFETY_RULE if rule is None else rule
    k = np.asarray(k, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    routed_to_pb = np.asarray(routed_to_pb, dtype=bool)

    # Evaluate the flip test only where it can possibly fire. A locus with
    # k < k_min is excluded by the rule outright, and a locus with k == 0 can
    # never produce a positive call at any epsilon (its LLR is negative for
    # every hypothesis), so neither can flip. Restricting the two logpmf passes
    # to the survivors is a pure evaluation-order change: the returned mask is
    # bit-identical to the dense computation, asserted in test_safety_layer.py.
    gate = (k >= max(float(rule["k_min"]), 1.0)) & ~routed_to_pb
    index = np.flatnonzero(gate)
    escalate = np.zeros(k.shape, dtype=bool)
    if index.size == 0:
        return escalate

    k_gated, n_gated = k[index], n[index]
    if "epsilon_lo" in rule and "epsilon_hi" in rule:
        flip = epsilon_flip(k_gated, n_gated, rule["epsilon_lo"], rule["epsilon_hi"])
    elif "quality_floor" in rule:
        if quality_sum is None:
            raise ValueError("a quality_floor rule needs quality_sum")
        flip = quality_flip(k_gated, n_gated,
                            np.asarray(quality_sum, dtype=np.float64)[index],
                            rule["quality_floor"])
    else:
        raise ValueError(f"rule carries no usable parameterisation: {sorted(rule)}")

    escalate[index[flip]] = True
    assert not np.any(escalate & routed_to_pb), "escalation must be disjoint from routing"
    return escalate


def guarded_route(binomial_llr: np.ndarray, k: np.ndarray, n: np.ndarray,
                  rule: dict | None = None,
                  quality_sum: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """The full Stage-1 + safety decision: who goes to PB, and why.

    Args:
        binomial_llr: ``[N]`` fixed-epsilon (0.01) binomial LLR -- the router's
            own input, unchanged.
        k: ``[N]`` alternate-allele support.
        n: ``[N]`` usable observations.
        rule: Rule parameters; defaults to :data:`FROZEN_SAFETY_RULE`.
        quality_sum: ``[N]`` summed Phred, required only by a
            ``quality_floor`` rule.

    Returns:
        ``(to_pb, escalated)``: ``to_pb`` is the union of the frozen router's
        mask and the safety escalations; ``escalated`` is the safety layer's
        contribution alone.
    """
    routed = route_stage1(binomial_llr)
    escalated = safety_escalate(k, n, routed, rule, quality_sum)
    return routed | escalated, escalated


__all__ = [
    "CANDIDATE_EPSILON_BRACKETS", "CANDIDATE_K_MIN",
    "FROZEN_SAFETY_RULE", "FROZEN_SAFETY_RULE_LEAN",
    "FROZEN_BINOMIAL_THRESHOLD", "FROZEN_ROUTER_CUTOFF",
    "binomial_llr_at_epsilon", "binomial_llr_per_locus_epsilon", "quality_epsilon",
    "epsilon_flip", "quality_flip", "safety_escalate", "guarded_route",
]
