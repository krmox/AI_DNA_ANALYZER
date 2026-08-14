"""Label-free routing between the Poisson-binomial caller and the PB-residual Mamba.

The architecture under test
---------------------------
::

    reads -> pileup evidence -> Poisson-binomial LLR
                                      |
                              uncertainty signal u
                                      |
                    +-----------------+-----------------+
              u < cutoff (easy)                  u >= cutoff (hard)
                    |                                   |
              PB call: llr >= T_pb          Mamba call: s_M >= T_M
                    +-----------------+-----------------+
                                      v
                                 final call

This is *Design A* of the experiment brief: the PB LLR is computed for every
locus and the router reads it. Design B (a cheap pre-router that fires before
the PB prior is paid for) is deliberately not built here; it is gated on Design
A showing a useful accuracy/compute trade-off.

Why the two arms need two thresholds
------------------------------------
PB and the residual model share an axis -- the model's score is
``PB_LLR + residual`` -- but each arm carries the threshold *its own*
validation selection chose (PB 10.5; the seeds 7.0/9.5/7.5). Mixing arms
therefore means mixing operating points, not scores. The hybrid never
interpolates scores; it picks which arm's *decision* to keep. That keeps each
arm exactly as validated and makes the hybrid a pure routing question.

Leakage boundary
----------------
Every function here takes evidence (features, counts, depth) and caller output
(the PB LLR, its frozen threshold) and nothing else. No label, no truth VCF, no
position, no model output from the arm being routed to. ``truth`` never appears
in a signature in this module, which is the property ``test_hybrid_router``
asserts by permuting labels and checking the routing is bit-identical.

The routing *cutoff* is a quantile of the validation uncertainty distribution.
That uses validation loci -- allowed, and it uses only their evidence, not
their labels. The *choice* of signal and coverage does consult validation
labels, which is the same budget every threshold in this project is selected
from, and is frozen before test scoring.
"""

from __future__ import annotations

import numpy as np

#: Pre-registered Mamba coverage grid, as fractions of loci routed to the
#: neural arm. Fixed before any test locus is scored.
COVERAGE_GRID: tuple[float, ...] = (
    0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.0,
)

#: Coverage at which the primary routing signal is chosen. Pre-registered so
#: the signal is not selected by scanning the whole curve for a winner.
SELECTION_COVERAGE: float = 0.10

#: A configuration is preferred over a larger one when it comes within this
#: much validation F1 of the best point on the curve. Pre-registered.
KNEE_TOLERANCE: float = 0.0002


def _column(features: np.ndarray, feature_names: list[str], name: str) -> np.ndarray:
    """One named evidence column.

    Args:
        features: ``[N, F]`` locus evidence.
        feature_names: Column names, in order.
        name: Column to extract.

    Returns:
        ``[N]`` float64 column.
    """
    return np.asarray(features[:, list(feature_names).index(name)], dtype=np.float64)


def pb_margin_uncertainty(llr: np.ndarray, pb_threshold: float) -> np.ndarray:
    """Signal A: closeness to the PB decision boundary.

    The hypothesis under test is that the loci the neural residual can move are
    the ones PB is nearly indifferent about. Returned as a *negated* distance
    so that larger always means "more uncertain", the convention every signal
    in this module follows.

    Args:
        llr: ``[N]`` Poisson-binomial log-likelihood ratios.
        pb_threshold: The validation-frozen PB operating point.

    Returns:
        ``[N]`` uncertainty, higher meaning closer to the boundary.
    """
    return -np.abs(np.asarray(llr, dtype=np.float64) - float(pb_threshold))


def pb_posterior_uncertainty(llr: np.ndarray, log_prior_odds: float = np.log(3.0)) -> np.ndarray:
    """Signal B: binary entropy of the caller's own posterior.

    Distinct from :func:`pb_margin_uncertainty` rather than a monotone
    relabelling of it: the entropy peaks where the caller's *probability* is
    0.5, i.e. at ``llr == log 3`` under the four-class prior of
    ``residual_prior``, whereas the margin peaks at the F1-maximizing threshold
    of 10.5. The two disagree about which loci are ambiguous, and the
    experiment measures which disagreement is right.

    Args:
        llr: ``[N]`` Poisson-binomial log-likelihood ratios.
        log_prior_odds: Offset turning the LLR into the SNP log-odds, matching
            ``softmax(prior)[SNP] = sigmoid(LLR - log 3)``.

    Returns:
        ``[N]`` binary entropy in nats, higher meaning less decided.
    """
    logit = np.asarray(llr, dtype=np.float64) - log_prior_odds
    probability = 1.0 / (1.0 + np.exp(-np.clip(logit, -60.0, 60.0)))
    p = np.clip(probability, 1e-12, 1.0 - 1e-12)
    return -(p * np.log(p) + (1.0 - p) * np.log(1.0 - p))


def evidence_count_uncertainty(features: np.ndarray, feature_names: list[str],
                               vaf_threshold: float) -> np.ndarray:
    """Signal C: closeness to the naive frequency caller's boundary.

    The count-only view of ambiguity: a locus whose VAF sits at the frequency
    caller's validated operating point is one where read counts alone cannot
    decide. Uses ``vaf`` from the evidence matrix, which is ``k / n`` with the
    integer semantics the 44-channel representation restored.

    Args:
        features: ``[N, F]`` locus evidence.
        feature_names: Column names, in order.
        vaf_threshold: Frequency-caller threshold frozen on validation.

    Returns:
        ``[N]`` uncertainty, higher meaning closer to the VAF boundary.
    """
    return -np.abs(_column(features, feature_names, "vaf") - float(vaf_threshold))


def low_depth_uncertainty(depth: np.ndarray) -> np.ndarray:
    """Signal C': shallow loci first.

    The crudest evidence-count router, included as a control: if simply sending
    the shallowest loci to the neural arm works as well as the score margin,
    the margin hypothesis has not earned anything.

    Args:
        depth: ``[N]`` per-locus depth.

    Returns:
        ``[N]`` uncertainty, higher meaning shallower.
    """
    return -np.asarray(depth, dtype=np.float64)


def alt_quality_uncertainty(features: np.ndarray, feature_names: list[str]) -> np.ndarray:
    """Signal D: ambiguity of the ALT-supporting reads' base quality.

    Experiment 8 established that ALT-read quality is the decisive variable for
    the statistical caller. The ambiguous case is not "bad quality" but *mixed*
    quality -- some ALT reads clearing Q30 and some not -- because a uniformly
    good or uniformly bad ALT pile is a locus the error model already handles.
    Mixedness is measured as ``1 - |2 f - 1|`` on the Q30 fraction, which is 1
    at a 50/50 split and 0 at either extreme, and is scaled up by the mean ALT
    error probability so that a coin-flip split of poor reads outranks a
    coin-flip split of excellent ones. Loci with no ALT read get 0: there is no
    quality evidence there to be ambiguous about.

    Args:
        features: ``[N, F]`` locus evidence.
        feature_names: Column names, in order.

    Returns:
        ``[N]`` uncertainty, higher meaning more mixed ALT quality evidence.
    """
    q30 = _column(features, feature_names, "alt_fraction_q30")
    error = _column(features, feature_names, "alt_mean_error")
    has_alt = _column(features, feature_names, "alt_has_reads")
    mixedness = 1.0 - np.abs(2.0 * q30 - 1.0)
    return has_alt * mixedness * (1.0 + error)


def _standardize(values: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """z-score ``values`` using ``reference``'s location and scale."""
    reference = np.asarray(reference, dtype=np.float64)
    centre = float(np.median(reference))
    scale = float(np.median(np.abs(reference - centre))) or 1.0
    return (np.asarray(values, dtype=np.float64) - centre) / scale


def composite_uncertainty(llr: np.ndarray, pb_threshold: float, features: np.ndarray,
                          feature_names: list[str],
                          reference: dict[str, np.ndarray] | None = None) -> np.ndarray:
    """Signal E: score margin plus ALT-quality ambiguity, equally weighted.

    Deliberately the smallest composite that can be written: two z-scored terms
    added, no fitted coefficients, so it stays interpretable and cannot quietly
    become a second model trained on validation labels.

    Args:
        llr: ``[N]`` Poisson-binomial LLRs.
        pb_threshold: Frozen PB operating point.
        features: ``[N, F]`` locus evidence.
        feature_names: Column names, in order.
        reference: Optional ``{"margin": ..., "quality": ...}`` arrays supplying
            the standardization statistics, so test loci are scaled by
            validation statistics rather than their own.

    Returns:
        ``[N]`` composite uncertainty.
    """
    margin = pb_margin_uncertainty(llr, pb_threshold)
    quality = alt_quality_uncertainty(features, feature_names)
    reference = reference or {"margin": margin, "quality": quality}
    return _standardize(margin, reference["margin"]) + _standardize(quality, reference["quality"])


def all_signals(llr: np.ndarray, pb_threshold: float, features: np.ndarray,
                feature_names: list[str], depth: np.ndarray, vaf_threshold: float,
                reference: dict[str, np.ndarray] | None = None) -> dict[str, np.ndarray]:
    """Every candidate routing signal, as a dict of ``[N]`` uncertainties.

    The candidate family is closed here, before any test evaluation, exactly as
    the brief requires.

    Args:
        llr: ``[N]`` Poisson-binomial LLRs.
        pb_threshold: Frozen PB operating point.
        features: ``[N, F]`` locus evidence.
        feature_names: Column names, in order.
        depth: ``[N]`` per-locus depth.
        vaf_threshold: Frozen frequency-caller operating point.
        reference: Optional standardization reference for the composite.

    Returns:
        Signal name -> ``[N]`` uncertainty, higher meaning route to the neural
        arm sooner.
    """
    return {
        "pb_margin": pb_margin_uncertainty(llr, pb_threshold),
        "pb_posterior_entropy": pb_posterior_uncertainty(llr),
        "vaf_margin": evidence_count_uncertainty(features, feature_names, vaf_threshold),
        "low_depth": low_depth_uncertainty(depth),
        "alt_quality": alt_quality_uncertainty(features, feature_names),
        "composite": composite_uncertainty(llr, pb_threshold, features, feature_names,
                                           reference),
    }


def cutoff_for_coverage(uncertainty: np.ndarray, coverage: float) -> float:
    """Uncertainty cutoff routing ``coverage`` of these loci to the neural arm.

    Args:
        uncertainty: ``[N]`` uncertainties from the split the cutoff is fitted
            on (validation).
        coverage: Target fraction in ``[0, 1]``.

    Returns:
        Cutoff ``c`` such that ``uncertainty >= c`` selects about ``coverage``
        of the loci. ``+inf`` at coverage 0 and ``-inf`` at coverage 1, so the
        endpoints route nothing and everything exactly, whatever ties exist.
    """
    if coverage <= 0.0:
        return float("inf")
    if coverage >= 1.0:
        return float("-inf")
    return float(np.quantile(np.asarray(uncertainty, dtype=np.float64), 1.0 - coverage))


def route(uncertainty: np.ndarray, cutoff: float) -> np.ndarray:
    """Boolean mask of loci sent to the neural arm.

    Args:
        uncertainty: ``[N]`` uncertainties.
        cutoff: Frozen cutoff from :func:`cutoff_for_coverage`.

    Returns:
        ``[N]`` boolean mask, True meaning "route to Mamba".
    """
    return np.asarray(uncertainty, dtype=np.float64) >= cutoff


def hybrid_calls(pb_score: np.ndarray, pb_threshold: float, neural_score: np.ndarray,
                 neural_threshold: float, routed: np.ndarray) -> np.ndarray:
    """Combine the two arms' decisions under a routing mask.

    Args:
        pb_score: ``[N]`` PB LLRs.
        pb_threshold: PB's frozen threshold.
        neural_score: ``[N]`` PB-residual model scores.
        neural_threshold: The checkpoint's own validation-frozen threshold.
        routed: ``[N]`` boolean routing mask.

    Returns:
        ``[N]`` boolean calls: the neural arm's decision where routed, the
        statistical arm's decision elsewhere.
    """
    routed = np.asarray(routed, dtype=bool)
    return np.where(routed, np.asarray(neural_score) >= neural_threshold,
                    np.asarray(pb_score) >= pb_threshold)


def routing_enrichment(routed: np.ndarray, pb_calls: np.ndarray,
                       truth: np.ndarray) -> dict:
    """How much of the statistical caller's error the router actually captures.

    Uses labels and is therefore only ever called on validation loci, or on
    test loci *after* the frozen evaluation, as reporting rather than
    selection. It answers the brief's second hypothesis: routing 10% of loci
    and catching 10% of the errors means the router does nothing.

    Args:
        routed: ``[N]`` routing mask.
        pb_calls: ``[N]`` boolean PB decisions at its frozen threshold.
        truth: ``[N]`` boolean truth.

    Returns:
        Routed fraction, captured error fractions, and the enrichment ratio
        ``P(PB error | routed) / P(PB error)``.
    """
    routed = np.asarray(routed, dtype=bool)
    truth = np.asarray(truth, dtype=bool)
    false_positive = pb_calls & ~truth
    false_negative = ~pb_calls & truth
    error = false_positive | false_negative

    routed_fraction = float(routed.mean())
    base_rate = float(error.mean())
    captured = float((error & routed).sum() / max(error.sum(), 1))
    conditional = float((error & routed).sum() / max(routed.sum(), 1))
    return {
        "routed_fraction": routed_fraction,
        "pb_errors": int(error.sum()),
        "pb_fp": int(false_positive.sum()),
        "pb_fn": int(false_negative.sum()),
        "fp_routed": int((false_positive & routed).sum()),
        "fn_routed": int((false_negative & routed).sum()),
        "fp_captured_fraction": float((false_positive & routed).sum()
                                      / max(false_positive.sum(), 1)),
        "fn_captured_fraction": float((false_negative & routed).sum()
                                      / max(false_negative.sum(), 1)),
        "error_captured_fraction": captured,
        "enrichment": float(conditional / base_rate) if base_rate > 0 else float("nan"),
        "snp_routed_fraction": float(routed[truth].mean()) if truth.any() else 0.0,
    }


def window_coverage(routed: np.ndarray, seq_len: int) -> dict:
    """Fraction of *windows* the neural arm must still run.

    The load-bearing systems fact of this experiment: the Mamba consumes
    windows of ``seq_len`` consecutive loci, so a locus routed in isolation
    costs a whole window of neural compute. Locus coverage is the accuracy
    variable; window coverage is the compute variable, and they are not the
    same number.

    Args:
        routed: ``[N]`` routing mask, N a multiple of ``seq_len``.
        seq_len: Window width in loci.

    Returns:
        Locus and window coverage plus the resulting amplification factor.
    """
    routed = np.asarray(routed, dtype=bool)
    windows = routed.reshape(-1, seq_len).any(axis=1)
    locus_fraction = float(routed.mean())
    window_fraction = float(windows.mean())
    return {
        "locus_fraction": locus_fraction,
        "window_fraction": window_fraction,
        "windows_total": int(windows.size),
        "windows_routed": int(windows.sum()),
        "amplification": float(window_fraction / locus_fraction)
        if locus_fraction > 0 else float("nan"),
    }


__all__ = [
    "COVERAGE_GRID", "SELECTION_COVERAGE", "KNEE_TOLERANCE",
    "pb_margin_uncertainty", "pb_posterior_uncertainty", "evidence_count_uncertainty",
    "low_depth_uncertainty", "alt_quality_uncertainty", "composite_uncertainty",
    "all_signals", "cutoff_for_coverage", "route", "hybrid_calls",
    "routing_enrichment", "window_coverage",
]
