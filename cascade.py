"""Three-stage cascade: binomial -> [frozen router] -> PB -> [margin router] -> Mamba.

The architecture under test
----------------------------
::

    binomial LLR --> |frozen cheap router| --> confident: binomial answer
                                             --> uncertain: PB LLR
                                                              |
                                                    PB margin = |PB_LLR - FROZEN_PB_THRESHOLD|
                                                              |
                                     margin >= pb_to_mamba_threshold: PB answer (confident)
                                     margin <  pb_to_mamba_threshold: Mamba answer (final)

Stage 1 is exactly the frozen pre-router from ``cheap_router``/``History/13_DEVLOG.md``
and ``History/14_DEVLOG.md``: the three constants below are never refit here. Stage 2
is new: a second margin threshold, this time on the *PB* LLR rather than the binomial
LLR, deciding which of the PB-routed loci are uncertain enough to hand to the neural
residual model.

Leakage boundary
-----------------
Every function here takes evidence and caller scores and nothing else. No truth
label, no VCF, no position, ever appears in a routing function's signature --
enforced by ``test_cascade.py`` via signature introspection, following the same
pattern ``test_cheap_router.py`` and ``test_hybrid_router.py`` use.
"""

from __future__ import annotations

import numpy as np

#: Frozen pre-router constants, from ``cheap_router``/devlog 13. Never refit here.
FROZEN_BINOMIAL_THRESHOLD = 7.0
FROZEN_PB_THRESHOLD = 10.5
FROZEN_ROUTER_CUTOFF = 5.411872376933351

#: Stage names, in cascade order. Every locus's final call is attributed to
#: exactly one of these.
STAGE_BINOMIAL = "binomial"
STAGE_PB = "pb"
STAGE_MAMBA = "mamba"


def route_stage1(binomial_llr: np.ndarray) -> np.ndarray:
    """Frozen pre-router: which loci are uncertain enough to deserve PB.

    Reproduces ``|binomial_LLR - FROZEN_BINOMIAL_THRESHOLD| <= FROZEN_ROUTER_CUTOFF``
    verbatim from ``cheap_router``/devlog 13-14. Nothing here is fitted.

    Args:
        binomial_llr: ``[N]`` fixed-epsilon binomial log-likelihood ratios.

    Returns:
        ``[N]`` boolean mask, True meaning "route to PB".
    """
    binomial_llr = np.asarray(binomial_llr, dtype=np.float64)
    return np.abs(binomial_llr - FROZEN_BINOMIAL_THRESHOLD) <= FROZEN_ROUTER_CUTOFF


def route_stage2(pb_llr: np.ndarray, pb_to_mamba_threshold: float) -> np.ndarray:
    """New PB-margin router: which PB-routed loci are uncertain enough for Mamba.

    Args:
        pb_llr: ``[N]`` Poisson-binomial log-likelihood ratios.
        pb_to_mamba_threshold: Margin cutoff. A locus whose PB score sits within
            this many LLR units of ``FROZEN_PB_THRESHOLD`` is handed to Mamba.

    Returns:
        ``[N]`` boolean mask, True meaning "route to Mamba" (among PB-routed loci;
        callers are responsible for intersecting with ``route_stage1``'s output
        where that matters).
    """
    pb_llr = np.asarray(pb_llr, dtype=np.float64)
    margin = np.abs(pb_llr - FROZEN_PB_THRESHOLD)
    return margin < float(pb_to_mamba_threshold)


def run_cascade(binomial_llr: np.ndarray, pb_llr: np.ndarray, mamba_score: np.ndarray,
                pb_to_mamba_threshold: float) -> tuple[np.ndarray, np.ndarray]:
    """Run the full three-stage cascade and attribute each call to its stage.

    Args:
        binomial_llr: ``[N]`` fixed-epsilon binomial LLR.
        pb_llr: ``[N]`` Poisson-binomial LLR.
        mamba_score: ``[N]`` PB-residual Mamba decision score (same axis as
            ``pb_llr``: ``score >= FROZEN_PB_THRESHOLD`` is the SNP call).
        pb_to_mamba_threshold: Frozen (or swept) stage-2 margin threshold.

    Returns:
        ``(final_calls, stage_used)``: ``final_calls`` is an ``[N]`` boolean SNP
        decision per locus; ``stage_used`` is an ``[N]`` array of strings in
        ``{"binomial", "pb", "mamba"}`` naming which stage produced that call.
    """
    binomial_llr = np.asarray(binomial_llr, dtype=np.float64)
    pb_llr = np.asarray(pb_llr, dtype=np.float64)
    mamba_score = np.asarray(mamba_score, dtype=np.float64)
    n = binomial_llr.size
    assert pb_llr.size == n and mamba_score.size == n, "input length mismatch"

    routed_to_pb = route_stage1(binomial_llr)
    routed_to_mamba = routed_to_pb & route_stage2(pb_llr, pb_to_mamba_threshold)
    pb_final = routed_to_pb & ~routed_to_mamba

    final_calls = np.where(
        routed_to_mamba, mamba_score >= FROZEN_PB_THRESHOLD,
        np.where(pb_final, pb_llr >= FROZEN_PB_THRESHOLD,
                 binomial_llr >= FROZEN_BINOMIAL_THRESHOLD))

    stage_used = np.full(n, STAGE_BINOMIAL, dtype=object)
    stage_used[pb_final] = STAGE_PB
    stage_used[routed_to_mamba] = STAGE_MAMBA

    # Every locus produced by exactly one stage.
    assert np.array_equal(~routed_to_pb, stage_used == STAGE_BINOMIAL)
    assert np.array_equal(pb_final, stage_used == STAGE_PB)
    assert np.array_equal(routed_to_mamba, stage_used == STAGE_MAMBA)
    counts = (stage_used == STAGE_BINOMIAL).sum() + (stage_used == STAGE_PB).sum() \
        + (stage_used == STAGE_MAMBA).sum()
    assert counts == n, "every locus must be produced by exactly one stage"

    return final_calls.astype(bool), stage_used


__all__ = [
    "FROZEN_BINOMIAL_THRESHOLD", "FROZEN_PB_THRESHOLD", "FROZEN_ROUTER_CUTOFF",
    "STAGE_BINOMIAL", "STAGE_PB", "STAGE_MAMBA",
    "route_stage1", "route_stage2", "run_cascade",
]
