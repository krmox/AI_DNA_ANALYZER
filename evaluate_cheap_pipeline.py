"""Two-stage caller: a cheap pre-router decides which loci deserve PB / Mamba.

The architecture under test
---------------------------
::

    BAM -> pileup counts            (paid by every arm, common cost)
             |
             +-- fixed-epsilon binomial LLR        1.33 M loci/s
             |
        cheap pre-router  (no PB, no Mamba, no labels)
             |
      +------+---------------------------+
      |                                  |
   confident                          uncertain
      |                                  |
   binomial call                     Poisson-binomial      8 k loci/s
   (llr >= T_b)                          |
                                    [optional] PB-margin -> Mamba
                                         |
                                    PB / Mamba call

Unlike ``evaluate_hybrid``, the expensive stage is genuinely skipped: a locus
the router calls confident never has its PB LLR computed. Accuracy is scored
using the cached PB values, but the compute model charges PB only where a
deployment would actually pay for it -- the routed loci, plus every locus of
any window the neural stage runs on, since the model reads the PB prior across
the whole window.

Arms compared
-------------
1. ``pb_only``                  PB everywhere
2. ``mamba_only``               PB everywhere + neural everywhere
3. ``hybrid_post_pb``           PB everywhere + neural on the frozen 3% (devlog 12)
4. ``cheap_router + PB``        PB only where the pre-router says so
5. ``cheap_router + PB + Mamba``  and neural inside that subset

Protocol
--------
Router trained on the **train** split against a distillation target computed
from PB on that split (no truth labels anywhere in router fitting). Coverage
selected on **validation** by a pre-registered rule. Frozen. Test evaluated
once.
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch

from analysis_30x import depth_bins_30x
from cheap_router import (
    alt_count_uncertainty,
    binomial_llr,
    binomial_margin_uncertainty,
    cheap_features,
    depth_uncertainty,
    fit_logistic_router,
    measure_throughput,
    oracle_target,
)
from config import LABEL_SNP
from evaluate_binomial_baseline import pick_threshold, prf
from evaluate_pb_residual import model_scores
from evaluate_quality_error import LLR_GRID, paired_bootstrap
from genomic_split import assign_window_roles, split_summary
from hybrid_router import (
    alt_quality_uncertainty,
    composite_uncertainty,
    cutoff_for_coverage,
    pb_margin_uncertainty,
    route,
    routing_enrichment,
)
from residual_metrics import bootstrap_f1
from train_raw_pileup import SEQ_LEN

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

#: Pre-registered coverage grid for the expensive path.
COVERAGE_GRID: tuple[float, ...] = (
    0.0, 0.001, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.0,
)

#: Pre-registered operating-point rule: the smallest coverage whose validation
#: F1 comes within this much of PB-everywhere. Chosen before any test scoring.
F1_TOLERANCE: float = 0.0005

#: Seeds, matching every previous experiment in this project.
SEEDS: tuple[int, ...] = (20260811, 424242, 13371337)


def two_stage_calls(routed: np.ndarray, cheap_score: np.ndarray, cheap_threshold: float,
                    pb_llr: np.ndarray, pb_threshold: float,
                    neural_score: np.ndarray | None = None,
                    neural_threshold: float | None = None,
                    neural_mask: np.ndarray | None = None) -> np.ndarray:
    """Decisions of the two-stage caller.

    Args:
        routed: ``[N]`` loci the pre-router sends to the expensive path.
        cheap_score: ``[N]`` binomial LLR, used on the confident majority.
        cheap_threshold: The binomial caller's frozen threshold.
        pb_llr: ``[N]`` Poisson-binomial LLR; consulted only where routed.
        pb_threshold: PB's frozen threshold.
        neural_score: ``[N]`` model score, optional third stage.
        neural_threshold: The checkpoint's frozen threshold.
        neural_mask: ``[N]`` loci handed to the neural stage; must be a subset
            of ``routed``.

    Returns:
        ``[N]`` boolean calls.
    """
    calls = np.where(routed, pb_llr >= pb_threshold, cheap_score >= cheap_threshold)
    if neural_score is not None:
        assert neural_mask is not None
        assert not np.any(neural_mask & ~routed), "neural stage outside the routed set"
        calls = np.where(neural_mask, neural_score >= neural_threshold, calls)
    return calls


@dataclass
class WindowIndex:
    """Maps each scored locus to the inference window that would contain it.

    Scoring is restricted to the SNP/Normal frame, which is *not* window
    aligned -- indel loci are dropped from the middle of windows. Reshaping the
    frame-restricted array into windows would therefore silently mix loci from
    different windows, so the window each locus belongs to is carried
    explicitly instead.

    Attributes:
        of_locus: ``[N]`` window id for each scored locus.
        total: Number of windows in the split, including any holding no scored
            locus at all -- the caller still pays for those.
    """

    of_locus: np.ndarray
    total: int


def build_window_index(split_mask: np.ndarray, frame_mask: np.ndarray,
                       seq_len: int = SEQ_LEN) -> WindowIndex:
    """Window ids for the scored loci of one split.

    Args:
        split_mask: ``[N_all]`` all loci of the split, window aligned.
        frame_mask: ``[N_all]`` the scored subset (split AND frame).
        seq_len: Window width.

    Returns:
        The window index for the scored loci.
    """
    global_window = np.arange(split_mask.size) // seq_len
    split_windows = np.unique(global_window[split_mask])
    return WindowIndex(
        of_locus=np.searchsorted(split_windows, global_window[frame_mask]),
        total=int(split_windows.size))


def compute_model(routed: np.ndarray, neural_mask: np.ndarray | None, rates: dict,
                  windows: WindowIndex, seq_len: int = SEQ_LEN) -> dict:
    """Estimated cost of one configuration, from measured per-stage rates.

    Args:
        routed: ``[N]`` mask of loci sent to PB.
        neural_mask: ``[N]`` mask of loci sent to the neural stage, or None.
        rates: Measured loci/second for each stage.
        windows: Window index for this split.
        seq_len: Window width, since the neural stage is billed per window.

    Returns:
        Per-stage seconds plus totals and coverage fractions.
    """
    total = routed.size
    if neural_mask is None or not neural_mask.any():
        touched = 0
        in_touched_window = np.zeros(total, dtype=bool)
    else:
        occupancy = np.bincount(windows.of_locus[neural_mask], minlength=windows.total)
        touched = int(np.count_nonzero(occupancy))
        in_touched_window = occupancy[windows.of_locus] > 0
    window_fraction = touched / max(windows.total, 1)

    # The neural stage reads the PB prior for *every* locus of a window it
    # runs on, not only for the loci that were routed. Charging PB for the
    # routed set alone would undercount the third stage's true cost, so the
    # PB bill is the whole of each touched window plus any routed locus that
    # falls outside them.
    pb_loci = touched * seq_len + int(np.count_nonzero(routed & ~in_touched_window))
    router_seconds = total / rates["router"]
    binomial_seconds = total / rates["binomial"]
    pb_seconds = pb_loci / rates["poisson_binomial"]
    neural_seconds = (touched * seq_len) / rates["mamba"]
    # The neural stage consumes the 44-channel evidence, which is only built
    # for the windows it actually runs on.
    evidence_seconds = (touched * seq_len) / rates["evidence"]
    return {
        "pb_locus_fraction": float(pb_loci / max(total, 1)),
        "routed_locus_fraction": float(routed.mean()),
        "neural_locus_fraction": 0.0 if neural_mask is None else float(neural_mask.mean()),
        "neural_window_fraction": window_fraction,
        "windows_touched": touched,
        "windows_total": windows.total,
        "window_amplification": (float(window_fraction / neural_mask.mean())
                                 if neural_mask is not None and neural_mask.any()
                                 else float("nan")),
        "router_seconds": router_seconds,
        "binomial_seconds": binomial_seconds,
        "poisson_binomial_seconds": pb_seconds,
        "evidence_seconds": evidence_seconds,
        "mamba_seconds": neural_seconds,
        "total_seconds": (router_seconds + binomial_seconds + pb_seconds
                          + evidence_seconds + neural_seconds),
        "loci_per_second": total / (router_seconds + binomial_seconds + pb_seconds
                                    + evidence_seconds + neural_seconds),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", default="cache/big15x_evidence.npz")
    parser.add_argument("--counts", default="cache/big15x_counts.npz")
    parser.add_argument("--positions", default="cache/big15x_positions.npy")
    parser.add_argument("--region-start", type=int, default=32000000)
    parser.add_argument("--checkpoint", action="append", default=[], metavar="NAME=PATH")
    parser.add_argument("--reads", default="data/giab_hg002_15x_readlevel/test_15x_reads.npz")
    parser.add_argument("--stat-loci", type=int, default=100_000)
    parser.add_argument("--out", default="results/cheap_router/cheap_results.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    blob = np.load(args.evidence, allow_pickle=True)
    features, labels, pb_llr, depth = (blob["features"], blob["labels"],
                                       blob["pb_llr"], blob["depth"])
    feature_names = list(blob["feature_names"])
    counts = np.load(args.counts)["counts"]
    positions = np.load(args.positions)

    masks = assign_window_roles(positions, args.region_start, seq_len=SEQ_LEN)
    snp = labels == LABEL_SNP
    frame = (labels == 0) | snp
    train = masks["train"] & frame
    validation = masks["validation"] & frame
    test = masks["test"] & frame

    report: dict = {
        "protocol": {
            "router_training": "distillation of the PB-margin oracle on the TRAIN split; "
                               "no truth label enters router fitting",
            "coverage_grid": list(COVERAGE_GRID),
            "f1_tolerance_rule": F1_TOLERANCE,
            "seeds": list(SEEDS),
        },
        "split": split_summary(positions, labels, args.region_start, seq_len=SEQ_LEN),
    }

    validation_windows = build_window_index(masks["validation"], validation)
    test_windows = build_window_index(masks["test"], test)

    # ---- cheap statistics, computed once for every locus -------------------
    logger.info("computing the cheap binomial LLR for all %d loci", labels.size)
    cheap_llr = binomial_llr(counts)

    # ---- thresholds, all selected on validation ----------------------------
    pb_threshold, _ = pick_threshold(pb_llr[validation], snp[validation], LLR_GRID)
    cheap_threshold, cheap_validation = pick_threshold(cheap_llr[validation],
                                                       snp[validation], LLR_GRID)
    logger.info("validation thresholds: PB=%.2f binomial=%.2f (binomial F1=%.4f)",
                pb_threshold, cheap_threshold, cheap_validation["f1"])

    neural: dict[str, np.ndarray] = {}
    neural_thresholds: dict[str, float] = {}
    for entry in args.checkpoint:
        name, _, path = entry.partition("=")
        score, diagnostics, _ = model_scores(Path(path), features, pb_llr, device)
        neural[name] = score
        neural_thresholds[name] = diagnostics["validation_threshold"]
    report["thresholds"] = {
        "poisson_binomial": pb_threshold, "binomial_cheap_path": cheap_threshold,
        **neural_thresholds}

    # ---- measured stage rates ---------------------------------------------
    logger.info("measuring stage throughput")
    sample = counts[test][:1_000_000]
    rates = {
        "binomial": measure_throughput(binomial_llr, sample)["loci_per_second"],
        "router": measure_throughput(
            lambda c: binomial_margin_uncertainty(c, cheap_threshold), sample
        )["loci_per_second"],
        "cheap_features": measure_throughput(cheap_features, sample)["loci_per_second"],
    }
    rates["poisson_binomial"], rates["evidence"] = measured_statistical_rates(
        args.reads, args.stat_loci)
    rates["mamba"] = measured_neural_rate(args.checkpoint, features[test], pb_llr[test],
                                          device)
    report["throughput_loci_per_second"] = rates
    logger.info("rates: %s", {k: round(v) for k, v in rates.items()})

    # ---- router family -----------------------------------------------------
    # Analytical rules are seed-independent; the learned routers are refitted
    # per seed so replication is meaningful.
    logger.info("fitting learned routers")
    target = oracle_target(pb_llr[train], pb_threshold)
    learned: dict[str, dict] = {}
    report["router_fits"] = {}
    for seed in SEEDS:
        for use_binomial in (False, True):
            name = f"learned_{'counts_binom' if use_binomial else 'counts'}_seed{seed}"
            router, diagnostics = fit_logistic_router(
                counts[train], target, use_binomial=use_binomial,
                binomial_threshold=cheap_threshold, seed=seed)
            learned[name] = {
                "validation": router.score(counts[validation], cheap_threshold,
                                           cheap_llr[validation]),
                "test": router.score(counts[test], cheap_threshold, cheap_llr[test]),
            }
            report["router_fits"][name] = diagnostics

    rng = np.random.default_rng(20260813)
    signals: dict[str, dict[str, np.ndarray]] = {
        # -- true pre-routers: computable before PB --------------------------
        "binomial_margin": {
            split_name: binomial_margin_uncertainty(counts[mask], cheap_threshold,
                                                    cheap_llr[mask])
            for split_name, mask in (("validation", validation), ("test", test))},
        "alt_count": {
            split_name: alt_count_uncertainty(counts[mask])
            for split_name, mask in (("validation", validation), ("test", test))},
        "depth_rule": {
            split_name: depth_uncertainty(counts[mask])
            for split_name, mask in (("validation", validation), ("test", test))},
        # -- controls --------------------------------------------------------
        "CONTROL_random": {
            split_name: rng.uniform(size=int(mask.sum()))
            for split_name, mask in (("validation", validation), ("test", test))},
        "CONTROL_reverse_binomial_margin": {
            split_name: -binomial_margin_uncertainty(counts[mask], cheap_threshold,
                                                     cheap_llr[mask])
            for split_name, mask in (("validation", validation), ("test", test))},
        # -- post-PB diagnostics: NOT pre-routers ----------------------------
        "ORACLE_pb_margin": {
            split_name: pb_margin_uncertainty(pb_llr[mask], pb_threshold)
            for split_name, mask in (("validation", validation), ("test", test))},
    }
    for name, scores in learned.items():
        signals[name] = scores

    report["signal_roles"] = {
        "true_pre_router": ["binomial_margin", "alt_count", "depth_rule"]
        + list(learned.keys()),
        "control": ["CONTROL_random", "CONTROL_reverse_binomial_margin"],
        "post_pb_diagnostic_only": ["ORACLE_pb_margin"],
    }

    # ---- curves ------------------------------------------------------------
    report["validation_curves"] = {}
    report["test_curves"] = {}
    for name, scores in signals.items():
        report["validation_curves"][name] = coverage_curve(
            scores["validation"], scores["validation"], cheap_llr[validation],
            cheap_threshold, pb_llr[validation], pb_threshold, snp[validation],
            {k: v[validation] for k, v in neural.items()}, neural_thresholds, rates,
            validation_windows)
        report["test_curves"][name] = coverage_curve(
            scores["test"], scores["validation"], cheap_llr[test], cheap_threshold,
            pb_llr[test], pb_threshold, snp[test],
            {k: v[test] for k, v in neural.items()}, neural_thresholds, rates,
            test_windows)
        logger.info("curve done: %s", name)

    # ---- frozen operating point -------------------------------------------
    pb_validation_f1 = prf(pb_llr[validation] >= pb_threshold, snp[validation])["f1"]
    curve = report["validation_curves"]["binomial_margin"]
    eligible = [float(k) for k, entry in curve.items()
                if entry["cheap_plus_pb"]["f1"] >= pb_validation_f1 - F1_TOLERANCE]
    frozen_coverage = min(eligible) if eligible else 1.0
    report["frozen_configuration"] = {
        "signal": "binomial_margin",
        "rule": "smallest coverage whose validation F1 is within F1_TOLERANCE of "
                "PB-everywhere; signal fixed a priori as the cheapest analytical rule "
                "that ranks loci",
        "coverage": frozen_coverage,
        "cutoff": cutoff_for_coverage(signals["binomial_margin"]["validation"],
                                      frozen_coverage),
        "pb_validation_f1": pb_validation_f1,
        "eligible_coverages": sorted(eligible),
    }
    logger.info("FROZEN: binomial_margin at coverage %.4f", frozen_coverage)

    # ---- reference arms ----------------------------------------------------
    report["arms"] = reference_arms(pb_llr, pb_threshold, cheap_llr, cheap_threshold,
                                    neural, neural_thresholds, features, feature_names,
                                    depth, snp, validation, test, rates, signals,
                                    frozen_coverage, test_windows)

    # ---- depth stratification at the frozen point --------------------------
    cutoff = report["frozen_configuration"]["cutoff"]
    routed_test = route(signals["binomial_margin"]["test"], cutoff)
    report["by_depth"] = {}
    for label, bin_mask in depth_bins_30x(depth):
        stratum = bin_mask[test]
        if not stratum.any():
            continue
        entry = {
            "loci": int(stratum.sum()), "snp": int(snp[test][stratum].sum()),
            "routed_fraction": float(routed_test[stratum].mean()),
            "pb_only": prf(pb_llr[test][stratum] >= pb_threshold, snp[test][stratum]),
            "binomial_only": prf(cheap_llr[test][stratum] >= cheap_threshold,
                                 snp[test][stratum]),
            "cheap_plus_pb": prf(
                two_stage_calls(routed_test, cheap_llr[test], cheap_threshold,
                                pb_llr[test], pb_threshold)[stratum], snp[test][stratum]),
        }
        for name, score in neural.items():
            entry[f"mamba_only_{name}"] = prf(
                score[test][stratum] >= neural_thresholds[name], snp[test][stratum])
        report["by_depth"][label] = entry

    report["dataset"] = {
        "loci": int(labels.size), "train_frame_loci": int(train.sum()),
        "validation_frame_loci": int(validation.sum()), "test_frame_loci": int(test.sum()),
        "test_snp": int(snp[test].sum()), "validation_snp": int(snp[validation].sum())}

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2, default=float))
    logger.info("Wrote %s", args.out)
    print_summary(report, frozen_coverage)


def measured_statistical_rates(reads_npz: str, loci: int) -> tuple[float, float]:
    """Measured PB-prior and 44-channel-evidence rates, on a real read tensor.

    Args:
        reads_npz: Cached read-level tensor for a 15x region.
        loci: Loci to time.

    Returns:
        ``(pb_loci_per_second, evidence_loci_per_second)``.
    """
    from benchmark_hybrid_runtime import time_statistical_stages
    measured = time_statistical_stages(reads_npz, loci)
    return (measured["poisson_binomial_prior"]["loci_per_second"],
            measured["feature_extraction"]["loci_per_second"])


@torch.no_grad()
def measured_neural_rate(checkpoints: list[str], features: np.ndarray,
                         pb_llr: np.ndarray, device: torch.device) -> float:
    """Measured Mamba throughput on this machine, in loci/second."""
    if not checkpoints:
        return float("inf")
    from benchmark_hybrid_runtime import time_neural
    from locus_evidence import FEATURE_DIM
    from model_pb_residual import PoissonBinomialResidualMamba

    _, _, path = checkpoints[0].partition("=")
    state = torch.load(path, map_location=device, weights_only=False)
    model = PoissonBinomialResidualMamba(
        feature_dim=FEATURE_DIM, use_gate=bool(state.get("use_gate", False))).to(device)
    model.load_state_dict(state["model_state_dict"])
    model.eval()
    subset = min(features.shape[0] // SEQ_LEN * SEQ_LEN, 1_000_000 // SEQ_LEN * SEQ_LEN)
    routed = np.ones(subset, dtype=bool)
    timing = time_neural(model, features[:subset], pb_llr[:subset], routed, device)
    return float(timing["loci_per_second"])


def coverage_curve(signal: np.ndarray, fitting_signal: np.ndarray, cheap_llr: np.ndarray,
                   cheap_threshold: float, pb_llr: np.ndarray, pb_threshold: float,
                   truth: np.ndarray, neural: dict[str, np.ndarray],
                   neural_thresholds: dict[str, float], rates: dict,
                   windows: WindowIndex) -> dict:
    """Accuracy and compute across the coverage grid for one routing signal.

    Cutoffs always come from ``fitting_signal`` (the validation distribution),
    so the test curve uses frozen cutoffs rather than test quantiles.

    Args:
        signal: Uncertainty on the split being scored.
        fitting_signal: Uncertainty on validation, supplying the cutoffs.
        cheap_llr: Binomial LLR on the scored split.
        cheap_threshold: Binomial threshold.
        pb_llr: PB LLR on the scored split.
        pb_threshold: PB threshold.
        truth: Boolean truth on the scored split.
        neural: Arm name -> model score on the scored split.
        neural_thresholds: Arm name -> threshold.
        rates: Measured stage throughputs.
        windows: Window index for the split being scored.

    Returns:
        Coverage key -> metrics.
    """
    pb_calls = pb_llr >= pb_threshold
    pb_stats = prf(pb_calls, truth)
    out: dict[str, dict] = {}
    for coverage in COVERAGE_GRID:
        cutoff = cutoff_for_coverage(fitting_signal, coverage)
        routed = route(signal, cutoff)
        calls = two_stage_calls(routed, cheap_llr, cheap_threshold, pb_llr, pb_threshold)
        stats = prf(calls, truth)
        entry = {
            "requested_coverage": coverage,
            "achieved_coverage": float(routed.mean()),
            "cheap_plus_pb": stats,
            "delta_f1_vs_pb": stats["f1"] - pb_stats["f1"],
            "routing": routing_enrichment(routed, pb_calls, truth),
            "compute": compute_model(routed, None, rates, windows),
            "with_mamba": {},
        }
        # Third stage: the neural model runs on exactly the routed set. The
        # pre-router has already selected the uncertain loci, so adding a second
        # constant here would be an unfrozen magic number for no stated reason.
        neural_mask = routed
        for name, score in neural.items():
            neural_calls = two_stage_calls(routed, cheap_llr, cheap_threshold, pb_llr,
                                           pb_threshold, score, neural_thresholds[name],
                                           neural_mask)
            neural_stats = prf(neural_calls, truth)
            entry["with_mamba"][name] = {**neural_stats,
                                         "delta_f1_vs_pb": neural_stats["f1"] - pb_stats["f1"]}
        entry["with_mamba_mean_f1"] = (
            float(np.mean([s["f1"] for s in entry["with_mamba"].values()]))
            if entry["with_mamba"] else float("nan"))
        entry["compute_with_mamba"] = compute_model(routed, neural_mask, rates, windows)
        out[f"{coverage:g}"] = entry
    return out


def reference_arms(pb_llr, pb_threshold, cheap_llr, cheap_threshold, neural,
                   neural_thresholds, features, feature_names, depth, snp,
                   validation, test, rates, signals, frozen_coverage,
                   test_windows) -> dict:
    """The five architectures, scored on test with everything frozen."""
    total = int(test.sum())
    arms: dict[str, dict] = {}

    pb_calls = pb_llr[test] >= pb_threshold
    everything = np.ones(total, dtype=bool)
    arms["pb_only"] = {
        **prf(pb_calls, snp[test]), "f1_ci95": bootstrap_f1(pb_calls, snp[test]),
        "compute": compute_model(everything, None, rates, test_windows)}
    arms["binomial_only"] = {
        **prf(cheap_llr[test] >= cheap_threshold, snp[test]),
        "compute": compute_model(np.zeros(total, dtype=bool), None, rates, test_windows)}

    for name, score in neural.items():
        calls = score[test] >= neural_thresholds[name]
        arms[f"mamba_only_{name}"] = {
            **prf(calls, snp[test]), "f1_ci95": bootstrap_f1(calls, snp[test]),
            "vs_pb": paired_bootstrap(calls, pb_calls, snp[test]),
            "compute": compute_model(everything, everything, rates, test_windows)}

    # Arm 3: the previous hybrid, PB everywhere plus the frozen 3% neural subset.
    reference = {"margin": pb_margin_uncertainty(pb_llr[validation], pb_threshold),
                 "quality": alt_quality_uncertainty(features[validation], feature_names)}
    previous = composite_uncertainty(pb_llr[test], pb_threshold, features[test],
                                     feature_names, reference)
    previous_mask = route(previous, 3.628494794922591)
    for name, score in neural.items():
        calls = np.where(previous_mask, score[test] >= neural_thresholds[name], pb_calls)
        arms[f"hybrid_post_pb_{name}"] = {
            **prf(calls, snp[test]), "vs_pb": paired_bootstrap(calls, pb_calls, snp[test]),
            "compute": compute_model(everything, previous_mask, rates, test_windows)}

    # Arms 4 and 5: the frozen cheap pre-router.
    routed = route(signals["binomial_margin"]["test"],
                   cutoff_for_coverage(signals["binomial_margin"]["validation"],
                                       frozen_coverage))
    calls = two_stage_calls(routed, cheap_llr[test], cheap_threshold, pb_llr[test],
                            pb_threshold)
    arms["cheap_router_pb"] = {
        **prf(calls, snp[test]), "f1_ci95": bootstrap_f1(calls, snp[test]),
        "disagreements_vs_pb_only": int(np.count_nonzero(calls != pb_calls)),
        "vs_pb": paired_bootstrap(calls, pb_calls, snp[test]),
        "routing": routing_enrichment(routed, pb_calls, snp[test]),
        "compute": compute_model(routed, None, rates, test_windows)}

    neural_mask = routed
    for name, score in neural.items():
        combined = two_stage_calls(routed, cheap_llr[test], cheap_threshold, pb_llr[test],
                                   pb_threshold, score[test], neural_thresholds[name],
                                   neural_mask)
        arms[f"cheap_router_pb_mamba_{name}"] = {
            **prf(combined, snp[test]),
            "f1_ci95": bootstrap_f1(combined, snp[test]),
            "vs_pb": paired_bootstrap(combined, pb_calls, snp[test]),
            "vs_mamba_only": paired_bootstrap(
                combined, score[test] >= neural_thresholds[name], snp[test]),
            "disagreements_vs_mamba_only": int(np.count_nonzero(
                combined != (score[test] >= neural_thresholds[name]))),
            "compute": compute_model(routed, neural_mask, rates, test_windows)}
    return arms


def print_summary(report: dict, frozen_coverage: float) -> None:
    """Console report."""
    print("\n=== MEASURED STAGE THROUGHPUT (loci/s) ===")
    for name, rate in report["throughput_loci_per_second"].items():
        print(f"  {name:<20} {rate:>12,.0f}")

    print("\n=== ROUTING SIGNALS: PB errors captured on TEST ===")
    print(f"  {'signal':<34} " + " ".join(f"{c:>7g}" for c in (0.001, 0.005, 0.01, 0.05)))
    for name, curve in report["test_curves"].items():
        row = " ".join(f"{curve[f'{c:g}']['routing']['error_captured_fraction']:>7.3f}"
                       for c in (0.001, 0.005, 0.01, 0.05))
        print(f"  {name:<34} {row}")

    print("\n=== CHEAP PRE-ROUTER + PB, test accuracy/compute curve "
          "(signal=binomial_margin) ===")
    print(f"  {'cov':>8} {'PB loci%':>9} {'F1':>8} {'dF1 vs PB':>10} {'errCaught':>10} "
          f"{'enrich':>8} {'e2e s':>8} {'loci/s':>10}")
    for key, entry in report["test_curves"]["binomial_margin"].items():
        compute = entry["compute"]
        print(f"  {key:>8} {entry['achieved_coverage']:>9.4f} "
              f"{entry['cheap_plus_pb']['f1']:>8.4f} {entry['delta_f1_vs_pb']:>+10.4f} "
              f"{entry['routing']['error_captured_fraction']:>10.3f} "
              f"{entry['routing']['enrichment']:>8.1f} "
              f"{compute['total_seconds']:>8.1f} {compute['loci_per_second']:>10,.0f}")

    print(f"\n=== ARCHITECTURES (test, frozen; pre-router at coverage {frozen_coverage:g}) ===")
    print(f"  {'arm':<32} {'F1':>8} {'TP':>6} {'FP':>5} {'FN':>5} "
          f"{'PB%':>7} {'win%':>7} {'e2e s':>8} {'speedup':>8}")
    baseline = report["arms"]["pb_only"]["compute"]["total_seconds"]
    for name, arm in report["arms"].items():
        compute = arm["compute"]
        print(f"  {name:<32} {arm['f1']:>8.4f} {arm['tp']:>6} {arm['fp']:>5} "
              f"{arm['fn']:>5} {compute['pb_locus_fraction']:>7.4f} "
              f"{compute['neural_window_fraction']:>7.4f} "
              f"{compute['total_seconds']:>8.1f} "
              f"{baseline / compute['total_seconds']:>7.2f}x")

    print("\n=== DEPTH STRATA at the frozen point (F1) ===")
    print(f"  {'depth':<8} {'loci':>9} {'SNP':>6} {'routed':>8} {'PB':>8} "
          f"{'binomial':>9} {'cheap+PB':>9}")
    for label, entry in report["by_depth"].items():
        print(f"  {label:<8} {entry['loci']:>9,} {entry['snp']:>6} "
              f"{entry['routed_fraction']:>8.4f} {entry['pb_only']['f1']:>8.4f} "
              f"{entry['binomial_only']['f1']:>9.4f} {entry['cheap_plus_pb']['f1']:>9.4f}")


if __name__ == "__main__":
    main()
