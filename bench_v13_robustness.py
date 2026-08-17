"""bench_v13 — frozen cheap-router → PB after the PB depth-clipping fix.

Protocol is pre-registered in ``History/13_DEVLOG.md`` §3–§6 and is not
re-derived here. Nothing in this file fits, selects or modifies a threshold:
the three frozen constants are imported from ``cascade`` and only ever read,
and the cutoff sweep exists solely to draw the accuracy-vs-compute curve.

Two data families are scored, both with the **repaired** Poisson-binomial:

* **Group A** — six freshly extracted region/depth realisations
  (``cache/bench_v13/*.npz``), each carrying the legacy (pre-fix) PB LLR
  alongside the repaired one, computed on identical reads in the same pass, so
  the fix can be A/B compared exactly.
* **Group B** — the seven already-cached regions of ``results/robustness_v11``
  (19.44 M loci). Their cached PB is reused unchanged; ``pb_fix_scope`` records
  per region whether any locus could have been touched by the defect, which is
  what licenses the reuse.

Every metric helper is imported from the existing benchmark modules rather than
reimplemented, so the numbers are on the same definitions as devlogs 11 and 12.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import numpy as np

from bench_v12_stage1 import arms_on_mask, strata_masks
from binomial_baseline import BinomialVariantCaller
from cascade import (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD,
                     FROZEN_ROUTER_CUTOFF)
from config import LABEL_SNP
from evaluate_binomial_baseline import prf
from evaluate_quality_error import paired_bootstrap
from residual_metrics import vaf_from_counts
from robustness_benchmark import (accuracy_block, controls, cutoff_sweep,
                                  depth_stratified, evaluate_arms,
                                  load_all_regions, quality_stratified, route_mask)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

#: Group A, exactly as pre-registered. ``tainted`` marks the regions that fall
#: inside the router's cutoff-selection span (chr21:32-40 Mb) and therefore
#: cannot support a *generalisation* claim about the router. They are still
#: scored, and are fully valid for the PB-defect question.
NEW_REGIONS: dict[str, dict] = {
    "new32_33M_15x": {"path": "cache/bench_v13/new32_33M_15x.npz", "tainted": True},
    "new32_33M_30x": {"path": "cache/bench_v13/new32_33M_30x.npz", "tainted": True},
    "new32_33M_full": {"path": "cache/bench_v13/new32_33M_full.npz", "tainted": True},
    "r30M_full": {"path": "cache/bench_v13/r30M_full.npz", "tainted": False},
    "re31_32M_full": {"path": "cache/bench_v13/re31_32M_full.npz", "tainted": False},
    "re31_32M_30x": {"path": "cache/bench_v13/re31_32M_30x.npz", "tainted": False},
}

#: Finer than ``robustness_benchmark.DEPTH_BINS``, whose top bin is "30x+" and
#: would hide the 48x boundary this benchmark exists to probe. The existing
#: bins are still reported unchanged via ``depth_stratified``.
FINE_DEPTH_BINS = [(1, 9), (10, 29), (30, 47), (48, 59), (60, 79), (80, 10_000)]

VAF_BINS = [(0.0, 0.15), (0.15, 0.25), (0.25, 0.35), (0.35, 0.45),
            (0.45, 0.55), (0.55, 0.75), (0.75, 1.01)]

OUT_DIR = Path("results/bench_v13")
RESULTS = OUT_DIR / "robustness_results.json"
DISAGREEMENTS = OUT_DIR / "disagreements.json"
PB_FIX_AB = OUT_DIR / "pb_fix_ab.json"


# -- loading ----------------------------------------------------------------


def load_new(name: str, path: str) -> dict:
    """Load one Group-A extraction into the region dict the helpers expect."""
    blob = np.load(path, allow_pickle=True)
    counts = blob["counts"]
    depth = blob["depth"].astype(np.float64)
    region = {
        "name": name,
        "counts": counts,
        "binomial_llr": blob["binomial_llr"].astype(np.float64),
        "pb_llr": blob["pb_llr"].astype(np.float64),
        "labels": blob["labels"],
        "depth": depth,
        "mean_base_quality": np.where(depth > 0, counts[:, 7] / np.maximum(depth, 1), 0.0),
        "mean_mapq": np.where(depth > 0, counts[:, 8] / np.maximum(depth, 1), 0.0),
        "region": [int(x) for x in blob["region"]],
        "counts_time_s": float(blob["timing"][0]),
        "reads_time_s": float(blob["timing"][1]),
        "pb_time_s": float(blob["timing"][2]),
        "bam": str(blob["bam"][0]),
    }
    for optional in ("pb_llr_legacy", "k_full", "k_retained", "n_counted"):
        if optional in blob:
            region[optional] = blob[optional]
    return region


# -- stratifications not already provided by robustness_benchmark ------------


def _arms(region: dict, frame: np.ndarray) -> tuple[np.ndarray, ...]:
    """(snp, binomial calls, pb calls, router calls, routed mask) on a frame."""
    labels = region["labels"]
    snp = labels[frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    return snp, binomial_calls, pb_calls, np.where(routed, pb_calls, binomial_calls), routed


def scoring_frame(region: dict) -> np.ndarray:
    labels = region["labels"]
    return (labels == 0) | (labels == LABEL_SNP)


def _binned(region: dict, values: np.ndarray, bins, label) -> dict:
    """Metrics on inclusive ``[lo, hi]`` bins of a per-locus evidence value."""
    frame = scoring_frame(region)
    snp, binomial_calls, pb_calls, router_calls, routed = _arms(region, frame)
    values = values[frame]
    out = {}
    for lo, hi in bins:
        mask = (values >= lo) & (values <= hi)
        n = int(mask.sum())
        key = label(lo, hi)
        if n == 0:
            out[key] = {"loci": 0}
            continue
        out[key] = {
            "loci": n,
            "snp": int(snp[mask].sum()),
            "binomial_only": accuracy_block(binomial_calls[mask], snp[mask]),
            "pb_only": accuracy_block(pb_calls[mask], snp[mask]),
            "cheap_router_pb": accuracy_block(router_calls[mask], snp[mask]),
            "fraction_routed_to_pb": float(routed[mask].mean()),
            "disagreement_vs_pb": int(np.count_nonzero(router_calls[mask] != pb_calls[mask])),
            "disagreement_rate": float(np.mean(router_calls[mask] != pb_calls[mask])),
        }
        out[key]["delta_f1_vs_pb"] = (out[key]["cheap_router_pb"]["f1"]
                                      - out[key]["pb_only"]["f1"])
    return out


def fine_depth_stratified(region: dict) -> dict:
    """Depth strata that resolve the 48x tensor-width boundary."""
    return _binned(region, region["depth"], FINE_DEPTH_BINS, lambda lo, hi: f"{lo}-{hi}x")


def vaf_stratified(region: dict) -> dict:
    """VAF strata over loci that have at least one non-reference read."""
    vaf = vaf_from_counts(region["counts"])
    frame = scoring_frame(region)
    snp, binomial_calls, pb_calls, router_calls, routed = _arms(region, frame)
    vaf_frame = vaf[frame]
    out = {}
    for lo, hi in VAF_BINS:
        mask = (vaf_frame > lo) & (vaf_frame <= hi)
        n = int(mask.sum())
        key = f"vaf{lo:.2f}-{hi:.2f}"
        if n == 0:
            out[key] = {"loci": 0}
            continue
        out[key] = {
            "loci": n, "snp": int(snp[mask].sum()),
            "binomial_only": accuracy_block(binomial_calls[mask], snp[mask]),
            "pb_only": accuracy_block(pb_calls[mask], snp[mask]),
            "cheap_router_pb": accuracy_block(router_calls[mask], snp[mask]),
            "fraction_routed_to_pb": float(routed[mask].mean()),
            "disagreement_vs_pb": int(np.count_nonzero(router_calls[mask] != pb_calls[mask])),
        }
        out[key]["delta_f1_vs_pb"] = (out[key]["cheap_router_pb"]["f1"]
                                      - out[key]["pb_only"]["f1"])
    return out


# -- the PB fix, measured -----------------------------------------------------


def pb_fix_scope(region: dict) -> dict:
    """How much of this region the depth-clipping defect could have touched.

    A locus is exposed only if the ALT count taken from the uncapped pileup
    exceeds the number of reads the tensor retained. Where ``k_full`` /
    ``n_counted`` were not cached (Group B), depth vs tensor width is the
    conservative bound: at or below the tensor width no read is dropped, so no
    locus can be exposed.
    """
    depth = region["depth"]
    scope = {
        "max_depth": float(depth.max()),
        "loci_depth_gt_48": int(np.count_nonzero(depth > 48)),
        "fraction_depth_gt_48": float(np.mean(depth > 48)),
    }
    if "k_full" in region and "n_counted" in region:
        exposed = region["k_full"].astype(np.int64) > region["n_counted"].astype(np.int64)
        truncated = region["k_full"].astype(np.int64) != region["k_retained"].astype(np.int64)
        scope["loci_k_out_of_support"] = int(exposed.sum())
        scope["fraction_k_out_of_support"] = float(exposed.mean())
        scope["loci_k_rescaled"] = int(truncated.sum())
        scope["fraction_k_rescaled"] = float(truncated.mean())
        scope["exposed_snp"] = int(np.count_nonzero(exposed & (region["labels"] == LABEL_SNP)))
    else:
        scope["note"] = ("k/n_counted not cached; depth <= tensor width is the "
                         "sufficient condition for 'unaffected'")
    return scope


def pb_fix_ab(region: dict) -> dict | None:
    """PB-only and cascade metrics under the legacy and the repaired PB.

    Same reads, same counts, same thresholds: the only thing that moves is the
    defect. Everything here is about PB versus truth, not about the router.
    """
    if "pb_llr_legacy" not in region:
        return None
    frame = scoring_frame(region)
    labels = region["labels"]
    snp = labels[frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    depth = region["depth"][frame]

    out = {"loci": int(frame.sum()), "snp": int(snp.sum())}
    for arm, llr in (("legacy", region["pb_llr_legacy"][frame]),
                     ("fixed", region["pb_llr"][frame])):
        pb_calls = llr >= FROZEN_PB_THRESHOLD
        router_calls = np.where(routed, pb_calls, binomial_calls)
        zero = llr == 0.0
        out[arm] = {
            "pb_only": accuracy_block(pb_calls, snp),
            "cheap_router_pb": accuracy_block(router_calls, snp),
            "pb_llr_zero_loci": int(zero.sum()),
            "pb_llr_zero_snp": int(np.count_nonzero(zero & snp)),
            "pb_llr_zero_rate_depth_gt_48": (
                float(np.mean(zero[depth > 48])) if np.any(depth > 48) else None),
            "pb_llr_zero_snp_rate_depth_gt_48": (
                float(np.mean(zero[(depth > 48) & snp])) if np.any((depth > 48) & snp) else None),
        }
    legacy_pb = region["pb_llr_legacy"][frame] >= FROZEN_PB_THRESHOLD
    fixed_pb = region["pb_llr"][frame] >= FROZEN_PB_THRESHOLD
    out["pb_call_changes"] = int(np.count_nonzero(legacy_pb != fixed_pb))
    out["pb_llr_changes"] = int(np.count_nonzero(
        region["pb_llr_legacy"][frame] != region["pb_llr"][frame]))
    out["delta_f1_fixed_minus_legacy"] = (out["fixed"]["pb_only"]["f1"]
                                          - out["legacy"]["pb_only"]["f1"])
    out["fixed_vs_legacy_paired_bootstrap"] = paired_bootstrap(fixed_pb, legacy_pb, snp)
    return out


# -- failure boundary ---------------------------------------------------------


def disagreement_records(region: dict, limit: int = 2000) -> dict:
    """Every locus where ``router→PB`` differs from ``PB-only``, with its features.

    ``limit`` caps only the per-locus dump; the counts and the stratified
    summaries below it are computed over all disagreements.
    """
    frame = scoring_frame(region)
    labels = region["labels"]
    snp = labels[frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]
    depth = region["depth"][frame]
    quality = region["mean_base_quality"][frame]
    vaf = vaf_from_counts(region["counts"])[frame]
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    router_calls = np.where(routed, pb_calls, binomial_calls)

    disagree = np.nonzero(router_calls != pb_calls)[0]
    records = []
    for index in disagree[:limit]:
        records.append({
            "frame_index": int(index),
            "depth": float(depth[index]),
            "vaf": float(vaf[index]),
            "mean_base_quality": float(quality[index]),
            "binomial_llr": float(binomial_llr[index]),
            "pb_llr": float(pb_llr[index]),
            "router_margin": float(abs(binomial_llr[index] - FROZEN_BINOMIAL_THRESHOLD)),
            "pb_margin": float(abs(pb_llr[index] - FROZEN_PB_THRESHOLD)),
            "routed_to_pb": bool(routed[index]),
            "truth_is_snp": bool(snp[index]),
            "router_call": bool(router_calls[index]),
            "pb_call": bool(pb_calls[index]),
            # Direction: does the disagreement cost the cascade a TP, or spare
            # it an FP? Both are disagreements; only the first is a regression.
            "effect": ("router_fn" if snp[index] and pb_calls[index] and not router_calls[index]
                       else "router_fp" if not snp[index] and router_calls[index] and not pb_calls[index]
                       else "router_recovers_fn" if snp[index] and router_calls[index] and not pb_calls[index]
                       else "router_avoids_fp" if not snp[index] and pb_calls[index] and not router_calls[index]
                       else "other"),
        })

    def summarise(mask: np.ndarray) -> dict:
        return {"loci": int(mask.sum()),
                "disagreements": int(np.count_nonzero((router_calls != pb_calls) & mask)),
                "rate": float(np.mean((router_calls != pb_calls)[mask])) if mask.any() else None}

    return {
        "region": region["name"],
        "loci_scored": int(frame.sum()),
        "disagreements": int(disagree.size),
        "disagreement_rate": float(disagree.size / max(int(frame.sum()), 1)),
        "records": records,
        "records_truncated": bool(disagree.size > limit),
        "by_depth": {f"{lo}-{hi}x": summarise((depth >= lo) & (depth <= hi))
                     for lo, hi in FINE_DEPTH_BINS},
        "by_quality": {f"q{lo}-{hi}": summarise((depth > 0) & (quality >= lo) & (quality < hi))
                       for lo, hi in ((0, 20), (20, 25), (25, 30), (30, 35), (35, 60))},
        "by_vaf": {f"vaf{lo:.2f}-{hi:.2f}": summarise((vaf > lo) & (vaf <= hi))
                   for lo, hi in VAF_BINS},
        "by_truth": {"snp": summarise(snp), "non_snp": summarise(~snp)},
    }


def pooled_failure_boundary(disagreements: dict[str, dict]) -> dict:
    """Disagreement rate per difficulty bin, pooled over every scored region.

    Per-region counts are tiny by construction, so the boundary question --
    "does the disagreement rate rise in harder regimes?" -- can only be asked
    of the pool. Rates are per-locus within the bin, so a bin holding few loci
    can show a high rate on one event; ``disagreements`` is reported beside
    every rate so that case stays visible instead of hiding behind a ratio.
    """
    axes = ("by_depth", "by_quality", "by_vaf", "by_truth")
    pooled: dict[str, dict[str, dict[str, int]]] = {axis: {} for axis in axes}
    for entry in disagreements.values():
        for axis in axes:
            for key, block in entry[axis].items():
                cell = pooled[axis].setdefault(key, {"loci": 0, "disagreements": 0})
                cell["loci"] += int(block["loci"])
                cell["disagreements"] += int(block["disagreements"])
    for axis in axes:
        for cell in pooled[axis].values():
            cell["rate"] = cell["disagreements"] / cell["loci"] if cell["loci"] else None
    pooled["total"] = {
        "loci": sum(entry["loci_scored"] for entry in disagreements.values()),
        "disagreements": sum(entry["disagreements"] for entry in disagreements.values()),
    }
    pooled["total"]["rate"] = (pooled["total"]["disagreements"]
                               / max(pooled["total"]["loci"], 1))
    return pooled


# -- pooling ------------------------------------------------------------------


def pooled_block(snp: np.ndarray, binomial_llr: np.ndarray, pb_llr: np.ndarray) -> dict:
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    router_calls = np.where(routed, pb_calls, binomial_calls)
    entry = {
        "loci": int(snp.size), "snp": int(snp.sum()),
        "binomial_only": accuracy_block(binomial_calls, snp),
        "pb_only": accuracy_block(pb_calls, snp),
        "cheap_router_pb": accuracy_block(router_calls, snp),
        "fraction_routed_to_pb": float(routed.mean()),
        "pb_compute_loci": int(routed.sum()),
        "disagreement_vs_pb": int(np.count_nonzero(router_calls != pb_calls)),
        "disagreement_rate": float(np.mean(router_calls != pb_calls)),
        "vs_pb_paired_bootstrap": paired_bootstrap(router_calls, pb_calls, snp),
    }
    entry["delta_f1_vs_pb"] = (prf(router_calls, snp)["f1"] - prf(pb_calls, snp)["f1"])
    return entry


# -- wall clock ---------------------------------------------------------------


def measure_throughput(bam: str, vcf: str, bed: str, span: tuple[int, int]) -> dict:
    """Binomial and repaired-PB throughput, both measured in this process.

    Run after all extractions have finished so no other job is competing for
    the CPU; the caller records whether that held.
    """
    from extract_bench_v12 import poisson_binomial_for_chunk
    from pileup_counts import load_counts
    from read_level_pileup import load_reads

    fasta = "data/reference/chr21_full.fa"
    counts, _ = load_counts(fasta, bam, vcf, bed, span)
    reads, _ = load_reads(fasta, bam, vcf, bed, span)

    caller = BinomialVariantCaller(error_rate=0.01)
    start = time.perf_counter()
    caller.score_counts(counts[:, 0:4], counts[:, 9].astype(int))
    binomial_seconds = time.perf_counter() - start

    start = time.perf_counter()
    poisson_binomial_for_chunk(reads, counts)
    pb_seconds = time.perf_counter() - start

    n = counts.shape[0]
    return {
        "bam": bam, "span": list(span), "sample_loci": int(n),
        "mean_depth": float(counts[:, 6].mean()),
        "binomial_seconds": binomial_seconds,
        "binomial_loci_per_second": n / max(binomial_seconds, 1e-9),
        "pb_seconds": pb_seconds,
        "pb_loci_per_second": n / max(pb_seconds, 1e-9),
    }


def projected_wallclock(region: dict, throughput: dict) -> dict:
    """Cascade vs PB-only runtime from measured rates.

    Locus-level throughout: both callers here are per-locus, so no window
    amplification applies. The extraction cost (counts + reads) is reported
    beside the projection rather than folded into it, because both arms pay the
    counts pass and only the routed loci need the read pass.
    """
    frame = scoring_frame(region)
    routed = route_mask(region["binomial_llr"][frame], FROZEN_ROUTER_CUTOFF)
    n, k = int(frame.sum()), int(routed.sum())
    b_rate = throughput["binomial_loci_per_second"]
    pb_rate = throughput["pb_loci_per_second"]
    pb_only = n / pb_rate
    cascade = n / b_rate + k / pb_rate
    return {
        "loci": n, "pb_loci": k,
        "pb_compute_fraction": k / max(n, 1),
        "pb_only_seconds": pb_only,
        "cascade_seconds": cascade,
        "speedup": pb_only / max(cascade, 1e-12),
        "measured_extraction_counts_seconds": region.get("counts_time_s"),
        "measured_extraction_reads_seconds": region.get("reads_time_s"),
        "measured_extraction_pb_seconds": region.get("pb_time_s"),
        "note": ("extraction timings were produced with up to four extractions "
                 "running concurrently and are therefore upper bounds, not "
                 "clean single-process measurements"),
    }


# -- driver -------------------------------------------------------------------


def accuracy_compute_curve(region: dict, throughput: dict | None) -> list:
    """The accuracy-vs-compute tradeoff, cutoff by cutoff.

    **Diagnostic only.** The production cutoff stays frozen at
    ``FROZEN_ROUTER_CUTOFF``; this curve exists to show what the frozen point
    buys and costs relative to its neighbours, not to choose a new one, and no
    row of it may be promoted to a production threshold.
    """
    rows = cutoff_sweep(region)
    if throughput is None:
        return rows
    b_rate = throughput["binomial_loci_per_second"]
    pb_rate = throughput["pb_loci_per_second"]
    n = int(scoring_frame(region).sum())
    pb_only_seconds = n / pb_rate
    for row in rows:
        pb_loci = row["fraction_routed_to_pb"] * n
        row["pb_compute_fraction"] = row["fraction_routed_to_pb"]
        row["projected_cascade_seconds"] = n / b_rate + pb_loci / pb_rate
        row["projected_pb_only_seconds"] = pb_only_seconds
        row["projected_speedup"] = pb_only_seconds / max(row["projected_cascade_seconds"], 1e-12)
    return rows


#: Throughput probes, measured once each in a quiet process after extraction.
THROUGHPUT_PROBES = {
    "full_depth": {
        "bam": "data/giab_hg002_real_train/hg002_chr21_31_33M.bam",
        "vcf": "data/giab_hg002_real_train/hg002_chr21_31_33M.vcf.gz",
        "bed": "data/giab_hg002_real_train/hg002_chr21_31_33M_highconf.bed",
        "span": (31_000_000, 31_128_000)},
    "30x": {
        "bam": "data/giab_hg002_real_30x/hg002_chr21_31_33M_30x.bam",
        "vcf": "data/giab_hg002_real_train/hg002_chr21_31_33M.vcf.gz",
        "bed": "data/giab_hg002_real_train/hg002_chr21_31_33M_highconf.bed",
        "span": (31_000_000, 31_128_000)},
    "15x": {
        "bam": "data/giab_hg002_real_15x/hg002_chr21_31_33M_15x.bam",
        "vcf": "data/giab_hg002_real_train/hg002_chr21_31_33M.vcf.gz",
        "bed": "data/giab_hg002_real_train/hg002_chr21_31_33M_highconf.bed",
        "span": (31_000_000, 31_128_000)},
}


def probe_for(region_name: str) -> str:
    """Match a region to the throughput probe of its own depth regime."""
    if region_name.endswith("full"):
        return "full_depth"
    return "30x" if region_name.endswith("30x") else "15x"


def score_region(region: dict, tainted: bool, throughput: dict | None = None) -> dict:
    entry = evaluate_arms(region)
    entry["router_tuning_region_overlap"] = tainted
    entry["by_depth"] = depth_stratified(region)
    entry["by_depth_fine"] = fine_depth_stratified(region)
    entry["by_quality"] = quality_stratified(region)
    entry["by_vaf"] = vaf_stratified(region)
    entry["controls"] = controls(region)
    entry["strata"] = {label: arms_on_mask(region, mask)
                       for label, mask in strata_masks(region).items()}
    entry["accuracy_compute_curve"] = accuracy_compute_curve(region, throughput)
    if throughput is not None:
        entry["projected_wallclock"] = projected_wallclock(region, throughput)
    entry["pb_fix_scope"] = pb_fix_scope(region)
    ab = pb_fix_ab(region)
    if ab is not None:
        entry["pb_fix_ab"] = ab
    return entry


def main(measure: bool = True) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    report = {
        "frozen": {"binomial_threshold": FROZEN_BINOMIAL_THRESHOLD,
                   "pb_threshold": FROZEN_PB_THRESHOLD,
                   "router_cutoff": FROZEN_ROUTER_CUTOFF},
        "preregistration": "History/13_DEVLOG.md sections 1-7",
        "new_regions": {}, "cached_regions": {}, "compute": {},
    }
    throughput: dict[str, dict] = {}
    if measure:
        for probe, spec in THROUGHPUT_PROBES.items():
            logger.info("throughput probe: %s", probe)
            throughput[probe] = measure_throughput(spec["bam"], spec["vcf"], spec["bed"],
                                                   spec["span"])
    report["compute"]["throughput"] = throughput
    report["compute"]["throughput_note"] = (
        "measured in this process; the benchmark is run only after all "
        "extractions have exited, so no sibling job competes for the CPU")

    disagreements: dict[str, dict] = {}
    ab_report: dict[str, dict] = {}
    pools: dict[str, dict[str, list]] = {}

    def add_to_pool(pool: str, region: dict) -> None:
        frame = scoring_frame(region)
        bucket = pools.setdefault(pool, {"snp": [], "binomial": [], "pb": []})
        bucket["snp"].append(region["labels"][frame] == LABEL_SNP)
        bucket["binomial"].append(region["binomial_llr"][frame])
        bucket["pb"].append(region["pb_llr"][frame])

    for name, spec in NEW_REGIONS.items():
        if not Path(spec["path"]).exists():
            logger.warning("MISSING extraction %s (%s) -- skipped", name, spec["path"])
            continue
        logger.info("scoring %s", name)
        region = load_new(name, spec["path"])
        report["new_regions"][name] = score_region(
            region, spec["tainted"], throughput.get(probe_for(name)))
        disagreements[name] = disagreement_records(region)
        if "pb_fix_ab" in report["new_regions"][name]:
            ab_report[name] = report["new_regions"][name]["pb_fix_ab"]
        add_to_pool("group_a_all", region)
        add_to_pool("all_regions", region)
        if not spec["tainted"]:
            add_to_pool("group_a_untainted", region)
        if region["depth"].mean() > 48:
            add_to_pool("high_depth", region)
        del region

    for name, region in load_all_regions().items():
        logger.info("scoring cached %s", name)
        region["name"] = name
        report["cached_regions"][name] = score_region(
            region, tainted=False, throughput=throughput.get(probe_for(name)))
        disagreements[name] = disagreement_records(region)
        add_to_pool("group_b_cached", region)
        add_to_pool("all_regions", region)
        del region

    report["pooled"] = {}
    for pool, bucket in pools.items():
        report["pooled"][pool] = pooled_block(np.concatenate(bucket["snp"]),
                                              np.concatenate(bucket["binomial"]),
                                              np.concatenate(bucket["pb"]))
    report["failure_boundary"] = pooled_failure_boundary(disagreements)
    report["disagreement_summary"] = {
        name: {"loci_scored": entry["loci_scored"],
               "disagreements": entry["disagreements"],
               "disagreement_rate": entry["disagreement_rate"]}
        for name, entry in disagreements.items()}

    RESULTS.write_text(json.dumps(report, indent=2, default=float))
    DISAGREEMENTS.write_text(json.dumps(disagreements, indent=2, default=float))
    PB_FIX_AB.write_text(json.dumps(ab_report, indent=2, default=float))
    logger.info("Wrote %s, %s, %s", RESULTS, DISAGREEMENTS, PB_FIX_AB)


if __name__ == "__main__":
    main()
