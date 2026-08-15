"""bench_v12 Stage 1 — frozen cheap-router -> PB on new data and new hard strata.

Protocol is pre-registered in ``History/12_DEVLOG.md`` §6.1 and is not
re-derived here. Nothing in this file fits, selects or modifies a threshold:
the three frozen constants are imported from ``cascade`` and only ever read.

Two data families are evaluated:

* **new** — chr21:31.0-32.0 Mb at three depth realisations (15x/30x/full),
  extracted by ``extract_bench_v12.py``. Never scored by any previous devlog
  and disjoint from the router's 32-40 Mb cutoff-selection region.
* **hard strata** — four pre-registered subsets (low depth, PB-uncertain, low
  VAF, low base quality) cut from the seven already-cached regions of
  ``results/robustness_v11``. New *cuts*, not new regions.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import numpy as np

from binomial_baseline import BinomialVariantCaller
from cascade import (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD,
                     FROZEN_ROUTER_CUTOFF)
from config import LABEL_SNP
from evaluate_quality_error import paired_bootstrap
from residual_metrics import prf, vaf_from_counts
from robustness_benchmark import (CUTOFF_SWEEP, accuracy_block, controls,
                                  cutoff_sweep, depth_stratified, evaluate_arms,
                                  load_all_regions, quality_stratified, route_mask)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

NEW_REGIONS = {
    "new_31_32M_15x": "cache/bench_v12/new31_32M_15x.npz",
    "new_31_32M_30x": "cache/bench_v12/new31_32M_30x.npz",
    "new_31_32M_full": "cache/bench_v12/new31_32M_full.npz",
}

OUT = "results/bench_v12/stage1_results.json"


def load_new(name: str, path: str) -> dict:
    blob = np.load(path, allow_pickle=True)
    counts = blob["counts"]
    depth = blob["depth"].astype(np.float64)
    return {
        "name": name,
        "counts": counts,
        "binomial_llr": blob["binomial_llr"].astype(np.float64),
        "pb_llr": blob["pb_llr"].astype(np.float64),
        "labels": blob["labels"],
        "depth": depth,
        "mean_base_quality": np.where(depth > 0, counts[:, 7] / np.maximum(depth, 1), 0.0),
        "mean_mapq": np.where(depth > 0, counts[:, 8] / np.maximum(depth, 1), 0.0),
        "region": [int(x) for x in blob["region"]],
        "pb_time_s": float(blob["timing"][2]),
        "counts_time_s": float(blob["timing"][0]),
        "reads_time_s": float(blob["timing"][1]),
    }


def strata_masks(region: dict) -> dict:
    """Pre-registered hard-stratum masks. Defined on evidence only, never labels."""
    counts = region["counts"]
    depth = region["depth"]
    vaf = vaf_from_counts(counts)
    quality = region["mean_base_quality"]
    return {
        "hard_low_depth": depth <= 9,
        "hard_pb_uncertain": np.abs(region["pb_llr"] - FROZEN_PB_THRESHOLD) <= 13.5,
        "hard_low_vaf": (vaf > 0) & (vaf <= 0.35),
        "hard_low_quality": (depth > 0) & (quality < 30),
    }


def arms_on_mask(region: dict, mask: np.ndarray) -> dict:
    """All three arms restricted to a boolean subset of the scoring frame."""
    labels = region["labels"]
    frame = ((labels == 0) | (labels == LABEL_SNP)) & mask
    snp = labels[frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]

    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    router_calls = np.where(routed, pb_calls, binomial_calls)

    entry = {
        "loci": int(frame.sum()), "snp": int(snp.sum()),
        "binomial_only": accuracy_block(binomial_calls, snp),
        "pb_only": accuracy_block(pb_calls, snp),
        "cheap_router_pb": accuracy_block(router_calls, snp),
        "fraction_routed_to_pb": float(routed.mean()) if frame.sum() else 0.0,
        "pb_compute_loci": int(routed.sum()),
        "disagreement_vs_pb": int(np.count_nonzero(router_calls != pb_calls)),
    }
    entry["delta_f1_vs_pb"] = entry["cheap_router_pb"]["f1"] - entry["pb_only"]["f1"]
    if frame.sum():
        entry["vs_pb_paired_bootstrap"] = paired_bootstrap(router_calls, pb_calls, snp)
    return entry


def measure_throughput(region: dict) -> dict:
    """Same-session binomial and Poisson-binomial throughput on the new region.

    Devlog 11 §13 had to mix a fresh binomial rate with a cross-session PB rate.
    Here both are measured in this process, on this machine, on a slice of the
    very region being scored.
    """
    from extract_bench_v12 import poisson_binomial_for_chunk
    from read_level_pileup import load_reads
    from pileup_counts import load_counts

    fasta = "data/reference/chr21_full.fa"
    bam = ("data/giab_hg002_real_15x/hg002_chr21_31_33M_15x.bam"
           if region["name"].endswith("15x") else
           "data/giab_hg002_real_30x/hg002_chr21_31_33M_30x.bam"
           if region["name"].endswith("30x") else
           "data/giab_hg002_real_train/hg002_chr21_31_33M.bam")
    vcf = "data/giab_hg002_real_train/hg002_chr21_31_33M.vcf.gz"
    bed = "data/giab_hg002_real_train/hg002_chr21_31_33M_highconf.bed"
    span = (31_000_000, 31_128_000)

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
        "sample_loci": int(n),
        "binomial_seconds": binomial_seconds,
        "binomial_loci_per_second": n / max(binomial_seconds, 1e-9),
        "pb_seconds": pb_seconds,
        "pb_loci_per_second": n / max(pb_seconds, 1e-9),
        "note": "both measured in the same process/session on the new region",
    }


def projected_wallclock(region: dict, throughput: dict) -> dict:
    labels = region["labels"]
    frame = (labels == 0) | (labels == LABEL_SNP)
    routed = route_mask(region["binomial_llr"][frame], FROZEN_ROUTER_CUTOFF)
    n = int(frame.sum())
    k = int(routed.sum())
    b_rate = throughput["binomial_loci_per_second"]
    pb_rate = throughput["pb_loci_per_second"]
    pb_only = n / pb_rate
    cascade = n / b_rate + k / pb_rate
    return {
        "loci": n, "pb_loci": k,
        "pb_only_seconds": pb_only,
        "cascade_seconds": cascade,
        "speedup": pb_only / max(cascade, 1e-12),
        "measured_extraction_pb_seconds": region["pb_time_s"],
    }


def main() -> None:
    report = {
        "frozen": {"binomial_threshold": FROZEN_BINOMIAL_THRESHOLD,
                   "pb_threshold": FROZEN_PB_THRESHOLD,
                   "router_cutoff": FROZEN_ROUTER_CUTOFF},
        "new_regions": {}, "hard_strata": {}, "compute": {},
    }

    # ---- new region, three depth realisations -------------------------------
    pooled = {"labels": [], "binomial": [], "pb": []}
    for name, path in NEW_REGIONS.items():
        logger.info("Stage 1: %s", name)
        region = load_new(name, path)
        entry = evaluate_arms(region)
        entry["by_depth"] = depth_stratified(region)
        entry["by_quality"] = quality_stratified(region)
        entry["cutoff_sweep"] = cutoff_sweep(region)
        entry["controls"] = controls(region)
        entry["strata"] = {label: arms_on_mask(region, mask)
                           for label, mask in strata_masks(region).items()}
        entry["extraction_timing_seconds"] = {
            "counts": region["counts_time_s"], "reads": region["reads_time_s"],
            "pb": region["pb_time_s"]}
        report["new_regions"][name] = entry

        throughput = measure_throughput(region)
        report["compute"][name] = {"throughput": throughput,
                                   "projection": projected_wallclock(region, throughput)}

        labels = region["labels"]
        frame = (labels == 0) | (labels == LABEL_SNP)
        pooled["labels"].append(labels[frame] == LABEL_SNP)
        pooled["binomial"].append(region["binomial_llr"][frame])
        pooled["pb"].append(region["pb_llr"][frame])
        del region

    snp = np.concatenate(pooled["labels"])
    binomial_llr = np.concatenate(pooled["binomial"])
    pb_llr = np.concatenate(pooled["pb"])
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    router_calls = np.where(routed, pb_calls, binomial_calls)
    report["new_pooled"] = {
        "loci": int(snp.size), "snp": int(snp.sum()),
        "binomial_only": accuracy_block(binomial_calls, snp),
        "pb_only": accuracy_block(pb_calls, snp),
        "cheap_router_pb": accuracy_block(router_calls, snp),
        "fraction_routed_to_pb": float(routed.mean()),
        "disagreement_vs_pb": int(np.count_nonzero(router_calls != pb_calls)),
        "vs_pb_paired_bootstrap": paired_bootstrap(router_calls, pb_calls, snp),
        "delta_f1_vs_pb": prf(router_calls, snp)["f1"] - prf(pb_calls, snp)["f1"],
    }

    # ---- new hard strata over the seven devlog-11 regions --------------------
    old = load_all_regions()
    strata_pool: dict[str, dict] = {}
    for name, region in old.items():
        logger.info("Stage 1 strata: %s", name)
        masks = strata_masks(region)
        report["hard_strata"][name] = {label: arms_on_mask(region, mask)
                                       for label, mask in masks.items()}
        labels = region["labels"]
        frame = (labels == 0) | (labels == LABEL_SNP)
        for label, mask in masks.items():
            bucket = strata_pool.setdefault(label, {"snp": [], "binomial": [], "pb": []})
            sel = frame & mask
            bucket["snp"].append(labels[sel] == LABEL_SNP)
            bucket["binomial"].append(region["binomial_llr"][sel])
            bucket["pb"].append(region["pb_llr"][sel])
        del region

    report["hard_strata_pooled"] = {}
    for label, bucket in strata_pool.items():
        snp = np.concatenate(bucket["snp"])
        binomial_llr = np.concatenate(bucket["binomial"])
        pb_llr = np.concatenate(bucket["pb"])
        binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
        pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
        routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
        router_calls = np.where(routed, pb_calls, binomial_calls)
        report["hard_strata_pooled"][label] = {
            "loci": int(snp.size), "snp": int(snp.sum()),
            "binomial_only": accuracy_block(binomial_calls, snp),
            "pb_only": accuracy_block(pb_calls, snp),
            "cheap_router_pb": accuracy_block(router_calls, snp),
            "fraction_routed_to_pb": float(routed.mean()),
            "disagreement_vs_pb": int(np.count_nonzero(router_calls != pb_calls)),
            "vs_pb_paired_bootstrap": paired_bootstrap(router_calls, pb_calls, snp),
        }

    Path(OUT).parent.mkdir(parents=True, exist_ok=True)
    Path(OUT).write_text(json.dumps(report, indent=2, default=float))
    logger.info("Wrote %s", OUT)


if __name__ == "__main__":
    main()
