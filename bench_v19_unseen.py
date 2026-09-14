"""bench_v19 — Validation Round 2: frozen cascade on HG004, three new contigs.

Protocol is pre-registered in ``VALIDATION_ROUND_2.md`` and is not re-derived
here. Nothing in this file fits, selects, or modifies a threshold: the three
frozen constants are imported from ``cascade`` and only ever read. Every
metric helper (arms_on_mask, scoring_frame, pooled_block, disagreement
machinery, controls, accuracy-vs-compute curve) is imported unchanged from
v12/v13/v14, so v19's numbers are produced by byte-identical scoring code to
every earlier devlog. Two things are genuinely new in this file, both
pre-registered in VALIDATION_ROUND_2.md §7 and neither touching a frozen
file: a genomic-block bootstrap (as opposed to locus-level) and McNemar's
test on the paired cascade/PB disagreement table.

Three regions (chr2 ordinary, chr3 difficult, chr5 stress), three depths each,
one sample (HG004), never previously used by this repository.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

from bench_v12_stage1 import arms_on_mask, strata_masks
from bench_v13_robustness import (disagreement_records, fine_depth_stratified,
                                  pooled_block, pooled_failure_boundary,
                                  scoring_frame, vaf_stratified)
from bench_v14_crosschrom import (STRATIFICATIONS, accuracy_compute_curve, classify,
                                  load_cell, projected_wallclock, structure_masks,
                                  structure_stratified)
from binomial_baseline import BinomialVariantCaller
from cascade import (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD,
                     FROZEN_ROUTER_CUTOFF)
from config import LABEL_SNP
from robustness_benchmark import (accuracy_block, controls, depth_stratified,
                                  evaluate_arms, quality_stratified, route_mask)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

CHUNK_BP = 500_000
SEQ_LEN = 64
BLOCK_BP = 50_000  # block-bootstrap block size, fixed in VALIDATION_ROUND_2.md §7.2

BASE_REGIONS = {
    "chr2_ordinary": {
        "contig": "chr2", "start": 210_000_000, "stop": 213_000_000,
        "regime": "Rule R19-ordinary: window closest to contig-wide alldifficult_hc_fraction",
        "role": "ordinary",
    },
    "chr3_difficult": {
        "contig": "chr3", "start": 75_000_000, "stop": 78_000_000,
        "regime": "Rule R19-difficult: window maximising lowmap_segdup_hc_fraction",
        "role": "difficult",
    },
    "chr5_stress": {
        "contig": "chr5", "start": 27_000_000, "stop": 30_000_000,
        "regime": "Rule R19-stress: window maximising tandemrepeats_hc_fraction (new axis)",
        "role": "stress",
    },
}
DEPTH_TAGS = ("full", "30x", "15x")

PREVIOUSLY_USED_CONTIGS = {"chr21", "chr20", "chr19", "chr4", "chr1",
                          "chr16", "chr15", "chr7", "chr14", "chr13", "chr17"}
EXPECTED_SAMPLE = "HG004"
EXPECTED_BAM_TOKEN = "HG004.GRCh38.2x250"

DATA_DIR = Path("data/giab_hg004_v19")
CACHE_DIR = Path("cache/bench_v19")
OUT_DIR = Path("results/bench_v19")
RESULTS = OUT_DIR / "benchmark_results.json"
DISAGREEMENTS = OUT_DIR / "disagreements.json"
AUDIT = OUT_DIR / "independence_audit.json"


def cells() -> dict[str, dict]:
    return {f"{name}_{tag}": {**base, "base": name, "depth_tag": tag,
                              "path": str(CACHE_DIR / f"{name}_{tag}.npz")}
            for name, base in BASE_REGIONS.items() for tag in DEPTH_TAGS}


# -- coordinates --------------------------------------------------------------


def genomic_positions(spec: dict) -> np.ndarray:
    from pileup_counts import PileupCountsProvider
    positions = []
    for chunk_start in range(spec["start"], spec["stop"], CHUNK_BP):
        chunk_stop = min(chunk_start + CHUNK_BP, spec["stop"])
        provider = PileupCountsProvider(
            fasta_path=f"data/reference/{spec['contig']}_full.fa",
            bam_path=str(DATA_DIR / f"{spec['base']}_{spec['depth_tag']}.bam"),
            vcf_path=str(DATA_DIR / f"{spec['base']}.vcf.gz"),
            contig=spec["contig"], region=(chunk_start, chunk_stop), seq_len=SEQ_LEN,
            high_confidence_bed=str(DATA_DIR / f"{spec['base']}_highconf.bed"),
        )
        with provider:
            for window in provider._build_windows():
                positions.append(np.arange(window.start, window.end, dtype=np.int64))
    return np.concatenate(positions) if positions else np.zeros(0, dtype=np.int64)


# -- independence audit --------------------------------------------------------


def _sample_is_hg004(bam_path: str, base_name: str) -> bool:
    """Two-signal identity check, same pattern as devlog 18 (upstream GIAB BAMs
    carry no @RG SM tag): truth VCF sample column + BAM @PG slice provenance."""
    import subprocess

    import pysam

    vcf_path = str(Path(bam_path).parent / f"{base_name}.vcf.gz")
    try:
        out = subprocess.run(["bcftools", "query", "-l", vcf_path],
                             capture_output=True, text=True, check=True)
        vcf_ok = out.stdout.strip() == EXPECTED_SAMPLE
    except Exception:
        vcf_ok = False
    try:
        with pysam.AlignmentFile(bam_path, "rb") as fh:
            pg_lines = fh.header.get("PG", [])
        bam_ok = any(EXPECTED_BAM_TOKEN in pg.get("CL", "") for pg in pg_lines)
    except Exception:
        bam_ok = False
    return bool(vcf_ok and bam_ok)


def verify_caches(loaded: dict[str, dict]) -> dict:
    audit = {"previously_used_contigs": sorted(PREVIOUSLY_USED_CONTIGS), "cells": {}}
    for key, region in loaded.items():
        spec = region["spec"]
        labels = region["labels"]
        positions = region["positions"]
        entry = {
            "contig": region["contig"], "region": region["region"],
            "expected_region": [spec["start"], spec["stop"]],
            "loci_cached": int(labels.size),
            "loci_from_window_tiling": int(positions.size),
            "window_aligned": bool(labels.size % SEQ_LEN == 0),
            "snp": int((labels == LABEL_SNP).sum()),
            "non_finite_binomial": int((~np.isfinite(region["binomial_llr"])).sum()),
            "non_finite_pb": int((~np.isfinite(region["pb_llr"])).sum()),
            "positions_inside_region": bool(
                positions.size and positions.min() >= spec["start"]
                and positions.max() < spec["stop"]),
            "contig_previously_used": region["contig"] in PREVIOUSLY_USED_CONTIGS,
            "bam": region["bam"],
            "sample_is_hg004": _sample_is_hg004(region["bam"], spec["base"]),
        }
        entry["ok"] = bool(
            entry["contig"] == spec["contig"] and entry["region"] == entry["expected_region"]
            and entry["window_aligned"] and entry["loci_cached"] == entry["loci_from_window_tiling"]
            and entry["non_finite_binomial"] == 0 and entry["non_finite_pb"] == 0
            and entry["positions_inside_region"] and not entry["contig_previously_used"]
            and entry["sample_is_hg004"])
        audit["cells"][key] = entry
    audit["all_ok"] = all(c["ok"] for c in audit["cells"].values())
    audit["n_cells"] = len(audit["cells"])
    return audit


# -- new statistics: block bootstrap and McNemar -------------------------------


def block_bootstrap(cascade_calls: np.ndarray, pb_calls: np.ndarray, truth: np.ndarray,
                    positions: np.ndarray, block_bp: int = BLOCK_BP,
                    resamples: int = 2000, seed: int = 20260904) -> dict:
    """Block bootstrap for ΔF1 = F1(cascade) - F1(PB), resampling contiguous
    genomic blocks (not individual loci) with replacement.

    Only "interesting" loci (either arm positive, or truth positive) are kept
    per block, mirroring ``evaluate_quality_error.paired_bootstrap``'s
    inert-locus optimisation, but grouped by block so within-block correlation
    is preserved by construction: an entire block is included or excluded as a
    unit in each resample, never split.
    """
    interesting = cascade_calls | pb_calls | truth
    idx = np.nonzero(interesting)[0]
    block_id = (positions[idx] // block_bp).astype(np.int64)
    unique_blocks = np.unique(block_id)
    n_blocks = unique_blocks.size

    # Per-block (tp, fp, fn) for each arm, precomputed once.
    per_block = {}
    for b in unique_blocks:
        sel = idx[block_id == b]
        c, p, t = cascade_calls[sel], pb_calls[sel], truth[sel]
        per_block[b] = (
            int((c & t).sum()), int((c & ~t).sum()), int((~c & t).sum()),
            int((p & t).sum()), int((p & ~t).sum()), int((~p & t).sum()),
        )
    table = np.array([per_block[b] for b in unique_blocks])  # [n_blocks, 6]

    def f1(tp, fp, fn):
        denom = 2 * tp + fp + fn
        return np.where(denom > 0, 2 * tp / np.maximum(denom, 1), 0.0)

    rng = np.random.default_rng(seed)
    deltas = np.empty(resamples)
    for r in range(resamples):
        draw = rng.integers(0, n_blocks, size=n_blocks)
        summed = table[draw].sum(axis=0)
        tp_c, fp_c, fn_c, tp_p, fp_p, fn_p = summed
        deltas[r] = f1(tp_c, fp_c, fn_c) - f1(tp_p, fp_p, fn_p)

    point_tp_c, point_fp_c, point_fn_c = table[:, 0].sum(), table[:, 1].sum(), table[:, 2].sum()
    point_tp_p, point_fp_p, point_fn_p = table[:, 3].sum(), table[:, 4].sum(), table[:, 5].sum()
    point_delta = float(f1(point_tp_c, point_fp_c, point_fn_c) - f1(point_tp_p, point_fp_p, point_fn_p))
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    return {
        "n_blocks": int(n_blocks), "block_bp": block_bp,
        "delta_f1": point_delta, "ci95": [float(lo), float(hi)],
        "ci_half_width": float((hi - lo) / 2), "resamples": resamples,
    }


def mcnemar_test(cascade_calls: np.ndarray, pb_calls: np.ndarray, truth: np.ndarray) -> dict:
    """McNemar's test on the paired cascade-vs-PB disagreement table.

    Cell b: cascade wrong, PB right. Cell c: cascade right, PB wrong (both
    relative to truth). Only loci where the two arms disagree contribute.
    Uses the exact binomial form (appropriate at the small b+c counts this
    project has always seen), not the chi-square approximation.
    """
    cascade_correct = cascade_calls == truth
    pb_correct = pb_calls == truth
    b = int(np.sum(~cascade_correct & pb_correct))   # cascade wrong, PB right
    c = int(np.sum(cascade_correct & ~pb_correct))    # cascade right, PB wrong
    n = b + c
    if n == 0:
        return {"b_cascade_wrong_pb_right": 0, "c_cascade_right_pb_wrong": 0,
                "n_discordant": 0, "p_value": 1.0, "note": "no discordant pairs"}
    result = binomtest(min(b, c), n, 0.5, alternative="two-sided")
    return {"b_cascade_wrong_pb_right": b, "c_cascade_right_pb_wrong": c,
            "n_discordant": n, "p_value": float(result.pvalue)}


# -- disagreement enrichment (verbatim pattern from v17/v18) ------------------


def enrich_disagreements(region: dict, block: dict) -> dict:
    frame = scoring_frame(region)
    mapq = region["mean_mapq"][frame]
    positions = region["positions"][frame]
    labels = {name: mask[frame]
              for name, mask in structure_masks(region["spec"], region["positions"]).items()
              if name.startswith("in_")}
    for record in block["records"]:
        index = record["frame_index"]
        record["mean_mapq"] = float(mapq[index])
        record["position"] = int(positions[index])
        record["annotations"] = sorted(name[3:] for name, mask in labels.items()
                                       if bool(mask[index]))
    counts: dict[str, int] = {}
    for record in block["records"]:
        counts[record["effect"]] = counts.get(record["effect"], 0) + 1
        for name in record["annotations"] or ["none"]:
            counts[f"in_{name}"] = counts.get(f"in_{name}", 0) + 1
    block["record_summary"] = counts
    return block


# -- wall clock (caller-only throughput, for the projected figure) -----------


def measure_throughput(spec: dict, span_bp: int = 128_000) -> dict:
    from extract_bench_v12 import poisson_binomial_for_chunk
    from pileup_counts import load_counts
    from read_level_pileup import load_reads

    fasta = f"data/reference/{spec['contig']}_full.fa"
    bam = str(DATA_DIR / f"{spec['base']}_{spec['depth_tag']}.bam")
    vcf = str(DATA_DIR / f"{spec['base']}.vcf.gz")
    bed = str(DATA_DIR / f"{spec['base']}_highconf.bed")
    span = (spec["start"], spec["start"] + span_bp)

    counts, _ = load_counts(fasta, bam, vcf, bed, span, contig=spec["contig"])
    reads, _ = load_reads(fasta, bam, vcf, bed, span, contig=spec["contig"])

    caller = BinomialVariantCaller(error_rate=0.01)
    start = time.perf_counter()
    caller.score_counts(counts[:, 0:4], counts[:, 9].astype(int))
    binomial_seconds = time.perf_counter() - start

    start = time.perf_counter()
    poisson_binomial_for_chunk(reads, counts)
    pb_seconds = time.perf_counter() - start

    n = counts.shape[0]
    return {"bam": bam, "span": list(span), "sample_loci": int(n),
            "mean_depth": float(counts[:, 6].mean()),
            "binomial_seconds": binomial_seconds,
            "binomial_loci_per_second": n / max(binomial_seconds, 1e-9),
            "pb_seconds": pb_seconds,
            "pb_loci_per_second": n / max(pb_seconds, 1e-9)}


# -- sensitivity sweep: ±10/20/30% of the frozen cutoff -----------------------


def cutoff_sensitivity(region: dict) -> dict:
    frame = scoring_frame(region)
    snp = region["labels"][frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    out = {}
    for pct in (-0.30, -0.20, -0.10, 0.0, 0.10, 0.20, 0.30):
        cutoff = FROZEN_ROUTER_CUTOFF * (1.0 + pct)
        routed = route_mask(binomial_llr, cutoff)
        calls = np.where(routed, pb_calls, binomial_calls)
        stats = accuracy_block(calls, snp)
        f1 = stats["f1"] if "f1" in stats else (
            2 * stats["tp"] / max(2 * stats["tp"] + stats["fp"] + stats["fn"], 1))
        out[f"{pct:+.0%}"] = {
            "cutoff": cutoff, "f1": float(f1), "fp": int(stats["fp"]), "fn": int(stats["fn"]),
            "routed_fraction": float(routed.mean()), "routed_loci": int(routed.sum()),
        }
    return out


# -- driver -------------------------------------------------------------------


def score_cell(region: dict, throughput: dict | None) -> dict:
    entry = evaluate_arms(region)
    entry["contig"] = region["contig"]
    entry["depth_tag"] = region["spec"]["depth_tag"]
    entry["role"] = region["spec"]["role"]
    entry["mean_depth"] = float(region["depth"].mean())
    entry["max_depth"] = float(region["depth"].max())
    entry["by_depth"] = depth_stratified(region)
    entry["by_depth_fine"] = fine_depth_stratified(region)
    entry["by_quality"] = quality_stratified(region)
    entry["by_vaf"] = vaf_stratified(region)
    entry["by_structure"] = structure_stratified(region)
    entry["strata"] = {label: arms_on_mask(region, mask)
                       for label, mask in strata_masks(region).items()}
    entry["controls"] = controls(region)
    entry["accuracy_compute_curve"] = accuracy_compute_curve(region, throughput)
    entry["cutoff_sensitivity"] = cutoff_sensitivity(region)
    if throughput is not None:
        entry["projected_wallclock"] = projected_wallclock(region, throughput)
    entry["verdict"] = classify(entry)

    frame = scoring_frame(region)
    snp = region["labels"][frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    cascade_calls = np.where(routed, pb_llr >= FROZEN_PB_THRESHOLD, binomial_llr >= FROZEN_BINOMIAL_THRESHOLD)
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    positions = region["positions"][frame]
    entry["block_bootstrap"] = block_bootstrap(cascade_calls, pb_calls, snp, positions)
    entry["mcnemar"] = mcnemar_test(cascade_calls, pb_calls, snp)
    return entry


def _pool(regions: list[dict]) -> dict:
    frames = [(r, scoring_frame(r)) for r in regions]
    snp = np.concatenate([r["labels"][f] == LABEL_SNP for r, f in frames])
    binomial = np.concatenate([r["binomial_llr"][f] for r, f in frames])
    pb = np.concatenate([r["pb_llr"][f] for r, f in frames])
    entry = pooled_block(snp, binomial, pb)
    entry["verdict"] = classify(entry)

    routed = route_mask(binomial, FROZEN_ROUTER_CUTOFF)
    cascade_calls = np.where(routed, pb >= FROZEN_PB_THRESHOLD, binomial >= FROZEN_BINOMIAL_THRESHOLD)
    pb_calls = pb >= FROZEN_PB_THRESHOLD
    # Pooled block bootstrap needs a single monotonically-meaningful position
    # axis; offset each region's positions by a large constant per contig so
    # blocks never collide across regions (contigs are never adjacent).
    offset_positions = []
    for i, (r, f) in enumerate(frames):
        offset_positions.append(r["positions"][f] + i * 1_000_000_000)
    positions = np.concatenate(offset_positions)
    entry["block_bootstrap"] = block_bootstrap(cascade_calls, pb_calls, snp, positions)
    entry["mcnemar"] = mcnemar_test(cascade_calls, pb_calls, snp)
    return entry


def main(measure: bool = True) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    specs = cells()
    missing = [key for key, spec in specs.items() if not Path(spec["path"]).exists()]
    if missing:
        raise SystemExit(f"missing extractions, refusing to score a partial benchmark: {missing}")

    loaded = {}
    for key, spec in specs.items():
        logger.info("loading %s", key)
        region = load_cell(key, spec)
        region["positions"] = genomic_positions(spec)
        loaded[key] = region

    audit = verify_caches(loaded)
    AUDIT.write_text(json.dumps(audit, indent=2))
    logger.info("independence audit: all_ok=%s over %d cells", audit["all_ok"], audit["n_cells"])
    if not audit["all_ok"]:
        raise SystemExit(f"cache validation failed; see {AUDIT}")

    throughput: dict[str, dict] = {}
    if measure:
        for base_name in BASE_REGIONS:
            for tag in DEPTH_TAGS:
                key = f"{base_name}_{tag}"
                logger.info("throughput probe: %s", key)
                throughput[key] = measure_throughput(specs[key])

    report = {
        "frozen": {"binomial_threshold": FROZEN_BINOMIAL_THRESHOLD,
                   "pb_threshold": FROZEN_PB_THRESHOLD,
                   "router_cutoff": FROZEN_ROUTER_CUTOFF},
        "preregistration": "VALIDATION_ROUND_2.md",
        "region": BASE_REGIONS,
        "independence_audit": {"all_ok": audit["all_ok"], "path": str(AUDIT)},
        "compute": {"throughput": throughput,
                    "throughput_note": "measured in this process on the v19 BAMs, after every extraction had exited"},
        "cells": {}, "pooled": {},
    }

    disagreements: dict[str, dict] = {}
    for key, region in loaded.items():
        logger.info("scoring %s", key)
        report["cells"][key] = score_cell(region, throughput.get(key))
        disagreements[key] = enrich_disagreements(region, disagreement_records(region))

    report["pooled"]["all"] = _pool(list(loaded.values()))
    for tag in DEPTH_TAGS:
        report["pooled"][f"depth_{tag}"] = _pool([r for r in loaded.values() if r["spec"]["depth_tag"] == tag])
    for base_name in BASE_REGIONS:
        report["pooled"][f"region_{base_name}"] = _pool([r for r in loaded.values() if r["spec"]["base"] == base_name])
    report["failure_boundary"] = pooled_failure_boundary(disagreements)

    RESULTS.write_text(json.dumps(report, indent=2))
    DISAGREEMENTS.write_text(json.dumps(disagreements, indent=2))
    logger.info("wrote %s and %s", RESULTS, DISAGREEMENTS)


if __name__ == "__main__":
    main()
