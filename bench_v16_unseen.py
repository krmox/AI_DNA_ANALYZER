"""bench_v16 — the frozen cascade on one representative, completely unseen region.

Protocol is pre-registered in ``History/16_DEVLOG.md`` §1–§10 and is not
re-derived here. Nothing in this file fits, selects or modifies a threshold: the
three frozen constants are imported from ``cascade`` and only ever read, and the
cutoff sweep exists solely to draw the diagnostic accuracy-vs-compute curve.

Three cells are scored — one 3 Mb window on chr13, a contig no experiment in
this project has touched, at native / 30x / 15x depth. The window was chosen by
``select_v16_region.py`` from public annotation alone, before any caller ran and
before the chr13 FASTA existed.

Almost everything here is imported. What is specific to v16 is only:

* :func:`genomic_positions` / :func:`measure_throughput` / :func:`verify_caches`
  — the v14 versions of these are bound to v14's data directory and its
  ``PREVIOUSLY_USED_CONTIGS`` set, so they are re-expressed against v16's, and
* :func:`enrich_disagreements` — adds mapping quality, the true GRCh38
  coordinate and the external difficulty labels to each disagreement record,
  which devlog 16 §8 (phase 8) asks for and the shared helper does not carry.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import numpy as np

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
from robustness_benchmark import (controls, depth_stratified, evaluate_arms,
                                  quality_stratified, route_mask)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

#: Must match what ``run_v16_extract.sh`` used (``extract_bench_v12.CHUNK_BP``),
#: because the window tiling — and therefore which loci exist — is rebuilt per
#: chunk.
CHUNK_BP = 500_000
SEQ_LEN = 64

#: The one pre-registered region (devlog 16 §4.3), at three depths.
BASE_REGIONS = {
    "chr13_representative": {
        "contig": "chr13", "start": 70_000_000, "stop": 73_000_000,
        "regime": "representative chr13 sequence (median-difficulty contig, "
                  "contig-typical difficulty window)",
        "gc": 0.35249166666666665,
    },
}
DEPTH_TAGS = ("full", "30x", "15x")

#: Every contig any previous experiment used. The audit asserts no cell touches
#: one; Rule C16 could not have selected one in the first place.
PREVIOUSLY_USED_CONTIGS = {"chr21", "chr20", "chr19", "chr4", "chr1",
                           "chr16", "chr15", "chr7", "chr14"}

DATA_DIR = Path("data/giab_hg002_v16")
CACHE_DIR = Path("cache/bench_v16")
OUT_DIR = Path("results/bench_v16")
RESULTS = OUT_DIR / "benchmark_results.json"
DISAGREEMENTS = OUT_DIR / "disagreements.json"
AUDIT = OUT_DIR / "independence_audit.json"


def cells() -> dict[str, dict]:
    """The 3 pre-registered region/depth cells, keyed as they are cached."""
    return {f"{name}_{tag}": {**base, "base": name, "depth_tag": tag,
                              "path": str(CACHE_DIR / f"{name}_{tag}.npz")}
            for name, base in BASE_REGIONS.items() for tag in DEPTH_TAGS}


# -- coordinates --------------------------------------------------------------


def genomic_positions(spec: dict) -> np.ndarray:
    """True GRCh38 coordinate of every locus in a cell, in cache order.

    The caches store ``positions`` as ``chunk_start + arange(...)``, which is not
    a coordinate whenever a window was dropped as uncallable (devlog 12 §16.5).
    The tiling is therefore rebuilt with the same provider, chunking and BED that
    produced the cache. Tiling reads the FASTA and the BED only — no pileup.
    """
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


# -- validation and leakage audit ---------------------------------------------


def verify_caches(loaded: dict[str, dict]) -> dict:
    """Mechanical validation and independence audit, run before any metric."""
    audit = {"previously_used_contigs": sorted(PREVIOUSLY_USED_CONTIGS),
             "historical_regions": HISTORICAL_REGIONS, "cells": {}}
    for key, region in loaded.items():
        spec = region["spec"]
        labels = region["labels"]
        positions = region["positions"]
        overlaps = [f"{c}:{lo}-{hi}" for c, lo, hi in HISTORICAL_REGIONS
                    if c == spec["contig"] and lo < spec["stop"] and hi > spec["start"]]
        entry = {
            "contig": region["contig"],
            "region": region["region"],
            "expected_region": [spec["start"], spec["stop"]],
            "loci_cached": int(labels.size),
            "loci_from_window_tiling": int(positions.size),
            "window_aligned": bool(labels.size % SEQ_LEN == 0),
            "snp": int((labels == LABEL_SNP).sum()),
            "non_finite_binomial": int((~np.isfinite(region["binomial_llr"])).sum()),
            "non_finite_pb": int((~np.isfinite(region["pb_llr"])).sum()),
            "position_min": int(positions.min()) if positions.size else None,
            "position_max": int(positions.max()) if positions.size else None,
            "positions_inside_region": bool(
                positions.size and positions.min() >= spec["start"]
                and positions.max() < spec["stop"]),
            "overlaps_previously_used_contig": region["contig"] in PREVIOUSLY_USED_CONTIGS,
            "coordinate_overlaps_historical_regions": overlaps,
            "bam": region["bam"],
        }
        entry["ok"] = bool(
            entry["contig"] == spec["contig"]
            and entry["region"] == entry["expected_region"]
            and entry["window_aligned"]
            and entry["loci_cached"] == entry["loci_from_window_tiling"]
            and entry["non_finite_binomial"] == 0 and entry["non_finite_pb"] == 0
            and entry["positions_inside_region"]
            and not entry["overlaps_previously_used_contig"]
            and not overlaps)
        audit["cells"][key] = entry
    audit["all_ok"] = all(cell["ok"] for cell in audit["cells"].values())
    audit["n_cells"] = len(audit["cells"])
    return audit


#: Every region any previous devlog extracted, for the coordinate-level half of
#: the audit. Contig-level exclusion already makes overlap impossible; this is
#: the mechanical check that says so rather than assuming it.
HISTORICAL_REGIONS = [
    ("chr21", 5_010_000, 5_510_000), ("chr21", 10_000_000, 15_000_000),
    ("chr21", 29_999_999, 30_470_000), ("chr21", 31_000_000, 32_000_000),
    ("chr21", 32_000_000, 44_000_000),
    ("chr20", 33_000_000, 34_000_000), ("chr19", 1_000_000, 2_000_000),
    ("chr4", 101_000_000, 102_000_000), ("chr1", 120_000_000, 121_000_000),
    ("chr16", 29_000_000, 30_000_000), ("chr15", 30_000_000, 31_000_000),
    ("chr7", 57_000_000, 58_000_000), ("chr14", 53_000_000, 54_000_000),
]


# -- disagreement enrichment --------------------------------------------------


def enrich_disagreements(region: dict, block: dict) -> dict:
    """Add mapping quality, GRCh38 coordinate and difficulty labels per record.

    ``disagreement_records`` indexes the *scoring frame*; the structure masks and
    ``mean_mapq`` are defined over all cached loci, so the frame is applied to
    both before indexing. No metric is recomputed here.
    """
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


# -- wall clock ---------------------------------------------------------------


def measure_throughput(spec: dict, span_bp: int = 128_000) -> dict:
    """Binomial and repaired-PB throughput on this experiment's own data.

    Run only after every extraction has exited, so nothing competes for the CPU.
    """
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
    return {
        "bam": bam, "span": list(span), "sample_loci": int(n),
        "mean_depth": float(counts[:, 6].mean()),
        "binomial_seconds": binomial_seconds,
        "binomial_loci_per_second": n / max(binomial_seconds, 1e-9),
        "pb_seconds": pb_seconds,
        "pb_loci_per_second": n / max(pb_seconds, 1e-9),
    }


# -- driver -------------------------------------------------------------------


def score_cell(region: dict, throughput: dict | None) -> dict:
    entry = evaluate_arms(region)
    entry["contig"] = region["contig"]
    entry["depth_tag"] = region["spec"]["depth_tag"]
    entry["regime"] = region["spec"]["regime"]
    entry["gc_fraction"] = region["spec"]["gc"]
    entry["mean_depth"] = float(region["depth"].mean())
    entry["max_depth"] = float(region["depth"].max())
    entry["by_depth"] = depth_stratified(region)
    entry["by_depth_fine"] = fine_depth_stratified(region)
    entry["by_quality"] = quality_stratified(region)
    entry["by_vaf"] = vaf_stratified(region)
    entry["by_mapq"] = mapq_stratified(region)
    entry["by_structure"] = structure_stratified(region)
    entry["strata"] = {label: arms_on_mask(region, mask)
                       for label, mask in strata_masks(region).items()}
    entry["controls"] = controls(region)
    entry["accuracy_compute_curve"] = accuracy_compute_curve(region, throughput)
    if throughput is not None:
        entry["projected_wallclock"] = projected_wallclock(region, throughput)
    entry["verdict"] = classify(entry)
    return entry


def mapq_stratified(region: dict) -> dict:
    """All three arms binned by mean mapping quality.

    ``robustness_benchmark`` ships depth, base-quality and VAF stratifiers but no
    mapping-quality one, and devlog 16 §6 lists mapping quality among the axes to
    report.

    The edges originally mirrored the base-quality stratifier's and topped out at
    61, which was simply wrong for this data: novoalign emits MAPQ 70, so those
    bins covered 0.02% of loci. The top edge is 71 here. This is a defect fix in
    a **diagnostic** axis — no arm, threshold, endpoint or verdict depends on it
    — and it is disclosed in devlog 16 §17.
    """
    edges = [(0, 20), (20, 40), (40, 60), (60, 71)]
    frame = scoring_frame(region)
    mapq = region["mean_mapq"][frame]
    depth = region["depth"][frame]
    return {f"mq{lo}-{hi}": arms_on_mask(
        region, _expand(frame, (depth > 0) & (mapq >= lo) & (mapq < hi)))
        for lo, hi in edges}


def _expand(frame: np.ndarray, inner: np.ndarray) -> np.ndarray:
    """Lift a mask defined on the scoring frame back to all cached loci."""
    out = np.zeros(frame.size, dtype=bool)
    out[np.nonzero(frame)[0][inner]] = True
    return out


def _pool(regions: list[dict]) -> dict:
    frames = [(r, scoring_frame(r)) for r in regions]
    snp = np.concatenate([r["labels"][f] == LABEL_SNP for r, f in frames])
    binomial = np.concatenate([r["binomial_llr"][f] for r, f in frames])
    pb = np.concatenate([r["pb_llr"][f] for r, f in frames])
    entry = pooled_block(snp, binomial, pb)
    entry["verdict"] = classify(entry)
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
        for tag in DEPTH_TAGS:
            logger.info("throughput probe: %s", tag)
            throughput[tag] = measure_throughput(specs[f"chr13_representative_{tag}"])

    report = {
        "frozen": {"binomial_threshold": FROZEN_BINOMIAL_THRESHOLD,
                   "pb_threshold": FROZEN_PB_THRESHOLD,
                   "router_cutoff": FROZEN_ROUTER_CUTOFF},
        "preregistration": "History/16_DEVLOG.md sections 1-10",
        "region": BASE_REGIONS,
        "independence_audit": {"all_ok": audit["all_ok"], "path": str(AUDIT)},
        "compute": {"throughput": throughput,
                    "throughput_note": ("measured in this process on the v16 BAMs, after "
                                        "every extraction had exited")},
        "cells": {}, "pooled": {},
    }

    disagreements: dict[str, dict] = {}
    for key, region in loaded.items():
        logger.info("scoring %s", key)
        report["cells"][key] = score_cell(region, throughput.get(region["spec"]["depth_tag"]))
        disagreements[key] = enrich_disagreements(region, disagreement_records(region))

    report["pooled"]["all"] = _pool(list(loaded.values()))
    for tag in DEPTH_TAGS:
        report["pooled"][f"depth_{tag}"] = _pool(
            [r for r in loaded.values() if r["spec"]["depth_tag"] == tag])
    report["failure_boundary"] = pooled_failure_boundary(disagreements)

    RESULTS.write_text(json.dumps(report, indent=2))
    DISAGREEMENTS.write_text(json.dumps(disagreements, indent=2))
    logger.info("wrote %s and %s", RESULTS, DISAGREEMENTS)


if __name__ == "__main__":
    main()
