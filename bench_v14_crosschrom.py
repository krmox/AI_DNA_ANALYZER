"""bench_v14 — does the frozen cheap-router → PB cascade leave chr21?

Protocol is pre-registered in ``History/14_DEVLOG.md`` §1–§8 and is not
re-derived here. Nothing in this file fits, selects or modifies a threshold: the
three frozen constants are imported from ``cascade`` and only ever read, and the
cutoff sweep exists solely to draw the diagnostic accuracy-vs-compute curve.

Twelve cells are scored — four 1 Mb regions on four chromosomes the project has
never touched (chr20, chr19, chr4, chr1), each at native / 30x / 15x depth. The
regions were chosen by ``select_v14_regions.py`` from public annotation alone,
before any caller ran.

Every metric helper is imported from the existing benchmark modules rather than
reimplemented, so the numbers sit on the same definitions as devlogs 11–13. What
is genuinely new here is only:

* :func:`genomic_positions` — reconstructs each locus' true GRCh38 coordinate,
  which the caches do not store (devlog 12 §16.5: ``positions`` are cache-order
  indices, not coordinates), so that
* :func:`structure_masks` can label every locus with the external GIAB v3.1
  difficulty annotations, and
* :func:`verify_caches` can audit contig/coordinate independence mechanically.
"""

from __future__ import annotations

import gzip
import json
import logging
import time
from pathlib import Path

import numpy as np

from bench_v12_stage1 import arms_on_mask, strata_masks
from bench_v13_robustness import (VAF_BINS, _arms, _binned, disagreement_records,
                                  fine_depth_stratified, pooled_block,
                                  pooled_failure_boundary, scoring_frame,
                                  vaf_stratified)
from binomial_baseline import BinomialVariantCaller
from cascade import (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD,
                     FROZEN_ROUTER_CUTOFF)
from config import LABEL_SNP
from evaluate_binomial_baseline import prf
from robustness_benchmark import (accuracy_block, controls, cutoff_sweep,
                                  depth_stratified, evaluate_arms,
                                  quality_stratified, route_mask)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

#: Extraction chunk width. Must match what ``run_v14_extract.sh`` actually used
#: (``extract_bench_v12.CHUNK_BP``), because the window tiling — and therefore
#: which loci exist at all — is rebuilt per chunk.
CHUNK_BP = 500_000
SEQ_LEN = 64

#: The four pre-registered regions (devlog 14 §5.3), each at three depths.
BASE_REGIONS = {
    "chr20_neutral": {"contig": "chr20", "start": 33_000_000, "stop": 34_000_000,
                      "regime": "ordinary euchromatin (neutral control)", "gc": 0.477},
    "chr19_gcrich": {"contig": "chr19", "start": 1_000_000, "stop": 2_000_000,
                     "regime": "GC-rich, Alu/repeat-dense subtelomeric", "gc": 0.591},
    "chr4_atrich": {"contig": "chr4", "start": 101_000_000, "stop": 102_000_000,
                    "regime": "AT-rich, gene-poor isochore", "gc": 0.363},
    "chr1_segdup": {"contig": "chr1", "start": 120_000_000, "stop": 121_000_000,
                    "regime": "pericentromeric segmental duplication / low mappability",
                    "gc": 0.406},
}
DEPTH_TAGS = ("full", "30x", "15x")

#: External difficulty annotation, GIAB genome-stratifications v3.1 (devlog 14
#: §5.5). Fixed public BED files; not derived from anything measured here.
STRATIFICATIONS = {
    "lowmap_segdup": "data/strat_v31/lowmap_segdup.bed.gz",
    "tandemrepeats": "data/strat_v31/tandemrepeats.bed.gz",
    "alldifficult": "data/strat_v31/alldifficult.bed.gz",
}

#: Every region this project used before devlog 14 lay on chr21. The audit in
#: :func:`verify_caches` asserts no cell touches it.
PREVIOUSLY_USED_CONTIGS = {"chr21"}

DATA_DIR = Path("data/giab_hg002_v14")
CACHE_DIR = Path("cache/bench_v14")
OUT_DIR = Path("results/bench_v14")
RESULTS = OUT_DIR / "robustness_results.json"
DISAGREEMENTS = OUT_DIR / "disagreements.json"
AUDIT = OUT_DIR / "independence_audit.json"


def cells() -> dict[str, dict]:
    """The 12 pre-registered region/depth cells, keyed as they are cached."""
    out = {}
    for name, base in BASE_REGIONS.items():
        for tag in DEPTH_TAGS:
            out[f"{name}_{tag}"] = {**base, "base": name, "depth_tag": tag,
                                    "path": str(CACHE_DIR / f"{name}_{tag}.npz")}
    return out


# -- coordinates --------------------------------------------------------------


def genomic_positions(spec: dict) -> np.ndarray:
    """True GRCh38 coordinate of every locus in a cell, in cache order.

    The caches store ``positions`` as ``chunk_start + arange(n_loci_in_chunk)``,
    which is *not* a coordinate whenever a window was dropped as uncallable
    (devlog 12 §16.5). Rather than trust it, the window tiling is rebuilt here
    with the same provider, the same chunking and the same BED that produced the
    cache, and the loci are laid out from the surviving windows.

    Tiling reads the FASTA and the BED only — no pileup — so this is cheap, and
    the result is checked against the cache's own locus count by the caller.
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


def _bed_mask(bed_path: str, contig: str, start: int, stop: int) -> np.ndarray:
    """Boolean mask over ``[start, stop)`` of bases covered by a BED's intervals."""
    mask = np.zeros(stop - start, dtype=bool)
    wanted = {contig, contig[3:] if contig.startswith("chr") else "chr" + contig}
    opener = gzip.open if bed_path.endswith(".gz") else open
    with opener(bed_path, "rt") as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if fields[0] not in wanted:
                continue
            lo, hi = int(fields[1]), int(fields[2])
            if hi <= start or lo >= stop:
                continue
            mask[max(lo, start) - start:min(hi, stop) - start] = True
    return mask


def structure_masks(spec: dict, positions: np.ndarray) -> dict[str, np.ndarray]:
    """Per-locus external difficulty labels, and their complements.

    Defined entirely by public annotation and the locus' coordinate. No caller
    score, no truth label and no measurement from this benchmark participates.
    """
    out = {}
    offset = positions - spec["start"]
    for key, path in STRATIFICATIONS.items():
        covered = _bed_mask(path, spec["contig"], spec["start"], spec["stop"])
        inside = covered[offset]
        out[f"in_{key}"] = inside
        out[f"not_in_{key}"] = ~inside
    return out


# -- loading ------------------------------------------------------------------


def load_cell(key: str, spec: dict) -> dict:
    """Load one extracted cell into the region dict the shared helpers expect."""
    blob = np.load(spec["path"], allow_pickle=True)
    counts = blob["counts"]
    depth = blob["depth"].astype(np.float64)
    region = {
        "name": key,
        "counts": counts,
        "binomial_llr": blob["binomial_llr"].astype(np.float64),
        "pb_llr": blob["pb_llr"].astype(np.float64),
        "labels": blob["labels"],
        "depth": depth,
        "mean_base_quality": np.where(depth > 0, counts[:, 7] / np.maximum(depth, 1), 0.0),
        "mean_mapq": np.where(depth > 0, counts[:, 8] / np.maximum(depth, 1), 0.0),
        "region": [int(x) for x in blob["region"]],
        "contig": str(blob["contig"][0]) if "contig" in blob else None,
        "counts_time_s": float(blob["timing"][0]),
        "reads_time_s": float(blob["timing"][1]),
        "pb_time_s": float(blob["timing"][2]),
        "bam": str(blob["bam"][0]),
        "spec": spec,
    }
    return region


# -- validation and leakage audit ---------------------------------------------


def verify_caches(loaded: dict[str, dict]) -> dict:
    """Mechanical validation and independence audit, run before any metric.

    Asserts, per cell: the stored contig and coordinates match the
    pre-registration; the contig is not one this project has used before; the
    locus count is window-aligned; no LLR is non-finite; and the reconstructed
    genomic coordinates line up with the cached locus count and stay inside the
    region.
    """
    audit = {"previously_used_contigs": sorted(PREVIOUSLY_USED_CONTIGS), "cells": {}}
    for key, region in loaded.items():
        spec = region["spec"]
        labels = region["labels"]
        positions = region["positions"]
        entry = {
            "contig": region["contig"],
            "region": region["region"],
            "expected_region": [spec["start"], spec["stop"]],
            "loci_cached": int(labels.size),
            "loci_from_window_tiling": int(positions.size),
            "window_aligned": bool(labels.size % SEQ_LEN == 0),
            "non_finite_binomial": int((~np.isfinite(region["binomial_llr"])).sum()),
            "non_finite_pb": int((~np.isfinite(region["pb_llr"])).sum()),
            "position_min": int(positions.min()) if positions.size else None,
            "position_max": int(positions.max()) if positions.size else None,
            "positions_inside_region": bool(
                positions.size and positions.min() >= spec["start"]
                and positions.max() < spec["stop"]),
            "overlaps_previously_used_contig": region["contig"] in PREVIOUSLY_USED_CONTIGS,
        }
        entry["ok"] = bool(
            entry["contig"] == spec["contig"]
            and entry["region"] == entry["expected_region"]
            and entry["window_aligned"]
            and entry["loci_cached"] == entry["loci_from_window_tiling"]
            and entry["non_finite_binomial"] == 0 and entry["non_finite_pb"] == 0
            and entry["positions_inside_region"]
            and not entry["overlaps_previously_used_contig"])
        audit["cells"][key] = entry
    audit["all_ok"] = all(cell["ok"] for cell in audit["cells"].values())
    audit["n_cells"] = len(audit["cells"])
    return audit


# -- structure stratification -------------------------------------------------


def structure_stratified(region: dict) -> dict:
    """All three arms inside and outside each external difficulty annotation."""
    return {label: arms_on_mask(region, mask)
            for label, mask in structure_masks(region["spec"], region["positions"]).items()}


# -- wall clock ---------------------------------------------------------------


def measure_throughput(spec: dict, span_bp: int = 128_000) -> dict:
    """Binomial and repaired-PB throughput on this experiment's own data.

    Devlog 13 measured throughput on chr21 BAMs. Re-measuring per depth regime
    on the v14 BAMs keeps the speedup projection on the same data the accuracy
    numbers come from. Run only after every extraction has exited, so nothing
    competes for the CPU.
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


def projected_wallclock(region: dict, throughput: dict) -> dict:
    """Cascade vs PB-only runtime from measured rates, plus measured extraction."""
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
        "note": ("extraction timings were produced with four extractions running "
                 "concurrently and are upper bounds; the throughput probes used "
                 "for every speedup claim were measured in a quiet single process"),
    }


def accuracy_compute_curve(region: dict, throughput: dict | None) -> list:
    """The accuracy-vs-compute tradeoff, cutoff by cutoff. **Diagnostic only.**

    The production cutoff stays frozen at ``FROZEN_ROUTER_CUTOFF``. This curve
    shows what the frozen point buys and costs relative to its neighbours; no row
    of it may be promoted to a production threshold.
    """
    rows = cutoff_sweep(region)
    if throughput is None:
        return rows
    b_rate = throughput["binomial_loci_per_second"]
    pb_rate = throughput["pb_loci_per_second"]
    n = int(scoring_frame(region).sum())
    pb_only_seconds = n / pb_rate
    for row in rows:
        row["pb_compute_fraction"] = row["fraction_routed_to_pb"]
        row["projected_cascade_seconds"] = (n / b_rate
                                            + row["fraction_routed_to_pb"] * n / pb_rate)
        row["projected_pb_only_seconds"] = pb_only_seconds
        row["projected_speedup"] = pb_only_seconds / max(row["projected_cascade_seconds"], 1e-12)
    return rows


# -- verdict ------------------------------------------------------------------


def classify(entry: dict) -> str:
    """The devlog 13 §3.4 / devlog 14 §3.4 decision rule, applied verbatim.

    Order matters: UNDERPOWERED is checked first, so a cell that cannot support
    a claim is never rounded to PRESERVED.
    """
    delta = entry["delta_f1_vs_pb"]
    boot = entry.get("vs_pb_paired_bootstrap") or {}
    interval = boot.get("delta_f1_ci95")
    snp = int(entry.get("snp", 0))
    if not interval:
        return "UNDERPOWERED"
    low, high = float(interval[0]), float(interval[1])
    half_width = (high - low) / 2.0
    if snp < 100 or half_width > 0.01:
        return "UNDERPOWERED"
    contains_zero = low <= 0.0 <= high
    if contains_zero and abs(delta) < 0.001:
        return "PRESERVED"
    if not contains_zero and delta > 0:
        return "IMPROVED"
    return "DEGRADED"


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
    entry["by_structure"] = structure_stratified(region)
    entry["strata"] = {label: arms_on_mask(region, mask)
                       for label, mask in strata_masks(region).items()}
    entry["controls"] = controls(region)
    entry["accuracy_compute_curve"] = accuracy_compute_curve(region, throughput)
    if throughput is not None:
        entry["projected_wallclock"] = projected_wallclock(region, throughput)
    entry["verdict"] = classify(entry)
    return entry


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
            spec = specs[f"chr20_neutral_{tag}"]
            logger.info("throughput probe: %s", tag)
            throughput[tag] = measure_throughput(spec)

    report = {
        "frozen": {"binomial_threshold": FROZEN_BINOMIAL_THRESHOLD,
                   "pb_threshold": FROZEN_PB_THRESHOLD,
                   "router_cutoff": FROZEN_ROUTER_CUTOFF},
        "preregistration": "History/14_DEVLOG.md sections 1-8",
        "independence_audit": {"all_ok": audit["all_ok"], "path": str(AUDIT)},
        "compute": {"throughput": throughput,
                    "throughput_note": ("measured in this process on the v14 BAMs, after "
                                        "every extraction had exited")},
        "cells": {}, "pooled": {},
    }

    disagreements: dict[str, dict] = {}
    for key, region in loaded.items():
        logger.info("scoring %s", key)
        report["cells"][key] = score_cell(region, throughput.get(region["spec"]["depth_tag"]))
        disagreements[key] = disagreement_records(region)

    report["pooled"]["all"] = _pool(list(loaded.values()))
    for name in BASE_REGIONS:
        report["pooled"][f"region_{name}"] = _pool(
            [r for k, r in loaded.items() if r["spec"]["base"] == name])
    for tag in DEPTH_TAGS:
        report["pooled"][f"depth_{tag}"] = _pool(
            [r for k, r in loaded.items() if r["spec"]["depth_tag"] == tag])
    # Native-depth cells are the only ones that are not downsamples of a
    # downsample; pooling them separately keeps the independent-library claim
    # honest (devlog 14 §5.4).
    report["pooled"]["ordinary_structure_only"] = _pool(
        [r for k, r in loaded.items() if r["spec"]["base"] != "chr1_segdup"])

    report["failure_boundary"] = pooled_failure_boundary(disagreements)

    RESULTS.write_text(json.dumps(report, indent=2))
    DISAGREEMENTS.write_text(json.dumps(disagreements, indent=2))
    logger.info("wrote %s and %s", RESULTS, DISAGREEMENTS)


if __name__ == "__main__":
    main()
