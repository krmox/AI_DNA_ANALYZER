"""bench_v19 §11 — real end-to-end wall-clock timing, PB-only vs frozen cascade.

Not a frozen file: this is new orchestration glue for Validation Round 2,
pre-registered in VALIDATION_ROUND_2.md §11. It calls only frozen functions
(``pileup_counts.load_counts``, ``read_level_pileup.load_reads``,
``binomial_baseline.BinomialVariantCaller.score_counts``,
``extract_bench_v12.poisson_binomial_for_chunk``, ``cascade.route_stage1``)
and does not modify any of them.

Both arms cover the complete path BAM -> pileup/read extraction -> caller ->
final calls; nothing is excluded from the timer. The cascade arm still loads
the full read tensor for the span (the router's own decision depends on the
binomial pass, which itself needs pileup counts first, so BAM I/O cannot be
skipped ahead of routing) -- only the expensive Poisson-binomial computation
is restricted to the routed subset. This is disclosed as the realistic
scenario, not an optimistic one: if I/O dominates, it will show up as both
arms taking similar wall-clock despite very different PB-locus counts.

Sub-window: chr2:210,000,000-210,500,000 (500 kb), a compute-budget
concession fixed in the pre-registration, not used for any accuracy metric.
"""

from __future__ import annotations

import json
import resource
import time
from pathlib import Path

import numpy as np

from binomial_baseline import BinomialVariantCaller
from cascade import route_stage1
from extract_bench_v12 import poisson_binomial_for_chunk
from pileup_counts import load_counts
from read_level_pileup import load_reads

FASTA = "data/reference/chr2_full.fa"
BAM = "data/giab_hg004_v19/chr2_timing_full.bam"
VCF = "data/giab_hg004_v19/chr2_timing.vcf.gz"
BED = "data/giab_hg004_v19/chr2_timing_highconf.bed"
CONTIG = "chr2"
SPAN = (210_000_000, 210_500_000)
REPETITIONS = 3
OUT = Path("results/bench_v19/timing.json")


def peak_rss_kb() -> int:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


def run_pb_only() -> dict:
    rss0 = peak_rss_kb()
    t0 = time.perf_counter()
    counts, _ = load_counts(FASTA, BAM, VCF, BED, SPAN, contig=CONTIG)
    reads, _ = load_reads(FASTA, BAM, VCF, BED, SPAN, contig=CONTIG)
    pb_llr = poisson_binomial_for_chunk(reads, counts)
    wall = time.perf_counter() - t0
    return {"wall_seconds": wall, "loci": int(counts.shape[0]),
            "pb_loci": int(counts.shape[0]), "peak_rss_kb": peak_rss_kb() - rss0}


def run_cascade() -> dict:
    rss0 = peak_rss_kb()
    t0 = time.perf_counter()
    counts, _ = load_counts(FASTA, BAM, VCF, BED, SPAN, contig=CONTIG)
    reads, _ = load_reads(FASTA, BAM, VCF, BED, SPAN, contig=CONTIG)
    caller = BinomialVariantCaller(error_rate=0.01)
    binomial_score = caller.score_counts(counts[:, 0:4], counts[:, 9].astype(int))
    routed = route_stage1(binomial_score.llr)
    if routed.any():
        poisson_binomial_for_chunk(reads[routed], counts[routed])
    wall = time.perf_counter() - t0
    return {"wall_seconds": wall, "loci": int(counts.shape[0]),
            "pb_loci": int(routed.sum()), "peak_rss_kb": peak_rss_kb() - rss0}


def main() -> None:
    results = {"pb_only": [], "cascade": []}
    for rep in range(REPETITIONS):
        r = run_pb_only()
        print(f"pb_only rep {rep}: {r}")
        results["pb_only"].append(r)
    for rep in range(REPETITIONS):
        r = run_cascade()
        print(f"cascade rep {rep}: {r}")
        results["cascade"].append(r)

    for arm in ("pb_only", "cascade"):
        times = [r["wall_seconds"] for r in results[arm]]
        results[f"{arm}_summary"] = {
            "median_s": float(np.median(times)), "min_s": float(np.min(times)),
            "max_s": float(np.max(times)), "std_s": float(np.std(times)),
        }
    median_pb = results["pb_only_summary"]["median_s"]
    median_cascade = results["cascade_summary"]["median_s"]
    results["measured_end_to_end_speedup"] = median_pb / max(median_cascade, 1e-9)
    results["span"] = list(SPAN)
    results["contig"] = CONTIG
    results["note"] = ("Real wall-clock, full BAM->pileup->caller path, both arms. "
                       "NOT the caller-only projected speedup reported elsewhere "
                       "in this project -- see VALIDATION_ROUND_2.md section 11.")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
