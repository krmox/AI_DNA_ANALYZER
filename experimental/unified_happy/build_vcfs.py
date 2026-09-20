"""Step 2/3: build ai_pb_only.vcf.gz and ai_cascade.vcf.gz from recovered
true positions + the existing frozen per-locus calls in
experimental/head_to_head/ai_cascade_per_locus.npz.

Reuses (imports, does not modify) frozen decision arrays already computed by
experimental/head_to_head/analyze_cascade.py from cascade.py / robustness_benchmark.py
logic: `pb_calls` (PB-only, threshold FROZEN_PB_THRESHOLD) and `cascade_calls`
(route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF) -> pb_calls else binomial_calls),
exactly as cascade.py's run_cascade / robustness_benchmark.py's route_mask define
(no mamba/neural stage2 was scored for this study -- `--features none` was used
for the underlying extraction, so the cascade here is the frozen binomial<->PB
two-stage router, not the 3-stage mamba cascade; this matches
experimental/head_to_head/COMMANDS.md section 5/7 exactly).
"""
import gzip
import subprocess
import sys

import numpy as np

import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
REPO = _ROOT
OUT = f"{REPO}/experimental/unified_happy"

recovered = np.load(f"{OUT}/recovered_full.npz")
labels_full = recovered["labels"]
positions_full = recovered["positions"]  # 0-based
counts_full = recovered["counts"]

frame = (labels_full == 0) | (labels_full == 1)
positions = positions_full[frame]
counts = counts_full[frame]
labels = labels_full[frame]

print(f"recovered_full: {labels_full.size} loci total, {frame.sum()} in SNP-eval frame")

stored = np.load(f"{REPO}/experimental/head_to_head/ai_cascade_per_locus.npz")
stored_labels = stored["labels"]
stored_snp = stored["snp"]
pb_calls = stored["pb_calls"]
cascade_calls = stored["cascade_calls"]

print(f"stored (ai_cascade_per_locus.npz): {stored_labels.size} loci in frame")

# --- Consistency check: recomputed frame-filtered labels must match stored labels
# exactly (same size AND same values), which proves the chunk/window ordering of
# this recomputation matches extract_bench_v12.py's original pass byte-for-byte,
# so `positions` (recovered here) can be safely zipped with `pb_calls` /
# `cascade_calls` (from the stored file) index-for-index.
if labels.size != stored_labels.size:
    print(f"FATAL: size mismatch {labels.size} vs {stored_labels.size} -- "
          f"cannot safely align positions with stored calls.", file=sys.stderr)
    sys.exit(1)

n_mismatch = int(np.count_nonzero(labels != stored_labels))
match_frac = 1.0 - n_mismatch / labels.size
print(f"label agreement: {match_frac:.6f} ({n_mismatch} mismatches out of {labels.size})")
if match_frac < 0.999:
    print("FATAL: label arrays disagree beyond tolerance -- ordering assumption invalid.",
          file=sys.stderr)
    sys.exit(1)

snp_recomputed = labels == 1
snp_match = np.array_equal(snp_recomputed, stored_snp)
print(f"snp array exact match: {snp_match}")

REF_TOKENS = ("N", "A", "C", "G", "T")


def alt_allele(count_row: np.ndarray, ref_base: str) -> str | None:
    """Majority non-reference base actually observed in the count matrix."""
    base_counts = {"A": count_row[0], "C": count_row[1], "G": count_row[2], "T": count_row[3]}
    base_counts.pop(ref_base, None)
    if not base_counts:
        return None
    best_base, best_count = max(base_counts.items(), key=lambda kv: kv[1])
    if best_count <= 0:
        return None
    return best_base


def build_vcf(call_mask: np.ndarray, out_path: str, label: str) -> int:
    header = f"""##fileformat=VCFv4.2
##source=ai_dna_analyzer_unified_happy_{label}
##contig=<ID=chr21,length=46709983>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth (passing-filter reads)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (heuristic: always 0/1 -- see README limitation)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
"""
    n_written = 0
    n_skipped_no_alt = 0
    lines = []
    idxs = np.nonzero(call_mask)[0]
    for i in idxs:
        pos0 = int(positions[i])
        ref_idx = int(counts[i][9])
        ref_base = REF_TOKENS[ref_idx] if 0 <= ref_idx < 5 else "N"
        if ref_base == "N":
            n_skipped_no_alt += 1
            continue
        alt = alt_allele(counts[i], ref_base)
        if alt is None:
            n_skipped_no_alt += 1
            continue
        depth = int(counts[i][6])
        pos1 = pos0 + 1  # VCF is 1-based
        lines.append(f"chr21\t{pos1}\t.\t{ref_base}\t{alt}\t.\tPASS\tDP={depth}\tGT\t0/1")
        n_written += 1

    with open(out_path, "w") as fh:
        fh.write(header)
        fh.write("\n".join(lines))
        if lines:
            fh.write("\n")

    print(f"{label}: {n_written} records written ({n_skipped_no_alt} calls skipped: "
          f"ref=N or no alt-supporting reads in count matrix)")
    return n_written


n_pb = build_vcf(pb_calls, f"{OUT}/ai_pb_only.vcf", "ai_pb_only")
n_casc = build_vcf(cascade_calls, f"{OUT}/ai_cascade.vcf", "ai_cascade")

print("Done. Run sort/bgzip/tabix next.")
