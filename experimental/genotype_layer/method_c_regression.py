"""Stage-2 revalidation: build an OPTIMIZED (vectorized) Method C implementation,
prove it is numerically equivalent to the REFERENCE (per-site scipy.stats.binom
loop, verbatim from build_vcfs_validation_C.py / build_vcfs.py) Method C, on the
EXACT SAME inputs used for the Stage-1 unified hap.py head-to-head benchmark
(experimental/unified_happy): recovered_full.npz (counts/positions) + the frozen
pb_calls/cascade_calls masks from experimental/head_to_head/ai_cascade_per_locus.npz.

This does NOT reuse experimental/genotype_layer/cache/validation_region.npz --
that file's pb_calls/cascade_calls masks are from a DIFFERENT extraction pass
(16121/16124 calls, non-identical positions/counts) and must not be silently
substituted for the Stage-1 unified_happy scope (16094/16097 calls). This script
only touches the unified_happy-scope arrays so the resulting VCFs are on IDENTICAL
input to Stage 1.

Two Method C implementations, mathematically identical formula:
  REFERENCE: python per-site loop, scipy.stats.binom.logpmf called once per
             hypothesis per site (verbatim from build_vcfs_validation_C.py).
  OPTIMIZED: fully vectorized numpy; closed-form binomial log-pmf via
             scipy.special.gammaln, computed once over all sites/hypotheses as
             array ops (no python-level per-site loop, no scipy.stats overhead).

Writes:
  experimental/genotype_layer/cache/method_c_regression_report.json
  experimental/unified_happy/ai_pb_only_C.vcf(.gz/.tbi)
  experimental/unified_happy/ai_cascade_C.vcf(.gz/.tbi)
"""
import json
import os
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

import numpy as np
from scipy.stats import binom
from scipy.special import gammaln

UNIFIED = f"{REPO}/experimental/unified_happy"
GTL = f"{REPO}/experimental/genotype_layer"

REF_TOKENS = ("N", "A", "C", "G", "T")
EPS = 0.01
GQ_CAP = 99.0

recovered = np.load(f"{UNIFIED}/recovered_full.npz")
labels_recovered_full = recovered["labels"]
positions_recovered_full = recovered["positions"]
counts_recovered_full = recovered["counts"]

# Reproduce experimental/unified_happy/build_vcfs.py's frame filter exactly:
# SNP-eval frame = labels in {0, 1}. This is what aligns recovered_full.npz's
# per-locus arrays with the frozen pb_calls/cascade_calls masks below.
frame = (labels_recovered_full == 0) | (labels_recovered_full == 1)
positions_full = positions_recovered_full[frame]
counts_full = counts_recovered_full[frame]

stored = np.load(f"{REPO}/experimental/head_to_head/ai_cascade_per_locus.npz")
stored_labels = stored["labels"]
pb_calls = stored["pb_calls"]
cascade_calls = stored["cascade_calls"]

assert labels_recovered_full[frame].size == stored_labels.size, (
    "size mismatch between recovered_full.npz (frame-filtered) and "
    "ai_cascade_per_locus.npz -- cannot safely zip positions with masks")
n_label_mismatch = int(np.count_nonzero(labels_recovered_full[frame] != stored_labels))
match_frac = 1.0 - n_label_mismatch / stored_labels.size
print(f"label agreement (frame-filtered recovered_full vs stored): {match_frac:.6f} "
      f"({n_label_mismatch} mismatches)")
assert match_frac >= 0.999, (
    "label arrays disagree beyond tolerance -- ordering assumption invalid, "
    "same check as build_vcfs.py")
assert positions_full.shape[0] == pb_calls.shape[0] == cascade_calls.shape[0], (
    "size mismatch between recovered_full.npz and ai_cascade_per_locus.npz -- "
    "cannot safely zip positions with masks")


def alt_allele(count_row, ref_base):
    base_counts = {"A": count_row[0], "C": count_row[1], "G": count_row[2], "T": count_row[3]}
    base_counts.pop(ref_base, None)
    if not base_counts:
        return None
    best_base, best_count = max(base_counts.items(), key=lambda kv: kv[1])
    if best_count <= 0:
        return None
    return best_base


def extract_records(call_mask, label):
    idxs = np.nonzero(call_mask)[0]
    records = []
    n_skipped = 0
    for i in idxs:
        pos0 = int(positions_full[i])
        ref_idx = int(counts_full[i][9])
        ref_base = REF_TOKENS[ref_idx] if 0 <= ref_idx < 5 else "N"
        if ref_base == "N":
            n_skipped += 1
            continue
        alt = alt_allele(counts_full[i], ref_base)
        if alt is None:
            n_skipped += 1
            continue
        base_counts = {"A": counts_full[i][0], "C": counts_full[i][1],
                        "G": counts_full[i][2], "T": counts_full[i][3]}
        k = float(base_counts[alt])
        ref_n = float(base_counts[ref_base])
        n = k + ref_n
        depth = int(counts_full[i][6])
        records.append(dict(pos0=pos0, pos1=pos0 + 1, ref=ref_base, alt=alt, depth=depth,
                             k=k, n=n))
    print(f"[{label}] {len(records)} records ({n_skipped} skipped: ref=N or no alt support)")
    return records


def reference_method_c(records):
    """Verbatim per-site scipy.stats.binom loop (build_vcfs_validation_C.py logic)."""
    t0 = time.perf_counter()
    p_00 = np.clip(0.0 * (1 - EPS) + 1.0 * EPS, 1e-12, 1 - 1e-12)
    p_01 = np.clip(0.5 * (1 - EPS) + 0.5 * EPS, 1e-12, 1 - 1e-12)
    p_11 = np.clip(1.0 * (1 - EPS) + 0.0 * EPS, 1e-12, 1 - 1e-12)
    gts, gqs, flags = [], [], []
    for r in records:
        k, n = r["k"], r["n"]
        logp = {
            "0/0": binom.logpmf(k, n, p_00) if n > 0 else np.log(1.0),
            "0/1": binom.logpmf(k, n, p_01) if n > 0 else np.log(1.0),
            "1/1": binom.logpmf(k, n, p_11) if n > 0 else np.log(1.0),
        }
        order = sorted(logp.items(), key=lambda kv: kv[1], reverse=True)
        best_class, best_logp = order[0]
        second_logp = order[1][1]
        gq = min(10.0 * (best_logp - second_logp) / np.log(10.0), GQ_CAP)
        gq = max(gq, 0.0)
        if best_class == "0/0":
            gts.append("0/1")
            gqs.append(round(min(gq, 5.0), 1))
            flags.append(True)
        else:
            gts.append(best_class)
            gqs.append(round(gq, 1))
            flags.append(False)
    runtime = time.perf_counter() - t0
    return np.array(gts), np.array(gqs), np.array(flags), runtime


def optimized_method_c(records):
    """Fully vectorized: closed-form binomial log-pmf via gammaln, array ops only."""
    t0 = time.perf_counter()
    k = np.array([r["k"] for r in records], dtype=np.float64)
    n = np.array([r["n"] for r in records], dtype=np.float64)
    p = np.array([
        np.clip(0.0 * (1 - EPS) + 1.0 * EPS, 1e-12, 1 - 1e-12),
        np.clip(0.5 * (1 - EPS) + 0.5 * EPS, 1e-12, 1 - 1e-12),
        np.clip(1.0 * (1 - EPS) + 0.0 * EPS, 1e-12, 1 - 1e-12),
    ])  # order: 0/0, 0/1, 1/1

    # log C(n,k) + k*log(p) + (n-k)*log(1-p), vectorized over (site, hypothesis)
    log_binom_coeff = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    logp = (log_binom_coeff[:, None]
            + k[:, None] * np.log(p)[None, :]
            + (n - k)[:, None] * np.log(1 - p)[None, :])
    # n==0 sites: reference used log(1.0)=0.0 for all three hypotheses
    zero_n = (n == 0)
    logp[zero_n, :] = 0.0

    order = np.argsort(-logp, axis=1)  # descending
    best_idx = order[:, 0]
    second_idx = order[:, 1]
    best_logp = logp[np.arange(len(records)), best_idx]
    second_logp = logp[np.arange(len(records)), second_idx]

    gq_raw = 10.0 * (best_logp - second_logp) / np.log(10.0)
    gq_raw = np.minimum(gq_raw, GQ_CAP)
    gq_raw = np.maximum(gq_raw, 0.0)

    classes = np.array(["0/0", "0/1", "1/1"])
    best_class = classes[best_idx]

    is_theta0 = best_class == "0/0"
    gt_out = np.where(is_theta0, "0/1", best_class)
    gq_capped = np.where(is_theta0, np.minimum(gq_raw, 5.0), gq_raw)
    gq_out = np.round(gq_capped, 1)
    flags = is_theta0

    runtime = time.perf_counter() - t0
    return gt_out, gq_out, flags, runtime


HEADER_TMPL = """##fileformat=VCFv4.2
##source=ai_dna_analyzer_unified_happy_{label}_C_optimized
##contig=<ID=chr21,length=46709983>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth (passing-filter reads)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (binomial genotype posterior argmax over {{0/0,0/1,1/1}}, eps=0.01 fixed, optimized vectorized implementation)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality (10*log10(P_best/P_second), capped at 99)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
"""


def write_vcf(records, gt, gq, out_path, label):
    order = np.argsort([r["pos1"] for r in records])
    lines = []
    for i in order:
        r = records[i]
        lines.append(f"chr21\t{r['pos1']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\t"
                      f"DP={r['depth']}\tGT:GQ\t{gt[i]}:{int(round(float(gq[i])))}")
    with open(out_path, "w") as fh:
        fh.write(HEADER_TMPL.format(label=label))
        fh.write("\n".join(lines))
        if lines:
            fh.write("\n")
    print(f"wrote {out_path}: {len(lines)} records")


report = {"scope": "unified_happy (Stage 1 identical inputs)", "datasets": {}}

for mask, label, out_name in [
    (pb_calls, "ai_pb_only", "ai_pb_only_C.vcf"),
    (cascade_calls, "ai_cascade", "ai_cascade_C.vcf"),
]:
    records = extract_records(mask, label)

    ref_gt, ref_gq, ref_flag, ref_rt = reference_method_c(records)
    opt_gt, opt_gq, opt_flag, opt_rt = optimized_method_c(records)

    gt_equal = bool(np.array_equal(ref_gt, opt_gt))
    gq_equal = bool(np.array_equal(ref_gq, opt_gq))
    flag_equal = bool(np.array_equal(ref_flag, opt_flag))
    n_gt_mismatch = int(np.count_nonzero(ref_gt != opt_gt))
    n_gq_mismatch = int(np.count_nonzero(ref_gq != opt_gq))
    max_gq_abs_diff = float(np.max(np.abs(ref_gq.astype(float) - opt_gq.astype(float)))) if len(records) else 0.0

    print(f"[{label}] GT equal={gt_equal} ({n_gt_mismatch} mismatches), "
          f"GQ equal={gq_equal} ({n_gq_mismatch} mismatches, max|diff|={max_gq_abs_diff}), "
          f"flags equal={flag_equal}")
    print(f"[{label}] reference runtime={ref_rt*1000:.2f}ms "
          f"({ref_rt/len(records)*1e6:.2f} us/site), "
          f"optimized runtime={opt_rt*1000:.2f}ms "
          f"({opt_rt/len(records)*1e6:.2f} us/site), "
          f"speedup={ref_rt/opt_rt:.2f}x")

    report["datasets"][label] = dict(
        n_records=len(records),
        gt_equal=gt_equal, gq_equal=gq_equal, flag_equal=flag_equal,
        n_gt_mismatch=n_gt_mismatch, n_gq_mismatch=n_gq_mismatch,
        max_gq_abs_diff=max_gq_abs_diff,
        reference_runtime_ms=ref_rt * 1000, optimized_runtime_ms=opt_rt * 1000,
        speedup=ref_rt / opt_rt if opt_rt > 0 else None,
    )

    if not (gt_equal and flag_equal):
        print(f"FATAL: optimized Method C changed at least one genotype call vs "
              f"reference on [{label}]. STOPPING per protocol -- no VCF written, "
              f"no hap.py run.", file=sys.stderr)
        with open(f"{GTL}/cache/method_c_regression_report.json", "w") as fh:
            json.dump(report, fh, indent=2)
        sys.exit(1)

    # Use OPTIMIZED output (proven identical) for the actual benchmark VCF.
    write_vcf(records, opt_gt, opt_gq, f"{UNIFIED}/{out_name}", label)

with open(f"{GTL}/cache/method_c_regression_report.json", "w") as fh:
    json.dump(report, fh, indent=2)

print("\nAll genotype arrays identical between reference and optimized Method C.")
print(json.dumps(report, indent=2))
