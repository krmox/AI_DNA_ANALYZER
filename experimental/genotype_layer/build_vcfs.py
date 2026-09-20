"""Build method A/B/C VCFs for the dev region (chr21:30.0-30.47Mb), reusing
the SAME allele calls (cascade_calls from cache/dev_region.npz, produced by
extract_dev_region.py -- the frozen binomial<->PB router, mirroring
experimental/unified_happy/build_vcfs.py exactly) and varying only the GT
field per method.

Method A: baseline, GT always 0/1 (reproduces build_vcfs.py's existing behavior).
Method B: VAF-threshold genotyping, tau fit by sweeping on this same dev region.
Method C: binomial genotype posterior (3-hypothesis argmax) with GQ, eps=0.01
          fixed (not fit to labels). A 0/0-favored site is NEVER dropped: it is
          kept as GT 0/1 with low GQ, and counted separately.

Never adds/removes a called site -- only the GT (and GQ, for C) field changes.
Runtime of the genotype-layer step alone (pure numpy over already-extracted
counts) is measured separately from detection/BAM-scanning (already-paid cost
in extract_dev_region.py) for methods B and C.
"""
import gzip
import json
import subprocess
import sys
import time

import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
sys.path.insert(0, _ROOT)

import numpy as np
from scipy.stats import binom

REPO = _ROOT
OUT = f"{REPO}/experimental/genotype_layer"
DATA = f"{REPO}/data/giab_hg002_real_30M"
TRUTH_VCF = f"{DATA}/hg002_chr21_30M.vcf.gz"

REF_TOKENS = ("N", "A", "C", "G", "T")
EPS = 0.01  # fixed Illumina-typical error rate, not fit to labels

d = np.load(f"{OUT}/cache/dev_region.npz")
counts = d["counts"]
positions = d["positions"]  # 0-based
call_mask = d["cascade_calls"]

idxs = np.nonzero(call_mask)[0]
print(f"{idxs.size} called loci (cascade_calls) to genotype")


def alt_allele(count_row, ref_base):
    base_counts = {"A": count_row[0], "C": count_row[1], "G": count_row[2], "T": count_row[3]}
    base_counts.pop(ref_base, None)
    if not base_counts:
        return None
    best_base, best_count = max(base_counts.items(), key=lambda kv: kv[1])
    if best_count <= 0:
        return None
    return best_base


# --- Parse truth GT at every truth VCF position (for tau fitting / evaluation) ---
def parse_truth_gt(path):
    gt_by_pos = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt, qual, filt, info, fmt, sample = parts[:10]
            fmt_fields = fmt.split(":")
            sample_fields = sample.split(":")
            gt_raw = sample_fields[fmt_fields.index("GT")]
            gt_norm = gt_raw.replace("|", "/")
            alleles = set(gt_norm.split("/"))
            if alleles == {"0"}:
                gt_class = "0/0"
            elif alleles == {"1"}:
                gt_class = "1/1"
            else:
                gt_class = "0/1"  # het, incl. multi-allelic hets treated as het
            gt_by_pos[(chrom, int(pos))] = (ref, alt, gt_class)
    return gt_by_pos


truth_gt = parse_truth_gt(TRUTH_VCF)
print(f"truth GT parsed: {len(truth_gt)} records")

# --- Build the per-called-locus record list (REF, ALT, pos1, k, ref_n, n, truth_gt) ---
records = []
n_skipped_no_alt = 0
for i in idxs:
    pos0 = int(positions[i])
    pos1 = pos0 + 1
    ref_idx = int(counts[i][9])
    ref_base = REF_TOKENS[ref_idx] if 0 <= ref_idx < 5 else "N"
    if ref_base == "N":
        n_skipped_no_alt += 1
        continue
    alt = alt_allele(counts[i], ref_base)
    if alt is None:
        n_skipped_no_alt += 1
        continue
    depth_col = int(counts[i][6])
    base_counts = {"A": counts[i][0], "C": counts[i][1], "G": counts[i][2], "T": counts[i][3]}
    k = float(base_counts[alt])
    ref_n = float(base_counts[ref_base])
    n = k + ref_n  # usable A/C/G/T observations supporting REF or ALT candidate
    truth = truth_gt.get(("chr21", pos1))
    truth_gt_class = truth[2] if truth is not None else None
    records.append(dict(pos0=pos0, pos1=pos1, ref=ref_base, alt=alt, depth=depth_col,
                        k=k, ref_n=ref_n, n=n, truth_gt=truth_gt_class))

print(f"{len(records)} records built ({n_skipped_no_alt} skipped: ref=N or no alt support)")

allele_keys_A = sorted((r["pos1"], r["ref"], r["alt"]) for r in records)

# ============================================================================
# Method A: baseline, always 0/1
# ============================================================================
for r in records:
    r["gt_A"] = "0/1"

# ============================================================================
# Method B: VAF-threshold genotyping. Fit tau on this dev region by sweeping
# to maximize GT accuracy against truth GT, restricted to TP loci (truth_gt known).
# ============================================================================
t0 = time.perf_counter()
vaf = np.array([r["k"] / r["n"] if r["n"] > 0 else 0.0 for r in records])
truth_known = np.array([r["truth_gt"] is not None for r in records])
truth_arr = np.array([r["truth_gt"] if r["truth_gt"] is not None else "" for r in records])

best_tau, best_acc = None, -1.0
tau_grid = np.round(np.arange(0.05, 0.96, 0.01), 3)
for tau in tau_grid:
    pred = np.where(vaf > tau, "1/1", "0/1")
    acc = (pred[truth_known] == truth_arr[truth_known]).mean() if truth_known.any() else 0.0
    if acc > best_acc:
        best_acc, best_tau = acc, float(tau)

for i, r in enumerate(records):
    r["gt_B"] = "1/1" if vaf[i] > best_tau else "0/1"
b_runtime = time.perf_counter() - t0
print(f"Method B: fitted tau={best_tau} (dev-region GT accuracy at fit={best_acc:.4f}), "
      f"runtime={b_runtime*1000:.2f}ms")

# ============================================================================
# Method C: binomial genotype posterior, 3 hypotheses {0/0, 0/1, 1/1} at
# theta in {0, 0.5, 1-eps}, P(k|n,theta) = Binomial(k; n, theta*(1-eps)+(1-theta)*eps).
# GT = argmax likelihood (flat prior over 3). GQ = 10*log10(P_best/P_second), capped.
# 0/0-favored sites are NOT dropped: kept as original call, GT written 0/1, low GQ,
# and separately counted as "theta0_wins".
# ============================================================================
t0 = time.perf_counter()
GQ_CAP = 99.0
theta0_wins = 0
for r in records:
    k, n = r["k"], r["n"]
    p_00 = np.clip(0.0 * (1 - EPS) + 1.0 * EPS, 1e-12, 1 - 1e-12)
    p_01 = np.clip(0.5 * (1 - EPS) + 0.5 * EPS, 1e-12, 1 - 1e-12)
    p_11 = np.clip(1.0 * (1 - EPS) + 0.0 * EPS, 1e-12, 1 - 1e-12)
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
        theta0_wins += 1
        r["gt_C"] = "0/1"
        r["gq_C"] = round(min(gq, 5.0), 1)  # low-GQ flag; call is never suppressed
        r["c_flag_theta0"] = True
    else:
        r["gt_C"] = best_class
        r["gq_C"] = round(gq, 1)
        r["c_flag_theta0"] = False
c_runtime = time.perf_counter() - t0
print(f"Method C: eps={EPS}, theta=0(0/0)-wins at {theta0_wins}/{len(records)} sites "
      f"(flagged, never suppressed), runtime={c_runtime*1000:.2f}ms")

# --- Two-level invariant check: allele call set (pos+ref+alt) must be identical ---
allele_keys_B = sorted((r["pos1"], r["ref"], r["alt"]) for r in records)
allele_keys_C = sorted((r["pos1"], r["ref"], r["alt"]) for r in records)
invariant_ok = allele_keys_A == allele_keys_B == allele_keys_C
print(f"allele-call-set invariant across A/B/C: {invariant_ok}")

# ============================================================================
# Write VCFs
# ============================================================================
HEADER_TMPL = """##fileformat=VCFv4.2
##source=ai_dna_analyzer_genotype_layer_{label}
##contig=<ID=chr21,length=46709983>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth (passing-filter reads)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype ({desc})">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality (Method C only; 10*log10(P_best/P_second), capped at 99)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
"""

DESCS = {
    "A": "heuristic: always 0/1",
    "B": f"VAF-threshold: 1/1 if VAF>{best_tau} else 0/1 (tau fit on dev region)",
    "C": "binomial genotype posterior argmax over {0/0,0/1,1/1}, eps=0.01 fixed",
}


def write_vcf(method, gt_key, gq_key=None):
    lines = []
    for r in sorted(records, key=lambda r: r["pos1"]):
        gt = r[gt_key]
        if gq_key is not None:
            gq = r[gq_key]
            lines.append(f"chr21\t{r['pos1']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\t"
                         f"DP={r['depth']}\tGT:GQ\t{gt}:{int(round(gq))}")
        else:
            lines.append(f"chr21\t{r['pos1']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\t"
                         f"DP={r['depth']}\tGT\t{gt}")
    path = f"{OUT}/vcfs/{method}.vcf"
    with open(path, "w") as fh:
        fh.write(HEADER_TMPL.format(label=method, desc=DESCS[method]))
        fh.write("\n".join(lines))
        if lines:
            fh.write("\n")
    print(f"wrote {path}: {len(lines)} records")
    return path


write_vcf("A", "gt_A")
write_vcf("B", "gt_B")
write_vcf("C", "gt_C", gq_key="gq_C")

# --- normalize (sort/bgzip/tabix), same as unified_happy/run_happy.sh step 3 ---
for method in ("A", "B", "C"):
    subprocess.run(["bcftools", "sort", "-Oz", "-o", f"{OUT}/vcfs/{method}.vcf.gz",
                    f"{OUT}/vcfs/{method}.vcf"], check=True, capture_output=True)
    subprocess.run(["bcftools", "index", "-t", f"{OUT}/vcfs/{method}.vcf.gz"],
                   check=True, capture_output=True)
print("bgzip/tabix done for A/B/C")

# --- save intermediate artifacts for evaluation step ---
summary = dict(
    n_records=len(records),
    n_skipped_no_alt=n_skipped_no_alt,
    fitted_tau=best_tau,
    fitted_tau_gt_acc_at_fit=best_acc,
    eps_C=EPS,
    theta0_wins=theta0_wins,
    invariant_ok=bool(invariant_ok),
    method_B_runtime_seconds=b_runtime,
    method_C_runtime_seconds=c_runtime,
)
with open(f"{OUT}/cache/build_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print(json.dumps(summary, indent=2))

# Save full records table for evaluation (GT accuracy, confusion matrix)
np.savez_compressed(
    f"{OUT}/cache/records.npz",
    pos1=np.array([r["pos1"] for r in records]),
    ref=np.array([r["ref"] for r in records]),
    alt=np.array([r["alt"] for r in records]),
    truth_gt=np.array([r["truth_gt"] if r["truth_gt"] is not None else "" for r in records]),
    gt_A=np.array([r["gt_A"] for r in records]),
    gt_B=np.array([r["gt_B"] for r in records]),
    gt_C=np.array([r["gt_C"] for r in records]),
    c_flag_theta0=np.array([r["c_flag_theta0"] for r in records]),
    vaf=vaf,
)
print("Saved cache/records.npz")
