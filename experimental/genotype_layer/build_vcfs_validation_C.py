"""Build method C VCF for the TRUE hold-out validation region (chr21:32-44Mb),
reusing the SAME allele calls (cascade_calls from cache/validation_region.npz --
identical 16,124 calls used for A/B) and applying EXACTLY the Method C binomial
genotype posterior logic from build_vcfs.py (dev region), verbatim: eps=0.01
fixed, theta in {0.0, 0.5, 1-eps}, GT=argmax, GQ=10*log10(P_best/P_second) capped
at 99, 0/0-favored sites kept (never dropped) as GT 0/1 with low GQ and flagged.

Does NOT touch validation_A.vcf / validation_B.vcf or tau=0.69.
Writes experimental/genotype_layer/vcfs/validation_C.vcf(.gz/.tbi) -- new file.
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
DATA = f"{REPO}/data/giab_hg002_chr21_12Mb"
TRUTH_VCF = f"{DATA}/hg002_chr21_32_44M.vcf.gz"

REF_TOKENS = ("N", "A", "C", "G", "T")
EPS = 0.01  # fixed, not fit -- same as dev-region Method C

d = np.load(f"{OUT}/cache/validation_region.npz")
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
                gt_class = "0/1"
            gt_by_pos[(chrom, int(pos))] = (ref, alt, gt_class)
    return gt_by_pos


truth_gt = parse_truth_gt(TRUTH_VCF)
print(f"truth GT parsed: {len(truth_gt)} records")

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
    n = k + ref_n
    truth = truth_gt.get(("chr21", pos1))
    truth_gt_class = truth[2] if truth is not None else None
    records.append(dict(pos0=pos0, pos1=pos1, ref=ref_base, alt=alt, depth=depth_col,
                        k=k, ref_n=ref_n, n=n, truth_gt=truth_gt_class))

print(f"{len(records)} records built ({n_skipped_no_alt} skipped: ref=N or no alt support)")

allele_keys_A = sorted((r["pos1"], r["ref"], r["alt"]) for r in records)

# ============================================================================
# Method C: binomial genotype posterior -- EXACT SAME LOGIC as build_vcfs.py
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
        r["gq_C"] = round(min(gq, 5.0), 1)
        r["c_flag_theta0"] = True
    else:
        r["gt_C"] = best_class
        r["gq_C"] = round(gq, 1)
        r["c_flag_theta0"] = False
c_runtime = time.perf_counter() - t0
print(f"Method C: eps={EPS}, theta=0(0/0)-wins at {theta0_wins}/{len(records)} sites "
      f"(flagged, never suppressed), runtime={c_runtime*1000:.2f}ms")

vaf = np.array([r["k"] / r["n"] if r["n"] > 0 else 0.0 for r in records])

allele_keys_C = sorted((r["pos1"], r["ref"], r["alt"]) for r in records)
invariant_ok = allele_keys_A == allele_keys_C
print(f"allele-call-set invariant vs A/B (same call set by construction): {invariant_ok}")
setA, setC = set(allele_keys_A), set(allele_keys_C)
print(f"|A|={len(setA)} |C|={len(setC)} A-C={len(setA - setC)} C-A={len(setC - setA)}")

HEADER_TMPL = """##fileformat=VCFv4.2
##source=ai_dna_analyzer_genotype_layer_validation_C
##contig=<ID=chr21,length=46709983>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth (passing-filter reads)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (binomial genotype posterior argmax over {{0/0,0/1,1/1}}, eps=0.01 fixed)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality (10*log10(P_best/P_second), capped at 99)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
"""


def write_vcf():
    lines = []
    for r in sorted(records, key=lambda r: r["pos1"]):
        gt = r["gt_C"]
        gq = r["gq_C"]
        lines.append(f"chr21\t{r['pos1']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\t"
                     f"DP={r['depth']}\tGT:GQ\t{gt}:{int(round(gq))}")
    path = f"{OUT}/vcfs/validation_C.vcf"
    with open(path, "w") as fh:
        fh.write(HEADER_TMPL)
        fh.write("\n".join(lines))
        if lines:
            fh.write("\n")
    print(f"wrote {path}: {len(lines)} records")
    return path


write_vcf()

subprocess.run(["bcftools", "sort", "-Oz", "-o", f"{OUT}/vcfs/validation_C.vcf.gz",
                f"{OUT}/vcfs/validation_C.vcf"], check=True, capture_output=True)
subprocess.run(["bcftools", "index", "-t", f"{OUT}/vcfs/validation_C.vcf.gz"],
               check=True, capture_output=True)
print("bgzip/tabix done for validation C")

summary = dict(
    n_records=len(records),
    n_skipped_no_alt=n_skipped_no_alt,
    eps_C=EPS,
    theta0_wins=theta0_wins,
    invariant_ok=bool(invariant_ok),
    invariant_setA_minus_setC=len(setA - setC),
    invariant_setC_minus_setA=len(setC - setA),
    method_C_runtime_seconds=c_runtime,
)
with open(f"{OUT}/cache/build_summary_validation_C.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print(json.dumps(summary, indent=2))

np.savez_compressed(
    f"{OUT}/cache/records_validation_C.npz",
    pos1=np.array([r["pos1"] for r in records]),
    ref=np.array([r["ref"] for r in records]),
    alt=np.array([r["alt"] for r in records]),
    depth=np.array([r["depth"] for r in records]),
    truth_gt=np.array([r["truth_gt"] if r["truth_gt"] is not None else "" for r in records]),
    gt_C=np.array([r["gt_C"] for r in records]),
    gq_C=np.array([r["gq_C"] for r in records]),
    c_flag_theta0=np.array([r["c_flag_theta0"] for r in records]),
    vaf=vaf,
    n=np.array([r["n"] for r in records]),
)
print("Saved cache/records_validation_C.npz")
