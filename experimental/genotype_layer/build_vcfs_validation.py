"""Build method A/B VCFs for the TRUE hold-out validation region
(chr21:32-44Mb), reusing the SAME allele calls (cascade_calls from
cache/validation_region.npz) and varying only the GT field per method.

tau=0.69 is FROZEN (fit on the dev region, chr21:30.0-30.47Mb) and is NOT
refit here -- this is a strict hold-out validation, per task instructions.

Method A: baseline, GT always 0/1.
Method B: VAF-threshold genotyping, GT = 1/1 if VAF > tau else 0/1, tau=0.69 fixed.
Method C is intentionally NOT run for this task.

Never adds/removes a called site -- only GT changes. Writes to
experimental/genotype_layer/vcfs/validation_{A,B}.vcf(.gz/.tbi) -- new
filenames, does not touch vcfs/A.vcf / B.vcf from the dev region.
"""
import gzip
import json
import subprocess
import sys
import time

sys.path.insert(0, "/home/mark/Documents/Projects/AI_DNA_ANALYZER")

import numpy as np

REPO = "/home/mark/Documents/Projects/AI_DNA_ANALYZER"
OUT = f"{REPO}/experimental/genotype_layer"
DATA = f"{REPO}/data/giab_hg002_chr21_12Mb"
TRUTH_VCF = f"{DATA}/hg002_chr21_32_44M.vcf.gz"

REF_TOKENS = ("N", "A", "C", "G", "T")
TAU = 0.69  # FROZEN, fit on dev region -- do not touch

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

for r in records:
    r["gt_A"] = "0/1"

t0 = time.perf_counter()
vaf = np.array([r["k"] / r["n"] if r["n"] > 0 else 0.0 for r in records])
for i, r in enumerate(records):
    r["gt_B"] = "1/1" if vaf[i] > TAU else "0/1"
b_runtime = time.perf_counter() - t0
print(f"Method B: tau={TAU} (FROZEN, not fit on this region), runtime={b_runtime*1000:.2f}ms")

allele_keys_B = sorted((r["pos1"], r["ref"], r["alt"]) for r in records)
invariant_ok = allele_keys_A == allele_keys_B
print(f"allele-call-set invariant across A/B: {invariant_ok}")
# independent set comparison, printed explicitly
setA, setB = set(allele_keys_A), set(allele_keys_B)
print(f"|A|={len(setA)} |B|={len(setB)} A-B={len(setA - setB)} B-A={len(setB - setA)}")

HEADER_TMPL = """##fileformat=VCFv4.2
##source=ai_dna_analyzer_genotype_layer_validation_{label}
##contig=<ID=chr21,length=46709983>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth (passing-filter reads)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype ({desc})">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002
"""

DESCS = {
    "A": "heuristic: always 0/1",
    "B": f"VAF-threshold: 1/1 if VAF>{TAU} else 0/1 (tau FROZEN from dev region, not refit)",
}


def write_vcf(method, gt_key):
    lines = []
    for r in sorted(records, key=lambda r: r["pos1"]):
        gt = r[gt_key]
        lines.append(f"chr21\t{r['pos1']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\t"
                     f"DP={r['depth']}\tGT\t{gt}")
    path = f"{OUT}/vcfs/validation_{method}.vcf"
    with open(path, "w") as fh:
        fh.write(HEADER_TMPL.format(label=method, desc=DESCS[method]))
        fh.write("\n".join(lines))
        if lines:
            fh.write("\n")
    print(f"wrote {path}: {len(lines)} records")
    return path


write_vcf("A", "gt_A")
write_vcf("B", "gt_B")

for method in ("A", "B"):
    subprocess.run(["bcftools", "sort", "-Oz", "-o", f"{OUT}/vcfs/validation_{method}.vcf.gz",
                    f"{OUT}/vcfs/validation_{method}.vcf"], check=True, capture_output=True)
    subprocess.run(["bcftools", "index", "-t", f"{OUT}/vcfs/validation_{method}.vcf.gz"],
                   check=True, capture_output=True)
print("bgzip/tabix done for validation A/B")

summary = dict(
    n_records=len(records),
    n_skipped_no_alt=n_skipped_no_alt,
    tau=TAU,
    tau_source="FROZEN from dev region chr21:30.0-30.47Mb, NOT refit here",
    invariant_ok=bool(invariant_ok),
    invariant_setA_minus_setB=len(setA - setB),
    invariant_setB_minus_setA=len(setB - setA),
    method_B_runtime_seconds=b_runtime,
)
with open(f"{OUT}/cache/build_summary_validation.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print(json.dumps(summary, indent=2))

np.savez_compressed(
    f"{OUT}/cache/records_validation.npz",
    pos1=np.array([r["pos1"] for r in records]),
    ref=np.array([r["ref"] for r in records]),
    alt=np.array([r["alt"] for r in records]),
    depth=np.array([r["depth"] for r in records]),
    truth_gt=np.array([r["truth_gt"] if r["truth_gt"] is not None else "" for r in records]),
    gt_A=np.array([r["gt_A"] for r in records]),
    gt_B=np.array([r["gt_B"] for r in records]),
    vaf=vaf,
    n=np.array([r["n"] for r in records]),
)
print("Saved cache/records_validation.npz")
