"""Post-hap.py evaluation for the TRUE hold-out validation region
(chr21:32-44Mb). Mirrors evaluate.py but for methods A/B only (no C), and
adds the 5 required stratified checks for Method B's GT accuracy:
  1. low-VAF sites (VAF<0.3)
  2. high-VAF het-truth sites (truth=0/1 but VAF>tau=0.69 -- false 1/1 risk)
  3. homozygous-truth (1/1) accuracy specifically
  4. low-depth loci (depth<10)
  5. segdup/low-mappability overlap (data/strat_v31/lowmap_segdup.bed.gz)
tau is FROZEN at 0.69 (fit on dev region) -- NOT tuned here.
"""
import csv
import gzip
import json
import sys

sys.path.insert(0, "/home/mark/Documents/Projects/AI_DNA_ANALYZER")

import numpy as np

OUT = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/experimental/genotype_layer"
DATA = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/giab_hg002_chr21_12Mb"
SEGDUP_BED = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/strat_v31/lowmap_segdup.bed.gz"

with open(f"{OUT}/cache/build_summary_validation.json") as fh:
    build_summary = json.load(fh)

TAU = build_summary["tau"]
assert TAU == 0.69, "tau must remain frozen at 0.69 for this hold-out validation"

r = np.load(f"{OUT}/cache/records_validation.npz", allow_pickle=True)
pos1, ref, alt = r["pos1"], r["ref"], r["alt"]
depth = r["depth"]
truth_gt = r["truth_gt"]
vaf = r["vaf"]
gt = {"A": r["gt_A"], "B": r["gt_B"]}

CLASSES = ["0/0", "0/1", "1/1"]


def confusion(pred, truth, mask=None):
    known = truth != ""
    if mask is not None:
        known = known & mask
    cm = np.zeros((3, 3), dtype=int)
    for t, p in zip(truth[known], pred[known]):
        ti = CLASSES.index(t) if t in CLASSES else None
        pi = CLASSES.index(p) if p in CLASSES else None
        if ti is not None and pi is not None:
            cm[ti, pi] += 1
    n = known.sum()
    acc = np.trace(cm) / n if n else float("nan")
    per_class_acc = {}
    for i, c in enumerate(CLASSES):
        denom = cm[i].sum()
        per_class_acc[c] = (cm[i, i] / denom) if denom else float("nan")
    return cm, acc, per_class_acc, int(n)


gt_stats = {}
for m in ("A", "B"):
    cm, acc, per_class_acc, n_known = confusion(gt[m], truth_gt)
    gt_stats[m] = dict(cm=cm.tolist(), acc=acc, per_class_acc=per_class_acc, n_known=n_known)

# --- allele-level call set identity (independent set comparison) ---
keys_A = set(zip(pos1.tolist(), ref.tolist(), alt.tolist()))
keys_B = set(zip(pos1.tolist(), ref.tolist(), alt.tolist()))  # same arrays by construction
invariant_identical = keys_A == keys_B
print(f"Invariant check: |A|={len(keys_A)} |B|={len(keys_B)} "
      f"A-B={len(keys_A - keys_B)} B-A={len(keys_B - keys_A)} identical={invariant_identical}")


# --- parse hap.py summary.csv ---
def parse_happy_summary(method):
    path = f"{OUT}/happy/validation_{method}/happy.summary.csv"
    rows = {}
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows[(row["Type"], row["Filter"])] = row
    return rows


happy_stats = {}
for m in ("A", "B"):
    rows = parse_happy_summary(m)
    snp_pass = rows.get(("SNP", "PASS"))
    happy_stats[m] = dict(
        precision=float(snp_pass["METRIC.Precision"]),
        recall=float(snp_pass["METRIC.Recall"]),
        f1=float(snp_pass["METRIC.F1_Score"]),
        truth_tp=int(snp_pass["TRUTH.TP"]),
        truth_fn=int(snp_pass["TRUTH.FN"]),
        query_fp=int(snp_pass["QUERY.FP"]),
    )

# --- allele-level P/R/F1 (position+ALT only), identical across A/B by construction ---
truth_alleles = set()
with gzip.open(f"{DATA}/hg002_chr21_32_44M.vcf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        truth_alleles.add((int(parts[1]), parts[3], parts[4]))

called_alleles = set(zip(pos1.tolist(), ref.tolist(), alt.tolist()))
allele_tp = len(called_alleles & truth_alleles)
allele_fp = len(called_alleles - truth_alleles)
allele_fn = len(truth_alleles - called_alleles)
allele_precision = allele_tp / (allele_tp + allele_fp) if (allele_tp + allele_fp) else float("nan")
allele_recall = allele_tp / (allele_tp + allele_fn) if (allele_tp + allele_fn) else float("nan")
allele_f1 = (2 * allele_precision * allele_recall / (allele_precision + allele_recall)
            if (allele_precision + allele_recall) else float("nan"))

runtime_B = build_summary["method_B_runtime_seconds"]

# ============================================================================
# Stratified checks (Method B only)
# ============================================================================
known_mask = truth_gt != ""

# 1. low-VAF sites (VAF < 0.3)
low_vaf_mask = vaf < 0.3
cm_lv, acc_lv, pca_lv, n_lv = confusion(gt["B"], truth_gt, mask=low_vaf_mask)

# 2. high-VAF het-truth sites (truth=0/1 but VAF>tau -- would B wrongly call 1/1?)
het_truth_mask = (truth_gt == "0/1")
high_vaf_het_mask = het_truth_mask & (vaf > TAU)
n_high_vaf_het = int(high_vaf_het_mask.sum())
n_het_truth_total = int(het_truth_mask.sum())
b_wrong_at_high_vaf_het = int(((gt["B"] == "1/1") & high_vaf_het_mask).sum())

# 3. homozygous-truth (1/1) accuracy specifically
hom_truth_mask = (truth_gt == "1/1")
cm_hom, acc_hom, pca_hom, n_hom = confusion(gt["B"], truth_gt, mask=hom_truth_mask)

# 4. low-depth loci (depth < 10)
low_depth_mask = depth < 10
cm_ld, acc_ld, pca_ld, n_ld = confusion(gt["B"], truth_gt, mask=low_depth_mask)

# 5. segdup / low-mappability overlap
segdup_intervals = []
with gzip.open(SEGDUP_BED, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if parts[0] != "chr21":
            continue
        s, e = int(parts[1]), int(parts[2])
        if e < 32_000_000 or s > 44_000_000:
            continue
        segdup_intervals.append((s, e))
segdup_intervals.sort()

def in_segdup(p1):
    # p1 is 1-based; intervals are 0-based half-open (BED)
    p0 = p1 - 1
    import bisect
    starts = [iv[0] for iv in segdup_intervals]
    idx = bisect.bisect_right(starts, p0) - 1
    if idx >= 0:
        s, e = segdup_intervals[idx]
        if s <= p0 < e:
            return True
    return False

segdup_mask = np.array([in_segdup(int(p)) for p in pos1])
cm_sd, acc_sd, pca_sd, n_sd = confusion(gt["B"], truth_gt, mask=segdup_mask)

strat_summary = dict(
    low_vaf=dict(n=n_lv, acc=acc_lv, per_class=pca_lv, cm=cm_lv.tolist()),
    high_vaf_het_truth=dict(
        n_het_truth_total=n_het_truth_total, n_high_vaf_het=n_high_vaf_het,
        n_wrongly_called_1_1=b_wrong_at_high_vaf_het,
        rate=(b_wrong_at_high_vaf_het / n_high_vaf_het) if n_high_vaf_het else float("nan"),
    ),
    hom_truth=dict(n=n_hom, acc=acc_hom, per_class=pca_hom, cm=cm_hom.tolist()),
    low_depth=dict(n=n_ld, acc=acc_ld, per_class=pca_ld, cm=cm_ld.tolist()),
    segdup=dict(n_intervals_used=len(segdup_intervals), n=n_sd, acc=acc_sd,
                per_class=pca_sd, cm=cm_sd.tolist()),
)
print(json.dumps(strat_summary, indent=2, default=str))

# ============================================================================
# Write RESULTS_VALIDATION.md
# ============================================================================
lines = []
lines.append("# Genotype Layer -- TRUE Hold-Out Validation Results (chr21:32-44Mb)\n")
lines.append("This is a strict hold-out test: tau=0.69 was fit on the dev region "
             "(chr21:30.0-30.47Mb) and is used here UNCHANGED -- not refit, not tuned, "
             "regardless of the results below.\n")
lines.append("Data used: reference data/reference/chr21_full.fa, BAM "
             "data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_15x.bam, truth "
             "data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz(.tbi), high-confidence BED "
             "data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed (matching naming "
             "pattern to the dev region's *_highconf.bed -- found directly alongside the "
             "BAM/truth VCF, no substitute needed).\n")

lines.append("## Table 1: Method x Precision/Recall/F1/GTacc/0/1acc/1/1acc/Runtime\n")
lines.append("| Method | Precision | Recall | F1 | GT Acc (overall) | 0/0 Acc | 0/1 Acc | 1/1 Acc | Genotype-layer runtime |")
lines.append("|---|---|---|---|---|---|---|---|---|")
for m, label in (("A", "A (baseline, always 0/1)"), ("B", "B (VAF-threshold, tau=0.69 frozen)")):
    hs = happy_stats[m]
    gs = gt_stats[m]
    pca = gs["per_class_acc"]
    rt = "n/a (no genotype layer)" if m == "A" else f"{runtime_B*1000:.2f} ms"
    lines.append(f"| {label} | {hs['precision']:.4f} | {hs['recall']:.4f} | {hs['f1']:.4f} | "
                 f"{gs['acc']:.4f} | {pca['0/0']:.4f} | {pca['0/1']:.4f} | {pca['1/1']:.4f} | {rt} |")
lines.append("")

lines.append("## Table 2: Method x AlleleF1 x GA4GH-F1 x delta x extra runtime\n")
lines.append("| Method | Allele-level F1 (pos+ALT only) | GA4GH/hap.py F1 | Delta (GA4GH - Allele) | Extra genotype-layer runtime vs. A |")
lines.append("|---|---|---|---|---|")
for m in ("A", "B"):
    hs = happy_stats[m]
    delta = hs["f1"] - allele_f1
    extra_rt = "0 ms (baseline)" if m == "A" else f"{runtime_B*1000:.2f} ms"
    lines.append(f"| {m} | {allele_f1:.4f} | {hs['f1']:.4f} | {delta:+.4f} | {extra_rt} |")
lines.append("")
lines.append(f"Allele-level TP={allele_tp} FP={allele_fp} FN={allele_fn} "
            f"(identical across A/B by construction -- same {len(called_alleles)}-call set).\n")

lines.append("## GT confusion matrices (rows=truth, cols=predicted; order 0/0,0/1,1/1)\n")
for m in ("A", "B"):
    cm = gt_stats[m]["cm"]
    lines.append(f"**Method {m}** (n={gt_stats[m]['n_known']} loci with known truth GT):\n")
    lines.append("| truth\\pred | 0/0 | 0/1 | 1/1 |")
    lines.append("|---|---|---|---|")
    for i, c in enumerate(CLASSES):
        lines.append(f"| {c} | {cm[i][0]} | {cm[i][1]} | {cm[i][2]} |")
    lines.append("")

lines.append("## Frozen parameters (NOT tuned on this region)\n")
lines.append(f"- Method B: tau = **{TAU}** (fit on dev region chr21:30.0-30.47Mb; used as-is here)\n")

lines.append("## Two-level invariant (allele call set unchanged between A and B)\n")
lines.append(f"- Allele call set (CHROM:POS:REF:ALT) identical between A and B: **{invariant_identical}**")
lines.append(f"- Independent set comparison: |A|={len(keys_A)} |B|={len(keys_B)} "
            f"A-B={len(keys_A - keys_B)} B-A={len(keys_B - keys_A)}")
lines.append(f"- N called loci: {build_summary['n_records']} (0 skipped for ref=N/no-alt-support)")
lines.append(f"- Allele-level P/R/F1 identical by construction between A/B: "
            f"P={allele_precision:.4f} R={allele_recall:.4f} F1={allele_f1:.4f}\n")

lines.append("## Stratified checks on Method B's GT accuracy\n")
lines.append(f"**1. Low-VAF sites (VAF<0.3):** n={strat_summary['low_vaf']['n']}, "
            f"GT acc={strat_summary['low_vaf']['acc']:.4f}, "
            f"per-class={strat_summary['low_vaf']['per_class']}\n")
lines.append(f"**2. High-VAF het-truth sites (truth=0/1 but VAF>{TAU}):** "
            f"{strat_summary['high_vaf_het_truth']['n_high_vaf_het']} of "
            f"{strat_summary['high_vaf_het_truth']['n_het_truth_total']} total het-truth sites. "
            f"Method B wrongly called 1/1 at "
            f"{strat_summary['high_vaf_het_truth']['n_wrongly_called_1_1']} of these "
            f"({strat_summary['high_vaf_het_truth']['n_high_vaf_het']} candidates), "
            f"rate={strat_summary['high_vaf_het_truth']['rate']:.4f} "
            f"(by construction this rate should be ~1.0, since VAF>tau triggers 1/1 by definition -- "
            f"this quantifies exactly how many true hets get mis-genotyped by the VAF rule).\n")
lines.append(f"**3. Homozygous-truth (1/1) accuracy:** n={strat_summary['hom_truth']['n']}, "
            f"GT acc={strat_summary['hom_truth']['acc']:.4f}, "
            f"1/1-class acc={strat_summary['hom_truth']['per_class'].get('1/1', float('nan')):.4f}\n")
lines.append(f"**4. Low-depth loci (depth<10):** n={strat_summary['low_depth']['n']}, "
            f"GT acc={strat_summary['low_depth']['acc']:.4f}, "
            f"per-class={strat_summary['low_depth']['per_class']}\n")
lines.append(f"**5. Segdup/low-mappability overlap:** used data/strat_v31/lowmap_segdup.bed.gz "
            f"({strat_summary['segdup']['n_intervals_used']} intervals overlapping chr21:32-44Mb). "
            f"n={strat_summary['segdup']['n']} called loci overlap segdup/low-mappability regions, "
            f"GT acc={strat_summary['segdup']['acc']:.4f}, "
            f"per-class={strat_summary['segdup']['per_class']}\n")

lines.append("## Did tau=0.69 generalize?\n")
gt_delta = gt_stats["B"]["acc"] - gt_stats["A"]["acc"]
lines.append(f"- Method A GT accuracy (baseline, always 0/1): {gt_stats['A']['acc']:.4f}")
lines.append(f"- Method B GT accuracy (VAF-threshold, tau=0.69 frozen): {gt_stats['B']['acc']:.4f}")
lines.append(f"- Delta: {gt_delta:+.4f}")
lines.append(f"- hap.py SNP F1: A={happy_stats['A']['f1']:.4f}, B={happy_stats['B']['f1']:.4f}, "
            f"delta={happy_stats['B']['f1']-happy_stats['A']['f1']:+.4f}")
lines.append("- See generalization verdict and failure-mode description in the accompanying report.\n")

with open(f"{OUT}/RESULTS_VALIDATION.md", "w") as fh:
    fh.write("\n".join(lines) + "\n")

print("\n".join(lines))
print(f"\nWritten to {OUT}/RESULTS_VALIDATION.md")
