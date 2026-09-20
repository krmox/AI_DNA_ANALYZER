"""Evaluation of Method C (binomial genotype posterior) on the TRUE hold-out
validation region (chr21:32-44Mb), computed identically to evaluate_validation.py
(methods A/B), plus the error-type breakdown table requested for B vs C.
Pulls A/B numbers from RESULTS_VALIDATION.md verbatim (does not recompute).
"""
import bisect
import gzip
import json
import sys

import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
sys.path.insert(0, _ROOT)

import numpy as np

OUT = _ROOT + "/experimental/genotype_layer"
DATA = _ROOT + "/data/giab_hg002_chr21_12Mb"
SEGDUP_BED = _ROOT + "/data/strat_v31/lowmap_segdup.bed.gz"
TANDEM_BED = _ROOT + "/data/strat_v31/tandemrepeats.bed.gz"

with open(f"{OUT}/cache/build_summary_validation.json") as fh:
    build_summary_AB = json.load(fh)
with open(f"{OUT}/cache/build_summary_validation_C.json") as fh:
    build_summary_C = json.load(fh)

r_ab = np.load(f"{OUT}/cache/records_validation.npz", allow_pickle=True)
r_c = np.load(f"{OUT}/cache/records_validation_C.npz", allow_pickle=True)

# Sanity: same pos1/ref/alt ordering issue -- align by (pos1,ref,alt) not by index,
# since both scripts iterate the SAME idxs from cascade_calls in the same order,
# but let's not assume order identical; build dict keyed by (pos1, ref, alt).
def to_dict(r, gt_key, extra=None):
    d = {}
    for i in range(len(r["pos1"])):
        key = (int(r["pos1"][i]), str(r["ref"][i]), str(r["alt"][i]))
        val = dict(truth_gt=str(r["truth_gt"][i]), gt=str(r[gt_key][i]),
                   depth=int(r["depth"][i]), vaf=float(r["vaf"][i]))
        if extra:
            for k in extra:
                val[k] = r[k][i]
        d[key] = val
    return d

dict_B = to_dict(r_ab, "gt_B")
dict_C = to_dict(r_c, "gt_C", extra=["gq_C", "c_flag_theta0"])

keys_B = set(dict_B.keys())
keys_C = set(dict_C.keys())
print(f"Invariant: |B|={len(keys_B)} |C|={len(keys_C)} B-C={len(keys_B-keys_C)} C-B={len(keys_C-keys_B)}")
assert keys_B == keys_C, "Method C call set must be identical to Method B's"

CLASSES = ["0/0", "0/1", "1/1"]


def confusion_from_dict(d, mask_fn=None):
    cm = np.zeros((3, 3), dtype=int)
    n = 0
    for key, v in d.items():
        if v["truth_gt"] == "":
            continue
        if mask_fn is not None and not mask_fn(v):
            continue
        t, p = v["truth_gt"], v["gt"]
        if t in CLASSES and p in CLASSES:
            cm[CLASSES.index(t), CLASSES.index(p)] += 1
            n += 1
    acc = np.trace(cm) / n if n else float("nan")
    per_class_acc = {}
    for i, c in enumerate(CLASSES):
        denom = cm[i].sum()
        per_class_acc[c] = (cm[i, i] / denom) if denom else float("nan")
    return cm, acc, per_class_acc, n


cm_C, acc_C, pca_C, n_known_C = confusion_from_dict(dict_C)
cm_B, acc_B, pca_B, n_known_B = confusion_from_dict(dict_B)
print(f"Method C GT acc={acc_C:.4f} n={n_known_C}; Method B GT acc={acc_B:.4f} n={n_known_B}")

# --- hap.py summary for C ---
def parse_happy_summary(method):
    path = f"{OUT}/happy/validation_{method}/happy.summary.csv"
    rows = {}
    import csv
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows[(row["Type"], row["Filter"])] = row
    return rows

rows_C = parse_happy_summary("C")
snp_pass_C = rows_C[("SNP", "PASS")]
happy_C = dict(
    precision=float(snp_pass_C["METRIC.Precision"]),
    recall=float(snp_pass_C["METRIC.Recall"]),
    f1=float(snp_pass_C["METRIC.F1_Score"]),
    truth_tp=int(snp_pass_C["TRUTH.TP"]),
    truth_fn=int(snp_pass_C["TRUTH.FN"]),
    query_fp=int(snp_pass_C["QUERY.FP"]),
)

# --- allele-level P/R/F1 for C (identical to A/B by construction: same call set) ---
truth_alleles = set()
with gzip.open(f"{DATA}/hg002_chr21_32_44M.vcf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        truth_alleles.add((int(parts[1]), parts[3], parts[4]))

called_alleles_C = keys_C
allele_tp_C = len(called_alleles_C & truth_alleles)
allele_fp_C = len(called_alleles_C - truth_alleles)
allele_fn_C = len(truth_alleles - called_alleles_C)
allele_p_C = allele_tp_C / (allele_tp_C + allele_fp_C) if (allele_tp_C + allele_fp_C) else float("nan")
allele_r_C = allele_tp_C / (allele_tp_C + allele_fn_C) if (allele_tp_C + allele_fn_C) else float("nan")
allele_f1_C = (2 * allele_p_C * allele_r_C / (allele_p_C + allele_r_C)) if (allele_p_C + allele_r_C) else float("nan")

runtime_C = build_summary_C["method_C_runtime_seconds"]

# ============================================================================
# Stratified checks for C -- mirrors evaluate_validation.py methodology exactly
# ============================================================================
def mask_low_vaf(v):
    return v["vaf"] < 0.3

def mask_het_03_05(v):
    return v["truth_gt"] == "0/1" and 0.3 <= v["vaf"] < 0.5

def mask_het_05_069(v):
    return v["truth_gt"] == "0/1" and 0.5 <= v["vaf"] <= 0.69

def mask_het_gt069(v):
    return v["truth_gt"] == "0/1" and v["vaf"] > 0.69

def mask_hom(v):
    return v["truth_gt"] == "1/1"

def mask_hom_nearvaf1(v):
    return v["truth_gt"] == "1/1" and v["vaf"] >= 0.95

def mask_hom_lowcov(v):
    return v["truth_gt"] == "1/1" and v["depth"] < 10

def mask_depth_lt10(v):
    return v["depth"] < 10

def mask_depth_10_19(v):
    return 10 <= v["depth"] < 20

def mask_depth_20_29(v):
    return 20 <= v["depth"] < 30

def mask_depth_30plus(v):
    return v["depth"] >= 30

strat_masks = dict(
    low_vaf=mask_low_vaf,
    het_vaf_lt03=mask_low_vaf,
    het_vaf_03_05=mask_het_03_05,
    het_vaf_05_069=mask_het_05_069,
    het_vaf_gt069=mask_het_gt069,
    hom_truth=mask_hom,
    hom_nearvaf1=mask_hom_nearvaf1,
    hom_lowcov=mask_hom_lowcov,
    depth_lt10=mask_depth_lt10,
    depth_10_19=mask_depth_10_19,
    depth_20_29=mask_depth_20_29,
    depth_30plus=mask_depth_30plus,
)

strat_results_C = {}
strat_results_B = {}
for name, fn in strat_masks.items():
    cm_c, acc_c, pca_c, n_c = confusion_from_dict(dict_C, mask_fn=fn)
    cm_b, acc_b, pca_b, n_b = confusion_from_dict(dict_B, mask_fn=fn)
    strat_results_C[name] = dict(n=n_c, acc=acc_c, per_class=pca_c, cm=cm_c.tolist())
    strat_results_B[name] = dict(n=n_b, acc=acc_b, per_class=pca_b, cm=cm_b.tolist())

# --- segdup / low-mappability ---
def load_intervals(bed_path):
    intervals = []
    with gzip.open(bed_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] != "chr21":
                continue
            s, e = int(parts[1]), int(parts[2])
            if e < 32_000_000 or s > 44_000_000:
                continue
            intervals.append((s, e))
    intervals.sort()
    return intervals

segdup_intervals = load_intervals(SEGDUP_BED)
tandem_intervals = load_intervals(TANDEM_BED)


def make_in_region(intervals):
    starts = [iv[0] for iv in intervals]
    def in_region(p1):
        p0 = p1 - 1
        idx = bisect.bisect_right(starts, p0) - 1
        if idx >= 0:
            s, e = intervals[idx]
            if s <= p0 < e:
                return True
        return False
    return in_region

in_segdup = make_in_region(segdup_intervals)
in_tandem = make_in_region(tandem_intervals)

def mask_segdup(v, pos1):
    return in_segdup(pos1)

def mask_tandem(v, pos1):
    return in_tandem(pos1)

def confusion_region(d, region_fn):
    cm = np.zeros((3, 3), dtype=int)
    n = 0
    for key, v in d.items():
        if v["truth_gt"] == "":
            continue
        if not region_fn(v, key[0]):
            continue
        t, p = v["truth_gt"], v["gt"]
        if t in CLASSES and p in CLASSES:
            cm[CLASSES.index(t), CLASSES.index(p)] += 1
            n += 1
    acc = np.trace(cm) / n if n else float("nan")
    per_class_acc = {}
    for i, c in enumerate(CLASSES):
        denom = cm[i].sum()
        per_class_acc[c] = (cm[i, i] / denom) if denom else float("nan")
    return cm, acc, per_class_acc, n

cm_sd_C, acc_sd_C, pca_sd_C, n_sd_C = confusion_region(dict_C, mask_segdup)
cm_sd_B, acc_sd_B, pca_sd_B, n_sd_B = confusion_region(dict_B, mask_segdup)
cm_td_C, acc_td_C, pca_td_C, n_td_C = confusion_region(dict_C, mask_tandem)
cm_td_B, acc_td_B, pca_td_B, n_td_B = confusion_region(dict_B, mask_tandem)

# ============================================================================
# Error-type breakdown: 0/1->1/1, 1/1->0/1, other (for B and C)
# ============================================================================
def error_breakdown(d):
    counts = {"0/1->1/1": 0, "1/1->0/1": 0, "other": 0}
    for key, v in d.items():
        t, p = v["truth_gt"], v["gt"]
        if t == "" or t == p:
            continue
        if t not in CLASSES or p not in CLASSES:
            continue
        if t == "0/1" and p == "1/1":
            counts["0/1->1/1"] += 1
        elif t == "1/1" and p == "0/1":
            counts["1/1->0/1"] += 1
        else:
            counts["other"] += 1
    return counts

err_B = error_breakdown(dict_B)
err_C = error_breakdown(dict_C)
print("Error breakdown B:", err_B)
print("Error breakdown C:", err_C)

# ============================================================================
# theta0_wins flag stats (C-specific)
# ============================================================================
theta0_flagged = [k for k, v in dict_C.items() if v["c_flag_theta0"]]
theta0_correct = sum(1 for k in theta0_flagged if dict_C[k]["truth_gt"] == "0/1")
theta0_truth_00 = sum(1 for k in theta0_flagged if dict_C[k]["truth_gt"] == "0/0")
print(f"theta0-flagged sites: {len(theta0_flagged)}; of these truth=0/1: {theta0_correct}, truth=0/0: {theta0_truth_00}")

# ============================================================================
# Known A/B numbers pulled VERBATIM from RESULTS_VALIDATION.md (not recomputed)
# ============================================================================
KNOWN = dict(
    A=dict(precision=0.6382, recall=0.6089, f1=0.6232, gt_acc=0.6464,
           acc00=float("nan"), acc01=1.0000, acc11=0.0000,
           allele_f1=0.8644, runtime_ms=0.0),
    B=dict(precision=0.9135, recall=0.8716, f1=0.8921, gt_acc=0.9250,
           acc00=float("nan"), acc01=0.8840, acc11=0.9998,
           allele_f1=0.8644, runtime_ms=15.68),
)

runtime_C_ms = runtime_C * 1000

# ============================================================================
# Write RESULTS_VALIDATION_C.md
# ============================================================================
lines = []
lines.append("# Genotype Layer -- Method C (Binomial Genotype Posterior) on TRUE Hold-Out Validation (chr21:32-44Mb)\n")
lines.append("Method C reuses the EXACT SAME 16,124 cascade calls used for A/B "
             "(cache/validation_region.npz) and applies the identical Method C logic "
             "from build_vcfs.py verbatim: eps=0.01 fixed (not fit), theta in "
             "{0.0, 0.5, 1-eps} representing {0/0, 0/1, 1/1}, GT=argmax likelihood, "
             "GQ=10*log10(P_best/P_second) capped at 99. Sites where theta=0 (0/0) wins "
             "are NEVER dropped -- kept as GT 0/1 with low GQ and flagged "
             "(never-suppress-a-detection invariant).\n")
lines.append("A and B numbers below are pulled verbatim from RESULTS_VALIDATION.md "
             "(not recomputed here).\n")

lines.append("## Table 1: A/B/C x GT Acc x 0/1 Acc x 1/1 Acc x Allele F1 x hap.py F1 x Runtime\n")
lines.append("| Method | GT Acc (overall) | 0/1 Acc | 1/1 Acc | Allele-level F1 | hap.py SNP F1 | Genotype-layer runtime |")
lines.append("|---|---|---|---|---|---|---|")
lines.append(f"| A (baseline, always 0/1) | {KNOWN['A']['gt_acc']:.4f} | {KNOWN['A']['acc01']:.4f} | {KNOWN['A']['acc11']:.4f} | {KNOWN['A']['allele_f1']:.4f} | {KNOWN['A']['f1']:.4f} | n/a |")
lines.append(f"| B (VAF-threshold, tau=0.69 frozen) | {KNOWN['B']['gt_acc']:.4f} | {KNOWN['B']['acc01']:.4f} | {KNOWN['B']['acc11']:.4f} | {KNOWN['B']['allele_f1']:.4f} | {KNOWN['B']['f1']:.4f} | {KNOWN['B']['runtime_ms']:.2f} ms |")
lines.append(f"| C (binomial posterior, eps=0.01 fixed) | {acc_C:.4f} | {pca_C['0/1']:.4f} | {pca_C['1/1']:.4f} | {allele_f1_C:.4f} | {happy_C['f1']:.4f} | {runtime_C_ms:.2f} ms |")
lines.append("")

lines.append("## Error-type breakdown (truth->pred mismatches, B vs C)\n")
lines.append("| Error type | B | C |")
lines.append("|---|---|---|")
lines.append(f"| 0/1 -> 1/1 | {err_B['0/1->1/1']} | {err_C['0/1->1/1']} |")
lines.append(f"| 1/1 -> 0/1 | {err_B['1/1->0/1']} | {err_C['1/1->0/1']} |")
lines.append(f"| other | {err_B['other']} | {err_C['other']} |")
lines.append(f"| **total errors** | {sum(err_B.values())} | {sum(err_C.values())} |")
lines.append("")

lines.append("## GT confusion matrix, Method C (rows=truth, cols=predicted; order 0/0,0/1,1/1)\n")
lines.append(f"n={n_known_C} loci with known truth GT:\n")
lines.append("| truth\\pred | 0/0 | 0/1 | 1/1 |")
lines.append("|---|---|---|---|")
for i, c in enumerate(CLASSES):
    lines.append(f"| {c} | {cm_C[i][0]} | {cm_C[i][1]} | {cm_C[i][2]} |")
lines.append("")
lines.append("## GT confusion matrix, Method B (for direct comparison; from records_validation.npz)\n")
lines.append(f"n={n_known_B} loci with known truth GT:\n")
lines.append("| truth\\pred | 0/0 | 0/1 | 1/1 |")
lines.append("|---|---|---|---|")
for i, c in enumerate(CLASSES):
    lines.append(f"| {c} | {cm_B[i][0]} | {cm_B[i][1]} | {cm_B[i][2]} |")
lines.append("")

lines.append("## Direct comparison: what changed between B and C\n")
lines.append(f"- GT accuracy: B={acc_B:.4f}, C={acc_C:.4f}, delta={acc_C-acc_B:+.4f}")
lines.append(f"- 0/1 accuracy: B={pca_B['0/1']:.4f}, C={pca_C['0/1']:.4f}, delta={pca_C['0/1']-pca_B['0/1']:+.4f}")
lines.append(f"- 1/1 accuracy: B={pca_B['1/1']:.4f}, C={pca_C['1/1']:.4f}, delta={pca_C['1/1']-pca_B['1/1']:+.4f}")
lines.append(f"- hap.py SNP F1: B={KNOWN['B']['f1']:.4f}, C={happy_C['f1']:.4f}, delta={happy_C['f1']-KNOWN['B']['f1']:+.4f}")
lines.append(f"- hap.py Precision: B={KNOWN['B']['precision']:.4f}, C={happy_C['precision']:.4f}, delta={happy_C['precision']-KNOWN['B']['precision']:+.4f}")
lines.append(f"- hap.py Recall: B={KNOWN['B']['recall']:.4f}, C={happy_C['recall']:.4f}, delta={happy_C['recall']-KNOWN['B']['recall']:+.4f}")
lines.append(f"- Error count: B={sum(err_B.values())}, C={sum(err_C.values())}, delta={sum(err_C.values())-sum(err_B.values()):+d}")
lines.append(f"- theta=0(0/0)-favored sites (C only): {build_summary_C['theta0_wins']} of {build_summary_C['n_records']} "
             f"({100*build_summary_C['theta0_wins']/build_summary_C['n_records']:.3f}%); flagged, kept as GT 0/1, low GQ, never dropped. "
             f"Of these, truth=0/1 at {theta0_correct}, truth=0/0 at {theta0_truth_00} (truth=0/0 not in call set otherwise by construction -- these are candidate false detections upstream of the genotype layer).")
lines.append("")

lines.append("## Stratified GT accuracy, C vs B\n")
strat_labels = dict(
    het_vaf_lt03="Het truth, VAF<0.3",
    het_vaf_03_05="Het truth, VAF 0.3-0.5",
    het_vaf_05_069="Het truth, VAF 0.5-0.69",
    het_vaf_gt069="Het truth, VAF>0.69",
    hom_nearvaf1="Hom truth (1/1), VAF>=0.95 (near-VAF-1)",
    hom_lowcov="Hom truth (1/1), depth<10 (low-coverage)",
    depth_lt10="Depth<10",
    depth_10_19="Depth 10-19",
    depth_20_29="Depth 20-29",
    depth_30plus="Depth 30+",
)
lines.append("| Stratum | n | B GT acc | C GT acc | Delta |")
lines.append("|---|---|---|---|---|")
for key, label in strat_labels.items():
    b = strat_results_B[key]
    c = strat_results_C[key]
    b_acc = b["acc"]
    c_acc = c["acc"]
    delta = (c_acc - b_acc) if not (np.isnan(b_acc) or np.isnan(c_acc)) else float("nan")
    lines.append(f"| {label} | {c['n']} | {b_acc:.4f} | {c_acc:.4f} | {delta:+.4f} |")
lines.append(f"| Segdup/low-mappability overlap | {n_sd_C} | {acc_sd_B:.4f} | {acc_sd_C:.4f} | {(acc_sd_C-acc_sd_B):+.4f} |")
lines.append(f"| Tandem repeat overlap | {n_td_C} | {acc_td_B:.4f} | {acc_td_C:.4f} | {(acc_td_C-acc_td_B):+.4f} |")
lines.append("")
lines.append(f"Difficult-region strat used: data/strat_v31/lowmap_segdup.bed.gz "
             f"({len(segdup_intervals)} intervals overlapping chr21:32-44Mb) and "
             f"data/strat_v31/tandemrepeats.bed.gz ({len(tandem_intervals)} intervals overlapping chr21:32-44Mb). "
             f"No local GIAB HG004 stratification BED applicable to this chr21/HG002 region was found under "
             f"data/giab_hg004_v19/ (that directory holds chr2/chr3/chr5 HG004 test data unrelated to this analysis).\n")

lines.append("## Efficiency: C vs B\n")
lines.append(f"- Method B total genotype-layer runtime: {KNOWN['B']['runtime_ms']:.2f} ms "
             f"({KNOWN['B']['runtime_ms']/build_summary_C['n_records']*1000:.4f} us/site)")
lines.append(f"- Method C total genotype-layer runtime: {runtime_C_ms:.2f} ms "
             f"({runtime_C_ms/build_summary_C['n_records']*1000:.4f} us/site)")
lines.append(f"- Runtime ratio C/B: {runtime_C_ms/KNOWN['B']['runtime_ms']:.1f}x")
lines.append(f"- Both measured identically: wall-clock time.perf_counter() around the pure-numpy/scipy "
             f"genotyping loop over already-extracted (k, n) counts, per-site, excluding BAM-scanning/detection "
             f"(a cost already paid identically for A/B/C by the shared cascade extraction).\n")

lines.append("## Invariant confirmation\n")
lines.append(f"- Call-set size before (B): {len(keys_B)}")
lines.append(f"- Call-set size after (C): {len(keys_C)}")
lines.append(f"- Symmetric difference B vs C: |B-C|={len(keys_B-keys_C)}, |C-B|={len(keys_C-keys_B)} (both zero -> identical call sets)")
lines.append(f"- n_records in C build: {build_summary_C['n_records']}, n_skipped_no_alt: {build_summary_C['n_skipped_no_alt']}")
lines.append("")

with open(f"{OUT}/RESULTS_VALIDATION_C.md", "w") as fh:
    fh.write("\n".join(lines) + "\n")

print(f"\nWritten to {OUT}/RESULTS_VALIDATION_C.md")
print("\n".join(lines[:40]))
