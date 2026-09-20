"""Genotype-level and allele-level metrics for optimized Method C applied to the
Stage-1 unified_happy scope (ai_pb_only_C.vcf, ai_cascade_C.vcf), plus the
call-set-invariant check (set of alleles before vs after genotyping).
"""
import gzip
import sys

import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
sys.path.insert(0, _ROOT)
import numpy as np

REPO = _ROOT
UNIFIED = f"{REPO}/experimental/unified_happy"
DATA = f"{REPO}/data/giab_hg002_chr21_12Mb"
TRUTH_VCF = f"{DATA}/hg002_chr21_32_44M.vcf.gz"

CLASSES = ["0/0", "0/1", "1/1"]


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


def parse_query_vcf(path):
    records = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, _id, ref, alt, qual, filt, info, fmt, sample = line.rstrip("\n").split("\t")
            fmt_fields = fmt.split(":")
            sample_fields = sample.split(":")
            gt = sample_fields[fmt_fields.index("GT")]
            records.append((chrom, int(pos), ref, alt, gt))
    return records


truth_gt = parse_truth_gt(TRUTH_VCF)
print(f"truth GT parsed: {len(truth_gt)} records")


def evaluate(vcf_path, forced_vcf_path, label):
    records = parse_query_vcf(vcf_path)
    forced_records = parse_query_vcf(forced_vcf_path)

    allele_keys_after = sorted((pos, ref, alt) for (_c, pos, ref, alt, _gt) in records)
    allele_keys_before = sorted((pos, ref, alt) for (_c, pos, ref, alt, _gt) in forced_records)
    set_before, set_after = set(allele_keys_before), set(allele_keys_after)
    invariant_ok = set_before == set_after
    print(f"[{label}] call-set invariant (before genotyping vs after Method C genotyping): "
          f"{invariant_ok} (|before|={len(set_before)} |after|={len(set_after)} "
          f"before-after={len(set_before - set_after)} after-before={len(set_after - set_before)})")

    cm = np.zeros((3, 3), dtype=int)
    n_known = 0
    n_no_truth = 0
    for chrom, pos, ref, alt, gt in records:
        gt_norm = gt.replace("|", "/")
        parts = set(gt_norm.split("/"))
        if parts == {"0"}:
            gt_class = "0/0"
        elif parts == {"1"}:
            gt_class = "1/1"
        else:
            gt_class = "0/1"
        t = truth_gt.get((chrom, pos))
        if t is None:
            n_no_truth += 1
            continue
        truth_class = t[2]
        cm[CLASSES.index(truth_class), CLASSES.index(gt_class)] += 1
        n_known += 1

    acc = np.trace(cm) / n_known if n_known else float("nan")
    per_class = {}
    for i, c in enumerate(CLASSES):
        denom = cm[i].sum()
        per_class[c] = (cm[i, i] / denom) if denom else float("nan")

    print(f"[{label}] n_known_truth={n_known} n_no_truth_match={n_no_truth}")
    print(f"[{label}] GT accuracy (overall, known truth)={acc:.4f}")
    print(f"[{label}] 0/1 accuracy={per_class['0/1']:.4f}  1/1 accuracy={per_class['1/1']:.4f}  "
          f"0/0 accuracy={per_class['0/0']:.4f}")
    print(f"[{label}] confusion matrix (rows=truth, cols=pred; order 0/0,0/1,1/1):")
    print(cm)

    return dict(label=label, n_known=n_known, n_no_truth=n_no_truth, gt_acc=float(acc),
                per_class_acc={k: float(v) for k, v in per_class.items()},
                confusion_matrix=cm.tolist(), invariant_ok=invariant_ok,
                n_before=len(set_before), n_after=len(set_after))


results = {}
results["ai_pb_only_C"] = evaluate(f"{UNIFIED}/ai_pb_only_C.vcf", f"{UNIFIED}/ai_pb_only.vcf", "ai_pb_only_C")
results["ai_cascade_C"] = evaluate(f"{UNIFIED}/ai_cascade_C.vcf", f"{UNIFIED}/ai_cascade.vcf", "ai_cascade_C")

import json
with open(f"{REPO}/experimental/genotype_layer/cache/unified_C_genotype_metrics.json", "w") as fh:
    json.dump(results, fh, indent=2)
print("\nWrote cache/unified_C_genotype_metrics.json")
