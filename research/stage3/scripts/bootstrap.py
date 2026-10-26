#!/usr/bin/env python3
"""2-Mb block bootstrap of F1/recall/precision per condition and the paired Varion - bcftools difference.

Per-record TP/FP/FN come from the hap.py annotated VCF (BD/BVT); SNP records only; blocks = floor((pos-1)/2e6).
4,000 resamples, seed 20261027. Blocks are not independent replicates; the interval reflects spatial sampling variability within one contig.
Output: BOOTSTRAP.tsv (one row per condition x quantity).
"""
import csv, subprocess, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import stage3 as s3  # noqa: E402

BLOCK = 2_000_000
RESAMPLES = 4000
OUT = Path("/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage3/BOOTSTRAP.tsv")


def block_counts(vcf):
    text = subprocess.run(["bcftools", "query", "-f", "%POS\t[%BD,%BVT\t]\n", str(vcf)], capture_output=True, text=True, check=True).stdout
    blocks = {}
    for line in text.splitlines():
        pos, truth, query = line.rstrip("\t").split("\t")
        acc = blocks.setdefault((int(pos) - 1) // BLOCK, [0, 0, 0])
        bd, vt = query.split(",")
        if vt == "SNP" and bd == "TP":
            acc[0] += 1
        elif vt == "SNP" and bd == "FP":
            acc[1] += 1
        bd, vt = truth.split(",")
        if vt == "SNP" and bd == "FN":
            acc[2] += 1
    return blocks


def metrics(a):
    tp, fp, fn = a.sum(axis=0)
    return np.array([2 * tp / (2 * tp + fp + fn), tp / (tp + fn), tp / (tp + fp)])


def main():
    rng = np.random.default_rng(20261027)
    rows = []
    for run_dir in sorted((s3.S3 / "runs").iterdir()):
        v = run_dir / "happy_varion" / "h.vcf.gz"
        b = run_dir / "happy_bcftools" / "h.vcf.gz"
        if not (v.exists() and b.exists()):
            continue
        cond = run_dir.name
        va, ba = block_counts(v), block_counts(b)
        keys = sorted(set(va) | set(ba))
        X = np.array([va.get(k, [0, 0, 0]) for k in keys])
        Y = np.array([ba.get(k, [0, 0, 0]) for k in keys])
        n = len(keys)
        idx = [rng.integers(0, n, n) for _ in range(RESAMPLES)]
        mv = np.array([metrics(X[i]) for i in idx])
        mb = np.array([metrics(Y[i]) for i in idx])
        point_v, point_b = metrics(X), metrics(Y)
        names = ("F1", "recall", "precision")
        for j, q in enumerate(names):
            for label, point, dist in (("Varion", point_v[j], mv[:, j]), ("bcftools", point_b[j], mb[:, j]), ("Varion-bcftools", point_v[j] - point_b[j], mv[:, j] - mb[:, j])):
                lo, hi = np.percentile(dist, [2.5, 97.5])
                rows.append(dict(condition=cond, quantity=q, series=label, estimate=f"{point:.6f}", ci_lo=f"{lo:.6f}", ci_hi=f"{hi:.6f}", n_blocks=n, resamples=RESAMPLES))
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, ["condition", "quantity", "series", "estimate", "ci_lo", "ci_hi", "n_blocks", "resamples"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(len(rows), "rows ->", OUT)


if __name__ == "__main__":
    main()
