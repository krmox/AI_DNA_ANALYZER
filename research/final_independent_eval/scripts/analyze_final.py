#!/usr/bin/env python3
"""Descriptive metrics for one hap.py run + block bootstrap (whole 2-Mb blocks, 4000 resamples, seed 20260926).
usage: analyze_final.py PREFIX_OF_HAPPY_OUTPUT OUT_JSON   (PREFIX like .../ai_A_C_chr8)"""
import csv, json, sys, gzip
import numpy as np, pysam
pre, outj = sys.argv[1], sys.argv[2]
BLOCK = 2_000_000; NBOOT = 4000; SEED = 20260926
summ = {}
for r in csv.DictReader(open(pre + ".summary.csv")):
    summ[(r["Type"], r["Filter"])] = r
snp = {f: summ[("SNP", f)] for f in ("ALL", "PASS")}
blocks = {}
def add(b, k):
    blocks.setdefault(b, dict(tpT=0, tpQ=0, fp=0, fn=0))[k] += 1
vf = pysam.VariantFile(pre + ".vcf.gz")
n_rec = 0
for rec in vf:
    n_rec += 1
    b = rec.pos // BLOCK
    for sname, s in rec.samples.items():
        if s.get("BVT") != "SNP":
            continue
        bd = s.get("BD")
        if sname == "TRUTH":
            if bd == "TP": add(b, "tpT")
            elif bd == "FN": add(b, "fn")
        elif sname == "QUERY":
            if bd == "TP": add(b, "tpQ")
            elif bd == "FP": add(b, "fp")
keys = sorted(blocks)
A = np.array([[blocks[k][c] for c in ("tpT", "tpQ", "fp", "fn")] for k in keys], float)
def prf(m):
    tpT, tpQ, fp, fn = m.sum(0)
    P = tpQ / (tpQ + fp); R = tpT / (tpT + fn); return P, R, 2 * P * R / (P + R)
rng = np.random.default_rng(SEED)
res = np.array([prf(A[rng.integers(0, len(A), len(A))]) for _ in range(NBOOT)])
ci = {n: [float(np.percentile(res[:, i], 2.5)), float(np.percentile(res[:, i], 97.5))] for i, n in enumerate(("precision", "recall", "f1"))}
P, R, F = prf(A)
out = dict(happy_summary_SNP=snp, n_vcf_records=n_rec, n_blocks=len(A), block_size=BLOCK, n_boot=NBOOT, seed=SEED,
           recomputed_from_vcf=dict(precision=P, recall=R, f1=F, TP_truth=int(A[:, 0].sum()), TP_query=int(A[:, 1].sum()), FP=int(A[:, 2].sum()), FN=int(A[:, 3].sum())),
           block_bootstrap_95ci=ci)
json.dump(out, open(outj, "w"), indent=1)
print(json.dumps(out, indent=1))
