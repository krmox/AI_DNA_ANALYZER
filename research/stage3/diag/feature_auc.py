#!/usr/bin/env python3
"""Single-feature AUC (true variant vs not) among k>=3 candidates, overall and inside GIAB alldifficultregions. Diagnostic only."""
import csv, gzip, sys, bisect
import numpy as np
D = sys.argv[1]
S = "/mnt/archive/AI_DNA_ANALYZER_stage3/strat/"
truth = {}
for r in csv.DictReader(open(D + "/records.tsv"), delimiter="\t"):
    if r["truth_bd"] in ("TP", "FN"): truth[int(r["pos"]) - 1] = 1
ivs = []
for l in gzip.open(S + "OtherDifficult/GRCh38_alldifficultregions.bed.gz", "rt") if False else gzip.open(next(__import__("glob").iglob(S + "*/GRCh38_alldifficultregions.bed.gz")), "rt"):
    p = l.split("\t")
    if p[0] == "chr22": ivs.append((int(p[1]), int(p[2])))
ivs.sort(); starts = [s for s, _ in ivs]
def in_diff(pos):
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and pos < ivs[i][1]
used = {l.split()[1] for l in open("/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage3/frozen/R1.model") if l.startswith("feat")}
with open(D + "/features.tsv") as f:
    hdr = f.readline().rstrip("\n").split("\t"); X = []; pos = []
    for line in f:
        p = line.rstrip("\n").split("\t"); pos.append(int(p[0])); X.append([float(v) if v != "nan" else np.nan for v in p])
X = np.array(X); pos = np.array(pos)
y = np.array([truth.get(int(p), 0) for p in pos]); dif = np.array([in_diff(int(p)) for p in pos])
def auc(v, y):
    m = ~np.isnan(v); v, y = v[m], y[m]
    if y.sum() == 0 or y.sum() == len(y): return np.nan
    r = np.argsort(np.argsort(v, kind="stable"), kind="stable") + 1.0
    # average ranks for ties
    order = np.argsort(v, kind="stable"); sv = v[order]; ranks = np.empty(len(v)); i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and sv[j + 1] == sv[i]: j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1; i = j + 1
    n1 = y.sum(); n0 = len(y) - n1
    return (ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
print(f"candidates {len(y)}  true variants {int(y.sum())}  in alldifficult {int(dif.sum())} (true {int(y[dif].sum())}) ; outside {int((~dif).sum())} (true {int(y[~dif].sum())})")
print(f"{'feature':20s} {'inR1':5s} {'AUC_all':>8s} {'AUC_diff':>9s} {'AUC_easy':>9s}")
res = []
for j, h in enumerate(hdr):
    if j < 6 or h.startswith("pb") or h in ("label", "routed", "n_adm", "alt", "ref_idx", "pos0"): continue
    a, b, c = auc(X[:, j], y), auc(X[dif, j], y[dif]), auc(X[~dif, j], y[~dif])
    res.append((max(abs(b - .5) if b == b else 0, 0), h, a, b, c))
for _, h, a, b, c in sorted(res, reverse=True)[:30]:
    print(f"{h:20s} {'yes' if h in used else 'NO':5s} {a:8.3f} {b:9.3f} {c:9.3f}")
