#!/usr/bin/env python3
"""Diagnostic (read-only) error anatomy of one Varion R1.1 run. Does not alter any prediction.

usage: diag_errors.py DIAGDIR BAM
DIAGDIR contains records.tsv (hap.py per-record), loci_dump.tsv (--dump-loci at FN/FP loci), features.tsv (--dump-features).
Writes DIAGDIR/summary.txt and DIAGDIR/fn_fp_annotated.tsv.
"""
import collections, csv, math, sys
from pathlib import Path
import numpy as np
import pysam

D = Path(sys.argv[1]); BAM = sys.argv[2]
MODEL = Path("/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage3/frozen/R1.model")
BASE = "ACGT"
out = open(D / "summary.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.write(s + "\n")

recs = {}
for r in csv.DictReader(open(D / "records.tsv"), delimiter="\t"):
    recs[int(r["pos"]) - 1] = r
dump = {int(r["pos0"]): r for r in csv.DictReader(open(D / "loci_dump.tsv"), delimiter="\t")}

# ---------- raw read evidence at each FN/FP locus (no filters) ----------
bam = pysam.AlignmentFile(BAM)
def raw_evidence(pos0, alt_char):
    ev = collections.Counter()
    for col in bam.pileup("chr22", pos0, pos0 + 1, truncate=True, stepper="nofilter", min_base_quality=0, ignore_overlaps=False,
                          ignore_orphans=False, max_depth=100000, min_mapping_quality=0, fastaFile=None):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            a = pr.alignment
            base = a.query_sequence[pr.query_position]
            ev["raw_depth"] += 1
            if base != alt_char:
                continue
            ev["raw_alt"] += 1
            reason = None
            if a.is_unmapped: reason = "unmapped"
            elif a.is_duplicate: reason = "dup"
            elif a.is_qcfail: reason = "qcfail"
            elif a.is_secondary or a.is_supplementary: reason = "sec/supp"
            elif a.is_paired and not a.is_proper_pair: reason = "not_proper_pair"
            elif a.mapping_quality < 20: reason = "mapq<20"
            elif a.query_qualities[pr.query_position] < 13: reason = "bq<13"
            ev["alt_kept" if reason is None else "alt_drop_" + reason] += 1
    return ev

rows = []
for pos0, r in sorted(recs.items()):
    is_fn = r["truth_bd"] == "FN"
    is_fp = r["query_bd"] == "FP"
    if not (is_fn or is_fp) or pos0 not in dump:
        continue
    d = dump[pos0]
    c = [int(d[f"c{i}"]) for i in range(4)]
    depth = int(d["c6"]); ref_idx = int(d["ref_idx"])
    truth_alt = r["alt"].split(",")[0] if is_fn else None
    q_alt = r["alt"].split(",")[0]
    alt_char = truth_alt if is_fn else q_alt
    k_truth_alt = c["ACGT".index(alt_char)] if alt_char in BASE else -1
    ev = raw_evidence(pos0, alt_char) if alt_char in BASE else collections.Counter()
    rows.append(dict(pos0=pos0, kind=("FN+FP" if is_fn and is_fp else "FN" if is_fn else "FP"), ref=r["ref"], alt=alt_char,
                     truth_gt=r["truth_gt"], query_gt=r["query_gt"], truth_bk=r["truth_bk"], query_bk=r["query_bk"], frame=int(d["frame"]),
                     depth_f=depth, counts="/".join(map(str, c)), k_alt_f=k_truth_alt, k_best=int(d["k"]), best_alt=BASE[int(d["alt"])],
                     vaf_alt=round(k_truth_alt / max(1, sum(c)), 3), bllr=float.fromhex(d["binomial_llr_hex"]), called=int(d["cascade_call"]),
                     routed=int(d["routed"]), regions=r["regions"], **{k: ev.get(k, 0) for k in ("raw_depth", "raw_alt", "alt_kept")},
                     drops=";".join(f"{k[9:]}={v}" for k, v in sorted(ev.items()) if k.startswith("alt_drop_"))))

cols = list(rows[0].keys())
with open(D / "fn_fp_annotated.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, cols, delimiter="\t"); w.writeheader(); w.writerows(rows)

# ---------- model score on every candidate (k>=3) ----------
feat = {}
model = {"feats": []}
for line in open(MODEL):
    p = line.split()
    if p[0] == "intercept": model["b"] = float(p[1])
    elif p[0] == "zthr": model["thr"] = float(p[1])
    elif p[0] == "feat": model["feats"].append((p[1], float(p[2]), float(p[3]), float(p[4])))
LOG1P = {"depth", "n", "k", "nref", "k2", "alt_lowbq"}
cand = {}
with open(D / "features.tsv") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {h: i for i, h in enumerate(hdr)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        z = model["b"]
        for name, mean, scale, coef in model["feats"]:
            x = float(p[ix[name]]);
            if name in LOG1P: x = math.log1p(x)
            z += coef * (x - mean) / scale
        cand[int(p[0])] = dict(z=z, k=float(p[ix["k"]]), vaf=float(p[ix["vaf"]]), depth=float(p[ix["depth"]]), routed=int(p[ix["routed"]]),
                               alt_mq=float(p[ix["alt_mq_mean"]]), clip=float(p[ix["alt_clip_frac"]]), nm=float(p[ix["alt_nm_mean"]]),
                               bq=float(p[ix["alt_bq_mean"]]), strand=float(p[ix["strand_diff"]]), hp=float(p[ix["hp_len"]]))

say("== RUN:", D)
say("hap.py SNP totals: TP", sum(r["truth_bd"] == "TP" for r in recs.values()), "FN", sum(r["truth_bd"] == "FN" for r in recs.values()),
    "FP", sum(r["query_bd"] == "FP" for r in recs.values()))
say("candidates (k>=3, in frame):", len(cand), " model zthr:", model["thr"])

# ---------- FN anatomy ----------
fn = [r for r in rows if r["kind"] in ("FN", "FN+FP")]
def fn_cat(r):
    if r["kind"] == "FN+FP": return "E genotype/allele mismatch (called, wrong GT/allele)"
    if not r["frame"]: return "A outside call frame (BED edge)"
    if r["k_alt_f"] == 0:
        if r["raw_alt"] == 0: return "B1 no ALT read in BAM at all (data-limited)"
        return "B2 ALT reads exist but all removed by read filters"
    if r["k_alt_f"] < 3: return "C 1-2 ALT reads after filters (below candidate rule k>=3; binomial<7)"
    if r["k_best"] > 0 and r["best_alt"] != r["alt"]: return "D2 k>=3 but best non-ref allele is a different base"
    return "D1 k>=3 candidate rejected by model (z<thr)"
cat = collections.Counter(fn_cat(r) for r in fn)
say("\n== FN anatomy (n=%d)" % len(fn))
for k, v in sorted(cat.items()): say(f"  {k:75s} {v:5d}  {100*v/len(fn):5.1f}%")
b2 = collections.Counter()
for r in fn:
    if fn_cat(r).startswith("B2"):
        for item in r["drops"].split(";"):
            if item: b2[item.split("=")[0]] += int(item.split("=")[1])
say("  B2 reasons (ALT reads dropped, summed over loci):", dict(b2))
dep = [r["depth_f"] for r in fn if r["frame"]]
say("  FN filtered depth: median %.0f, p10 %.0f, p90 %.0f ; (genome median ~30)" % (np.median(dep), np.percentile(dep, 10), np.percentile(dep, 90)))
d1 = [r for r in fn if fn_cat(r).startswith("D1")]
if d1:
    zs = [cand[r["pos0"]]["z"] for r in d1 if r["pos0"] in cand]
    say("  D1 z-scores (thr %.3f): median %.2f, min %.2f, max %.2f; within 1.0 below thr: %d / %d" % (
        model["thr"], np.median(zs), min(zs), max(zs), sum(model["thr"] - 1.0 <= z < model["thr"] for z in zs), len(zs)))
    say("  D1 mean alt MAPQ %.1f, alt clip frac %.2f, alt NM %.2f, vaf %.2f" % tuple(
        np.mean([cand[r["pos0"]][k] for r in d1 if r["pos0"] in cand]) for k in ("alt_mq", "clip", "nm", "vaf")))
rg = collections.Counter()
for r in fn:
    rg["in alldifficultregions" if "alldifficultregions" in r["regions"].split(",") and "notinalldifficult" not in r["regions"] else "outside difficult"] += 1
say("  FN in/out of GIAB alldifficultregions:", dict(rg))

# ---------- FP anatomy ----------
fp = [r for r in rows if r["kind"] in ("FP", "FN+FP")]
say("\n== FP anatomy (n=%d; %d are genotype/allele mismatches at a true variant)" % (len(fp), sum(r["kind"] == "FN+FP" for r in fp)))
pure = [r for r in fp if r["kind"] == "FP"]
via_model = [r for r in pure if r["pos0"] in cand]
say("  pure FP called via model (k>=3 candidate): %d ; via Binomial rule (k<3): %d" % (len(via_model), len(pure) - len(via_model)))
if via_model:
    for name, key in (("VAF", "vaf"), ("alt MAPQ mean", "alt_mq"), ("alt clip frac", "clip"), ("alt NM mean", "nm"), ("alt BQ mean", "bq"), ("strand diff", "strand"), ("z", "z"), ("homopolymer len", "hp")):
        say(f"    {name:16s} median {np.median([cand[r['pos0']][key] for r in via_model]):.3f}")
    vafs = np.array([cand[r["pos0"]]["vaf"] for r in via_model])
    say("    VAF<0.2: %d ; 0.2-0.35: %d ; 0.35-0.65: %d ; >0.65: %d" % ((vafs < .2).sum(), ((vafs >= .2) & (vafs < .35)).sum(), ((vafs >= .35) & (vafs < .65)).sum(), (vafs >= .65).sum()))
    zs = np.array([cand[r["pos0"]]["z"] for r in via_model])
    say("    z of FPs: within 1.0 above thr: %d / %d" % (((zs >= model["thr"]) & (zs < model["thr"] + 1)).sum(), len(zs)))
rgp = collections.Counter()
for r in pure:
    rgp["in alldifficultregions" if "alldifficultregions" in r["regions"].split(",") and "notinalldifficult" not in r["regions"] else "outside difficult"] += 1
say("  pure FP in/out of GIAB alldifficultregions:", dict(rgp))
ra = [r["raw_alt"] for r in pure]; kf = [r["k_alt_f"] for r in pure]
say("  pure FP: raw ALT reads median %.0f; filtered ALT k median %.0f; filtered depth median %.0f" % (np.median(ra), np.median(kf), np.median([r["depth_f"] for r in pure])))

# ---------- oracle z-threshold sweep on candidates (diagnostic only; NOT applied) ----------
tp0 = sum(r["truth_bd"] == "TP" for r in recs.values()); fn0 = sum(r["truth_bd"] == "FN" for r in recs.values()); fp0 = sum(r["query_bd"] == "FP" for r in recs.values())
def contribution(pos0, called):
    r = recs.get(pos0)
    if r is None: return (0, 1, 0) if called else (0, 0, 0)       # no truth: FP if called
    truth = r["truth_bd"] in ("TP", "FN")
    if truth and called: return (1, 0, 0) if (r["truth_bd"] == "TP") else (0, 1, 1)   # mismatch: FP and FN
    if truth and not called: return (0, 0, 1)
    return (0, 1, 0) if called else (0, 0, 0)
cur = np.zeros(3, int); thr = model["thr"]
cur_called = {p: c["z"] >= thr for p, c in cand.items()}
for p, called in cur_called.items(): cur += contribution(p, called)
fixed = np.array([tp0, fp0, fn0]) - cur
say("\n== ORACLE sweep of z threshold on the %d candidates (diagnostic; uses truth; shows ceiling of recalibrating only that threshold)" % len(cand))
zs_sorted = sorted(set(round(c["z"], 2) for c in cand.values()))
best = (0, None)
res = []
zarr = np.array([c["z"] for c in cand.values()]); pos = list(cand.keys())
contrib_called = np.array([contribution(p, True) for p in pos]); contrib_not = np.array([contribution(p, False) for p in pos])
for t in np.arange(-3, 6.01, 0.25):
    m = zarr >= t
    tot = fixed + contrib_called[m].sum(0) + contrib_not[~m].sum(0)   # (tp, fp, fn) with convention: columns tp, fp, fn
    f1 = 2 * tot[0] / (2 * tot[0] + tot[1] + tot[2])
    res.append((t, tot[0], tot[1], tot[2], f1))
    if f1 > best[0]: best = (f1, t)
for t, a, b, c_, f in res:
    if abs(t - thr) < 0.13 or abs(t - best[1]) < 1e-9 or t in (-1, 0, 1, 2):
        say(f"  z>={t:5.2f}: TP {int(a)} FP {int(b)} FN {int(c_)} F1 {f:.5f}" + ("   <- frozen threshold" if abs(t - thr) < 0.13 else "") + ("   <- oracle best" if abs(t - best[1]) < 1e-9 else ""))
say(f"  frozen-threshold F1 (this run) = {2*tp0/(2*tp0+fp0+fn0):.5f}; oracle-best single-threshold F1 = {best[0]:.5f} at z>={best[1]}")
out.close()
