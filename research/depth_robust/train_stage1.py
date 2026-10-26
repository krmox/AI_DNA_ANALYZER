"""Stage 1 fit (see PREREGISTRATION.md).  DEV blocks (even 2-Mb) of HG002 chr20 only, four depths pooled.
Writes R1.model / R2.model and FINALIST.txt.  Nothing here reads an odd block's label for fitting."""
import sys, json, hashlib
import numpy as np, pandas as pd, pysam
from scipy.stats import binom
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import GroupKFold

R = "/mnt/archive/AI_DNA_ANALYZER_depth_robust"
OUT = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/depth_robust"
B = "/mnt/archive/AI_DNA_ANALYZER_benchmark"
TRUTH = f"{B}/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
BED = f"{B}/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
BLOCK = 2_000_000
DEPTHS = ["f010", "f020", "f033", "f100"]
LOGF = {"depth", "n", "k", "nref", "k2", "alt_lowbq"}   # same set model.hpp log1p-transforms

A_SCALEFREE = ["alt_bq_mean", "ref_bq_mean", "bq_diff", "alt_bq_hi_frac", "alt_mq_mean", "ref_mq_mean", "mq_diff",
               "alt_mq_lt40_frac", "alt_fwd_frac", "ref_fwd_frac", "strand_diff", "alt_min_strand", "alt_end_mean",
               "ref_end_mean", "alt_near_end_frac", "alt_nm_mean", "ref_nm_mean", "alt_clip_frac", "ref_clip_frac",
               "alt_uniq_frac", "gap_frac", "ins_frac", "del_next_frac"]
R1 = ["vaf", "bllr_pr", "kcap", "lowbq_frac"] + A_SCALEFREE
R2 = R1 + ["depth"]
BASES = "ACGT"

truth = {}
for r in pysam.VariantFile(TRUTH).fetch("chr20"):
    alts = r.alts or ()
    if len(r.ref) == 1 and alts and all(len(a) == 1 for a in alts):
        gt = r.samples[0]["GT"]
        truth[r.pos - 1] = (r.ref, tuple(alts), "1/1" if (gt[0] == gt[1] and gt[0] != 0) else "0/1")
bed = np.array(sorted((int(a[1]), int(a[2])) for a in (l.split() for l in open(BED)) if a[0] == "chr20"))
inbed_mask = np.zeros(64444168, bool)
for a, b in bed: inbed_mask[a:b] = True
tpos = np.array(sorted(p for p in truth if inbed_mask[p]))
truth_total_dev = int(((tpos // BLOCK) % 2 == 0).sum())

P = np.clip(np.array([0.01, 0.5, 0.99]), 1e-12, 1 - 1e-12)
def methodc_gt(k, n):
    lp = np.stack([binom.logpmf(k, n, p) for p in P], 1)
    return np.where(lp.argmax(1) == 2, "1/1", "0/1")

frames = []
for d in DEPTHS:
    df = pd.read_csv(f"{R}/feat/features_{d}.tsv", sep="\t", float_precision="round_trip")
    df["depth_tag"] = d
    df["block"] = (df.pos0 // BLOCK).astype(int)
    df = df[(df.block % 2 == 0) & inbed_mask[df.pos0.values]].copy()          # DEV blocks, high-confidence BED
    gtc = methodc_gt(df.k.values, (df.k + df.nref).values)
    df["y"] = [int(p in truth and BASES[a] in truth[p][1] and truth[p][2] == g) for p, a, g in zip(df.pos0, df.alt, gtc)]
    df["ctrl"] = np.where(df.routed == 1, df.pb48 >= 10.5, df.bllr >= 7.0)
    frames.append(df)
dev = pd.concat(frames, ignore_index=True)

def design(df, feats):
    X = df[feats].copy()
    for c in feats:
        if c in LOGF: X[c] = np.log1p(X[c])
    return X.values

def f1_candidate_only(call, y):
    tp = int((call & (y == 1)).sum()); fp = int(call.sum()) - tp
    fn = truth_total_dev - tp          # truth total approximated by dev truth count (same for every depth)
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0

def mean_depth_f1(call):
    return float(np.mean([f1_candidate_only(call[(dev.depth_tag == d).values], dev.y.values[(dev.depth_tag == d).values]) for d in DEPTHS]))

def fit(feats):
    X, y, g = design(dev, feats), dev.y.values, dev.block.values
    best = None
    for C in (0.03, 0.1, 0.3, 1.0, 3.0):
        oof = np.zeros(len(y))
        for tr, te in GroupKFold(8).split(X, y, g):
            oof[te] = make_pipeline(StandardScaler(), LogisticRegression(C=C, max_iter=3000)).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
        for t in np.linspace(0.05, 0.95, 91):
            s = mean_depth_f1(oof >= t)
            if best is None or s > best[0]: best = (s, C, t, oof.copy())
    s, C, t, oof = best
    m = make_pipeline(StandardScaler(), LogisticRegression(C=C, max_iter=3000)).fit(X, y)
    return dict(score=s, C=C, thr=t, model=m, oof=oof)

def write_model(name, feats, fit_res):
    sc = fit_res["model"].named_steps["standardscaler"]; lr = fit_res["model"].named_steps["logisticregression"]
    t = fit_res["thr"]; lg = float(np.log(t / (1 - t)))
    lines = [f"name {name}", "kind logistic", "kmin 3", "routed_variant 1", f"intercept {float(lr.intercept_[0])!r}", f"zthr {lg!r}"]
    for j, f in enumerate(feats):
        lines.append(f"feat {f} {float(sc.mean_[j])!r} {float(sc.scale_[j])!r} {float(lr.coef_[0][j])!r}")
    path = f"{OUT}/models/{name}.model"
    open(path, "w").write("\n".join(lines) + "\n")
    return hashlib.sha256(open(path, "rb").read()).hexdigest()

import os; os.makedirs(f"{OUT}/models", exist_ok=True)
report = {"n_dev_rows": {d: int((dev.depth_tag == d).sum()) for d in DEPTHS}, "truth_total_dev": truth_total_dev}
ctrl_scores = {d: f1_candidate_only(dev.ctrl.values[(dev.depth_tag == d).values], dev.y.values[(dev.depth_tag == d).values]) for d in DEPTHS}
report["control_dev_f1_candidate_only"] = ctrl_scores
for name, feats in (("R1", R1), ("R2", R2)):
    res = fit(feats)
    per = {d: f1_candidate_only((res["oof"] >= res["thr"])[(dev.depth_tag == d).values], dev.y.values[(dev.depth_tag == d).values]) for d in DEPTHS}
    report[name] = dict(mean_dev_oof_f1=res["score"], C=res["C"], thr=res["thr"], per_depth_oof_f1=per, sha256=write_model(name, feats, res))
    print(name, json.dumps(report[name]))
finalist = "R1" if report["R1"]["mean_dev_oof_f1"] >= report["R2"]["mean_dev_oof_f1"] else "R2"
report["finalist"] = finalist
open(f"{OUT}/FINALIST.txt", "w").write(finalist + "\n")
json.dump(report, open(f"{OUT}/stage1_fit_report.json", "w"), indent=1)
print("control", ctrl_scores, "FINALIST", finalist)
