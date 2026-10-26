#!/usr/bin/env python3
"""R2 offline ablation (preregistered arms, PREREG_R2.md). Reads feature dumps from R2 binary runs (--r2 --dump-features) and fits:
  A0  R1.model as is (control)            A0T R1.model + 3..4 parameter tail (gate fires on difficult loci)
  A1  R1's 29 features retrained          A2  A1 + sample-relative/context features
  A2R A2 with absolute MAPQ/BQ replaced by sample-relative ones (domain-robust variant)
  A3  A2 + three-scale pattern layer      A3R A3 with the robust substitution
Candidate-level F1 approximation: TP = true candidates called, FP = false candidates called, FN = truth SNPs in BED - TP
(ignores calls at k<3 and genotype errors; final numbers come from the real binary + hap.py).
Blocked CV (2-Mb blocks, all depths of a block in one fold) on TRAIN only; threshold = argmax of mean per-dump F1 on OOF scores.
VALID dumps are only scored, never used for fitting.
usage: fit_arms.py DUMPDIR OUTDIR
"""
import csv, glob, json, math, os, subprocess, sys
import numpy as np
from scipy.optimize import minimize
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold

DUMP, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
R1_MODEL = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage3/frozen/R1.model"
TRUTH = {"HG002": "/mnt/archive/AI_DNA_ANALYZER_benchmark/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
         "HG003": "/mnt/archive/AI_DNA_ANALYZER_final_eval/truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
         "HG004": "/mnt/archive/AI_DNA_ANALYZER_stage2/truth/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"}
LOG1P_R1 = {"depth", "n", "k", "nref", "k2", "alt_lowbq"}
# (name, log1p, cap) for the extra features
EXTRA = {"depth_rel": (0, 4.0), "nc25": (1, 10.0), "lmr25": (0, 1e300), "nearest": (0, 26.0), "depth_ratio25": (0, 4.0),
         "neigh_max_vaf25": (0, 1e300), "hp_len": (1, 30.0), "mq_rel": (0, 1e300), "bq_rel": (0, 1e300),
         "dens100": (1, 20.0), "dens5k": (1, 20.0), "lowvaf5k": (0, 1e300)}
ABS_TO_REL = ["alt_mq_mean", "ref_mq_mean", "mq_diff", "alt_bq_mean", "ref_bq_mean", "bq_diff"]
GATE = dict(depth_rel=1.25, dens5k=3.0, lowvaf5k=0.40)   # fixed a priori, not tuned
TAIL_FEATS = ["depth_rel", "dens5k", "lowvaf5k", "vaf"]


def parse_r1_model():
    m = {"feats": []}
    for line in open(R1_MODEL):
        p = line.split()
        if p[0] == "intercept": m["b"] = float(p[1])
        elif p[0] == "zthr": m["thr"] = float(p[1])
        elif p[0] == "feat": m["feats"].append((p[1], float(p[2]), float(p[3]), float(p[4])))
    return m


R1 = parse_r1_model()
R1_NAMES = [f[0] for f in R1["feats"]]


def load_dump(path):
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f]
    ix = {h: i for i, h in enumerate(hdr)}
    cols = {h: np.array([float(r[ix[h]]) if r[ix[h]] != "nan" else np.nan for r in rows]) for h in hdr if not h.startswith("pb")}
    return cols


def truth_snps(sample, bed):
    key = f"{OUT}/truth_{sample}_{os.path.basename(os.path.dirname(bed))}.txt"
    if os.path.exists(key):
        return int(open(key).read())
    txt = subprocess.run(["bcftools", "view", "-H", "-v", "snps", "-R", bed, TRUTH[sample]], capture_output=True, text=True).stdout
    n = sum(1 for l in txt.splitlines() if l)
    open(key, "w").write(str(n))
    return n


sets = []
for d in sorted(glob.glob(f"{DUMP}/*/features.tsv")):
    name = os.path.basename(os.path.dirname(d))
    role = name[:2]
    sample = {"hg002": "HG002", "hg003": "HG003", "hg004": "HG004"}[name.split("_")[1]]
    contig = name.split("_")[2]
    cols = load_dump(d)
    bed = os.path.join(os.path.dirname(d), "conf.bed")
    meta = json.load(open(os.path.join(os.path.dirname(d), "stdout.json")))
    sets.append(dict(name=name, role=role, sample=sample, contig=contig, cols=cols, y=(cols["label"] == 1).astype(int),
                     truth=truth_snps(sample, bed), cal=meta.get("r2", {}), pos=cols["pos0"],
                     group=np.array([f"{sample}_{contig}_{int(p) // 2_000_000}" for p in cols["pos0"]])))
print("dumps:", [(s["name"], len(s["y"]), s["truth"]) for s in sets])


def xform(name, v):
    if name in EXTRA:
        lg, cap = EXTRA[name]
        v = np.minimum(v, cap)
        return np.log1p(v) if lg else v
    return np.log1p(v) if name in LOG1P_R1 else v


def design(s, names):
    return np.column_stack([xform(n, s["cols"][n]) for n in names])


def r1_score(s):
    z = np.full(len(s["y"]), R1["b"])
    for n, mean, scale, coef in R1["feats"]:
        z += coef * (xform(n, s["cols"][n]) - mean) / scale
    return z


def f1_at(s, score, thr):
    called = score >= thr
    tp = int((called & (s["y"] == 1)).sum()); fp = int((called & (s["y"] == 0)).sum()); fn = s["truth"] - tp
    return 2 * tp / (2 * tp + fp + fn), tp, fp, fn


def best_thr(sc_list, y_sets, truths, grid):
    best = (-1, None)
    for t in grid:
        f = []
        for sc, s in zip(sc_list, y_sets):
            f.append(f1_at(s, sc, t)[0])
        m = float(np.mean(f))
        if m > best[0]: best = (m, t)
    return best[1], best[0]


train = [s for s in sets if s["role"] == "tr"]
valid = [s for s in sets if s["role"] == "va"]
rows_out = []


def report(arm, thr, scores_valid, scores_train=None):
    line = {"arm": arm, "thr": round(float(thr), 4)}
    fv = []
    for s, sc in zip(valid, scores_valid):
        f, tp, fp, fn = f1_at(s, sc, thr)
        line[s["name"]] = round(f, 5); fv.append(f)
    line["valid_mean_F1"] = round(float(np.mean(fv)), 5)
    if scores_train is not None:
        line["train_oof_mean_F1"] = round(float(np.mean([f1_at(s, sc, thr)[0] for s, sc in zip(train, scores_train)])), 5)
    rows_out.append(line)
    print(line)


grid = np.arange(-3, 6.01, 0.1)
# ---- A0: R1 as is
a0_tr = [r1_score(s) for s in train]; a0_va = [r1_score(s) for s in valid]
report("A0 R1.model (frozen thr)", R1["thr"], a0_va, a0_tr)

# ---- A0T: frozen z + tail
def gate(s):
    c = s["cols"]
    return ((c["depth_rel"] > GATE["depth_rel"]) | (c["dens5k"] > GATE["dens5k"]) | (c["lowvaf5k"] > GATE["lowvaf5k"])).astype(float)

def tail_design(s):
    return np.column_stack([xform(n, s["cols"][n]) for n in TAIL_FEATS])

Xt = np.vstack([tail_design(s) for s in train]); mu = Xt.mean(0); sd = Xt.std(0) + 1e-9
Gt = np.concatenate([gate(s) for s in train]); Z0 = np.concatenate(a0_tr); Yt = np.concatenate([s["y"] for s in train])
def nll(w):
    adj = Gt * (((Xt - mu) / sd) @ w[:-1] + w[-1])
    z = Z0 + adj
    return float(np.mean(np.logaddexp(0, z) - Yt * z)) + 1e-3 * float(np.sum(w[:-1] ** 2))
w = minimize(nll, np.zeros(len(TAIL_FEATS) + 1), method="L-BFGS-B").x
print("tail weights", dict(zip(TAIL_FEATS + ["bias"], np.round(w, 3))), "gated fraction", round(float(Gt.mean()), 3))
def a0t_score(s, a0):
    return a0 + gate(s) * (((tail_design(s) - mu) / sd) @ w[:-1] + w[-1])
a0t_tr = [a0t_score(s, z) for s, z in zip(train, a0_tr)]; a0t_va = [a0t_score(s, z) for s, z in zip(valid, a0_va)]
report("A0T frozen z + tail (frozen thr)", R1["thr"], a0t_va, a0t_tr)
t0t, _ = best_thr(a0t_tr, train, None, grid)
report("A0T frozen z + tail (thr re-fit on train)", t0t, a0t_va, a0t_tr)

# ---- retrained arms
def fit_arm(arm, names):
    X = np.vstack([design(s, names) for s in train]); y = np.concatenate([s["y"] for s in train]); g = np.concatenate([s["group"] for s in train])
    mean = X.mean(0); scale = X.std(0) + 1e-9
    Xs = (X - mean) / scale
    best = None
    for C in (0.03, 0.1, 0.3, 1.0, 3.0):
        oof = np.zeros(len(y))
        for tr, te in GroupKFold(n_splits=5).split(Xs, y, g):
            m = LogisticRegression(C=C, max_iter=3000).fit(Xs[tr], y[tr])
            oof[te] = m.decision_function(Xs[te])
        # per-dump OOF F1 at best common threshold
        parts = np.cumsum([0] + [len(s["y"]) for s in train])
        oofs = [oof[parts[i]:parts[i + 1]] for i in range(len(train))]
        t, f = best_thr(oofs, train, None, grid)
        if best is None or f > best[0]: best = (f, C, t, oofs)
    f, C, t, oofs = best
    m = LogisticRegression(C=C, max_iter=3000).fit(Xs, y)
    sc_va = [m.decision_function((design(s, names) - mean) / scale) for s in valid]
    report(f"{arm} (C={C})", t, sc_va, oofs)
    return dict(arm=arm, names=names, mean=mean, scale=scale, coef=m.coef_[0], intercept=float(m.intercept_[0]), thr=float(t), C=C)

A1 = list(R1_NAMES)
A2 = A1 + ["depth_rel", "nc25", "lmr25", "nearest", "depth_ratio25", "neigh_max_vaf25", "hp_len", "mq_rel", "bq_rel"]
A2R = [n for n in A2 if n not in ABS_TO_REL]
A3 = A2 + ["dens100", "dens5k", "lowvaf5k"]
A3R = [n for n in A3 if n not in ABS_TO_REL]
fits = {a: fit_arm(a, n) for a, n in (("A1", A1), ("A2", A2), ("A2R", A2R), ("A3", A3), ("A3R", A3R))}

with open(f"{OUT}/arms_results.csv", "w", newline="") as f:
    cols = sorted({k for r in rows_out for k in r}, key=lambda k: (k not in ("arm", "thr"), k))
    wr = csv.DictWriter(f, cols); wr.writeheader(); wr.writerows(rows_out)

# ---- export model files
def envelope():
    eps, mq, bq, dp = [], [], [], []
    for s in sets:
        c = s["cal"]
        if not c.get("calibration_valid"): continue
        eps.append(c["eps_hat"]); mq.append(c["mq_p90"]); bq.append(c["bq_p50"]); dp.append(c["depth_med"])
    return dict(eps=(0.5 * min(eps), 3.0 * max(eps)), mq=(min(mq) - 10, max(mq) + 10), bq=(min(bq) - 10, max(bq) + 10), depth=(4.0, 800.0))
ENV = envelope()
print("OOD envelope", ENV)

def write_model(path, name, thr, intercept, feats_lines, tail_lines=(), tail_bias=0.0):
    with open(path, "w") as f:
        f.write(f"name {name}\nkind logistic\nkmin 3\nrouted_variant 1\nintercept {intercept!r}\nzthr {thr!r}\n")
        for l in feats_lines: f.write(l + "\n")
        for l in tail_lines: f.write(l + "\n")
        if tail_lines: f.write(f"tail_bias {tail_bias!r}\n")
        f.write(f"gate_depth_rel {GATE['depth_rel']}\ngate_dens5k {GATE['dens5k']}\ngate_lowvaf5k {GATE['lowvaf5k']}\n")
        f.write(f"ood_eps {ENV['eps'][0]:.6g} {ENV['eps'][1]:.6g}\nood_mq {ENV['mq'][0]:.4g} {ENV['mq'][1]:.4g}\n")
        f.write(f"ood_bq {ENV['bq'][0]:.4g} {ENV['bq'][1]:.4g}\nood_depth {ENV['depth'][0]} {ENV['depth'][1]}\n")

def feat_lines(fit):
    out = []
    for n, mean, scale, coef in zip(fit["names"], fit["mean"], fit["scale"], fit["coef"]):
        lg, cap = (EXTRA[n] if n in EXTRA else ((1 if n in LOG1P_R1 else 0), 1e300))
        out.append(f"feat2 {n} {float(mean)!r} {float(scale)!r} {float(coef)!r} {lg} {'inf' if cap >= 1e299 else cap}")
    return out

os.makedirs(f"{OUT}/models", exist_ok=True)
r1_lines = [f"feat2 {n} {mean!r} {scale!r} {coef!r} {1 if n in LOG1P_R1 else 0} inf" for n, mean, scale, coef in R1["feats"]]
a0t_lines = []
for n, mi, si, wi in zip(TAIL_FEATS, mu, sd, w[:-1]):
    lg, cap = EXTRA.get(n, (0, 1e300))
    a0t_lines.append(f"tail_feat2 {n} {float(mi)!r} {float(si)!r} {float(wi)!r} {lg} {'inf' if cap >= 1e299 else cap}")
write_model(f"{OUT}/models/A0T.model", "A0T", R1["thr"], R1["b"], r1_lines, a0t_lines, float(w[-1]))
write_model(f"{OUT}/models/A0.model", "A0", R1["thr"], R1["b"], r1_lines)
for a, fit in fits.items():
    write_model(f"{OUT}/models/{a}.model", a, fit["thr"], fit["intercept"], feat_lines(fit))
print("wrote", OUT)
