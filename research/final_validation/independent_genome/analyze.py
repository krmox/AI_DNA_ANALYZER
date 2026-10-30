#!/usr/bin/env python3
"""RESULTS.jsonl -> SUMMARY.md, BOOTSTRAP.tsv, ENDPOINTS.md, FAILED_RUNS.md (HG001 independent-genome campaign).

Bootstrap = preregistered stage3/PLAN.md + stage4/PREREG_R2.md method: 2-Mb blocks (floor((pos-1)/2e6)), 4,000 resamples, seed 20261027,
percentile 95% CI, SNP records only (hap.py BD/BVT), paired (same resampled blocks for both callers).
Aggregation over contigs x seeds (NOT specified numerically in the PREREG, documented here): in each replicate, blocks are resampled
independently PER CONTIG and the same block draw is applied to all seeds and callers of that contig; delta-F1 is computed per (contig, seed)
and averaged. Seeds share one BAM, so they are NOT independent samples; blocks within a contig are not independent either.
"""
import collections, json, statistics as st
from pathlib import Path
import numpy as np

H = Path(__file__).parent
R = [json.loads(l) for l in open(H / "RESULTS.jsonl")]
OK = [r for r in R if r.get("status") == "OK"]
RES, SEED = 4000, 20261027
CALLERS = ("r2", "r11", "bcftools")
NAME = {"r2": "R2", "r11": "R1.1", "bcftools": "bcftools"}


def f1c(tp, fp, fn):
    d = 2 * tp + fp + fn
    return 2 * tp / d if d else float("nan")


def arr(r, keys):
    d = {b: (tp, fp, fn) for b, tp, fp, fn in r["snp_blocks_2Mb"]}
    return np.array([d.get(k, (0, 0, 0)) for k in keys], dtype=float)


def f1_blocks(X, idx):  # X: nblocks x 3; idx: resample indices -> F1
    tp, fp, fn = X[idx].sum(axis=0)
    return 2 * tp / (2 * tp + fp + fn)


# group: (contig, target) -> seed -> caller -> rec
G = collections.defaultdict(lambda: collections.defaultdict(dict))
for r in OK:
    G[(r["contig"], r["target"])][r["seed"]][r["caller"]] = r
complete = {k: {s: c for s, c in v.items() if all(x in c for x in CALLERS)} for k, v in G.items()}
contigs = sorted({k[0] for k in complete}); depths = sorted({int(k[1]) for k in complete})

# ---------------- bootstrap
rng = np.random.default_rng(SEED)
boot_rows = []
agg = {}   # depth -> dict(delta_r11=(point, lo, hi), delta_bcf=...)
percell = {}
for t in depths:
    cells = []  # (contig, seed, X_r2, X_r11, X_bcf, blocks idx space)
    per_contig = {}
    for c in contigs:
        seeds = complete.get((c, str(t)), {})
        if not seeds:
            continue
        keys = sorted({b for s in seeds.values() for cal in CALLERS for b, *_ in s[cal]["snp_blocks_2Mb"]})
        per_contig[c] = (keys, {s: [arr(seeds[s][cal], keys) for cal in CALLERS] for s in sorted(seeds)})
    if not per_contig:
        continue
    # per-contig draws shared by all seeds/callers
    draws = {c: [rng.integers(0, len(per_contig[c][0]), len(per_contig[c][0])) for _ in range(RES)] for c in per_contig}
    # per (contig, seed): F1 for each caller point + CI; paired deltas
    d11 = np.zeros((RES, 0)); dbf = np.zeros((RES, 0)); P11 = []; PBF = []
    for c, (keys, seeds) in per_contig.items():
        for s, (X2, X11, XB) in seeds.items():
            idx_all = np.arange(len(keys))
            f2 = np.array([f1_blocks(X2, i) for i in draws[c]]); f11 = np.array([f1_blocks(X11, i) for i in draws[c]]); fb = np.array([f1_blocks(XB, i) for i in draws[c]])
            p2, p11, pb = f1_blocks(X2, idx_all), f1_blocks(X11, idx_all), f1_blocks(XB, idx_all)
            for lab, point, dist in (("R2", p2, f2), ("R1.1", p11, f11), ("bcftools", pb, fb), ("R2-R1.1", p2 - p11, f2 - f11), ("R2-bcftools", p2 - pb, f2 - fb)):
                lo, hi = np.percentile(dist, [2.5, 97.5])
                boot_rows.append(dict(contig=c, depth=t, seed=s, series=lab, quantity="F1", estimate=point, ci_lo=lo, ci_hi=hi, n_blocks=len(keys), resamples=RES, boot_seed=SEED))
            d11 = np.column_stack([d11, f2 - f11]); dbf = np.column_stack([dbf, f2 - fb]); P11.append(p2 - p11); PBF.append(p2 - pb)
    for lab, D, P in (("R2-R1.1", d11, P11), ("R2-bcftools", dbf, PBF)):
        m = D.mean(axis=1); lo, hi = np.percentile(m, [2.5, 97.5]); pt = float(np.mean(P))
        agg[(t, lab)] = (pt, float(lo), float(hi), D.shape[1], len(per_contig))
        boot_rows.append(dict(contig="ALL(mean over contigs x seeds)", depth=t, seed="all", series=lab, quantity="F1", estimate=pt, ci_lo=lo, ci_hi=hi,
                              n_blocks=sum(len(per_contig[c][0]) for c in per_contig), resamples=RES, boot_seed=SEED))
with open(H / "BOOTSTRAP.tsv", "w") as f:
    cols = ["contig", "depth", "seed", "series", "quantity", "estimate", "ci_lo", "ci_hi", "n_blocks", "resamples", "boot_seed"]
    f.write("\t".join(cols) + "\n")
    for r in boot_rows:
        f.write("\t".join(f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

# ---------------- consistency check TP/FP/FN vs reported P/R/F1
bad = []
for r in OK:
    tp, fn, fp, qtp = r["SNP_TP"], r["SNP_FN"], r["SNP_FP"], r["SNP_QTP"]
    chk = dict(precision=qtp / (qtp + fp), recall=tp / (tp + fn), F1=f1c(tp, fp, fn))
    for k, v in chk.items():
        if abs(v - r[f"SNP_{k}"]) > 2e-4:
            bad.append((r["cond"], r["caller"], k, v, r[f"SNP_{k}"]))
    s = sum(x[1] for x in r["snp_blocks_2Mb"])
    if abs(s - qtp) > max(5, 0.001 * qtp):
        bad.append((r["cond"], r["caller"], "blockTP_vs_QTP", s, qtp))

# ---------------- SUMMARY.md
L = ["# HG001 (NA12878) independent-genome results — hap.py 0.3.15, SNP PASS, confident-region BED", "",
     "Per-run rows: RESULTS.jsonl. All F1 values are point estimates from one contig and one downsampled BAM (seed); see BOOTSTRAP.tsv / ENDPOINTS.md for intervals.", ""]
L += ["## Per run", "", "| contig | depth (achieved) | seed | caller | TP | FP | FN | precision | recall | F1 | wall s | CPU s | peak RSS MB | exit |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in sorted(R, key=lambda r: (r["contig"], str(r["target"]).zfill(4), r.get("seed") or 0, CALLERS.index(r["caller"]))):
    if r.get("status") == "OK":
        L.append(f"| {r['contig']} | {r['target']}x ({r.get('achieved_bed_depth', 0):.1f}) | {r['seed']} | {NAME[r['caller']]} | {r['SNP_TP']} | {r['SNP_FP']} | {r['SNP_FN']} | "
                 f"{r['SNP_precision']:.5f} | {r['SNP_recall']:.5f} | {r['SNP_F1']:.5f} | {r['wall_s']:.1f} | {r['cpu_s']:.1f} | {r['peak_rss_kb']/1024:.0f} | {r['rc']} |")
    else:
        L.append(f"| {r['contig']} | {r['target']}x | {r.get('seed')} | {NAME.get(r['caller'], r['caller'])} | — | — | — | — | — | — | {r.get('wall_s', 0):.1f} | {r.get('cpu_s', 0):.1f} | {r.get('peak_rss_kb', 0)/1024:.0f} | {r.get('rc')} ({r.get('status')}) |")
L += ["", "## Mean over seeds (range = min–max over the 3 seeds of one BAM; NOT a CI)", "", "| contig | depth | R2 F1 | R1.1 F1 | bcftools F1 | R2−R1.1 | R2−bcftools | R2 wall/CPU/RSS | R1.1 wall/CPU/RSS | bcftools wall/CPU/RSS |", "|---|---|---|---|---|---|---|---|---|---|"]
for c in contigs:
    for t in depths:
        s = complete.get((c, str(t)))
        if not s:
            continue
        g = lambda cal, k: [x[cal][k] for x in s.values()]
        m = lambda v: f"{st.mean(v):.5f} ({min(v):.5f}–{max(v):.5f})" if len(v) > 1 else f"{v[0]:.5f}"
        res = lambda cal: f"{st.mean(g(cal,'wall_s')):.1f}s / {st.mean(g(cal,'cpu_s')):.0f}s / {max(g(cal,'peak_rss_kb'))/1024:.0f}MB"
        d1 = st.mean(a - b for a, b in zip(g("r2", "SNP_F1"), g("r11", "SNP_F1"))); d2 = st.mean(a - b for a, b in zip(g("r2", "SNP_F1"), g("bcftools", "SNP_F1")))
        L.append(f"| {c} | {t}x | {m(g('r2','SNP_F1'))} | {m(g('r11','SNP_F1'))} | {m(g('bcftools','SNP_F1'))} | {d1:+.5f} | {d2:+.5f} | {res('r2')} | {res('r11')} | {res('bcftools')} |")
L += ["", "## INDEL (reported separately; Varion emits SNPs only, so INDEL recall is not a Varion capability measure)", "", "| contig | depth | seed | caller | INDEL TP | FN | FP |", "|---|---|---|---|---|---|---|"]
for r in sorted(OK, key=lambda r: (r["contig"], r["target"], r["seed"], CALLERS.index(r["caller"]))):
    if r["seed"] == 101:
        L.append(f"| {r['contig']} | {r['target']}x | {r['seed']} | {NAME[r['caller']]} | {r.get('INDEL_TP')} | {r.get('INDEL_FN')} | {r.get('INDEL_FP')} |")
L += ["", f"Metric consistency check (precision=QTP/(QTP+FP), recall=TP/(TP+FN), F1=2TP/(2TP+FP+FN), sum of block TP≈QTP): {len(bad)} inconsistencies"] + [f"- {b}" for b in bad]
(H / "SUMMARY.md").write_text("\n".join(L) + "\n")

# ---------------- ENDPOINTS.md
E = ["# PREREG_R2 endpoints on HG001 (point estimates + preregistered block-bootstrap CIs)", "",
     "Thresholds are the preregistered ones, unchanged. Depths 15x and 50x of P2/P4 were NOT run here (matrix is 10/20/30x), so P2/P4 are evaluated on the tested depths only.",
     "'inconclusive' is a reporting convention added here (point estimate below threshold but CI upper bound reaches it); it is not preregistered.", ""]
a = agg.get((30, "R2-R1.1"))
if a:
    pt, lo, hi, n, nc = a
    verdict = "MET" if (pt >= 0.003 and lo > 0) else ("INCONCLUSIVE" if hi >= 0.003 else "NOT MET")
    E += [f"* **P1** (30x: mean ΔF1 R2−R1.1 ≥ +0.003 and CI lower bound > 0; mean over {nc} contig(s) × seeds, n={n} cells): ΔF1 = {pt:+.5f}, 95% block-bootstrap CI [{lo:+.5f}, {hi:+.5f}] → **{verdict}**"]
E += ["", "* **P2** (paired ΔF1 R2−bcftools, lower 95% bound ≥ −0.0005, per depth tested):"]
for t in depths:
    a = agg.get((t, "R2-bcftools"))
    if a:
        E.append(f"  * {t}x: ΔF1 = {a[0]:+.5f}, CI [{a[1]:+.5f}, {a[2]:+.5f}] → **{'MET' if a[1] >= -0.0005 else 'NOT MET'}**")
E += ["", "* **P4** (F1(R2) ≥ F1(R1.1) − 0.002 at 10x; 15x not run), per contig and seed mean:"]
for c in contigs:
    s = complete.get((c, "10"))
    if s:
        d = st.mean(x["r2"]["SNP_F1"] - x["r11"]["SNP_F1"] for x in s.values())
        E.append(f"  * {c}: ΔF1 (R2−R1.1) = {d:+.5f} → **{'MET' if d >= -0.002 else 'NOT MET'}**")
E += ["", "* **P3** (F1(R2) ≥ 0.990 at ≥30x; claimed only where bcftools also reaches it, else 'data-limited'), mean over seeds:"]
for c in contigs:
    s = complete.get((c, "30"))
    if s:
        f2 = st.mean(x["r2"]["SNP_F1"] for x in s.values()); fb = st.mean(x["bcftools"]["SNP_F1"] for x in s.values())
        E.append(f"  * {c} 30x: R2 {f2:.5f}, bcftools {fb:.5f} → {'R2 reaches 0.990' if f2 >= 0.990 else 'R2 below 0.990'}; bcftools {'also reaches' if fb >= 0.990 else 'does not reach'} it"
                 + (" (ceiling claim allowed)" if f2 >= 0.990 and fb >= 0.990 else " (target not claimable as a data ceiling)" if f2 >= 0.990 else " (data-limited if bcftools also below)" if fb < 0.990 else ""))
E += ["", "* **Guardrails vs R1.1** (wall ≤ +10 % or ≤ +2 s, CPU ≤ +10 %, RSS ≤ +15 %). Preregistered design was 5 interleaved reps; here runs are sequential single runs per seed, so this is indicative:"]
for c in contigs:
    for t in depths:
        s = complete.get((c, str(t)))
        if s:
            ratio = lambda k: st.mean(x["r2"][k] for x in s.values()) / st.mean(x["r11"][k] for x in s.values()) - 1
            wd = st.mean(x["r2"]["wall_s"] - x["r11"]["wall_s"] for x in s.values())
            ok = (ratio("wall_s") <= 0.10 or wd <= 2) and ratio("cpu_s") <= 0.10 and ratio("peak_rss_kb") <= 0.15
            E.append(f"  * {c} {t}x: wall {ratio('wall_s'):+.1%} ({wd:+.1f}s), CPU {ratio('cpu_s'):+.1%}, RSS {ratio('peak_rss_kb'):+.1%} → {'within guardrail' if ok else 'OUTSIDE guardrail'}")
(H / "ENDPOINTS.md").write_text("\n".join(E) + "\n")

# ---------------- FAILED_RUNS.md
F = ["# Failed or skipped runs", ""]
fails = [r for r in R if r.get("status") != "OK"]
F += [f"* {r['cond']} {r['caller']}: {r.get('status')} rc={r.get('rc')} cmd: `{r.get('command','')}` tail: {r.get('stderr_tail','')!r}" for r in fails] or ["No failed runs recorded in RESULTS.jsonl."]
F += ["", "Skipped by design (see REPORT.md): native 300x depth (impractical: ~3x the 100x cost for every caller and hap.py; the protocol targets 10–30x); 15x and 50x depths (PREREG P2/P4 depths) not in the requested matrix."]
(H / "FAILED_RUNS.md").write_text("\n".join(F) + "\n")
print("\n".join(E)); print("inconsistencies:", len(bad), "failures:", len(fails), "ok runs:", len(OK))
