#!/usr/bin/env python3
"""Expanded benchmark analysis (HG001 + HG007, both PREVIOUSLY EVALUATED; not an independent-genome test).
Implements expanded_benchmark_preregistration.md exactly: hap.py SNP PASS, 2-Mb block bootstrap (4000 resamples, seed 20261027, percentile 95%),
resampling per cell (sample x contig) with one block draw shared by all seeds and callers, delta averaged over seeds then over cells.
Outputs (evidence_package/): expanded_results.csv, expanded_runtime_memory.csv, BOOTSTRAP_EXPANDED.tsv, ENDPOINTS_EXPANDED.md, DEEPVARIANT_TABLE.md
Run with ~/varion_mamba/envs/vt/bin/python (numpy).
"""
import collections, csv, json, statistics as st
from pathlib import Path
import numpy as np

F = Path("/home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER/research/final_validation")
OUT = F / "evidence_package"
RES, SEED = 4000, 20261027
CELLS = [("HG001", "chr17"), ("HG001", "chr19"), ("HG007", "chr17"), ("HG007", "chr19")]
DEPTHS = [10, 15, 20, 30, 50]
NAME = {"r2": "R2", "r11": "R1.1", "bcftools": "bcftools 1.21", "deepvariant": "DeepVariant 1.10.0"}
MARGIN = -0.0005


def load():
    rows = []
    for l in open(F / "independent_genome/RESULTS.jsonl"):
        r = json.loads(l); r["prov"] = "EXISTING"; rows.append(r)
    blk = {(b["cond"], b["caller"]): b["snp_blocks_2Mb"] for b in map(json.loads, open(F / "expanded/hg007_existing_blocks.jsonl"))}
    for l in open(F / "hg007/RESULTS.jsonl"):
        r = json.loads(l); r["prov"] = "EXISTING"; r["snp_blocks_2Mb"] = blk[(r["cond"], r["caller"])]; rows.append(r)
    p = F / "expanded/RESULTS.jsonl"
    if p.exists():
        for l in open(p):
            r = json.loads(l); r["prov"] = "NEW"; rows.append(r)
    return rows


def f1c(tp, fp, fn):
    d = 2 * tp + fp + fn
    return 2 * tp / d if d else float("nan")


def arr(r, keys):
    d = {b: (tp, fp, fn) for b, tp, fp, fn in r["snp_blocks_2Mb"]}
    return np.array([d.get(k, (0, 0, 0)) for k in keys], dtype=float)


def f1_draws(X, idx):  # X nb x 3 ; idx RES x nb -> F1 per draw
    s = X[idx].sum(axis=1)
    return 2 * s[:, 0] / (2 * s[:, 0] + s[:, 1] + s[:, 2])


def f1_point(X):
    s = X.sum(axis=0)
    return 2 * s[0] / (2 * s[0] + s[1] + s[2])


rows = load()
OK = [r for r in rows if r.get("status") == "OK"]
G = collections.defaultdict(dict)  # (sample,contig,depth) -> {(seed): {caller: rec}}
for r in OK:
    if r["target"] == "native":
        continue
    G[(r["sample"], r["contig"], int(r["target"]))].setdefault(r["seed"], {})[r["caller"]] = r


def seeds_complete(cell, d, callers):
    return {s: c for s, c in G.get((cell[0], cell[1], d), {}).items() if all(x in c for x in callers)}


def boot(pairs, callers_a, callers_b, depth, rng, cells):
    """pairs: list of (a, b) caller keys to compare (delta a-b). Returns dict level -> (point, lo, hi, ncell) or None if UNASSESSED, per cell draws."""
    out = {}
    cell_delta = {}
    for cell in cells:
        sc = seeds_complete(cell, depth, set([callers_a, callers_b]))
        if not sc:
            cell_delta[cell] = None; continue
        keys = sorted({b for s in sc.values() for cal in (callers_a, callers_b) for b, *_ in s[cal]["snp_blocks_2Mb"]})
        idx = rng.integers(0, len(keys), (RES, len(keys)))
        dd = np.zeros(RES); pt = 0.0
        for s, c in sc.items():
            Xa, Xb = arr(c[callers_a], keys), arr(c[callers_b], keys)
            dd += (f1_draws(Xa, idx) - f1_draws(Xb, idx)) / len(sc); pt += (f1_point(Xa) - f1_point(Xb)) / len(sc)
        cell_delta[cell] = (pt, dd, len(sc))
    for lvl, cs in (("HG001", [c for c in cells if c[0] == "HG001"]), ("HG007", [c for c in cells if c[0] == "HG007"]), ("POOLED(4 cells)", list(cells))):
        if not cs or any(cell_delta[c] is None for c in cs):
            out[lvl] = None; continue
        m = np.mean([cell_delta[c][1] for c in cs], axis=0); pt = float(np.mean([cell_delta[c][0] for c in cs]))
        lo, hi = np.percentile(m, [2.5, 97.5])
        out[lvl] = (pt, float(lo), float(hi), len(cs), sum(cell_delta[c][2] for c in cs))
    return out


# ---------------- bootstraps
rng = np.random.default_rng(SEED)
boot_rows = []; B = {}
for d in DEPTHS:
    for a, b in (("r2", "r11"), ("r2", "bcftools")):
        res = boot(None, a, b, d, rng, CELLS)
        for lvl, v in res.items():
            B[(d, a, b, lvl)] = v
            boot_rows.append(dict(depth=d, comparison=f"{NAME[a]} - {NAME[b]}", level=lvl, state="UNASSESSED" if v is None else "ok",
                                  delta_F1=None if v is None else v[0], ci_lo=None if v is None else v[1], ci_hi=None if v is None else v[2],
                                  n_cells=None if v is None else v[3], n_cell_seeds=None if v is None else v[4], resamples=RES, boot_seed=SEED))
with open(OUT / "BOOTSTRAP_EXPANDED.tsv", "w") as f:
    cols = list(boot_rows[0]); f.write("\t".join(cols) + "\n")
    for r in boot_rows:
        f.write("\t".join("" if r[c] is None else (f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c])) for c in cols) + "\n")

# ---------------- results.csv / runtime csv
bad = []
with open(OUT / "expanded_results.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["provenance", "sample", "contig", "depth", "seed", "caller", "status", "achieved_bed_depth", "SNP_TP", "SNP_FP", "SNP_FN", "SNP_precision", "SNP_recall", "SNP_F1_happy", "SNP_F1_from_counts", "INDEL_TP", "INDEL_FP", "INDEL_FN", "INDEL_F1", "bam_sha256", "vcf_sha256"])
    for r in sorted(rows, key=lambda r: (r["prov"], r["sample"], r["contig"], str(r["target"]).zfill(6), r.get("seed") or 0, r["caller"])):
        if r.get("status") != "OK":
            w.writerow([r["prov"], r["sample"], r["contig"], r["target"], r.get("seed"), NAME.get(r["caller"], r["caller"]), r.get("status")]); continue
        tp, fp, fn = r["SNP_TP"], r["SNP_FP"], r["SNP_FN"]; fc = f1c(tp, fp, fn)
        if abs(fc - r["SNP_F1"]) > 2e-4:
            bad.append((r["cond"], r["caller"], fc, r["SNP_F1"]))
        w.writerow([r["prov"], r["sample"], r["contig"], r["target"], r.get("seed"), NAME.get(r["caller"], r["caller"]), "OK", f"{r.get('achieved_bed_depth', 0):.2f}" if r.get("achieved_bed_depth") else "",
                    tp, fp, fn, r["SNP_precision"], r["SNP_recall"], r["SNP_F1"], f"{fc:.6f}", r.get("INDEL_TP"), r.get("INDEL_FP"), r.get("INDEL_FN"), r.get("INDEL_F1"), r.get("bam_sha256", ""), r.get("vcf_sha256", "")])
with open(OUT / "expanded_runtime_memory.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["provenance", "sample", "contig", "depth", "seed", "caller", "wall_s", "cpu_s", "peak_rss_mb", "loadavg_before", "timing_note"])
    for r in sorted(rows, key=lambda r: (r["prov"], r["sample"], r["contig"], str(r["target"]).zfill(6), r.get("seed") or 0, r["caller"])):
        if "wall_s" in r:
            note = "EXISTING: no loadavg recorded; HG001 chr19 overlapped a download" if r["prov"] == "EXISTING" and r["sample"] == "HG001" and r["contig"] == "chr19" else ("EXISTING: no loadavg recorded" if r["prov"] == "EXISTING" else "")
            w.writerow([r["prov"], r["sample"], r["contig"], r["target"], r.get("seed"), NAME.get(r["caller"], r["caller"]), f"{r['wall_s']:.2f}", f"{r['cpu_s']:.2f}", f"{r['peak_rss_kb']/1024:.0f}", f"{r['loadavg_before']:.2f}" if "loadavg_before" in r else "", note])

# ---------------- endpoints
def st_(v):
    return "UNASSESSED" if v is None else None


def fmt(v):
    return "UNASSESSED (data missing)" if v is None else f"Δ={v[0]:+.5f}, 95% CI [{v[1]:+.5f}, {v[2]:+.5f}] ({v[3]} cells, {v[4]} cell-seeds)"


E = ["# Expanded benchmark — endpoint evaluation (rules: expanded_benchmark_preregistration.md)", "",
     "HG001 and HG007 are **previously evaluated** samples; this is NOT an independent-genome validation. Δ = absolute F1 difference. CIs = 2-Mb block bootstrap within 2 contigs/sample; pooled = 4 cells of 2 samples (not a population interval).", ""]
E += ["## P1 — gain over R1.1 at 30× (PASS needs pooled Δ ≥ +0.003 and CI lower > 0)", ""]
for lvl in ("HG001", "HG007", "POOLED(4 cells)"):
    v = B[(30, "r2", "r11", lvl)]
    state = "UNASSESSED" if v is None else ("PASS" if v[0] >= 0.003 and v[1] > 0 else "FAIL")
    note = "" if v is None else (" — any improvement (P1b): " + ("PASS" if v[1] > 0 else "FAIL"))
    E.append(f"* {lvl}: {fmt(v)} → **P1 {state}**{note}")
E += ["", "## P2 — non-inferiority to bcftools (PASS if CI lower ≥ −0.0005); superiority if lower > 0", ""]
for d in DEPTHS:
    E.append(f"### {d}×")
    for lvl in ("HG001", "HG007", "POOLED(4 cells)"):
        v = B[(d, "r2", "bcftools", lvl)]
        state = "UNASSESSED" if v is None else ("PASS" if v[1] >= MARGIN else "FAIL")
        sup = "" if v is None else (", superior (lower>0)" if v[1] > 0 else "")
        E.append(f"* {lvl}: {fmt(v)} → **P2 {state}**{sup}")
E += ["", "## P4 — no low-depth regression vs R1.1 (point Δ ≥ −0.002)", ""]
for d in (10, 15):
    for lvl in ("HG001", "HG007", "POOLED(4 cells)"):
        v = B[(d, "r2", "r11", lvl)]
        E.append(f"* {d}× {lvl}: {fmt(v)} → **P4 {'UNASSESSED' if v is None else ('PASS' if v[0] >= -0.002 else 'FAIL')}**")
E += ["", "## P3 — absolute level (mean over seeds per cell; R2 ≥ 0.990; bcftools ceiling check)", "", "| depth | sample | contig | R2 F1 | R1.1 F1 | bcftools F1 | R2 ≥ 0.990 | bcftools ≥ 0.990 | claim |", "|---|---|---|---|---|---|---|---|---|"]
for d in (30, 50):
    for cell in CELLS:
        sc = seeds_complete(cell, d, ["r2", "r11", "bcftools"])
        if not sc:
            E.append(f"| {d}× | {cell[0]} | {cell[1]} | UNASSESSED | | | | | UNASSESSED |"); continue
        m = {c: st.mean(x[c]["SNP_F1"] for x in sc.values()) for c in ("r2", "r11", "bcftools")}
        claim = "PASS (ceiling check ok)" if m["r2"] >= .99 and m["bcftools"] >= .99 else ("PASS (reached by R2 only; data-limited claim not allowed)" if m["r2"] >= .99 else "FAIL")
        E.append(f"| {d}× | {cell[0]} | {cell[1]} | {m['r2']:.5f} | {m['r11']:.5f} | {m['bcftools']:.5f} | {m['r2'] >= .99} | {m['bcftools'] >= .99} | {claim} |")
E += ["", "## Mean SNP F1 per cell/depth (seeds averaged; n seeds)", "", "| sample | contig | depth | n seeds | R2 | R1.1 | bcftools | provenance |", "|---|---|---|---|---|---|---|---|"]
for cell in CELLS:
    for d in DEPTHS:
        sc = seeds_complete(cell, d, ["r2", "r11", "bcftools"])
        if sc:
            m = {c: st.mean(x[c]["SNP_F1"] for x in sc.values()) for c in ("r2", "r11", "bcftools")}
            pv = sorted({x["r2"]["prov"] for x in sc.values()})
            E.append(f"| {cell[0]} | {cell[1]} | {d}× | {len(sc)} | {m['r2']:.5f} | {m['r11']:.5f} | {m['bcftools']:.5f} | {'/'.join(pv)} |")
        else:
            E.append(f"| {cell[0]} | {cell[1]} | {d}× | 0 | UNASSESSED | | | |")
E += ["", "## Resource guardrails (R2 vs R1.1; median over seeds). Verdicts use NEW rows only (sequential, load recorded); EXISTING rows shown for context", "",
      "| sample | contig | depth | prov | wall Δ% (Δs) | CPU Δ% | RSS Δ% | verdict |", "|---|---|---|---|---|---|---|---|"]
for cell in CELLS:
    for d in DEPTHS:
        sc = seeds_complete(cell, d, ["r2", "r11"])
        if not sc:
            continue
        med = lambda c, k: st.median(x[c][k] for x in sc.values())
        wr = med("r2", "wall_s") / med("r11", "wall_s") - 1; wd = med("r2", "wall_s") - med("r11", "wall_s")
        cr = med("r2", "cpu_s") / med("r11", "cpu_s") - 1; rr = med("r2", "peak_rss_kb") / med("r11", "peak_rss_kb") - 1
        pv = sorted({x["r2"]["prov"] for x in sc.values()})[0] if len({x["r2"]["prov"] for x in sc.values()}) == 1 else "MIXED"
        ok = (wr <= .10 or wd <= 2) and cr <= .10 and rr <= .15
        verdict = ("PASS" if ok else "FAIL") if pv == "NEW" else ("within" if ok else "outside") + " (context only)"
        E.append(f"| {cell[0]} | {cell[1]} | {d}× | {pv} | {wr:+.1%} ({wd:+.1f}s) | {cr:+.1%} | {rr:+.1%} | {verdict} |")
E += ["", f"Metric consistency check (F1 from TP/FP/FN vs hap.py F1, tol 2e-4): {len(bad)} inconsistencies"] + [f"- {b}" for b in bad]
(OUT / "ENDPOINTS_EXPANDED.md").write_text("\n".join(E) + "\n")

# ---------------- DeepVariant separate table (30x seed 101; same hap.py/truth/BED/ref)
D = ["# DeepVariant comparison (separate table; descriptive only; same hap.py 0.3.15 / truth / BED / reference / BAM / seed 101 for all callers)", "",
     "No non-inferiority or superiority claim vs DeepVariant is preregistered. CIs: 2-Mb block bootstrap, single seed. HG001/HG007 are previously evaluated samples.", ""]
dv = [(c, d) for c in CELLS for d in DEPTHS if 101 in G.get((c[0], c[1], d), {}) and "deepvariant" in G[(c[0], c[1], d)][101]]
if not dv:
    D.append("**UNASSESSED** — no successful DeepVariant run (see deepvariant_feasibility.md / FAILED_RUNS_EXPANDED.md).")
else:
    D += ["| sample | contig | depth | caller | TP | FP | FN | precision | recall | F1 | F1 95% CI | wall s | CPU s | peak RSS MB |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    rng2 = np.random.default_rng(SEED + 1)
    for cell, d in dv:
        recs = G[(cell[0], cell[1], d)][101]
        keys = sorted({b for r in recs.values() for b, *_ in r["snp_blocks_2Mb"]})
        idx = rng2.integers(0, len(keys), (RES, len(keys)))
        for cal in ("r2", "r11", "bcftools", "deepvariant"):
            if cal not in recs:
                continue
            r = recs[cal]; lo, hi = np.percentile(f1_draws(arr(r, keys), idx), [2.5, 97.5])
            D.append(f"| {cell[0]} | {cell[1]} | {d}× | {NAME[cal]} | {r['SNP_TP']} | {r['SNP_FP']} | {r['SNP_FN']} | {r['SNP_precision']:.5f} | {r['SNP_recall']:.5f} | {r['SNP_F1']:.5f} | [{lo:.5f}, {hi:.5f}] | {r['wall_s']:.0f} | {r['cpu_s']:.0f} | {r['peak_rss_kb']/1024:.0f} |")
(OUT / "DEEPVARIANT_TABLE.md").write_text("\n".join(D) + "\n")
fails = [r for r in rows if r.get("status") != "OK"]
print("\n".join(E[:60])); print("rows", len(rows), "ok", len(OK), "failed", len(fails), "inconsistencies", len(bad))
