"""Figures 1, 3, 5 and 6 of PAPER_MANUSCRIPT_V2.md (Figures 2 and 4: make_v2_figures.py).

Numbers are copied from research/TABLES_FINAL.md, PAPER_MANUSCRIPT_V2.md tables, the hap.py /
HEAD_TO_HEAD result files and LOAD_READS_OPTIMIZATION_REPORT.md. Figure 4a is read directly from the frozen
results/bench_v14|16|17|18|19|20 disagreement JSON files (no recomputation). Nothing is simulated.
Run:  python3 make_v2_figures_more.py   (matplotlib)
"""
import json
import os
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
from matplotlib.lines import Line2D

import make_v2_figures as base
from make_v2_figures import BLUE, VERM, GREY, INK, MUTED, OUT

ROOT = os.path.abspath(os.path.join(OUT, "..", ".."))
GREEN, SKY, YEL, PINK = "#009E73", "#56B4E9", "#E69F00", "#CC79A7"


# ------------------------------------------------------------------ Figure 1
def box(ax, x, y, w, h, title, sub="", fc="#eff6ff", ec=BLUE, tc=INK, ls="-", fs=8.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06", fc=fc, ec=ec, lw=1.2, ls=ls, zorder=2))
    ax.text(x + w / 2, y + h * (0.66 if sub else 0.5), title, ha="center", va="center", fontsize=fs, fontweight="bold", color=tc, zorder=3)
    if sub:
        ax.text(x + w / 2, y + h * 0.28, sub, ha="center", va="center", fontsize=7.3, color=MUTED, zorder=3, linespacing=1.25)


def arrow(ax, p, q, color=INK, ls="-", rad=0.0, lw=1.2, text=None, tpos=None, tc=MUTED):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=11, color=color, lw=lw, ls=ls,
                                 connectionstyle=f"arc3,rad={rad}", zorder=1))
    if text:
        ax.text(*(tpos or ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2 + 0.1)), text, fontsize=7.2, color=tc, ha="center", va="bottom")


def fig1():
    fig, ax = plt.subplots(figsize=(13.4, 6.0))
    ax.set_xlim(0, 14); ax.set_ylim(0, 6.4); ax.axis("off")
    # inputs
    for k, (t, s) in enumerate([("Aligned reads", "BAM, MAPQ $\\geq$ 20"), ("Reference", "FASTA (GRCh38)"), ("Confident regions", "GIAB v4.2.1 BED")]):
        box(ax, 0.1, 4.75 - k * 1.55, 1.9, 1.05, t, s, fc="#f9fafb", ec="#9ca3af")
        arrow(ax, (2.0, 5.27 - k * 1.55), (2.75, 3.75 if k == 1 else (4.35 if k == 0 else 3.15)), color="#9ca3af")
    box(ax, 2.75, 2.55, 2.35, 2.2, "Window extraction",
        "64-bp windows; a window is dropped\nif any base is outside the BED\nnine count channels per locus\nbackend: native C/htslib\n(pysam = reference / fallback)", fc="#f3f4f6", ec="#6b7280")
    box(ax, 5.55, 3.15, 2.05, 1.6, "Binomial screen", "fixed error $\\varepsilon$ = 0.01\n$\\mathrm{LLR_{bin}}$; call if $\\geq$ 7.0\nfrozen threshold 7.0", fc="#e0f2fe")
    arrow(ax, (5.1, 3.95), (5.55, 3.95))
    box(ax, 8.05, 3.15, 2.2, 1.6, "Router", "$|\\mathrm{LLR_{bin}} - 7.0| \\leq$ 5.4119\nband [1.588, 12.412]\nmeasures uncertainty, not error", fc="#fff7ed", ec=YEL)
    arrow(ax, (7.6, 3.95), (8.05, 3.95))
    # not routed
    box(ax, 8.05, 1.05, 2.2, 1.15, "Not routed", "binomial call\n(over 99.8 % of loci)", fc="#f9fafb", ec="#9ca3af")
    arrow(ax, (9.15, 3.15), (9.15, 2.2), text="outside band", tpos=(9.75, 2.55))
    # routed
    box(ax, 10.7, 4.55, 2.1, 1.45, "Read tensor", "$\\leq$ 48 reads per locus\nnative C or pure Python\nBLAKE2b read ordering", fc="#f3f4f6", ec="#6b7280")
    arrow(ax, (10.25, 4.35), (10.7, 5.0), text="routed:\n0.04-0.16 %", tpos=(10.05, 5.05))
    box(ax, 10.7, 2.85, 2.1, 1.35, "Poisson-binomial", "per-read error $q_i$\n$\\mathrm{LLR_{PB}}$; call if $\\geq$ 10.5\nfrozen threshold 10.5", fc="#dbeafe")
    arrow(ax, (11.75, 4.55), (11.75, 4.2))
    # merge -> calls -> genotype
    box(ax, 10.7, 0.5, 2.1, 1.3, "Allele calls", "cascade output\nSNP-only", fc="#ecfdf5", ec=GREEN)
    arrow(ax, (11.75, 2.85), (11.75, 1.8))
    arrow(ax, (10.25, 1.62), (10.7, 1.15))
    box(ax, 13.05, 0.5, 0.9, 1.3, "GT / GQ", "Method C\n(allele set\nunchanged)", fc="#ecfdf5", ec=GREEN, fs=7.8)
    ax.texts[-1].set_fontsize(6.2)
    arrow(ax, (12.8, 1.15), (13.05, 1.15))
    box(ax, 12.85, 3.25, 1.1, 2.2, "Scoring", "internal locus\nevaluator\nand\nhap.py 0.3.15", fc="#faf5ff", ec=PINK, fs=8.2)
    ax.texts[-1].set_fontsize(6.6)
    arrow(ax, (13.5, 1.8), (13.5, 3.25), color=PINK)
    # PB-only control
    box(ax, 5.55, 0.35, 2.05, 1.3, "PB-only control", "PB at every locus,\nthreshold 10.5\n(benchmark extractor)", fc="white", ec=VERM, ls="--")
    arrow(ax, (4.0, 2.55), (5.55, 1.2), color=VERM, ls="--", rad=0.15)
    arrow(ax, (7.6, 0.6), (10.7, 0.6), color=VERM, ls="--", text="paired comparison: $\\Delta F_1$ = cascade $-$ PB-only", tpos=(9.15, 0.12))
    ax.text(0.1, 6.15, "Frozen: 7.0 / 10.5 / 5.4119. Not routed loci never reach PB in the cascade; the benchmark computes PB everywhere only to build the control.",
            fontsize=8, color=MUTED, va="center")
    fig.suptitle("Figure 1. Cascade architecture and data flow", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(OUT, "fig1_v2_architecture.png"), dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------ Figure 4
def load_records():
    recs = []
    for e in ("v14", "v16", "v17", "v18", "v19", "v20"):
        f = os.path.join(ROOT, "results", f"bench_{e}", "disagreements.json" if e != "v20" else "chr20_disagreements.json")
        d = json.load(open(f))
        cells = d.items() if e != "v20" else [("chr20", d)]
        for _, v in cells:
            assert not v.get("records_truncated")
            for r in v["records"]:
                recs.append((e, r["effect"], r["binomial_llr"], r["vaf"]))
    return recs


def fig4():
    recs = load_records()
    cnt = Counter((e, k) for e, k, *_ in recs)
    # sanity gates against the manuscript's Table 5 (raw records only)
    assert cnt[("v14", "router_fn")] == 5 and cnt[("v17", "router_fn")] == 2 and cnt[("v19", "router_fn")] == 1 and cnt[("v20", "router_fn")] == 1
    assert sum(v for (e, k), v in cnt.items() if k == "router_avoids_fp") == 154
    assert sum(v for (e, k), v in cnt.items() if k == "router_fp" and e != "v20") == 13 and cnt[("v20", "router_fp")] == 8

    fig, (a, b) = plt.subplots(1, 2, figsize=(13.2, 5.6), gridspec_kw={"width_ratios": [1.6, 1]})
    a.axvspan(1.588, 12.412, color="#e5e7eb", alpha=.8, lw=0, zorder=0)
    a.axvline(7.0, color="#9ca3af", lw=.9, ls="--", zorder=1)
    style = {"router_avoids_fp": (SKY, "o", 22, "PB false positive avoided (n = 154)"),
             "router_fp": (BLUE, "^", 46, "cascade-introduced false positive (n = 21 with records)"),
             "router_fn": (VERM, "D", 62, "lost true SNP (n = 9 with covariates)")}
    for k in ("router_avoids_fp", "router_fp", "router_fn"):
        col, mk, sz, lab = style[k]
        pts = [(x, y, e) for e, kk, x, y in recs if kk == k]
        a.scatter([p[0] for p in pts], [p[1] for p in pts], s=sz, marker=mk, c=col, alpha=.75 if k == "router_avoids_fp" else .95,
                  edgecolors=["black" if p[2] == "v20" else "white" for p in pts], linewidths=[1.3 if p[2] == "v20" else .5 for p in pts], zorder=3, label=lab)
    a.set_xlabel("binomial LLR at the locus"); a.set_ylabel("variant allele fraction")
    a.set_title("a  Where each class of departure sits relative to the router band", loc="left", fontsize=9.5, fontweight="bold")
    a.text(7.0, 0.655, "router band [1.588, 12.412]", ha="center", fontsize=7.6, color=MUTED)
    a.set_ylim(0, 0.68); a.set_xlim(-23, 29)
    a.grid(color="#eef0f3", lw=.5, zorder=0)
    for s in ("top", "right"): a.spines[s].set_visible(False)
    a.legend(loc="upper left", fontsize=7.4, frameon=False, bbox_to_anchor=(0.0, 0.93))
    a.text(-22.5, 0.02, "black ring = v20 (chr20); the v15 lost SNP has no stored covariates; HG005 events are record-only and not plotted",
           fontsize=6.9, color=MUTED)

    # (b) composition of classified departures: pre-chr20 vs v20 (Table 5)
    groups = [("v14-v19\n(pre-chr20)", 15, 153, 9), ("v20\nchr20 whole", 8, 1, 1)]
    for i, (lab, fp, av, lo) in enumerate(groups):
        tot = fp + av + lo
        left = 0
        for v, col, nm in ((av, SKY, "PB FP avoided"), (fp, BLUE, "cascade-introduced FP"), (lo, VERM, "lost true SNP")):
            b.barh(i, v / tot * 100, left=left, color=col, edgecolor="white", height=.5, lw=1.2)
            if v / tot > .06:
                b.text(left + v / tot * 50, i, f"{v}", ha="center", va="center", fontsize=8, color="white", fontweight="bold")
            left += v / tot * 100
        b.text(101, i, f"n = {tot}", va="center", fontsize=8, color=INK)
    b.set_yticks([0, 1]); b.set_yticklabels([g[0] for g in groups], fontsize=8.4); b.invert_yaxis()
    b.set_xlim(0, 112); b.set_xlabel("share of classified departures (%)")
    b.set_title("b  Composition reverses on chr20", loc="left", fontsize=9.5, fontweight="bold")
    for s in ("top", "right"): b.spines[s].set_visible(False)
    b.legend(handles=[Patch(color=SKY, label="PB FP avoided (cascade gains)"), Patch(color=BLUE, label="cascade-introduced FP"),
                      Patch(color=VERM, label="lost true SNP")], loc="lower center", bbox_to_anchor=(0.45, -0.36), fontsize=7.6, frameon=False, ncol=1)
    fig.text(0.60, 0.012, "Pre-chr20 counts include 2 record-only HG005 false positives and the v15 lost SNP.\n"
             "v14's 2 router_recovers_fn events and the 80 v15 departures are not classified here.", fontsize=6.8, color=MUTED, va="bottom")
    fig.suptitle("Figure 3. Cascade-specific failure analysis", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(os.path.join(OUT, "fig3_v2_failure_analysis.png"), dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------ Figure 5
CALL = [  # name, precision, recall, F1, TP, FP, FN, colour, marker  (Table 9; hap.py 0.3.15)
    ("V1 PB-only + Method C", 0.9785, 0.9319, 0.9547, 15748, 346, 1150, BLUE, "o"),
    ("V1 cascade + Method C", 0.9783, 0.9319, 0.9546, 15748, 349, 1150, SKY, "o"),
    ("DeepVariant 1.6.1", 0.9340, 0.9654, 0.9494, 16313, 1153, 585, GREEN, "s"),
    ("Clair3 1.0.10", 0.8965, 0.9638, 0.9290, 16287, 1880, 611, YEL, "s"),
    ("GATK HC 4.5.0.0 (no BQSR)", 0.9877, 0.9493, 0.9682, 16042, 199, 856, PINK, "s"),
]


def fig5():
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.0, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
    r = np.linspace(0.90, 1.0, 400)
    for f1 in (0.93, 0.95, 0.97):
        p = f1 * r / (2 * r - f1)
        a.plot(r, p, color="#d1d5db", lw=.9, ls="--", zorder=1)
        idx = np.argmin(abs(r - 0.975)); a.text(0.9775, p[idx], f"F1 {f1}", fontsize=7, color=MUTED, va="bottom")
    for n, p, rc, f1, tp, fp, fn, c, m in CALL:
        a.scatter(rc, p, s=90, c=c, marker=m, edgecolors="white", linewidths=1, zorder=3)
    lab = {"V1 PB-only + Method C": (-0.0015, 0.0055, "right"), "V1 cascade + Method C": (0.0, -0.0075, "center"),
           "DeepVariant 1.6.1": (0.002, -0.004, "left"), "Clair3 1.0.10": (0.002, -0.004, "left"), "GATK HC 4.5.0.0 (no BQSR)": (0.002, 0.003, "left")}
    for n, p, rc, f1, *_ in CALL:
        if n.startswith("V1"):
            continue
        dx, dy, ha = lab[n]
        a.text(rc + dx, p + dy, f"{n}\nF1 {f1:.4f}", fontsize=7.6, ha=ha, va="center")
    a.text(0.9345, 0.9745, "V1 + Method C\nPB-only F1 0.9547\ncascade F1 0.9546\n(the two markers overlap)", fontsize=7.6, ha="left", va="top")
    a.axvline(0.982, color=BLUE, lw=.9, ls=":", zorder=1)
    a.text(0.9825, 0.8925, "V1 recall\nceiling 0.982\n(window filter)", fontsize=7, color=BLUE, va="bottom")
    a.set_xlim(0.925, 0.995); a.set_ylim(0.885, 1.0)
    a.set_xlabel("recall"); a.set_ylabel("precision")
    a.set_title("a  hap.py, HG002 chr21:32-44 Mb, 15x, SNP, PASS", loc="left", fontsize=9.5, fontweight="bold")
    a.grid(color="#eef0f3", lw=.5, zorder=0)
    for s in ("top", "right"): a.spines[s].set_visible(False)

    names = [c[0] for c in CALL][::-1]
    y = np.arange(len(CALL))
    fps = [c[5] for c in CALL][::-1]; fns = [c[6] for c in CALL][::-1]
    b.barh(y + .2, fps, height=.36, color=VERM, label="false positives")
    b.barh(y - .2, fns, height=.36, color=BLUE, label="false negatives")
    for yi, v in zip(y + .2, fps): b.text(v + 25, yi, f"{v:,}", va="center", fontsize=7.6)
    for yi, v in zip(y - .2, fns): b.text(v + 25, yi, f"{v:,}", va="center", fontsize=7.6)
    b.set_yticks(y); b.set_yticklabels(names, fontsize=8)
    b.set_xlim(0, 2300); b.set_xlabel("records (TRUTH.TOTAL = 16,898 for every caller)")
    b.set_title("b  Shape of the errors", loc="left", fontsize=9.5, fontweight="bold")
    b.legend(frameon=False, fontsize=8, loc="lower right")
    for s in ("top", "right"): b.spines[s].set_visible(False)
    fig.text(0.01, 0.015, "One region, one sample, one depth. The region contains the chr21:32-40 Mb span used to select the V1 constants (not held out). "
             "GATK ran without BQSR. Single runs. No ranking is implied and no runtime is compared.", fontsize=7.4, color=MUTED)
    fig.suptitle("Figure 7. External caller comparison", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(os.path.join(OUT, "fig6_v2_external_callers.png"), dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------ Figure 6
SPEED = [  # label, value, level, record-only
    ("load_counts, HG005 3 Mb", 21.03, "function", False), ("load_counts, HG002 12 Mb (earlier run)", 65.1, "function", False),
    ("load_counts, HG002 0.2 Mb", 29.9, "function", False), ("load_reads, HG005 30x", 26.9, "function", False),
    ("load_reads, HG002 15x", 35.1, "function", False), ("load_reads, Python threads x4", 0.73, "function", False),
    ("pipeline, native serial (800 kb)", 2.72, "pipeline", False), ("pipeline, native + 4 procs vs pure-Python serial", 7.11, "pipeline", False),
    ("pipeline, native 1 -> 4 workers (3 Mb)", 309 / 146, "pipeline", False),
    ("pipeline, counts-only native (3 Mb)", 1.22, "pipeline", True),
    ("BAM->calls, cascade vs PB-only (0.5 Mb)", 1.65, "measured", False),
]
LCOL = {"function": SKY, "pipeline": BLUE, "measured": GREEN}


def fig6():
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.4, 6.0), gridspec_kw={"width_ratios": [1.25, 1]})
    n = len(SPEED)
    for i, (lab, v, lv, ro) in enumerate(SPEED):
        yy = n - 1 - i
        a.barh(yy, v, height=.6, color="white" if ro else LCOL[lv], edgecolor=LCOL[lv], hatch="///" if ro else None, lw=1, zorder=2)
        a.text(v * 1.08, yy, f"{v:.3g}x", va="center", fontsize=8)
    a.axvline(1, color=INK, lw=.9, zorder=3)
    a.set_yticks(range(n)); a.set_yticklabels([s[0] for s in SPEED][::-1], fontsize=8)
    a.set_xscale("log"); a.set_xlim(0.5, 140)
    a.set_xticks([0.5, 1, 2, 5, 10, 20, 50, 100]); a.set_xticklabels(["0.5", "1", "2", "5", "10", "20", "50", "100"])
    a.set_xlabel("speedup versus the pure-Python / pysam reference (log scale)")
    a.set_title("a  Speedups by level: function-level figures are not pipeline speedups", loc="left", fontsize=9.5, fontweight="bold")
    a.grid(axis="x", color="#eef0f3", lw=.5, zorder=0)
    for s in ("top", "right"): a.spines[s].set_visible(False)
    a.legend(handles=[Patch(color=SKY, label="function level"), Patch(color=BLUE, label="pipeline level"), Patch(color=GREEN, label="measured BAM->calls"),
                      Patch(fc="white", ec=BLUE, hatch="///", label="historical, record-only")], loc="lower right", fontsize=7.6, frameon=False)
    fig.text(0.01, 0.012, "Projected caller-stage speedups (171-396x) and whole-genome estimates are projections, not results, and are not plotted. "
             "\nHatched = historical, record-only. Stage shares use seconds summed over stages and workers; chr20 PB time is inflated by 6 workers on 4 physical cores; the binomial stage was not timed separately.",
             fontsize=6.9, color=MUTED, va="bottom")

    # (b) stage seconds -> shares
    cfg = [("pure Python\n3 Mb (record-only)", (344, 696, 351), True), ("counts-only native\n3 Mb (record-only)", (15, 775, 349), True),
           ("native, 1 worker\n3 Mb", (5, 15, 287), False), ("native, 6 workers\nchr20 (56 M loci)", (194, 444, 15036), False)]
    cols = (SKY, YEL, VERM)
    for i, (lab, st, ro) in enumerate(cfg):
        tot = sum(st); left = 0
        for v, c in zip(st, cols):
            share = v / tot * 100
            b.barh(i, share, left=left, color=c, edgecolor="white", height=.55, lw=1.2, hatch="///" if ro else None, alpha=.55 if ro else 1)
            if share > 7:
                b.text(left + share / 2, i, f"{share:.0f} %", ha="center", va="center", fontsize=7.6, color="white" if c == VERM else INK, fontweight="bold")
            left += share
        b.text(101, i, f"{tot:,} s", va="center", fontsize=7.6)
    b.set_yticks(range(len(cfg))); b.set_yticklabels([c[0] for c in cfg], fontsize=8); b.invert_yaxis()
    b.set_xlim(0, 118); b.set_xlabel("share of summed stage time (%)")
    b.set_title("b  The bottleneck moved to the Poisson-binomial stage", loc="left", fontsize=9.5, fontweight="bold")
    for s in ("top", "right"): b.spines[s].set_visible(False)
    b.legend(handles=[Patch(color=SKY, label="load_counts"), Patch(color=YEL, label="load_reads"), Patch(color=VERM, label="Poisson-binomial LLR")],
             loc="lower center", bbox_to_anchor=(0.45, -0.20), ncol=3, fontsize=7.6, frameon=False)
    fig.suptitle("Figure 5. Where the time goes", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    fig.savefig(os.path.join(OUT, "fig5_v2_performance.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    fig1(); fig4(); fig5(); fig6()
    print("wrote fig1_v2_architecture, fig3_v2_failure_analysis, fig6_v2_external_callers, fig5_v2_performance")
