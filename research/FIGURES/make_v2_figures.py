"""Figures 2 and 4 of PAPER_MANUSCRIPT_V2.md.

Every plotted number is copied from research/TABLES_FINAL.md (Tables 2, S1, S2, 3c),
research/CHR20_RESULTS.csv / results/bench_v20/chr20_results.json, or the hap.py
summaries named below. Nothing is recomputed, simulated or fitted.
Run from anywhere:  python3 make_v2_figures.py   (needs matplotlib)
Outputs: fig2_v2_accuracy_delta_f1.png, fig4_v2_router_behavior.png
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

OUT = os.path.dirname(os.path.abspath(__file__))
BLUE, VERM, GREY, INK, MUTED = "#0072B2", "#D55E00", "#6b7280", "#1f2937", "#6b7280"
COL = {"IMPROVED": BLUE, "PRESERVED": GREY, "DEGRADED": VERM}
MK = {"IMPROVED": "o", "PRESERVED": "s", "DEGRADED": "D"}
plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": "#9ca3af", "axes.linewidth": 0.6,
                     "text.color": INK, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK})

# ------------------------------------------------------------------ Figure 2
# Panel a: experiment level, internal evaluator, depths pooled (TABLES_FINAL.md Table 2).
# (label, dF1, lo, hi, block_lo, block_hi, pooled verdict, experiment-level note)
A = [
 ("v13  HG002 chr21, 13 regions",            +0.000103, 0.0,       +0.00022,  None,      None,      "PRESERVED", "PRESERVED (lower CI printed as 0.0000)"),
 ("v14  HG002, 4 x 1 Mb",                    +0.000581, +0.000075, +0.001089, None,      None,      "IMPROVED",  "pooled IMPROVED; experiment NEGATIVE"),
 ("v15  HG002, 4 held-out contigs",          +0.002686, +0.001968, +0.003414, None,      None,      "IMPROVED",  "IMPROVED"),
 ("v16  HG002 chr13:70-73 Mb",               -0.000035, -0.000171, +0.000070, None,      None,      "PRESERVED", "PRESERVED"),
 ("v17  HG002 chr17:18-21 Mb",               +0.001892, +0.001250, +0.002615, None,      None,      "IMPROVED",  "IMPROVED"),
 ("v18  HG003 chr13:70-73 Mb",               +0.000262, +0.000076, +0.000483, None,      None,      "IMPROVED",  "IMPROVED"),
 ("v19  HG004 chr2/3/5",                     +0.001014, +0.000787, +0.001254, +0.000520, +0.001601, "IMPROVED",  "IMPROVED"),
 ("v20  HG002 chr20 whole, 15x only",        -0.0000564, -0.000105, -0.0000141, -0.000108, -0.0000075, "DEGRADED", "DEGRADED; 10 discordant loci in 56.2 M"),
]
# Internal evaluator on the U-H2 region (TABLES_FINAL.md Table 2, row H2H).
H2H = ("U-H2  HG002 chr21:32-44 Mb, 15x (internal)", -0.0000282, -0.000122, +0.0000647)
# hap.py points, no CI exists. dF1 = PB-only F1 - cascade F1 from happy.summary.csv
# (chr21: 0.954656/0.954569; chr20: 0.962851/0.962796); HG005 from the project record (0.952098/0.951862).
HAPPY = [
 ("hap.py  HG002 chr21:32-44 Mb, 15x",  -0.000087, False),
 ("hap.py  HG002 chr20 whole, 15x",     -0.000055, False),
 ("hap.py  HG005 chr1:1-4 Mb  [record-only]", -0.000236, True),
]

# Panel b: per-cell, internal evaluator (TABLES_FINAL.md S1 = v14, S2 = v19).
# (cell, depth, dF1, lo, hi, verdict)
B14 = [
 ("chr20:33-34 Mb", "full", -0.001262, -0.003849, +0.001228, "DEGRADED"),
 ("chr20:33-34 Mb", "30x", 0, 0, 0, "PRESERVED"), ("chr20:33-34 Mb", "15x", 0, 0, 0, "PRESERVED"),
 ("chr19:1-2 Mb", "full", -0.000875, -0.002668, +0.000856, "PRESERVED"),
 ("chr19:1-2 Mb", "30x", 0, 0, 0, "PRESERVED"), ("chr19:1-2 Mb", "15x", -0.000454, -0.001415, 0, "PRESERVED"),
 ("chr4:101-102 Mb", "full", -0.000447, -0.001393, 0, "PRESERVED"),
 ("chr4:101-102 Mb", "30x", 0, 0, 0, "PRESERVED"), ("chr4:101-102 Mb", "15x", 0, 0, 0, "PRESERVED"),
 ("chr1:120-121 Mb (segdup)", "full", +0.005218, +0.001418, +0.008972, "IMPROVED"),
 ("chr1:120-121 Mb (segdup)", "30x", +0.003056, +0.000211, +0.006140, "IMPROVED"),
 ("chr1:120-121 Mb (segdup)", "15x", 0, 0, 0, "PRESERVED"),
]
B19 = [
 ("chr2:210-213 Mb", "full", +0.000493, +0.000121, +0.001005, "IMPROVED"),
 ("chr2:210-213 Mb", "30x", 0, 0, 0, "PRESERVED"), ("chr2:210-213 Mb", "15x", 0, 0, 0, "PRESERVED"),
 ("chr3:75-78 Mb", "full", +0.005582, +0.004019, +0.007273, "IMPROVED"),
 ("chr3:75-78 Mb", "30x", +0.002972, +0.001853, +0.004211, "IMPROVED"),
 ("chr3:75-78 Mb", "15x", +0.000124, 0, +0.000382, "PRESERVED"),
 ("chr5:27-30 Mb", "full", +0.000507, -0.000106, +0.001152, "PRESERVED"),
 ("chr5:27-30 Mb", "30x", +0.000102, 0, +0.000311, "PRESERVED"),
 ("chr5:27-30 Mb", "15x", -0.000102, -0.000313, 0, "PRESERVED"),
]


def fig2():
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.6, 7.6), gridspec_kw={"width_ratios": [1.0, 1.0]})
    k = 1e3
    rows = A + [(H2H[0], H2H[1], H2H[2], H2H[3], None, None, "PRESERVED", "PRESERVED")]
    n = len(rows) + len(HAPPY) + 1
    y = n
    ticks, labels = [], []
    for lab, d, lo, hi, blo, bhi, v, note in rows:
        y -= 1
        a.plot([lo * k, hi * k], [y, y], color=COL[v], lw=1.6, solid_capstyle="butt", zorder=3)
        if blo is not None:
            a.plot([blo * k, bhi * k], [y - 0.22, y - 0.22], color=COL[v], lw=1.0, alpha=.6, zorder=2)
        a.plot(d * k, y, MK[v], color=COL[v], ms=6, mec="white", mew=.8, zorder=4)
        ticks.append(y); labels.append(lab + "\n" + note)
    y -= 1.0
    for lab, d, ro in HAPPY:
        y -= 1
        a.plot(d * k, y, "^", ms=6.5, mfc="white" if ro else INK, mec=INK, mew=1.1, zorder=4)
        ticks.append(y); labels.append(lab + "\nno CI exists" + ("; raw output unavailable" if ro else ""))
    a.axvspan(-1, 1, color="#e5e7eb", alpha=.7, zorder=0, lw=0)
    a.axvline(0, color="#9ca3af", lw=.8, zorder=1)
    a.set_yticks(ticks); a.set_yticklabels(labels, fontsize=7.4, linespacing=1.15)
    a.set_xlim(-1.6, 3.6); a.set_ylim(y - .8, n - .2)
    a.set_xlabel(r"$\Delta F_1$ = cascade $-$ PB-only  ($\times10^{-3}$)")
    a.set_title("a  Experiment level (depths pooled)", loc="left", fontsize=9.5, fontweight="bold")
    a.text(0, n - .35, "pre-specified margin |dF1| < 0.001", ha="center", va="bottom", fontsize=7, color=MUTED)
    a.grid(axis="x", color="#e5e7eb", lw=.5, zorder=0)
    for s in ("top", "right"): a.spines[s].set_visible(False)
    a.set_clip_on(False)

    # panel b
    cells = [("v14 HG002", c) for c in B14] + [("v19 HG004", c) for c in B19]
    m = len(cells) + 1
    yy = m; tk, lb = [], []
    prev = None
    for tag, (cell, dep, d, lo, hi, v) in cells:
        if prev is not None and tag != prev: yy -= 0.8
        prev = tag; yy -= 1
        if lo != hi:
            b.plot([lo * k, hi * k], [yy, yy], color=COL[v], lw=1.4, solid_capstyle="butt", zorder=3)
        b.plot(d * k, yy, MK[v], color=COL[v], ms=5.2, mec="white", mew=.7, zorder=4)
        tk.append(yy); lb.append(f"{tag[:3]}  {cell}, {dep}")
    b.axvspan(-1, 1, color="#e5e7eb", alpha=.7, zorder=0, lw=0)
    b.axvline(0, color="#9ca3af", lw=.8, zorder=1)
    b.set_yticks(tk); b.set_yticklabels(lb, fontsize=7.4)
    b.set_ylim(yy - .8, m - .3)
    b.set_xlim(-4.6, 9.4)
    b.set_xlabel(r"$\Delta F_1$ per cell ($\times10^{-3}$); bars: 95 % locus bootstrap")
    b.set_title("b  v14 and v19 cells (region, depth)", loc="left", fontsize=9.5, fontweight="bold")
    b.grid(axis="x", color="#e5e7eb", lw=.5, zorder=0)
    for s in ("top", "right"): b.spines[s].set_visible(False)
    b.annotate("v14 chr20:33-34 Mb, full depth:\nDEGRADED by the |dF1| >= 0.001 clause\n(CI contains 0; 2 extra FPs)",
               xy=(-1.262, tk[0]), xytext=(2.4, tk[0] - 1.2), fontsize=7.4, color=VERM,
               arrowprops=dict(arrowstyle="-", color=VERM, lw=.7), va="center")

    h = [Line2D([], [], marker=MK[v], color=COL[v], ls="", ms=6, label=v) for v in ("IMPROVED", "PRESERVED", "DEGRADED")]
    h += [Line2D([], [], marker="^", color=INK, ls="", ms=6, label="hap.py point (no CI)"),
          Line2D([], [], marker="^", mfc="white", mec=INK, color=INK, ls="", ms=6, label="hap.py, record-only"),
          Line2D([], [], color=GREY, lw=1.0, alpha=.6, label="50-kb block bootstrap CI (v19, v20)"),
          Patch(color="#e5e7eb", label="|dF1| < 0.001")]
    fig.legend(handles=h, loc="lower center", ncol=7, frameon=False, fontsize=7.6, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Figure 2. Cascade versus PB-only: no global verdict; v14 NEGATIVE, v20 DEGRADED", x=0.01, ha="left",
                 fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(os.path.join(OUT, "fig2_v2_accuracy_delta_f1.png"), dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------ Figure 3
# Panel a: routed fraction, percent of scored loci (TABLES_FINAL.md Table 2; S1, S2 for cells).
ROUTED = [  # label, pooled %, hatched(record-only)
 ("v13 HG002 chr21 (13 reg.)", 0.1084, False), ("v14 HG002", 0.0704, False), ("v15 HG002", 0.0581, False),
 ("v16 HG002 chr13", 0.0392, False), ("v17 HG002 chr17", 0.0582, False), ("v18 HG003 chr13", 0.0409, False),
 ("v19 HG004 chr2/3/5", 0.0409, False), ("v20 HG002 chr20 whole (15x)", 0.13735, False),
 ("U-H2 HG002 chr21:32-44 Mb (15x)", 0.155, False), ("U-H5 HG005 chr1:1-4 Mb (~30x)  [record-only]", 0.0666, True),
]
# per-depth-cell routed % (S1 for v14, S2 for v19): full, 30x, 15x
CELL14 = {"full": [0.0053, 0.0255, 0.0030, 0.0507], "30x": [0.0325, 0.0599, 0.0133, 0.0597], "15x": [0.1743, 0.2184, 0.0936, 0.1233]}
CELL19 = {"full": [0.0030, 0.0042, 0.0036], "30x": [0.0143, 0.0146, 0.0142], "15x": [0.1070, 0.1032, 0.1037]}
DEPCOL = {"full": "#0072B2", "30x": "#E69F00", "15x": "#CC79A7"}

# Panel b: router capture of PB-rescuable loci (results/bench_v20/chr20_results.json rescue_composition;
# HG005 from the project record, raw unavailable). (label, routed, total, record_only)
CAP = [
 ("chr20  true-SNP rescuable\n(PB calls a real SNP the binomial misses)", 1327, 1328, False),
 ("chr20  FP-avoidance\n(binomial FP that PB rejects)",                    1490, 1498, False),
 ("HG005  true-SNP rescuable  [record-only]",                              13,   13,   True),
 ("HG005  FP-avoidance  [record-only]",                                    127,  129,  True),
]


def fig3():
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.0, 5.6), gridspec_kw={"width_ratios": [1.15, 1]})
    n = len(ROUTED)
    for i, (lab, v, ro) in enumerate(ROUTED):
        y = n - 1 - i
        a.barh(y, v, height=.55, color="white" if ro else "#93c5fd", edgecolor=BLUE, hatch="///" if ro else None, lw=.9, zorder=2)
        a.text(0.66, y, f"{v:.3g} %", va="center", ha="right", fontsize=7.6)
        if lab.startswith("v14"):
            for dep, vals in CELL14.items(): a.plot(vals, [y] * 4, "o", ms=3.6, color=DEPCOL[dep], alpha=.9, zorder=3)
        if lab.startswith("v19"):
            for dep, vals in CELL19.items(): a.plot(vals, [y] * 3, "o", ms=3.6, color=DEPCOL[dep], alpha=.9, zorder=3)
    a.set_yticks(range(n)); a.set_yticklabels([r[0] for r in ROUTED][::-1], fontsize=8)
    a.set_xscale("log"); a.set_xlim(0.002, 0.7)
    a.set_xlabel("loci routed to PB (% of scored loci; log scale)")
    a.set_title("a  Routing fraction", loc="left", fontsize=9.5, fontweight="bold")
    a.grid(axis="x", color="#e5e7eb", lw=.5, zorder=0)
    for s in ("top", "right"): a.spines[s].set_visible(False)
    ha = [Patch(fc="#93c5fd", ec=BLUE, label="pooled per experiment"), Patch(fc="white", ec=BLUE, hatch="///", label="record-only")]
    ha += [Line2D([], [], marker="o", ls="", color=c, ms=4.5, label=f"v14 / v19 cells, {d} depth") for d, c in DEPCOL.items()]
    fig.legend(handles=ha, loc="lower left", ncol=5, fontsize=7.2, frameon=False, bbox_to_anchor=(0.01, 0.0))

    for i, (lab, r, t, ro) in enumerate(CAP):
        y = len(CAP) - 1 - i
        f = r / t
        b.plot(f, y, "o", ms=7, mfc="white" if ro else BLUE, mec=BLUE, mew=1.3, zorder=3)
        b.text(f, y + 0.28, f"{r:,}/{t:,} routed  ({f:.4f})" + ("  n = 13" if lab.startswith("HG005  true") else ""),
               ha="center", fontsize=7.6, color=MUTED if ro else INK)
    b.set_yticks(range(len(CAP))); b.set_yticklabels([c[0] for c in CAP][::-1], fontsize=7.8)
    b.set_xlim(0.975, 1.0015); b.set_ylim(-0.6, len(CAP) - 0.3)
    b.axhline(1.5, color="#e5e7eb", lw=.8)
    b.set_xlabel("fraction of PB-rescuable loci routed to PB (axis from 0.975)")
    b.set_title("b  Router capture, by population", loc="left", fontsize=9.5, fontweight="bold")
    b.grid(axis="x", color="#e5e7eb", lw=.5, zorder=0)
    for s in ("top", "right"): b.spines[s].set_visible(False)
    fig.text(0.55, 0.045, "PB-rescuable = PB-only correct, binomial-only wrong. Capture is relative to PB, not truth recall.\n"
             "The mixed HG005 figure (140/142) is not shown as a true-SNP statistic.", fontsize=7, color=MUTED, va="bottom")
    fig.suptitle("Figure 4. Router behaviour", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0.09, 1, 0.95))
    fig.savefig(os.path.join(OUT, "fig4_v2_router_behavior.png"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    fig2(); fig3()
    print("wrote fig2_v2_accuracy_delta_f1.png, fig4_v2_router_behavior.png")
