"""Regenerates the manuscript figures (research/PAPER_MANUSCRIPT.md) from recorded results.

Run from research/FIGURES/:  python3 make_manuscript_figures.py

Every number below is copied from a named source file (see FIGURE_PLAN.md); nothing is
simulated. Fig 7 reads the HG005 post-fix cache directly (path below) to split the
'PB-rescuable' loci into true-SNP rescues vs false-positive avoidance.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.rcParams.update({"font.size": 9, "figure.dpi": 150, "axes.spines.top": False,
                     "axes.spines.right": False})
BLUE, ORANGE, GREEN, RED, PURPLE, GREY = "#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b2", "#8c8c8c"
HERE = os.path.dirname(os.path.abspath(__file__))
HG005_CACHE = os.path.join(HERE, "..", "..", "cache", "hg005_stress_postfix", "chr1_1000000_4000000.npz")


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, name))
    plt.close(fig)


# ------------------------------------------------------------------ Fig 1
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.axis("off")
ax.set_xlim(0, 12)
ax.set_ylim(-1.0, 7)


def box(x, y, w, h, text, fc="#eef3f8", ec="#333", fs=8):
    ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", fc=fc, ec=ec, lw=1.1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs)


def arrow(x0, y0, x1, y1, txt=None):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", color="#333"))
    if txt:
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.15, txt, fontsize=7, ha="center", color="#333")


box(0.2, 5.3, 2.0, 1.2, "BAM (Illumina,\nGRCh38-aligned)")
box(3.0, 5.9, 2.9, 0.9, "Native C/htslib pileup\n(bam_mplp_*; default)", fc="#e3f1e6")
box(3.0, 4.6, 2.9, 0.9, "Python/pysam pileup\n(reference / fallback)", fc="#f3f3f3")
box(6.8, 5.3, 2.4, 1.2, "Per-locus counts\n+ 64-bp window tiling\n(BED / N filter)")
box(9.9, 5.3, 2.0, 1.2, "Binomial LLR\n(thr 7.0, eps=0.01)")
arrow(2.2, 5.9, 3.0, 6.3)
arrow(2.2, 5.9, 3.0, 5.05)
arrow(5.9, 6.3, 6.8, 6.0)
arrow(5.9, 5.05, 6.8, 5.7)
arrow(9.2, 5.9, 9.9, 5.9)
box(9.6, 3.3, 2.6, 1.2, "Router\n|LLR - 7.0| <= 5.4119", fc="#fbeede")
arrow(10.9, 5.3, 10.9, 4.5)
box(6.6, 1.4, 3.0, 1.3, "Poisson-binomial LLR\n(thr 10.5; <=48 reads/locus;\nload_reads extraction)", fc="#f5e6e6")
arrow(10.2, 3.3, 8.6, 2.7, "routed")
box(9.9, 1.4, 2.0, 1.3, "Binomial call\n(LLR >= 7.0)")
arrow(11.3, 3.3, 11.1, 2.7, "not routed")
box(6.6, -0.3, 5.3, 1.0, "Cascade call set (allele level)", fc="#eef3f8")
arrow(8.1, 1.4, 8.1, 0.7)
arrow(10.9, 1.4, 10.9, 0.7)
box(0.2, -0.3, 5.5, 1.0, "Method C genotype posterior\n(0/0, 0/1, 1/1; eps=0.01) -> VCF", fc="#e9e4f4")
arrow(6.6, 0.2, 5.7, 0.2)
ax.text(0.2, 3.3, "Frozen constants: 7.0 / 10.5 / 5.4119\nFit once (HG002 chr21 validation\nblocks); never refit.", fontsize=8, va="top")
ax.set_title("Figure 1. Cascade architecture and extraction backends", fontsize=10)
save(fig, "fig1_architecture_v2.png")

# ------------------------------------------------------------------ Fig 2
# (label, sample, evaluator, loci_scored, truth_snp)  -- research/BENCHMARK_MASTER.csv
exp = [
    ("v13 chr21 pool", "HG002", "internal", 24583412, 33893),
    ("v14 cross-chrom (4x1Mb)", "HG002", "internal", 10368078, 12261),
    ("v15 held-out (chr16/15/7/14)", "HG002", "internal", 9498705, 10863),
    ("v16 chr13:70-73Mb", "HG002", "internal", 8726751, 14274),
    ("v17 chr17:18-21Mb", "HG002", "internal", 7400175, 8484),
    ("v18 chr13:70-73Mb", "HG003", "internal", 8699208, 13266),
    ("v19 chr2/3/5", "HG004", "internal", 25336365, 38340),
    ("chr21:32-44Mb (hap.py)", "HG002", "hap.py", 10975654, 16898),
    ("chr1:1.0-4.0Mb (hap.py)", "HG005", "hap.py", 2402190, 4090),
]
colors = {"HG002": BLUE, "HG003": GREEN, "HG004": ORANGE, "HG005": RED}
fig, ax = plt.subplots(figsize=(7.4, 4.5))
ypos = np.arange(len(exp))[::-1]
for y, (lab, s, ev, loci, snp) in zip(ypos, exp):
    ax.barh(y, loci / 1e6, color=colors[s], hatch="//" if ev == "hap.py" else None,
            edgecolor="white" if ev != "hap.py" else "black", lw=0.6)
    ax.text(loci / 1e6 + 0.3, y, f"{snp:,} truth SNPs", va="center", fontsize=7)
ax.set_yticks(ypos)
ax.set_yticklabels([e[0] for e in exp], fontsize=8)
ax.set_xlabel("Scored loci (millions; depth cells pooled where applicable)")
ax.set_xlim(0, 33)
handles = [mpatches.Patch(color=c, label=s) for s, c in colors.items()]
handles += [mpatches.Patch(fc="white", ec="black", hatch="//", label="hap.py (genotype-aware)"),
            mpatches.Patch(fc="white", ec="#999", label="internal row-index evaluator")]
ax.legend(handles=handles, fontsize=7, loc="upper center", bbox_to_anchor=(0.45, -0.18), ncol=3, frameon=False)
ax.set_title("Figure 2. Validation experiments (HG002-4 form a related trio; HG005 is unrelated)", fontsize=9)
save(fig, "fig2_validation_design.png")

# ------------------------------------------------------------------ Fig 3
# (label, dF1, lo, hi, kind)  -- CIs: paired locus bootstrap, 10,000 resamples
# sources: research/BENCHMARK_MASTER.csv, results/bench_v13..v19, VALIDATION_ROUND_2.md, head_to_head report
fx = [
    ("v13 chr21 pool (3 of 13 regions inside router-fit span)", 0.000103, 0.0, 0.00022),
    ("v14 cross-chromosome (pooled; pre-reg. verdict NEGATIVE)", 0.000581, 0.000075, 0.001089),
    ("v15 held-out test (chr16/15/7/14)", 0.002686, 0.001968, 0.003414),
    ("v16 HG002 chr13", -0.0000349, -0.000171, 0.0000703),
    ("v17 HG002 chr17 segdup", 0.001892, 0.001250, 0.002615),
    ("v18 HG003 chr13", 0.000262, 0.000076, 0.000483),
    ("v19 HG004 chr2/3/5 (locus CI)", 0.001014, 0.000787, 0.001254),
    ("v19 HG004 chr2/3/5 (50-kb block CI)", 0.001014, 0.000520, 0.001601),
    ("HG002 chr21:32-44Mb, 15x (internal eval.)", -0.0000282, -0.000122, 0.0000647),
]
pts_only = [  # hap.py Delta F1, no CI computed
    ("HG002 chr21:32-44Mb hap.py + Method C (no CI)", 0.954569 - 0.954656),
    ("HG005 chr1:1-4Mb hap.py + Method C (no CI)", 0.951862 - 0.952098),
]
fig, ax = plt.subplots(figsize=(7.6, 4.4))
n = len(fx) + len(pts_only)
y = np.arange(n)[::-1]
for yi, (lab, d, lo, hi) in zip(y[:len(fx)], fx):
    ax.errorbar(d * 1e3, yi, xerr=[[ (d - lo) * 1e3], [(hi - d) * 1e3]], fmt="o", color=BLUE, capsize=2)
for yi, (lab, d) in zip(y[len(fx):], pts_only):
    ax.plot(d * 1e3, yi, "D", color=RED)
ax.axvline(0, color="black", lw=0.8)
ax.axvspan(-1, 1, color="#dddddd", alpha=0.5, lw=0)
ax.set_yticks(y)
ax.set_yticklabels([f[0] for f in fx] + [p[0] for p in pts_only], fontsize=7)
ax.set_xlabel("Delta F1 (cascade - PB-only) x 10^-3\n(grey band: pre-registered |dF1| < 0.001 margin)")
ax.set_title("Figure 3. Cascade vs PB-only across validation experiments", fontsize=9)
save(fig, "fig3_cascade_vs_pb_delta_f1.png")

# ------------------------------------------------------------------ Fig 4
fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 3.9), gridspec_kw={"width_ratios": [1.15, 1]})
# lost true-SNP events: (VAF, binomial LLR, label) -- results/bench_v14/disagreements.json, bench_v17, bench_v19
ev = [(0.125, 0.887, "v14"), (0.125, 0.665, "v14"), (0.1099, -5.297, "v14"), (0.120, -0.018, "v14"),
      (0.1224, 0.648, "v14"), (0.129, 1.57, "v17"), (0.115, -0.70, "v17"), (0.130, 1.349, "v19")]
cmap = {"v14": BLUE, "v17": ORANGE, "v19": GREEN}
for v, l, k in ev:
    a.scatter(v, l, color=cmap[k], s=28, zorder=3)
a.axhspan(7.0 - 5.4119, 7.0 + 5.4119, color="#fbeede", alpha=0.8, label="router band (PB invoked)")
a.axhline(7.0, color="k", ls="--", lw=0.8, label="binomial threshold 7.0")
a.set_xlim(0.0, 0.5)
a.set_ylim(-8, 15)
a.set_xlabel("Alt allele fraction (VAF)")
a.set_ylabel("Binomial LLR")
a.legend(handles=[mpatches.Patch(color=cmap[k], label=k) for k in cmap]
         + [mpatches.Patch(color="#fbeede", label="router band"),
            plt.Line2D([], [], color="k", ls="--", label="threshold 7.0")], fontsize=7, loc="upper right")
a.set_title("(a) 8 cascade-specific true-SNP losses (unrouted;\nPB LLR 11.8-24.1 would have called each)\n"
            "v14 = 5 events at 4 distinct loci", fontsize=8)
# HG005 FN taxonomy (n=244), research/ERROR_ANALYSIS.csv (flags overlap)
labs = ["window dropped\n(no candidate)", "alldifficult", "lowmap_segdup", "low VAF", "low depth", "low MAPQ", "low BQ", "tandem rep.\n(n=6)"]
vals = [142, 170, 123, 66, 46, 30, 19, 6]
cols = [GREY, RED, RED, RED, RED, RED, RED, "#e6b0b0"]
b.barh(np.arange(len(labs))[::-1], vals, color=cols)
b.set_yticks(np.arange(len(labs))[::-1])
b.set_yticklabels(labs, fontsize=7)
b.set_xlabel("HG005 false negatives (of 244; flags overlap)")
b.set_title("(b) HG005 chr1:1-4Mb false-negative flags\n(grey = candidate-generation artifact, M-1)", fontsize=8)
save(fig, "fig4_failure_taxonomy.png")

# ------------------------------------------------------------------ Fig 5
fig, ax = plt.subplots(figsize=(6.6, 3.9))
cases = [("HG002 chr21\n0.2 Mb", 16.14, 0.54), ("HG005 chr1\n3 Mb", 373.11, 17.74), ("HG002 chr21\n12 Mb", 1092.07, 16.76)]
x = np.arange(len(cases))
w = 0.36
ax.bar(x - w / 2, [c[1] for c in cases], w, color=GREY, label="Python/pysam (reference)")
ax.bar(x + w / 2, [c[2] for c in cases], w, color=GREEN, label="Native C/htslib")
for i, c in enumerate(cases):
    ax.text(i, c[1] * 1.25, ["29.9x", "21.0x", "65.1x"][i], ha="center", fontsize=9, fontweight="bold")
ax.set_yscale("log")
ax.set_ylim(0.3, 4000)
ax.set_xticks(x)
ax.set_xticklabels([c[0] for c in cases])
ax.set_ylabel("load_counts wall time (s, log scale)")
ax.legend(fontsize=8, loc="upper left")
ax.set_title("Figure 5. Extraction-stage (load_counts) speedup only.\n12-Mb pair: original measurement (source worktree), not re-run", fontsize=9)
save(fig, "fig5_native_extraction_speedup.png")

# ------------------------------------------------------------------ Fig 6
fig, ax = plt.subplots(figsize=(5.6, 4.2))
stages = [("load_counts (pileup)", 344, 15, GREEN), ("load_reads (PB read tensor)", 696, 775, BLUE), ("Poisson-binomial LLR", 351, 349, ORANGE)]
for j, (lab, old, new, c) in enumerate(stages):
    bo = sum(s[1] for s in stages[:j])
    bn = sum(s[2] for s in stages[:j])
    ax.bar(0, old, bottom=bo, color=c, label=lab)
    ax.bar(1, new, bottom=bn, color=c)
    ax.text(0, bo + old / 2, f"{old} s", ha="center", va="center", color="white", fontsize=8)
    ax.text(1, bn + new / 2 if new > 100 else bn + new + 40, f"{new} s", ha="center", va="center",
            color="white" if new > 100 else "black", fontsize=8)
ax.set_xticks([0, 1])
ax.set_xticklabels(["pysam extraction\n1392 s total", "native extraction\n1140 s total"])
ax.set_ylabel("Wall time, HG005 chr1:1.0-4.0Mb (s)")
ax.legend(fontsize=7, loc="upper right")
ax.set_title("Figure 6. Bottleneck migration: 1.22x end-to-end\n(single runs; native load_counts 22.9x on this stage)", fontsize=9)
save(fig, "fig6_bottleneck_migration.png")

# ------------------------------------------------------------------ Fig 7
if os.path.exists(HG005_CACHE):
    d = np.load(HG005_CACHE)
    lab, bl, pl = d["labels"], d["binomial_llr"], d["pb_llr"]
    fr = (lab == 0) | (lab == 1)
    t, bl, pl = (lab[fr] == 1), bl[fr], pl[fr]
    bcall, pcall = bl >= 7.0, pl >= 10.5
    routed = np.abs(bl - 7.0) <= 5.411872376933351
    resc = (pcall == t) & (bcall != t)
    tp_type, fp_type = resc & t, resc & ~t
    fig, (a, b) = plt.subplots(1, 2, figsize=(8.2, 3.6))
    cats = ["true-SNP rescue\n(PB calls, binomial misses)", "false-positive avoidance\n(binomial calls, PB rejects)"]
    routed_ct = [int((tp_type & routed).sum()), int((fp_type & routed).sum())]
    missed_ct = [int((tp_type & ~routed).sum()), int((fp_type & ~routed).sum())]
    a.bar(cats, routed_ct, color=GREEN, label="routed to PB")
    a.bar(cats, missed_ct, bottom=routed_ct, color=RED, label="not routed (missed)")
    for i in range(2):
        a.text(i, routed_ct[i] + missed_ct[i] + 2, f"{routed_ct[i] + missed_ct[i]}", ha="center", fontsize=8)
    a.set_ylabel("Loci")
    a.legend(fontsize=7)
    a.set_title(f"(a) {int(resc.sum())} PB-rescuable loci, HG005 chr1:1-4Mb", fontsize=8)
    n_r = int(routed.sum())
    needed_correct = int((routed & (bcall == t)).sum())
    b.bar(["routed loci"], [needed_correct], color=GREY, label="binomial already correct")
    b.bar(["routed loci"], [n_r - needed_correct], bottom=[needed_correct], color=BLUE, label="binomial wrong")
    b.text(0, n_r + 20, f"{n_r} routed ({100 * n_r / fr.sum():.3f}% of {int(fr.sum()):,} loci)", ha="center", fontsize=8)
    b.legend(fontsize=7, loc="center right")
    b.set_title("(b) what the router sends to PB", fontsize=8)
    fig.suptitle("Figure 7. Router rescue analysis (diagnostic; router not modified)", fontsize=9)
    save(fig, "fig7_router_rescue.png")
    print("fig7 counts", routed_ct, missed_ct, n_r, needed_correct)
print("done")
