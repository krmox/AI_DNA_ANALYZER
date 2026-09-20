"""Generates the 5 dossier figures from BENCHMARK_MASTER.csv / ERROR_ANALYSIS.csv
/ PERFORMANCE_RESULTS.csv. Run from research/FIGURES/: python3 make_figures.py
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams.update({"font.size": 10, "figure.dpi": 150})

# ---------------------------------------------------------------- Figure 1
fig, ax = plt.subplots(figsize=(5, 7))
ax.axis("off")
stages = ["BAM/CRAM", "Pileup extraction\n(Python/pysam)", "Binomial LLR\n(threshold=7.0)",
          "Router\n(|LLR-7.0|<=5.41)", "Poisson-Binomial LLR\n(threshold=10.5, routed loci only)",
          "Cascade decision", "Method C genotype\n(posterior over 0/0,0/1,1/1)", "VCF"]
y = list(range(len(stages), 0, -1))
for yi, s in zip(y, stages):
    box = mpatches.FancyBboxPatch((0.15, yi - 0.4), 0.7, 0.7, boxstyle="round,pad=0.02",
                                   linewidth=1.2, edgecolor="#333", facecolor="#eef3f8")
    ax.add_patch(box)
    ax.text(0.5, yi - 0.05, s, ha="center", va="center", fontsize=9)
for i in range(len(y) - 1):
    ax.annotate("", xy=(0.5, y[i + 1] + 0.3), xytext=(0.5, y[i] - 0.4),
                arrowprops=dict(arrowstyle="->", color="#333"))
ax.set_xlim(0, 1)
ax.set_ylim(0, len(stages) + 1)
ax.set_title("Figure 1: Cascade architecture (frozen thresholds)", fontsize=11)
fig.tight_layout()
fig.savefig("fig1_architecture.png")
plt.close(fig)

# ---------------------------------------------------------------- Figure 2
rows = list(csv.DictReader(open("../BENCHMARK_MASTER.csv")))
pairs = {}
for r in rows:
    if r["arm"] not in ("pb_only", "pb_only+MethodC", "cascade", "cheap_router_pb (cascade)", "cascade+MethodC"):
        continue
    if not r["f1"]:
        continue
    key = r["experiment_id"]
    pairs.setdefault(key, {})
    if "pb_only" in r["arm"]:
        pairs[key]["pb"] = (float(r["f1"]), r["evaluator"])
    else:
        pairs[key]["cascade"] = (float(r["f1"]), r["evaluator"])

labels, pb_vals, casc_vals, evaluators = [], [], [], []
for k, v in pairs.items():
    if "pb" in v and "cascade" in v:
        labels.append(k.replace("_", "\n"))
        pb_vals.append(v["pb"][0])
        casc_vals.append(v["cascade"][0])
        evaluators.append(v["pb"][1])

fig, ax = plt.subplots(figsize=(11, 5))
x = range(len(labels))
w = 0.35
b1 = ax.bar([i - w / 2 for i in x], pb_vals, width=w, label="PB-only", color="#4c72b0")
b2 = ax.bar([i + w / 2 for i in x], casc_vals, width=w, label="Cascade", color="#dd8452")
ax.set_xticks(list(x))
ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel("F1")
ax.set_ylim(0.9, 1.0)
ax.set_title("Figure 2: PB-only vs Cascade F1 across experiments\n"
             "(bar color = arm, NOT a ranking across evaluators -- see evaluator labels below)", fontsize=10)
for i, ev in enumerate(evaluators):
    ax.text(i, 0.905, ev.replace("_", "\n"), ha="center", fontsize=6, rotation=0, color="#555")
ax.legend()
fig.tight_layout()
fig.savefig("fig2_pb_vs_cascade_f1.png")
plt.close(fig)

# ---------------------------------------------------------------- Figure 3
rescue_rows = list(csv.DictReader(open("../RESCUE_ANALYSIS.csv")))
hg005 = {r["metric"]: r["value"] for r in rescue_rows if r["experiment_id"] == "HG005_postfix"}
labels = ["PB-rescuable\nloci", "Rescued\nby router", "Missed\nrescue", "Routed\ntotal", "False\nrescues"]
values = [float(hg005["pb_rescuable_loci"]), float(hg005["rescued_by_router"]),
          float(hg005["missed_rescue"]), float(hg005["routed_total"]), float(hg005["false_rescue_count"])]
fig, ax = plt.subplots(figsize=(7, 4.5))
bars = ax.bar(labels, values, color=["#4c72b0", "#55a868", "#c44e52", "#8172b2", "#c44e52"])
for b, v in zip(bars, values):
    ax.text(b.get_x() + b.get_width() / 2, v + max(values) * 0.01, f"{int(v)}", ha="center", fontsize=9)
ax.set_ylabel("Loci")
ax.set_title("Figure 3: Router rescue analysis, HG005 chr1:1,000,001-4,000,000\n"
             f"Rescue recall = {float(hg005['rescue_recall']):.3f}, "
             f"false rescue rate = {float(hg005['false_rescue_rate']):.4f} (diagnostic only, router not tuned)",
             fontsize=9)
fig.tight_layout()
fig.savefig("fig3_router_rescue.png")
plt.close(fig)

# ---------------------------------------------------------------- Figure 4
strata = ["lowmap_segdup", "alldifficult", "tandemrepeats"]
fn_enrich = [7.73, 1.66, 0.57]
fp_enrich = [6.35, 1.98, 7.46]
underpowered = [False, False, True]  # tandemrepeats FN n=6
fig, ax = plt.subplots(figsize=(7, 4.5))
x = range(len(strata))
w = 0.35
ax.bar([i - w / 2 for i in x], fn_enrich, width=w, label="FN enrichment", color="#c44e52")
ax.bar([i + w / 2 for i in x], fp_enrich, width=w, label="FP enrichment", color="#dd8452")
ax.axhline(1.0, color="#333", linewidth=0.8, linestyle="--", label="background (1.0x)")
ax.set_xticks(list(x))
ax.set_xticklabels(strata)
ax.set_ylabel("Enrichment over background rate")
ax.set_title("Figure 4: GIAB stratification error enrichment, HG005 chr1:1,000,001-4,000,000\n"
              "tandemrepeats FN (n=6) is UNDERPOWERED -- point estimate shown, do not over-interpret", fontsize=9)
for i, u in enumerate(underpowered):
    if u:
        ax.text(i - w / 2, fn_enrich[i] + 0.15, "n=6\n(low power)", ha="center", fontsize=7, color="#c44e52")
ax.legend()
fig.tight_layout()
fig.savefig("fig4_stratification_enrichment.png")
plt.close(fig)

# ---------------------------------------------------------------- Figure 5
perf_rows = list(csv.DictReader(open("../PERFORMANCE_RESULTS.csv")))
fig, ax = plt.subplots(figsize=(8, 5))
labels, vals, colors, tags = [], [], [], []
for r in perf_rows:
    if r["measured_or_extrapolated"] == "MEASURED" and r["wall_seconds"]:
        labels.append(f"{r['experiment_id']}\n{r['configuration']}"[:40])
        vals.append(float(r["wall_seconds"]))
        colors.append("#4c72b0" if "htslib_thread" in r["experiment_id"] else
                       ("#55a868" if "multiprocess" in r["experiment_id"] else "#8172b2"))
        tags.append(r["region"])
y = range(len(labels))
ax.barh(list(y), vals, color=colors)
ax.set_yticks(list(y))
ax.set_yticklabels(labels, fontsize=7)
ax.set_xlabel("Wall-clock seconds (MEASURED only -- see PERFORMANCE_RESULTS.csv for EXTRAPOLATION rows, excluded here)")
ax.set_title("Figure 5: Measured runtime -- htslib thread scaling (blue) shows NO speedup;\n"
              "process-level chunking (green) shows real, sub-linear speedup. No C/htslib\n"
              "extension exists in this codebase (see PIPELINE_INTEGRITY_AUDIT.md H-1).", fontsize=9)
fig.tight_layout()
fig.savefig("fig5_runtime_python_vs_threading.png")
plt.close(fig)

print("Wrote 5 figures.")
