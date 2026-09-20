"""Regenerates Figure 5 from the corrected PERFORMANCE_RESULTS.csv."""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 10, "figure.dpi": 150})

rows = list(csv.DictReader(open("../PERFORMANCE_RESULTS.csv")))

bars = []  # (label, seconds, color, group)
color_map = {
    "chtslib_extraction_stage": "#55a868",
    "full_pipeline_extract_bench_v12": "#8172b2",
    "htslib_thread_scaling": "#c44e52",
    "multiprocess_chunking": "#4c72b0",
}
for r in rows:
    if r["measured_or_extrapolated"] != "MEASURED":
        continue
    if r["wall_seconds"] in ("", "PENDING"):
        continue
    label = f"{r['experiment_id']}\n{r['configuration']} ({r['region_mb']}Mb)"
    bars.append((label, float(r["wall_seconds"]), color_map.get(r["experiment_id"], "#777"), r["experiment_id"]))

fig, ax = plt.subplots(figsize=(9, 6))
labels = [b[0] for b in bars]
vals = [b[1] for b in bars]
colors = [b[2] for b in bars]
y = range(len(labels))
ax.barh(list(y), vals, color=colors)
ax.set_yticks(list(y))
ax.set_yticklabels(labels, fontsize=7)
ax.set_xscale("log")
ax.set_xlabel("Wall-clock seconds (log scale) -- MEASURED only")
ax.set_title(
    "Figure 5: Extraction runtime, three DISTINCT mechanisms\n"
    "GREEN = C/htslib extension (integrated, this session): real 21-65x extraction-stage speedup.\n"
    "PURPLE = full extract_bench_v12.py pipeline (load_counts+load_reads+PB): smaller real speedup,\n"
    "         because load_reads/PB are unaccelerated and dominate.\n"
    "RED = htslib-internal thread scaling: genuinely ~0x (a different, separate, unintegrated mechanism).\n"
    "BLUE = process-level chunking: real but modest (~3.1x@8 workers), unintegrated.",
    fontsize=8, ha="center",
)
fig.tight_layout()
fig.savefig("fig5_runtime_python_vs_threading.png")
print("Wrote fig5_runtime_python_vs_threading.png")
