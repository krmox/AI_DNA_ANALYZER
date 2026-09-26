"""Figures 6 and 8 of PAPER_MANUSCRIPT_V2.md (V2.x evidence tiers; chr20 300x external callers).

Every plotted number is read from an existing artefact:
  research/V2X_ACCURACY_PERFORMANCE_BENCHMARK.csv   (V2 control and V2.x A+C, HG002 chr20 dev / holdout / full)
  research/final_independent_eval/results_chr8_analysis.json  (HG003 chr8 hap.py summary + block-bootstrap intervals)
  research/benchmark_chr20_300x/happy/*/*.summary.csv          (V1/V2 control, DeepVariant, GATK; SNP PASS rows)
Nothing is simulated, fitted or edited by hand.  Run:  python3 make_v2x_figures.py   (needs matplotlib)
"""
import csv, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.abspath(os.path.join(HERE, ".."))
BLUE, VERM, GREY, INK, MUTED, GREEN = "#0072B2", "#D55E00", "#6b7280", "#1f2937", "#6b7280", "#009E73"
plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": "#9ca3af", "axes.linewidth": 0.6, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": INK})

def rows(path):
    return list(csv.DictReader(open(path)))

acc = rows(f"{R}/V2X_ACCURACY_PERFORMANCE_BENCHMARK.csv")
def pick(model, tag):
    for r in acc:
        if r["model"] == model and r["development_or_holdout"].startswith(tag):
            return {k: float(r[k]) for k in ("precision", "recall", "F1", "FP", "FN", "TP")}
    raise KeyError((model, tag))
ctrl = {t: pick("V2_CONTROL", t) for t in ("development", "holdout")}
axc = {t: pick("V2.x_A+C", t) for t in ("development", "holdout")}
c8 = json.load(open(f"{R}/final_independent_eval/results_chr8_analysis.json"))
s8 = c8["happy_summary_SNP"]["PASS"]
ci = c8["block_bootstrap_95ci"]
chr8 = dict(precision=float(s8["METRIC.Precision"]), recall=float(s8["METRIC.Recall"]), F1=float(s8["METRIC.F1_Score"]))

# ---------------------------------------------------------------- Figure 6
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.9), gridspec_kw={"width_ratios": [1.15, 1]})
x = [0, 1, 2]
labels = ["HG002 chr20\ndevelopment blocks\n(in-sample)", "HG002 chr20\nholdout blocks", "HG003 chr8\nfinal evaluation\n(whole chromosome)"]
a = ax[0]
a.plot(x[:2], [ctrl["development"]["F1"], ctrl["holdout"]["F1"]], "o", mfc="white", mec=GREY, ms=8, label="V2 control (same blocks)")
a.plot(x, [axc["development"]["F1"], axc["holdout"]["F1"], chr8["F1"]], "o", color=BLUE, ms=8, label="V2.x A+C")
a.errorbar([2], [chr8["F1"]], yerr=[[chr8["F1"] - ci["f1"][0]], [ci["f1"][1] - chr8["F1"]]], color=BLUE, capsize=4, lw=1.4)
for xi, v in zip(x, [axc["development"]["F1"], axc["holdout"]["F1"], chr8["F1"]]):
    a.annotate(f"{v:.5f}", (xi, v), textcoords="offset points", xytext=(10, 4), fontsize=7.5, color=BLUE)
for xi, v in zip(x[:2], [ctrl["development"]["F1"], ctrl["holdout"]["F1"]]):
    a.annotate(f"{v:.5f}", (xi, v), textcoords="offset points", xytext=(10, -10), fontsize=7.5, color=GREY)
a.annotate("no V2 control run", (2, 0.9729), ha="center", fontsize=7.5, color=MUTED)
a.set_xticks(x); a.set_xticklabels(labels, fontsize=7.6); a.set_xlim(-0.5, 2.6); a.set_ylim(0.965, 0.995)
a.set_ylabel("hap.py SNP F1"); a.set_title("(a) F1 by evidence tier", fontsize=9, loc="left")
a.legend(loc="center left", fontsize=7.5, frameon=False, bbox_to_anchor=(0.0, 0.62))
a.text(2.55, 0.9655, "bar: 95 % block-bootstrap CI (chr8 only)", ha="right", fontsize=7, color=MUTED)
b = ax[1]
for j, (nm, key, off) in enumerate((("Precision", "precision", -0.12), ("Recall", "recall", 0.12))):
    vals = [axc["development"][key], axc["holdout"][key], chr8[key]]
    b.plot([xi + off for xi in x], vals, "o" if key == "precision" else "s", color=BLUE if key == "precision" else GREEN, ms=7, label=f"V2.x A+C {nm.lower()}")
    b.errorbar([2 + off], [chr8[key]], yerr=[[chr8[key] - ci[key][0]], [ci[key][1] - chr8[key]]], color=BLUE if key == "precision" else GREEN, capsize=3, lw=1.2)
    ctl = [ctrl["development"][key], ctrl["holdout"][key]]
    b.plot([0 + off, 1 + off], ctl, "o" if key == "precision" else "s", mfc="white", mec=GREY, ms=6, label=f"V2 control {nm.lower()}" if True else None)
b.set_xticks(x); b.set_xticklabels(labels, fontsize=7.6); b.set_xlim(-0.5, 2.5); b.set_ylim(0.955, 1.0)
b.set_title("(b) Precision and recall", fontsize=9, loc="left"); b.set_ylabel("value")
b.legend(fontsize=7, frameon=False, loc="lower left", ncol=1)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig6_v2x_evidence_tiers.png"), dpi=200); plt.close(fig)

# ---------------------------------------------------------------- Figure 8
def happy_snp(path):
    for r in csv.DictReader(open(path)):
        if r["Type"] == "SNP" and r["Filter"] == "PASS":
            return float(r["METRIC.Precision"]), float(r["METRIC.Recall"]), float(r["METRIC.F1_Score"])
B = f"{R}/benchmark_chr20_300x/happy"
ext = [("V1 / V2 control", happy_snp(f"{B}/ai_dna_analyzer/ai_cascade.summary.csv"), GREY, "s"),
       ("DeepVariant 1.10.0", happy_snp(f"{B}/deepvariant/deepvariant.summary.csv"), VERM, "D"),
       ("GATK 4.6.2.0", happy_snp(f"{B}/gatk/gatk.summary.csv"), GREEN, "^")]

fig, a = plt.subplots(figsize=(5.6, 4.0))
import numpy as np
P, Rr = np.meshgrid(np.linspace(0.95, 1.0, 200), np.linspace(0.95, 1.0, 200)); F = 2 * P * Rr / (P + Rr)
cs = a.contour(Rr, P, F, levels=[0.96, 0.97, 0.98, 0.99, 0.995], colors="#d1d5db", linewidths=0.7)
a.clabel(cs, fmt="F1=%.3f", fontsize=6.5)
for nm, (p, r, f), col, mk in ext:
    a.plot(r, p, mk, color=col, ms=8, label=f"{nm} (F1 {f:.4f})")
p2 = {r["model"]: r for r in acc}
full_row = [r for r in acc if r["model"] == "V2.x_A+C" and r["development_or_holdout"].startswith("full")][0]
a.plot(float(full_row["recall"]), float(full_row["precision"]), "o", color=BLUE, ms=8, label=f"V2.x A+C, whole chr20 (F1 {float(full_row['F1']):.4f})")
a.plot(axc["holdout"]["recall"], axc["holdout"]["precision"], "o", mfc="white", mec=BLUE, ms=8, label=f"V2.x A+C, holdout blocks (F1 {axc['holdout']['F1']:.4f})")
a.set_xlim(0.975, 1.0); a.set_ylim(0.955, 1.0005); a.set_xlabel("Recall (hap.py SNP, PASS)"); a.set_ylabel("Precision")
a.legend(fontsize=7, frameon=True, framealpha=0.95, edgecolor="#e5e7eb", loc="lower right")
a.set_title("HG002 chr20, 300x: descriptive, not a ranking", fontsize=9, loc="left")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig8_v2x_chr20_external.png"), dpi=200); plt.close(fig)
print("ok")
