"""Render the bench_v19 (Validation Round 2) results JSON as markdown tables.

Pure formatting, computes nothing new: reuses summarize_v14's table renderers,
re-pointed at the v19 artefacts, plus new tables for the two statistics that
are genuinely new this round (block bootstrap, McNemar) and the cutoff
sensitivity sweep.
"""

from __future__ import annotations

import json
from pathlib import Path

import summarize_v14 as v14

RESULTS = Path("results/bench_v19/benchmark_results.json")
DISAGREEMENTS = Path("results/bench_v19/disagreements.json")

v14.ORDER = ["chr2_ordinary", "chr3_difficult", "chr5_stress"]
v14.TAGS = ["full", "30x", "15x"]
v14.DISAGREEMENTS = DISAGREEMENTS


def block_bootstrap_table(report: dict) -> str:
    lines = ["| Cell | n_blocks | block bp | ΔF1 (block) | 95% CI (block) | half-width | ΔF1 (locus) | 95% CI (locus) |",
             "|---|---:|---:|---:|---|---:|---:|---|"]
    for tag in v14.TAGS:
        for base in v14.ORDER:
            key = f"{base}_{tag}"
            cell = report["cells"].get(key)
            if cell is None:
                continue
            bb = cell["block_bootstrap"]
            locus = cell["vs_pb_paired_bootstrap"]
            lines.append(
                f"| {key} | {bb['n_blocks']} | {bb['block_bp']:,} | {bb['delta_f1']:+.6f} | "
                f"[{bb['ci95'][0]:+.6f}, {bb['ci95'][1]:+.6f}] | {bb['ci_half_width']:.6f} | "
                f"{locus['delta_f1']:+.6f} | [{locus['delta_f1_ci95'][0]:+.6f}, {locus['delta_f1_ci95'][1]:+.6f}] |")
    for key, entry in report.get("pooled", {}).items():
        bb = entry.get("block_bootstrap")
        locus = entry.get("vs_pb_paired_bootstrap")
        if not bb or not locus:
            continue
        lines.append(
            f"| **pooled/{key}** | {bb['n_blocks']} | {bb['block_bp']:,} | {bb['delta_f1']:+.6f} | "
            f"[{bb['ci95'][0]:+.6f}, {bb['ci95'][1]:+.6f}] | {bb['ci_half_width']:.6f} | "
            f"{locus['delta_f1']:+.6f} | [{locus['delta_f1_ci95'][0]:+.6f}, {locus['delta_f1_ci95'][1]:+.6f}] |")
    return "\n".join(lines)


def mcnemar_table(report: dict) -> str:
    lines = ["| Cell | b (cascade wrong, PB right) | c (cascade right, PB wrong) | n discordant | p-value |",
             "|---|---:|---:|---:|---:|"]
    for tag in v14.TAGS:
        for base in v14.ORDER:
            key = f"{base}_{tag}"
            cell = report["cells"].get(key)
            if cell is None:
                continue
            m = cell["mcnemar"]
            lines.append(f"| {key} | {m['b_cascade_wrong_pb_right']} | {m['c_cascade_right_pb_wrong']} | "
                        f"{m['n_discordant']} | {m['p_value']:.4g} |")
    for key, entry in report.get("pooled", {}).items():
        m = entry.get("mcnemar")
        if not m:
            continue
        lines.append(f"| **pooled/{key}** | {m['b_cascade_wrong_pb_right']} | {m['c_cascade_right_pb_wrong']} | "
                    f"{m['n_discordant']} | {m['p_value']:.4g} |")
    return "\n".join(lines)


def sensitivity_table(report: dict) -> str:
    lines = ["| Cell | pct | cutoff | F1 | FP | FN | routed fraction |",
             "|---|---|---:|---:|---:|---:|---:|"]
    for tag in v14.TAGS:
        for base in v14.ORDER:
            key = f"{base}_{tag}"
            cell = report["cells"].get(key)
            if cell is None:
                continue
            for pct, row in cell["cutoff_sensitivity"].items():
                lines.append(f"| {key} | {pct} | {row['cutoff']:.4f} | {row['f1']:.6f} | "
                            f"{row['fp']} | {row['fn']} | {row['routed_fraction']*100:.4f}% |")
    return "\n".join(lines)


def main() -> None:
    report = json.loads(RESULTS.read_text())
    print("### Main result table\n")
    print(v14.main_table(report))
    print("\n### Pooled\n")
    print(v14.pooled_table(report))
    print("\n### All arms\n")
    print(v14.arms_table(report))
    print("\n### Compute\n")
    print(v14.compute_table(report))
    print("\n### Depth-stratified (fine bins)\n")
    print(v14.depth_table(report))
    print("\n### Structure-stratified\n")
    print(v14.structure_table(report))
    print("\n### Controls\n")
    print(v14.controls_table(report))
    print("\n### Failure boundary\n")
    print(v14.failure_table(report))
    print("\n### Block bootstrap vs locus bootstrap\n")
    print(block_bootstrap_table(report))
    print("\n### McNemar's test\n")
    print(mcnemar_table(report))
    print("\n### Router cutoff sensitivity (diagnostic only)\n")
    print(sensitivity_table(report))
    print("\n### Disagreement effects\n")
    print(v14.disagreement_effects())


if __name__ == "__main__":
    main()
