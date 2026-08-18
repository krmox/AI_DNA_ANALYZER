"""Render the bench_v16 results JSON as the markdown tables devlog 16 needs.

Pure formatting: it computes nothing the benchmark did not already compute, so
the devlog and the artefact can never disagree. The table renderers are
``summarize_v14``'s, re-pointed at the v16 artefacts — only the generic
stratum table (VAF, base quality, mapping quality, Stage-1 hard strata) is new,
because devlog 16 §6 reports axes devlog 14 stored but never rendered.
"""

from __future__ import annotations

import json
from pathlib import Path

import summarize_v14 as v14

RESULTS = Path("results/bench_v16/benchmark_results.json")
DISAGREEMENTS = Path("results/bench_v16/disagreements.json")

v14.ORDER = ["chr13_representative"]
v14.TAGS = ["full", "30x", "15x"]
v14.DISAGREEMENTS = DISAGREEMENTS


def stratum_table(report: dict, axis: str, title: str) -> str:
    lines = [f"**{title}**", "",
             "| Cell | stratum | loci | SNP | F1 PB | F1 cascade | ΔF1 | routed→PB | disagree |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for tag in v14.TAGS:
        cell = report["cells"].get(f"chr13_representative_{tag}")
        if cell is None:
            continue
        for name, block in cell[axis].items():
            if not block.get("loci"):
                continue
            # ``quality_stratified`` stores the two F1 values but no delta; the
            # difference of two numbers already in the artefact is formatting,
            # not a new metric.
            delta = block.get("delta_f1_vs_pb",
                              block["cheap_router_pb"]["f1"] - block["pb_only"]["f1"])
            lines.append(
                f"| {tag} | {name} | {block['loci']:,} | {block['snp']} | "
                f"{block['pb_only']['f1']:.6f} | {block['cheap_router_pb']['f1']:.6f} | "
                f"{delta:+.6f} | "
                f"{block['fraction_routed_to_pb'] * 100:.4f}% | "
                f"{block.get('disagreement_vs_pb', '—')} |")
    return "\n".join(lines)


def disagreement_detail() -> str:
    if not DISAGREEMENTS.exists():
        return "(no disagreement dump)"
    data = json.loads(DISAGREEMENTS.read_text())
    lines = ["| cell | pos | truth SNP | routed | effect | depth | VAF | baseQ | MAPQ | binom LLR | PB LLR | annotations |",
             "|---|---:|---|---|---|---:|---:|---:|---:|---:|---:|---|"]
    for key, entry in data.items():
        for record in entry["records"]:
            lines.append(
                f"| {key.replace('chr13_representative_', '')} | {record['position']:,} | "
                f"{'yes' if record['truth_is_snp'] else 'no'} | "
                f"{'yes' if record['routed_to_pb'] else 'no'} | `{record['effect']}` | "
                f"{record['depth']:.0f} | {record['vaf']:.3f} | {record['mean_base_quality']:.1f} | "
                f"{record['mean_mapq']:.1f} | {record['binomial_llr']:+.2f} | "
                f"{record['pb_llr']:+.2f} | {', '.join(record['annotations']) or '—'} |")
    return "\n".join(lines) if len(lines) > 2 else "No disagreement between the frozen cascade and PB-only."


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
    # The corrected mapping-quality axis lives in its own artefact (devlog 16
    # §17); splice it over the report's stale by_mapq block before rendering.
    corrected = Path("results/bench_v16/mapq_stratified.json")
    if corrected.exists():
        for key, block in json.loads(corrected.read_text())["cells"].items():
            report["cells"][key]["by_mapq"] = block
    for axis, title in (("by_vaf", "VAF-stratified"),
                        ("by_quality", "Base-quality-stratified"),
                        ("by_mapq", "Mapping-quality-stratified (corrected bins)"),
                        ("strata", "Stage-1 hard strata")):
        print(f"\n### {title}\n")
        print(stratum_table(report, axis, title))
    print("\n### Controls\n")
    print(v14.controls_table(report))
    print("\n### Failure boundary\n")
    print(v14.failure_table(report))
    print("\n### Disagreement effects\n")
    print(v14.disagreement_effects())
    print("\n### Every disagreement\n")
    print(disagreement_detail())
    for tag in v14.TAGS:
        key = f"chr13_representative_{tag}"
        if key in report["cells"]:
            print(f"\n### Operating curve — {tag}\n")
            print(v14.curve_table(report, key))


if __name__ == "__main__":
    main()
