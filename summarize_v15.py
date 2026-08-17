"""Render the bench_v15 results JSON as the markdown tables devlog 15 needs.

Pure formatting. It reads ``results/bench_v15/test_results.json`` and computes
nothing the benchmark did not already compute, so the devlog and the artefact
can never disagree.
"""

from __future__ import annotations

import json
from pathlib import Path

RESULTS = Path("results/bench_v15/test_results.json")
SUMMARY = Path("results/bench_v15/summary.md")

ORDER = ["chr16_segdup", "chr15_segdup", "chr7_segdup", "chr14_neutral"]
TAGS = ["full", "30x", "15x"]
POOLS = ["all", "segdup_regions", "chr16_segdup", "chr15_segdup", "chr7_segdup",
         "chr14_neutral", "depth_full", "depth_30x", "depth_15x"]


def _ci(block: dict) -> str:
    boot = block.get("vs_pb_paired_bootstrap") or {}
    ci = boot.get("delta_f1_ci95")
    return f"[{ci[0]:+.6f}, {ci[1]:+.6f}]" if ci else "—"


def main_table(report: dict) -> str:
    tag = "frozen_full"
    lines = ["### Per-cell: frozen cascade (C) vs safety layer (D) vs PB-only (B)", "",
             "| cell | depth | loci | SNP | F1 PB | F1 C | F1 D | ΔF1(D−PB) | ΔF1(D−C) | "
             "lost SNP C | lost SNP D | routed | escalated | PB budget ×C |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for base in ORDER:
        for depth in TAGS:
            cell = report["cells"].get(f"{base}_{depth}")
            if cell is None:
                continue
            arm = cell["arms"][tag]
            safety = cell["safety"][tag]
            lines.append(
                f"| `{base}` | {depth} | {cell['loci_scored']:,} | {cell['snp']:,} | "
                f"{cell['pb_only']['f1']:.6f} | {cell['cheap_router_pb']['f1']:.6f} | "
                f"{arm['f1']:.6f} | {arm['delta_f1_vs_pb']:+.6f} | "
                f"{arm['delta_f1_vs_frozen_cascade']:+.6f} | "
                f"{cell['cheap_router_pb']['lost_true_snp_vs_pb']} | "
                f"{arm['lost_true_snp_vs_pb']} | "
                f"{cell['fraction_routed_to_pb'] * 100:.4f}% | "
                f"{safety['escalated']} | {safety['pb_fraction_multiple_vs_frozen']:.2f} |")
    return "\n".join(lines)


def pool_table(report: dict, tag: str = "frozen_full") -> str:
    lines = [f"### Pooled results ({tag})", "",
             "| pool | loci | SNP | F1 B | F1 C | F1 D | ΔF1(D−PB) | 95% CI | "
             "lost C | lost D | PB frac C | PB frac D |",
             "|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|"]
    for name in POOLS:
        pool = report["pools"].get(name)
        if pool is None:
            continue
        arm = pool["arms"][tag]
        safety = pool["safety"][tag]
        lines.append(
            f"| `{name}` | {pool['loci_scored']:,} | {pool['snp']:,} | "
            f"{pool['pb_only']['f1']:.6f} | {pool['cheap_router_pb']['f1']:.6f} | "
            f"{arm['f1']:.6f} | {arm['delta_f1_vs_pb']:+.6f} | {_ci(arm)} | "
            f"{pool['cheap_router_pb']['lost_true_snp_vs_pb']} | "
            f"{arm['lost_true_snp_vs_pb']} | "
            f"{pool['fraction_routed_to_pb']:.3e} | {safety['total_pb_fraction']:.3e} |")
    return "\n".join(lines)


def all_arms_table(report: dict) -> str:
    pool = report["pools"]["all"]
    lines = ["### All arms, pooled test set", "",
             "| arm | F1 | precision | recall | TP | FP | FN | TN | accuracy | "
             "balanced acc. | ΔF1 vs PB | lost true SNP | PB fraction |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    rows = [("A. binomial-only", pool["binomial_only"], 0.0),
            ("B. PB-only", pool["pb_only"], 1.0),
            ("C. frozen router → PB", pool["cheap_router_pb"], pool["fraction_routed_to_pb"])]
    for tag, arm in pool["arms"].items():
        rows.append((f"D. safety ({tag})", arm, pool["safety"][tag]["total_pb_fraction"]))
    for name, block, fraction in rows:
        lines.append(
            f"| {name} | {block['f1']:.6f} | {block['precision']:.6f} | {block['recall']:.6f} | "
            f"{block['tp']} | {block['fp']} | {block['fn']} | {block['tn']:,} | "
            f"{block['accuracy']:.8f} | {block['balanced_accuracy']:.6f} | "
            f"{block['delta_f1_vs_pb']:+.6f} | {block['lost_true_snp_vs_pb']} | {fraction:.3e} |")
    return "\n".join(lines)


def safety_table(report: dict, tag: str = "frozen_full") -> str:
    pool = report["pools"]["all"]["safety"][tag]
    lines = [f"### Safety-layer behaviour, pooled ({tag})", "", "| quantity | value |", "|---|---:|"]
    for key, value in pool.items():
        lines.append(f"| `{key}` | {value:,.6g} |" if isinstance(value, float)
                     else f"| `{key}` | {value:,} |")
    return "\n".join(lines)


def controls_table(report: dict, tag: str = "frozen_full") -> str:
    controls = report["pools"]["all"]["controls"][tag]
    pool = report["pools"]["all"]
    lines = [f"### Controls at matched extra PB coverage ({controls['matched_extra_pb_coverage']:,} loci)",
             "", "| arm | F1 | ΔF1 vs PB | TP | FP | FN | lost true SNP |", "|---|---:|---:|---:|---:|---:|---:|"]
    rows = [("frozen cascade (no extra PB)", {
        "f1": pool["cheap_router_pb"]["f1"],
        "delta_f1_vs_pb": pool["cheap_router_pb"]["delta_f1_vs_pb"],
        "tp": pool["cheap_router_pb"]["tp"], "fp": pool["cheap_router_pb"]["fp"],
        "fn": pool["cheap_router_pb"]["fn"],
        "lost_true_snp_vs_pb": pool["cheap_router_pb"]["lost_true_snp_vs_pb"]}),
        ("safety layer", controls["safety_layer"]),
        ("reverse-confidence escalation", controls["reverse_confidence_escalation"])]
    for name, block in rows:
        lines.append(f"| {name} | {block['f1']:.6f} | {block['delta_f1_vs_pb']:+.6f} | "
                     f"{block['tp']} | {block['fp']} | {block['fn']} | "
                     f"{block['lost_true_snp_vs_pb']} |")
    random = controls["random_escalation"]
    for index, run in enumerate(random["runs"]):
        lines.append(f"| random escalation, seed {index} | {run['f1']:.6f} | "
                     f"{run['delta_f1_vs_pb']:+.6f} | {run['tp']} | {run['fp']} | {run['fn']} | "
                     f"{run['lost_true_snp_vs_pb']} |")
    lines += ["", f"random mean F1 {random['mean_f1']:.6f} ± {random['std_f1']:.6f}; "
                  f"mean lost true SNP {random['mean_lost_true_snp']:.1f} "
                  f"(best of 5 seeds: {random['min_lost_true_snp']})"]
    return "\n".join(lines)


def curve_table(report: dict) -> str:
    rows = report["pools"]["all"]["accuracy_compute_curve"]
    rows = sorted(rows, key=lambda r: r["escalated_fraction"])
    lines = ["### Accuracy-vs-compute curve, pooled test set (diagnostic only)", "",
             "| candidate | escalated | escalated fraction | PB budget × frozen | F1 | "
             "ΔF1 vs PB | lost true SNP |", "|---|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        mark = " **(FROZEN-FULL)**" if row.get("name") == \
            report["frozen_rules"]["frozen_full"]["name"] else ""
        mark = " **(FROZEN-LEAN)**" if row.get("name") == \
            report["frozen_rules"].get("frozen_lean", {}).get("name") else mark
        lines.append(
            f"| `{row.get('name', '?')}`{mark} | {row['escalated']:,} | "
            f"{row['escalated_fraction']:.3e} | {row['pb_fraction_multiple_vs_frozen']:.2f} | "
            f"{row['f1']:.6f} | {row['delta_f1_vs_pb']:+.6f} | {row['lost_true_snp_vs_pb']} |")
    return "\n".join(lines)


def structure_table(report: dict, tag: str = "frozen_full") -> str:
    lines = ["### Structure-stratified (external GIAB v3.1 annotation, reporting only)", "",
             "| stratum | loci | SNP | escalated | F1 PB | F1 C | F1 D | lost C | lost D |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    totals: dict[str, dict] = {}
    for cell in report["cells"].values():
        for label, arms in cell.get("safety_by_structure", {}).items():
            block = arms[tag]
            bucket = totals.setdefault(label, {key: 0 for key in
                                               ("loci", "snp", "escalated",
                                                "lost_true_snp_vs_pb",
                                                "frozen_cascade_lost_true_snp_vs_pb")})
            for key in bucket:
                bucket[key] += block[key]
    for label, bucket in sorted(totals.items()):
        lines.append(f"| `{label}` | {bucket['loci']:,} | {bucket['snp']:,} | "
                     f"{bucket['escalated']:,} | — | — | — | "
                     f"{bucket['frozen_cascade_lost_true_snp_vs_pb']} | "
                     f"{bucket['lost_true_snp_vs_pb']} |")
    return "\n".join(lines)


def compute_table(report: dict) -> str:
    total = report["pooled_wallclock"]
    lines = ["### Compute and projected wall clock (pooled test set)", ""]
    for unit, title in (("caller_only_seconds", "Caller-only (devlogs 11–14 convention)"),
                        ("end_to_end_seconds", "End-to-end (incl. pileup construction)")):
        if unit not in total["pb_only"]:
            continue
        label = unit.replace("_seconds", "")
        lines += [f"**{title}**", "",
                  "| arm | seconds | speedup vs PB-only | × frozen cascade |",
                  "|---|---:|---:|---:|",
                  f"| PB-only | {total['pb_only'][unit]:,.1f} | 1.00× | — |",
                  f"| frozen cascade | {total['frozen_cascade'][unit]:,.1f} | "
                  f"{total['frozen_cascade'][f'{label}_speedup_vs_pb_only']:.1f}× | 1.000× |"]
        for tag, arm in total["arms"].items():
            lines.append(f"| safety ({tag}) | {arm[unit]:,.1f} | "
                         f"{arm[f'{label}_speedup_vs_pb_only']:.1f}× | "
                         f"{arm[f'{label}_speedup_relative_to_frozen_cascade']:.3f}× |")
        lines.append("")
    probe = next(iter(report["throughput"].values()), None)
    if probe:
        lines += ["", "Measured per-locus rates (single quiet process, one probe per cell; "
                      "first probe shown):", "",
                  f"* binomial caller: {probe['binomial_loci_per_second']:,.0f} loci/s",
                  f"* safety layer (production, gated): {probe['safety_loci_per_second']:,.0f} "
                  f"loci/s ({probe['safety_cost_relative_to_binomial']:.2f}× the binomial caller)",
                  f"* safety layer (naive dense form, not used): "
                  f"{probe['sample_loci'] / probe['dense_safety_seconds']:,.0f} loci/s "
                  f"({probe['dense_safety_cost_relative_to_binomial']:.2f}× the binomial caller)",
                  f"* PB: {probe['pb_loci_per_second']:,.0f} loci/s "
                  f"({probe['pb_cost_relative_to_safety']:,.0f}× the safety layer per locus)"]
    return "\n".join(lines)


def acceptance_table(report: dict) -> str:
    lines = ["### Acceptance criteria (devlog 15 §7.5) and verdict (§7.6)", ""]
    for tag, criteria in report["acceptance"].items():
        lines += [f"**{tag}** — frozen cascade lost {criteria['frozen_cascade_lost_true_snp']} "
                  f"true SNP; this arm lost {criteria['safety_arm_lost_true_snp']}.", "",
                  "| criterion | result | evidence |", "|---|---|---|",
                  f"| A1 failure repaired | {'PASS' if criteria['A1_failure_repaired'] else 'FAIL'}"
                  f"{'' if criteria['A1_informative'] else ' (UNINFORMATIVE)'} | "
                  f"{criteria['frozen_cascade_lost_true_snp']} → "
                  f"{criteria['safety_arm_lost_true_snp']} lost true SNP |",
                  f"| A2 non-inferior to PB | {'PASS' if criteria['A2_non_inferior_to_pb'] else 'FAIL'} | "
                  f"ΔF1 95% CI [{criteria['A2_delta_f1_ci95'][0]:+.6f}, "
                  f"{criteria['A2_delta_f1_ci95'][1]:+.6f}], bound −0.001 |",
                  f"| A3 compute preserved | {'PASS' if criteria['A3_compute_preserved'] else 'FAIL'} | "
                  f"PB budget {criteria['A3_pb_fraction_multiple']:.2f}× (≤2); "
                  f"{criteria['A3_accounting']} speedup "
                  f"{criteria['A3_speedup_relative_to_frozen_cascade']:.3f}× of cascade (≥0.8); "
                  f"caller-only "
                  f"{criteria['A3_caller_only_speedup_relative_to_frozen_cascade']:.3f}× → "
                  f"{'would also pass' if criteria['A3_would_pass_on_caller_only_accounting'] else 'would FAIL'} |",
                  f"| A4 beats matched-budget controls | "
                  f"{'PASS' if criteria['A4_beats_matched_budget_controls'] else 'FAIL'} | "
                  f"random best {criteria['A4_random_min_lost']}, reverse "
                  f"{criteria['A4_reverse_confidence_lost']}, safety "
                  f"{criteria['safety_arm_lost_true_snp']} lost |",
                  "", f"**Verdict ({tag}): {report['verdict'][tag]}**", ""]
    return "\n".join(lines)


def main() -> None:
    report = json.loads(RESULTS.read_text())
    blocks = [main_table(report), pool_table(report), all_arms_table(report),
              safety_table(report), controls_table(report), curve_table(report),
              structure_table(report), compute_table(report), acceptance_table(report)]
    if "frozen_lean" in report["frozen_rules"]:
        blocks.insert(2, pool_table(report, "frozen_lean"))
        blocks.append(controls_table(report, "frozen_lean"))
    text = "\n\n".join(blocks) + "\n"
    SUMMARY.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
