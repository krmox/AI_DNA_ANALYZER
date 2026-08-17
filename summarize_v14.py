"""Render the bench_v14 results JSON as the markdown tables devlog 14 needs.

Pure formatting. It reads ``results/bench_v14/robustness_results.json`` and
computes nothing that the benchmark did not already compute, so the devlog and
the artefact can never disagree.
"""

from __future__ import annotations

import json
from pathlib import Path

RESULTS = Path("results/bench_v14/robustness_results.json")
DISAGREEMENTS = Path("results/bench_v14/disagreements.json")

ORDER = ["chr20_neutral", "chr19_gcrich", "chr4_atrich", "chr1_segdup"]
TAGS = ["full", "30x", "15x"]


def _f(value, digits=6):
    return "—" if value is None else f"{value:.{digits}f}"


def main_table(report: dict) -> str:
    lines = ["| Cell | depth | loci | SNP | F1 PB-only | F1 cascade | ΔF1 | ΔF1 95% CI | routed→PB | disagree | verdict |",
             "|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---|"]
    for base in ORDER:
        for tag in TAGS:
            key = f"{base}_{tag}"
            cell = report["cells"].get(key)
            if cell is None:
                continue
            boot = cell.get("vs_pb_paired_bootstrap", {})
            ci = boot.get("delta_f1_ci95")
            ci_text = f"[{ci[0]:+.6f}, {ci[1]:+.6f}]" if ci else "—"
            lines.append(
                f"| `{base}` | {tag} | {cell['loci_scored']:,} | {cell['snp']:,} | "
                f"{cell['pb_only']['f1']:.6f} | {cell['cheap_router_pb']['f1']:.6f} | "
                f"{cell['delta_f1_vs_pb']:+.6f} | {ci_text} | "
                f"{cell['fraction_routed_to_pb'] * 100:.4f}% | {cell['disagreement_vs_pb']} | "
                f"**{cell['verdict']}** |")
    return "\n".join(lines)


def arms_table(report: dict) -> str:
    lines = ["| Cell | arm | F1 | precision | recall | TP | FP | FN | accuracy | balanced acc |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for base in ORDER:
        for tag in TAGS:
            cell = report["cells"].get(f"{base}_{tag}")
            if cell is None:
                continue
            for arm, label in (("binomial_only", "binomial"), ("pb_only", "PB"),
                               ("cheap_router_pb", "cascade")):
                block = cell[arm]
                lines.append(
                    f"| `{base}` {tag} | {label} | {block['f1']:.6f} | {block['precision']:.6f} | "
                    f"{block['recall']:.6f} | {block['tp']} | {block['fp']} | {block['fn']} | "
                    f"{block['accuracy']:.8f} | {block.get('balanced_accuracy', float('nan')):.6f} |")
    return "\n".join(lines)


def pooled_table(report: dict) -> str:
    lines = ["| Pool | loci | SNP | F1 PB-only | F1 cascade | ΔF1 | ΔF1 95% CI | routed→PB | disagree | verdict |",
             "|---|---:|---:|---:|---:|---:|---|---:|---:|---|"]
    for name, pool in report["pooled"].items():
        ci = pool.get("vs_pb_paired_bootstrap", {}).get("delta_f1_ci95")
        ci_text = f"[{ci[0]:+.6f}, {ci[1]:+.6f}]" if ci else "—"
        lines.append(
            f"| `{name}` | {pool['loci']:,} | {pool['snp']:,} | {pool['pb_only']['f1']:.6f} | "
            f"{pool['cheap_router_pb']['f1']:.6f} | {pool['delta_f1_vs_pb']:+.6f} | {ci_text} | "
            f"{pool['fraction_routed_to_pb'] * 100:.4f}% | {pool['disagreement_vs_pb']} | "
            f"**{pool['verdict']}** |")
    return "\n".join(lines)


def compute_table(report: dict) -> str:
    lines = ["| Cell | loci | PB loci | PB compute fraction | PB-only s | cascade s | speedup |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for base in ORDER:
        for tag in TAGS:
            cell = report["cells"].get(f"{base}_{tag}")
            if cell is None or "projected_wallclock" not in cell:
                continue
            w = cell["projected_wallclock"]
            lines.append(
                f"| `{base}` {tag} | {w['loci']:,} | {w['pb_loci']:,} | "
                f"{w['pb_compute_fraction'] * 100:.4f}% | {w['pb_only_seconds']:.1f} | "
                f"{w['cascade_seconds']:.1f} | {w['speedup']:.1f}x |")
    return "\n".join(lines)


def curve_table(report: dict, key: str) -> str:
    cell = report["cells"][key]
    lines = [f"Accuracy-vs-compute, `{key}` (diagnostic only; the cutoff stays frozen)", "",
             "| cutoff | routed→PB | F1 | ΔF1 vs PB | disagreements | projected speedup |",
             "|---:|---:|---:|---:|---:|---:|"]
    for row in cell["accuracy_compute_curve"]:
        marker = " **←frozen**" if row["is_frozen_operating_point"] else ""
        speed = f"{row['projected_speedup']:.1f}x" if "projected_speedup" in row else "—"
        lines.append(
            f"| {row['cutoff']:.3f}{marker} | {row['fraction_routed_to_pb'] * 100:.4f}% | "
            f"{row['f1']:.6f} | {row['delta_f1_vs_pb']:+.6f} | {row['disagreement_vs_pb']} | {speed} |")
    return "\n".join(lines)


def depth_table(report: dict) -> str:
    lines = ["| Cell | depth bin | loci | SNP | F1 PB | F1 cascade | ΔF1 | routed→PB | disagree |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for base in ORDER:
        for tag in TAGS:
            cell = report["cells"].get(f"{base}_{tag}")
            if cell is None:
                continue
            for name, block in cell["by_depth_fine"].items():
                if not block.get("loci"):
                    continue
                lines.append(
                    f"| `{base}` {tag} | {name} | {block['loci']:,} | {block['snp']} | "
                    f"{block['pb_only']['f1']:.6f} | {block['cheap_router_pb']['f1']:.6f} | "
                    f"{block['delta_f1_vs_pb']:+.6f} | "
                    f"{block['fraction_routed_to_pb'] * 100:.4f}% | {block['disagreement_vs_pb']} |")
    return "\n".join(lines)


def structure_table(report: dict) -> str:
    lines = ["| Cell | stratum | loci | SNP | F1 PB | F1 cascade | ΔF1 | routed→PB | disagree |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for base in ORDER:
        for tag in TAGS:
            cell = report["cells"].get(f"{base}_{tag}")
            if cell is None:
                continue
            for name, block in cell["by_structure"].items():
                if not block.get("loci"):
                    continue
                lines.append(
                    f"| `{base}` {tag} | {name} | {block['loci']:,} | {block['snp']} | "
                    f"{block['pb_only']['f1']:.6f} | {block['cheap_router_pb']['f1']:.6f} | "
                    f"{block['delta_f1_vs_pb']:+.6f} | "
                    f"{block['fraction_routed_to_pb'] * 100:.4f}% | {block['disagreement_vs_pb']} |")
    return "\n".join(lines)


def controls_table(report: dict) -> str:
    lines = ["| Cell | PB budget | frozen router ΔF1 | random routing ΔF1 (mean±sd) | reverse-confidence ΔF1 |",
             "|---|---:|---:|---:|---:|"]
    for base in ORDER:
        for tag in TAGS:
            cell = report["cells"].get(f"{base}_{tag}")
            if cell is None:
                continue
            c = cell["controls"]
            random_deltas = [r["delta_f1_vs_pb"] for r in c["random_routing"]["runs"]]
            mean = sum(random_deltas) / len(random_deltas)
            sd = (sum((d - mean) ** 2 for d in random_deltas) / len(random_deltas)) ** 0.5
            lines.append(
                f"| `{base}` {tag} | {c['matched_pb_coverage']:,} | "
                f"{c['frozen_router']['delta_f1_vs_pb']:+.6f} | {mean:+.6f} ± {sd:.6f} | "
                f"{c['reverse_confidence_routing']['delta_f1_vs_pb']:+.6f} |")
    return "\n".join(lines)


def failure_table(report: dict) -> str:
    boundary = report["failure_boundary"]
    out = []
    for axis in ("by_depth", "by_quality", "by_vaf", "by_truth"):
        out.append(f"\n**{axis}**\n")
        out.append("| bin | loci | disagreements | rate |")
        out.append("|---|---:|---:|---:|")
        for name, cell in boundary[axis].items():
            rate = "—" if cell["rate"] is None else f"{cell['rate']:.3e}"
            out.append(f"| {name} | {cell['loci']:,} | {cell['disagreements']} | {rate} |")
    total = boundary["total"]
    out.append(f"\nTotal: {total['disagreements']} disagreements in {total['loci']:,} loci "
               f"(rate {total['rate']:.3e}).")
    return "\n".join(out)


def disagreement_effects() -> str:
    if not DISAGREEMENTS.exists():
        return "(no disagreement dump)"
    data = json.loads(DISAGREEMENTS.read_text())
    tally: dict[str, int] = {}
    for entry in data.values():
        for record in entry["records"]:
            tally[record["effect"]] = tally.get(record["effect"], 0) + 1
    if not tally:
        return "No disagreement between the frozen cascade and PB-only anywhere."
    lines = ["| effect | count |", "|---|---:|"]
    for name, count in sorted(tally.items(), key=lambda kv: -kv[1]):
        lines.append(f"| `{name}` | {count} |")
    return "\n".join(lines)


def main() -> None:
    report = json.loads(RESULTS.read_text())
    print("### Main result table\n")
    print(main_table(report))
    print("\n### Pooled\n")
    print(pooled_table(report))
    print("\n### All arms\n")
    print(arms_table(report))
    print("\n### Compute\n")
    print(compute_table(report))
    print("\n### Depth-stratified\n")
    print(depth_table(report))
    print("\n### Structure-stratified\n")
    print(structure_table(report))
    print("\n### Controls\n")
    print(controls_table(report))
    print("\n### Failure boundary\n")
    print(failure_table(report))
    print("\n### Disagreement effects\n")
    print(disagreement_effects())
    for key in ("chr20_neutral_full", "chr1_segdup_full", "chr19_gcrich_15x"):
        if key in report["cells"]:
            print("\n### " + key + "\n")
            print(curve_table(report, key))


if __name__ == "__main__":
    main()
