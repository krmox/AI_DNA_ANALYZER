"""bench_v17 -- mechanical, results-blind selection of ONE hard unseen test region.

Runs before any v17 caller executes and never reads a caller score, a router
decision or a benchmark metric. Its only inputs are the GIAB v4.2.1
high-confidence BED, the GIAB v3.1 stratifications, and (only after Rule C17 has
already committed to a contig) the reference FASTA for GC/N content.

Direct descendant of select_v14_regions.py / select_v15_regions.py /
select_v16_region.py; BED machinery reused unchanged.

Two rules, both fixed and written down before execution:

Rule C17 (which contig)
    Consider the GRCh38 autosomes no previous experiment has touched -- every
    autosome except chr21 (devlogs 3-13), chr20/chr19/chr4/chr1 (devlog 14),
    chr16/chr15/chr7/chr14 (devlog 15) and chr13 (devlog 16). Rank the survivors
    by segdup_hc_fraction = |HC intersect lowmap_segdup| / |HC| over the whole
    contig, descending, and take rank 0 -- the HARDEST remaining contig by
    segmental-duplication burden. Ties broken by smaller contig number.

    Rationale, recorded before the rule runs: this experiment is explicitly a
    hard-data stress test (devlog 16 took the median contig; devlog 14's control
    took an easy window). Segmental duplication is the one regime where devlog 14
    found the cascade's only known failure mode, so the remaining contig with the
    largest segdup burden is the honest place to look for it again. The ranking
    is fully determined by two public BED files and cannot be gamed.

Rule R17 (which window, "hardest")
    Among all 3 Mb windows at whole-Mb offsets on the chosen contig with
    hc_fraction >= 0.50 and n_fraction == 0, take the window MAXIMISING
    segdup_hc_fraction. Ties broken by the smaller coordinate.

    hc_fraction >= 0.50 exists only so the window carries enough truth to score
    and is verbatim from devlog 15's Rule H and devlog 16's Rule R16. Width 3 Mb
    is fixed here for statistical power, copied from devlog 16, not chosen from
    any property of this data.

No rule can select chr21, chr20, chr19, chr4, chr1, chr16, chr15, chr7, chr14 or
chr13, so no previously used coordinate can be selected.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pysam

from select_v15_regions import CONTIG_LENGTHS, HC_BED, STRAT, bed_mask

WINDOW = 3_000_000

PREVIOUSLY_USED_CONTIGS = ("chr21", "chr20", "chr19", "chr4", "chr1",
                           "chr16", "chr15", "chr7", "chr14", "chr13")

CANDIDATE_CONTIGS = tuple(c for c in (f"chr{i}" for i in range(1, 23))
                          if c not in PREVIOUSLY_USED_CONTIGS)

FASTA_TEMPLATE = "data/reference/{contig}_full.fa"
OUT = Path("results/bench_v17/region_selection.json")


def rank_contigs() -> list[dict]:
    """Rule C17 ranking. Uses only the two public BED files -- no FASTA needed."""
    rows = []
    for contig in CANDIDATE_CONTIGS:
        length = CONTIG_LENGTHS[contig]
        hc = bed_mask(HC_BED, contig, 0, length)
        segdup = bed_mask(STRAT["lowmap_segdup"], contig, 0, length)
        difficult = bed_mask(STRAT["alldifficult"], contig, 0, length)
        hc_bases = int(hc.sum())
        rows.append({
            "contig": contig,
            "length": length,
            "hc_bases": hc_bases,
            "hc_fraction": hc_bases / length,
            "segdup_hc_fraction": float((hc & segdup).sum() / max(hc_bases, 1)),
            "alldifficult_hc_fraction": float((hc & difficult).sum() / max(hc_bases, 1)),
        })
    order = sorted(rows, key=lambda r: (-r["segdup_hc_fraction"], int(r["contig"][3:])))
    for rank, row in enumerate(order):
        row["rank_by_segdup"] = rank
    return order


def select_contig(order: list[dict]) -> dict:
    """Rule C17: rank 0 -- the hardest remaining contig."""
    return order[0]


def window_stats(contig: str) -> list[dict]:
    """Annotation statistics for every 3 Mb window at whole-Mb offsets."""
    handle = pysam.FastaFile(FASTA_TEMPLATE.format(contig=contig))
    names = {contig, contig[3:]}
    name = next(iter(names & set(handle.references)))
    length = handle.get_reference_length(name)
    hc = bed_mask(HC_BED, contig, 0, length)
    strat = {key: bed_mask(path, contig, 0, length) for key, path in STRAT.items()}

    rows = []
    for start in range(0, length - WINDOW + 1, 1_000_000):
        stop = start + WINDOW
        sequence = np.frombuffer(handle.fetch(name, start, stop).upper().encode(), dtype="S1")
        hc_win = hc[start:stop]
        hc_bases = int(hc_win.sum())
        row = {
            "contig": contig, "start": start, "stop": stop,
            "hc_fraction": hc_bases / WINDOW,
            "n_fraction": float((sequence == b"N").sum() / WINDOW),
            "gc_fraction": float(np.isin(sequence, [b"G", b"C"]).sum() / WINDOW),
        }
        for key, mask in strat.items():
            win = mask[start:stop]
            row[f"{key}_fraction"] = float(win.sum() / WINDOW)
            row[f"{key}_hc_fraction"] = float((win & hc_win).sum() / max(hc_bases, 1))
        rows.append(row)
    handle.close()
    return rows


def select_window(rows: list[dict]) -> dict:
    """Rule R17: eligible window maximising segdup_hc_fraction; ties -> lower start."""
    eligible = [r for r in rows if r["hc_fraction"] >= 0.50 and r["n_fraction"] == 0.0]
    if not eligible:
        raise SystemExit("no eligible window under Rule R17")
    return sorted(eligible, key=lambda r: (-r["lowmap_segdup_hc_fraction"], r["start"]))[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("contig", "window"), default="window")
    args = parser.parse_args()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    order = rank_contigs()
    chosen = select_contig(order)
    report = {
        "rule_c17": "rank-0 contig of the descending segdup_hc_fraction ranking (hardest remaining)",
        "rule_r17": ("3 Mb window, hc_fraction >= 0.50, n_fraction == 0, "
                     "maximising lowmap_segdup_hc_fraction; ties -> lower start"),
        "window_bp": WINDOW,
        "previously_used_contigs": list(PREVIOUSLY_USED_CONTIGS),
        "candidate_contigs": list(CANDIDATE_CONTIGS),
        "contig_ranking": order,
        "selected_contig": chosen,
    }
    if args.stage == "window":
        rows = window_stats(chosen["contig"])
        report["n_windows"] = len(rows)
        report["n_eligible_windows"] = sum(
            1 for r in rows if r["hc_fraction"] >= 0.50 and r["n_fraction"] == 0.0)
        report["selected_window"] = select_window(rows)
        report["window_stats"] = rows
    OUT.write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ("contig_ranking", "window_stats")}, indent=2))


if __name__ == "__main__":
    main()
