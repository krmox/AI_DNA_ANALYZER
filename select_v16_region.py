"""bench_v16 — mechanical, results-blind selection of the ONE unseen test region.

This script runs before any v16 caller is executed and it never reads a caller
score, a router decision or a benchmark metric. Its only inputs are

* the GIAB v4.2.1 high-confidence BED (which loci carry usable truth),
* the GIAB v3.1 genome stratifications (which loci are annotated difficult),
* the reference FASTA (GC content, N content), fetched only *after* Rule C16 has
  already committed to a contig.

It is a direct descendant of ``select_v14_regions.py`` / ``select_v15_regions.py``
and reuses their BED machinery unchanged. Two rules, both fixed before execution:

Rule C16 (which contig)
    Consider the GRCh38 autosomes no previous experiment has touched —
    every autosome except ``chr21`` (devlogs 3-13), ``chr20``, ``chr19``,
    ``chr4``, ``chr1`` (devlog 14) and ``chr16``, ``chr15``, ``chr7``, ``chr14``
    (devlog 15). Rank the survivors by

        segdup_hc_fraction = |HC ∩ lowmap_segdup| / |HC|

    over the whole contig, descending, and take the **median-rank** contig
    (index ``len // 2`` of the descending ranking; for an even-sized candidate
    set the later of the two central ranks). Ties broken by smaller contig
    number.

    Rationale, recorded before the rule runs: devlog 15 deliberately took the
    three *hardest* remaining contigs, and devlog 14's neutral control took an
    easy window. Neither is representative. The median of the remaining
    candidates is the contig whose segmental-duplication burden is typical of
    what is left, and it cannot be gamed because the ranking is fully determined
    by two public BED files.

Rule R16 (which window, "representative")
    Among all 3 Mb windows at whole-Mb offsets on the chosen contig with
    ``hc_fraction >= 0.50`` and ``n_fraction == 0``, take the window whose
    ``alldifficult_hc_fraction`` is **closest to the contig-wide**
    ``alldifficult_hc_fraction``. Ties broken by the smaller coordinate.

    Rationale: this selects a window of typical difficulty for the contig
    rather than an easy one (devlog 14 Rule O required ``difficult <= 0.15``,
    which is explicitly easy) or a maximally hard one (devlog 15 Rule H, which
    maximised segdup). The ``hc_fraction >= 0.50`` floor exists only so the
    window carries enough truth to score, and is verbatim from Rule H.

Window width, 3 Mb, is fixed here for statistical power and not for any property
of the data: devlog 14's 1 Mb cells carried ~1,000 SNP and gave a paired
bootstrap ΔF1 half-width of ~0.0025, so 3 Mb should give ~0.0015, comfortably
inside the pre-registered ``<= 0.01`` power criterion.

No rule can select chr21, chr20, chr19, chr4, chr1, chr16, chr15, chr7 or chr14,
so no previously used coordinate can be selected.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pysam

from select_v15_regions import CONTIG_LENGTHS, HC_BED, STRAT, bed_mask

WINDOW = 3_000_000

#: Contigs any previous experiment used. Rule C16 may not consider them.
PREVIOUSLY_USED_CONTIGS = ("chr21", "chr20", "chr19", "chr4", "chr1",
                           "chr16", "chr15", "chr7", "chr14")

#: The autosomes Rule C16 ranks: every autosome except the ones already used.
CANDIDATE_CONTIGS = tuple(c for c in (f"chr{i}" for i in range(1, 23))
                          if c not in PREVIOUSLY_USED_CONTIGS)

FASTA_TEMPLATE = "data/reference/{contig}_full.fa"
OUT = Path("results/bench_v16/region_selection.json")


def rank_contigs() -> list[dict]:
    """Rule C16 ranking. Uses only the two public BED files -- no FASTA needed."""
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
    """Rule C16: the median-rank contig of the descending segdup ranking."""
    return order[len(order) // 2]


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
        hc_window = hc[start:stop]
        row = {
            "contig": contig, "start": start, "stop": stop,
            "gc": float(np.mean((sequence == b"G") | (sequence == b"C"))),
            "n_fraction": float(np.mean(sequence == b"N")),
            "hc_fraction": float(hc_window.mean()),
        }
        for key, mask in strat.items():
            row[f"{key}_fraction"] = float(mask[start:stop].mean())
            row[f"{key}_hc_fraction"] = float(
                (mask[start:stop] & hc_window).sum() / max(hc_window.sum(), 1))
        rows.append(row)
    handle.close()
    return rows


def pick_representative(rows: list[dict], contig_difficult_hc: float) -> dict:
    """Rule R16, applied to the 3 Mb windows of the chosen contig."""
    eligible = [r for r in rows if r["hc_fraction"] >= 0.50 and r["n_fraction"] == 0.0]
    if not eligible:
        raise SystemExit(f"Rule R16 found no eligible window on {rows[0]['contig']}")
    for row in eligible:
        row["difficulty_gap"] = abs(row["alldifficult_hc_fraction"] - contig_difficult_hc)
    return min(eligible, key=lambda r: (r["difficulty_gap"], r["start"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contig-only", action="store_true",
                        help="run Rule C16 only (no FASTA needed) and stop")
    args = parser.parse_args()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    order = rank_contigs()
    chosen = select_contig(order)
    payload = {
        "rule_c16": "median-rank contig of the descending segdup_hc_fraction ranking",
        "rule_r16": ("3 Mb window, hc_fraction >= 0.50, n_fraction == 0, minimising "
                     "|alldifficult_hc_fraction - contig alldifficult_hc_fraction|"),
        "window_bp": WINDOW,
        "previously_used_contigs": list(PREVIOUSLY_USED_CONTIGS),
        "candidate_contigs": list(CANDIDATE_CONTIGS),
        "contig_ranking": order,
        "chosen_contig": chosen,
    }
    if args.contig_only:
        OUT.write_text(json.dumps(payload, indent=2))
        print(json.dumps({"chosen_contig": chosen["contig"],
                          "rank": chosen["rank_by_segdup"],
                          "of": len(order)}, indent=2))
        return

    rows = window_stats(chosen["contig"])
    window = pick_representative(rows, chosen["alldifficult_hc_fraction"])
    payload["windows_considered"] = len(rows)
    payload["chosen_window"] = window
    payload["all_windows"] = rows
    OUT.write_text(json.dumps(payload, indent=2))
    print(json.dumps({"contig": chosen["contig"], "window": window}, indent=2))


if __name__ == "__main__":
    main()
