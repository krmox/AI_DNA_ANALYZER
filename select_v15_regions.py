"""bench_v15 — mechanical, results-blind selection of the *test* windows.

This script runs before any v15 caller is executed and it never reads a caller
score, a router decision, a safety-layer decision or a benchmark metric. Its
only inputs are

* the GIAB v4.2.1 high-confidence BED (which loci carry usable truth),
* the GIAB v3.1 genome stratifications (which loci are annotated difficult),
* the reference FASTA (GC content, N content).

It is a direct descendant of ``select_v14_regions.py``: the window rules are
that script's Rule H and Rule O, quoted below unchanged. What is new is a
*contig* rule, because devlog 15 needs test chromosomes that no previous devlog
has touched and it must pick them without looking at a result either.

Rule C (which contigs)
    Consider the GRCh38 autosomes ``chr2 ... chr18, chr22`` — every autosome
    except those any previous experiment used (``chr21`` for devlogs 3–13;
    ``chr20``, ``chr19``, ``chr4``, ``chr1`` for devlog 14) . Rank them by

        segdup_hc_fraction = |HC ∩ lowmap_segdup| / |HC|

    computed over the whole contig from the two public BED files. Take the
    **top three** as the hard test contigs and the **bottom one** as the
    ordinary-sequence control contig. Ties broken by the smaller contig number.

Rule H (hard window), applied to each of the three hard contigs
    Among all 1 Mb windows at whole-Mb offsets with ``hc_fraction >= 0.50`` and
    ``n_fraction == 0``, take the one maximising the fraction of
    high-confidence bases inside ``GRCh38_alllowmapandsegdupregions``; ties
    broken by the smaller coordinate. Verbatim from devlog 14 §5.2.

Rule O (neutral window), applied to the control contig
    Among 1 Mb windows at whole-Mb offsets with ``hc_fraction >= 0.90``,
    ``difficult_fraction <= 0.15`` and ``n_fraction == 0``, take the window
    whose start is closest to the contig midpoint; ties broken by the smaller
    coordinate. Verbatim from devlog 14 §5.2.

Rule H is applied to three contigs rather than one because devlog 14 found the
frozen cascade's only true-variant losses inside segmental duplication, and a
test set that cannot contain that failure mode cannot test a repair for it. The
control window exists so the safety layer's false-escalation cost is measured on
ordinary sequence, which is where almost all of the genome lives.

No rule can select chr21, chr20, chr19, chr4 or chr1, so no previously used
coordinate can be selected.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pysam

WINDOW = 1_000_000

#: Contigs any previous experiment used. Rule C may not consider them.
PREVIOUSLY_USED_CONTIGS = ("chr21", "chr20", "chr19", "chr4", "chr1")

#: The autosomes Rule C ranks: every autosome except the ones already used.
CANDIDATE_CONTIGS = tuple(c for c in
                          tuple(f"chr{i}" for i in range(2, 19)) + ("chr22",)
                          if c not in PREVIOUSLY_USED_CONTIGS)

HC_BED = "data/giab_hg002_v14/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
STRAT = {
    "lowmap_segdup": "data/strat_v31/lowmap_segdup.bed.gz",
    "alldifficult": "data/strat_v31/alldifficult.bed.gz",
    "tandemrepeats": "data/strat_v31/tandemrepeats.bed.gz",
}
FASTA_TEMPLATE = "data/reference/{contig}_full.fa"
OUT = Path("results/bench_v15/region_selection.json")

#: GRCh38 autosome lengths, so contig-level fractions can be computed from the
#: BED files alone -- before any reference FASTA is downloaded.
CONTIG_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468,
}


def _names(contig: str) -> set[str]:
    """Both spellings of a contig, so BEDs with or without ``chr`` both match."""
    return {contig, contig[3:] if contig.startswith("chr") else "chr" + contig}


def bed_mask(bed_path: str, contig: str, start: int, stop: int) -> np.ndarray:
    """Boolean mask over ``[start, stop)`` of bases covered by a BED's intervals."""
    mask = np.zeros(stop - start, dtype=bool)
    wanted = _names(contig)
    opener = gzip.open if bed_path.endswith(".gz") else open
    with opener(bed_path, "rt") as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if fields[0] not in wanted:
                continue
            lo, hi = int(fields[1]), int(fields[2])
            if hi <= start or lo >= stop:
                continue
            mask[max(lo, start) - start:min(hi, stop) - start] = True
    return mask


def rank_contigs() -> list[dict]:
    """Rule C: rank the candidate autosomes by segdup fraction of their HC bases.

    Uses only the two public BED files, so this runs before any reference FASTA
    is fetched and cannot be informed by anything measured later.
    """
    rows = []
    for contig in CANDIDATE_CONTIGS:
        length = CONTIG_LENGTHS[contig]
        hc = bed_mask(HC_BED, contig, 0, length)
        segdup = bed_mask(STRAT["lowmap_segdup"], contig, 0, length)
        hc_bases = int(hc.sum())
        rows.append({
            "contig": contig,
            "hc_bases": hc_bases,
            "hc_fraction": hc_bases / length,
            "segdup_hc_bases": int((hc & segdup).sum()),
            "segdup_hc_fraction": float((hc & segdup).sum() / max(hc_bases, 1)),
        })
    order = sorted(rows, key=lambda r: (-r["segdup_hc_fraction"],
                                        int(r["contig"][3:])))
    for rank, row in enumerate(order):
        row["rank_by_segdup"] = rank
    return rows


def select_contigs(rows: list[dict]) -> dict:
    """Top three by Rule C are the hard contigs; the bottom one is the control."""
    order = sorted(rows, key=lambda r: (-r["segdup_hc_fraction"], int(r["contig"][3:])))
    return {"hard": [r["contig"] for r in order[:3]], "neutral": order[-1]["contig"]}


def window_stats(contig: str) -> list[dict]:
    """Annotation statistics for every 1 Mb window at whole-Mb offsets."""
    fasta_path = FASTA_TEMPLATE.format(contig=contig)
    handle = pysam.FastaFile(fasta_path)
    name = next(iter(_names(contig) & set(handle.references)))
    length = handle.get_reference_length(name)
    hc = bed_mask(HC_BED, contig, 0, length)
    strat = {key: bed_mask(path, contig, 0, length) for key, path in STRAT.items()}

    rows = []
    for start in range(0, length - WINDOW + 1, WINDOW):
        stop = start + WINDOW
        sequence = np.frombuffer(handle.fetch(name, start, stop).upper().encode(), dtype="S1")
        n_fraction = float(np.mean(sequence == b"N"))
        gc = float(np.mean((sequence == b"G") | (sequence == b"C")))
        hc_window = hc[start:stop]
        hc_fraction = float(hc_window.mean())
        row = {"contig": contig, "start": start, "stop": stop, "gc": gc,
               "n_fraction": n_fraction, "hc_fraction": hc_fraction}
        for key, mask in strat.items():
            row[f"{key}_fraction"] = float(mask[start:stop].mean())
            row[f"{key}_hc_fraction"] = float(
                (mask[start:stop] & hc_window).sum() / max(hc_window.sum(), 1))
        rows.append(row)
    handle.close()
    return rows


def pick_hard(rows: list[dict]) -> dict:
    """Rule H, verbatim from devlog 14 §5.2."""
    eligible = [r for r in rows if r["hc_fraction"] >= 0.50 and r["n_fraction"] == 0.0]
    if not eligible:
        raise SystemExit(f"Rule H found no eligible window on {rows[0]['contig']}")
    return min(eligible, key=lambda r: (-r["lowmap_segdup_hc_fraction"], r["start"]))


def pick_neutral(rows: list[dict]) -> dict:
    """Rule O, verbatim from devlog 14 §5.2."""
    eligible = [r for r in rows if r["hc_fraction"] >= 0.90
                and r["alldifficult_fraction"] <= 0.15 and r["n_fraction"] == 0.0]
    if not eligible:
        raise SystemExit(f"Rule O found no eligible window on {rows[0]['contig']}")
    midpoint = rows[-1]["stop"] / 2.0
    return min(eligible, key=lambda r: (abs(r["start"] - midpoint), r["start"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs-only", action="store_true",
                        help="run Rule C only (no FASTA needed) and stop")
    args = parser.parse_args()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    ranking = rank_contigs()
    chosen = select_contigs(ranking)
    report = {
        "previously_used_contigs": list(PREVIOUSLY_USED_CONTIGS),
        "candidate_contigs": list(CANDIDATE_CONTIGS),
        "contig_ranking": sorted(ranking, key=lambda r: r["rank_by_segdup"]),
        "chosen_contigs": chosen,
    }
    if args.contigs_only:
        print(json.dumps(report, indent=2))
        OUT.with_name("contig_selection.json").write_text(json.dumps(report, indent=2))
        return

    selected = {}
    for contig in chosen["hard"]:
        rows = window_stats(contig)
        selected[f"{contig}_segdup"] = {**pick_hard(rows), "rule": "H",
                                        "regime": "segmental duplication / low mappability"}
    rows = window_stats(chosen["neutral"])
    selected[f"{chosen['neutral']}_neutral"] = {**pick_neutral(rows), "rule": "O",
                                                "regime": "ordinary euchromatin (control)"}
    report["selected"] = selected
    OUT.write_text(json.dumps(report, indent=2))
    print(json.dumps(selected, indent=2))


if __name__ == "__main__":
    main()
