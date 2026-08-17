"""bench_v14 — mechanical, results-blind selection of the benchmark windows.

This script runs *before* any caller is executed, and it never reads a caller
score, a router decision or a benchmark metric. Its only inputs are

* the GIAB v4.2.1 high-confidence BED (which loci carry usable truth),
* the GIAB v3.1 genome stratifications (which loci are annotated difficult),
* the reference FASTA (GC content, N content).

The selection rules are fixed in ``History/14_DEVLOG.md`` §5 and implemented
verbatim below, so the chosen coordinates are a deterministic function of public
annotation rather than of anything this benchmark measures.

All three "ordinary" rules share a callability filter, **F**: 1 Mb windows at
whole-Mb offsets with ``hc_fraction >= 0.90``, ``difficult_fraction <= 0.15``
and ``n_fraction == 0``.

Rule O (neutral window, chr20)
    Among F, take the window whose start is closest to the chromosome midpoint;
    ties broken by the smaller coordinate.

Rule G+ (GC-rich window, chr19)
    Among windows with ``hc_fraction >= 0.90`` and ``n_fraction == 0`` — filter F
    **without** its difficulty cap — take the window with the highest GC
    fraction. The cap is dropped for this arm alone because GC-rich human
    sequence is Alu-dense by construction: under the full filter F the most
    GC-rich survivor on chr19 is 45.3%, below chr20's neutral window, so the cap
    silently deletes the very regime this arm exists to test. The selected
    window is repeat-rich but *not* segdup-rich (1.8% of its HC bases are
    low-mappability/segdup), so it stays a distinct axis from the chr1 arm.

Rule G- (AT-rich window, chr4)
    Among F, take the window with the lowest GC fraction.

Rule H (hard window, chr1)
    Among all 1 Mb windows at whole-Mb offsets with ``hc_fraction >= 0.50``
    and ``n_fraction == 0``, take the one maximising the fraction of
    high-confidence bases that fall inside
    ``GRCh38_alllowmapandsegdupregions``; ties broken by the smaller coordinate.

Revision note, recorded because it matters for pre-registration hygiene: Rules
G+/G- replaced an earlier "closest to the chromosome midpoint" rule applied to
all three ordinary contigs. That rule put the chr19 window at 28–29 Mb, which is
pericentromeric and GC 43.0% — the opposite of the GC-rich regime chr19 was
chosen to represent, so the rule failed to express the design. The revision was
made **before any caller, router or metric was run on any of these windows**;
the only outputs in existence at that moment were the annotation statistics
computed by this script. No benchmark result informed it. See devlog 14 §5.

No rule can see chr21, so no previously used coordinate can be selected.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pysam

WINDOW = 1_000_000

#: (contig, reference FASTA, selection rule).
TARGETS = [
    ("chr20", "data/reference/chr20_full.fa", "neutral"),
    ("chr19", "data/reference/chr19_full.fa", "gc_high"),
    ("chr4", "data/reference/chr4_full.fa", "gc_low"),
    ("chr1", "data/reference/chr1_full.fa", "hard"),
]

HC_BED = "data/giab_hg002_v14/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
STRAT = {
    "lowmap_segdup": "data/strat_v31/lowmap_segdup.bed.gz",
    "alldifficult": "data/strat_v31/alldifficult.bed.gz",
    "tandemrepeats": "data/strat_v31/tandemrepeats.bed.gz",
}
OUT = Path("results/bench_v14/region_selection.json")


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def coverage_mask(bed_path: str, contig: str, length: int) -> np.ndarray:
    """Boolean ``[length]`` mask of bases covered by ``contig``'s BED intervals."""
    mask = np.zeros(length, dtype=bool)
    wanted = {contig, contig[3:] if contig.startswith("chr") else "chr" + contig}
    with _open(bed_path) as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if fields[0] not in wanted:
                continue
            start, stop = int(fields[1]), int(fields[2])
            if start >= length:
                continue
            mask[start:min(stop, length)] = True
    return mask


def characterise(contig: str, fasta_path: str) -> list[dict]:
    """Per-window annotation statistics for one contig. No caller is involved."""
    fasta = pysam.FastaFile(fasta_path)
    name = fasta.references[0]
    length = fasta.get_reference_length(name)
    hc = coverage_mask(HC_BED, contig, length)
    strat = {key: coverage_mask(path, contig, length) for key, path in STRAT.items()}

    rows = []
    for start in range(0, length - WINDOW + 1, WINDOW):
        stop = start + WINDOW
        sequence = np.frombuffer(fasta.fetch(name, start, stop).upper().encode(), dtype="S1")
        n_frac = float(np.mean(sequence == b"N"))
        gc = float(np.mean(np.isin(sequence, [b"G", b"C"])) / max(1e-12, 1.0 - n_frac))
        window_hc = hc[start:stop]
        hc_frac = float(window_hc.mean())
        row = {
            "contig": contig, "start": start, "stop": stop,
            "hc_fraction": hc_frac, "n_fraction": n_frac, "gc_fraction": gc,
        }
        for key, mask in strat.items():
            window_mask = mask[start:stop]
            row[f"{key}_fraction"] = float(window_mask.mean())
            row[f"{key}_in_hc"] = float((window_mask & window_hc).sum() / max(1, window_hc.sum()))
        row["difficult_fraction"] = row["alldifficult_fraction"]
        rows.append(row)
    fasta.close()
    return rows


def apply_rule(rows: list[dict], rule: str, midpoint: int) -> dict:
    """Rule O / Rule H from the module docstring."""
    if rule in ("neutral", "gc_high", "gc_low"):
        # Rule G+ drops the difficulty cap; see the module docstring.
        cap = 1.01 if rule == "gc_high" else 0.15
        pool = [r for r in rows if r["hc_fraction"] >= 0.90
                and r["difficult_fraction"] <= cap and r["n_fraction"] == 0.0]
        assert pool, f"no window satisfies filter F for rule {rule}"
        if rule == "neutral":
            return min(pool, key=lambda r: (abs(r["start"] - midpoint), r["start"]))
        sign = -1.0 if rule == "gc_high" else 1.0
        return min(pool, key=lambda r: (sign * r["gc_fraction"], r["start"]))
    pool = [r for r in rows if r["hc_fraction"] >= 0.50 and r["n_fraction"] == 0.0]
    assert pool, "no window satisfies Rule H"
    return min(pool, key=lambda r: (-r["lowmap_segdup_in_hc"], r["start"]))


def main() -> None:
    selection = {}
    for contig, fasta_path, rule in TARGETS:
        rows = characterise(contig, fasta_path)
        length = max(r["stop"] for r in rows)
        chosen = apply_rule(rows, rule, midpoint=length // 2)
        chosen["rule"] = rule
        chosen["n_candidate_windows"] = len(rows)
        selection[contig] = chosen
        print(f"{contig} [{rule}] -> {chosen['start']:,}-{chosen['stop']:,} "
              f"hc={chosen['hc_fraction']:.3f} gc={chosen['gc_fraction']:.3f} "
              f"difficult={chosen['difficult_fraction']:.3f} "
              f"segdup_in_hc={chosen['lowmap_segdup_in_hc']:.3f} "
              f"tr_in_hc={chosen['tandemrepeats_in_hc']:.3f}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(selection, indent=2))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
