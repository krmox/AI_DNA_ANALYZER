"""bench_v19 (Validation Round 2) — mechanical, results-blind selection of
three new benchmark windows on three contigs no prior devlog has touched.

Runs before any caller, router, or metric exists for HG004. Inputs are public
annotation only: HG004's own GIAB v4.2.1 high-confidence BED, the GIAB v3.1
stratification BEDs, and the reference FASTA. No caller score, router
decision or benchmark metric is read by this script.

Rule C19 (which contig gets which role) — mechanical, not annotation-driven:
untouched GRCh38 autosomes (every autosome except chr21, chr20, chr19, chr4,
chr1, chr16, chr15, chr7, chr14, chr13, chr17 — used by devlogs 11-18), sorted
ascending by chromosome number: chr2, chr3, chr5, chr6, chr8, chr9, chr10,
chr11, chr12, chr18, chr22. The first three in that fixed order are assigned,
in the order the three roles are listed in the Round-2 pre-registration
(ordinary, difficult, stress): chr2 -> ordinary, chr3 -> difficult,
chr5 -> stress. Role assignment reads no annotation of any kind.

Rule R19-ordinary (chr2, "typical difficulty") — verbatim devlog 16 Rule R16:
among 3 Mb windows at whole-Mb offsets with hc_fraction >= 0.50 and
n_fraction == 0, take the window whose alldifficult_hc_fraction is closest to
the contig-wide alldifficult_hc_fraction. Ties -> smaller coordinate.

Rule R19-difficult (chr3, "hardest scorable") — verbatim devlog 17 Rule R17:
among the same window pool, take the window maximising
lowmap_segdup_hc_fraction. Ties -> smaller coordinate.

Rule R19-stress (chr5, "worst-case stress") — NEW axis: among the same window
pool, take the window maximising tandemrepeats_hc_fraction. Every prior devlog
used tandemrepeats only as a stratification axis, never as the primary
selection criterion; this rule exists so Round 2 exercises a genome-context
failure surface (paralogous/repetitive low-complexity sequence) the project
has never deliberately selected FOR before, distinct from the
already-characterised lowmap_segdup stress case.

Window width, 3 Mb: matches devlog 16/17/18 for direct statistical
comparability (paired bootstrap CI half-width ~0.0015 at ~4,000+ SNP).
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pysam

WINDOW = 3_000_000

TARGETS = [
    ("chr2", "data/reference/chr2_full.fa", "ordinary"),
    ("chr3", "data/reference/chr3_full.fa", "difficult"),
    ("chr5", "data/reference/chr5_full.fa", "stress"),
]

PREVIOUSLY_USED_CONTIGS = {"chr21", "chr20", "chr19", "chr4", "chr1",
                          "chr16", "chr15", "chr7", "chr14", "chr13", "chr17"}

HC_BED = "data/giab_hg004_v19/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
STRAT = {
    "lowmap_segdup": "data/strat_v31/lowmap_segdup.bed.gz",
    "alldifficult": "data/strat_v31/alldifficult.bed.gz",
    "tandemrepeats": "data/strat_v31/tandemrepeats.bed.gz",
}
OUT = Path("results/bench_v19/region_selection.json")


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def coverage_mask(bed_path: str, contig: str, length: int) -> np.ndarray:
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
    fasta = pysam.FastaFile(fasta_path)
    name = fasta.references[0]
    length = fasta.get_reference_length(name)
    hc = coverage_mask(HC_BED, contig, length)
    strat = {key: coverage_mask(path, contig, length) for key, path in STRAT.items()}
    contig_hc_bases = max(1, hc.sum())
    contig_alldifficult_hc = float((strat["alldifficult"] & hc).sum() / contig_hc_bases)

    rows = []
    for start in range(0, length - WINDOW + 1, WINDOW):
        stop = start + WINDOW
        sequence = np.frombuffer(fasta.fetch(name, start, stop).upper().encode(), dtype="S1")
        n_frac = float(np.mean(sequence == b"N"))
        gc = float(np.mean(np.isin(sequence, [b"G", b"C"])) / max(1e-12, 1.0 - n_frac))
        window_hc = hc[start:stop]
        hc_frac = float(window_hc.mean())
        row = {"contig": contig, "start": start, "stop": stop,
               "hc_fraction": hc_frac, "n_fraction": n_frac, "gc_fraction": gc}
        for key, mask in strat.items():
            window_mask = mask[start:stop]
            row[f"{key}_fraction"] = float(window_mask.mean())
            row[f"{key}_hc_fraction"] = float((window_mask & window_hc).sum() / max(1, window_hc.sum()))
        rows.append(row)
    fasta.close()
    return rows, contig_alldifficult_hc


def apply_rule(rows: list[dict], role: str, contig_alldifficult_hc: float) -> dict:
    pool = [r for r in rows if r["hc_fraction"] >= 0.50 and r["n_fraction"] == 0.0]
    assert pool, f"no window satisfies the callability floor for role {role}"
    if role == "ordinary":
        return min(pool, key=lambda r: (abs(r["alldifficult_hc_fraction"] - contig_alldifficult_hc), r["start"]))
    if role == "difficult":
        return min(pool, key=lambda r: (-r["lowmap_segdup_hc_fraction"], r["start"]))
    if role == "stress":
        return min(pool, key=lambda r: (-r["tandemrepeats_hc_fraction"], r["start"]))
    raise ValueError(role)


def main() -> None:
    selection = {}
    for contig, fasta_path, role in TARGETS:
        assert contig not in PREVIOUSLY_USED_CONTIGS, contig
        rows, contig_alldifficult_hc = characterise(contig, fasta_path)
        chosen = apply_rule(rows, role, contig_alldifficult_hc)
        chosen["role"] = role
        chosen["n_candidate_windows"] = len(rows)
        chosen["contig_alldifficult_hc_fraction"] = contig_alldifficult_hc
        selection[contig] = chosen
        print(f"{contig} [{role}] -> {chosen['start']:,}-{chosen['stop']:,} "
              f"hc={chosen['hc_fraction']:.4f} gc={chosen['gc_fraction']:.4f} "
              f"alldifficult_hc={chosen['alldifficult_hc_fraction']:.4f} "
              f"segdup_hc={chosen['lowmap_segdup_hc_fraction']:.4f} "
              f"tr_hc={chosen['tandemrepeats_hc_fraction']:.4f} "
              f"(n_windows={chosen['n_candidate_windows']})")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(selection, indent=2))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
