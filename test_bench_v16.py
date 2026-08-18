"""Tests for the bench_v16 unseen-region benchmark driver.

The v14 suite already covers the shared helpers this driver imports. What is
tested here is only what v16 adds, and it is the set of things that could
silently corrupt *this* experiment:

1. **The frozen constants are read, never written**, and the pre-registered
   region is the one the selection script actually chose.
2. **The independence audit really fails** when a cell's contig, coordinates,
   locus count or LLRs are wrong — an audit that passes unconditionally is worse
   than no audit.
3. **The new mapping-quality stratifier and the disagreement enrichment index
   the scoring frame correctly**, since both mix frame-space and cache-space
   arrays and an off-by-one there would mislabel every record.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pytest

import bench_v16_unseen as bench
from cascade import (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD,
                     FROZEN_ROUTER_CUTOFF)
from config import LABEL_SNP


# -- frozen configuration and pre-registration --------------------------------


def test_frozen_constants_are_the_devlog_13_values():
    assert FROZEN_BINOMIAL_THRESHOLD == 7.0
    assert FROZEN_PB_THRESHOLD == 10.5
    assert FROZEN_ROUTER_CUTOFF == 5.411872376933351


def test_driver_never_assigns_to_a_frozen_constant():
    source = Path("bench_v16_unseen.py").read_text()
    for name in ("FROZEN_BINOMIAL_THRESHOLD", "FROZEN_PB_THRESHOLD", "FROZEN_ROUTER_CUTOFF"):
        assert f"{name} =" not in source


def test_region_matches_the_results_blind_selection():
    """The driver's hard-coded region must equal what select_v16_region.py chose."""
    import json
    chosen = json.loads(Path("results/bench_v16/region_selection.json").read_text())["chosen_window"]
    spec = bench.BASE_REGIONS["chr13_representative"]
    assert (spec["contig"], spec["start"], spec["stop"]) == (
        chosen["contig"], chosen["start"], chosen["stop"])
    assert spec["gc"] == pytest.approx(chosen["gc"])


def test_region_is_on_no_previously_used_contig():
    for spec in bench.BASE_REGIONS.values():
        assert spec["contig"] not in bench.PREVIOUSLY_USED_CONTIGS
    used = {contig for contig, _, _ in bench.HISTORICAL_REGIONS}
    assert used == bench.PREVIOUSLY_USED_CONTIGS


def test_three_cells_one_region_three_depths():
    specs = bench.cells()
    assert len(specs) == 3
    assert {s["depth_tag"] for s in specs.values()} == {"full", "30x", "15x"}


def test_no_routing_helper_can_see_truth():
    for function in (bench.genomic_positions, bench.measure_throughput):
        parameters = set(inspect.signature(function).parameters)
        assert "labels" not in parameters and "truth" not in parameters


# -- synthetic cell ------------------------------------------------------------


def _cell(n_loci: int = 128, contig: str = "chr13", start: int = 70_000_000):
    rng = np.random.default_rng(0)
    labels = np.zeros(n_loci, dtype=np.int64)
    labels[::16] = LABEL_SNP
    counts = np.zeros((n_loci, 10), dtype=np.float32)
    counts[:, 6] = 30.0
    counts[:, 0] = 20.0
    counts[:, 1] = 10.0
    counts[:, 9] = 30.0
    depth = counts[:, 6].astype(np.float64)
    return {
        "name": "synthetic",
        "counts": counts,
        "binomial_llr": rng.normal(7.0, 6.0, n_loci),
        "pb_llr": rng.normal(10.5, 6.0, n_loci),
        "labels": labels,
        "depth": depth,
        "mean_base_quality": np.full(n_loci, 32.0),
        "mean_mapq": np.linspace(0.0, 60.0, n_loci),
        "region": [start, start + n_loci],
        "contig": contig,
        "bam": "synthetic.bam",
        "positions": np.arange(start, start + n_loci, dtype=np.int64),
        "spec": {"contig": contig, "start": start, "stop": start + n_loci,
                 "base": "chr13_representative", "depth_tag": "full", "gc": 0.35},
    }


# -- independence audit --------------------------------------------------------


def test_audit_passes_on_a_consistent_cell():
    cell = _cell()
    audit = bench.verify_caches({"c": cell})
    assert audit["all_ok"] is True
    assert audit["cells"]["c"]["coordinate_overlaps_historical_regions"] == []


@pytest.mark.parametrize("corrupt", ["contig", "region", "tiling", "nonfinite", "outside"])
def test_audit_fails_on_every_kind_of_corruption(corrupt):
    cell = _cell()
    if corrupt == "contig":
        cell["contig"] = "chr21"
    elif corrupt == "region":
        cell["region"] = [0, 10]
    elif corrupt == "tiling":
        cell["positions"] = cell["positions"][:-64]
    elif corrupt == "nonfinite":
        cell["pb_llr"][3] = np.nan
    elif corrupt == "outside":
        cell["positions"] = cell["positions"] + 10_000_000
    audit = bench.verify_caches({"c": cell})
    assert audit["all_ok"] is False


def test_audit_flags_a_coordinate_overlap_with_a_historical_region():
    """Even on an allowed contig, an overlapping span must fail the audit."""
    cell = _cell(contig="chr13", start=70_000_000)
    cell["spec"] = {**cell["spec"], "contig": "chr13"}
    bench.HISTORICAL_REGIONS.append(("chr13", 70_000_000, 70_001_000))
    try:
        audit = bench.verify_caches({"c": cell})
        assert audit["all_ok"] is False
        assert audit["cells"]["c"]["coordinate_overlaps_historical_regions"]
    finally:
        bench.HISTORICAL_REGIONS.pop()


# -- new stratifier and enrichment --------------------------------------------


def test_expand_lifts_frame_space_masks_back_to_cache_space():
    frame = np.array([True, False, True, True])
    inner = np.array([False, True, False])
    lifted = bench._expand(frame, inner)
    assert lifted.tolist() == [False, False, True, False]


def test_mapq_strata_partition_the_scored_loci():
    cell = _cell()
    strata = bench.mapq_stratified(cell)
    assert set(strata) == {"mq0-20", "mq20-40", "mq40-60", "mq60-71"}
    total = sum(block["loci"] for block in strata.values())
    assert total == int(((cell["labels"] == 0) | (cell["labels"] == LABEL_SNP)).sum())


def test_mapq_strata_route_with_the_frozen_cutoff():
    """A stratum's routed set must equal the frozen rule applied to that stratum."""
    cell = _cell()
    strata = bench.mapq_stratified(cell)
    frame = (cell["labels"] == 0) | (cell["labels"] == LABEL_SNP)
    mapq, depth = cell["mean_mapq"][frame], cell["depth"][frame]
    llr = cell["binomial_llr"][frame]
    for (lo, hi), key in zip([(0, 20), (20, 40), (40, 60), (60, 71)], strata):
        mask = (depth > 0) & (mapq >= lo) & (mapq < hi)
        expected = int(np.count_nonzero(
            np.abs(llr[mask] - FROZEN_BINOMIAL_THRESHOLD) <= FROZEN_ROUTER_CUTOFF))
        assert strata[key]["pb_compute_loci"] == expected


def test_enrichment_attaches_frame_correct_mapq_and_position(monkeypatch):
    from bench_v13_robustness import disagreement_records

    cell = _cell()
    monkeypatch.setattr(bench, "structure_masks",
                        lambda spec, positions: {"in_alldifficult": positions % 2 == 0,
                                                 "not_in_alldifficult": positions % 2 == 1})
    block = bench.enrich_disagreements(cell, disagreement_records(cell))
    frame = (cell["labels"] == 0) | (cell["labels"] == LABEL_SNP)
    mapq, positions = cell["mean_mapq"][frame], cell["positions"][frame]
    assert block["records"], "synthetic cell produced no disagreements to check"
    for record in block["records"]:
        index = record["frame_index"]
        assert record["mean_mapq"] == pytest.approx(float(mapq[index]))
        assert record["position"] == int(positions[index])
        assert ("alldifficult" in record["annotations"]) == (positions[index] % 2 == 0)
    assert sum(v for k, v in block["record_summary"].items() if k.startswith("router_")) \
        == len(block["records"])
