"""Tests for the bench_v14 cross-chromosome benchmark driver.

Three things are worth testing here, and they are the three things that could
silently corrupt the experiment:

1. **The frozen constants are read, never written.** The whole experiment is
   worthless if the driver can retune the router.
2. **The decision rule matches the pre-registration**, including the ordering
   that makes UNDERPOWERED win over PRESERVED.
3. **Coordinates and difficulty masks are correct**, because the structure
   strata and the leakage audit both rest on them, and the caches do not store
   real coordinates.
"""

from __future__ import annotations

import gzip
import inspect
from pathlib import Path

import numpy as np
import pytest

import bench_v14_crosschrom as bench
from cascade import (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD,
                     FROZEN_ROUTER_CUTOFF)


# -- the frozen configuration -------------------------------------------------


def test_frozen_constants_are_the_devlog_13_values():
    assert FROZEN_BINOMIAL_THRESHOLD == 7.0
    assert FROZEN_PB_THRESHOLD == 10.5
    assert FROZEN_ROUTER_CUTOFF == 5.411872376933351


def test_driver_never_assigns_to_a_frozen_constant():
    source = Path("bench_v14_crosschrom.py").read_text()
    for name in ("FROZEN_BINOMIAL_THRESHOLD", "FROZEN_PB_THRESHOLD", "FROZEN_ROUTER_CUTOFF"):
        assert f"{name} =" not in source, f"{name} is assigned in the benchmark driver"


def test_frozen_cutoff_is_in_the_diagnostic_sweep():
    from robustness_benchmark import CUTOFF_SWEEP
    assert FROZEN_ROUTER_CUTOFF in CUTOFF_SWEEP


def test_no_routing_function_can_see_truth():
    """Leakage boundary: routing depends on caller scores, never on labels."""
    for function in (bench.projected_wallclock, bench.accuracy_compute_curve):
        parameters = set(inspect.signature(function).parameters)
        assert "labels" not in parameters and "truth" not in parameters


# -- the pre-registered decision rule -----------------------------------------


def _entry(delta, low, high, snp):
    return {"delta_f1_vs_pb": delta, "snp": snp,
            "vs_pb_paired_bootstrap": {"delta_f1_ci95": [low, high]}}


def test_preserved_requires_both_ci_containing_zero_and_small_delta():
    assert bench.classify(_entry(0.0001, -0.0005, 0.0007, 500)) == "PRESERVED"


def test_degraded_when_ci_excludes_zero_below():
    assert bench.classify(_entry(-0.002, -0.004, -0.001, 500)) == "DEGRADED"


def test_degraded_when_delta_is_large_even_if_ci_contains_zero():
    """|ΔF1| >= 0.001 is DEGRADED by the rule, regardless of the interval."""
    assert bench.classify(_entry(-0.005, -0.009, 0.001, 500)) == "DEGRADED"


def test_improved_when_ci_excludes_zero_above():
    assert bench.classify(_entry(0.003, 0.001, 0.006, 500)) == "IMPROVED"


def test_underpowered_wins_over_preserved_when_snp_count_is_small():
    """A cell with < 100 SNP is never rounded to PRESERVED."""
    assert bench.classify(_entry(0.0, -0.0002, 0.0002, 40)) == "UNDERPOWERED"


def test_underpowered_when_interval_is_wide():
    assert bench.classify(_entry(0.0, -0.05, 0.05, 5000)) == "UNDERPOWERED"


def test_underpowered_when_bootstrap_is_missing():
    assert bench.classify({"delta_f1_vs_pb": 0.0, "snp": 500}) == "UNDERPOWERED"


# -- cell inventory -----------------------------------------------------------


def test_twelve_cells_four_regions_three_depths():
    specs = bench.cells()
    assert len(specs) == 12
    assert {s["base"] for s in specs.values()} == set(bench.BASE_REGIONS)
    assert {s["depth_tag"] for s in specs.values()} == set(bench.DEPTH_TAGS)


def test_no_cell_is_on_a_previously_used_chromosome():
    """The independence claim, asserted rather than asserted-in-prose."""
    for spec in bench.cells().values():
        assert spec["contig"] not in bench.PREVIOUSLY_USED_CONTIGS


def test_every_region_is_window_aligned():
    for spec in bench.cells().values():
        assert (spec["stop"] - spec["start"]) % bench.SEQ_LEN == 0


# -- coordinates and difficulty masks -----------------------------------------


def test_bed_mask_marks_exactly_the_covered_bases(tmp_path):
    bed = tmp_path / "s.bed.gz"
    with gzip.open(bed, "wt") as handle:
        handle.write("#comment\nchr7\t100\t110\nchr7\t200\t205\nchr8\t100\t110\n")
    mask = bench._bed_mask(str(bed), "chr7", 90, 210)
    assert mask.sum() == 15
    assert mask[10:20].all() and mask[110:115].all()
    assert not mask[0:10].any() and not mask[20:110].any()


def test_bed_mask_accepts_either_contig_spelling(tmp_path):
    bed = tmp_path / "s.bed"
    bed.write_text("7\t100\t110\n")
    assert bench._bed_mask(str(bed), "chr7", 90, 120).sum() == 10


def test_bed_mask_clips_intervals_to_the_window(tmp_path):
    bed = tmp_path / "s.bed"
    bed.write_text("chr7\t0\t1000\n")
    mask = bench._bed_mask(str(bed), "chr7", 100, 150)
    assert mask.all() and mask.size == 50


def test_structure_masks_are_complementary_pairs():
    spec = {"contig": "chr7", "start": 0, "stop": 100, "base": "x", "depth_tag": "full"}
    positions = np.arange(0, 100, dtype=np.int64)
    original = bench.STRATIFICATIONS
    try:
        bench.STRATIFICATIONS = {}
        masks = bench.structure_masks(spec, positions)
        assert masks == {}
    finally:
        bench.STRATIFICATIONS = original


def test_structure_masks_partition_the_loci(tmp_path):
    bed = tmp_path / "s.bed"
    bed.write_text("chr7\t10\t20\n")
    spec = {"contig": "chr7", "start": 0, "stop": 100, "base": "x", "depth_tag": "full"}
    positions = np.arange(0, 100, dtype=np.int64)
    original = bench.STRATIFICATIONS
    try:
        bench.STRATIFICATIONS = {"t": str(bed)}
        masks = bench.structure_masks(spec, positions)
        assert masks["in_t"].sum() == 10
        assert np.array_equal(masks["in_t"], ~masks["not_in_t"])
    finally:
        bench.STRATIFICATIONS = original


def test_structure_masks_follow_positions_not_locus_order(tmp_path):
    """Loci are labelled by their coordinate, so a gap in the tiling must shift
    the labels with it -- this is the bug that cache-order indices would cause."""
    bed = tmp_path / "s.bed"
    bed.write_text("chr7\t50\t60\n")
    spec = {"contig": "chr7", "start": 0, "stop": 100, "base": "x", "depth_tag": "full"}
    positions = np.array([0, 1, 55, 56, 90], dtype=np.int64)
    original = bench.STRATIFICATIONS
    try:
        bench.STRATIFICATIONS = {"t": str(bed)}
        masks = bench.structure_masks(spec, positions)
        assert list(masks["in_t"]) == [False, False, True, True, False]
    finally:
        bench.STRATIFICATIONS = original


# -- the audit ----------------------------------------------------------------


def _region(contig="chr20", start=0, stop=128, n=128, positions=None, contig_stored=None):
    labels = np.zeros(n, dtype=np.int64)
    return {
        "contig": contig_stored if contig_stored is not None else contig,
        "region": [start, stop],
        "labels": labels,
        "binomial_llr": np.zeros(n), "pb_llr": np.zeros(n),
        "positions": np.arange(start, start + n) if positions is None else positions,
        "spec": {"contig": contig, "start": start, "stop": stop,
                 "base": "x", "depth_tag": "full"},
    }


def test_audit_passes_a_consistent_cell():
    audit = bench.verify_caches({"a": _region()})
    assert audit["all_ok"] and audit["cells"]["a"]["ok"]


def test_audit_fails_on_contig_mismatch():
    audit = bench.verify_caches({"a": _region(contig_stored="chr9")})
    assert not audit["all_ok"]


def test_audit_fails_on_a_previously_used_contig():
    audit = bench.verify_caches({"a": _region(contig="chr21")})
    assert not audit["all_ok"]
    assert audit["cells"]["a"]["overlaps_previously_used_contig"]


def test_audit_fails_when_locus_count_disagrees_with_the_tiling():
    region = _region()
    region["positions"] = region["positions"][:-1]
    audit = bench.verify_caches({"a": region})
    assert not audit["all_ok"]


def test_audit_fails_on_non_finite_scores():
    region = _region()
    region["pb_llr"] = region["pb_llr"].copy()
    region["pb_llr"][0] = np.nan
    audit = bench.verify_caches({"a": region})
    assert not audit["all_ok"]


def test_audit_fails_when_positions_escape_the_region():
    region = _region()
    region["positions"] = region["positions"] + 10_000
    audit = bench.verify_caches({"a": region})
    assert not audit["all_ok"]


def test_audit_fails_when_loci_are_not_window_aligned():
    audit = bench.verify_caches({"a": _region(n=100)})
    assert not audit["all_ok"]


# -- refusal to score a partial benchmark -------------------------------------


def test_main_refuses_when_an_extraction_is_missing(monkeypatch, tmp_path):
    monkeypatch.setattr(bench, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(bench, "OUT_DIR", tmp_path)
    with pytest.raises(SystemExit) as excinfo:
        bench.main(measure=False)
    assert "missing extractions" in str(excinfo.value)


# -- integration: coordinates vs truth ----------------------------------------


@pytest.mark.skipif(not Path("cache/bench_v14/chr20_neutral_15x.npz").exists(),
                    reason="requires the v14 extraction")
def test_reconstructed_coordinates_land_on_real_vcf_variants():
    """The structure strata and the audit both rest on these coordinates.

    If the tiling were reconstructed wrongly -- an off-by-one, a missed dropped
    window, the wrong chunk width -- SNP labels would land on positions the truth
    VCF says are reference. Requiring exact coincidence in both directions is the
    strongest available check that locus i really is base ``positions[i]``.
    """
    import pysam

    spec = bench.cells()["chr20_neutral_15x"]
    positions = bench.genomic_positions(spec)
    labels = np.load(spec["path"], allow_pickle=True)["labels"]
    assert positions.size == labels.size

    labelled = set(positions[labels == 1].tolist())
    tiled = set(positions.tolist())
    with pysam.VariantFile(str(bench.DATA_DIR / "chr20_neutral.vcf.gz")) as vcf:
        snv = {record.start for record in vcf.fetch(spec["contig"], spec["start"], spec["stop"])
               if len(record.ref) == 1 and record.alts
               and all(len(alt) == 1 for alt in record.alts)}
    assert labelled, "no SNP labels at all -- the truth join is broken"
    assert labelled - snv == set(), "a SNP label sits where the VCF has no SNV"
    assert (snv & tiled) - labelled == set(), "a tiled VCF SNV was not labelled"
