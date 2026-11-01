"""Unit tests for the bench_v13 benchmark infrastructure.

These test the *new* code only. The frozen router/cascade/PB implementations
have their own suites (``test_cascade.py``, ``test_cheap_router.py``,
``test_quality_error.py``) and are not re-tested here.
"""

from __future__ import annotations

import inspect

import numpy as np

import bench_v13_robustness as v13
from cascade import (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD,
                     FROZEN_ROUTER_CUTOFF)


def _region(n: int = 640, seed: int = 0, high_depth: bool = False) -> dict:
    rng = np.random.default_rng(seed)
    counts = np.zeros((n, 10), dtype=np.float32)
    counts[:, 0:4] = rng.integers(0, 30 if high_depth else 12, size=(n, 4))
    counts[:, 6] = counts[:, 0:4].sum(axis=1)
    counts[:, 7] = counts[:, 6] * rng.uniform(20, 38, size=n)
    counts[:, 8] = counts[:, 6] * 60
    counts[:, 9] = rng.integers(1, 5, size=n)
    depth = counts[:, 6].astype(np.float64)
    return {
        "name": "synthetic",
        "counts": counts,
        "binomial_llr": rng.normal(7.0, 6.0, size=n),
        "pb_llr": rng.normal(10.5, 9.0, size=n),
        "labels": rng.integers(0, 2, size=n),
        "depth": depth,
        "mean_base_quality": np.where(depth > 0, counts[:, 7] / np.maximum(depth, 1), 0.0),
        "mean_mapq": np.where(depth > 0, counts[:, 8] / np.maximum(depth, 1), 0.0),
        "region": [0, n],
    }


class TestFrozenConfiguration:
    def test_driver_never_defines_its_own_thresholds(self):
        """Every frozen constant must be imported, never re-declared."""
        source = inspect.getsource(v13)
        code = "".join(part for index, part in enumerate(source.split('"""')) if index % 2 == 0)
        for forbidden in ("FROZEN_ROUTER_CUTOFF =", "FROZEN_PB_THRESHOLD =",
                          "FROZEN_BINOMIAL_THRESHOLD =", "5.411872376933351"):
            assert forbidden not in code, forbidden

    def test_frozen_cutoff_is_the_devlog_11_value(self):
        assert FROZEN_ROUTER_CUTOFF == 5.411872376933351
        assert (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD) == (7.0, 10.5)


class TestArms:
    def test_router_equals_pb_when_everything_is_routed(self):
        region = _region()
        region["binomial_llr"] = np.full(region["labels"].size, FROZEN_BINOMIAL_THRESHOLD)
        entry = v13.pooled_block(region["labels"] == 1, region["binomial_llr"], region["pb_llr"])
        assert entry["disagreement_vs_pb"] == 0
        assert entry["cheap_router_pb"] == entry["pb_only"]
        assert entry["fraction_routed_to_pb"] == 1.0

    def test_router_equals_binomial_when_nothing_is_routed(self):
        region = _region()
        far = np.full(region["labels"].size, FROZEN_BINOMIAL_THRESHOLD + 10 * FROZEN_ROUTER_CUTOFF)
        entry = v13.pooled_block(region["labels"] == 1, far, region["pb_llr"])
        assert entry["fraction_routed_to_pb"] == 0.0
        assert entry["cheap_router_pb"] == entry["binomial_only"]

    def test_delta_f1_is_router_minus_pb(self):
        region = _region(seed=3)
        entry = v13.pooled_block(region["labels"] == 1, region["binomial_llr"], region["pb_llr"])
        assert entry["delta_f1_vs_pb"] == (entry["cheap_router_pb"]["f1"]
                                           - entry["pb_only"]["f1"])


class TestStratification:
    def test_fine_depth_bins_partition_the_scoring_frame(self):
        region = _region(high_depth=True)
        strata = v13.fine_depth_stratified(region)
        counted = sum(entry["loci"] for entry in strata.values())
        frame = v13.scoring_frame(region)
        zero_depth = int(np.count_nonzero(region["depth"][frame] == 0))
        assert counted + zero_depth == int(frame.sum())

    def test_fine_depth_bins_resolve_the_48x_boundary(self):
        """The tensor-width boundary must not be hidden inside a '30x+' bin."""
        edges = {lo for lo, _ in v13.FINE_DEPTH_BINS} | {hi for _, hi in v13.FINE_DEPTH_BINS}
        assert 48 in edges and 47 in edges

    def test_vaf_strata_are_label_free(self):
        region = _region(seed=5)
        flipped = dict(region, labels=1 - region["labels"])
        a, b = v13.vaf_stratified(region), v13.vaf_stratified(flipped)
        assert [entry["loci"] for entry in a.values()] == [entry["loci"] for entry in b.values()]


class TestPbFix:
    def test_scope_uses_k_versus_retained_reads_when_available(self):
        region = _region()
        n = region["labels"].size
        region["k_full"] = np.full(n, 54, dtype=np.int32)
        region["n_counted"] = np.full(n, 48, dtype=np.int32)
        region["k_retained"] = np.full(n, 30, dtype=np.int32)
        scope = v13.pb_fix_scope(region)
        assert scope["loci_k_out_of_support"] == n
        assert scope["loci_k_rescaled"] == n

    def test_scope_falls_back_to_depth_when_k_is_absent(self):
        scope = v13.pb_fix_scope(_region())
        assert "loci_k_out_of_support" not in scope
        assert scope["loci_depth_gt_48"] == 0

    def test_ab_returns_none_without_a_legacy_arm(self):
        assert v13.pb_fix_ab(_region()) is None

    def test_ab_recovers_snp_that_the_legacy_path_zeroed(self):
        """The defect's signature: true SNP scored exactly 0 by the legacy arm."""
        region = _region(seed=11)
        n = region["labels"].size
        region["labels"] = np.zeros(n, dtype=np.int64)
        region["labels"][:50] = 1
        region["pb_llr"] = np.where(region["labels"] == 1, 300.0, -5.0)
        region["pb_llr_legacy"] = np.where(region["labels"] == 1, 0.0, -5.0)
        ab = v13.pb_fix_ab(region)
        assert ab["legacy"]["pb_only"]["tp"] == 0
        assert ab["fixed"]["pb_only"]["tp"] == 50
        assert ab["delta_f1_fixed_minus_legacy"] > 0
        assert ab["legacy"]["pb_llr_zero_snp"] == 50
        assert ab["fixed"]["pb_llr_zero_snp"] == 0


class TestDisagreementAnalysis:
    def test_counts_match_a_manual_computation(self):
        region = _region(seed=7)
        result = v13.disagreement_records(region)
        frame = v13.scoring_frame(region)
        binomial_llr = region["binomial_llr"][frame]
        pb_calls = region["pb_llr"][frame] >= FROZEN_PB_THRESHOLD
        routed = np.abs(binomial_llr - FROZEN_BINOMIAL_THRESHOLD) <= FROZEN_ROUTER_CUTOFF
        router_calls = np.where(routed, pb_calls, binomial_llr >= FROZEN_BINOMIAL_THRESHOLD)
        assert result["disagreements"] == int(np.count_nonzero(router_calls != pb_calls))

    def test_no_disagreement_can_be_inside_the_routed_band(self):
        """A routed locus takes PB's answer by construction, so it cannot differ."""
        region = _region(seed=8)
        for record in v13.disagreement_records(region)["records"]:
            assert record["routed_to_pb"] is False

    def test_effect_labels_are_exhaustive_and_directional(self):
        region = _region(seed=9)
        effects = {record["effect"] for record in v13.disagreement_records(region)["records"]}
        assert effects <= {"router_fn", "router_fp", "router_recovers_fn",
                           "router_avoids_fp", "other"}
        assert "other" not in effects

    def test_stratified_disagreement_counts_sum_to_the_total(self):
        region = _region(seed=10)
        result = v13.disagreement_records(region)
        by_truth = result["by_truth"]["snp"]["disagreements"] \
            + result["by_truth"]["non_snp"]["disagreements"]
        assert by_truth == result["disagreements"]

    def test_record_dump_is_capped_but_counts_are_not(self):
        region = _region(n=6400, seed=12)
        result = v13.disagreement_records(region, limit=3)
        assert len(result["records"]) <= 3
        assert result["disagreements"] >= len(result["records"])
        assert result["records_truncated"] == (result["disagreements"] > 3)


class TestCompute:
    def test_projection_uses_measured_rates_and_routed_fraction(self):
        region = _region(seed=13)
        throughput = {"binomial_loci_per_second": 1e6, "pb_loci_per_second": 1e4}
        projection = v13.projected_wallclock(region, throughput)
        expected = projection["loci"] / 1e6 + projection["pb_loci"] / 1e4
        assert projection["cascade_seconds"] == expected
        assert projection["speedup"] == projection["pb_only_seconds"] / expected

    def test_pb_compute_fraction_equals_routed_fraction(self):
        region = _region(seed=14)
        throughput = {"binomial_loci_per_second": 1e6, "pb_loci_per_second": 1e4}
        projection = v13.projected_wallclock(region, throughput)
        frame = v13.scoring_frame(region)
        routed = np.abs(region["binomial_llr"][frame] - FROZEN_BINOMIAL_THRESHOLD) \
            <= FROZEN_ROUTER_CUTOFF
        assert projection["pb_compute_fraction"] == float(routed.mean())


class TestLeakage:
    def test_no_analysis_function_routes_on_truth(self):
        """Routing is a function of the binomial LLR alone, everywhere in this file."""
        source = inspect.getsource(v13)
        code = "".join(part for index, part in enumerate(source.split('"""')) if index % 2 == 0)
        assert "route_mask(binomial_llr" in code or "route_mask(region[\"binomial_llr\"]" in code
        assert "route_mask(snp" not in code and "route_mask(labels" not in code


class TestPooledFailureBoundary:
    def test_pools_counts_across_regions(self):
        region_a = v13.disagreement_records(_region(seed=21))
        region_b = v13.disagreement_records(_region(seed=22))
        pooled = v13.pooled_failure_boundary({"a": region_a, "b": region_b})
        assert pooled["total"]["disagreements"] == (region_a["disagreements"]
                                                    + region_b["disagreements"])
        assert pooled["total"]["loci"] == region_a["loci_scored"] + region_b["loci_scored"]

    def test_empty_bins_report_no_rate_rather_than_zero(self):
        region = v13.disagreement_records(_region(seed=23))
        pooled = v13.pooled_failure_boundary({"a": region})
        for cell in pooled["by_depth"].values():
            assert (cell["rate"] is None) == (cell["loci"] == 0)


class TestAccuracyComputeCurve:
    def test_curve_contains_the_frozen_operating_point_exactly_once(self):
        rows = v13.accuracy_compute_curve(_region(seed=24), None)
        assert sum(row["is_frozen_operating_point"] for row in rows) == 1

    def test_curve_speedup_falls_as_more_loci_are_routed(self):
        throughput = {"binomial_loci_per_second": 1e6, "pb_loci_per_second": 1e4}
        rows = v13.accuracy_compute_curve(_region(seed=25), throughput)
        by_coverage = sorted(rows, key=lambda row: row["pb_compute_fraction"])
        speedups = [row["projected_speedup"] for row in by_coverage]
        assert all(later <= earlier for earlier, later in zip(speedups, speedups[1:]))

    def test_probe_matches_the_regions_own_depth_regime(self):
        assert v13.probe_for("re31_32M_full") == "full_depth"
        assert v13.probe_for("new32_33M_30x") == "30x"
        assert v13.probe_for("gen14_far_13_17M") == "15x"
