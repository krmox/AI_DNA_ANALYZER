"""Regression tests for the bench_v15 safety-layer benchmark driver.

Mirrors ``test_bench_v14.py``: the frozen constants must stay frozen, the test
cells must stay independent of every previously used contig, the acceptance
criteria and verdict mapping must implement devlog 15 §7.5–§7.6 exactly as
written, and no metric may be computed on a partial extraction.
"""

from __future__ import annotations

import inspect
import json

import numpy as np
import pytest

import bench_v15_safety as bench
import cascade
import safety_layer


# -- frozen constants ---------------------------------------------------------


def test_frozen_constants_are_the_devlog_13_values():
    assert bench.FROZEN_BINOMIAL_THRESHOLD == 7.0
    assert bench.FROZEN_PB_THRESHOLD == 10.5
    assert bench.FROZEN_ROUTER_CUTOFF == 5.411872376933351


def test_driver_never_assigns_to_a_frozen_constant():
    source = inspect.getsource(bench)
    for name in ("FROZEN_BINOMIAL_THRESHOLD", "FROZEN_PB_THRESHOLD", "FROZEN_ROUTER_CUTOFF"):
        assert f"{name} =" not in source, f"{name} must be imported and read only"


def test_driver_routes_through_the_frozen_router_function():
    """Stage 1 must be ``cascade.route_stage1``, not a local reimplementation."""
    assert bench.route_stage1 is cascade.route_stage1


# -- test-set shape and independence ------------------------------------------


def test_twelve_cells_four_regions_three_depths():
    cells = bench.test_cells()
    assert len(cells) == 12
    assert len({spec["base"] for spec in cells.values()}) == 4
    assert len({spec["depth_tag"] for spec in cells.values()}) == 3


def test_no_test_cell_is_on_a_previously_used_chromosome():
    for key, spec in bench.test_cells().items():
        assert spec["contig"] not in bench.PREVIOUSLY_USED_CONTIGS, key


def test_validation_and_test_contigs_are_disjoint():
    validation = {spec["contig"] for spec in bench.validation_cells().values()}
    test = {spec["contig"] for spec in bench.test_cells().values()}
    assert not (validation & test)


def test_test_regions_match_the_selection_record():
    """The coordinates scored must be the ones select_v15_regions.py chose."""
    from pathlib import Path
    record = Path("results/bench_v15/region_selection.json")
    if not record.exists():
        pytest.skip("region_selection.json not present in this checkout")
    selected = json.loads(record.read_text())["selected"]
    for name, spec in bench.TEST_REGIONS.items():
        assert selected[name]["contig"] == spec["contig"]
        assert selected[name]["start"] == spec["start"]
        assert selected[name]["stop"] == spec["stop"]


def test_every_region_is_a_whole_megabase_window():
    for name, spec in bench.TEST_REGIONS.items():
        assert spec["stop"] - spec["start"] == 1_000_000, name
        assert spec["start"] % 1_000_000 == 0, name


# -- the candidate grid -------------------------------------------------------


def test_candidate_grid_covers_every_pre_registered_family():
    families = {rule["family"] for rule in bench.candidate_rules()}
    assert families == {"S", "Sq", "K", "KV", "SK"}


def test_candidate_names_are_unique():
    names = [rule["name"] for rule in bench.candidate_rules()]
    assert len(names) == len(set(names))


def test_frozen_rules_are_in_the_candidate_grid():
    from pathlib import Path
    if not Path("results/bench_v15/rule_selection.json").exists():
        pytest.skip("rule_selection.json not present in this checkout")
    grid = {rule["name"] for rule in bench.candidate_rules()}
    for rule in bench.frozen_rules().values():
        assert rule["name"] in grid


# -- escalation is always disjoint from the frozen router ---------------------


@pytest.fixture
def synthetic_frame():
    rng = np.random.default_rng(15)
    n = rng.integers(0, 90, size=3000).astype(np.float64)
    k = np.floor(n * rng.random(3000) ** 3).astype(np.float64)
    binomial_llr = safety_layer.binomial_llr_at_epsilon(k, n, 0.01)
    return {
        "snp": rng.random(3000) < 0.002,
        "binomial_llr": binomial_llr,
        "pb_llr": binomial_llr + rng.normal(0, 4, size=3000),
        "k": k, "n": n, "vaf": np.where(n > 0, k / np.maximum(n, 1), 0.0),
        "quality_sum": n * rng.uniform(15, 40, size=3000),
        "depth": n,
        "mean_base_quality": rng.uniform(15, 40, size=3000),
        "mean_mapq": np.full(3000, 60.0),
    }


@pytest.mark.parametrize("rule", bench.candidate_rules())
def test_every_candidate_escalation_is_disjoint_from_routing(rule, synthetic_frame):
    routed = cascade.route_stage1(synthetic_frame["binomial_llr"])
    mask = bench.escalation_mask(rule, synthetic_frame, routed)
    assert not np.any(mask & routed), rule["name"]


def test_unknown_family_is_rejected(synthetic_frame):
    routed = cascade.route_stage1(synthetic_frame["binomial_llr"])
    with pytest.raises(ValueError):
        bench.escalation_mask({"family": "nonsense", "k_min": 1}, synthetic_frame, routed)


def test_lost_true_snp_counts_only_pb_positives_the_arm_missed():
    snp = np.array([True, True, True, False])
    pb = np.array([True, True, False, True])
    calls = np.array([True, False, False, True])
    assert bench.lost_true_snp(calls, pb, snp) == 1


# -- acceptance criteria (§7.5) and verdict mapping (§7.6) --------------------


def _pooled(cascade_lost, arm_lost, ci, pb_multiple, random_lost, reverse_lost):
    return {
        "cheap_router_pb": {"lost_true_snp_vs_pb": cascade_lost},
        "arms": {"frozen_full": {"lost_true_snp_vs_pb": arm_lost,
                                 "vs_pb_paired_bootstrap": {"delta_f1_ci95": ci}}},
        "safety": {"frozen_full": {"pb_fraction_multiple_vs_frozen": pb_multiple}},
        "controls": {"frozen_full": {
            "random_escalation": {"min_lost_true_snp": random_lost},
            "reverse_confidence_escalation": {"lost_true_snp_vs_pb": reverse_lost}}},
    }


def _wallclock(relative_speedup, caller_only=None):
    """A3 reads the end-to-end ratio; the caller-only one is reported alongside."""
    return {"arms": {"frozen_full": {
        "end_to_end_speedup_relative_to_frozen_cascade": relative_speedup,
        "caller_only_speedup_relative_to_frozen_cascade":
            relative_speedup if caller_only is None else caller_only}}}


def test_a1_requires_at_least_halving_the_losses():
    good = bench.acceptance(_pooled(6, 3, [-0.0001, 0.001], 1.2, 6, 6), _wallclock(0.99))
    assert good["A1_failure_repaired"]
    weak = bench.acceptance(_pooled(6, 4, [-0.0001, 0.001], 1.2, 6, 6), _wallclock(0.99))
    assert not weak["A1_failure_repaired"]


def test_a1_is_uninformative_when_the_cascade_lost_nothing():
    criteria = bench.acceptance(_pooled(0, 0, [-0.0001, 0.001], 1.2, 0, 0), _wallclock(0.99))
    assert not criteria["A1_informative"]
    assert bench.verdict_from(criteria, 0) == "C. NULL / INCONCLUSIVE"


def test_a2_uses_the_lower_bound_of_the_ci():
    assert bench.acceptance(_pooled(4, 0, [-0.0009, 0.002], 1.2, 4, 4),
                            _wallclock(0.99))["A2_non_inferior_to_pb"]
    assert not bench.acceptance(_pooled(4, 0, [-0.0011, 0.002], 1.2, 4, 4),
                                _wallclock(0.99))["A2_non_inferior_to_pb"]


def test_a3_fails_when_pb_budget_more_than_doubles():
    criteria = bench.acceptance(_pooled(4, 0, [-0.0001, 0.001], 2.4, 4, 4), _wallclock(0.99))
    assert not criteria["A3_compute_preserved"]
    assert bench.verdict_from(criteria, 4) == "E. ARCHITECTURAL FAILURE"


def test_a3_fails_when_the_speedup_collapses():
    criteria = bench.acceptance(_pooled(4, 0, [-0.0001, 0.001], 1.1, 4, 4), _wallclock(0.5))
    assert not criteria["A3_compute_preserved"]


def test_a4_requires_beating_both_controls():
    assert bench.acceptance(_pooled(4, 0, [-0.0001, 0.001], 1.2, 2, 3),
                            _wallclock(0.99))["A4_beats_matched_budget_controls"]
    assert not bench.acceptance(_pooled(4, 1, [-0.0001, 0.001], 1.2, 1, 3),
                                _wallclock(0.99))["A4_beats_matched_budget_controls"]


def test_verdict_strong_positive_needs_all_four_and_enough_events():
    criteria = bench.acceptance(_pooled(6, 0, [-0.0001, 0.001], 1.2, 5, 6), _wallclock(0.99))
    assert bench.verdict_from(criteria, 6) == "A. STRONG POSITIVE"


def test_verdict_weak_positive_when_controls_do_not_separate():
    criteria = bench.acceptance(_pooled(6, 0, [-0.0001, 0.001], 1.2, 0, 0), _wallclock(0.99))
    assert bench.verdict_from(criteria, 6) == "B. WEAK POSITIVE"


def test_verdict_weak_positive_when_too_few_events():
    criteria = bench.acceptance(_pooled(2, 0, [-0.0001, 0.001], 1.2, 2, 2), _wallclock(0.99))
    assert bench.verdict_from(criteria, 2) == "B. WEAK POSITIVE"


def test_verdict_negative_when_the_failure_is_not_repaired():
    criteria = bench.acceptance(_pooled(6, 6, [-0.0001, 0.001], 1.2, 6, 6), _wallclock(0.99))
    assert bench.verdict_from(criteria, 6) == "D. NEGATIVE"


def test_verdict_negative_when_accuracy_drops_below_the_bound():
    criteria = bench.acceptance(_pooled(6, 0, [-0.01, -0.005], 1.2, 6, 6), _wallclock(0.99))
    assert bench.verdict_from(criteria, 6) == "D. NEGATIVE"


def test_architectural_failure_outranks_negative():
    """A repair that destroys the compute advantage is E, not D."""
    criteria = bench.acceptance(_pooled(6, 6, [-0.0001, 0.001], 5.0, 6, 6), _wallclock(0.2))
    assert bench.verdict_from(criteria, 6) == "E. ARCHITECTURAL FAILURE"


# -- refuse to score a partial benchmark --------------------------------------


def test_run_test_refuses_when_an_extraction_is_missing(monkeypatch, tmp_path):
    monkeypatch.setattr(bench, "TEST_CACHE", tmp_path)
    monkeypatch.setattr(bench, "test_cells", lambda: {
        "missing_full": {"contig": "chr16", "start": 0, "stop": 1, "base": "missing",
                         "depth_tag": "full", "path": str(tmp_path / "missing_full.npz")}})
    with pytest.raises(SystemExit):
        bench.run_test(measure=False)


# -- leakage boundary ---------------------------------------------------------


def test_no_escalation_function_can_see_truth():
    forbidden = ("label", "truth", "snp", "vcf", "position")
    for function in (bench.escalation_mask, bench.candidate_rules):
        for name in inspect.signature(function).parameters:
            assert not any(token in name.lower() for token in forbidden), function.__name__


def test_structure_annotation_is_never_a_rule_input():
    """GIAB stratification BEDs may label results; they may not drive routing."""
    source = inspect.getsource(bench.escalation_mask)
    for banned in ("segdup", "difficult", "tandem", "STRATIFICATIONS", "_bed_mask"):
        assert banned not in source


def test_a3_is_evaluated_on_the_end_to_end_accounting():
    """§7.5 A3 says 'end-to-end'; the caller-only figure is reported, not applied."""
    criteria = bench.acceptance(_pooled(4, 0, [-0.0001, 0.001], 1.2, 4, 4),
                                _wallclock(0.95, caller_only=0.45))
    assert criteria["A3_accounting"] == "end_to_end"
    assert criteria["A3_compute_preserved"]
    assert not criteria["A3_would_pass_on_caller_only_accounting"]
    assert criteria["A3_caller_only_speedup_relative_to_frozen_cascade"] == 0.45


def test_a3_falls_back_to_caller_only_when_no_extraction_timing_exists():
    wallclock = {"arms": {"frozen_full": {
        "caller_only_speedup_relative_to_frozen_cascade": 0.9}}}
    criteria = bench.acceptance(_pooled(4, 0, [-0.0001, 0.001], 1.2, 4, 4), wallclock)
    assert criteria["A3_accounting"] == "caller_only"
    assert criteria["A3_compute_preserved"]
