"""bench_v15 — can the frozen cascade's one failure mode be detected cheaply?

Protocol is pre-registered in ``History/15_DEVLOG.md`` §1–§8 and is not
re-derived here. Nothing in this file fits, selects or modifies a *router*
threshold: the three frozen cascade constants are imported from ``cascade`` and
only ever read. The one thing this file does select is the **safety rule**, and
it does so exactly once, on the validation cells only, by the §7.4 criterion —
``--select`` writes the frozen rule to disk and ``--test`` consumes it.

Two disjoint sets (devlog 15 §5):

* **validation** — the 12 devlog 14 cells (chr20/chr19/chr4/chr1 × 3 depths).
  Discovery data; contains the 5 known failures; used only to choose the rule.
* **test** — 12 new cells (chr16/chr15/chr7 segdup + chr14 neutral × 3 depths)
  on contigs no previous devlog has touched. Untouched until the rule is frozen.

Every metric helper is imported from the existing benchmark modules rather than
reimplemented, so the numbers sit on the same definitions as devlogs 11–14. What
is genuinely new here is only the safety arm itself, its controls at matched
*additional* PB coverage, and the ``lost_true_snp`` endpoint that names the
failure mode being repaired.
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np

from bench_v12_stage1 import arms_on_mask, strata_masks
from bench_v13_robustness import (fine_depth_stratified, scoring_frame,
                                  vaf_stratified)
from bench_v14_crosschrom import (BASE_REGIONS as V14_BASE_REGIONS, CHUNK_BP,
                                  DEPTH_TAGS, SEQ_LEN, STRATIFICATIONS,
                                  _bed_mask, classify)
from binomial_baseline import BinomialVariantCaller
from cascade import (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD,
                     FROZEN_ROUTER_CUTOFF, route_stage1)
from config import LABEL_SNP
from evaluate_binomial_baseline import prf
from evaluate_quality_error import paired_bootstrap
from robustness_benchmark import accuracy_block, depth_stratified, quality_stratified
from safety_layer import (CANDIDATE_EPSILON_BRACKETS, CANDIDATE_K_MIN,
                          FROZEN_SAFETY_RULE, binomial_llr_at_epsilon,
                          epsilon_flip, safety_escalate)  # noqa: F401

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

#: The four pre-registered test regions (devlog 15 §5.3), each at three depths.
TEST_REGIONS = {
    "chr16_segdup": {"contig": "chr16", "start": 29_000_000, "stop": 30_000_000,
                     "regime": "segmental duplication / low mappability", "gc": 0.490},
    "chr15_segdup": {"contig": "chr15", "start": 30_000_000, "stop": 31_000_000,
                     "regime": "segmental duplication / low mappability", "gc": 0.427},
    "chr7_segdup": {"contig": "chr7", "start": 57_000_000, "stop": 58_000_000,
                    "regime": "segmental duplication / low mappability", "gc": 0.419},
    "chr14_neutral": {"contig": "chr14", "start": 53_000_000, "stop": 54_000_000,
                      "regime": "ordinary euchromatin (control)", "gc": 0.393},
}

#: Contigs used by devlogs 3–14. No test cell may sit on one.
PREVIOUSLY_USED_CONTIGS = {"chr21", "chr20", "chr19", "chr4", "chr1"}

VALIDATION_CACHE = Path("cache/bench_v14")
VALIDATION_DATA = Path("data/giab_hg002_v14")
TEST_CACHE = Path("cache/bench_v15")
TEST_DATA = Path("data/giab_hg002_v15")
OUT_DIR = Path("results/bench_v15")
RULE_FILE = OUT_DIR / "rule_selection.json"
TEST_RESULTS = OUT_DIR / "test_results.json"
AUDIT = OUT_DIR / "independence_audit.json"


def validation_cells() -> dict[str, dict]:
    """The 12 devlog 14 cells, addressed in their existing cache."""
    return {f"{name}_{tag}": {**base, "base": name, "depth_tag": tag, "split": "validation",
                              "cache": VALIDATION_CACHE, "data": VALIDATION_DATA,
                              "path": str(VALIDATION_CACHE / f"{name}_{tag}.npz")}
            for name, base in V14_BASE_REGIONS.items() for tag in DEPTH_TAGS}


def test_cells() -> dict[str, dict]:
    """The 12 pre-registered test cells (devlog 15 §5.3)."""
    return {f"{name}_{tag}": {**base, "base": name, "depth_tag": tag, "split": "test",
                              "cache": TEST_CACHE, "data": TEST_DATA,
                              "path": str(TEST_CACHE / f"{name}_{tag}.npz")}
            for name, base in TEST_REGIONS.items() for tag in DEPTH_TAGS}


# -- loading ------------------------------------------------------------------


def load_cell(key: str, spec: dict) -> dict:
    """Load one extracted cell, adding the ``(k, n)`` the safety layer needs.

    ``k`` and ``n`` are recomputed with the *same* ``BinomialVariantCaller`` the
    cache's ``binomial_llr`` came from, and the recomputed LLR is asserted equal
    to the cached one — so the safety layer is provably reading the caller's own
    candidate-allele choice rather than a re-derivation that might differ.
    """
    blob = np.load(spec["path"], allow_pickle=True)
    counts = blob["counts"]
    depth = blob["depth"].astype(np.float64)
    caller = BinomialVariantCaller(error_rate=0.01)
    score = caller.score_counts(counts[:, 0:4], counts[:, 9].astype(int))
    cached_llr = blob["binomial_llr"].astype(np.float64)
    assert np.allclose(score.llr, cached_llr, atol=1e-8), \
        f"{key}: recomputed binomial LLR disagrees with the cache"
    return {
        "name": key,
        "counts": counts,
        "binomial_llr": cached_llr,
        "pb_llr": blob["pb_llr"].astype(np.float64),
        "labels": blob["labels"],
        "depth": depth,
        "k": score.k,
        "n": score.n,
        "vaf": score.vaf,
        "mean_base_quality": np.where(depth > 0, counts[:, 7] / np.maximum(depth, 1), 0.0),
        "mean_mapq": np.where(depth > 0, counts[:, 8] / np.maximum(depth, 1), 0.0),
        "region": [int(x) for x in blob["region"]],
        "contig": str(blob["contig"][0]) if "contig" in blob else None,
        "counts_time_s": float(blob["timing"][0]),
        "reads_time_s": float(blob["timing"][1]),
        "pb_time_s": float(blob["timing"][2]),
        "bam": str(blob["bam"][0]),
        "spec": spec,
    }


def frame_of(region: dict) -> dict[str, np.ndarray]:
    """Every array the arms need, restricted to the scoring frame, in one dict."""
    f = scoring_frame(region)
    return {
        "snp": region["labels"][f] == LABEL_SNP,
        "binomial_llr": region["binomial_llr"][f],
        "pb_llr": region["pb_llr"][f],
        "k": region["k"][f], "n": region["n"][f], "vaf": region["vaf"][f],
        "quality_sum": region["counts"][:, 7].astype(np.float64)[f],
        "depth": region["depth"][f],
        "mean_base_quality": region["mean_base_quality"][f],
        "mean_mapq": region["mean_mapq"][f],
        "frame": f,
    }


# -- the arms -----------------------------------------------------------------


def arm_masks(d: dict, escalated: np.ndarray | None = None) -> dict[str, np.ndarray]:
    """Boolean SNP calls for every arm, plus the routing masks behind them."""
    binomial_calls = d["binomial_llr"] >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = d["pb_llr"] >= FROZEN_PB_THRESHOLD
    routed = route_stage1(d["binomial_llr"])
    out = {
        "binomial_calls": binomial_calls,
        "pb_calls": pb_calls,
        "routed": routed,
        "cascade_calls": np.where(routed, pb_calls, binomial_calls),
    }
    if escalated is not None:
        to_pb = routed | escalated
        out["escalated"] = escalated
        out["to_pb"] = to_pb
        out["safety_calls"] = np.where(to_pb, pb_calls, binomial_calls)
    return out


def lost_true_snp(calls: np.ndarray, pb_calls: np.ndarray, snp: np.ndarray) -> int:
    """True SNP that PB-only finds and this arm does not — the failure mode."""
    return int(np.count_nonzero(snp & pb_calls & ~calls))


def arm_block(calls: np.ndarray, d: dict, pb_calls: np.ndarray, pb_f1: float,
              bootstrap: bool = True) -> dict:
    """Full metric block for one arm against truth, with PB as the reference."""
    entry = accuracy_block(calls, d["snp"])
    entry["delta_f1_vs_pb"] = entry["f1"] - pb_f1
    entry["lost_true_snp_vs_pb"] = lost_true_snp(calls, pb_calls, d["snp"])
    entry["gained_true_snp_vs_pb"] = int(np.count_nonzero(d["snp"] & ~pb_calls & calls))
    entry["disagreement_vs_pb"] = int(np.count_nonzero(calls != pb_calls))
    entry["disagreement_rate_vs_pb"] = float(np.mean(calls != pb_calls))
    if bootstrap:
        entry["vs_pb_paired_bootstrap"] = paired_bootstrap(calls, pb_calls, d["snp"])
    return entry


# -- candidate rule grid (devlog 15 §7.2) -------------------------------------


def candidate_rules() -> list[dict]:
    """Every candidate in the pre-registered grid, as evaluable specs."""
    rules: list[dict] = []
    for lo, hi in CANDIDATE_EPSILON_BRACKETS:
        rules.append({"family": "S", "name": f"S[{lo:g},{hi:g}]",
                      "epsilon_lo": lo, "epsilon_hi": hi, "k_min": 1})
    for floor in (1e-4, 1e-3, 3e-3):
        rules.append({"family": "Sq", "name": f"Sq(floor={floor:g})",
                      "quality_floor": floor, "k_min": 1})
    for t in (2, 3, 4, 5, 6):
        rules.append({"family": "K", "name": f"K(k>={t})", "k_min": t})
    for t in (3, 4):
        for lo, hi in ((0.05, 0.35), (0.05, 0.25), (0.08, 0.30)):
            rules.append({"family": "KV", "name": f"KV(k>={t},vaf[{lo},{hi}])",
                          "k_min": t, "vaf_lo": lo, "vaf_hi": hi})
    for lo, hi in CANDIDATE_EPSILON_BRACKETS:
        for t in CANDIDATE_K_MIN:
            if t == 1:
                continue  # identical to family S
            rules.append({"family": "SK", "name": f"SK[{lo:g},{hi:g}]&k>={t}",
                          "epsilon_lo": lo, "epsilon_hi": hi, "k_min": t})
    return rules


def escalation_mask(rule: dict, d: dict, routed: np.ndarray) -> np.ndarray:
    """Escalation mask for any candidate family, always disjoint from ``routed``."""
    k, n = d["k"], d["n"]
    if rule["family"] in ("S", "SK"):
        mask = safety_escalate(k, n, routed, rule)
        return mask
    if rule["family"] == "Sq":
        return safety_escalate(k, n, routed, rule, d["quality_sum"])
    if rule["family"] == "K":
        return (k >= rule["k_min"]) & ~routed
    if rule["family"] == "KV":
        return ((k >= rule["k_min"]) & (d["vaf"] >= rule["vaf_lo"])
                & (d["vaf"] <= rule["vaf_hi"]) & ~routed)
    raise ValueError(f"unknown rule family {rule['family']!r}")


# -- rule selection (devlog 15 §7.4) ------------------------------------------

#: Pre-registered feasibility bound on the escalated fraction (devlog 15 §7.4c).
MAX_ESCALATED_FRACTION = 5e-4


def pooled_validation() -> dict:
    """Every validation cell concatenated onto one scoring frame."""
    parts: dict[str, list] = {}
    cells = validation_cells()
    for key, spec in cells.items():
        region = load_cell(key, spec)
        d = frame_of(region)
        d["cell"] = np.full(int(d["frame"].sum()), key, dtype=object)
        for field, value in d.items():
            if field == "frame":
                continue
            parts.setdefault(field, []).append(value)
    return {field: np.concatenate(values) for field, values in parts.items()}


def evaluate_candidate(rule: dict, d: dict, masks: dict, pb_f1: float) -> dict:
    """One candidate's budget and accuracy on the pooled validation frame."""
    escalated = escalation_mask(rule, d, masks["routed"])
    calls = np.where(masks["routed"] | escalated, masks["pb_calls"], masks["binomial_calls"])
    stats = prf(calls, d["snp"])
    n = d["snp"].size
    return {
        **{key: value for key, value in rule.items()},
        "escalated": int(escalated.sum()),
        "escalated_fraction": float(escalated.mean()),
        "total_pb_fraction": float((masks["routed"] | escalated).mean()),
        "pb_fraction_multiple_vs_frozen": float((masks["routed"] | escalated).sum()
                                                / max(masks["routed"].sum(), 1)),
        "f1": stats["f1"], "precision": stats["precision"], "recall": stats["recall"],
        "tp": stats["tp"], "fp": stats["fp"], "fn": stats["fn"],
        "delta_f1_vs_pb": stats["f1"] - pb_f1,
        "lost_true_snp_vs_pb": lost_true_snp(calls, masks["pb_calls"], d["snp"]),
        "loci": n,
    }


def _tie_key(row: dict) -> tuple:
    """§7.4 tie-break: narrower bracket, then larger k_min, then grid order."""
    lo = row.get("epsilon_lo")
    hi = row.get("epsilon_hi")
    width = (hi / lo) if (lo and hi) else float("inf")
    return (row["escalated_fraction"], width, -row.get("k_min", 1))


def select_rule() -> dict:
    """Run the §7.4 selection on validation data and write the frozen rule."""
    logger.info("loading validation cells")
    d = pooled_validation()
    masks = arm_masks(d)
    pb_f1 = prf(masks["pb_calls"], d["snp"])["f1"]
    cascade_stats = prf(masks["cascade_calls"], d["snp"])
    frozen_lost = lost_true_snp(masks["cascade_calls"], masks["pb_calls"], d["snp"])
    logger.info("validation: %d loci, %d SNP, routed %.4e, frozen-cascade lost %d true SNP",
                d["snp"].size, int(d["snp"].sum()), float(masks["routed"].mean()), frozen_lost)

    rows = [evaluate_candidate(rule, d, masks, pb_f1) for rule in candidate_rules()]

    def feasible(max_lost: int) -> list[dict]:
        return [r for r in rows
                if r["lost_true_snp_vs_pb"] <= max_lost
                and r["delta_f1_vs_pb"] >= 0.0
                and r["escalated_fraction"] <= MAX_ESCALATED_FRACTION]

    full_set, lean_set = feasible(0), feasible(1)
    report = {
        "protocol": "History/15_DEVLOG.md §7.4",
        "validation_cells": sorted(validation_cells()),
        "validation_loci": int(d["snp"].size),
        "validation_snp": int(d["snp"].sum()),
        "frozen_router_routed_fraction": float(masks["routed"].mean()),
        "pb_only_f1": pb_f1,
        "frozen_cascade_f1": cascade_stats["f1"],
        "frozen_cascade_delta_f1_vs_pb": cascade_stats["f1"] - pb_f1,
        "frozen_cascade_lost_true_snp": frozen_lost,
        "max_escalated_fraction": MAX_ESCALATED_FRACTION,
        "candidates": rows,
        "n_feasible_full": len(full_set),
        "n_feasible_lean": len(lean_set),
    }
    report["frozen_full"] = min(full_set, key=_tie_key) if full_set else None
    report["frozen_lean"] = min(lean_set, key=_tie_key) if lean_set else None
    if report["frozen_lean"] is None:
        report["outcome"] = "NO_RULE_FROZEN"
    else:
        report["outcome"] = "RULE_FROZEN"
    report["frozen_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RULE_FILE.write_text(json.dumps(report, indent=2))
    return report


def frozen_rules() -> dict[str, dict]:
    """The frozen rule(s), read back from disk. Never recomputed at test time."""
    report = json.loads(RULE_FILE.read_text())
    out = {}
    for tag in ("frozen_full", "frozen_lean"):
        rule = report.get(tag)
        if rule is None:
            continue
        out[tag] = {key: rule[key] for key in
                    ("family", "name", "k_min", "epsilon_lo", "epsilon_hi",
                     "quality_floor", "vaf_lo", "vaf_hi") if key in rule}
    if "frozen_full" in out and out.get("frozen_lean") == out.get("frozen_full"):
        out.pop("frozen_lean")
    return out


# -- coordinates and structure labels -----------------------------------------


def genomic_positions(spec: dict) -> np.ndarray:
    """True GRCh38 coordinate of every locus in a cell, in cache order.

    Same reconstruction as ``bench_v14_crosschrom.genomic_positions`` (the
    caches store chunk-relative indices, not coordinates — devlog 12 §16.5),
    generalised over the data directory so validation and test cells both work.
    """
    from pileup_counts import PileupCountsProvider

    data = Path(spec["data"])
    positions = []
    for chunk_start in range(spec["start"], spec["stop"], CHUNK_BP):
        chunk_stop = min(chunk_start + CHUNK_BP, spec["stop"])
        provider = PileupCountsProvider(
            fasta_path=f"data/reference/{spec['contig']}_full.fa",
            bam_path=str(data / f"{spec['base']}_{spec['depth_tag']}.bam"),
            vcf_path=str(data / f"{spec['base']}.vcf.gz"),
            contig=spec["contig"], region=(chunk_start, chunk_stop), seq_len=SEQ_LEN,
            high_confidence_bed=str(data / f"{spec['base']}_highconf.bed"),
        )
        with provider:
            for window in provider._build_windows():
                positions.append(np.arange(window.start, window.end, dtype=np.int64))
    return np.concatenate(positions) if positions else np.zeros(0, dtype=np.int64)


def structure_masks(spec: dict, positions: np.ndarray) -> dict[str, np.ndarray]:
    """Per-locus external difficulty labels, and their complements.

    Defined entirely by public annotation and the locus' coordinate. No caller
    score, no safety decision and no truth label participates. These are used
    for *reporting* only: devlog 15 §7.1 bans annotation from every rule input.
    """
    out = {}
    offset = positions - spec["start"]
    for key, path in STRATIFICATIONS.items():
        inside = _bed_mask(path, spec["contig"], spec["start"], spec["stop"])[offset]
        out[f"in_{key}"] = inside
        out[f"not_in_{key}"] = ~inside
    return out


# -- independence audit (devlog 15 §6) ----------------------------------------


def verify_caches(loaded: dict[str, dict]) -> dict:
    """Mechanical validation and independence audit, run before any test metric."""
    audit = {
        "previously_used_contigs": sorted(PREVIOUSLY_USED_CONTIGS),
        "rule_frozen_at": json.loads(RULE_FILE.read_text())["frozen_at"],
        "rule_file_mtime": time.strftime("%Y-%m-%dT%H:%M:%S%z",
                                         time.localtime(RULE_FILE.stat().st_mtime)),
        "cells": {},
    }
    for key, region in loaded.items():
        spec = region["spec"]
        positions = region["positions"]
        labels = region["labels"]
        cache_mtime = Path(spec["path"]).stat().st_mtime
        entry = {
            "contig": region["contig"],
            "region": region["region"],
            "expected_region": [spec["start"], spec["stop"]],
            "loci_cached": int(labels.size),
            "loci_from_window_tiling": int(positions.size),
            "window_aligned": bool(labels.size % SEQ_LEN == 0),
            "non_finite_binomial": int((~np.isfinite(region["binomial_llr"])).sum()),
            "non_finite_pb": int((~np.isfinite(region["pb_llr"])).sum()),
            "position_min": int(positions.min()) if positions.size else None,
            "position_max": int(positions.max()) if positions.size else None,
            "positions_inside_region": bool(
                positions.size and positions.min() >= spec["start"]
                and positions.max() < spec["stop"]),
            "overlaps_previously_used_contig": region["contig"] in PREVIOUSLY_USED_CONTIGS,
            "cache_mtime": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(cache_mtime)),
            "extracted_after_rule_was_frozen": bool(cache_mtime > RULE_FILE.stat().st_mtime),
            "bam": region["bam"],
        }
        entry["ok"] = bool(
            entry["contig"] == spec["contig"]
            and entry["region"] == entry["expected_region"]
            and entry["window_aligned"]
            and entry["loci_cached"] == entry["loci_from_window_tiling"]
            and entry["non_finite_binomial"] == 0 and entry["non_finite_pb"] == 0
            and entry["positions_inside_region"]
            and not entry["overlaps_previously_used_contig"]
            and entry["extracted_after_rule_was_frozen"])
        audit["cells"][key] = entry
    audit["n_cells"] = len(audit["cells"])
    audit["all_ok"] = all(cell["ok"] for cell in audit["cells"].values())
    return audit


# -- scoring one cell ---------------------------------------------------------


def escalations_for(d: dict, rules: dict[str, dict],
                    routed: np.ndarray) -> dict[str, np.ndarray]:
    """Escalation masks for every frozen rule, on one frame."""
    return {tag: escalation_mask(rule, d, routed) for tag, rule in rules.items()}


def safety_diagnostics(d: dict, masks: dict, escalated: np.ndarray) -> dict:
    """What the escalations actually did, as a routing mechanism."""
    cascade_calls, pb_calls = masks["cascade_calls"], masks["pb_calls"]
    safety_calls = np.where(masks["routed"] | escalated, pb_calls, masks["binomial_calls"])
    changed = safety_calls != cascade_calls
    cascade_disagrees = cascade_calls != pb_calls
    unrouted_disagreements = cascade_disagrees & ~masks["routed"]
    cascade_lost = d["snp"] & pb_calls & ~cascade_calls
    n_escalated = int(escalated.sum())
    return {
        "escalated": n_escalated,
        "escalated_fraction": float(escalated.mean()),
        "escalated_fraction_of_unrouted": float(
            n_escalated / max(int((~masks["routed"]).sum()), 1)),
        "total_pb_fraction": float((masks["routed"] | escalated).mean()),
        "pb_fraction_multiple_vs_frozen": float(
            (masks["routed"] | escalated).sum() / max(masks["routed"].sum(), 1)),
        "calls_changed_vs_frozen_cascade": int(changed.sum()),
        "capture_of_unrouted_disagreements": int((escalated & unrouted_disagreements).sum()),
        "unrouted_disagreements": int(unrouted_disagreements.sum()),
        "capture_rate": float((escalated & unrouted_disagreements).sum()
                              / max(int(unrouted_disagreements.sum()), 1)),
        "capture_of_cascade_lost_true_snp": int((escalated & cascade_lost).sum()),
        "cascade_lost_true_snp": int(cascade_lost.sum()),
        "false_escalations": int((escalated & ~cascade_disagrees).sum()),
        "false_escalation_rate": float((escalated & ~cascade_disagrees).sum()
                                       / max(n_escalated, 1)),
        "escalations_that_fixed_a_call": int((changed & (safety_calls == d["snp"])
                                              & (cascade_calls != d["snp"])).sum()),
        "escalations_that_broke_a_call": int((changed & (cascade_calls == d["snp"])
                                              & (safety_calls != d["snp"])).sum()),
    }


def controls_at_matched_budget(d: dict, masks: dict, escalated: np.ndarray,
                               seeds=(0, 1, 2, 3, 4)) -> dict:
    """Random and reverse-confidence escalation at the *same* extra PB budget.

    The frozen router's routed set is held fixed in every control; only the
    *additional* loci differ. That isolates the question the controls exist to
    answer: is the benefit the rule's structure, or just the extra PB compute?
    """
    routed, pb_calls = masks["routed"], masks["pb_calls"]
    binomial_calls = masks["binomial_calls"]
    unrouted = np.flatnonzero(~routed)
    budget = int(escalated.sum())
    pb_f1 = prf(pb_calls, d["snp"])["f1"]

    def score(mask: np.ndarray) -> dict:
        calls = np.where(routed | mask, pb_calls, binomial_calls)
        stats = prf(calls, d["snp"])
        return {"f1": stats["f1"], "delta_f1_vs_pb": stats["f1"] - pb_f1,
                "tp": stats["tp"], "fp": stats["fp"], "fn": stats["fn"],
                "lost_true_snp_vs_pb": lost_true_snp(calls, pb_calls, d["snp"]),
                "extra_pb_loci": int(mask.sum())}

    random_runs = []
    for seed in seeds:
        mask = np.zeros(d["snp"].size, dtype=bool)
        if budget and unrouted.size:
            rng = np.random.default_rng(seed)
            mask[rng.choice(unrouted, size=min(budget, unrouted.size), replace=False)] = True
        random_runs.append(score(mask))

    # reverse-confidence: the unrouted loci FARTHEST from the binomial threshold
    margin = np.abs(d["binomial_llr"] - FROZEN_BINOMIAL_THRESHOLD)
    reverse = np.zeros(d["snp"].size, dtype=bool)
    if budget and unrouted.size:
        order = unrouted[np.argsort(margin[unrouted])[::-1][:budget]]
        reverse[order] = True

    return {
        "matched_extra_pb_coverage": budget,
        "safety_layer": score(escalated),
        "random_escalation": {
            "runs": random_runs,
            "mean_f1": float(np.mean([r["f1"] for r in random_runs])),
            "std_f1": float(np.std([r["f1"] for r in random_runs])),
            "mean_lost_true_snp": float(np.mean([r["lost_true_snp_vs_pb"] for r in random_runs])),
            "min_lost_true_snp": int(min(r["lost_true_snp_vs_pb"] for r in random_runs)),
        },
        "reverse_confidence_escalation": score(reverse),
    }


def accuracy_compute_curve(d: dict, masks: dict) -> list:
    """Every candidate in the §7.2 grid, as an (extra PB, accuracy) point.

    **Diagnostic only.** The frozen rule is the one in ``rule_selection.json``;
    no row of this curve may be promoted after the fact.
    """
    pb_f1 = prf(masks["pb_calls"], d["snp"])["f1"]
    rows = [{"name": "no safety layer (frozen cascade)", "family": "-",
             **{key: value for key, value in
                evaluate_candidate({"family": "K", "k_min": 10 ** 9}, d, masks, pb_f1).items()
                if key not in ("family", "k_min")}}]
    for rule in candidate_rules():
        rows.append(evaluate_candidate(rule, d, masks, pb_f1))
    return rows


def score_cell(region: dict, rules: dict[str, dict], throughput: dict | None) -> dict:
    """All arms, all strata, controls and curves for one cell."""
    d = frame_of(region)
    routed = route_stage1(d["binomial_llr"])
    escalations = escalations_for(d, rules, routed)
    masks = arm_masks(d)
    pb_f1 = prf(masks["pb_calls"], d["snp"])["f1"]

    entry = {
        "contig": region["contig"],
        "depth_tag": region["spec"]["depth_tag"],
        "regime": region["spec"]["regime"],
        "gc_fraction": region["spec"]["gc"],
        "region": region["region"],
        "loci_total": int(region["labels"].size),
        "loci_scored": int(d["snp"].size),
        "snp": int(d["snp"].sum()),
        "mean_depth": float(region["depth"].mean()),
        "max_depth": float(region["depth"].max()),
        "fraction_routed_to_pb": float(routed.mean()),
        "pb_compute_loci": int(routed.sum()),
        "binomial_only": arm_block(masks["binomial_calls"], d, masks["pb_calls"], pb_f1),
        "pb_only": arm_block(masks["pb_calls"], d, masks["pb_calls"], pb_f1, bootstrap=False),
        "cheap_router_pb": arm_block(masks["cascade_calls"], d, masks["pb_calls"], pb_f1),
        "arms": {},
        "safety": {},
        "controls": {},
    }
    entry["delta_f1_vs_pb"] = entry["cheap_router_pb"]["delta_f1_vs_pb"]
    entry["vs_pb_paired_bootstrap"] = entry["cheap_router_pb"]["vs_pb_paired_bootstrap"]
    entry["verdict_frozen_cascade"] = classify(entry)

    for tag, escalated in escalations.items():
        calls = np.where(routed | escalated, masks["pb_calls"], masks["binomial_calls"])
        block = arm_block(calls, d, masks["pb_calls"], pb_f1)
        block["delta_f1_vs_frozen_cascade"] = block["f1"] - entry["cheap_router_pb"]["f1"]
        entry["arms"][tag] = block
        entry["safety"][tag] = safety_diagnostics(d, masks, escalated)
        entry["controls"][tag] = controls_at_matched_budget(d, masks, escalated)

    entry["by_depth"] = depth_stratified(region)
    entry["by_depth_fine"] = fine_depth_stratified(region)
    entry["by_quality"] = quality_stratified(region)
    entry["by_vaf"] = vaf_stratified(region)
    entry["strata"] = {label: arms_on_mask(region, mask)
                       for label, mask in strata_masks(region).items()}
    entry["by_structure"] = {
        label: arms_on_mask(region, mask)
        for label, mask in structure_masks(region["spec"], region["positions"]).items()}
    entry["safety_by_structure"] = {}
    for label, mask in structure_masks(region["spec"], region["positions"]).items():
        sub = mask[d["frame"]]
        if not sub.any():
            continue
        cell = {}
        for tag, escalated in escalations.items():
            calls = np.where(routed | escalated, masks["pb_calls"], masks["binomial_calls"])
            cell[tag] = {
                "loci": int(sub.sum()),
                "snp": int((d["snp"] & sub).sum()),
                "escalated": int((escalated & sub).sum()),
                "f1": prf(calls[sub], d["snp"][sub])["f1"],
                "frozen_cascade_f1": prf(masks["cascade_calls"][sub], d["snp"][sub])["f1"],
                "pb_only_f1": prf(masks["pb_calls"][sub], d["snp"][sub])["f1"],
                "lost_true_snp_vs_pb": lost_true_snp(calls[sub], masks["pb_calls"][sub],
                                                     d["snp"][sub]),
                "frozen_cascade_lost_true_snp_vs_pb": lost_true_snp(
                    masks["cascade_calls"][sub], masks["pb_calls"][sub], d["snp"][sub]),
            }
        entry["safety_by_structure"][label] = cell

    entry["accuracy_compute_curve"] = accuracy_compute_curve(d, masks)
    if throughput is not None:
        entry["extraction_rates"] = extraction_rates(region)
        entry["projected_wallclock"] = projected_wallclock(
            d, routed, escalations, throughput, entry["extraction_rates"])
    return entry


# -- compute accounting -------------------------------------------------------


def measure_throughput(spec: dict, span_bp: int = 128_000) -> dict:
    """Binomial, safety-layer and PB throughput, measured on this experiment's data.

    The safety layer is timed as the production path would pay for it: over
    **every** locus, because the flip test must be evaluated everywhere the
    router declined, not only where it fires.
    """
    from extract_bench_v12 import poisson_binomial_for_chunk
    from pileup_counts import load_counts
    from read_level_pileup import load_reads
    data = Path(spec["data"])
    fasta = f"data/reference/{spec['contig']}_full.fa"
    bam = str(data / f"{spec['base']}_{spec['depth_tag']}.bam")
    vcf = str(data / f"{spec['base']}.vcf.gz")
    bed = str(data / f"{spec['base']}_highconf.bed")
    span = (spec["start"], spec["start"] + span_bp)

    counts, _ = load_counts(fasta, bam, vcf, bed, span, contig=spec["contig"])
    reads, _ = load_reads(fasta, bam, vcf, bed, span, contig=spec["contig"])

    caller = BinomialVariantCaller(error_rate=0.01)
    clock = time.perf_counter()
    score = caller.score_counts(counts[:, 0:4], counts[:, 9].astype(int))
    binomial_seconds = time.perf_counter() - clock

    # Timed as production pays for it: the full escalation decision over every
    # locus, through the same sparse path the benchmark scores with.
    routed = route_stage1(score.llr)
    clock = time.perf_counter()
    safety_escalate(score.k, score.n, routed, FROZEN_SAFETY_RULE)
    safety_seconds = time.perf_counter() - clock

    # The dense form, for the record: it is what a naive implementation of the
    # same rule would cost, and the gap is reported rather than hidden.
    clock = time.perf_counter()
    epsilon_flip(score.k, score.n,
                 FROZEN_SAFETY_RULE["epsilon_lo"], FROZEN_SAFETY_RULE["epsilon_hi"])
    dense_safety_seconds = time.perf_counter() - clock

    clock = time.perf_counter()
    poisson_binomial_for_chunk(reads, counts)
    pb_seconds = time.perf_counter() - clock

    n = counts.shape[0]
    return {
        "bam": bam, "span": list(span), "sample_loci": int(n),
        "mean_depth": float(counts[:, 6].mean()),
        "binomial_seconds": binomial_seconds,
        "binomial_loci_per_second": n / max(binomial_seconds, 1e-9),
        "safety_seconds": safety_seconds,
        "safety_loci_per_second": n / max(safety_seconds, 1e-9),
        "safety_cost_relative_to_binomial": safety_seconds / max(binomial_seconds, 1e-9),
        "dense_safety_seconds": dense_safety_seconds,
        "dense_safety_cost_relative_to_binomial":
            dense_safety_seconds / max(binomial_seconds, 1e-9),
        "counts_extraction_loci_per_second": None,
        "pb_seconds": pb_seconds,
        "pb_loci_per_second": n / max(pb_seconds, 1e-9),
        "pb_cost_relative_to_safety": pb_seconds / max(safety_seconds, 1e-9),
    }


def projected_wallclock(d: dict, routed: np.ndarray, escalations: dict[str, np.ndarray],
                        throughput: dict, extraction: dict | None = None) -> dict:
    """Cascade, safety-cascade and PB-only runtime, under two accountings.

    **caller-only** — binomial + safety + PB, the convention devlogs 11–14 used
    for every "170–400x" claim. It excludes pileup construction, which both arms
    pay, and therefore flatters PB-only's competitors the least: with the PB
    fraction at ~1e-3, the binomial pass over *every* locus already dominates
    the cascade, so any per-locus stage added next to it is expensive in these
    units even when it is trivial in absolute terms.

    **end-to-end** — the same plus the measured cost of actually building the
    evidence: the count pileup over every locus (both arms), and the read tensor
    only where PB runs. This is what a user waits for. Devlog 15 §7.5's criterion
    A3 says "end-to-end", so A3 is evaluated on this one and the caller-only
    figure is reported next to it rather than instead of it.

    Args:
        d: One scoring frame.
        routed: The frozen router's mask on that frame.
        escalations: Escalation mask per frozen rule.
        throughput: A :func:`measure_throughput` result.
        extraction: Optional measured extraction rates (``counts_loci_per_second``,
            ``reads_loci_per_second``); end-to-end is omitted without it.

    Returns:
        Nested timings and speedups per arm and per accounting.
    """
    n = int(d["snp"].size)
    b_rate = throughput["binomial_loci_per_second"]
    s_rate = throughput["safety_loci_per_second"]
    pb_rate = throughput["pb_loci_per_second"]

    def timings(to_pb: int, with_safety: bool) -> dict:
        caller = n / b_rate + to_pb / pb_rate + (n / s_rate if with_safety else 0.0)
        out = {"pb_loci": to_pb, "pb_compute_fraction": to_pb / max(n, 1),
               "safety_overhead_seconds": (n / s_rate) if with_safety else 0.0,
               "caller_only_seconds": caller}
        if extraction:
            out["end_to_end_seconds"] = (caller + n / extraction["counts_loci_per_second"]
                                         + to_pb / extraction["reads_loci_per_second"])
        return out

    pb_only = {"caller_only_seconds": n / pb_rate}
    if extraction:
        pb_only["end_to_end_seconds"] = (n / pb_rate + n / extraction["counts_loci_per_second"]
                                         + n / extraction["reads_loci_per_second"])
    cascade_t = timings(int(routed.sum()), with_safety=False)

    def speedups(entry: dict, reference: dict) -> dict:
        out = {}
        for key in ("caller_only_seconds", "end_to_end_seconds"):
            if key in entry and key in pb_only:
                unit = key.replace("_seconds", "")
                out[f"{unit}_speedup_vs_pb_only"] = pb_only[key] / max(entry[key], 1e-12)
                out[f"{unit}_speedup_relative_to_frozen_cascade"] = \
                    reference[key] / max(entry[key], 1e-12)
        return out

    out = {"loci": n, "pb_only": pb_only,
           "frozen_cascade": {**cascade_t, "pb_compute_fraction": float(routed.mean()),
                              **speedups(cascade_t, cascade_t)},
           "arms": {}}
    for tag, escalated in escalations.items():
        entry = timings(int((routed | escalated).sum()), with_safety=True)
        out["arms"][tag] = {**entry, **speedups(entry, cascade_t)}
    return out


def extraction_rates(region: dict) -> dict:
    """Measured pileup-construction throughput for one cell, from its cache.

    These timings were produced with four extractions running concurrently, so
    they are *upper bounds* on the per-locus cost and therefore a conservative
    input to the end-to-end speedup (they inflate the shared term both arms pay,
    which pulls every speedup ratio toward 1).
    """
    loci = float(region["labels"].size)
    return {
        "counts_loci_per_second": loci / max(region["counts_time_s"], 1e-9),
        "reads_loci_per_second": loci / max(region["reads_time_s"], 1e-9),
        "measured_counts_seconds": region["counts_time_s"],
        "measured_reads_seconds": region["reads_time_s"],
        "measured_pb_seconds": region["pb_time_s"],
        "note": ("extraction timings come from four concurrent extractions and are "
                 "upper bounds; caller and safety rates come from quiet single-process "
                 "probes"),
    }


# -- pooling ------------------------------------------------------------------


def pool(frames: list[dict], rules: dict[str, dict]) -> dict:
    """Every arm on the concatenation of several cells' scoring frames."""
    keys = ("snp", "binomial_llr", "pb_llr", "k", "n", "vaf", "quality_sum",
            "depth", "mean_base_quality", "mean_mapq")
    d = {key: np.concatenate([f[key] for f in frames]) for key in keys}
    routed = route_stage1(d["binomial_llr"])
    masks = arm_masks(d)
    pb_f1 = prf(masks["pb_calls"], d["snp"])["f1"]
    entry = {
        "loci_scored": int(d["snp"].size),
        "snp": int(d["snp"].sum()),
        "fraction_routed_to_pb": float(routed.mean()),
        "pb_compute_loci": int(routed.sum()),
        "binomial_only": arm_block(masks["binomial_calls"], d, masks["pb_calls"], pb_f1),
        "pb_only": arm_block(masks["pb_calls"], d, masks["pb_calls"], pb_f1, bootstrap=False),
        "cheap_router_pb": arm_block(masks["cascade_calls"], d, masks["pb_calls"], pb_f1),
        "arms": {}, "safety": {}, "controls": {},
    }
    entry["delta_f1_vs_pb"] = entry["cheap_router_pb"]["delta_f1_vs_pb"]
    entry["vs_pb_paired_bootstrap"] = entry["cheap_router_pb"]["vs_pb_paired_bootstrap"]
    entry["verdict_frozen_cascade"] = classify(entry)
    for tag, rule in rules.items():
        escalated = escalation_mask(rule, d, routed)
        calls = np.where(routed | escalated, masks["pb_calls"], masks["binomial_calls"])
        block = arm_block(calls, d, masks["pb_calls"], pb_f1)
        block["delta_f1_vs_frozen_cascade"] = block["f1"] - entry["cheap_router_pb"]["f1"]
        entry["arms"][tag] = block
        entry["safety"][tag] = safety_diagnostics(d, masks, escalated)
        entry["controls"][tag] = controls_at_matched_budget(d, masks, escalated)
        entry[f"verdict_{tag}"] = classify({"delta_f1_vs_pb": block["delta_f1_vs_pb"],
                                            "vs_pb_paired_bootstrap": block["vs_pb_paired_bootstrap"],
                                            "snp": entry["snp"]})
    entry["accuracy_compute_curve"] = accuracy_compute_curve(d, masks)
    return entry


def acceptance(pooled: dict, wallclock: dict, tag: str = "frozen_full") -> dict:
    """The §7.5 criteria A1–A4, evaluated as written. No judgement calls."""
    cascade_lost = pooled["cheap_router_pb"]["lost_true_snp_vs_pb"]
    arm = pooled["arms"][tag]
    safety = pooled["safety"][tag]
    controls = pooled["controls"][tag]
    ci = arm["vs_pb_paired_bootstrap"]["delta_f1_ci95"]

    a1_informative = cascade_lost > 0
    a1 = bool(a1_informative
              and arm["lost_true_snp_vs_pb"] < cascade_lost
              and arm["lost_true_snp_vs_pb"] <= cascade_lost / 2.0)
    a2 = bool(float(ci[0]) > -0.001)
    # §7.5 A3 says "end-to-end", so that is the accounting it is applied to. The
    # stricter caller-only figure (devlogs 11-14's convention, which excludes the
    # pileup construction both arms pay) is recorded next to it, and §14 reports
    # both, because the two do not agree and the disagreement is the point.
    arm_clock = wallclock["arms"][tag]
    end_to_end = arm_clock.get("end_to_end_speedup_relative_to_frozen_cascade")
    caller_only = arm_clock.get("caller_only_speedup_relative_to_frozen_cascade")
    speedup_ratio = end_to_end if end_to_end is not None else caller_only
    a3 = bool(safety["pb_fraction_multiple_vs_frozen"] <= 2.0 and speedup_ratio >= 0.8)
    random_lost = controls["random_escalation"]["min_lost_true_snp"]
    reverse_lost = controls["reverse_confidence_escalation"]["lost_true_snp_vs_pb"]
    a4 = bool(arm["lost_true_snp_vs_pb"] < random_lost
              and arm["lost_true_snp_vs_pb"] < reverse_lost)
    return {
        "arm": tag,
        "frozen_cascade_lost_true_snp": cascade_lost,
        "safety_arm_lost_true_snp": arm["lost_true_snp_vs_pb"],
        "A1_failure_repaired": a1,
        "A1_informative": a1_informative,
        "A2_non_inferior_to_pb": a2,
        "A2_delta_f1_ci95": [float(ci[0]), float(ci[1])],
        "A3_compute_preserved": a3,
        "A3_pb_fraction_multiple": safety["pb_fraction_multiple_vs_frozen"],
        "A3_accounting": "end_to_end" if end_to_end is not None else "caller_only",
        "A3_speedup_relative_to_frozen_cascade": speedup_ratio,
        "A3_end_to_end_speedup_relative_to_frozen_cascade": end_to_end,
        "A3_caller_only_speedup_relative_to_frozen_cascade": caller_only,
        "A3_would_pass_on_caller_only_accounting": bool(
            safety["pb_fraction_multiple_vs_frozen"] <= 2.0
            and (caller_only is None or caller_only >= 0.8)),
        "A4_beats_matched_budget_controls": a4,
        "A4_random_min_lost": random_lost,
        "A4_reverse_confidence_lost": reverse_lost,
    }


def verdict_from(criteria: dict, n_cascade_failures: int) -> str:
    """The §7.6 mapping, applied verbatim."""
    if not criteria["A1_informative"]:
        return "C. NULL / INCONCLUSIVE"
    if not criteria["A3_compute_preserved"]:
        return "E. ARCHITECTURAL FAILURE"
    if not criteria["A1_failure_repaired"] or not criteria["A2_non_inferior_to_pb"]:
        return "D. NEGATIVE"
    if criteria["A4_beats_matched_budget_controls"] and n_cascade_failures >= 3:
        return "A. STRONG POSITIVE"
    return "B. WEAK POSITIVE"


# -- driver -------------------------------------------------------------------


def run_test(measure: bool = True) -> dict:
    specs = test_cells()
    missing = [key for key, spec in specs.items() if not Path(spec["path"]).exists()]
    if missing:
        raise SystemExit(f"missing extractions, refusing to score a partial benchmark: {missing}")
    rules = frozen_rules()
    logger.info("frozen rules: %s", json.dumps(rules))

    loaded, frames = {}, {}
    for key, spec in specs.items():
        logger.info("loading %s", key)
        region = load_cell(key, spec)
        region["positions"] = genomic_positions(spec)
        loaded[key] = region

    audit = verify_caches(loaded)
    AUDIT.write_text(json.dumps(audit, indent=2))
    if not audit["all_ok"]:
        bad = [key for key, cell in audit["cells"].items() if not cell["ok"]]
        raise SystemExit(f"independence audit failed for {bad}; refusing to score")
    logger.info("independence audit passed for %d cells", audit["n_cells"])

    throughputs = {}
    if measure:
        for base in TEST_REGIONS:
            for tag in DEPTH_TAGS:
                spec = specs[f"{base}_{tag}"]
                logger.info("throughput probe %s_%s", base, tag)
                throughputs[f"{base}_{tag}"] = measure_throughput(spec)

    report = {
        "protocol": "History/15_DEVLOG.md",
        "frozen_cascade": {"binomial_threshold": FROZEN_BINOMIAL_THRESHOLD,
                           "pb_threshold": FROZEN_PB_THRESHOLD,
                           "router_cutoff": FROZEN_ROUTER_CUTOFF},
        "frozen_rules": rules,
        "independence_audit": audit,
        "throughput": throughputs,
        "cells": {},
        "pools": {},
    }
    for key, region in loaded.items():
        logger.info("scoring %s", key)
        report["cells"][key] = score_cell(region, rules, throughputs.get(key))
        frames[key] = frame_of(region)

    pools = {"all": list(frames)}
    for base in TEST_REGIONS:
        pools[base] = [f"{base}_{tag}" for tag in DEPTH_TAGS]
    for tag in DEPTH_TAGS:
        pools[f"depth_{tag}"] = [f"{base}_{tag}" for base in TEST_REGIONS]
    pools["segdup_regions"] = [f"{base}_{tag}" for base in TEST_REGIONS
                               for tag in DEPTH_TAGS if base.endswith("_segdup")]
    for name, members in pools.items():
        logger.info("pooling %s", name)
        report["pools"][name] = pool([frames[key] for key in members], rules)
        report["pools"][name]["members"] = members

    if throughputs:
        units = ("caller_only_seconds", "end_to_end_seconds")
        total = {"loci": sum(report["cells"][key]["loci_scored"] for key in frames),
                 "pb_only": {}, "frozen_cascade": {}, "arms": {}}
        for unit in units:
            if not all(unit in report["cells"][key]["projected_wallclock"]["pb_only"]
                       for key in frames):
                continue
            pb_only = sum(report["cells"][key]["projected_wallclock"]["pb_only"][unit]
                          for key in frames)
            cascade_s = sum(
                report["cells"][key]["projected_wallclock"]["frozen_cascade"][unit]
                for key in frames)
            label = unit.replace("_seconds", "")
            total["pb_only"][unit] = pb_only
            total["frozen_cascade"][unit] = cascade_s
            total["frozen_cascade"][f"{label}_speedup_vs_pb_only"] = pb_only / cascade_s
            total["frozen_cascade"][f"{label}_speedup_relative_to_frozen_cascade"] = 1.0
            for rule_tag in rules:
                seconds = sum(
                    report["cells"][key]["projected_wallclock"]["arms"][rule_tag][unit]
                    for key in frames)
                arm = total["arms"].setdefault(rule_tag, {})
                arm[unit] = seconds
                arm[f"{label}_speedup_vs_pb_only"] = pb_only / seconds
                arm[f"{label}_speedup_relative_to_frozen_cascade"] = cascade_s / seconds
        report["pooled_wallclock"] = total
        report["acceptance"] = {
            tag: acceptance(report["pools"]["all"], total, tag) for tag in rules}
        report["verdict"] = {
            tag: verdict_from(criteria,
                              report["pools"]["all"]["cheap_router_pb"]["lost_true_snp_vs_pb"])
            for tag, criteria in report["acceptance"].items()}

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TEST_RESULTS.write_text(json.dumps(report, indent=2))
    return report


def characterise_validation() -> dict:
    """§9: re-derive the failure regime from the cached evidence, independently.

    Devlog 14 §17 characterised the 41 router/PB disagreements and concluded the
    mechanism is epsilon misspecification. This experiment's entire design rests
    on that conclusion, so it is recomputed here from the per-locus evidence
    rather than taken on trust — including the one quantity devlog 14 did not
    report, the epsilon-sensitivity shift that the safety layer keys on.
    """
    d = pooled_validation()
    masks = arm_masks(d)
    cascade_calls, pb_calls = masks["cascade_calls"], masks["pb_calls"]
    disagree = cascade_calls != pb_calls
    snp = d["snp"]

    populations = {
        "router_fn": disagree & snp & pb_calls & ~cascade_calls,
        "router_recovers_fn": disagree & snp & ~pb_calls & cascade_calls,
        "router_fp": disagree & ~snp & cascade_calls,
        "router_avoids_fp": disagree & ~snp & ~cascade_calls,
    }
    shift = (binomial_llr_at_epsilon(d["k"], d["n"], 1e-4)
             - binomial_llr_at_epsilon(d["k"], d["n"], 5e-2))
    fields = {"binomial_llr": d["binomial_llr"], "pb_llr": d["pb_llr"], "vaf": d["vaf"],
              "k": d["k"], "n": d["n"], "depth": d["depth"],
              "mean_base_quality": d["mean_base_quality"], "mean_mapq": d["mean_mapq"],
              "epsilon_sensitivity_shift": shift,
              "router_margin": np.abs(d["binomial_llr"] - FROZEN_BINOMIAL_THRESHOLD)}

    def describe(mask: np.ndarray) -> dict:
        out = {"n": int(mask.sum())}
        if not mask.any():
            return out
        for name, values in fields.items():
            selected = values[mask]
            out[name] = {"min": float(selected.min()),
                         "median": float(np.median(selected)),
                         "max": float(selected.max())}
        out["cells"] = {str(cell): int(count) for cell, count in
                        zip(*np.unique(d["cell"][mask], return_counts=True))}
        return out

    report = {
        "loci": int(snp.size), "snp": int(snp.sum()),
        "routed": int(masks["routed"].sum()),
        "disagreements": int(disagree.sum()),
        "disagreement_rate": float(disagree.mean()),
        "populations": {name: describe(mask) for name, mask in populations.items()},
        "background_unrouted": describe(~masks["routed"]),
    }

    # Is "confidently wrong" visible in the router's own statistic? It must not
    # be, or the router could have caught it by widening (devlog 14 §13).
    unfavourable = populations["router_fn"] | populations["router_fp"]
    report["failures_inside_router_band"] = int((unfavourable & masks["routed"]).sum())
    report["failures_outside_router_band"] = int((unfavourable & ~masks["routed"]).sum())
    report["router_margin_at_failures"] = describe(unfavourable).get("router_margin")

    # Does the epsilon-sensitivity signal separate them from the background?
    for lo, hi in ((1e-3, 5e-2), (1e-4, 5e-2)):
        flip = epsilon_flip(d["k"], d["n"], lo, hi) & ~masks["routed"]
        report[f"flip_{lo:g}_{hi:g}"] = {
            "loci": int(flip.sum()),
            "fraction": float(flip.mean()),
            "captures_router_fn": int((flip & populations["router_fn"]).sum()),
            "captures_router_fp": int((flip & populations["router_fp"]).sum()),
            "captures_router_avoids_fp": int((flip & populations["router_avoids_fp"]).sum()),
            "enrichment_vs_background": float(
                ((flip & unfavourable).sum() / max(int(flip.sum()), 1))
                / max(float((unfavourable & ~masks["routed"]).mean()), 1e-12)),
        }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "validation_characterisation.json").write_text(json.dumps(report, indent=2))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--select", action="store_true",
                        help="run the §7.4 rule selection on validation data only")
    parser.add_argument("--characterise", action="store_true",
                        help="§9 failure-regime characterisation on validation data")
    parser.add_argument("--test", action="store_true",
                        help="score the frozen rule on the untouched test cells")
    parser.add_argument("--no-measure", action="store_true",
                        help="skip the throughput probes (no wall-clock claims)")
    args = parser.parse_args()
    if args.select:
        report = select_rule()
        print(json.dumps({key: report[key] for key in
                          ("outcome", "frozen_cascade_lost_true_snp", "n_feasible_full",
                           "n_feasible_lean", "frozen_full", "frozen_lean")}, indent=2))
        return
    if args.characterise:
        report = characterise_validation()
        print(json.dumps({key: value for key, value in report.items()
                          if key != "populations"}, indent=2)[:4000])
        return
    if args.test:
        report = run_test(measure=not args.no_measure)
        print(json.dumps(report.get("verdict", {"verdict": "not computed (--no-measure)"}),
                         indent=2))
        return
    parser.error("nothing to do: pass --select or --test")


if __name__ == "__main__":
    main()
