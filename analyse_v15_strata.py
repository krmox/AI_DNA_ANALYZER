"""Stratified reporting for the bench_v15 safety arm (devlog 15 §17).

The benchmark's inherited stratifiers (`vaf_stratified`, `quality_stratified`,
`depth_stratified`) score the three *pre-existing* arms only, because they were
written before a fourth arm existed. This script adds the safety arm to the same
bins, on the same cells, with the same frozen rule read from disk.

It is post-hoc **reporting**, not selection: the rule is frozen, the arms are
fixed, and nothing here can change a decision. It exists so the devlog can show
where in (VAF, base quality, depth) space the layer fires and what it does there.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from bench_v13_robustness import VAF_BINS
from bench_v15_safety import (DEPTH_TAGS, TEST_REGIONS, arm_masks, escalation_mask,
                              frame_of, frozen_rules, load_cell, lost_true_snp,
                              test_cells)
from cascade import route_stage1
from evaluate_binomial_baseline import prf
from robustness_benchmark import DEPTH_BINS

OUT = Path("results/bench_v15/stratified.json")
QUALITY_BINS = [(0, 20), (20, 25), (25, 30), (30, 35), (35, 60)]


def _bins(values: np.ndarray, bins, prefix: str) -> dict[str, np.ndarray]:
    """Boolean mask per half-open bin, labelled for the devlog table."""
    return {f"{prefix}{lo:g}-{hi:g}": (values >= lo) & (values < hi) for lo, hi in bins}


def strata(d: dict) -> dict[str, np.ndarray]:
    return {
        **_bins(d["vaf"], VAF_BINS, "vaf "),
        **_bins(d["mean_base_quality"], QUALITY_BINS, "meanQ "),
        **_bins(d["depth"], [(lo, hi + 1) for lo, hi in DEPTH_BINS], "depth "),
    }


def main() -> None:
    rules = frozen_rules()
    specs = test_cells()
    parts: dict[str, list] = {}
    for key, spec in specs.items():
        d = frame_of(load_cell(key, spec))
        for field, value in d.items():
            if field != "frame":
                parts.setdefault(field, []).append(value)
    d = {field: np.concatenate(values) for field, values in parts.items()}

    routed = route_stage1(d["binomial_llr"])
    masks = arm_masks(d)
    escalations = {tag: escalation_mask(rule, d, routed) for tag, rule in rules.items()}
    calls = {tag: np.where(routed | mask, masks["pb_calls"], masks["binomial_calls"])
             for tag, mask in escalations.items()}

    report = {"frozen_rules": rules, "loci": int(d["snp"].size), "strata": {}}
    for label, mask in strata(d).items():
        if not mask.any():
            continue
        snp = d["snp"][mask]
        entry = {
            "loci": int(mask.sum()),
            "snp": int(snp.sum()),
            "routed": int(routed[mask].sum()),
            "pb_only_f1": prf(masks["pb_calls"][mask], snp)["f1"],
            "frozen_cascade_f1": prf(masks["cascade_calls"][mask], snp)["f1"],
            "frozen_cascade_lost_true_snp": lost_true_snp(
                masks["cascade_calls"][mask], masks["pb_calls"][mask], snp),
            "arms": {},
        }
        for tag in rules:
            sub = calls[tag][mask]
            entry["arms"][tag] = {
                "escalated": int(escalations[tag][mask].sum()),
                "f1": prf(sub, snp)["f1"],
                "lost_true_snp": lost_true_snp(sub, masks["pb_calls"][mask], snp),
                "calls_changed_vs_cascade": int(
                    np.count_nonzero(sub != masks["cascade_calls"][mask])),
            }
        report["strata"][label] = entry

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2))

    tag = "frozen_full"
    print(f"{'stratum':18s} {'loci':>10s} {'snp':>6s} {'esc':>5s} {'F1 B':>9s} "
          f"{'F1 C':>9s} {'F1 D':>9s} {'lost C':>6s} {'lost D':>6s} {'changed':>7s}")
    for label, entry in report["strata"].items():
        arm = entry["arms"][tag]
        print(f"{label:18s} {entry['loci']:10,d} {entry['snp']:6d} {arm['escalated']:5d} "
              f"{entry['pb_only_f1']:9.6f} {entry['frozen_cascade_f1']:9.6f} {arm['f1']:9.6f} "
              f"{entry['frozen_cascade_lost_true_snp']:6d} {arm['lost_true_snp']:6d} "
              f"{arm['calls_changed_vs_cascade']:7d}")


if __name__ == "__main__":
    main()
