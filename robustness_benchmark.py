"""Final robustness benchmark: frozen cheap-router -> PB cascade, broader/harder data.

This is a VALIDATION experiment, not an optimization exercise. It reuses the
frozen production router (``cascade.py`` / ``cheap_router.py``, both copied
byte-for-byte from the validated worktree commit ``c06deec`` on branch
``worktree-agent-aaff295c27efb0e79``, see History/11_DEVLOG.md sec 2) and
evaluates it, unmodified, on genomic regions that were never used to select
the router's cutoff.

Nothing in this file fits, tunes, or adjusts the router. The only "learned"
artifacts touched are the frozen constants imported from ``cascade``:

    FROZEN_BINOMIAL_THRESHOLD = 7.0
    FROZEN_PB_THRESHOLD       = 10.5
    FROZEN_ROUTER_CUTOFF      = 5.411872376933351

Regions (see History/11_DEVLOG.md sec 4-5 for the leakage audit and
selection rule, pre-registered before any result in this file was read):

    gen14_near_40_44M        chr21:40,000,000-44,000,000   (3.37 Mb)
    gen14_far_13_17M         chr21:13,000,000-17,000,000   (3.57 Mb)
    gen14_far_17_21M         chr21:17,000,000-21,000,000   (3.87 Mb)
    gen14_far_21_25M         chr21:21,000,000-25,000,000   (3.88 Mb)
    gen14_far_25_29M         chr21:25,000,000-29,000,000   (3.84 Mb)
    test_15x                 chr21:30,000,000-30,469,999   (0.47 Mb, ~15x)
    test_30x                 chr21:30,000,000-30,469,999   (0.47 Mb, ~30x, same
                              loci as test_15x at different downsample depth)

The router-tuning region (devlog 13) was chr21:32,000,000-40,000,000 (big15x
train/validation windows fall at 32-38 Mb within it). None of the seven
regions above overlap that span.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np

from binomial_baseline import BinomialVariantCaller
from cascade import FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD, FROZEN_ROUTER_CUTOFF
from cheap_router import binomial_llr as compute_binomial_llr
from cheap_router import measure_throughput
from config import LABEL_SNP
from evaluate_binomial_baseline import prf
from evaluate_quality_error import paired_bootstrap

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

RESULTS_DIR = Path("results/robustness_v11")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CASCADE_SOURCE = {
    "branch": "worktree-agent-aaff295c27efb0e79",
    "commit": "c06deec",
    "file": "cascade.py (copied verbatim into main working tree for this benchmark)",
    "frozen_binomial_threshold": FROZEN_BINOMIAL_THRESHOLD,
    "frozen_pb_threshold": FROZEN_PB_THRESHOLD,
    "frozen_router_cutoff": FROZEN_ROUTER_CUTOFF,
}

DEPTH_BINS = [(1, 4), (5, 9), (10, 14), (15, 19), (20, 29), (30, 10_000)]
CUTOFF_SWEEP = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
               5.0, 5.411872376933351, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0]

GEN14_REGIONS = {
    "gen14_near_40_44M": "cache/gen14_near_40_44M.npz",
    "gen14_far_13_17M": "cache/gen14_far_13000000_17000000.npz",
    "gen14_far_17_21M": "cache/gen14_far_17000000_21000000.npz",
    "gen14_far_21_25M": "cache/gen14_far_21000000_25000000.npz",
    "gen14_far_25_29M": "cache/gen14_far_25000000_29000000.npz",
}


def accuracy_block(calls: np.ndarray, truth: np.ndarray) -> dict:
    stats = prf(calls, truth)
    true_negative = int(np.sum(~calls & ~truth))
    total = truth.size
    sensitivity = stats["tp"] / max(stats["tp"] + stats["fn"], 1)
    specificity = true_negative / max(true_negative + stats["fp"], 1)
    stats["tn"] = true_negative
    stats["accuracy"] = (stats["tp"] + true_negative) / max(total, 1)
    stats["balanced_accuracy"] = 0.5 * (sensitivity + specificity)
    return stats


def load_gen14(name: str, path: str) -> dict:
    blob = np.load(path, allow_pickle=True)
    counts = blob["counts"]
    return {
        "name": name,
        "counts": counts,
        "binomial_llr": blob["binomial_llr"].astype(np.float64),
        "pb_llr": blob["pb_llr"].astype(np.float64),
        "labels": blob["labels"],
        "depth": blob["depth"].astype(np.float64),
        "mean_base_quality": np.where(blob["depth"] > 0, counts[:, 7] / np.maximum(blob["depth"], 1), 0.0),
        "mean_mapq": np.where(blob["depth"] > 0, counts[:, 8] / np.maximum(blob["depth"], 1), 0.0),
        "region": [int(x) for x in blob["region"]],
        "pb_time_s": float(blob["timing"][2]) if "timing" in blob else None,
        "counts_time_s": float(blob["timing"][0]) if "timing" in blob else None,
    }


def load_15x30x(name: str, counts_path: str, evidence_path: str) -> dict:
    counts_blob = np.load(counts_path, allow_pickle=True)
    evidence_blob = np.load(evidence_path, allow_pickle=True)
    counts = counts_blob["counts"]
    depth = evidence_blob["depth"].astype(np.float64)
    caller = BinomialVariantCaller(error_rate=0.01)
    binomial_llr = compute_binomial_llr(counts, caller)
    labels = counts_blob["labels"]
    assert np.array_equal(labels, evidence_blob["labels"]), "label mismatch between counts/evidence caches"
    return {
        "name": name,
        "counts": counts,
        "binomial_llr": binomial_llr,
        "pb_llr": evidence_blob["pb_llr"].astype(np.float64),
        "labels": labels,
        "depth": depth,
        "mean_base_quality": np.where(depth > 0, counts[:, 7] / np.maximum(depth, 1), 0.0),
        "mean_mapq": np.where(depth > 0, counts[:, 8] / np.maximum(depth, 1), 0.0),
        "region": [int(x) for x in counts_blob["region"]],
        "pb_time_s": None,
        "counts_time_s": None,
    }


def load_all_regions() -> dict:
    regions = {}
    for name, path in GEN14_REGIONS.items():
        regions[name] = load_gen14(name, path)
    regions["test_15x"] = load_15x30x("test_15x", "cache/test_15x_counts.npz",
                                      "cache/test_15x_evidence.npz")
    regions["test_30x"] = load_15x30x("test_30x", "cache/test_30x_counts.npz",
                                      "cache/test_30x_evidence.npz")
    return regions


def route_mask(binomial_llr: np.ndarray, cutoff: float) -> np.ndarray:
    return np.abs(binomial_llr - FROZEN_BINOMIAL_THRESHOLD) <= cutoff


def evaluate_arms(region: dict) -> dict:
    labels = region["labels"]
    frame = (labels == 0) | (labels == LABEL_SNP)
    snp = labels[frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]
    depth = region["depth"][frame]
    quality = region["mean_base_quality"][frame]

    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    router_pb_calls = np.where(routed, pb_calls, binomial_calls)

    entry = {
        "loci_total": int(labels.size),
        "loci_scored": int(frame.sum()),
        "snp": int(snp.sum()),
        "region": region["region"],
        "binomial_only": accuracy_block(binomial_calls, snp),
        "pb_only": accuracy_block(pb_calls, snp),
        "cheap_router_pb": accuracy_block(router_pb_calls, snp),
        "fraction_routed_to_pb": float(routed.mean()),
        "pb_compute_loci": int(routed.sum()),
        "disagreement_vs_pb": int(np.count_nonzero(router_pb_calls != pb_calls)),
        "disagreement_rate": float(np.mean(router_pb_calls != pb_calls)),
        "vs_pb_paired_bootstrap": paired_bootstrap(router_pb_calls, pb_calls, snp),
        "vs_pb_binomial_only_paired_bootstrap": paired_bootstrap(binomial_calls, pb_calls, snp),
    }
    entry["delta_f1_vs_pb"] = entry["cheap_router_pb"]["f1"] - entry["pb_only"]["f1"]
    entry["delta_accuracy_vs_pb"] = (entry["cheap_router_pb"]["accuracy"]
                                     - entry["pb_only"]["accuracy"])

    # PB disagreement margins: |PB_LLR - FROZEN_PB_THRESHOLD| at disagreement loci
    disagree_mask = router_pb_calls != pb_calls
    entry["disagreement_pb_margins"] = [float(x) for x in
                                        np.sort(np.abs(pb_llr[disagree_mask] - FROZEN_PB_THRESHOLD))]

    # Router coverage distribution stats
    margin = np.abs(binomial_llr - FROZEN_BINOMIAL_THRESHOLD)
    entry["router_margin_percentiles"] = {
        str(p): float(np.percentile(margin, p)) for p in (1, 5, 25, 50, 75, 95, 99)
    }
    return entry


def depth_stratified(region: dict) -> dict:
    labels = region["labels"]
    frame = (labels == 0) | (labels == LABEL_SNP)
    snp = labels[frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]
    depth = region["depth"][frame]
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    router_pb_calls = np.where(routed, pb_calls, binomial_calls)

    out = {}
    for lo, hi in DEPTH_BINS:
        mask = (depth >= lo) & (depth <= hi)
        n = int(mask.sum())
        key = f"{lo}-{hi}x"
        if n == 0:
            out[key] = {"loci": 0}
            continue
        out[key] = {
            "loci": n,
            "snp": int(snp[mask].sum()),
            "binomial_only": accuracy_block(binomial_calls[mask], snp[mask]),
            "pb_only": accuracy_block(pb_calls[mask], snp[mask]),
            "cheap_router_pb": accuracy_block(router_pb_calls[mask], snp[mask]),
            "fraction_routed_to_pb": float(routed[mask].mean()),
        }
    return out


def quality_stratified(region: dict) -> dict:
    labels = region["labels"]
    frame = (labels == 0) | (labels == LABEL_SNP)
    snp = labels[frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]
    quality = region["mean_base_quality"][frame]
    depth = region["depth"][frame]
    has_reads = depth > 0
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    router_pb_calls = np.where(routed, pb_calls, binomial_calls)

    edges = [(0, 20), (20, 25), (25, 30), (30, 35), (35, 60)]
    out = {}
    for lo, hi in edges:
        mask = has_reads & (quality >= lo) & (quality < hi)
        n = int(mask.sum())
        key = f"q{lo}-{hi}"
        if n == 0:
            out[key] = {"loci": 0}
            continue
        out[key] = {
            "loci": n,
            "snp": int(snp[mask].sum()),
            "binomial_only": accuracy_block(binomial_calls[mask], snp[mask]),
            "pb_only": accuracy_block(pb_calls[mask], snp[mask]),
            "cheap_router_pb": accuracy_block(router_pb_calls[mask], snp[mask]),
            "fraction_routed_to_pb": float(routed[mask].mean()),
        }
    return out


def cutoff_sweep(region: dict) -> list:
    labels = region["labels"]
    frame = (labels == 0) | (labels == LABEL_SNP)
    snp = labels[frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    pb_f1 = prf(pb_calls, snp)["f1"]

    rows = []
    for cutoff in CUTOFF_SWEEP:
        routed = route_mask(binomial_llr, cutoff)
        calls = np.where(routed, pb_calls, binomial_calls)
        stats = prf(calls, snp)
        rows.append({
            "cutoff": cutoff,
            "is_frozen_operating_point": bool(abs(cutoff - FROZEN_ROUTER_CUTOFF) < 1e-9),
            "fraction_routed_to_pb": float(routed.mean()),
            "f1": stats["f1"], "precision": stats["precision"], "recall": stats["recall"],
            "tp": stats["tp"], "fp": stats["fp"], "fn": stats["fn"],
            "delta_f1_vs_pb": stats["f1"] - pb_f1,
            "disagreement_vs_pb": int(np.count_nonzero(calls != pb_calls)),
        })
    return rows


def controls(region: dict, seeds=(0, 1, 2, 3, 4)) -> dict:
    """Random and reverse-confidence routing at matched PB coverage."""
    labels = region["labels"]
    frame = (labels == 0) | (labels == LABEL_SNP)
    snp = labels[frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    pb_f1 = prf(pb_calls, snp)["f1"]

    routed_frozen = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    k = int(routed_frozen.sum())
    n = binomial_llr.size

    def score_mask(mask):
        calls = np.where(mask, pb_calls, binomial_calls)
        stats = prf(calls, snp)
        return {"f1": stats["f1"], "delta_f1_vs_pb": stats["f1"] - pb_f1,
                "tp": stats["tp"], "fp": stats["fp"], "fn": stats["fn"],
                "pb_coverage": int(mask.sum())}

    random_runs = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        idx = rng.choice(n, size=k, replace=False)
        mask = np.zeros(n, dtype=bool)
        mask[idx] = True
        random_runs.append(score_mask(mask))

    # reverse-confidence: route the k loci FARTHEST from the binomial decision
    # boundary (the router's ranking inverted).
    margin = np.abs(binomial_llr - FROZEN_BINOMIAL_THRESHOLD)
    reverse_idx = np.argsort(margin)[::-1][:k]
    reverse_mask = np.zeros(n, dtype=bool)
    reverse_mask[reverse_idx] = True

    return {
        "matched_pb_coverage": k,
        "frozen_router": score_mask(routed_frozen),
        "random_routing": {
            "runs": random_runs,
            "mean_f1": float(np.mean([r["f1"] for r in random_runs])),
            "std_f1": float(np.std([r["f1"] for r in random_runs])),
        },
        "reverse_confidence_routing": score_mask(reverse_mask),
    }


def main() -> None:
    regions = load_all_regions()
    report = {"cascade_source": CASCADE_SOURCE, "regions": {}}

    pooled_labels, pooled_binomial, pooled_pb, pooled_router = [], [], [], []

    for name, region in regions.items():
        logger.info("evaluating region %s", name)
        primary = evaluate_arms(region)
        depth_strat = depth_stratified(region)
        quality_strat = quality_stratified(region)
        sweep = cutoff_sweep(region)
        ctrl = controls(region)

        # throughput measurement on this hardware, this run
        counts_sample = region["counts"][:200_000] if region["counts"].shape[0] > 200_000 else region["counts"]
        binomial_throughput = measure_throughput(lambda c: compute_binomial_llr(c), counts_sample)

        report["regions"][name] = {
            "primary": primary,
            "depth_stratified": depth_strat,
            "quality_stratified": quality_strat,
            "cutoff_sweep": sweep,
            "controls": ctrl,
            "pb_extraction_time_s_from_original_run": region["pb_time_s"],
            "binomial_throughput_this_machine": binomial_throughput,
        }

        labels = region["labels"]
        frame = (labels == 0) | (labels == LABEL_SNP)
        pooled_labels.append(labels[frame] == LABEL_SNP)
        pooled_binomial.append(region["binomial_llr"][frame] >= FROZEN_BINOMIAL_THRESHOLD)
        pooled_pb.append(region["pb_llr"][frame] >= FROZEN_PB_THRESHOLD)
        routed = route_mask(region["binomial_llr"][frame], FROZEN_ROUTER_CUTOFF)
        pooled_router.append(np.where(routed, region["pb_llr"][frame] >= FROZEN_PB_THRESHOLD,
                                      region["binomial_llr"][frame] >= FROZEN_BINOMIAL_THRESHOLD))

    pooled_labels_a = np.concatenate(pooled_labels)
    pooled_binomial_a = np.concatenate(pooled_binomial)
    pooled_pb_a = np.concatenate(pooled_pb)
    pooled_router_a = np.concatenate(pooled_router)

    report["pooled"] = {
        "loci_scored": int(pooled_labels_a.size),
        "snp": int(pooled_labels_a.sum()),
        "binomial_only": accuracy_block(pooled_binomial_a, pooled_labels_a),
        "pb_only": accuracy_block(pooled_pb_a, pooled_labels_a),
        "cheap_router_pb": accuracy_block(pooled_router_a, pooled_labels_a),
        "vs_pb_paired_bootstrap": paired_bootstrap(pooled_router_a, pooled_pb_a, pooled_labels_a),
        "disagreement_vs_pb": int(np.count_nonzero(pooled_router_a != pooled_pb_a)),
        "disagreement_rate": float(np.mean(pooled_router_a != pooled_pb_a)),
    }
    report["pooled"]["delta_f1_vs_pb"] = (report["pooled"]["cheap_router_pb"]["f1"]
                                          - report["pooled"]["pb_only"]["f1"])

    out_path = RESULTS_DIR / "robustness_results.json"
    out_path.write_text(json.dumps(report, indent=2, default=float))
    logger.info("wrote %s", out_path)

    print("\n=== POOLED RESULT ===")
    print(json.dumps(report["pooled"], indent=2, default=float))
    for name, entry in report["regions"].items():
        p = entry["primary"]
        print(f"\n-- {name}: loci={p['loci_scored']:,} snp={p['snp']} "
              f"pb_only.F1={p['pb_only']['f1']:.4f} "
              f"router_pb.F1={p['cheap_router_pb']['f1']:.4f} "
              f"dF1={p['delta_f1_vs_pb']:+.5f} "
              f"routed={p['fraction_routed_to_pb']*100:.3f}% "
              f"disagree={p['disagreement_vs_pb']}")


if __name__ == "__main__":
    main()
