"""bench_v12 Stages 2 and 3 — PB -> Mamba, and the complete three-stage cascade.

Protocol pre-registered in ``History/12_DEVLOG.md`` §6.2-6.3. Nothing here is
fitted: the PB threshold, the router cutoff, the Mamba decision rule
(``score >= FROZEN_PB_THRESHOLD``), the candidate-pool margin (13.5) and the
Stage-2 sweep grid are all imported or copied verbatim from frozen sources.

Both stages share one expensive artefact -- the Mamba score for every locus in
every retained window -- so they run in one process and one output file.

Leakage note: the candidate pool and every routing mask are functions of caller
scores only. No truth label enters any routing decision.
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import torch

from cascade import (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD,
                     FROZEN_ROUTER_CUTOFF, run_cascade)
from config import LABEL_SNP
from evaluate_pb_residual import model_scores
from evaluate_quality_error import paired_bootstrap
from residual_metrics import depth_bin_masks, prf
from robustness_benchmark import accuracy_block, route_mask
from train_raw_pileup import SEQ_LEN

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

#: Default dataset; overridable with --region-npz so the identical, frozen
#: protocol can be applied to a second region without editing any threshold.
REGION_NPZ = "cache/bench_v12/new31_32M_15x.npz"
OUT = "results/bench_v12/stage23_results.json"
SCORE_CACHE = "results/bench_v12/mamba_scores_new31_32M_15x.npz"

#: Adopted unchanged from ``History/16_DEVLOG.md``; not re-searched here.
CANDIDATE_MARGIN = 13.5
THRESHOLD_SWEEP = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 13.5]

WORKTREE = Path(".claude/worktrees/agent-aaff295c27efb0e79")
CHECKPOINTS = {
    "big15x_seedA": WORKTREE / "checkpoints_big15x/best_big15x_seedA.pt",
    "big15x_seedB": WORKTREE / "checkpoints_big15x/best_big15x_seedB.pt",
    "big15x_seedC": WORKTREE / "checkpoints_big15x/best_big15x_seedC.pt",
    "big15x_shuffled_control": WORKTREE / "checkpoints_big15x/best_big15x_shuffled.pt",
}
#: Trained on chr21:31-33 Mb, i.e. ON this benchmark region. Reported only in a
#: separate, explicitly-marked block; excluded from every verdict.
LEAKED_CHECKPOINTS = {
    "LEAKED_pb_residual_v1_seedA": Path("checkpoints_pb_residual/best_v1_repr_seedA.pt"),
    "LEAKED_pb_residual_v1_seedB": Path("checkpoints_pb_residual/best_v1_repr_seedB.pt"),
    "LEAKED_pb_residual_v2_gate_seedA": Path("checkpoints_pb_residual/best_v2_gate_seedA.pt"),
}


def load_region(path: str = REGION_NPZ) -> dict:
    blob = np.load(path, allow_pickle=True)
    positions = blob["positions"].astype(np.int64)
    window_positions = blob["pool_window_positions"].astype(np.int64)
    index = np.searchsorted(positions, window_positions)
    assert np.array_equal(positions[index], window_positions), "window/locus index mismatch"
    assert index.size % SEQ_LEN == 0, "retained features are not window aligned"
    stored_scores = {k.split("score__", 1)[1]: blob[k].astype(np.float64)
                     for k in blob.files if k.startswith("score__")}
    stored_thresholds = {k.split("own_threshold__", 1)[1]: float(blob[k][0])
                         for k in blob.files if k.startswith("own_threshold__")}
    return {
        "labels": blob["labels"],
        "binomial_llr": blob["binomial_llr"].astype(np.float64),
        "pb_llr": blob["pb_llr"].astype(np.float64),
        "depth": blob["depth"].astype(np.float64),
        "features": blob["pool_window_features"] if "pool_window_features" in blob.files else None,
        "stored_scores": stored_scores, "stored_thresholds": stored_thresholds,
        "feature_index": index,
        "n_loci": int(blob["labels"].size),
    }


def mamba_scores(region: dict, cache_path: str = SCORE_CACHE) -> dict:
    """Score every retained window with every checkpoint (cached to disk)."""
    if region["stored_scores"]:
        logger.info("using per-locus Mamba scores stored at extraction time")
        out = {}
        for name, score in region["stored_scores"].items():
            out[name] = score
            out[name + "__own_threshold"] = np.asarray([region["stored_thresholds"][name]])
            out[name + "__seconds"] = np.asarray([float("nan"), float("nan")])
        return out

    cache = Path(cache_path)
    if cache.exists():
        blob = np.load(cache)
        logger.info("reusing cached Mamba scores %s", cache)
        return {k: blob[k] for k in blob.files}

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    features = region["features"]
    llr = region["pb_llr"][region["feature_index"]]
    out = {}
    for name, path in {**CHECKPOINTS, **LEAKED_CHECKPOINTS}.items():
        start = time.perf_counter()
        score, diagnostics, seconds = model_scores(Path(path), features, llr, device)
        out[name] = score
        out[name + "__seconds"] = np.asarray([seconds, time.perf_counter() - start])
        out[name + "__own_threshold"] = np.asarray([diagnostics["validation_threshold"]])
        logger.info("scored %s in %.1fs (%d loci)", name, seconds, score.size)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, **out)
    return out


def breakdown(pb_calls, final_calls, truth, routed) -> dict:
    """Paired change accounting restricted to the Mamba-routed subset."""
    r = routed
    pb_ok = pb_calls[r] == truth[r]
    new_ok = final_calls[r] == truth[r]
    return {
        "routed": int(r.sum()),
        "decisions_changed": int(np.count_nonzero(pb_calls[r] != final_calls[r])),
        "both_correct": int(np.count_nonzero(pb_ok & new_ok)),
        "mamba_fixes": int(np.count_nonzero(~pb_ok & new_ok)),
        "mamba_introduces": int(np.count_nonzero(pb_ok & ~new_ok)),
        "both_wrong": int(np.count_nonzero(~pb_ok & ~new_ok)),
    }


def stage2(region: dict, scores: dict) -> dict:
    """PB-only vs Mamba-on-pool vs PB->Mamba, on the PB-uncertain candidate pool."""
    labels = region["labels"][region["feature_index"]]
    pb_llr = region["pb_llr"][region["feature_index"]]
    depth = region["depth"][region["feature_index"]]
    frame = (labels == 0) | (labels == LABEL_SNP)
    pool = np.abs(pb_llr - FROZEN_PB_THRESHOLD) <= CANDIDATE_MARGIN
    sel = frame & pool

    truth = labels[sel] == LABEL_SNP
    pb = pb_llr[sel] >= FROZEN_PB_THRESHOLD
    margin = np.abs(pb_llr[sel] - FROZEN_PB_THRESHOLD)
    depth_sel = depth[sel]

    result = {
        "candidate_margin": CANDIDATE_MARGIN,
        "pool_loci_in_retained_windows": int(pool.sum()),
        "pool_loci_scored": int(sel.sum()),
        "pool_snp": int(truth.sum()),
        "retained_window_loci": int(labels.size),
        "region_loci": region["n_loci"],
        "pb_only": accuracy_block(pb, truth),
        "sweep": {}, "controls": {}, "by_depth": {},
    }

    result["own_threshold_diagnostic"] = {}
    for name in list(CHECKPOINTS) + list(LEAKED_CHECKPOINTS):
        score = scores[name][sel]
        mamba_calls = score >= FROZEN_PB_THRESHOLD

        # Secondary, pre-declared diagnostic: each checkpoint's own
        # validation-selected threshold (selected on its own validation split,
        # never on this region). Declared before any result was read; the
        # primary decision rule remains FROZEN_PB_THRESHOLD, per devlog 16.
        own = float(scores[name + "__own_threshold"][0])
        own_calls = score >= own
        own_block = accuracy_block(own_calls, truth)
        own_block["threshold"] = own
        own_block["delta_f1_vs_pb"] = own_block["f1"] - result["pb_only"]["f1"]
        own_block["paired_bootstrap"] = paired_bootstrap(own_calls, pb, truth)
        result["own_threshold_diagnostic"][name] = own_block

        entry = {"mamba_only_on_pool": accuracy_block(mamba_calls, truth)}
        entry["mamba_only_on_pool"]["delta_f1_vs_pb"] = (
            entry["mamba_only_on_pool"]["f1"] - result["pb_only"]["f1"])
        entry["mamba_only_on_pool"]["paired_bootstrap"] = paired_bootstrap(
            mamba_calls, pb, truth)
        for t in THRESHOLD_SWEEP:
            routed = margin < t
            final = np.where(routed, mamba_calls, pb)
            block = accuracy_block(final, truth)
            block["delta_f1_vs_pb"] = block["f1"] - result["pb_only"]["f1"]
            block["paired_bootstrap"] = paired_bootstrap(final, pb, truth)
            block["breakdown"] = breakdown(pb, final, truth, routed)
            block["fraction_of_pool_routed"] = float(routed.mean())
            entry[f"t={t}"] = block
        result["sweep"][name] = entry

        # depth stratification at the well-powered point (t = whole pool)
        by_depth = {}
        for label, mask in depth_bin_masks(depth_sel):
            if mask.sum() == 0:
                continue
            by_depth[label] = {
                "loci": int(mask.sum()), "snp": int(truth[mask].sum()),
                "pb_only_f1": prf(pb[mask], truth[mask])["f1"],
                "mamba_f1": prf(mamba_calls[mask], truth[mask])["f1"],
            }
        result["by_depth"][name] = by_depth

    # ---- controls at matched coverage, for the primary point t = 0.5 --------
    for t in (0.5, 13.5):
        routed = margin < t
        k = int(routed.sum())
        n = truth.size
        seed_runs = []
        for seed in range(5):
            rng = np.random.default_rng(seed)
            mask = np.zeros(n, dtype=bool)
            mask[rng.choice(n, size=k, replace=False)] = True
            calls = np.where(mask, scores["big15x_seedA"][sel] >= FROZEN_PB_THRESHOLD, pb)
            seed_runs.append(prf(calls, truth)["f1"])
        reverse = np.zeros(n, dtype=bool)
        reverse[np.argsort(margin)[::-1][:k]] = True
        reverse_calls = np.where(reverse, scores["big15x_seedA"][sel] >= FROZEN_PB_THRESHOLD, pb)
        shuffled_calls = np.where(routed,
                                  scores["big15x_shuffled_control"][sel] >= FROZEN_PB_THRESHOLD, pb)
        result["controls"][f"t={t}"] = {
            "matched_coverage": k,
            "pb_only_f1": result["pb_only"]["f1"],
            "frozen_routing_seedA_f1": result["sweep"]["big15x_seedA"][f"t={t}"]["f1"],
            "random_routing_mean_f1": float(np.mean(seed_runs)),
            "random_routing_std_f1": float(np.std(seed_runs)),
            "reverse_confidence_f1": prf(reverse_calls, truth)["f1"],
            "shuffled_checkpoint_f1": prf(shuffled_calls, truth)["f1"],
        }
    return result


def stage3(region: dict, scores: dict) -> dict:
    """Complete cascade on the whole region frame."""
    labels = region["labels"]
    frame = (labels == 0) | (labels == LABEL_SNP)
    truth = labels[frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]

    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    routed_pb = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    router_pb_calls = np.where(routed_pb, pb_calls, binomial_calls)

    # Mamba score defined only inside retained windows; PB elsewhere. Which loci
    # have a score is a function of PB LLR alone (the pool rule), never labels.
    has_score = np.zeros(labels.size, dtype=bool)
    has_score[region["feature_index"]] = True

    result = {
        "loci_scored": int(frame.sum()), "snp": int(truth.sum()),
        "arms": {
            "binomial_only": accuracy_block(binomial_calls, truth),
            "pb_only": accuracy_block(pb_calls, truth),
            "cheap_router_pb": accuracy_block(router_pb_calls, truth),
        },
        "stage_fractions": {}, "compute": {}, "sweep": {},
        "mamba_only_note": ("a whole-frame Mamba-only arm is NOT computable: the "
                            "44-channel features exist only for windows containing a "
                            "PB-uncertain locus. The 'mamba_on_pool_pb_elsewhere' arm "
                            "below is the closest computable stand-in and is labelled "
                            "as such."),
    }
    result["arms"]["cheap_router_pb"]["delta_f1_vs_pb"] = (
        result["arms"]["cheap_router_pb"]["f1"] - result["arms"]["pb_only"]["f1"])
    result["arms"]["cheap_router_pb"]["paired_bootstrap"] = paired_bootstrap(
        router_pb_calls, pb_calls, truth)

    for name in CHECKPOINTS:
        # PB LLR stands in wherever no Mamba score exists (no retained window).
        score_full = np.zeros(labels.size)
        score_full[region["feature_index"]] = scores[name]
        full_score = np.where(has_score[frame], score_full[frame], pb_llr)

        entry = {"mamba_on_pool_pb_elsewhere": accuracy_block(
            full_score >= FROZEN_PB_THRESHOLD, truth)}
        entry["mamba_on_pool_pb_elsewhere"]["delta_f1_vs_pb"] = (
            entry["mamba_on_pool_pb_elsewhere"]["f1"] - result["arms"]["pb_only"]["f1"])

        for t in THRESHOLD_SWEEP:
            final, stage_used = run_cascade(binomial_llr, pb_llr, full_score, t)
            block = accuracy_block(final, truth)
            block["delta_f1_vs_pb"] = block["f1"] - result["arms"]["pb_only"]["f1"]
            block["paired_bootstrap"] = paired_bootstrap(final, pb_calls, truth)
            block["disagreement_vs_pb"] = int(np.count_nonzero(final != pb_calls))
            n_mamba = int((stage_used == "mamba").sum())
            n_pb = int((stage_used == "pb").sum())
            block["stage_counts"] = {
                "binomial": int((stage_used == "binomial").sum()),
                "pb": n_pb, "mamba": n_mamba}
            block["fraction_reaching_pb"] = float((n_pb + n_mamba) / truth.size)
            block["fraction_reaching_mamba"] = float(n_mamba / truth.size)
            # window accounting: Mamba operates on whole SEQ_LEN windows
            mamba_mask = np.zeros(labels.size, dtype=bool)
            mamba_mask[np.nonzero(frame)[0][stage_used == "mamba"]] = True
            windows = mamba_mask.reshape(-1, SEQ_LEN).any(axis=1)
            block["windows_touched"] = int(windows.sum())
            block["window_amplified_loci"] = int(windows.sum() * SEQ_LEN)
            block["window_amplification_factor"] = (
                float(windows.sum() * SEQ_LEN / n_mamba) if n_mamba else None)
            entry[f"t={t}"] = block
        result["sweep"][name] = entry

    result["stage_fractions"]["fraction_routed_to_pb"] = float(routed_pb.mean())
    result["compute"]["mamba_scoring_seconds"] = {
        name: float(scores[name + "__seconds"][0]) for name in CHECKPOINTS}
    result["compute"]["retained_window_loci"] = int(region["feature_index"].size)
    result["compute"]["total_windows_in_region"] = int(region["n_loci"] / SEQ_LEN)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--region-npz", default=REGION_NPZ)
    parser.add_argument("--out", default=OUT)
    parser.add_argument("--score-cache", default=SCORE_CACHE)
    args = parser.parse_args()

    region = load_region(args.region_npz)
    scores = mamba_scores(region, args.score_cache)
    report = {
        "region_npz": args.region_npz,
        "checkpoints": {k: str(v) for k, v in CHECKPOINTS.items()},
        "leaked_checkpoints": {k: str(v) for k, v in LEAKED_CHECKPOINTS.items()},
        "threshold_sweep": THRESHOLD_SWEEP,
        "stage2": stage2(region, scores),
        "stage3": stage3(region, scores),
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2, default=float))
    logger.info("Wrote %s", args.out)


if __name__ == "__main__":
    main()
