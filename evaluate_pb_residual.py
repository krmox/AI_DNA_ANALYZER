"""Final evaluation of the PB-residual model against every established arm.

Comparison set, as required by the ablation principle:

1. Poisson-binomial (the current strong statistical baseline)
2. Binomial v1 (the previous statistical baseline)
3. RawPileupMamba, 14 channels (the previous neural model)
4. Residual binomial Mamba (the previous residual model)
5. The new PB-residual model, plus its controls and ablations

Every arm is scored on the same loci with thresholds frozen on validation. The
neural thresholds come from each checkpoint's own validation selection, made
before this script runs.

Runs only on checkpoints that were already selected on validation. This script
does not choose an architecture; it reports one.
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import torch

from binomial_baseline import BinomialVariantCaller
from config import LABEL_SNP
from evaluate_binomial_baseline import depth_bins, pick_threshold, prf
from evaluate_quality_error import LLR_GRID, paired_bootstrap
from locus_evidence import FEATURE_DIM, resolve_group_indices
from model_pb_residual import PoissonBinomialResidualMamba
from residual_prior import prior_logits_from_llr, snp_decision_score
from train_raw_pileup import SEQ_LEN

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


@torch.no_grad()
def model_scores(checkpoint: Path, features: np.ndarray, llr: np.ndarray,
                 device: torch.device, batch_size: int = 512) -> tuple[np.ndarray, dict, float]:
    """Decision scores from a trained PB-residual checkpoint.

    Args:
        checkpoint: Path to the checkpoint, opened read-only.
        features: ``[N, FEATURE_DIM]`` locus evidence.
        llr: ``[N]`` Poisson-binomial LLR.
        device: Inference device.
        batch_size: Windows per forward pass.

    Returns:
        ``(score, residual_diagnostics, seconds)``.
    """
    state = torch.load(checkpoint, map_location=device, weights_only=False)
    model = PoissonBinomialResidualMamba(
        feature_dim=FEATURE_DIM, use_gate=bool(state.get("use_gate", False))).to(device)
    model.load_state_dict(state["model_state_dict"])
    model.eval()

    prepared = np.array(features, dtype=np.float32, copy=True)
    for index in resolve_group_indices(tuple(state.get("drop_groups", []) or [])):
        prepared[:, index] = 0.0

    feature_windows = torch.tensor(prepared).reshape(-1, SEQ_LEN, FEATURE_DIM)
    prior_windows = torch.tensor(prior_logits_from_llr(llr),
                                 dtype=torch.float).reshape(-1, SEQ_LEN, 4)

    scores, gaps, gates = [], [], []
    start = time.perf_counter()
    for begin in range(0, feature_windows.shape[0], batch_size):
        logits, residual, gate = model(
            feature_windows[begin:begin + batch_size].to(device),
            prior_windows[begin:begin + batch_size].to(device), return_parts=True)
        scores.append(snp_decision_score(logits).reshape(-1).cpu())
        gaps.append(snp_decision_score(residual).reshape(-1).cpu())
        gates.append(gate.reshape(-1).cpu())
    seconds = time.perf_counter() - start

    gap = torch.cat(gaps)
    gate_values = torch.cat(gates)
    diagnostics = {
        "epoch": int(state.get("epoch", -1)),
        "validation_threshold": float(state["best_threshold"]),
        "use_gate": bool(state.get("use_gate", False)),
        "drop_groups": list(state.get("drop_groups", []) or []),
        "residual_gap_mean": float(gap.mean()), "residual_gap_std": float(gap.std()),
        "residual_gap_p1": float(torch.quantile(gap.float(), 0.01)),
        "residual_gap_p99": float(torch.quantile(gap.float(), 0.99)),
        "gate_mean": float(gate_values.mean()), "gate_std": float(gate_values.std()),
    }
    return torch.cat(scores).numpy().astype(np.float64), diagnostics, seconds


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-evidence", default="cache/test_15x_evidence.npz")
    parser.add_argument("--test-counts", default="cache/test_15x_counts.npz")
    parser.add_argument("--checkpoint", action="append", default=[],
                        metavar="NAME=PATH", help="Named checkpoint to evaluate.")
    parser.add_argument("--out", default="results/pb_residual_15x_results.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    blob = np.load(args.test_evidence, allow_pickle=True)
    features, labels, llr = blob["features"], blob["labels"], blob["pb_llr"]
    depth = blob["depth"]
    counts = np.load(args.test_counts)["counts"]

    snp = labels == LABEL_SNP
    frame = (labels == 0) | snp
    report: dict = {"primary": {}, "by_depth": {}, "thresholds": {},
                    "diagnostics": {}, "runtime": {}}

    scores: dict[str, np.ndarray] = {}
    thresholds: dict[str, float] = {}

    # -- statistical baselines, thresholds already frozen in earlier work ----
    binomial_v1 = BinomialVariantCaller(error_rate=0.01).score_counts(
        counts[:, 0:4], counts[:, 9].astype(int))
    scores["binomial_v1"] = binomial_v1.llr
    thresholds["binomial_v1"] = 10.5
    scores["poisson_binomial"] = llr
    thresholds["poisson_binomial"] = 12.5
    for name in ("binomial_v1", "poisson_binomial"):
        report["thresholds"][name] = {
            "threshold": thresholds[name],
            "source": "frozen on validation in the quality-error experiment"}

    # -- trained checkpoints, each carrying its own validation threshold ----
    for entry in args.checkpoint:
        name, _, path = entry.partition("=")
        score, diagnostics, seconds = model_scores(Path(path), features, llr, device)
        scores[name] = score
        thresholds[name] = diagnostics["validation_threshold"]
        report["diagnostics"][name] = diagnostics
        report["thresholds"][name] = {"threshold": diagnostics["validation_threshold"],
                                      "source": f"validation selection inside {path}"}
        report["runtime"][name] = {"inference_seconds": seconds,
                                   "loci_per_second": float(labels.size / max(seconds, 1e-9))}
        logger.info("%-28s epoch=%d thr=%.1f | residual gap %+.3f +- %.3f | gate %.3f",
                    name, diagnostics["epoch"], diagnostics["validation_threshold"],
                    diagnostics["residual_gap_mean"], diagnostics["residual_gap_std"],
                    diagnostics["gate_mean"])

    for name, score in scores.items():
        report["primary"][name] = prf(score[frame] >= thresholds[name], snp[frame])

    report["oracle_diagnostic"] = {
        "_note": "test-optimal thresholds; a diagnostic ceiling, NOT a result"}
    for name, score in scores.items():
        threshold, stats = pick_threshold(score[frame], snp[frame], LLR_GRID)
        report["oracle_diagnostic"][name] = {"threshold": threshold, **stats}

    for label, bin_mask in depth_bins(depth):
        combined = bin_mask & frame
        entry = {"loci": int(combined.sum()), "snp": int(snp[combined].sum())}
        for name, score in scores.items():
            entry[name] = prf(score[combined] >= thresholds[name], snp[combined])
        report["by_depth"][label] = entry

    baseline_calls = scores["poisson_binomial"][frame] >= thresholds["poisson_binomial"]
    report["bootstrap_vs_poisson_binomial"] = {
        name: paired_bootstrap(scores[name][frame] >= thresholds[name],
                               baseline_calls, snp[frame])
        for name in scores if name != "poisson_binomial"}

    report["dataset"] = {"loci": int(labels.size), "frame_loci": int(frame.sum()),
                         "snp": int(snp[frame].sum())}

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2, default=float))
    logger.info("Wrote %s", args.out)

    print("\n=== PRIMARY: SNP vs Normal, validation-frozen thresholds ===")
    for name in scores:
        stats = report["primary"][name]
        print(f"  {name:<28} thr={thresholds[name]:>7.2f} P={stats['precision']:.4f} "
              f"R={stats['recall']:.4f} F1={stats['f1']:.4f} "
              f"TP={stats['tp']} FP={stats['fp']} FN={stats['fn']}")
    print("\n=== Paired bootstrap vs Poisson-binomial (dF1, 95% CI) ===")
    for name, stats in report["bootstrap_vs_poisson_binomial"].items():
        low, high = stats["delta_f1_ci95"]
        print(f"  {name:<28} dF1={stats['delta_f1']:+.4f} [{low:+.4f}, {high:+.4f}] "
              f"dTP={stats['delta_tp']:+d} dFP={stats['delta_fp']:+d} "
              f"dFN={stats['delta_fn']:+d}")


if __name__ == "__main__":
    main()
