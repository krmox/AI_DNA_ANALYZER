"""Final comparison for the residual-binomial experiment on the 15x test loci.

Five methods, one cached count matrix, so every method sees byte-identical
input:

A. frequency (flat VAF rule)      B. binomial v1
C. RawPileupMamba 14ch (existing checkpoint, untouched)
D. RawPileupMambaLLR 15ch (existing checkpoint, untouched)
E. Residual Mamba  (+ the shuffled-prior control)

Scores live on different axes -- VAF in [0,1], LLR in logits, the two older
Mamba models in softmax probability -- so each method carries its own frozen,
validation-derived threshold. No threshold is ever chosen on test; the
test-optimal values appear only under ``oracle_diagnostic`` and are labelled
as an unreachable upper bound.

The residual model is reported in two selection modes, both validation-only:

* ``residual_selected`` -- the checkpoint validation actually prefers,
  *including* the untrained epoch 0. This is the honest deployment answer, and
  if it lands on epoch 0 the honest conclusion is "training the residual did
  not help";
* ``residual_best_trained`` -- the best epoch among 1..N, i.e. what happens if
  training is forced to move. This is what makes the diagnostics meaningful,
  because a residual pinned at exactly zero has nothing to analyse.

The diagnostics section is the point of the experiment: where the residual
moves, whether it moves *usefully*, and the four numbers that decide it --
binomial FNs recovered, binomial FPs corrected, new FPs introduced, binomial
TPs overturned.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import torch

from binomial_baseline import BinomialVariantCaller
from config import LABEL_SNP
from llr_features import LLRTransform, base_features_from_counts, raw_llr_from_counts
from model_raw_pileup import RawPileupMamba
from model_raw_pileup_llr import RawPileupMambaLLR
from model_residual_binomial import ResidualBinomialMamba
from residual_metrics import (
    bootstrap_f1,
    depth_bin_masks,
    paired_bootstrap_f1_delta,
    prf,
    ranking_metrics,
    vaf_from_counts,
)
from residual_prior import NUM_CLASSES, prior_logits_from_llr, snp_decision_score

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

SEQ_LEN = 64


@torch.no_grad()
def probability_model_score(features: np.ndarray, checkpoint: Path, model_class,
                            device: torch.device) -> tuple[np.ndarray, float]:
    """p(SNP) and frozen threshold for the two older softmax models."""
    state = torch.load(checkpoint, map_location=device, weights_only=False)
    model = model_class().to(device)
    model.load_state_dict(state["model_state_dict"])
    model.eval()
    windows = torch.tensor(features, dtype=torch.float).reshape(-1, SEQ_LEN,
                                                                features.shape[-1])
    out = []
    for begin in range(0, windows.shape[0], 256):
        logits = model(pileup_features=windows[begin:begin + 256].to(device))
        out.append(torch.softmax(logits, dim=-1)[..., LABEL_SNP].reshape(-1).cpu())
    return torch.cat(out).numpy().astype(np.float64), float(state["best_threshold"])


@torch.no_grad()
def residual_model_score(features: np.ndarray, llr: np.ndarray, checkpoint: Path,
                         device: torch.device) -> tuple[np.ndarray, np.ndarray, float]:
    """Decision score, per-locus SNP-channel residual, and frozen threshold.

    Args:
        features: ``[N, 14]`` raw pileup features.
        llr: ``[N]`` binomial LLRs forming the prior (already permuted for the
            shuffled control).
        checkpoint: Checkpoint to load.
        device: Inference device.

    Returns:
        ``(score, snp_residual, frozen_threshold)`` where ``score`` is the
        SNP-minus-Normal logit gap and ``snp_residual`` is that gap minus the
        prior LLR, i.e. exactly the correction the model applied.
    """
    state = torch.load(checkpoint, map_location=device, weights_only=False)
    model = ResidualBinomialMamba().to(device)
    model.load_state_dict(state["model_state_dict"])
    model.eval()
    windows = torch.tensor(features, dtype=torch.float).reshape(-1, SEQ_LEN,
                                                                features.shape[-1])
    prior = torch.tensor(prior_logits_from_llr(llr), dtype=torch.float).reshape(
        -1, SEQ_LEN, NUM_CLASSES)
    scores, residuals = [], []
    for begin in range(0, windows.shape[0], 256):
        logits, residual = model(pileup_features=windows[begin:begin + 256].to(device),
                                 prior_logits=prior[begin:begin + 256].to(device),
                                 return_residual=True)
        scores.append(snp_decision_score(logits).reshape(-1).cpu())
        residuals.append(snp_decision_score(residual).reshape(-1).cpu())
    return (torch.cat(scores).numpy().astype(np.float64),
            torch.cat(residuals).numpy().astype(np.float64),
            float(state["best_threshold"]))


def bin_profile(values: np.ndarray, by: np.ndarray, edges) -> dict:
    """Mean/median of ``values`` within bins of ``by``."""
    out = {}
    for low, high in edges:
        mask = (by >= low) & (by < high)
        label = f"{low}-{high}"
        out[label] = ({"n": int(mask.sum()),
                       "mean": float(values[mask].mean()),
                       "median": float(np.median(values[mask]))}
                      if mask.any() else {"n": 0})
    return out


def describe(values: np.ndarray) -> dict:
    """Distribution summary used throughout the diagnostics."""
    if values.size == 0:
        return {"n": 0}
    return {"n": int(values.size), "mean": float(values.mean()),
            "std": float(values.std()), "min": float(values.min()),
            "p01": float(np.percentile(values, 1)), "p25": float(np.percentile(values, 25)),
            "median": float(np.median(values)), "p75": float(np.percentile(values, 75)),
            "p99": float(np.percentile(values, 99)), "max": float(values.max()),
            "fraction_nonzero": float(np.mean(np.abs(values) > 1e-6)),
            "mean_abs": float(np.abs(values).mean())}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-counts", required=True)
    parser.add_argument("--freq-threshold", type=float, required=True)
    parser.add_argument("--binom-v1-threshold", type=float, required=True)
    parser.add_argument("--old-checkpoint", required=True)
    parser.add_argument("--llr-checkpoint", required=True)
    parser.add_argument("--residual-selected-checkpoint", required=True)
    parser.add_argument("--residual-best-trained-checkpoint", required=True)
    parser.add_argument("--shuffled-checkpoint", default="")
    parser.add_argument("--shuffle-seed", type=int, default=777)
    parser.add_argument("--out", default="results/residual_experiment_results.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    blob = np.load(args.test_counts, allow_pickle=True)
    counts, labels = blob["counts"], blob["labels"]
    snp = labels == LABEL_SNP
    primary = (labels == 0) | snp
    depth = counts[:, 6]
    vaf = vaf_from_counts(counts)
    base = base_features_from_counts(counts)

    caller = BinomialVariantCaller()
    llr = raw_llr_from_counts(counts, caller)

    old_score, old_threshold = probability_model_score(
        base, Path(args.old_checkpoint), RawPileupMamba, device)

    llr_state = torch.load(args.llr_checkpoint, map_location="cpu", weights_only=False)
    transform = LLRTransform(mean=llr_state["llr_transform"]["mean"],
                             scale=llr_state["llr_transform"]["scale"])
    features_llr = np.concatenate([base, transform.apply(llr)[:, None]], axis=1)
    llr_model_score, llr_model_threshold = probability_model_score(
        features_llr, Path(args.llr_checkpoint), RawPileupMambaLLR, device)

    selected_score, selected_residual, selected_threshold = residual_model_score(
        base, llr, Path(args.residual_selected_checkpoint), device)
    trained_score, trained_residual, trained_threshold = residual_model_score(
        base, llr, Path(args.residual_best_trained_checkpoint), device)

    methods = {
        "frequency": (vaf, args.freq_threshold),
        "binomial_v1": (llr, args.binom_v1_threshold),
        "mamba_14ch": (old_score, old_threshold),
        "mamba_llr_15ch": (llr_model_score, llr_model_threshold),
        "residual_selected": (selected_score, selected_threshold),
        "residual_best_trained": (trained_score, trained_threshold),
    }

    if args.shuffled_checkpoint:
        # The control model was trained against a permuted prior, so it must be
        # evaluated against the same permutation of the test LLRs.
        shuffled_llr = np.random.default_rng(args.shuffle_seed).permutation(llr)
        shuffled_score, shuffled_residual, shuffled_threshold = residual_model_score(
            base, shuffled_llr, Path(args.shuffled_checkpoint), device)
        methods["residual_shuffled_prior"] = (shuffled_score, shuffled_threshold)

    report: dict = {
        "test_counts": args.test_counts,
        "thresholds": {name: float(t) for name, (_, t) in methods.items()},
        "threshold_provenance": "all frozen from validation; test never used for selection",
        "primary": {}, "secondary": {}, "by_depth": {}, "oracle_diagnostic": {},
        "checkpoints": {
            "mamba_14ch": args.old_checkpoint, "mamba_llr_15ch": args.llr_checkpoint,
            "residual_selected": args.residual_selected_checkpoint,
            "residual_best_trained": args.residual_best_trained_checkpoint,
            "residual_shuffled_prior": args.shuffled_checkpoint,
        },
    }

    # --- primary: SNP vs Normal; secondary: SNP vs everything else ---
    for frame, mask in (("primary", primary), ("secondary", np.ones_like(snp, dtype=bool))):
        for name, (score, threshold) in methods.items():
            stats = prf(score[mask] >= threshold, snp[mask])
            stats.update(ranking_metrics(score[mask], snp[mask]))
            if frame == "primary":
                stats["f1_ci95"] = bootstrap_f1(score[mask] >= threshold, snp[mask])
            report[frame][name] = stats
        report[frame]["_loci"] = int(mask.sum())
        report[frame]["_snp"] = int(snp[mask].sum())

    # --- paired comparison against the binomial baseline ---
    binomial_call = llr >= args.binom_v1_threshold
    report["vs_binomial_paired_bootstrap"] = {
        name: paired_bootstrap_f1_delta(score[primary] >= threshold,
                                        binomial_call[primary], snp[primary])
        for name, (score, threshold) in methods.items() if name != "binomial_v1"
    }

    # --- oracle upper bound (diagnostic only, never used for any decision) ---
    for name, (score, _) in methods.items():
        grid = (np.arange(0.0, 1.0, 0.005) if name in ("frequency", "mamba_14ch",
                                                       "mamba_llr_15ch")
                else np.arange(-20.0, 200.0, 0.5))
        best, best_threshold = {"f1": -1.0}, 0.0
        for candidate in grid:
            stats = prf(score[primary] >= candidate, snp[primary])
            if stats["f1"] > best["f1"]:
                best, best_threshold = stats, float(candidate)
        report["oracle_diagnostic"][name] = {"threshold": best_threshold, **best}

    # --- depth stratification ---
    for label, bin_mask in depth_bin_masks(depth):
        combined = bin_mask & primary
        entry = {"loci": int(combined.sum()), "snp": int(snp[combined].sum())}
        for name, (score, threshold) in methods.items():
            entry[name] = prf(score[combined] >= threshold, snp[combined])
        report["by_depth"][label] = entry

    # --- THE DIAGNOSTIC: what did the residual actually do? ---
    diagnostics = {}
    for name, residual, score, threshold in (
        ("residual_selected", selected_residual, selected_score, selected_threshold),
        ("residual_best_trained", trained_residual, trained_score, trained_threshold),
    ):
        model_call = score >= threshold
        binomial_tp = primary & snp & binomial_call
        binomial_fn = primary & snp & ~binomial_call
        binomial_fp = primary & ~snp & binomial_call
        binomial_tn = primary & ~snp & ~binomial_call

        entry = {
            "distribution_all_loci": describe(residual[primary]),
            "distribution_snp": describe(residual[primary & snp]),
            "distribution_normal": describe(residual[primary & ~snp]),
            "by_depth": bin_profile(residual[primary], depth[primary],
                                    [(1, 5), (5, 10), (10, 15), (15, 20), (20, 30),
                                     (30, 10 ** 9)]),
            "by_vaf": bin_profile(residual[primary], vaf[primary],
                                  [(0.0, 0.05), (0.05, 0.15), (0.15, 0.3), (0.3, 0.45),
                                   (0.45, 0.6), (0.6, 0.85), (0.85, 1.01)]),
            "by_llr": bin_profile(residual[primary], llr[primary],
                                  [(-100, -10), (-10, 0), (0, 5), (5, 10.5), (10.5, 25),
                                   (25, 60), (60, 10 ** 9)]),
            "residual_on_binomial_fn": describe(residual[binomial_fn]),
            "residual_on_binomial_fp": describe(residual[binomial_fp]),
            "residual_on_binomial_tp": describe(residual[binomial_tp]),
            # The four numbers the experiment turns on.
            "binomial_fn_recovered": int((binomial_fn & model_call).sum()),
            "binomial_fn_total": int(binomial_fn.sum()),
            "binomial_fp_corrected": int((binomial_fp & ~model_call).sum()),
            "binomial_fp_total": int(binomial_fp.sum()),
            "new_fp_introduced": int((binomial_tn & model_call).sum()),
            "binomial_tp_overturned": int((binomial_tp & ~model_call).sum()),
            "net_call_changes": int((model_call[primary] != binomial_call[primary]).sum()),
        }
        entry["net_true_positive_gain"] = entry["binomial_fn_recovered"] - \
            entry["binomial_tp_overturned"]
        entry["net_false_positive_gain"] = entry["binomial_fp_corrected"] - \
            entry["new_fp_introduced"]
        diagnostics[name] = entry
    report["residual_diagnostics"] = diagnostics

    # --- is the residual anything more than a constant shift? ---
    # A residual that is a near-constant ``c`` is not new information: it just
    # slides the decision threshold to ``threshold - c``, which the binomial
    # caller can do on its own by turning its one knob. Three tests:
    #   1. the residual's own discriminative power (AUC as a standalone score);
    #   2. the constant-shift null, i.e. the binomial evaluated at the exactly
    #      equivalent shifted threshold;
    #   3. the plain binomial swept over nearby thresholds, to see whether the
    #      residual model's operating point is even the best one available to
    #      the baseline without any neural network at all.
    shift_analysis = {}
    for name, residual, score, threshold in (
        ("residual_selected", selected_residual, selected_score, selected_threshold),
        ("residual_best_trained", trained_residual, trained_score, trained_threshold),
    ):
        constant = float(residual.mean())
        equivalent = threshold - constant
        shift_analysis[name] = {
            "mean_residual": constant,
            "std_residual": float(residual.std()),
            "residual_alone_as_snp_score": ranking_metrics(residual[primary], snp[primary]),
            "corr_residual_llr": float(np.corrcoef(residual[primary], llr[primary])[0, 1])
            if residual.std() > 0 else 0.0,
            "corr_residual_depth": float(np.corrcoef(residual[primary], depth[primary])[0, 1])
            if residual.std() > 0 else 0.0,
            "corr_residual_vaf": float(np.corrcoef(residual[primary], vaf[primary])[0, 1])
            if residual.std() > 0 else 0.0,
            "model_at_frozen_threshold": prf(score[primary] >= threshold, snp[primary]),
            "constant_shift_null": {
                "equivalent_binomial_threshold": equivalent,
                **prf(llr[primary] >= equivalent, snp[primary]),
            },
            "plain_binomial_threshold_scan": {
                f"{t:.2f}": prf(llr[primary] >= t, snp[primary])
                for t in (equivalent - 0.5, equivalent - 0.25, equivalent,
                          threshold, args.binom_v1_threshold)
            },
        }
    report["constant_shift_analysis"] = shift_analysis

    report["dataset"] = {
        "loci": int(labels.size), "snp": int(snp.sum()), "normal": int((labels == 0).sum()),
        "insertion": int((labels == 2).sum()), "deletion": int((labels == 3).sum()),
        "depth_median": float(np.median(depth)),
    }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))

    print("\n=== PRIMARY: SNP vs Normal (frozen validation thresholds) ===")
    for name in methods:
        s = report["primary"][name]
        print(f"  {name:<24} P={s['precision']:.4f} R={s['recall']:.4f} F1={s['f1']:.4f} "
              f"[{s['f1_ci95'][0]:.4f},{s['f1_ci95'][1]:.4f}] "
              f"AUC={s['roc_auc']:.5f} AP={s['pr_auc']:.4f} "
              f"TP={s['tp']} FP={s['fp']} FN={s['fn']}")
    print("\n=== PAIRED delta F1 vs binomial v1 (95% CI) ===")
    for name, d in report["vs_binomial_paired_bootstrap"].items():
        print(f"  {name:<24} {d['delta_f1']:+.4f} [{d['ci95'][0]:+.4f},{d['ci95'][1]:+.4f}] "
              f"P(better)={d['fraction_of_draws_favouring_a']:.3f}")
    print("\n=== RESIDUAL DIAGNOSTICS ===")
    for name, d in diagnostics.items():
        print(f"  {name}: residual mean={d['distribution_all_loci']['mean']:+.4f} "
              f"sd={d['distribution_all_loci']['std']:.4f} | "
              f"FN recovered {d['binomial_fn_recovered']}/{d['binomial_fn_total']} | "
              f"FP corrected {d['binomial_fp_corrected']}/{d['binomial_fp_total']} | "
              f"new FP {d['new_fp_introduced']} | TP overturned {d['binomial_tp_overturned']}")
    print("\n=== IS THE RESIDUAL MORE THAN A CONSTANT SHIFT? ===")
    for name, s in shift_analysis.items():
        print(f"  {name}: mean={s['mean_residual']:+.4f} sd={s['std_residual']:.4f} | "
              f"residual-alone AUC={s['residual_alone_as_snp_score']['roc_auc']:.4f} "
              f"AP={s['residual_alone_as_snp_score']['pr_auc']:.5f} | "
              f"corr(llr)={s['corr_residual_llr']:+.4f}")
        print(f"      model F1={s['model_at_frozen_threshold']['f1']:.4f}  vs  "
              f"constant-shift null F1={s['constant_shift_null']['f1']:.4f} "
              f"(binomial @ {s['constant_shift_null']['equivalent_binomial_threshold']:.3f})")
    logger.info("Wrote %s", args.out)


if __name__ == "__main__":
    main()
