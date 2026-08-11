"""Final comparison for the pileup+LLR experiment: 5 methods, one test region.

Evaluates on cached 15x test counts so every method sees byte-identical input:

A. flat frequency        B. binomial v1        C. binomial v2
D. RawPileupMamba (14ch, existing checkpoint)
E. RawPileupMambaLLR (15ch)  + the shuffled-LLR control

PRIMARY metric is SNP vs Normal. Every threshold is validation-derived and
frozen; test-optimal values appear only under ``oracle_diagnostic``.
Bootstrap CIs are reported because 549 SNPs is small enough that a 0.001 F1
difference is noise.
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
from llr_features import (
    LLRTransform,
    base_features_from_counts,
    raw_llr_from_counts,
)
from model_raw_pileup import RawPileupMamba
from model_raw_pileup_llr import RawPileupMambaLLR

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

SEQ_LEN = 64


def prf(predicted: np.ndarray, truth: np.ndarray) -> dict:
    """Precision/recall/F1 plus confusion counts."""
    tp = int(np.sum(predicted & truth))
    fp = int(np.sum(predicted & ~truth))
    fn = int(np.sum(~predicted & truth))
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)
    return {"precision": precision, "recall": recall, "f1": f1, "tp": tp, "fp": fp, "fn": fn}


def bootstrap_f1(predicted: np.ndarray, truth: np.ndarray, draws: int = 2000,
                 seed: int = 20260811) -> tuple[float, float]:
    """Percentile bootstrap CI for F1, resampling loci with replacement."""
    rng = np.random.default_rng(seed)
    n = truth.size
    scores = np.empty(draws)
    for i in range(draws):
        idx = rng.integers(0, n, n)
        scores[i] = prf(predicted[idx], truth[idx])["f1"]
    return float(np.percentile(scores, 2.5)), float(np.percentile(scores, 97.5))


def model_probability(features: np.ndarray, checkpoint: Path, model_class,
                      device: torch.device) -> tuple[np.ndarray, float, float]:
    """p(SNP) per locus, frozen threshold, and pure inference seconds."""
    state = torch.load(checkpoint, map_location=device, weights_only=False)
    model = model_class().to(device)
    model.load_state_dict(state["model_state_dict"])
    model.eval()                                   # init excluded from the timer below
    windows = torch.tensor(features, dtype=torch.float).reshape(-1, SEQ_LEN,
                                                                features.shape[-1])
    out = []
    if device.type == "cuda":
        torch.cuda.synchronize()
    start = time.perf_counter()
    with torch.no_grad():
        for begin in range(0, windows.shape[0], 256):
            logits = model(pileup_features=windows[begin:begin + 256].to(device))
            out.append(torch.softmax(logits, dim=-1)[..., LABEL_SNP].reshape(-1).cpu())
    if device.type == "cuda":
        torch.cuda.synchronize()
    seconds = time.perf_counter() - start
    return torch.cat(out).numpy(), float(state["best_threshold"]), seconds


def depth_bins(depth: np.ndarray):
    """Standard 15x depth strata."""
    for low, high in [(1, 4), (5, 9), (10, 14), (15, 19), (20, 29), (30, 10 ** 9)]:
        label = f"{low}-{high}" if high < 10 ** 9 else f"{low}+"
        yield label, (depth >= low) & (depth <= high)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-counts", required=True)
    parser.add_argument("--test-counts", required=True)
    parser.add_argument("--freq-threshold", type=float, required=True)
    parser.add_argument("--binom-v1-threshold", type=float, required=True)
    parser.add_argument("--binom-v2-threshold", type=float, required=True)
    parser.add_argument("--old-checkpoint", required=True)
    parser.add_argument("--llr-checkpoint", required=True)
    parser.add_argument("--shuffled-checkpoint", default="")
    parser.add_argument("--out", default="results/llr_experiment_results.json")
    args = parser.parse_args()
    return args


def main() -> None:
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    train_blob = np.load(args.train_counts, allow_pickle=True)
    test_blob = np.load(args.test_counts, allow_pickle=True)
    counts, labels = test_blob["counts"], test_blob["labels"]
    snp = labels == LABEL_SNP
    primary_mask = (labels == 0) | snp
    depth = counts[:, 6]

    caller_v1 = BinomialVariantCaller(include_homozygous=True)
    caller_v2 = BinomialVariantCaller(quality_derived_epsilon=True, include_homozygous=True)

    start = time.perf_counter()
    llr_v1 = raw_llr_from_counts(counts, caller_v1)
    binomial_seconds = time.perf_counter() - start
    llr_v2 = raw_llr_from_counts(counts, caller_v2)

    base = base_features_from_counts(counts)
    vaf = np.where(counts[:, 0:4].sum(1) > 0,
                   np.max(np.where(
                       np.arange(4)[None, :] == (counts[:, 9].astype(int) - 1)[:, None],
                       -1.0, counts[:, 0:4]), axis=1) / np.maximum(counts[:, 0:4].sum(1), 1), 0.0)

    # LLR feature channel must be rebuilt with the SAME transform the model
    # trained under, which is stored in its checkpoint.
    llr_state = torch.load(args.llr_checkpoint, map_location="cpu", weights_only=False)
    stored = llr_state["llr_transform"]
    transform = LLRTransform(mean=stored["mean"], scale=stored["scale"])
    use_v2 = bool(llr_state.get("quality_epsilon", False))
    llr_channel = transform.apply(llr_v2 if use_v2 else llr_v1)
    features_llr = np.concatenate([base, llr_channel[:, None]], axis=1)

    old_probability, old_threshold, old_seconds = model_probability(
        base, Path(args.old_checkpoint), RawPileupMamba, device)
    llr_probability, llr_threshold, llr_seconds = model_probability(
        features_llr, Path(args.llr_checkpoint), RawPileupMambaLLR, device)

    methods = {
        "frequency": (vaf, args.freq_threshold),
        "binomial_v1": (llr_v1, args.binom_v1_threshold),
        "binomial_v2": (llr_v2, args.binom_v2_threshold),
        "mamba_14ch": (old_probability, old_threshold),
        "mamba_llr_15ch": (llr_probability, llr_threshold),
    }
    if args.shuffled_checkpoint:
        shuffled_state = torch.load(args.shuffled_checkpoint, map_location="cpu",
                                    weights_only=False)
        shuffled_transform = LLRTransform(mean=shuffled_state["llr_transform"]["mean"],
                                          scale=shuffled_state["llr_transform"]["scale"])
        rng = np.random.default_rng(777)
        shuffled_channel = rng.permutation(shuffled_transform.apply(llr_v1))
        features_shuffled = np.concatenate([base, shuffled_channel[:, None]], axis=1)
        shuffled_probability, shuffled_threshold, _ = model_probability(
            features_shuffled, Path(args.shuffled_checkpoint), RawPileupMambaLLR, device)
        methods["mamba_shuffled_llr"] = (shuffled_probability, shuffled_threshold)

    report: dict = {"primary": {}, "secondary": {}, "by_depth": {}, "oracle_diagnostic": {},
                    "thresholds": {name: float(t) for name, (_, t) in methods.items()}}

    for frame, mask in (("primary", primary_mask), ("secondary", np.ones_like(snp))):
        for name, (score, threshold) in methods.items():
            stats = prf(score[mask] >= threshold, snp[mask])
            if frame == "primary":
                low, high = bootstrap_f1(score[mask] >= threshold, snp[mask])
                stats["f1_ci95"] = [low, high]
            report[frame][name] = stats
        report[frame]["_loci"] = int(mask.sum())
        report[frame]["_snp"] = int(snp[mask].sum())

    for name, (score, _) in methods.items():
        grid = (np.arange(0.0, 1.0, 0.005) if name.startswith(("frequency", "mamba"))
                else np.arange(-20.0, 200.0, 0.5))
        best = {"f1": -1.0}
        best_threshold = 0.0
        for candidate in grid:
            stats = prf(score[primary_mask] >= candidate, snp[primary_mask])
            if stats["f1"] > best["f1"]:
                best, best_threshold = stats, float(candidate)
        report["oracle_diagnostic"][name] = {"threshold": best_threshold, **best}

    for label, bin_mask in depth_bins(depth):
        combined = bin_mask & primary_mask
        entry = {"loci": int(combined.sum()), "snp": int(snp[combined].sum())}
        for name, (score, threshold) in methods.items():
            entry[name] = prf(score[combined] >= threshold, snp[combined])
        report["by_depth"][label] = entry

    # --- is the model just a monotone transform of the LLR? ---
    order_llr = np.argsort(llr_v1[primary_mask])
    order_model = np.argsort(llr_probability[primary_mask])
    rank_llr = np.empty_like(order_llr); rank_llr[order_llr] = np.arange(order_llr.size)
    rank_model = np.empty_like(order_model); rank_model[order_model] = np.arange(order_model.size)
    spearman = float(np.corrcoef(rank_llr, rank_model)[0, 1])
    logit = np.log(np.clip(llr_probability[primary_mask], 1e-9, 1 - 1e-9) /
                   np.clip(1 - llr_probability[primary_mask], 1e-9, 1 - 1e-9))
    pearson = float(np.corrcoef(llr_v1[primary_mask], logit)[0, 1])
    report["llr_vs_model"] = {"spearman_rank_corr": spearman, "pearson_llr_vs_logit": pearson}

    # --- error analysis: does the model recover binomial false negatives? ---
    binomial_call = llr_v1 >= args.binom_v1_threshold
    model_call = llr_probability >= llr_threshold
    binomial_fn = primary_mask & snp & ~binomial_call
    model_fn = primary_mask & snp & ~model_call
    report["error_analysis"] = {
        "binomial_fn": int(binomial_fn.sum()),
        "model_fn": int(model_fn.sum()),
        "binomial_fn_recovered_by_model": int((binomial_fn & model_call).sum()),
        "model_fn_recovered_by_binomial": int((model_fn & binomial_call).sum()),
        "shared_fn": int((binomial_fn & model_fn).sum()),
        "binomial_fp": int((primary_mask & ~snp & binomial_call).sum()),
        "model_fp": int((primary_mask & ~snp & model_call).sum()),
        "model_new_fp_not_binomial": int((primary_mask & ~snp & model_call & ~binomial_call).sum()),
    }
    for tag, mask in (("binomial_fn", binomial_fn), ("model_fn", model_fn)):
        if int(mask.sum()):
            report["error_analysis"][f"{tag}_profile"] = {
                "depth_median": float(np.median(depth[mask])),
                "depth_min": float(depth[mask].min()), "depth_max": float(depth[mask].max()),
                "vaf_median": float(np.median(vaf[mask])),
                "llr_median": float(np.median(llr_v1[mask])),
                "mean_bq_median": float(np.median(base[mask, 7] * 40)),
                "mean_mq_median": float(np.median(base[mask, 8] * 60)),
            }

    report["runtime"] = {
        "binomial_scoring_seconds": binomial_seconds,
        "binomial_loci_per_second": float(labels.size / max(binomial_seconds, 1e-9)),
        "mamba_14ch_inference_seconds": old_seconds,
        "mamba_llr_inference_seconds": llr_seconds,
        "mamba_llr_loci_per_second": float(labels.size / max(llr_seconds, 1e-9)),
        "device": str(device),
        "note": "model init/checkpoint load excluded from inference timings",
    }
    report["dataset"] = {"loci": int(labels.size), "snp": int(snp.sum()),
                         "normal": int((labels == 0).sum()),
                         "insertion": int((labels == 2).sum()),
                         "deletion": int((labels == 3).sum()),
                         "depth_median": float(np.median(depth))}

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))

    print("\n=== PRIMARY: SNP vs Normal (frozen validation thresholds) ===")
    for name in methods:
        s = report["primary"][name]
        print(f"  {name:<20} P={s['precision']:.4f} R={s['recall']:.4f} F1={s['f1']:.4f} "
              f"[{s['f1_ci95'][0]:.4f},{s['f1_ci95'][1]:.4f}] "
              f"TP={s['tp']} FP={s['fp']} FN={s['fn']}")
    print("\n=== SECONDARY: SNP vs Rest ===")
    for name in methods:
        s = report["secondary"][name]
        print(f"  {name:<20} P={s['precision']:.4f} R={s['recall']:.4f} F1={s['f1']:.4f} "
              f"TP={s['tp']} FP={s['fp']} FN={s['fn']}")
    print(f"\nLLR vs model: spearman={spearman:.4f} pearson(LLR,logit)={pearson:.4f}")
    print("error analysis:", json.dumps(report["error_analysis"], indent=2))
    logger.info("Wrote %s", args.out)


if __name__ == "__main__":
    main()
