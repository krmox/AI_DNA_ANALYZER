"""Benchmark: flat frequency rule vs depth-aware binomial vs RawPileupMamba.

Runs on the 15x real HG002 data. All three methods are scored on exactly the
same loci, with every threshold frozen from train/validation data before the
held-out test region is touched once.

PRIMARY metric is SNP vs Normal -- the variant-detection question. Indel loci
are excluded from the primary because the project has already established that
the neural model can separate SNPs from indels, and letting them count as
negatives conflates that separate ability with detection. SNP vs Rest is
reported as a labelled secondary.
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import torch

from binomial_baseline import (
    BinomialVariantCaller,
    benjamini_hochberg,
    binomial_tail_pvalue,
)
from config import LABEL_SNP
from model_raw_pileup import RawPileupMamba
from pileup_counts import load_counts
from raw_pileup import FEATURE_DIM, RawPileupDataset, RawPileupProvider

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

TRAIN_FRACTION = 0.9  # matches run_train.py / train_raw_pileup.py


def prf(predicted: np.ndarray, truth: np.ndarray) -> dict:
    """Precision/recall/F1 and the confusion counts for a boolean call."""
    tp = int(np.sum(predicted & truth))
    fp = int(np.sum(predicted & ~truth))
    fn = int(np.sum(~predicted & truth))
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)
    return {"precision": precision, "recall": recall, "f1": f1, "tp": tp, "fp": fp, "fn": fn}


def pick_threshold(score: np.ndarray, truth: np.ndarray, grid: np.ndarray) -> tuple[float, dict]:
    """Choose the threshold maximizing F1 on the supplied (non-test) data."""
    best_threshold, best = float(grid[0]), {"f1": -1.0}
    for candidate in grid:
        stats = prf(score >= candidate, truth)
        if stats["f1"] > best["f1"]:
            best_threshold, best = float(candidate), stats
    return best_threshold, best


def model_snp_probability(features: np.ndarray, checkpoint: Path,
                          device: torch.device) -> np.ndarray:
    """p(SNP) per locus from the existing RawPileupMamba checkpoint.

    The checkpoint is loaded read-only and never written back.

    Args:
        features: ``[N, FEATURE_DIM]`` raw-pileup features, N a multiple of 64.
        checkpoint: Path to the trained checkpoint.
        device: Inference device.

    Returns:
        ``[N]`` softmax probability of the SNP class.
    """
    state = torch.load(checkpoint, map_location=device, weights_only=False)
    model = RawPileupMamba().to(device)
    model.load_state_dict(state["model_state_dict"])
    model.eval()
    windows = torch.tensor(features, dtype=torch.float).reshape(-1, 64, FEATURE_DIM)
    out = []
    with torch.no_grad():
        for start in range(0, windows.shape[0], 256):
            logits = model(pileup_features=windows[start:start + 256].to(device))
            out.append(torch.softmax(logits, dim=-1)[..., LABEL_SNP].reshape(-1).cpu())
    return torch.cat(out).numpy()


def load_features(fasta, bam, vcf, bed, region, window_slice=None) -> np.ndarray:
    """Materialize raw-pileup model features for a region (or window slice)."""
    provider = RawPileupProvider(
        fasta_path=fasta, bam_path=bam, vcf_path=vcf, contig="chr21",
        region=region, seq_len=64, high_confidence_bed=bed,
    )
    with provider:
        dataset = RawPileupDataset(provider)
        indices = range(len(dataset))[window_slice] if window_slice else range(len(dataset))
        return np.concatenate([dataset[i]["pileup_features"].numpy() for i in indices])


def depth_bins(depth: np.ndarray) -> list[tuple[str, np.ndarray]]:
    """Standard depth strata for the 15x regime."""
    edges = [(1, 4), (5, 9), (10, 14), (15, 19), (20, 29), (30, 10**9)]
    out = []
    for low, high in edges:
        label = f"{low}-{high}" if high < 10**9 else f"{low}+"
        out.append((label, (depth >= low) & (depth <= high)))
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--train-bam", required=True)
    parser.add_argument("--train-vcf", required=True)
    parser.add_argument("--train-bed", required=True)
    parser.add_argument("--train-region", type=int, nargs=2, required=True)
    parser.add_argument("--test-bam", required=True)
    parser.add_argument("--test-vcf", required=True)
    parser.add_argument("--test-bed", required=True)
    parser.add_argument("--test-region", type=int, nargs=2, required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--out", default="results/binomial_15x_results.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    report: dict = {"thresholds": {}, "primary": {}, "secondary": {}, "by_depth": {},
                    "runtime": {}, "validation": {}}

    # ---------- validation split: same last-10% of windows the model used ----------
    logger.info("Materializing VALIDATION split (last 10%% of the training region)...")
    probe = RawPileupProvider(
        fasta_path=args.fasta, bam_path=args.train_bam, vcf_path=args.train_vcf,
        contig="chr21", region=tuple(args.train_region), seq_len=64,
        high_confidence_bed=args.train_bed,
    )
    with probe:
        total_windows = len(probe)
    split = int(total_windows * TRAIN_FRACTION)
    val_slice = slice(split, total_windows)
    logger.info("Training region has %d windows; validation = windows %d..%d",
                total_windows, split, total_windows - 1)

    val_counts, val_labels = load_counts(
        args.fasta, args.train_bam, args.train_vcf, args.train_bed,
        tuple(args.train_region), window_slice=val_slice)
    val_snp = val_labels == LABEL_SNP
    val_keep = (val_labels == 0) | val_snp          # SNP-vs-Normal, the primary frame
    report["validation"]["loci"] = int(val_labels.size)
    report["validation"]["snp"] = int(val_snp.sum())

    caller_v1 = BinomialVariantCaller(error_rate=0.01, include_homozygous=True)
    caller_v2 = BinomialVariantCaller(quality_derived_epsilon=True, include_homozygous=True)

    val_base = val_counts[:, 0:4]
    val_ref = val_counts[:, 9].astype(int)
    val_qual = val_counts[:, 7]

    val_v1 = caller_v1.score_counts(val_base, val_ref)
    val_v2 = caller_v2.score_counts(val_base, val_ref, quality_sum=val_qual)
    val_vaf = val_v1.vaf

    llr_grid = np.arange(-20.0, 200.0, 0.5)
    vaf_grid = np.arange(0.01, 1.00, 0.005)

    thr_v1, val_stats_v1 = pick_threshold(val_v1.llr[val_keep], val_snp[val_keep], llr_grid)
    thr_v2, val_stats_v2 = pick_threshold(val_v2.llr[val_keep], val_snp[val_keep], llr_grid)
    thr_freq, val_stats_freq = pick_threshold(val_vaf[val_keep], val_snp[val_keep], vaf_grid)

    logger.info("FROZEN thresholds (validation-derived, %d loci / %d SNP):",
                int(val_keep.sum()), int(val_snp[val_keep].sum()))
    for name, threshold, stats in (("frequency", thr_freq, val_stats_freq),
                                   ("binomial_v1", thr_v1, val_stats_v1),
                                   ("binomial_v2_qual", thr_v2, val_stats_v2)):
        logger.info("  %-17s thr=%8.3f  val P=%.4f R=%.4f F1=%.4f",
                    name, threshold, stats["precision"], stats["recall"], stats["f1"])
        report["thresholds"][name] = {"threshold": threshold, "source": "validation split "
                                      "(last 10% of training region), test labels unused",
                                      "validation": stats}

    # ---------- held-out test region, touched once ----------
    logger.info("Materializing TEST region...")
    t0 = time.perf_counter()
    counts, labels = load_counts(args.fasta, args.test_bam, args.test_vcf, args.test_bed,
                                 tuple(args.test_region))
    extract_seconds = time.perf_counter() - t0

    base = counts[:, 0:4]
    reference_index = counts[:, 9].astype(int)
    quality_sum = counts[:, 7]
    depth = counts[:, 6]
    snp = labels == LABEL_SNP
    normal_or_snp = (labels == 0) | snp

    t0 = time.perf_counter()
    score_v1 = caller_v1.score_counts(base, reference_index)
    binomial_seconds = time.perf_counter() - t0
    score_v2 = caller_v2.score_counts(base, reference_index, quality_sum=quality_sum)
    vaf = score_v1.vaf

    logger.info("Running RawPileupMamba inference on the same loci...")
    features = load_features(args.fasta, args.test_bam, args.test_vcf, args.test_bed,
                             tuple(args.test_region))
    assert features.shape[0] == labels.size, "feature/count locus mismatch"
    t0 = time.perf_counter()
    model_probability = model_snp_probability(features, Path(args.checkpoint), device)
    model_seconds = time.perf_counter() - t0

    model_threshold = float(torch.load(args.checkpoint, map_location="cpu",
                                       weights_only=False)["best_threshold"])
    report["thresholds"]["rawpileup_mamba"] = {
        "threshold": model_threshold,
        "source": "validation split of the 15x training region (frozen at training time)",
    }

    methods = {
        "frequency": (vaf, thr_freq),
        "binomial_v1": (score_v1.llr, thr_v1),
        "binomial_v2_qual": (score_v2.llr, thr_v2),
        "rawpileup_mamba": (model_probability, model_threshold),
    }

    for frame, mask in (("primary", normal_or_snp), ("secondary", np.ones_like(snp))):
        for name, (score, threshold) in methods.items():
            report[frame][name] = prf(score[mask] >= threshold, snp[mask])
        report[frame]["_loci"] = int(mask.sum())
        report[frame]["_snp"] = int(snp[mask].sum())

    # Oracle diagnostics, clearly separated from the primary numbers.
    report["oracle_diagnostic"] = {}
    for name, (score, _) in methods.items():
        grid = vaf_grid if name == "frequency" else (
            np.arange(0.0, 1.0, 0.005) if name == "rawpileup_mamba" else llr_grid)
        oracle_threshold, stats = pick_threshold(score[normal_or_snp], snp[normal_or_snp], grid)
        report["oracle_diagnostic"][name] = {"threshold": oracle_threshold, **stats}

    # Multiple-testing diagnostic on the primary frame.
    pvalues = binomial_tail_pvalue(score_v1.k, score_v1.n, score_v1.epsilon)
    rejected = benjamini_hochberg(pvalues[normal_or_snp], alpha=0.05)
    report["multiple_testing_diagnostic"] = {
        "method": "binomial upper-tail p-value under H0 + Benjamini-Hochberg, alpha=0.05",
        **prf(rejected, snp[normal_or_snp]),
    }

    # ---------- depth-stratified, primary frame ----------
    for label, bin_mask in depth_bins(depth):
        combined = bin_mask & normal_or_snp
        entry = {"loci": int(combined.sum()), "snp": int(snp[combined].sum())}
        for name, (score, threshold) in methods.items():
            entry[name] = prf(score[combined] >= threshold, snp[combined])
        report["by_depth"][label] = entry

    report["runtime"] = {
        "count_extraction_seconds": extract_seconds,
        "binomial_scoring_seconds": binomial_seconds,
        "binomial_loci_per_second": float(labels.size / max(binomial_seconds, 1e-9)),
        "model_inference_seconds": model_seconds,
        "model_loci_per_second": float(labels.size / max(model_seconds, 1e-9)),
        "device": str(device),
    }
    report["dataset"] = {
        "test_region": args.test_region, "loci": int(labels.size),
        "snp": int(snp.sum()), "normal": int((labels == 0).sum()),
        "insertion": int((labels == 2).sum()), "deletion": int((labels == 3).sum()),
        "depth_median": float(np.median(depth)), "depth_mean": float(depth.mean()),
    }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    logger.info("Wrote %s", args.out)

    print("\n=== PRIMARY: SNP vs Normal ===")
    for name in methods:
        s = report["primary"][name]
        print(f"  {name:<18} P={s['precision']:.4f} R={s['recall']:.4f} F1={s['f1']:.4f} "
              f"TP={s['tp']} FP={s['fp']} FN={s['fn']}")
    print("\n=== SECONDARY: SNP vs Rest ===")
    for name in methods:
        s = report["secondary"][name]
        print(f"  {name:<18} P={s['precision']:.4f} R={s['recall']:.4f} F1={s['f1']:.4f} "
              f"TP={s['tp']} FP={s['fp']} FN={s['fn']}")


if __name__ == "__main__":
    main()
