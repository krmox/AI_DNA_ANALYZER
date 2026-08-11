"""Final comparison for the read-level experiment on the 15x test loci.

Five methods, one cached region, byte-identical inputs:

A. flat frequency      B. binomial v1      C. binomial v2 (quality-derived eps)
D. RawPileupMamba 14ch (existing checkpoint, untouched)
E. ReadLevelMamba      (+ any controls/ablations passed on the command line)

Every threshold is frozen from validation. Test-optimal values appear only
under ``oracle_diagnostic`` and are never used as a primary result.

The error analysis is the substantive part: for the variants the read-level
model recovers and the binomial misses (and vice versa), it profiles the
*read-level* evidence -- alternate-allele base qualities, strand balance,
read-position distribution -- to test the actual scientific hypothesis, which
is that alternate observations differ in ways an aggregate ``(n, k)`` statistic
cannot represent.
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
from llr_features import base_features_from_counts, raw_llr_from_counts
from model_raw_pileup import RawPileupMamba
from model_read_level import ReadLevelMamba
from read_level_controls import apply_control
from read_level_pileup import (
    FLAG_COUNTED,
    FLAG_NEAR_END,
    FLAG_VALID,
    build_features_torch,
    resolve_drop_indices,
)
from residual_metrics import (
    bootstrap_f1,
    depth_bin_masks,
    paired_bootstrap_f1_delta,
    prf,
    ranking_metrics,
    vaf_from_counts,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

SEQ_LEN = 64


@torch.no_grad()
def aggregate_model_score(features: np.ndarray, checkpoint: Path,
                          device: torch.device) -> tuple[np.ndarray, float, float]:
    """p(SNP), frozen threshold and inference seconds for the 14-channel model."""
    state = torch.load(checkpoint, map_location=device, weights_only=False)
    model = RawPileupMamba().to(device)
    model.load_state_dict(state["model_state_dict"])
    model.eval()
    windows = torch.tensor(features, dtype=torch.float).reshape(-1, SEQ_LEN,
                                                                features.shape[-1])
    out = []
    if device.type == "cuda":
        torch.cuda.synchronize()
    start = time.perf_counter()
    for begin in range(0, windows.shape[0], 256):
        logits = model(pileup_features=windows[begin:begin + 256].to(device))
        out.append(torch.softmax(logits, dim=-1)[..., LABEL_SNP].reshape(-1).cpu())
    if device.type == "cuda":
        torch.cuda.synchronize()
    return (torch.cat(out).numpy().astype(np.float64), float(state["best_threshold"]),
            time.perf_counter() - start)


@torch.no_grad()
def read_level_score(reads: np.ndarray, reference_index: np.ndarray, checkpoint: Path,
                     device: torch.device) -> tuple[np.ndarray, float, float]:
    """p(SNP), frozen threshold and inference seconds for the read-level model.

    The checkpoint records which control and which feature ablation it was
    trained under; both are re-applied here, so a control model is always
    evaluated on the same kind of input it was trained on.
    """
    state = torch.load(checkpoint, map_location=device, weights_only=False)
    model = ReadLevelMamba().to(device)
    model.load_state_dict(state["model_state_dict"])
    model.eval()

    control = state.get("control", "")
    if control:
        reads = apply_control(reads, control)
        logger.info("  re-applied control '%s' for %s", control, checkpoint.name)
    drop_indices = resolve_drop_indices(tuple(state.get("drop_groups", [])))

    windows = torch.from_numpy(np.ascontiguousarray(reads)).reshape(
        -1, SEQ_LEN, reads.shape[-2], reads.shape[-1])
    reference = torch.eye(5, dtype=torch.float32)[
        torch.tensor(reference_index.astype(np.int64)).clamp(0, 4)].reshape(-1, SEQ_LEN, 5)

    out = []
    if device.type == "cuda":
        torch.cuda.synchronize()
    start = time.perf_counter()
    for begin in range(0, windows.shape[0], 64):
        raw = windows[begin:begin + 64].to(device)
        valid = (raw[..., 6].to(torch.int64) & FLAG_VALID) > 0
        logits = model(read_features=build_features_torch(raw, drop_indices),
                       read_valid=valid,
                       reference_onehot=reference[begin:begin + 64].to(device))
        out.append(torch.softmax(logits, dim=-1)[..., LABEL_SNP].reshape(-1).cpu())
    if device.type == "cuda":
        torch.cuda.synchronize()
    return (torch.cat(out).numpy().astype(np.float64), float(state["best_threshold"]),
            time.perf_counter() - start)


def alternate_read_profile(reads: np.ndarray, counts: np.ndarray,
                           loci: np.ndarray) -> dict:
    """Profile the read-level evidence of the candidate alternate allele.

    The candidate alternate is the most-supported non-reference base, chosen
    from observed counts alone -- the same rule ``BinomialVariantCaller`` uses,
    and no truth is consulted.

    Args:
        reads: ``[N, R, C]`` uint8 read tensor.
        counts: ``[N, 10]`` aggregate counts (for the reference base).
        loci: Integer indices of the loci to profile.

    Returns:
        Median summaries of the alternate-supporting observations.
    """
    if loci.size == 0:
        return {"n_loci": 0}

    bases = counts[loci, 0:4].astype(np.float64)
    reference_column = counts[loci, 9].astype(int) - 1
    masked = np.where(np.arange(4)[None, :] == reference_column[:, None], -1.0, bases)
    alternate = masked.argmax(axis=1)

    subset = reads[loci]
    flags = subset[..., 6].astype(np.int64)
    counted = (flags & FLAG_COUNTED) > 0
    is_alternate = counted & (subset[..., 0].astype(np.int64) == alternate[:, None])
    is_reference = counted & (subset[..., 0].astype(np.int64) == reference_column[:, None])

    def summarize(selector: np.ndarray, column: int, scale: float = 1.0) -> dict:
        values = subset[..., column].astype(np.float64)[selector]
        if values.size == 0:
            return {"n": 0}
        return {"n": int(values.size), "mean": float(values.mean() * scale),
                "median": float(np.median(values) * scale)}

    alternate_count = is_alternate.sum(axis=1)
    forward = (is_alternate & (subset[..., 3] == 0)).sum(axis=1)
    with np.errstate(invalid="ignore"):
        strand_balance = np.where(alternate_count > 0,
                                  forward / np.maximum(alternate_count, 1), np.nan)
    near_end = (is_alternate & ((flags & FLAG_NEAR_END) > 0)).sum(axis=1)

    return {
        "n_loci": int(loci.size),
        "alt_base_quality": summarize(is_alternate, 1),
        "ref_base_quality": summarize(is_reference, 1),
        "alt_mapping_quality": summarize(is_alternate, 2),
        "alt_read_position": summarize(is_alternate, 4, 1 / 255.0),
        "alt_end_distance": summarize(is_alternate, 5, 1 / 255.0),
        "alt_forward_strand_fraction_median": float(np.nanmedian(strand_balance))
        if np.any(alternate_count > 0) else None,
        "alt_reads_near_end_fraction": float(
            near_end.sum() / max(alternate_count.sum(), 1)),
        "alt_reads_median": float(np.median(alternate_count)),
        "low_quality_alt_reads_median": float(np.median(
            ((subset[..., 0].astype(np.int64) == alternate[:, None])
             & ((flags & FLAG_VALID) > 0) & ~counted).sum(axis=1))),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-reads", required=True)
    parser.add_argument("--test-counts", required=True)
    parser.add_argument("--freq-threshold", type=float, required=True)
    parser.add_argument("--binom-v1-threshold", type=float, required=True)
    parser.add_argument("--binom-v2-threshold", type=float, required=True)
    parser.add_argument("--aggregate-checkpoint", required=True)
    parser.add_argument("--readlevel-checkpoint", required=True)
    parser.add_argument("--extra-checkpoint", nargs="*", default=[],
                        metavar="NAME=PATH", help="Control/ablation models to include.")
    parser.add_argument("--out", default="results/readlevel/readlevel_experiment_results.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    read_blob = np.load(args.test_reads, allow_pickle=True)
    count_blob = np.load(args.test_counts, allow_pickle=True)
    reads, labels = read_blob["reads"], read_blob["labels"]
    counts = count_blob["counts"]
    if not np.array_equal(count_blob["labels"], labels):
        raise SystemExit("read cache and count cache are not locus-aligned")

    snp = labels == LABEL_SNP
    primary = (labels == 0) | snp
    depth = counts[:, 6]
    vaf = vaf_from_counts(counts)
    reference_index = counts[:, 9]

    start = time.perf_counter()
    llr_v1 = raw_llr_from_counts(counts, BinomialVariantCaller())
    binomial_seconds = time.perf_counter() - start
    llr_v2 = raw_llr_from_counts(
        counts, BinomialVariantCaller(quality_derived_epsilon=True))

    aggregate_score, aggregate_threshold, aggregate_seconds = aggregate_model_score(
        base_features_from_counts(counts), Path(args.aggregate_checkpoint), device)
    readlevel, readlevel_threshold, readlevel_seconds = read_level_score(
        reads, reference_index, Path(args.readlevel_checkpoint), device)

    methods = {
        "frequency": (vaf, args.freq_threshold),
        "binomial_v1": (llr_v1, args.binom_v1_threshold),
        "binomial_v2": (llr_v2, args.binom_v2_threshold),
        "mamba_14ch": (aggregate_score, aggregate_threshold),
        "read_level": (readlevel, readlevel_threshold),
    }
    for entry in args.extra_checkpoint:
        name, path = entry.split("=", 1)
        score, threshold, _ = read_level_score(reads, reference_index, Path(path), device)
        methods[name] = (score, threshold)

    report: dict = {
        "test_reads": args.test_reads, "test_counts": args.test_counts,
        "thresholds": {name: float(t) for name, (_, t) in methods.items()},
        "threshold_provenance": "all frozen from validation; test never used for selection",
        "primary": {}, "secondary": {}, "by_depth": {}, "oracle_diagnostic": {},
        "checkpoints": {"mamba_14ch": args.aggregate_checkpoint,
                        "read_level": args.readlevel_checkpoint,
                        **dict(e.split("=", 1) for e in args.extra_checkpoint)},
    }

    for frame, mask in (("primary", primary), ("secondary", np.ones_like(snp, dtype=bool))):
        for name, (score, threshold) in methods.items():
            stats = prf(score[mask] >= threshold, snp[mask])
            stats.update(ranking_metrics(score[mask], snp[mask]))
            if frame == "primary":
                stats["f1_ci95"] = bootstrap_f1(score[mask] >= threshold, snp[mask])
            report[frame][name] = stats
        report[frame]["_loci"] = int(mask.sum())
        report[frame]["_snp"] = int(snp[mask].sum())

    binomial_call = llr_v1 >= args.binom_v1_threshold
    aggregate_call = aggregate_score >= aggregate_threshold
    report["vs_binomial_paired_bootstrap"] = {
        name: paired_bootstrap_f1_delta(score[primary] >= threshold,
                                        binomial_call[primary], snp[primary])
        for name, (score, threshold) in methods.items() if name != "binomial_v1"}
    report["vs_aggregate_mamba_paired_bootstrap"] = {
        name: paired_bootstrap_f1_delta(score[primary] >= threshold,
                                        aggregate_call[primary], snp[primary])
        for name, (score, threshold) in methods.items() if name != "mamba_14ch"}

    for name, (score, _) in methods.items():
        grid = (np.arange(-20.0, 200.0, 0.5) if name.startswith("binomial")
                else np.arange(0.0, 1.0, 0.005))
        best, best_threshold = {"f1": -1.0}, 0.0
        for candidate in grid:
            stats = prf(score[primary] >= candidate, snp[primary])
            if stats["f1"] > best["f1"]:
                best, best_threshold = stats, float(candidate)
        report["oracle_diagnostic"][name] = {"threshold": best_threshold, **best}

    for label, bin_mask in depth_bin_masks(depth):
        combined = bin_mask & primary
        entry = {"loci": int(combined.sum()), "snp": int(snp[combined].sum())}
        for name, (score, threshold) in methods.items():
            entry[name] = prf(score[combined] >= threshold, snp[combined])
        report["by_depth"][label] = entry

    # --- error analysis, including read-level evidence profiles ---
    model_call = readlevel >= readlevel_threshold
    binomial_fn = primary & snp & ~binomial_call
    model_fn = primary & snp & ~model_call
    binomial_fp = primary & ~snp & binomial_call
    model_fp = primary & ~snp & model_call
    model_only_recovered = binomial_fn & model_call
    model_only_fp = model_fp & ~binomial_call

    error_analysis = {
        "binomial_fn": int(binomial_fn.sum()), "model_fn": int(model_fn.sum()),
        "shared_fn": int((binomial_fn & model_fn).sum()),
        "model_only_fn": int((model_fn & binomial_call).sum()),
        "binomial_only_fn_recovered_by_model": int(model_only_recovered.sum()),
        "binomial_fp": int(binomial_fp.sum()), "model_fp": int(model_fp.sum()),
        "shared_fp": int((binomial_fp & model_fp).sum()),
        "model_only_fp": int(model_only_fp.sum()),
        "binomial_fp_corrected_by_model": int((binomial_fp & ~model_call).sum()),
        "binomial_tp_overturned_by_model": int(
            (primary & snp & binomial_call & ~model_call).sum()),
    }
    for tag, mask in (("model_only_recovered", model_only_recovered),
                      ("shared_fn", binomial_fn & model_fn),
                      ("model_only_fp", model_only_fp),
                      ("model_only_fn", model_fn & binomial_call),
                      ("true_positives_both", primary & snp & binomial_call & model_call)):
        loci = np.nonzero(mask)[0]
        error_analysis[f"{tag}_profile"] = {
            "count": int(loci.size),
            "depth_median": float(np.median(depth[loci])) if loci.size else None,
            "vaf_median": float(np.median(vaf[loci])) if loci.size else None,
            "llr_median": float(np.median(llr_v1[loci])) if loci.size else None,
            "read_level": alternate_read_profile(reads, counts, loci),
        }
    report["error_analysis"] = error_analysis

    report["runtime"] = {
        "binomial_scoring_seconds": binomial_seconds,
        "mamba_14ch_inference_seconds": aggregate_seconds,
        "read_level_inference_seconds": readlevel_seconds,
        "read_level_loci_per_second": float(labels.size / max(readlevel_seconds, 1e-9)),
        "device": str(device),
        "note": "model init/checkpoint load excluded from inference timings",
    }
    report["dataset"] = {
        "loci": int(labels.size), "snp": int(snp.sum()), "normal": int((labels == 0).sum()),
        "insertion": int((labels == 2).sum()), "deletion": int((labels == 3).sum()),
        "depth_median": float(np.median(depth))}

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))

    print("\n=== PRIMARY: SNP vs Normal (frozen validation thresholds) ===")
    for name in methods:
        s = report["primary"][name]
        print(f"  {name:<26} P={s['precision']:.4f} R={s['recall']:.4f} F1={s['f1']:.4f} "
              f"[{s['f1_ci95'][0]:.4f},{s['f1_ci95'][1]:.4f}] AUC={s['roc_auc']:.5f} "
              f"AP={s['pr_auc']:.4f} TP={s['tp']} FP={s['fp']} FN={s['fn']}")
    print("\n=== PAIRED delta F1 vs binomial v1 ===")
    for name, d in report["vs_binomial_paired_bootstrap"].items():
        print(f"  {name:<26} {d['delta_f1']:+.4f} [{d['ci95'][0]:+.4f},{d['ci95'][1]:+.4f}] "
              f"P(better)={d['fraction_of_draws_favouring_a']:.3f}")
    print("\n=== PAIRED delta F1 vs aggregate Mamba 14ch ===")
    for name, d in report["vs_aggregate_mamba_paired_bootstrap"].items():
        print(f"  {name:<26} {d['delta_f1']:+.4f} [{d['ci95'][0]:+.4f},{d['ci95'][1]:+.4f}] "
              f"P(better)={d['fraction_of_draws_favouring_a']:.3f}")
    print("\n=== ERROR ANALYSIS (read-level vs binomial) ===")
    print(json.dumps({k: v for k, v in error_analysis.items()
                      if not k.endswith("_profile")}, indent=2))
    logger.info("Wrote %s", args.out)


if __name__ == "__main__":
    main()
