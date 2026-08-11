"""Train RawPileupMambaLLR on cached 15x counts. One controlled variable: the LLR channel.

Every training hyperparameter is imported from ``train_raw_pileup`` rather than
restated, so the previous 15x run and this one cannot silently diverge: same
focal gamma, class-weight calibration, LR, batch size, seq_len, sequential
90/10 coordinate-order split, mutation-presence sampler, epoch count and
validation-threshold selection.

Two modes, selected by ``--shuffle-llr``:

* **real** -- channel 14 is the locus' own binomial LLR;
* **shuffled** -- channel 14 is a permutation of the same values across loci.
  The marginal distribution is preserved exactly while the association with
  the locus is destroyed. This is the control that distinguishes "the model
  uses the evidence" from "an extra input channel changed optimization".

The permutation seed is fixed and independent of any label or test statistic.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Subset, WeightedRandomSampler
from tqdm import tqdm

from binomial_baseline import BinomialVariantCaller
from config import LABEL_NORMAL
from llr_features import (
    FEATURE_DIM_LLR,
    LLRTransform,
    base_features_from_counts,
    raw_llr_from_counts,
)
from loss import FocalLoss, compute_calibrated_class_weights
from model_raw_pileup_llr import RawPileupMambaLLR
from threshold_selection import sweep_thresholds
from train_raw_pileup import (
    BATCH_SIZE,
    CLASS_WEIGHT_ALPHA,
    CLASS_WEIGHT_CAP,
    FOCAL_GAMMA,
    LR,
    NUM_CLASSES,
    SEQ_LEN,
    TRAIN_FRACTION,
    report,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


@torch.no_grad()
def collect_logits(model, batches, device):
    """Inference over a list of (features, labels) batches."""
    model.eval()
    logits_out, labels_out = [], []
    for features, labels in batches:
        logits = model(pileup_features=features.to(device))
        logits_out.append(logits.reshape(-1, logits.size(-1)).cpu())
        labels_out.append(labels.reshape(-1))
    return torch.cat(logits_out), torch.cat(labels_out)


def to_batches(features: np.ndarray, labels: np.ndarray, batch_size: int):
    """Reshape flat loci into window batches of ``[B, SEQ_LEN, F]``."""
    windows = torch.tensor(features, dtype=torch.float).reshape(-1, SEQ_LEN, features.shape[-1])
    window_labels = torch.tensor(labels, dtype=torch.long).reshape(-1, SEQ_LEN)
    return [(windows[i:i + batch_size], window_labels[i:i + batch_size])
            for i in range(0, windows.shape[0], batch_size)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-counts", required=True)
    parser.add_argument("--test-counts", required=True)
    parser.add_argument("--epochs", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--shuffle-llr", action="store_true",
                        help="Control run: permute the LLR channel across loci.")
    parser.add_argument("--shuffle-seed", type=int, default=777)
    parser.add_argument("--quality-epsilon", action="store_true",
                        help="Use binomial v2 (quality-derived epsilon) for the LLR channel.")
    parser.add_argument("--tag", required=True)
    parser.add_argument("--checkpoint-dir", default="checkpoints_rawpileup_llr")
    parser.add_argument("--out", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint_dir = Path(args.checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    train_blob = np.load(args.train_counts, allow_pickle=True)
    test_blob = np.load(args.test_counts, allow_pickle=True)
    train_counts, train_labels = train_blob["counts"], train_blob["labels"]
    test_counts, test_labels = test_blob["counts"], test_blob["labels"]

    caller = BinomialVariantCaller(
        quality_derived_epsilon=args.quality_epsilon, include_homozygous=True)
    logger.info("LLR channel from binomial %s", "v2 (quality-derived eps)"
                if args.quality_epsilon else "v1 (fixed eps=0.01)")

    # Sequential split by window, matching train_raw_pileup.py exactly.
    n_windows = train_labels.size // SEQ_LEN
    n_batches = (n_windows + BATCH_SIZE - 1) // BATCH_SIZE
    split_batch = int(n_batches * TRAIN_FRACTION)
    split_locus = split_batch * BATCH_SIZE * SEQ_LEN
    logger.info("Training region: %d loci / %d windows; fit=%d loci, val=%d loci",
                train_labels.size, n_windows, split_locus, train_labels.size - split_locus)

    fit_counts = train_counts[:split_locus]
    val_counts = train_counts[split_locus:]

    # Transform fitted on the FIT split only -- never validation, never test.
    fit_llr_raw = raw_llr_from_counts(fit_counts, caller)
    transform = LLRTransform.fit(fit_llr_raw)
    logger.info("LLR transform (fit split only): %s", transform.to_dict())
    logger.info("raw LLR on fit split: min=%.2f p50=%.2f p99=%.2f max=%.2f",
                fit_llr_raw.min(), np.median(fit_llr_raw),
                np.percentile(fit_llr_raw, 99), fit_llr_raw.max())

    def build(counts: np.ndarray, name: str) -> np.ndarray:
        base = base_features_from_counts(counts)
        llr = transform.apply(raw_llr_from_counts(counts, caller))
        if args.shuffle_llr:
            rng = np.random.default_rng(args.shuffle_seed)
            llr = rng.permutation(llr)
            logger.info("  [%s] LLR channel PERMUTED (marginal preserved, association destroyed)",
                        name)
        return np.concatenate([base, llr[:, None]], axis=1)

    fit_features = build(fit_counts, "fit")
    val_features = build(val_counts, "val")
    test_features = build(test_counts, "test")
    assert fit_features.shape[-1] == FEATURE_DIM_LLR

    fit_batches = to_batches(fit_features, train_labels[:split_locus], BATCH_SIZE)
    val_batches = to_batches(val_features, train_labels[split_locus:], BATCH_SIZE)
    test_batches = to_batches(test_features, test_labels, BATCH_SIZE)

    # Mutation-presence sampler, same rationale and construction as before.
    fit_windows = torch.tensor(fit_features, dtype=torch.float).reshape(-1, SEQ_LEN,
                                                                       FEATURE_DIM_LLR)
    fit_window_labels = torch.tensor(train_labels[:split_locus],
                                     dtype=torch.long).reshape(-1, SEQ_LEN)
    has_mutation = (fit_window_labels != LABEL_NORMAL).any(dim=1)
    n_mut = int(has_mutation.sum())
    n_bg = int(has_mutation.numel() - n_mut)
    logger.info("Fit windows with a mutation: %d | pure background: %d", n_mut, n_bg)
    weights = torch.where(has_mutation, 1.0 / max(n_mut, 1), 1.0 / max(n_bg, 1))

    window_ds = [{"pileup_features": fit_windows[i], "labels": fit_window_labels[i]}
                 for i in range(fit_windows.shape[0])]
    sampler = WeightedRandomSampler(weights, num_samples=len(window_ds), replacement=True)
    fit_loader = DataLoader(window_ds, batch_size=BATCH_SIZE, sampler=sampler)

    model = RawPileupMambaLLR().to(device)
    logger.info("Model parameters: %d", sum(p.numel() for p in model.parameters()))

    class_weights = compute_calibrated_class_weights(
        Subset(window_ds, range(min(200, len(window_ds)))), num_classes=NUM_CLASSES,
        sample_size=min(200, len(window_ds)), alpha=CLASS_WEIGHT_ALPHA,
        max_weight_cap=CLASS_WEIGHT_CAP).to(device)
    logger.info("Calibrated class weights: %s", class_weights.tolist())

    criterion = FocalLoss(alpha=class_weights, gamma=FOCAL_GAMMA, reduction="weighted_mean")
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR)

    history, best_score, best_threshold = [], float("-inf"), 0.5
    for epoch in range(1, args.epochs + 1):
        model.train()
        total, nb, grad_norm = 0.0, 0, 0.0
        for batch in tqdm(fit_loader, desc=f"epoch {epoch}", leave=False):
            features = batch["pileup_features"].to(device)
            labels = batch["labels"].to(device).reshape(-1)
            optimizer.zero_grad()
            logits = model(pileup_features=features)
            loss = criterion(logits.reshape(-1, NUM_CLASSES), labels)
            loss.backward()
            grad_norm = float(sum(p.grad.norm() ** 2 for p in model.parameters()
                                  if p.grad is not None) ** 0.5)
            optimizer.step()
            total += loss.item()
            nb += 1

        val_logits, val_labels_flat = collect_logits(model, val_batches, device)
        sweep = sweep_thresholds(val_logits, val_labels_flat)
        logger.info("epoch %d | train_loss=%.5f grad_norm=%.4f | val f1@opt=%.4f thr=%.2f",
                    epoch, total / max(nb, 1), grad_norm, sweep.f1_at_best, sweep.best_threshold)
        val_metrics = report(val_logits, val_labels_flat, sweep.best_threshold,
                             f"  VAL epoch{epoch}")
        history.append({"epoch": epoch, "train_loss": total / max(nb, 1),
                        "grad_norm": grad_norm, "val": val_metrics})

        payload = {"epoch": epoch, "model_state_dict": model.state_dict(),
                   "best_threshold": sweep.best_threshold, "val": val_metrics,
                   "representation": "raw_pileup_plus_llr_v1",
                   "llr_transform": transform.to_dict(),
                   "shuffled_llr": bool(args.shuffle_llr),
                   "quality_epsilon": bool(args.quality_epsilon)}
        torch.save(payload, checkpoint_dir / f"epoch{epoch}_{args.tag}.pt")
        if sweep.f1_at_best > best_score:
            best_score, best_threshold = sweep.f1_at_best, sweep.best_threshold
            torch.save(payload, checkpoint_dir / f"best_{args.tag}.pt")
    torch.save(payload, checkpoint_dir / f"last_{args.tag}.pt")

    best = torch.load(checkpoint_dir / f"best_{args.tag}.pt", map_location=device,
                      weights_only=False)
    model.load_state_dict(best["model_state_dict"])
    logger.info("Frozen validation threshold from epoch %d: %.2f", best["epoch"],
                best["best_threshold"])
    test_logits, test_labels_flat = collect_logits(model, test_batches, device)
    test_metrics = report(test_logits, test_labels_flat, best["best_threshold"],
                          "  TEST (frozen val thr)")
    oracle = sweep_thresholds(test_logits, test_labels_flat)
    logger.info("  [diagnostic upper bound only] test-optimal thr=%.2f f1=%.4f",
                oracle.best_threshold, oracle.f1_at_best)

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(
            {"history": history, "frozen_threshold": best["best_threshold"],
             "best_epoch": best["epoch"], "test": test_metrics,
             "llr_transform": transform.to_dict(), "shuffled_llr": bool(args.shuffle_llr),
             "diagnostic_test_optimal": {"threshold": oracle.best_threshold,
                                         "f1": oracle.f1_at_best}}, indent=2))


if __name__ == "__main__":
    main()
