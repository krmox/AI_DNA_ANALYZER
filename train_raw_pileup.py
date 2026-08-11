"""Train and evaluate the raw-pileup model on real GIAB HG002 data.

Mirrors ``run_train.py``'s scientific settings exactly -- same focal gamma,
class-weight calibration, LR, batch size, seq_len, sequential coordinate-order
splitting, mutation-presence sampler, and per-checkpoint optimal-threshold
selection -- so the only changed variable against the previous experiment is
the input representation.

Three disjoint coordinate ranges, all in absolute GRCh38 chr21:

* train / validation: one region, split sequentially (no shuffling, so
  neighbouring correlated windows cannot straddle the boundary);
* test: a separate region the model never sees during training or selection.

Thresholds are selected on validation only and frozen before touching test,
matching the discipline established for the token-model experiments.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import torch
from torch.utils.data import DataLoader, Subset, WeightedRandomSampler
from tqdm import tqdm

from config import LABEL_NAMES, LABEL_NORMAL, LABEL_SNP
from loss import FocalLoss, compute_calibrated_class_weights
from model_raw_pileup import RawPileupMamba
from raw_pileup import RawPileupDataset, RawPileupProvider
from run_train import (
    average_precision,
    binary_auc,
    per_class_precision_recall_f1,
    predict_with_threshold,
    torch_confusion_matrix,
)
from threshold_selection import sweep_thresholds

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

NUM_CLASSES = len(LABEL_NAMES)

# Unchanged from run_train.py.
SEQ_LEN = 64
BATCH_SIZE = 16
LR = 3e-4
TRAIN_FRACTION = 0.9
FOCAL_GAMMA = 3.0
CLASS_WEIGHT_ALPHA = 0.75
CLASS_WEIGHT_CAP = 20.0

CHECKPOINT_DIR = Path("checkpoints_rawpileup")


def materialize(provider_args: dict, batch_size: int) -> list[dict[str, torch.Tensor]]:
    """Build every window of a region once and keep the batches in memory.

    The pileup walk dominates runtime, so paying it once per region rather
    than once per epoch is what makes multi-epoch training practical here.

    Args:
        provider_args: Keyword arguments for :class:`RawPileupProvider`.
        batch_size: Batch size to collate at.

    Returns:
        List of collated batch dicts in coordinate order.
    """
    provider = RawPileupProvider(**provider_args)
    with provider:
        dataset = RawPileupDataset(provider)
        logger.info("Region %s -> %d windows", provider_args["region"], len(dataset))
        return list(DataLoader(dataset, batch_size=batch_size, shuffle=False))


def flatten(batches: list[dict[str, torch.Tensor]], key: str) -> torch.Tensor:
    """Concatenate one key across batches into a flat tensor."""
    if key == "pileup_features":
        return torch.cat([b[key].reshape(-1, b[key].shape[-1]) for b in batches])
    return torch.cat([b[key].reshape(-1) for b in batches])


@torch.no_grad()
def collect_logits(
    model: torch.nn.Module, batches: list[dict[str, torch.Tensor]], device: torch.device
) -> tuple[torch.Tensor, torch.Tensor]:
    """Run inference over batches and buffer flat logits and labels."""
    model.eval()
    logits_out, labels_out = [], []
    for batch in batches:
        features = batch["pileup_features"].to(device)
        logits = model(pileup_features=features)
        logits_out.append(logits.reshape(-1, logits.size(-1)).cpu())
        labels_out.append(batch["labels"].reshape(-1).cpu())
    return torch.cat(logits_out), torch.cat(labels_out)


def report(logits: torch.Tensor, labels: torch.Tensor, threshold: float, name: str) -> dict:
    """Compute and log the full metric bundle at one threshold."""
    preds = predict_with_threshold(logits, threshold)
    confusion = torch_confusion_matrix(labels, preds, NUM_CLASSES)
    precision, recall, f1 = per_class_precision_recall_f1(confusion)
    scores = torch.softmax(logits, dim=-1)[:, LABEL_SNP].float()
    positives = labels == LABEL_SNP
    snp_true, snp_pred = labels == LABEL_SNP, preds == LABEL_SNP
    metrics = {
        "threshold": threshold,
        "snp_auc": binary_auc(scores, positives),
        "snp_ap": average_precision(scores, positives),
        "snp_precision": float(precision[LABEL_SNP]),
        "snp_recall": float(recall[LABEL_SNP]),
        "snp_f1": float(f1[LABEL_SNP]),
        "mutation_macro_f1": float(f1[1:].mean()),
        "tp": int((snp_true & snp_pred).sum()),
        "fp": int((~snp_true & snp_pred).sum()),
        "fn": int((snp_true & ~snp_pred).sum()),
        "tn": int((~snp_true & ~snp_pred).sum()),
        "confusion": confusion.tolist(),
    }
    logger.info(
        "%s @thr=%.2f | auc=%.4f ap=%.4f P=%.4f R=%.4f F1=%.4f TP=%d FP=%d FN=%d",
        name, threshold, metrics["snp_auc"], metrics["snp_ap"], metrics["snp_precision"],
        metrics["snp_recall"], metrics["snp_f1"], metrics["tp"], metrics["fp"], metrics["fn"],
    )
    return metrics


def parse_args() -> argparse.Namespace:
    """CLI for the raw-pileup experiment."""
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
    parser.add_argument("--epochs", type=int, default=6)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--tag", default="rawpileup")
    parser.add_argument("--max-windows", type=int, default=None)
    parser.add_argument("--checkpoint-dir", default=str(CHECKPOINT_DIR),
                        help="Directory for this run's checkpoints; keeps regimes isolated.")
    parser.add_argument("--out", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    torch.manual_seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint_dir = Path(args.checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Device: %s | seed=%d | epochs=%d", device, args.seed, args.epochs)

    common = dict(fasta_path=args.fasta, contig="chr21", seq_len=SEQ_LEN,
                  max_windows=args.max_windows)
    train_batches = materialize(
        dict(common, bam_path=args.train_bam, vcf_path=args.train_vcf,
             high_confidence_bed=args.train_bed, region=tuple(args.train_region)), BATCH_SIZE)
    test_batches = materialize(
        dict(common, bam_path=args.test_bam, vcf_path=args.test_vcf,
             high_confidence_bed=args.test_bed, region=tuple(args.test_region)), BATCH_SIZE)

    # Sequential split of the training region: first 90% train, last 10% val.
    split = int(len(train_batches) * TRAIN_FRACTION)
    fit_batches, val_batches = train_batches[:split], train_batches[split:]
    logger.info("Batches -> fit %d | val %d | test %d",
                len(fit_batches), len(val_batches), len(test_batches))

    # Mutation-presence sampler weights, same rationale as run_train.py.
    window_labels = [w for b in fit_batches for w in b["labels"]]
    has_mutation = torch.tensor([bool((w != LABEL_NORMAL).any()) for w in window_labels])
    n_mut = int(has_mutation.sum())
    n_bg = len(has_mutation) - n_mut
    logger.info("Fit windows with a mutation: %d | pure background: %d", n_mut, n_bg)
    weights = torch.where(has_mutation, 1.0 / max(n_mut, 1), 1.0 / max(n_bg, 1))

    flat_features = torch.cat([b["pileup_features"] for b in fit_batches])
    flat_labels = torch.cat([b["labels"] for b in fit_batches])
    window_ds = [{"pileup_features": flat_features[i], "labels": flat_labels[i]}
                 for i in range(flat_labels.shape[0])]
    sampler = WeightedRandomSampler(weights, num_samples=len(window_ds), replacement=True)
    fit_loader = DataLoader(window_ds, batch_size=BATCH_SIZE, sampler=sampler)

    model = RawPileupMamba().to(device)
    logger.info("Model parameters: %d", sum(p.numel() for p in model.parameters()))

    class_weights = compute_calibrated_class_weights(
        Subset(window_ds, range(min(200, len(window_ds)))), num_classes=NUM_CLASSES,
        sample_size=min(200, len(window_ds)), alpha=CLASS_WEIGHT_ALPHA,
        max_weight_cap=CLASS_WEIGHT_CAP,
    ).to(device)
    logger.info("Calibrated class weights: %s", class_weights.tolist())

    criterion = FocalLoss(alpha=class_weights, gamma=FOCAL_GAMMA, reduction="weighted_mean")
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR)

    history = []
    best_score, best_threshold = float("-inf"), 0.5
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

        val_logits, val_labels = collect_logits(model, val_batches, device)
        sweep = sweep_thresholds(val_logits, val_labels)
        logger.info("epoch %d | train_loss=%.5f grad_norm=%.4f | val f1@0.5=%.4f "
                    "f1@opt=%.4f thr=%.2f", epoch, total / max(nb, 1), grad_norm,
                    sweep.f1_at_fixed, sweep.f1_at_best, sweep.best_threshold)
        val_metrics = report(val_logits, val_labels, sweep.best_threshold, f"  VAL epoch{epoch}")
        history.append({"epoch": epoch, "train_loss": total / max(nb, 1),
                        "grad_norm": grad_norm, "val": val_metrics})

        torch.save({"epoch": epoch, "model_state_dict": model.state_dict(),
                    "best_threshold": sweep.best_threshold, "val": val_metrics,
                    "representation": "raw_pileup_v1"},
                   checkpoint_dir / f"epoch{epoch}_{args.tag}.pt")

        if sweep.f1_at_best > best_score:
            best_score, best_threshold = sweep.f1_at_best, sweep.best_threshold
            torch.save({"epoch": epoch, "model_state_dict": model.state_dict(),
                        "best_threshold": best_threshold, "val": val_metrics,
                        "representation": "raw_pileup_v1"},
                       checkpoint_dir / f"best_{args.tag}.pt")

    # Test, once, with the frozen validation threshold.
    best = torch.load(checkpoint_dir / f"best_{args.tag}.pt", map_location=device,
                      weights_only=False)
    model.load_state_dict(best["model_state_dict"])
    frozen = best["best_threshold"]
    logger.info("Frozen validation threshold from epoch %d: %.2f", best["epoch"], frozen)
    test_logits, test_labels = collect_logits(model, test_batches, device)
    test_metrics = report(test_logits, test_labels, frozen, "  TEST (frozen val thr)")
    upper = sweep_thresholds(test_logits, test_labels)
    logger.info("  [diagnostic upper bound only] test-optimal thr=%.2f f1=%.4f",
                upper.best_threshold, upper.f1_at_best)

    if args.out:
        Path(args.out).write_text(json.dumps(
            {"history": history, "frozen_threshold": frozen, "test": test_metrics,
             "diagnostic_test_optimal": {"threshold": upper.best_threshold,
                                         "f1": upper.f1_at_best}}, indent=2))


if __name__ == "__main__":
    main()
