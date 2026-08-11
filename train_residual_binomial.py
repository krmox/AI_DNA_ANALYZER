"""Train the residual binomial Mamba on cached 15x counts.

    final_logits = binomial_prior_logits + residual(14-channel pileup)

Everything that can be inherited from the previous 15x experiments *is*
inherited rather than restated -- ``SEQ_LEN``, ``BATCH_SIZE``, ``LR``,
``TRAIN_FRACTION``, ``FOCAL_GAMMA``, ``CLASS_WEIGHT_ALPHA``,
``CLASS_WEIGHT_CAP`` all come from ``train_raw_pileup``, and the split,
sampler, optimizer, epoch count and seed are unchanged. The controlled
variable is the residual parameterization.

Two deliberate, documented protocol deviations
----------------------------------------------
1. **Threshold selection happens in logit space, not probability space.**
   ``softmax(prior)[SNP]`` saturates to exactly 1.0 in float32 above
   LLR ~= 17, which would tie together the large majority of true SNPs and
   make a probability threshold unable to express the binomial's operating
   point at all (LLR 10.5 is p = 0.99992, beyond the old 0.50-0.99 grid).
   This is a representational impossibility, not a preference, so the score
   is the SNP-minus-Normal logit gap and the grid is
   ``arange(-20, 200, 0.5)`` -- the *same* grid ``evaluate_binomial_baseline``
   used to pick the binomial's own 10.5. The residual model and the binomial
   baseline are therefore selected by an identical procedure on an identical
   axis, which is what makes the comparison fair.

2. **The selection objective is SNP F1, not mutation_macro_f1.** The binomial
   baseline's threshold was selected on SNP F1; matching it keeps the two
   methods' operating points comparable. SNP-vs-Normal is the primary metric
   in any case.

Neither deviation touches the loss, the data, the split or the optimizer.

The residual penalty ``lambda * mean(residual^2)`` is swept over
{0, 1e-4, 1e-3, 1e-2}; selection uses validation only.

Control mode ``--shuffle-prior`` permutes the LLR across loci before it
becomes the prior. The marginal distribution of the prior is preserved
exactly while its association with the locus is destroyed, so a residual model
that improves in the real condition but not the shuffled one is exploiting the
statistic rather than the extra degrees of freedom.
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
from config import LABEL_NORMAL, LABEL_SNP
from llr_features import base_features_from_counts, raw_llr_from_counts
from loss import FocalLoss, compute_calibrated_class_weights
from model_residual_binomial import ResidualBinomialMamba
from residual_metrics import prf, ranking_metrics
from residual_prior import NUM_CLASSES, prior_logits_from_llr, snp_decision_score
from train_raw_pileup import (
    BATCH_SIZE,
    CLASS_WEIGHT_ALPHA,
    CLASS_WEIGHT_CAP,
    FOCAL_GAMMA,
    LR,
    SEQ_LEN,
    TRAIN_FRACTION,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

#: Identical to the grid ``evaluate_binomial_baseline`` used for binomial v1.
LLR_THRESHOLD_GRID = np.arange(-20.0, 200.0, 0.5)


def select_threshold(score: np.ndarray, labels: np.ndarray) -> tuple[float, dict]:
    """Pick the SNP-F1-maximizing threshold on the SNP-vs-Normal frame.

    Args:
        score: ``[N]`` decision scores (logit gap).
        labels: ``[N]`` integer class labels.

    Returns:
        ``(threshold, stats_at_that_threshold)``.
    """
    mask = (labels == LABEL_NORMAL) | (labels == LABEL_SNP)
    snp = labels[mask] == LABEL_SNP
    masked_score = score[mask]
    best, best_threshold = {"f1": -1.0}, 0.0
    for candidate in LLR_THRESHOLD_GRID:
        stats = prf(masked_score >= candidate, snp)
        if stats["f1"] > best["f1"]:
            best, best_threshold = stats, float(candidate)
    return best_threshold, best


@torch.no_grad()
def collect_scores(model, features, prior, device, batch_size: int = 256):
    """Decision score and mean squared residual over a whole split."""
    model.eval()
    scores, squared = [], []
    for begin in range(0, features.shape[0], batch_size):
        logits, residual = model(
            pileup_features=features[begin:begin + batch_size].to(device),
            prior_logits=prior[begin:begin + batch_size].to(device),
            return_residual=True)
        scores.append(snp_decision_score(logits).reshape(-1).cpu())
        squared.append(residual.pow(2).mean().cpu())
    return torch.cat(scores).numpy().astype(np.float64), float(torch.stack(squared).mean())


def evaluate(model, features, prior, labels, device, threshold=None) -> tuple[dict, float]:
    """Full metric bundle on one split, optionally at a supplied threshold."""
    score, mean_sq_residual = collect_scores(model, features, prior, device)
    if threshold is None:
        threshold, _ = select_threshold(score, labels)
    mask = (labels == LABEL_NORMAL) | (labels == LABEL_SNP)
    snp = labels[mask] == LABEL_SNP
    metrics = {"threshold": float(threshold),
               **prf(score[mask] >= threshold, snp),
               **ranking_metrics(score[mask], snp),
               "mean_squared_residual": mean_sq_residual}
    return metrics, float(threshold)


def build_tensors(counts: np.ndarray, llr: np.ndarray):
    """Reshape flat loci into ``[W, SEQ_LEN, ...]`` feature and prior windows."""
    features = torch.tensor(base_features_from_counts(counts), dtype=torch.float)
    prior = torch.tensor(prior_logits_from_llr(llr), dtype=torch.float)
    return (features.reshape(-1, SEQ_LEN, features.shape[-1]),
            prior.reshape(-1, SEQ_LEN, NUM_CLASSES))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-counts", required=True)
    parser.add_argument("--epochs", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--residual-lambda", type=float, required=True)
    parser.add_argument("--shuffle-prior", action="store_true",
                        help="Control: permute the binomial LLR across loci.")
    parser.add_argument("--shuffle-seed", type=int, default=777)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--checkpoint-dir", default="checkpoints_residual_binomial")
    parser.add_argument("--out", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint_dir = Path(args.checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    blob = np.load(args.train_counts, allow_pickle=True)
    counts, labels = blob["counts"], blob["labels"]

    caller = BinomialVariantCaller()          # v1: fixed eps=0.01, max(het, hom)
    llr = raw_llr_from_counts(counts, caller)
    if args.shuffle_prior:
        llr = np.random.default_rng(args.shuffle_seed).permutation(llr)
        logger.info("CONTROL: binomial prior PERMUTED across loci "
                    "(marginal preserved, association destroyed)")

    # Sequential split by batch, byte-identical to train_raw_pileup_llr.py.
    n_windows = labels.size // SEQ_LEN
    n_batches = (n_windows + BATCH_SIZE - 1) // BATCH_SIZE
    split_locus = int(n_batches * TRAIN_FRACTION) * BATCH_SIZE * SEQ_LEN
    logger.info("Training region: %d loci / %d windows; fit=%d loci, val=%d loci",
                labels.size, n_windows, split_locus, labels.size - split_locus)

    fit_features, fit_prior = build_tensors(counts[:split_locus], llr[:split_locus])
    val_features, val_prior = build_tensors(counts[split_locus:], llr[split_locus:])
    fit_labels_flat, val_labels_flat = labels[:split_locus], labels[split_locus:]
    fit_window_labels = torch.tensor(fit_labels_flat, dtype=torch.long).reshape(-1, SEQ_LEN)

    # Mutation-presence sampler, same construction as the previous experiments.
    has_mutation = (fit_window_labels != LABEL_NORMAL).any(dim=1)
    n_mut = int(has_mutation.sum())
    n_bg = int(has_mutation.numel() - n_mut)
    logger.info("Fit windows with a mutation: %d | pure background: %d", n_mut, n_bg)
    weights = torch.where(has_mutation, 1.0 / max(n_mut, 1), 1.0 / max(n_bg, 1))

    window_ds = [{"pileup_features": fit_features[i], "prior_logits": fit_prior[i],
                  "labels": fit_window_labels[i]} for i in range(fit_features.shape[0])]
    sampler = WeightedRandomSampler(weights, num_samples=len(window_ds), replacement=True)
    fit_loader = DataLoader(window_ds, batch_size=BATCH_SIZE, sampler=sampler)

    model = ResidualBinomialMamba().to(device)
    logger.info("Model parameters: %d | residual lambda=%g | tag=%s",
                sum(p.numel() for p in model.parameters()), args.residual_lambda, args.tag)

    class_weights = compute_calibrated_class_weights(
        Subset(window_ds, range(min(200, len(window_ds)))), num_classes=NUM_CLASSES,
        sample_size=min(200, len(window_ds)), alpha=CLASS_WEIGHT_ALPHA,
        max_weight_cap=CLASS_WEIGHT_CAP).to(device)
    logger.info("Calibrated class weights: %s", class_weights.tolist())

    criterion = FocalLoss(alpha=class_weights, gamma=FOCAL_GAMMA, reduction="weighted_mean")
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR)

    # Epoch 0: the untrained model, which by construction is the binomial caller.
    epoch0, _ = evaluate(model, val_features, val_prior, val_labels_flat, device)
    logger.info("epoch 0 (untrained == binomial) | val F1=%.4f AUC=%.4f mean_res^2=%.3e",
                epoch0["f1"], epoch0["roc_auc"], epoch0["mean_squared_residual"])

    history, best_score, best_threshold, payload = [], epoch0["f1"], epoch0["threshold"], None
    torch.save({"epoch": 0, "model_state_dict": model.state_dict(),
                "best_threshold": epoch0["threshold"], "val": epoch0,
                "representation": "residual_binomial_v1",
                "residual_lambda": args.residual_lambda,
                "shuffled_prior": bool(args.shuffle_prior)},
               checkpoint_dir / f"best_{args.tag}.pt")

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss, total_task, total_penalty, nb, grad_norm = 0.0, 0.0, 0.0, 0, 0.0
        for batch in tqdm(fit_loader, desc=f"epoch {epoch}", leave=False):
            features = batch["pileup_features"].to(device)
            prior = batch["prior_logits"].to(device)
            batch_labels = batch["labels"].to(device).reshape(-1)
            optimizer.zero_grad()
            logits, residual = model(pileup_features=features, prior_logits=prior,
                                     return_residual=True)
            task_loss = criterion(logits.reshape(-1, NUM_CLASSES), batch_labels)
            penalty = args.residual_lambda * residual.pow(2).mean()
            loss = task_loss + penalty
            loss.backward()
            grad_norm = float(sum(p.grad.norm() ** 2 for p in model.parameters()
                                  if p.grad is not None) ** 0.5)
            optimizer.step()
            total_loss += float(loss.item())
            total_task += float(task_loss.item())
            total_penalty += float(penalty.item() if torch.is_tensor(penalty) else penalty)
            nb += 1

        val_metrics, epoch_threshold = evaluate(
            model, val_features, val_prior, val_labels_flat, device)
        logger.info(
            "epoch %d | loss=%.5f (task=%.5f pen=%.3e) grad=%.3f | "
            "VAL thr=%.1f F1=%.4f P=%.4f R=%.4f AUC=%.4f AP=%.4f mean_res^2=%.3e",
            epoch, total_loss / max(nb, 1), total_task / max(nb, 1),
            total_penalty / max(nb, 1), grad_norm, val_metrics["threshold"],
            val_metrics["f1"], val_metrics["precision"], val_metrics["recall"],
            val_metrics["roc_auc"], val_metrics["pr_auc"],
            val_metrics["mean_squared_residual"])
        history.append({"epoch": epoch, "train_loss": total_loss / max(nb, 1),
                        "task_loss": total_task / max(nb, 1),
                        "residual_penalty": total_penalty / max(nb, 1),
                        "grad_norm": grad_norm, "val": val_metrics})

        payload = {"epoch": epoch, "model_state_dict": model.state_dict(),
                   "best_threshold": epoch_threshold, "val": val_metrics,
                   "representation": "residual_binomial_v1",
                   "residual_lambda": args.residual_lambda,
                   "shuffled_prior": bool(args.shuffle_prior),
                   "shuffle_seed": args.shuffle_seed, "seed": args.seed}
        torch.save(payload, checkpoint_dir / f"epoch{epoch}_{args.tag}.pt")
        if val_metrics["f1"] > best_score:
            best_score, best_threshold = val_metrics["f1"], epoch_threshold
            torch.save(payload, checkpoint_dir / f"best_{args.tag}.pt")
    if payload is not None:
        torch.save(payload, checkpoint_dir / f"last_{args.tag}.pt")

    best = torch.load(checkpoint_dir / f"best_{args.tag}.pt", map_location=device,
                      weights_only=False)
    logger.info("Selected epoch %d on VALIDATION | val F1=%.4f | frozen threshold %.1f",
                best["epoch"], best_score, best["best_threshold"])

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps({
            "tag": args.tag, "residual_lambda": args.residual_lambda,
            "shuffled_prior": bool(args.shuffle_prior), "seed": args.seed,
            "epochs": args.epochs, "epoch0_untrained_binomial": epoch0,
            "history": history, "selected_epoch": best["epoch"],
            "selected_val_f1": best_score, "frozen_threshold": best["best_threshold"],
            "selection": "validation split only; test never consulted",
        }, indent=2))


if __name__ == "__main__":
    main()
