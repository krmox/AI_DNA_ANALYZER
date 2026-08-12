"""Train the Poisson-binomial residual Mamba on the 44-channel locus evidence.

    final_logits = PB_prior_logits + [gate(x) *] residual(44-channel evidence)

Everything that can be inherited from the previous residual experiment *is*
inherited rather than restated -- ``SEQ_LEN``, ``BATCH_SIZE``, ``LR``,
``TRAIN_FRACTION``, ``FOCAL_GAMMA``, ``CLASS_WEIGHT_ALPHA``,
``CLASS_WEIGHT_CAP``, the sequential split, the mutation-presence sampler, the
optimizer, the epoch count, the threshold grid and the selection objective all
come from ``train_raw_pileup`` / ``train_residual_binomial`` unchanged.

The controlled variable in a Variant 1 run is therefore **the input
representation alone**: 44 channels of count/quality/strand evidence in place
of 14 channels of fractions, and the Poisson-binomial LLR in place of the
fixed-epsilon binomial LLR as the prior. ``--use-gate`` adds the Variant 2
parameterization on top and is off by default.

Epoch 0 is evaluated before any gradient step. By construction the untrained
model is *exactly* the Poisson-binomial caller, so epoch 0 is the baseline
measured on this experiment's own validation split, with this experiment's own
threshold grid -- no cross-run comparison is needed to know what the model has
to beat.

Diagnostics that decide interpretation
--------------------------------------
The previous residual experiment failed by learning an almost-constant shift.
Two quantities are logged every epoch so that outcome is visible rather than
inferred:

* ``residual_gap_mean`` / ``residual_gap_std`` -- the mean and the across-loci
  standard deviation of the residual's SNP-minus-Normal component. A constant
  shift is exactly ``std ~ 0`` with ``mean != 0``. Locus-specific correction is
  ``std`` comparable to or larger than ``|mean|``.
* ``gate_mean`` / ``gate_std`` -- present only with ``--use-gate``.

Controls
--------
``--shuffle-prior`` permutes the PB LLR across loci: the marginal is preserved
and the association destroyed. ``--drop-groups`` zeroes named feature groups
(``locus_evidence.FEATURE_GROUPS``) while keeping the tensor shape, so an
ablation measures information rather than capacity.

One documented protocol change: selection metric
------------------------------------------------
Previous experiments selected the checkpoint on validation SNP F1. That is no
longer usable here. The validation split holds **251 SNPs**, and the
Poisson-binomial prior -- which is where this model *starts* -- already scores
249 TP / 0 FP / 2 FN on it. Validation F1 therefore has a total dynamic range
of two loci, and choosing between checkpoints or architectures on it would be
selecting on noise.

Selection is consequently keyed on validation **average precision**, which
ranks all ~188k validation loci instead of counting two decision flips, and is
threshold-free so it cannot be gamed by an operating-point shift. F1 is still
computed, logged and reported every epoch; it is simply not the selector. This
choice is a property of the validation split's size and was made without
reference to any test label.
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

from config import LABEL_NORMAL
from locus_evidence import FEATURE_DIM, FEATURE_GROUPS, resolve_group_indices
from loss import FocalLoss, compute_calibrated_class_weights
from model_pb_residual import PoissonBinomialResidualMamba
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
from train_residual_binomial import LLR_THRESHOLD_GRID, select_threshold

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


@torch.no_grad()
def collect(model, features, prior, device, batch_size: int = 256) -> dict:
    """Decision scores and residual/gate diagnostics over a whole split.

    Args:
        model: The residual model.
        features: ``[W, L, FEATURE_DIM]`` feature windows.
        prior: ``[W, L, NUM_CLASSES]`` prior windows.
        device: Inference device.
        batch_size: Windows per forward pass.

    Returns:
        Dict with ``score`` and the residual/gate summary statistics.
    """
    model.eval()
    scores, gaps, gates = [], [], []
    for begin in range(0, features.shape[0], batch_size):
        logits, residual, gate = model(
            pileup_features=features[begin:begin + batch_size].to(device),
            prior_logits=prior[begin:begin + batch_size].to(device),
            return_parts=True)
        scores.append(snp_decision_score(logits).reshape(-1).cpu())
        gaps.append(snp_decision_score(residual).reshape(-1).cpu())
        gates.append(gate.reshape(-1).cpu())
    gap = torch.cat(gaps)
    gate_values = torch.cat(gates)
    return {
        "score": torch.cat(scores).numpy().astype(np.float64),
        "residual_gap_mean": float(gap.mean()),
        "residual_gap_std": float(gap.std()),
        "residual_gap_abs_mean": float(gap.abs().mean()),
        "gate_mean": float(gate_values.mean()),
        "gate_std": float(gate_values.std()),
    }


def evaluate(model, features, prior, labels, device, threshold=None) -> tuple[dict, float]:
    """Full metric bundle on one split, optionally at a supplied threshold."""
    collected = collect(model, features, prior, device)
    score = collected["score"]
    if threshold is None:
        threshold, _ = select_threshold(score, labels)
    mask = (labels == LABEL_NORMAL) | (labels == 1)
    snp = labels[mask] == 1
    metrics = {"threshold": float(threshold),
               **prf(score[mask] >= threshold, snp),
               **ranking_metrics(score[mask], snp),
               **{key: value for key, value in collected.items() if key != "score"}}
    return metrics, float(threshold)


def build_tensors(features: np.ndarray, llr: np.ndarray,
                  drop_indices: tuple[int, ...] = ()):
    """Reshape flat loci into ``[W, SEQ_LEN, ...]`` feature and prior windows."""
    features = np.array(features, dtype=np.float32, copy=True)
    for index in drop_indices:
        features[:, index] = 0.0
    feature_tensor = torch.tensor(features, dtype=torch.float)
    prior_tensor = torch.tensor(prior_logits_from_llr(llr), dtype=torch.float)
    return (feature_tensor.reshape(-1, SEQ_LEN, feature_tensor.shape[-1]),
            prior_tensor.reshape(-1, SEQ_LEN, NUM_CLASSES))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-evidence", default="cache/train_15x_evidence.npz")
    parser.add_argument("--epochs", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--residual-lambda", type=float, default=0.0)
    parser.add_argument("--use-gate", action="store_true",
                        help="Variant 2: final = prior + gate(x) * residual(x)")
    parser.add_argument("--shuffle-prior", action="store_true",
                        help="Control: permute the PB LLR across loci.")
    parser.add_argument("--shuffle-seed", type=int, default=777)
    parser.add_argument("--drop-groups", nargs="*", default=[],
                        choices=sorted(FEATURE_GROUPS), help="Ablation: zero these groups.")
    parser.add_argument("--tag", required=True)
    parser.add_argument("--checkpoint-dir", default="checkpoints_pb_residual")
    parser.add_argument("--out", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint_dir = Path(args.checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    blob = np.load(args.train_evidence, allow_pickle=True)
    features, labels, llr = blob["features"], blob["labels"], blob["pb_llr"]

    if args.shuffle_prior:
        llr = np.random.default_rng(args.shuffle_seed).permutation(llr)
        logger.info("CONTROL: PB prior PERMUTED across loci "
                    "(marginal preserved, association destroyed)")

    drop_indices = resolve_group_indices(tuple(args.drop_groups))
    if drop_indices:
        logger.info("ABLATION: zeroing feature groups %s (%d channels)",
                    args.drop_groups, len(drop_indices))

    # Sequential split by batch, byte-identical to train_residual_binomial.py.
    n_windows = labels.size // SEQ_LEN
    n_batches = (n_windows + BATCH_SIZE - 1) // BATCH_SIZE
    split_locus = int(n_batches * TRAIN_FRACTION) * BATCH_SIZE * SEQ_LEN
    logger.info("Training region: %d loci / %d windows; fit=%d loci, val=%d loci",
                labels.size, n_windows, split_locus, labels.size - split_locus)

    fit_features, fit_prior = build_tensors(features[:split_locus], llr[:split_locus],
                                            drop_indices)
    val_features, val_prior = build_tensors(features[split_locus:], llr[split_locus:],
                                            drop_indices)
    fit_labels_flat, val_labels_flat = labels[:split_locus], labels[split_locus:]
    fit_window_labels = torch.tensor(fit_labels_flat, dtype=torch.long).reshape(-1, SEQ_LEN)

    has_mutation = (fit_window_labels != LABEL_NORMAL).any(dim=1)
    n_mut = int(has_mutation.sum())
    n_bg = int(has_mutation.numel() - n_mut)
    logger.info("Fit windows with a mutation: %d | pure background: %d", n_mut, n_bg)
    weights = torch.where(has_mutation, 1.0 / max(n_mut, 1), 1.0 / max(n_bg, 1))

    window_ds = [{"pileup_features": fit_features[i], "prior_logits": fit_prior[i],
                  "labels": fit_window_labels[i]} for i in range(fit_features.shape[0])]
    sampler = WeightedRandomSampler(weights, num_samples=len(window_ds), replacement=True)
    fit_loader = DataLoader(window_ds, batch_size=BATCH_SIZE, sampler=sampler)

    model = PoissonBinomialResidualMamba(feature_dim=FEATURE_DIM,
                                         use_gate=args.use_gate).to(device)
    logger.info("Model parameters: %d | gate=%s | lambda=%g | tag=%s",
                sum(p.numel() for p in model.parameters()), args.use_gate,
                args.residual_lambda, args.tag)

    class_weights = compute_calibrated_class_weights(
        Subset(window_ds, range(min(200, len(window_ds)))), num_classes=NUM_CLASSES,
        sample_size=min(200, len(window_ds)), alpha=CLASS_WEIGHT_ALPHA,
        max_weight_cap=CLASS_WEIGHT_CAP).to(device)
    logger.info("Calibrated class weights: %s", class_weights.tolist())

    criterion = FocalLoss(alpha=class_weights, gamma=FOCAL_GAMMA, reduction="weighted_mean")
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR)

    # Epoch 0: the untrained model, which by construction IS the PB caller.
    epoch0, _ = evaluate(model, val_features, val_prior, val_labels_flat, device)
    logger.info("epoch 0 (untrained == Poisson-binomial) | val thr=%.1f F1=%.4f "
                "P=%.4f R=%.4f AUC=%.4f AP=%.4f",
                epoch0["threshold"], epoch0["f1"], epoch0["precision"],
                epoch0["recall"], epoch0["roc_auc"], epoch0["pr_auc"])
    assert abs(epoch0["residual_gap_std"]) < 1e-9, \
        "untrained residual is not identically zero -- initialization is broken"

    # Selection metric: validation average precision. See the module docstring.
    history, best_score = [], epoch0["pr_auc"]
    payload = {"epoch": 0, "model_state_dict": model.state_dict(),
               "best_threshold": epoch0["threshold"], "val": epoch0,
               "representation": "locus_evidence_44ch", "use_gate": args.use_gate,
               "seed": args.seed}
    torch.save(payload, checkpoint_dir / f"best_{args.tag}.pt")

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss, total_task, total_penalty, nb, grad_norm = 0.0, 0.0, 0.0, 0, 0.0
        for batch in tqdm(fit_loader, desc=f"epoch {epoch}", leave=False):
            batch_features = batch["pileup_features"].to(device)
            prior = batch["prior_logits"].to(device)
            batch_labels = batch["labels"].to(device).reshape(-1)
            optimizer.zero_grad()
            logits, residual, _ = model(pileup_features=batch_features,
                                        prior_logits=prior, return_parts=True)
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
            "epoch %d | loss=%.5f grad=%.2f | VAL thr=%.1f F1=%.4f P=%.4f R=%.4f "
            "AUC=%.4f AP=%.4f | residual gap mean=%+.3f std=%.3f | gate %.3f+-%.3f",
            epoch, total_loss / max(nb, 1), grad_norm, val_metrics["threshold"],
            val_metrics["f1"], val_metrics["precision"], val_metrics["recall"],
            val_metrics["roc_auc"], val_metrics["pr_auc"],
            val_metrics["residual_gap_mean"], val_metrics["residual_gap_std"],
            val_metrics["gate_mean"], val_metrics["gate_std"])
        history.append({"epoch": epoch, "train_loss": total_loss / max(nb, 1),
                        "task_loss": total_task / max(nb, 1),
                        "residual_penalty": total_penalty / max(nb, 1),
                        "grad_norm": grad_norm, "val": val_metrics})

        payload = {"epoch": epoch, "model_state_dict": model.state_dict(),
                   "best_threshold": epoch_threshold, "val": val_metrics,
                   "representation": "locus_evidence_44ch", "use_gate": args.use_gate,
                   "drop_groups": args.drop_groups,
                   "shuffled_prior": bool(args.shuffle_prior), "seed": args.seed}
        torch.save(payload, checkpoint_dir / f"epoch{epoch}_{args.tag}.pt")
        if val_metrics["pr_auc"] > best_score:
            best_score = val_metrics["pr_auc"]
            torch.save(payload, checkpoint_dir / f"best_{args.tag}.pt")
    torch.save(payload, checkpoint_dir / f"last_{args.tag}.pt")

    best = torch.load(checkpoint_dir / f"best_{args.tag}.pt", map_location=device,
                      weights_only=False)
    logger.info("Selected epoch %d on VALIDATION AP=%.5f (PB baseline AP=%.5f) | "
                "its val F1=%.4f (PB %.4f) | frozen threshold %.1f",
                best["epoch"], best_score, epoch0["pr_auc"], best["val"]["f1"],
                epoch0["f1"], best["best_threshold"])

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps({
            "tag": args.tag, "use_gate": args.use_gate,
            "residual_lambda": args.residual_lambda,
            "drop_groups": args.drop_groups,
            "shuffled_prior": bool(args.shuffle_prior), "seed": args.seed,
            "epochs": args.epochs, "epoch0_untrained_poisson_binomial": epoch0,
            "history": history, "selected_epoch": best["epoch"],
            "selected_val_f1": best_score, "frozen_threshold": best["best_threshold"],
            "selection": "validation split only; test never consulted",
        }, indent=2, default=float))


if __name__ == "__main__":
    main()
