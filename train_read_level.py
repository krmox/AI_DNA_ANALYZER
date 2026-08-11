"""Train the read-level model on cached read tensors. One controlled variable.

Every training hyperparameter is imported from ``train_raw_pileup`` rather than
restated, so the 14-channel aggregate run and this one cannot silently diverge:
same focal gamma, class-weight calibration, LR, batch size, seq_len, sequential
90/10 coordinate-order split, mutation-presence sampler, epoch count, and the
same probability-space validation threshold sweep the aggregate model used.

The changed variable is the input representation: individual read observations
instead of per-locus summary statistics.

Threshold protocol is deliberately the *aggregate model's*, not the binomial's:
this model is a softmax classifier like ``RawPileupMamba``, so it is selected by
``sweep_thresholds`` on validation mutation_macro_f1 exactly as method D was.
Keeping the neural protocol fixed is what makes "read-level vs aggregate" a
one-variable comparison.

``--control`` and ``--drop-groups`` switch on the ablations. They transform the
cached tensor (or zero feature columns) and change nothing else -- not the
model, not the optimizer, not the tensor shape -- so an ablation measures
information rather than capacity.
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Subset, WeightedRandomSampler
from tqdm import tqdm

from config import LABEL_NORMAL
from loss import FocalLoss, compute_calibrated_class_weights
from model_read_level import ReadLevelMamba
from read_level_controls import apply_control
from read_level_pileup import FEATURE_GROUPS, build_features_torch, resolve_drop_indices
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

#: Bit 0 of the flags column; a row is a real read when set.
FLAG_VALID_BIT = 1


def to_windows(reads: np.ndarray) -> torch.Tensor:
    """Reshape flat loci into ``[W, SEQ_LEN, R, C]`` uint8 windows."""
    tensor = torch.from_numpy(np.ascontiguousarray(reads))
    return tensor.reshape(-1, SEQ_LEN, reads.shape[-2], reads.shape[-1])


def expand(batch_reads: torch.Tensor, drop_indices: tuple[int, ...],
           device: torch.device) -> tuple[torch.Tensor, torch.Tensor]:
    """Move one uint8 batch to the device and expand it to float features."""
    raw = batch_reads.to(device, non_blocking=True)
    valid = (raw[..., 6].to(torch.int64) & FLAG_VALID_BIT) > 0
    return build_features_torch(raw, drop_indices), valid


def reference_onehot(reference_index: np.ndarray) -> torch.Tensor:
    """One-hot the locus reference base over (N, A, C, G, T).

    Taken from the aggregate cache's ``reference_index`` column, which comes
    from the FASTA. The sanity gate proves the two caches are locus-aligned, so
    this cannot silently pair a locus with another locus' reference.
    """
    index = torch.tensor(np.asarray(reference_index, dtype=np.int64)).clamp(0, 4)
    return torch.eye(5, dtype=torch.float32)[index].reshape(-1, SEQ_LEN, 5)


@torch.no_grad()
def collect_logits(model, windows: torch.Tensor, reference: torch.Tensor,
                   labels: torch.Tensor, drop_indices: tuple[int, ...],
                   device: torch.device, batch_size: int = 64):
    """Inference over a split held as uint8 windows."""
    model.eval()
    logits_out = []
    for begin in range(0, windows.shape[0], batch_size):
        features, valid = expand(windows[begin:begin + batch_size], drop_indices, device)
        logits = model(read_features=features, read_valid=valid,
                       reference_onehot=reference[begin:begin + batch_size].to(device))
        logits_out.append(logits.reshape(-1, logits.size(-1)).cpu())
    return torch.cat(logits_out), labels.reshape(-1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-reads", required=True)
    parser.add_argument("--train-counts", required=True,
                        help="Aggregate cache supplying the locus reference base; "
                             "locus alignment is verified by sanity_read_level.py.")
    parser.add_argument("--epochs", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--control", default="",
                        choices=["", "permute_reads", "decouple", "aggregate_only"],
                        help="Tensor-level control applied to fit, val and test alike.")
    parser.add_argument("--control-seed", type=int, default=4321)
    parser.add_argument("--drop-groups", nargs="*", default=[],
                        choices=sorted(FEATURE_GROUPS), help="Feature ablations.")
    parser.add_argument("--tag", required=True)
    parser.add_argument("--checkpoint-dir", default="checkpoints_readlevel")
    parser.add_argument("--out", default="")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint_dir = Path(args.checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    blob = np.load(args.train_reads, allow_pickle=True)
    reads, labels = blob["reads"], blob["labels"]
    count_blob = np.load(args.train_counts, allow_pickle=True)
    if not np.array_equal(count_blob["labels"], labels):
        raise SystemExit("read cache and count cache are not locus-aligned")
    reference_index = count_blob["counts"][:, 9]

    if args.control:
        start = time.perf_counter()
        reads = apply_control(reads, args.control, args.control_seed)
        logger.info("CONTROL '%s' applied to the training tensor in %.1fs",
                    args.control, time.perf_counter() - start)
    drop_indices = resolve_drop_indices(tuple(args.drop_groups))
    if drop_indices:
        logger.info("ABLATION: zeroing feature groups %s -> indices %s",
                    args.drop_groups, drop_indices)

    # Sequential split by batch, byte-identical to the previous experiments.
    n_windows = labels.size // SEQ_LEN
    n_batches = (n_windows + BATCH_SIZE - 1) // BATCH_SIZE
    split_locus = int(n_batches * TRAIN_FRACTION) * BATCH_SIZE * SEQ_LEN
    logger.info("Training region: %d loci / %d windows; fit=%d loci, val=%d loci",
                labels.size, n_windows, split_locus, labels.size - split_locus)

    fit_windows = to_windows(reads[:split_locus])
    val_windows = to_windows(reads[split_locus:])
    fit_reference = reference_onehot(reference_index[:split_locus])
    val_reference = reference_onehot(reference_index[split_locus:])
    fit_labels = torch.tensor(labels[:split_locus], dtype=torch.long)
    val_labels = torch.tensor(labels[split_locus:], dtype=torch.long)
    fit_window_labels = fit_labels.reshape(-1, SEQ_LEN)

    # Mutation-presence sampler, same construction as the previous experiments.
    has_mutation = (fit_window_labels != LABEL_NORMAL).any(dim=1)
    n_mut = int(has_mutation.sum())
    n_bg = int(has_mutation.numel() - n_mut)
    logger.info("Fit windows with a mutation: %d | pure background: %d", n_mut, n_bg)
    weights = torch.where(has_mutation, 1.0 / max(n_mut, 1), 1.0 / max(n_bg, 1))

    window_ds = [{"reads": fit_windows[i], "reference": fit_reference[i],
                  "labels": fit_window_labels[i]}
                 for i in range(fit_windows.shape[0])]
    sampler = WeightedRandomSampler(weights, num_samples=len(window_ds), replacement=True)
    fit_loader = DataLoader(window_ds, batch_size=BATCH_SIZE, sampler=sampler)

    model = ReadLevelMamba().to(device)
    n_parameters = sum(p.numel() for p in model.parameters())
    logger.info("Model parameters: %d | control=%s | drop=%s | tag=%s",
                n_parameters, args.control or "none", args.drop_groups or "none", args.tag)

    # Class weights are calibrated from labels only, exactly as before; the
    # calibration helper reads the "labels" key and ignores the rest.
    class_weights = compute_calibrated_class_weights(
        Subset(window_ds, range(min(200, len(window_ds)))), num_classes=NUM_CLASSES,
        sample_size=min(200, len(window_ds)), alpha=CLASS_WEIGHT_ALPHA,
        max_weight_cap=CLASS_WEIGHT_CAP).to(device)
    logger.info("Calibrated class weights: %s", class_weights.tolist())

    criterion = FocalLoss(alpha=class_weights, gamma=FOCAL_GAMMA, reduction="weighted_mean")
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR)

    history, best_score, payload = [], float("-inf"), None
    for epoch in range(1, args.epochs + 1):
        model.train()
        total, nb, grad_norm = 0.0, 0, 0.0
        epoch_start = time.perf_counter()
        for batch in tqdm(fit_loader, desc=f"epoch {epoch}", leave=False):
            features, valid = expand(batch["reads"], drop_indices, device)
            batch_labels = batch["labels"].to(device).reshape(-1)
            optimizer.zero_grad()
            logits = model(read_features=features, read_valid=valid,
                           reference_onehot=batch["reference"].to(device))
            loss = criterion(logits.reshape(-1, NUM_CLASSES), batch_labels)
            loss.backward()
            grad_norm = float(sum(p.grad.norm() ** 2 for p in model.parameters()
                                  if p.grad is not None) ** 0.5)
            optimizer.step()
            total += loss.item()
            nb += 1
        train_seconds = time.perf_counter() - epoch_start

        val_logits, val_labels_flat = collect_logits(
            model, val_windows, val_reference, val_labels, drop_indices, device)
        sweep = sweep_thresholds(val_logits, val_labels_flat)
        logger.info("epoch %d | train_loss=%.5f grad_norm=%.4f %.0fs | val f1@opt=%.4f thr=%.2f",
                    epoch, total / max(nb, 1), grad_norm, train_seconds,
                    sweep.f1_at_best, sweep.best_threshold)
        val_metrics = report(val_logits, val_labels_flat, sweep.best_threshold,
                             f"  VAL epoch{epoch}")
        history.append({"epoch": epoch, "train_loss": total / max(nb, 1),
                        "grad_norm": grad_norm, "train_seconds": train_seconds,
                        "val": val_metrics})

        payload = {"epoch": epoch, "model_state_dict": model.state_dict(),
                   "best_threshold": sweep.best_threshold, "val": val_metrics,
                   "representation": "read_level_v1", "control": args.control,
                   "drop_groups": list(args.drop_groups), "seed": args.seed,
                   "parameters": n_parameters}
        torch.save(payload, checkpoint_dir / f"epoch{epoch}_{args.tag}.pt")
        if sweep.f1_at_best > best_score:
            best_score = sweep.f1_at_best
            torch.save(payload, checkpoint_dir / f"best_{args.tag}.pt")
    if payload is not None:
        torch.save(payload, checkpoint_dir / f"last_{args.tag}.pt")

    best = torch.load(checkpoint_dir / f"best_{args.tag}.pt", map_location=device,
                      weights_only=False)
    logger.info("Selected epoch %d on VALIDATION | val macro-F1=%.4f | frozen threshold %.2f",
                best["epoch"], best_score, best["best_threshold"])

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps({
            "tag": args.tag, "control": args.control, "drop_groups": list(args.drop_groups),
            "seed": args.seed, "epochs": args.epochs, "parameters": n_parameters,
            "history": history, "selected_epoch": best["epoch"],
            "selected_val_macro_f1": best_score,
            "frozen_threshold": best["best_threshold"],
            "selection": "validation split only; test never consulted",
        }, indent=2))


if __name__ == "__main__":
    main()
