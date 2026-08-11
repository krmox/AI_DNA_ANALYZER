"""Inference-only evaluation of trained checkpoints on a held-out region.

Answers one question: does a checkpoint's behaviour transfer from the
validation split it was selected on to genomic sequence it has never seen?

Nothing here trains, updates weights, builds an optimizer, or writes a
checkpoint. Every metric function, the model class, the provider and the
threshold sweep are imported from the existing modules rather than
reimplemented, so a number produced here is the same number the training
loop would have produced on the same tensors.

**Threshold discipline.** The primary result for every checkpoint uses a
threshold frozen on the *original validation split* — either the value
stored in the checkpoint by ``run_train.py`` (2026-08-10 Part B onward) or,
for checkpoints written before that field existed, the value recovered by
re-running the identical ``sweep_thresholds`` call against the original
validation data. Test labels never influence it. The test-optimal
threshold is also reported, but strictly as a clearly-labelled optimistic
upper bound, never as the headline number.

Coordinate note: the held-out FASTA is its own contig starting at local 0,
exactly like the original slice. Local 0-based ``x`` maps to absolute
GRCh38 chr21 1-based ``30,000,000 + x``. That offset is documentation only
-- no file records it, and the pipeline works purely in local coordinates,
as it already did for the original slice.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

import torch
from torch.utils.data import DataLoader, Subset

from config import LABEL_NAMES, LABEL_SNP
from dna_mamba_baseline_v5 import DNAMambaBaseline
from providers import GiabAlignmentProvider, ProviderDataset
from run_train import (
    BATCH_SIZE,
    NUM_WORKERS,
    SEQ_LEN,
    TRAIN_FRACTION,
    VOCAB_SIZE,
    average_precision,
    binary_auc,
    per_class_precision_recall_f1,
    predict_with_threshold,
    torch_confusion_matrix,
)
from threshold_selection import collect_logits, sweep_thresholds

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

NUM_CLASSES = len(LABEL_NAMES)

# The original dataset, and the split run_train.py applies to it. Repeated
# here (not imported) only because run_train defines them inside main().
VAL_FASTA = "data/giab_chr21/chr21_slice.fa"
VAL_BAM = "data/giab_chr21/chr21_slice.bam"
VAL_VCF = "data/giab_chr21/chr21_slice.vcf.gz"
VAL_REGION = (1, 5_000_000)


def materialize(dataset: Any, batch_size: int) -> list[dict[str, torch.Tensor]]:
    """Collate a dataset into a list of batches, once.

    ``ProviderDataset.__getitem__`` recomputes a pysam pileup on every
    access, so iterating a fresh ``DataLoader`` per checkpoint would repay
    that cost N times over for identical tensors. Materializing once turns
    an O(checkpoints x windows) pileup cost into O(windows).

    Args:
        dataset: Dataset (or ``Subset``) of ``VariantSample``-like dicts.
        batch_size: Batch size to collate at.

    Returns:
        A list of collated batch dicts, in dataset order. Accepted anywhere
        a ``DataLoader`` is (both are just iterables of batches).
    """
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=NUM_WORKERS)
    return list(loader)


def metrics_at_threshold(
    logits: torch.Tensor, labels: torch.Tensor, threshold: float
) -> dict[str, Any]:
    """Full metric bundle for one checkpoint at one decision threshold.

    Args:
        logits: Float tensor ``[N, num_classes]``.
        labels: Long tensor ``[N]`` of true class ids.
        threshold: Non-background decision threshold to apply.

    Returns:
        Dict of SNP-class precision/recall/F1 and TP/FP/TN/FN, the
        mutation macro F1 over classes 1..3, and the full confusion matrix.
    """
    preds = predict_with_threshold(logits, threshold)
    confusion = torch_confusion_matrix(labels, preds, NUM_CLASSES)
    precision, recall, f1 = per_class_precision_recall_f1(confusion)

    # SNP-vs-rest binary view, so precision/recall/F1/TP/FP/TN/FN are
    # unambiguous scalars rather than per-class vectors.
    snp_true = labels == LABEL_SNP
    snp_pred = preds == LABEL_SNP
    tp = int((snp_true & snp_pred).sum())
    fp = int((~snp_true & snp_pred).sum())
    fn = int((snp_true & ~snp_pred).sum())
    tn = int((~snp_true & ~snp_pred).sum())

    return {
        "threshold": threshold,
        "snp_precision": float(precision[LABEL_SNP]),
        "snp_recall": float(recall[LABEL_SNP]),
        "snp_f1": float(f1[LABEL_SNP]),
        "mutation_macro_f1": float(f1[1:].mean()),
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "confusion": confusion.tolist(),
    }


def threshold_free(logits: torch.Tensor, labels: torch.Tensor) -> dict[str, float]:
    """Ranking-quality metrics that no decision threshold can affect.

    Args:
        logits: Float tensor ``[N, num_classes]``.
        labels: Long tensor ``[N]`` of true class ids.

    Returns:
        ``{"snp_auc": ..., "snp_ap": ...}``.
    """
    scores = torch.softmax(logits, dim=-1)[:, LABEL_SNP].float()
    positives = labels == LABEL_SNP
    return {"snp_auc": binary_auc(scores, positives), "snp_ap": average_precision(scores, positives)}


def build_validation_batches(batch_size: int) -> list[dict[str, torch.Tensor]]:
    """Rebuild the exact validation split ``run_train.py`` selects on.

    Mirrors ``run_train.main`` precisely: same region, same seq_len, same
    sequential (coordinate-order) 90/10 split with no shuffling. Only the
    tail Subset is materialized, so the 90% training windows never incur a
    pileup.

    Returns:
        The validation split as a list of collated batches.
    """
    provider = GiabAlignmentProvider(
        fasta_path=VAL_FASTA, bam_path=VAL_BAM, vcf_path=VAL_VCF,
        contig="chr21", region=VAL_REGION, seq_len=SEQ_LEN,
    )
    with provider:
        full = ProviderDataset(provider)
        total = len(full)
        train_len = int(total * TRAIN_FRACTION)
        logger.info(
            "Original dataset: %d windows -> train %d | val %d (val = indices %d..%d)",
            total, train_len, total - train_len, train_len, total - 1,
        )
        return materialize(Subset(full, range(train_len, total)), batch_size)


def build_test_batches(
    fasta: str, bam: str, vcf: str, region: tuple[int, int], batch_size: int,
    bed: str | None = None,
) -> list[dict[str, torch.Tensor]]:
    """Tile the held-out region into batches with the same provider settings.

    Args:
        fasta: Held-out reference FASTA.
        bam: Held-out BAM.
        vcf: Held-out VCF.
        region: ``(start, end)`` in the coordinate system of the supplied
            FASTA — local for the synthetic slices, absolute chr21 for the
            real GIAB data (whose BAM/VCF/BED are all natively absolute, so
            nothing is rebased).
        batch_size: Batch size to collate at.
        bed: Optional high-confidence BED. Supplying it is what stops
            out-of-confidence positions being scored as reliable negatives:
            ``_build_windows`` drops any window not lying entirely inside a
            confident interval. Required for real GIAB, meaningless for the
            synthetic data (which has no notion of confidence).

    Returns:
        The held-out region as a list of collated batches.
    """
    provider = GiabAlignmentProvider(
        fasta_path=fasta, bam_path=bam, vcf_path=vcf,
        contig="chr21", region=region, seq_len=SEQ_LEN,
        high_confidence_bed=bed,
    )
    with provider:
        dataset = ProviderDataset(provider)
        logger.info("Held-out dataset: %d windows", len(dataset))
        return materialize(dataset, batch_size)


def load_model(path: Path, device: torch.device) -> tuple[torch.nn.Module, dict[str, Any]]:
    """Load one checkpoint into a fresh model, in eval mode.

    Args:
        path: Checkpoint file.
        device: Device to place the model on.

    Returns:
        ``(model, metadata)`` where metadata is every checkpoint key except
        the two state dicts.
    """
    checkpoint = torch.load(path, map_location=device, weights_only=False)
    model = DNAMambaBaseline(num_classes=NUM_CLASSES, vocab_size=VOCAB_SIZE).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    metadata = {
        key: value for key, value in checkpoint.items()
        if key not in ("model_state_dict", "optimizer_state_dict")
    }
    return model, metadata


def parse_args() -> argparse.Namespace:
    """CLI for the held-out evaluation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoints", nargs="+", required=True, help="Checkpoint paths.")
    parser.add_argument("--test-fasta", required=True)
    parser.add_argument("--test-bam", required=True)
    parser.add_argument("--test-vcf", required=True)
    parser.add_argument("--test-region-start", type=int, default=1)
    parser.add_argument("--test-region-end", type=int, required=True)
    parser.add_argument(
        "--test-bed", default=None,
        help="High-confidence BED restricting evaluated windows. Use for real GIAB.",
    )
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    parser.add_argument("--out", default="", help="Optional JSON results path.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info("Device: %s | inference only, no weights will be updated", device)

    val_batches = build_validation_batches(args.batch_size)
    test_batches = build_test_batches(
        args.test_fasta, args.test_bam, args.test_vcf,
        (args.test_region_start, args.test_region_end), args.batch_size, args.test_bed,
    )

    results = []
    for path in sorted(args.checkpoints):
        checkpoint_path = Path(path)
        model, metadata = load_model(checkpoint_path, device)

        # --- threshold frozen on ORIGINAL VALIDATION data only ---
        val_logits, val_labels = collect_logits(model, val_batches, device)
        val_sweep = sweep_thresholds(val_logits, val_labels)
        frozen_threshold = val_sweep.best_threshold
        stored = metadata.get("best_threshold")
        if stored is not None and abs(float(stored) - frozen_threshold) > 1e-9:
            logger.warning(
                "%s: recomputed validation threshold %.2f != stored %.2f",
                checkpoint_path.name, frozen_threshold, float(stored),
            )

        val_metrics = metrics_at_threshold(val_logits, val_labels, frozen_threshold)
        val_metrics.update(threshold_free(val_logits, val_labels))

        # --- held-out test, frozen threshold applied ---
        test_logits, test_labels = collect_logits(model, test_batches, device)
        test_metrics = metrics_at_threshold(test_logits, test_labels, frozen_threshold)
        test_metrics.update(threshold_free(test_logits, test_labels))

        # Diagnostic only. Selected ON TEST LABELS, therefore an optimistic
        # upper bound and never the headline result.
        test_sweep = sweep_thresholds(test_logits, test_labels)
        diagnostic = {
            "test_optimal_threshold": test_sweep.best_threshold,
            "test_mutation_macro_f1_at_test_optimal": test_sweep.f1_at_best,
        }

        results.append({
            "checkpoint": checkpoint_path.name,
            "metadata": {k: v for k, v in metadata.items() if not isinstance(v, torch.Tensor)},
            "frozen_validation_threshold": frozen_threshold,
            "stored_threshold": stored,
            "validation": val_metrics,
            "test": test_metrics,
            "diagnostic_upper_bound": diagnostic,
        })
        logger.info(
            "%-46s thr=%.2f | VAL auc=%.4f ap=%.4f f1=%.4f | TEST auc=%.4f ap=%.4f f1=%.4f",
            checkpoint_path.name, frozen_threshold,
            val_metrics["snp_auc"], val_metrics["snp_ap"], val_metrics["mutation_macro_f1"],
            test_metrics["snp_auc"], test_metrics["snp_ap"], test_metrics["mutation_macro_f1"],
        )

    if args.out:
        Path(args.out).write_text(json.dumps(results, indent=2))
        logger.info("Wrote %s", args.out)


if __name__ == "__main__":
    main()
