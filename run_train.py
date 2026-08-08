"""Production training loop for the reference-aware DNA-Mamba variant caller.

Trains on GIAB-derived alignment windows with a strictly sequential
train/validation split (windows are tiled along the genome in coordinate
order, so a random split would leak neighbouring, highly-correlated columns
across the boundary). The model's forward signature is inspected at runtime
so this loop keeps working if the model starts accepting more of the
provider's features (``base_quality``, ``depth``, ``vaf``) without any edit
here.

Before touching loss/sampler knobs, read ``_warn_if_bam_lacks_indel_evidence``
and the module docstring below it. On this project's current
``data/giab_chr21/chr21_slice.bam`` (produced by ``generate_mock_bam.py``),
every read's CIGAR is a flat ``(0, 150)`` match — the mock read generator
only ever injects the SNP branch of the VCF, never the indel branch. That
means every Insertion/Deletion-labelled column in the training data has
*zero token-level evidence*: ``input_ids`` and ``reference_ids`` at that
locus are identical, and ``vaf``/``depth`` show nothing unusual either. No
loss reweighting, sampler, threshold, or scheduler change can make classes
2/3 learnable while that holds — see the confirmed measurements in the
diagnosis this file's commit message / PR description links back to.
"""

from __future__ import annotations

import inspect
import logging
from pathlib import Path
from typing import Any, NamedTuple

import torch
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.utils.data import DataLoader, Dataset, Subset, WeightedRandomSampler
from tqdm import tqdm

from config import LABEL_NAMES, LABEL_NORMAL, LABEL_SNP, VOCAB
from dna_mamba_baseline_v5 import DNAMambaBaseline
from loss import FocalLoss, compute_calibrated_class_weights
from providers import GiabAlignmentProvider, ProviderDataset

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

FASTA_PATH = "data/giab_chr21/chr21_slice.fa"
VCF_PATH = "data/giab_chr21/chr21_slice.vcf.gz"
BAM_PATH = "data/giab_chr21/chr21_slice.bam"

CHECKPOINT_DIR = Path("checkpoints")
BEST_MODEL_PATH = CHECKPOINT_DIR / "best_model.pt"
LAST_CHECKPOINT_PATH = CHECKPOINT_DIR / "last_checkpoint.pt"

MAX_EPOCHS = 20
BATCH_SIZE = 16
LR = 3e-4
TRAIN_FRACTION = 0.9
NUM_WORKERS = 0
SEQ_LEN = 64
NUM_CLASSES = len(LABEL_NAMES)
VOCAB_SIZE = len(VOCAB)

# Focal-loss dynamics. gamma=3.0 pushes harder on hard-to-classify tokens
# than the gamma=2.0 default; "weighted_mean" reduction (see loss.py)
# divides by the sum of alpha weights rather than raw token count, so the
# ~99.9%-background distribution measured on this dataset (see
# summarise_label_distribution.py) can't dilute the minority-class gradient
# down to noise the way plain "mean" reduction did.
FOCAL_GAMMA = 3.0
CLASS_WEIGHT_ALPHA = 0.75
CLASS_WEIGHT_CAP = 20.0

# Non-background decision threshold used only for validation metrics: a
# locus is called as class 1/2/3 only if that class's softmax probability
# clears this bar, otherwise it falls back to class 0. Training itself still
# optimizes the raw per-token loss; only reported metrics use this.
VARIANT_PROB_THRESHOLD = 0.5

EARLY_STOPPING_PATIENCE = 5
EARLY_STOPPING_MIN_DELTA = 1e-4


def _warn_if_bam_lacks_indel_evidence(bam_path: str | Path, sample_reads: int = 5000) -> None:
    """Flag a BAM that cannot possibly support Insertion/Deletion learning.

    Scans up to ``sample_reads`` alignments for a CIGAR operation other than
    a plain match (op codes 1 = insertion, 2 = deletion). A BAM assembled by
    always emitting ``cigar=((0, read_len),)`` — which is what
    ``generate_mock_bam.py`` currently does for every read — never contains
    this evidence, so any Insertion/Deletion label in the paired VCF is
    unlearnable from the pileup no matter how the loss is tuned.

    Args:
        bam_path: Path to the BAM to inspect.
        sample_reads: Maximum alignments to scan before giving up.
    """
    import pysam

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    checked = 0
    found_indel_cigar = False
    try:
        for read in bam.fetch(until_eof=True):
            checked += 1
            if read.cigartuples and any(op in (1, 2) for op, _ in read.cigartuples):
                found_indel_cigar = True
                break
            if checked >= sample_reads:
                break
    finally:
        bam.close()

    if not found_indel_cigar:
        logger.warning(
            "No insertion/deletion CIGAR operation found in the first %d reads of %s. "
            "Insertion/Deletion labels cannot be learned from a BAM that never encodes "
            "an indel in any read — this is a data-generation bug, not a loss/imbalance "
            "problem. Fix the read generator (or use real aligned reads) before expecting "
            "Class 2/3 F1 above 0.",
            checked, bam_path,
        )


def filter_batch_for_model(model: torch.nn.Module, batch: dict[str, Any]) -> dict[str, Any]:
    """Keep only the batch keys that ``model.forward`` actually declares.

    Never widened into a blanket ``**batch`` unpack: the provider's batches
    carry fields (``base_quality``, ``depth``, ``vaf``, ``labels``) that
    ``DNAMambaBaseline.forward`` does not accept today, and a blind unpack
    would raise ``TypeError`` the moment any of those show up.

    Args:
        model: The model whose forward signature to inspect.
        batch: A collated batch dict, as produced by the ``DataLoader``.

    Returns:
        The subset of ``batch`` whose keys match a parameter of
        ``model.forward`` (excluding ``self``).
    """
    forward_params = inspect.signature(model.forward).parameters
    return {key: value for key, value in batch.items() if key in forward_params}


def move_batch_to_device(batch: dict[str, Any], device: torch.device) -> dict[str, Any]:
    """Move every tensor value in a batch dict to ``device``."""
    return {key: value.to(device) for key, value in batch.items()}


def predict_with_threshold(logits: torch.Tensor, threshold: float) -> torch.Tensor:
    """Call a non-background class only when its probability clears ``threshold``.

    Plain ``argmax`` lets class 0 win every locus where it is merely the
    *most* likely class, which on a ~99.9%-background distribution is
    nearly everywhere. Requiring the winning variant class to clear an
    absolute probability bar (falling back to class 0 otherwise) makes the
    decision boundary explicit and tunable independent of the loss.

    Args:
        logits: Float tensor ``[N, num_classes]`` of raw model outputs.
        threshold: Minimum softmax probability a non-background class needs
            to be predicted; below it, the prediction is class 0.

    Returns:
        Long tensor ``[N]`` of predicted class ids.
    """
    probs = torch.softmax(logits, dim=-1)
    variant_probs = probs[:, 1:]
    best_variant_prob, best_variant_class = variant_probs.max(dim=-1)
    predictions = torch.where(
        best_variant_prob >= threshold,
        best_variant_class + 1,
        torch.zeros_like(best_variant_class),
    )
    return predictions


def binary_auc(scores: torch.Tensor, positives: torch.Tensor) -> float:
    """Rank-based AUC (Mann-Whitney U statistic, normalised) for one class.

    Reported alongside — never instead of — the F1@``VARIANT_PROB_THRESHOLD``
    table, because the two answer different questions. F1 at a fixed absolute
    probability bar conflates *ranking quality* with *where the decision
    boundary happens to sit*, and on this dataset the model's ``p(SNP)`` on
    background tokens sits in a narrow band just under 0.5, so a drift of a
    few hundredths moves thousands of tokens across the bar and craters
    precision while the ranking is untouched. AUC is threshold-free and
    isolates the ranking half of that.

    Ties contribute 0.5, via average ranks within each tied group.

    Args:
        scores: Float tensor ``[N]`` of per-token scores for the class.
        positives: Bool tensor ``[N]``, True where the token's true label is
            the class in question.

    Returns:
        AUC in ``[0, 1]``, or ``float("nan")`` if either group is empty
        (undefined rather than silently 0 — a 0 here would read as
        "perfectly wrong ranking", which is a very different claim).
    """
    num_positive = int(positives.sum().item())
    num_negative = int(positives.numel() - num_positive)
    if num_positive == 0 or num_negative == 0:
        return float("nan")

    scores = scores.double()
    order = scores.argsort()
    sorted_scores = scores[order]
    _, inverse, counts = torch.unique(sorted_scores, return_inverse=True, return_counts=True)
    positions = torch.arange(1, scores.numel() + 1, dtype=torch.float64, device=scores.device)
    tie_group_sums = torch.zeros(counts.numel(), dtype=torch.float64, device=scores.device)
    tie_group_sums.index_add_(0, inverse, positions)
    mean_rank_per_group = tie_group_sums / counts.double()

    ranks = torch.empty_like(scores)
    ranks[order] = mean_rank_per_group[inverse]

    positive_rank_sum = ranks[positives].sum()
    u_statistic = positive_rank_sum - num_positive * (num_positive + 1) / 2
    return float((u_statistic / (num_positive * num_negative)).item())


def average_precision(scores: torch.Tensor, positives: torch.Tensor) -> float:
    """Area under the precision-recall curve for one class (step interpolation).

    The companion to :func:`binary_auc`, and the more informative of the two
    at this class balance: with ~374 positives against ~474k background
    tokens, AUC stays high for a model that still emits far more false
    positives than true ones, because the negatives it ranks correctly
    overwhelm the few it does not. AP is computed against the positive class
    directly and does not flatter that regime.

    Args:
        scores: Float tensor ``[N]`` of per-token scores for the class.
        positives: Bool tensor ``[N]``, True where the token's true label is
            the class in question.

    Returns:
        Average precision in ``[0, 1]``, or ``float("nan")`` if there are no
        positives.
    """
    num_positive = int(positives.sum().item())
    if num_positive == 0:
        return float("nan")

    order = scores.double().argsort(descending=True)
    positive_sorted = positives[order].double()
    true_positives = torch.cumsum(positive_sorted, dim=0)
    false_positives = torch.cumsum(1.0 - positive_sorted, dim=0)
    precision_at_k = true_positives / (true_positives + false_positives).clamp_min(1e-12)
    return float(((precision_at_k * positive_sorted).sum() / num_positive).item())


def torch_confusion_matrix(
    y_true: torch.Tensor, y_pred: torch.Tensor, num_classes: int
) -> torch.Tensor:
    """Compute a confusion matrix entirely with torch ops (no sklearn).

    Args:
        y_true: Long tensor ``[N]`` of ground-truth class ids.
        y_pred: Long tensor ``[N]`` of predicted class ids.
        num_classes: Number of classes.

    Returns:
        Long tensor ``[num_classes, num_classes]``; entry ``[i, j]`` is the
        count of true-class ``i`` tokens predicted as class ``j``.
    """
    indices = y_true * num_classes + y_pred
    return torch.bincount(indices, minlength=num_classes * num_classes).reshape(
        num_classes, num_classes
    )


def per_class_precision_recall_f1(
    confusion: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Derive per-class precision/recall/F1 from a confusion matrix.

    Args:
        confusion: ``[num_classes, num_classes]`` matrix as produced by
            :func:`torch_confusion_matrix`.

    Returns:
        Three float tensors ``[num_classes]``: precision, recall, F1.
        Classes with zero predicted (precision) or zero true (recall)
        instances get 0 rather than NaN.
    """
    confusion = confusion.float()
    true_positive = confusion.diag()
    predicted_totals = confusion.sum(dim=0)
    true_totals = confusion.sum(dim=1)

    precision = torch.where(
        predicted_totals > 0, true_positive / predicted_totals, torch.zeros_like(true_positive)
    )
    recall = torch.where(
        true_totals > 0, true_positive / true_totals, torch.zeros_like(true_positive)
    )
    denom = precision + recall
    f1 = torch.where(denom > 0, 2 * precision * recall / denom, torch.zeros_like(true_positive))
    return precision, recall, f1


def log_confusion_and_per_class_metrics(confusion: torch.Tensor) -> None:
    """Log a confusion matrix and per-class precision/recall/F1.

    Args:
        confusion: ``[num_classes, num_classes]`` matrix, rows=true, cols=pred.
    """
    matrix = confusion.tolist()
    header = "        " + "".join(f"{name:>12}" for name in LABEL_NAMES)
    rows = [
        f"{name:>8}" + "".join(f"{count:>12d}" for count in row)
        for name, row in zip(LABEL_NAMES, matrix)
    ]
    logger.info("Confusion matrix (rows=true, cols=pred):\n%s\n%s", header, "\n".join(rows))

    precision, recall, f1 = per_class_precision_recall_f1(confusion)
    support = confusion.sum(dim=1)
    for idx in range(1, NUM_CLASSES):
        logger.info(
            "  Class %d (%-9s) precision=%.4f recall=%.4f f1=%.4f support=%d",
            idx, LABEL_NAMES[idx], precision[idx].item(), recall[idx].item(),
            f1[idx].item(), int(support[idx].item()),
        )


def compute_mutation_presence_weights(dataset: Dataset, num_classes: int) -> torch.Tensor:
    """Build per-window sampler weights from a single pass over labels only.

    Each window gets one of two weights depending on whether it contains at
    least one non-background (SNP/Insertion/Deletion) token: ``1 /
    count_in_its_group``. That makes a :class:`WeightedRandomSampler` draw
    roughly half its samples from "has a mutation" windows and half from
    pure-background windows in expectation, regardless of how rare mutation
    windows actually are — the true fix for gradient starvation at ~99.9%
    background density, as opposed to per-token alpha weights alone, which
    still only get to act on however many mutation windows a random batch
    happens to contain.

    This never materializes the dataset: it iterates ``dataset[i]`` once,
    reads out ``"labels"`` into a Python bool, and lets the rest of the
    sample (``input_ids``, ``reference_ids``, ``depth``, ...) get garbage
    collected immediately. Peak extra memory is one ``bool`` tensor of
    length ``len(dataset)``.

    Args:
        dataset: A dataset (or ``Subset``) of ``VariantSample``-like dicts.
        num_classes: Unused directly but kept for signature symmetry with
            the rest of this module's metric helpers.

    Returns:
        Float tensor ``[len(dataset)]`` of per-window sampler weights.
    """
    del num_classes  # Not needed — presence is binary regardless of class count.
    has_mutation = torch.zeros(len(dataset), dtype=torch.bool)
    for index in range(len(dataset)):
        labels = dataset[index]["labels"]
        has_mutation[index] = bool((labels != LABEL_NORMAL).any())

    mutation_count = int(has_mutation.sum().item())
    background_count = len(dataset) - mutation_count
    logger.info(
        "Sampler weighting -> windows with a mutation: %d | pure background: %d",
        mutation_count, background_count,
    )

    weights = torch.empty(len(dataset), dtype=torch.float)
    weights[has_mutation] = 1.0 / max(mutation_count, 1)
    weights[~has_mutation] = 1.0 / max(background_count, 1)
    return weights


class EpochMetrics(NamedTuple):
    """Summary of one training or validation pass.

    ``snp_auc``/``snp_ap`` are threshold-free diagnostics appended after the
    2026-08-08 investigation; they are *reported only*. Checkpoint selection
    and early stopping still key off ``mutation_macro_f1`` (measured at
    ``VARIANT_PROB_THRESHOLD``) exactly as before, so scores stay comparable
    with every run already recorded in ``History/3_DEVLOG.md``.
    """

    loss: float
    macro_f1: float
    mutation_macro_f1: float
    confusion: torch.Tensor
    snp_auc: float = float("nan")
    snp_ap: float = float("nan")


def run_epoch(
    model: torch.nn.Module,
    loader: DataLoader,
    criterion: torch.nn.Module,
    device: torch.device,
    optimizer: torch.optim.Optimizer | None,
    epoch: int,
    total_epochs: int,
    variant_prob_threshold: float | None = None,
    report_diagnostics: bool = False,
) -> EpochMetrics:
    """Run one training or validation pass over ``loader``.

    Training vs. validation mode is entirely determined by ``optimizer``:
    pass an optimizer for a training pass (gradients flow, weights update),
    pass ``None`` for a validation pass (``model.eval()`` + ``no_grad``).

    Macro F1 is computed two ways. ``macro_f1`` averages over all 4 classes
    including Normal, whose F1 sits near 0.99 by construction at ~99.9%
    prevalence — that alone floors the metric around 0.25 even with zero
    variant-calling skill, so it reads as "mostly fine" when it is not.
    ``mutation_macro_f1`` averages only classes 1-3 and is what actually
    tracks variant-calling quality; it is what checkpointing and early
    stopping key off.

    Returns:
        :class:`EpochMetrics` for this pass.
    """
    is_train = optimizer is not None
    model.train(mode=is_train)

    phase = "Train" if is_train else "Val  "
    running_loss = 0.0
    num_batches = 0
    confusion = torch.zeros(NUM_CLASSES, NUM_CLASSES, dtype=torch.long, device=device)

    # Buffered on CPU for the threshold-free AUC/AP report. Both metrics are
    # rank-based, so they cannot be accumulated batch-wise the way the
    # confusion matrix can — the scores have to be ranked against each other
    # globally. One float + one bool per token (~5 MB per validation pass at
    # the current split); moved off-device immediately so this never competes
    # with the model for GPU memory.
    snp_scores: list[torch.Tensor] = []
    snp_positives: list[torch.Tensor] = []

    progress = tqdm(loader, desc=f"Epoch {epoch}/{total_epochs} [{phase}]", leave=False)

    context = torch.enable_grad() if is_train else torch.no_grad()
    with context:
        for batch in progress:
            batch = move_batch_to_device(batch, device)
            labels = batch["labels"].reshape(-1)

            if is_train:
                optimizer.zero_grad()

            model_inputs = filter_batch_for_model(model, batch)
            logits = model(**model_inputs)
            flat_logits = logits.reshape(-1, logits.size(-1))

            loss = criterion(flat_logits, labels)

            if is_train:
                loss.backward()
                optimizer.step()

            running_loss += loss.item()
            num_batches += 1

            with torch.no_grad():
                detached_logits = flat_logits.detach()
                if variant_prob_threshold is not None:
                    preds = predict_with_threshold(detached_logits, variant_prob_threshold)
                else:
                    preds = detached_logits.argmax(dim=-1)
                confusion += torch_confusion_matrix(labels, preds, NUM_CLASSES)

                # p(SNP) per token, independent of any decision threshold.
                snp_scores.append(
                    torch.softmax(detached_logits, dim=-1)[:, LABEL_SNP].float().cpu()
                )
                snp_positives.append((labels == LABEL_SNP).cpu())

            progress.set_postfix(loss=f"{running_loss / num_batches:.4f}")

    mean_loss = running_loss / max(num_batches, 1)
    _, _, f1 = per_class_precision_recall_f1(confusion)
    macro_f1 = float(f1.mean().item())
    mutation_macro_f1 = float(f1[1:].mean().item())

    if snp_scores:
        all_scores = torch.cat(snp_scores)
        all_positives = torch.cat(snp_positives)
        snp_auc = binary_auc(all_scores, all_positives)
        snp_ap = average_precision(all_scores, all_positives)
    else:
        snp_auc = snp_ap = float("nan")

    if report_diagnostics:
        log_confusion_and_per_class_metrics(confusion)

    return EpochMetrics(mean_loss, macro_f1, mutation_macro_f1, confusion, snp_auc, snp_ap)


def main() -> None:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info("Using device: %s", device)

    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    _warn_if_bam_lacks_indel_evidence(BAM_PATH)

    provider = GiabAlignmentProvider(
        fasta_path=FASTA_PATH,
        bam_path=BAM_PATH,
        vcf_path=VCF_PATH,
        contig="chr21",
        region=(1, 5_000_000),
        seq_len=SEQ_LEN,
    )

    with provider:
        full_dataset = ProviderDataset(provider)
        total_len = len(full_dataset)
        if total_len == 0:
            raise RuntimeError("Provider yielded zero windows; check region/contig/paths.")

        # Sequential (coordinate-order) split: no shuffling. Random splitting
        # here would leak overlapping/adjacent genomic context between train
        # and validation, inflating validation metrics.
        train_len = int(total_len * TRAIN_FRACTION)
        train_ds = Subset(full_dataset, range(0, train_len))
        val_ds = Subset(full_dataset, range(train_len, total_len))
        logger.info("Samples -> train: %d | val: %d", len(train_ds), len(val_ds))

        sampler_weights = compute_mutation_presence_weights(train_ds, NUM_CLASSES)
        train_sampler = WeightedRandomSampler(
            sampler_weights, num_samples=len(train_ds), replacement=True
        )
        train_loader = DataLoader(
            train_ds, batch_size=BATCH_SIZE, sampler=train_sampler, num_workers=NUM_WORKERS
        )
        val_loader = DataLoader(
            val_ds, batch_size=BATCH_SIZE, shuffle=False, num_workers=NUM_WORKERS
        )

        # vocab_size must be explicit, not the class default: the tokenizer's
        # real vocabulary (config.VOCAB) includes the gap token id (5), which
        # only started appearing in real windows once generate_mock_bam.py's
        # indel splicing fix landed. An undersized embedding table here is a
        # silent latent bug until then and a CUDA device-side assert after.
        model = DNAMambaBaseline(num_classes=NUM_CLASSES, vocab_size=VOCAB_SIZE).to(device)

        class_weights = compute_calibrated_class_weights(
            train_ds,
            num_classes=NUM_CLASSES,
            sample_size=min(200, len(train_ds)),
            alpha=CLASS_WEIGHT_ALPHA,
            max_weight_cap=CLASS_WEIGHT_CAP,
        ).to(device)
        logger.info("Calibrated class weights: %s", class_weights.tolist())

        # alpha is a fixed buffer, not an nn.Parameter: it is a data-derived
        # prior (measured class frequency), and letting the optimizer touch
        # it would hand it a trivial way to shrink the rare-class loss
        # contribution back toward zero instead of actually learning those
        # classes — exactly the collapse this loss is meant to prevent.
        criterion = FocalLoss(alpha=class_weights, gamma=FOCAL_GAMMA, reduction="weighted_mean")
        optimizer = torch.optim.AdamW(model.parameters(), lr=LR)
        scheduler = ReduceLROnPlateau(optimizer, mode="max", factor=0.5, patience=1)

        best_score = float("-inf")
        best_val_loss = float("inf")
        epochs_without_improvement = 0

        for epoch in range(1, MAX_EPOCHS + 1):
            train_metrics = run_epoch(
                model, train_loader, criterion, device, optimizer, epoch, MAX_EPOCHS
            )
            logger.info(
                "Epoch %d/%d [Train] loss=%.4f macro_f1=%.4f mutation_macro_f1=%.4f",
                epoch, MAX_EPOCHS, train_metrics.loss, train_metrics.macro_f1,
                train_metrics.mutation_macro_f1,
            )

            val_metrics = run_epoch(
                model,
                val_loader,
                criterion,
                device,
                None,
                epoch,
                MAX_EPOCHS,
                variant_prob_threshold=VARIANT_PROB_THRESHOLD,
                report_diagnostics=True,
            )
            logger.info(
                "Epoch %d/%d [Val]   loss=%.4f macro_f1=%.4f mutation_macro_f1=%.4f",
                epoch, MAX_EPOCHS, val_metrics.loss, val_metrics.macro_f1,
                val_metrics.mutation_macro_f1,
            )
            # Threshold-free companions to the F1@VARIANT_PROB_THRESHOLD table
            # above. A large gap between these two lines — AUC/AP flat while
            # mutation_macro_f1 falls — means the ranking is intact and only
            # the fixed decision bar has gone stale, not that the model got
            # worse. See History/3_DEVLOG.md, 2026-08-08.
            logger.info(
                "Epoch %d/%d [Val]   SNP-vs-rest (threshold-free): auc=%.4f ap=%.4f",
                epoch, MAX_EPOCHS, val_metrics.snp_auc, val_metrics.snp_ap,
            )

            current_lr = optimizer.param_groups[0]["lr"]
            scheduler.step(val_metrics.mutation_macro_f1)
            new_lr = optimizer.param_groups[0]["lr"]
            if new_lr != current_lr:
                logger.info("LR reduced: %.2e -> %.2e", current_lr, new_lr)

            improved = val_metrics.mutation_macro_f1 > best_score + EARLY_STOPPING_MIN_DELTA
            if improved:
                best_score = val_metrics.mutation_macro_f1
                best_val_loss = val_metrics.loss
                epochs_without_improvement = 0
                torch.save(
                    {
                        "epoch": epoch,
                        "model_state_dict": model.state_dict(),
                        "optimizer_state_dict": optimizer.state_dict(),
                        "best_score": best_score,
                        "val_loss": val_metrics.loss,
                        "val_macro_f1": val_metrics.macro_f1,
                        "val_mutation_macro_f1": val_metrics.mutation_macro_f1,
                    },
                    BEST_MODEL_PATH,
                )
                logger.info(
                    "New best model saved (mutation_macro_f1=%.4f) -> %s",
                    best_score, BEST_MODEL_PATH,
                )
            else:
                epochs_without_improvement += 1

            torch.save(
                {
                    "epoch": epoch,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "best_score": best_score,
                },
                LAST_CHECKPOINT_PATH,
            )

            if epochs_without_improvement >= EARLY_STOPPING_PATIENCE:
                logger.info(
                    "Early stopping: no mutation_macro_f1 improvement > %.1e for %d epochs.",
                    EARLY_STOPPING_MIN_DELTA, EARLY_STOPPING_PATIENCE,
                )
                break

        logger.info("Training complete. Best validation mutation macro F1: %.4f", best_score)


if __name__ == "__main__":
    main()
