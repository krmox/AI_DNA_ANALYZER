"""Command-line entry point for training and evaluating the variant caller.

Run from the project root::

    python -m bimamba_variant_caller.main
    python -m bimamba_variant_caller.main --epochs 15 --fusion sum
"""

from __future__ import annotations

import argparse
import random
from typing import List

import torch

from config import LABEL_NAMES, DataConfig, ExperimentConfig, ModelConfig, TrainingConfig
from dataset import DNATokenizer, SyntheticVariantDataset
from model import VariantCaller
from train import Trainer


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list, or ``None`` to read from ``sys.argv``.

    Returns:
        The parsed namespace.
    """
    parser = argparse.ArgumentParser(
        description="Train the reference-guided BiMamba variant caller on synthetic reads.",
    )
    parser.add_argument("--epochs", type=int, default=15, help="Maximum training epochs.")
    parser.add_argument("--batch-size", type=int, default=16, help="Mini-batch size.")
    parser.add_argument("--seq-len", type=int, default=64, help="Sequence length in tokens.")
    parser.add_argument("--train-samples", type=int, default=1000, help="Training sequences per epoch.")
    parser.add_argument("--val-samples", type=int, default=200, help="Validation sequences.")
    parser.add_argument("--mutation-rate", type=float, default=0.05, help="Per-position mutation rate.")
    parser.add_argument(
        "--noise-rate",
        type=float,
        default=0.0075,
        help="Per-base sequencing substitution error rate. Never labelled.",
    )
    parser.add_argument(
        "--max-insertion-len", type=int, default=6, help="Upper bound of Uniform(1, L) insertions."
    )
    parser.add_argument("--d-model", type=int, default=128, help="Hidden width.")
    parser.add_argument("--num-layers", type=int, default=6, help="Number of BiMamba blocks.")
    parser.add_argument("--dropout", type=float, default=0.15, help="Dropout probability.")
    parser.add_argument(
        "--fusion",
        choices=("concat", "sum"),
        default="concat",
        help="How reference and observed embeddings are combined.",
    )
    parser.add_argument("--lr", type=float, default=1e-3, help="Peak learning rate.")
    parser.add_argument("--warmup-epochs", type=int, default=3, help="Linear warmup epochs.")
    parser.add_argument("--patience", type=int, default=5, help="Early-stopping patience.")
    parser.add_argument("--seed", type=int, default=42, help="Global torch seed.")
    parser.add_argument(
        "--device",
        choices=("auto", "cpu", "cuda"),
        default="auto",
        help="Compute device.",
    )
    parser.add_argument(
        "--no-base-quality",
        action="store_true",
        help="Ablation: drop the Phred channel from the input embedding.",
    )
    parser.add_argument(
        "--no-depth",
        action="store_true",
        help="Ablation: drop the coverage channel from the input embedding.",
    )
    parser.add_argument(
        "--no-vaf",
        action="store_true",
        help="Ablation: drop the allele-fraction channel (inert on synthetic data).",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress per-epoch logging.")
    return parser.parse_args(argv)


def build_config(args: argparse.Namespace) -> ExperimentConfig:
    """Assemble an :class:`~.config.ExperimentConfig` from parsed arguments.

    Args:
        args: Parsed command-line arguments.

    Returns:
        The fully populated experiment configuration.
    """
    return ExperimentConfig(
        data=DataConfig(
            train_samples=args.train_samples,
            val_samples=args.val_samples,
            seq_len=args.seq_len,
            mutation_rate=args.mutation_rate,
            noise_rate=args.noise_rate,
            max_insertion_len=args.max_insertion_len,
            batch_size=args.batch_size,
        ),
        model=ModelConfig(
            d_model=args.d_model,
            num_layers=args.num_layers,
            dropout=args.dropout,
            fusion=args.fusion,
            use_base_quality=not args.no_base_quality,
            use_depth=not args.no_depth,
            use_vaf=not args.no_vaf,
        ),
        training=TrainingConfig(
            num_epochs=args.epochs,
            learning_rate=args.lr,
            warmup_epochs=args.warmup_epochs,
            early_stopping_patience=args.patience,
            seed=args.seed,
            device=args.device,
        ),
    )


def show_prediction(
    model: VariantCaller,
    dataset: SyntheticVariantDataset,
    device: torch.device,
    tokenizer: DNATokenizer | None = None,
) -> None:
    """Print a side-by-side prediction for one random validation sample.

    Aggregate F1 hides *which* positions fail; this shows the actual base
    context of every disagreement, which is where insertion/SNP confusions
    become visible.

    Args:
        model: Trained model, switched to eval mode internally.
        dataset: Validation dataset to draw a sample from.
        device: Device the model lives on.
        tokenizer: Tokenizer for decoding, or ``None`` to build one.
    """
    tokenizer = tokenizer or DNATokenizer()
    model.eval()

    sample = dataset[random.randrange(len(dataset))]
    input_ids = sample["input_ids"].unsqueeze(0).to(device)
    reference_ids = sample["reference_ids"].unsqueeze(0).to(device)
    base_quality = sample["base_quality"].unsqueeze(0).to(device)
    depth = sample["depth"].unsqueeze(0).to(device)
    vaf = sample["vaf"].unsqueeze(0).to(device)
    truth = sample["labels"]

    with torch.no_grad():
        prediction = (
            model(input_ids, reference_ids, base_quality, depth, vaf)
            .argmax(dim=-1)
            .squeeze(0)
            .cpu()
        )

    reference = tokenizer.decode(sample["reference_ids"])
    observed = tokenizer.decode(sample["input_ids"])

    print("\n" + "=" * 78)
    print("SAMPLE PREDICTION")
    print("=" * 78)
    print(f"Reference: {reference}")
    print(f"Observed:  {observed}")
    print(f"Truth:     {''.join(str(label.item()) for label in truth)}")
    print(f"Predicted: {''.join(str(label.item()) for label in prediction)}")

    print("\nNon-trivial positions (truth or prediction is a variant):")
    flagged = False
    for position in range(len(truth)):
        true_label = int(truth[position].item())
        predicted_label = int(prediction[position].item())
        if true_label == 0 and predicted_label == 0:
            continue
        flagged = True
        status = "ok" if true_label == predicted_label else "MISMATCH"
        print(
            f"  pos={position:3d} | ref={reference[position]} obs={observed[position]} | "
            f"Q={float(sample['base_quality'][position]):>5.1f} | "
            f"truth={LABEL_NAMES[true_label]:<10} pred={LABEL_NAMES[predicted_label]:<10} | {status}"
        )
    if not flagged:
        print("  (no variants and no false positives in this sample)")


def main(argv: List[str] | None = None) -> int:
    """Train the model, evaluate it and print a qualitative sample.

    Args:
        argv: Argument list, or ``None`` to read from ``sys.argv``.

    Returns:
        Process exit code; ``0`` on success.
    """
    args = parse_args(argv)
    config = build_config(args)

    random.seed(config.training.seed)

    trainer = Trainer(config, verbose=not args.quiet)
    trainer.fit()

    val_loss, metrics = trainer.evaluate()

    print("\n" + "=" * 78)
    print("FINAL VALIDATION (best checkpoint)")
    print("=" * 78)
    print(f"val_loss           : {val_loss:.6f}")
    print(f"macro_f1           : {metrics['macro_f1']:.4f}")
    print(f"f1_mutations_macro : {metrics['f1_mutations_macro']:.4f}")
    for class_id, name in enumerate(LABEL_NAMES):
        print(
            f"  {name:<10} precision={metrics[f'precision_class_{class_id}']:.4f} "
            f"recall={metrics[f'recall_class_{class_id}']:.4f} "
            f"f1={metrics[f'f1_class_{class_id}']:.4f}"
        )

    show_prediction(trainer.model, trainer.val_dataset, trainer.device)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
