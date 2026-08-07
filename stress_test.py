"""Ablation harness for the gap-padded generator and auxiliary channels.

Runs the same data and training schedule under four input configurations so
that any metric change is attributable to the channels rather than to
seeding or schedule differences. The point is to measure what quality and
depth each contribute, rather than assume that adding a channel helped.

Run with::

    python -m bimamba_variant_caller.stress_test --epochs 12
"""

from __future__ import annotations

import argparse
from typing import Dict, List, Tuple

import torch

from config import LABEL_NAMES, DataConfig, ExperimentConfig, ModelConfig, TrainingConfig
from dataset import SyntheticVariantDataset
from train import Trainer


def build_experiment(
    use_base_quality: bool,
    use_depth: bool,
    epochs: int,
    seq_len: int,
    mutation_rate: float,
    noise_rate: float,
) -> ExperimentConfig:
    """Assemble one ablation arm.

    Every arm shares seeds, sizes and schedule; only the input channels vary.

    Args:
        use_base_quality: Whether the Phred channel is consumed.
        use_depth: Whether the coverage channel is consumed.
        epochs: Training epochs.
        seq_len: Alignment columns per window.
        mutation_rate: Per-position true-variant rate.
        noise_rate: Per-base sequencing error rate.

    Returns:
        The configured experiment.
    """
    return ExperimentConfig(
        data=DataConfig(
            train_samples=1000,
            val_samples=200,
            seq_len=seq_len,
            mutation_rate=mutation_rate,
            noise_rate=noise_rate,
        ),
        model=ModelConfig(
            d_model=128,
            num_layers=6,
            dropout=0.15,
            use_base_quality=use_base_quality,
            use_depth=use_depth,
        ),
        training=TrainingConfig(num_epochs=epochs, warmup_epochs=3),
    )


def describe_separability(dataset: SyntheticVariantDataset, samples: int = 300) -> Dict[str, float]:
    """Measure how far the data alone determines the answer.

    Reports the ceiling a mismatch-only rule would reach, which is the
    reference point for judging whether the model is learning something or
    merely counting disagreements.

    Args:
        dataset: Dataset to inspect.
        samples: Number of samples to scan.

    Returns:
        Mapping of diagnostic name to value.
    """
    counted = min(samples, len(dataset))
    mismatch_normal = 0
    total_normal = 0
    snp_mismatch = 0
    total_mismatch = 0
    total_snp = 0

    for index in range(counted):
        sample = dataset[index]
        mismatch = sample["input_ids"] != sample["reference_ids"]
        labels = sample["labels"]

        normal = labels == 0
        snp = labels == 1

        total_normal += int(normal.sum())
        mismatch_normal += int((normal & mismatch).sum())
        total_snp += int(snp.sum())
        snp_mismatch += int((snp & mismatch).sum())
        total_mismatch += int(mismatch.sum())

    precision = snp_mismatch / total_mismatch if total_mismatch else 0.0
    recall = snp_mismatch / total_snp if total_snp else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0

    return {
        "normal_mismatch_fraction": mismatch_normal / max(total_normal, 1),
        "mismatch_rule_precision": precision,
        "mismatch_rule_recall": recall,
        "mismatch_rule_f1": f1,
    }


def run_arm(name: str, config: ExperimentConfig) -> Tuple[str, float, Dict[str, float]]:
    """Train and evaluate one ablation arm.

    Args:
        name: Human-readable arm label.
        config: The experiment to run.

    Returns:
        Tuple of ``(name, val_loss, metrics)``.
    """
    print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")
    trainer = Trainer(config, verbose=True)
    trainer.fit()
    val_loss, metrics = trainer.evaluate()
    return name, val_loss, metrics


def main(argv: List[str] | None = None) -> int:
    """Run all ablation arms and print a comparison table.

    Args:
        argv: Argument list, or ``None`` to read from ``sys.argv``.

    Returns:
        Process exit code.
    """
    parser = argparse.ArgumentParser(description="Ablate the auxiliary input channels.")
    parser.add_argument("--epochs", type=int, default=12)
    parser.add_argument("--seq-len", type=int, default=64)
    parser.add_argument("--mutation-rate", type=float, default=0.05)
    parser.add_argument("--noise-rate", type=float, default=0.0075)
    args = parser.parse_args(argv)

    probe = SyntheticVariantDataset(
        300,
        seq_len=args.seq_len,
        mutation_rate=args.mutation_rate,
        seed=123,
        noise_rate=args.noise_rate,
    )
    diagnostics = describe_separability(probe)

    print("=" * 78)
    print("DATA DIAGNOSTICS (gap-padded alignment)")
    print("=" * 78)
    print(
        f"Normal columns mismatching the reference : "
        f"{diagnostics['normal_mismatch_fraction']:.4f}"
    )
    print(
        f"'call every mismatch a SNP' ceiling      : "
        f"P={diagnostics['mismatch_rule_precision']:.3f} "
        f"R={diagnostics['mismatch_rule_recall']:.3f} "
        f"F1={diagnostics['mismatch_rule_f1']:.3f}"
    )

    arms = [
        ("A: bases only (no quality, no depth)", False, False),
        ("B: bases + quality", True, False),
        ("C: bases + depth", False, True),
        ("D: bases + quality + depth", True, True),
    ]

    results: List[Tuple[str, float, Dict[str, float]]] = []
    for name, quality, depth in arms:
        torch.manual_seed(42)
        config = build_experiment(
            quality, depth, args.epochs, args.seq_len, args.mutation_rate, args.noise_rate
        )
        results.append(run_arm(name, config))

    print("\n" + "=" * 78)
    print("ABLATION SUMMARY")
    print("=" * 78)
    header = f"{'arm':<38}{'val_loss':>10}{'macro':>8}{'mut':>8}" + "".join(
        f"{name[:6]:>9}" for name in LABEL_NAMES
    )
    print(header)
    print("-" * len(header))
    for name, val_loss, metrics in results:
        row = (
            f"{name:<38}{val_loss:>10.4f}"
            f"{metrics['macro_f1']:>8.3f}{metrics['f1_mutations_macro']:>8.3f}"
        )
        row += "".join(f"{metrics[f'f1_class_{index}']:>9.3f}" for index in range(len(LABEL_NAMES)))
        print(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
