"""Baselines for the raw-pileup experiment, on the same test region and BED.

Two baselines, answering different questions:

1. **Historical token baseline** ``input_id != reference_id`` -- the shortcut
   that made every previous result uninterpretable. Uses the *old* consensus
   provider, so it is only meaningful as a record of what the previous
   representation gave away for free.
2. **Raw-pileup statistical caller** -- a threshold on non-reference allele
   fraction computed from the new features. This is the bar the model must
   clear: if a one-parameter frequency rule matches the model, the model has
   learned nothing a caller could not do with arithmetic.

The frequency threshold for baseline 2 is selected on the *training* region,
never on test, mirroring the model's threshold discipline.
"""

from __future__ import annotations

import argparse

import torch
from torch.utils.data import DataLoader

from config import LABEL_SNP
from raw_pileup import RawPileupDataset, RawPileupProvider


def load_region(fasta, bam, vcf, bed, region, max_windows=None):
    """Materialize a region's raw-pileup features and labels."""
    provider = RawPileupProvider(
        fasta_path=fasta, bam_path=bam, vcf_path=vcf, contig="chr21",
        region=region, seq_len=64, high_confidence_bed=bed, max_windows=max_windows,
    )
    with provider:
        batches = list(DataLoader(RawPileupDataset(provider), batch_size=16, shuffle=False))
    features = torch.cat([b["pileup_features"].reshape(-1, b["pileup_features"].shape[-1])
                          for b in batches])
    labels = torch.cat([b["labels"].reshape(-1) for b in batches])
    return features, labels


def non_reference_fraction(features: torch.Tensor) -> torch.Tensor:
    """Fraction of the read stack disagreeing with the reference base.

    Derived here, in the baseline, from the same raw channels the model sees.
    It is deliberately NOT a model input -- computing it is exactly the
    comparison the model is supposed to learn.

    Args:
        features: ``[N, FEATURE_DIM]`` raw-pileup features.

    Returns:
        ``[N]`` non-reference allele fraction in [0, 1].
    """
    reference_index = features[:, 9:14].argmax(dim=1)   # 0=N,1=A,2=C,3=G,4=T
    reference_fraction = torch.zeros(features.shape[0])
    for slot in range(1, 5):
        mask = reference_index == slot
        reference_fraction[mask] = features[mask, slot - 1]
    covered = features[:, 0:5].sum(dim=1) > 0
    return torch.where(covered, 1.0 - reference_fraction, torch.zeros_like(reference_fraction))


def score(predicted: torch.Tensor, labels: torch.Tensor, name: str) -> dict:
    """Print and return SNP-vs-rest precision/recall/F1 for a boolean call."""
    truth = labels == LABEL_SNP
    tp = int((predicted & truth).sum())
    fp = int((predicted & ~truth).sum())
    fn = int((~predicted & truth).sum())
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)
    print(f"  {name:<52} P={precision:.4f} R={recall:.4f} F1={f1:.4f} "
          f"TP={tp} FP={fp} FN={fn}")
    return {"precision": precision, "recall": recall, "f1": f1, "tp": tp, "fp": fp, "fn": fn}


def main() -> None:
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
    parser.add_argument("--train-max-windows", type=int, default=6000)
    args = parser.parse_args()

    print("Selecting the frequency threshold on the TRAINING region only...")
    tr_features, tr_labels = load_region(
        args.fasta, args.train_bam, args.train_vcf, args.train_bed,
        tuple(args.train_region), args.train_max_windows)
    tr_fraction = non_reference_fraction(tr_features)

    best_threshold, best_f1 = 0.0, -1.0
    for candidate in [i / 100 for i in range(5, 96)]:
        stats = {"tp": int(((tr_fraction >= candidate) & (tr_labels == LABEL_SNP)).sum()),
                 "fp": int(((tr_fraction >= candidate) & (tr_labels != LABEL_SNP)).sum()),
                 "fn": int(((tr_fraction < candidate) & (tr_labels == LABEL_SNP)).sum())}
        precision = stats["tp"] / max(stats["tp"] + stats["fp"], 1)
        recall = stats["tp"] / max(stats["tp"] + stats["fn"], 1)
        f1 = 2 * precision * recall / max(precision + recall, 1e-12)
        if f1 > best_f1:
            best_threshold, best_f1 = candidate, f1
    print(f"  chosen non-reference-fraction threshold = {best_threshold:.2f} "
          f"(train F1={best_f1:.4f})")

    te_features, te_labels = load_region(
        args.fasta, args.test_bam, args.test_vcf, args.test_bed, tuple(args.test_region))
    te_fraction = non_reference_fraction(te_features)

    print(f"\nTEST region: {te_labels.numel():,} loci, "
          f"{int((te_labels == LABEL_SNP).sum())} SNP labels")
    score(te_fraction >= best_threshold, te_labels,
          f"raw-pileup frequency rule (thr={best_threshold:.2f}, train-derived)")
    for reference_threshold in (0.15, 0.25, 0.50):
        score(te_fraction >= reference_threshold, te_labels,
              f"raw-pileup frequency rule (thr={reference_threshold:.2f}, reference point)")


if __name__ == "__main__":
    main()
