"""Sanity checks for the raw-pileup representation. Read-only, no training.

Run before any training run. Verifies shapes, finiteness, distributions, and
-- most importantly -- that no feature channel is a proxy for "this locus
differs from the reference".
"""

from __future__ import annotations

import argparse

import torch
from torch.utils.data import DataLoader

from config import LABEL_NAMES
from raw_pileup import FEATURE_DIM, FEATURE_NAMES, RawPileupDataset, RawPileupProvider


def main() -> None:
    """Run every check and print PASS/FAIL per item."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--bed", default=None)
    parser.add_argument("--region-start", type=int, required=True)
    parser.add_argument("--region-end", type=int, required=True)
    parser.add_argument("--max-windows", type=int, default=2000)
    args = parser.parse_args()

    provider = RawPileupProvider(
        fasta_path=args.fasta, bam_path=args.bam, vcf_path=args.vcf,
        contig="chr21", region=(args.region_start, args.region_end),
        seq_len=64, high_confidence_bed=args.bed, max_windows=args.max_windows,
    )

    results: list[tuple[str, bool, str]] = []

    with provider:
        dataset = RawPileupDataset(provider)
        batches = list(DataLoader(dataset, batch_size=16, shuffle=False))
        features = torch.cat([b["pileup_features"].reshape(-1, FEATURE_DIM) for b in batches])
        labels = torch.cat([b["labels"].reshape(-1) for b in batches])
        depth = torch.cat([b["depth"].reshape(-1) for b in batches])

        print(f"windows={len(dataset)} tokens={labels.numel():,} feature_dim={FEATURE_DIM}")

        sample = batches[0]["pileup_features"]
        results.append((
            "tensor shape [B, L, FEATURE_DIM]",
            sample.ndim == 3 and sample.shape[1] == 64 and sample.shape[2] == FEATURE_DIM,
            str(tuple(sample.shape)),
        ))
        results.append((
            "all values finite",
            bool(torch.isfinite(features).all()),
            f"min={features.min():.4f} max={features.max():.4f}",
        ))
        results.append((
            "labels align 1:1 with loci",
            labels.numel() == features.shape[0],
            f"{labels.numel()} labels vs {features.shape[0]} loci",
        ))

        counts = torch.bincount(labels, minlength=4)
        present = {LABEL_NAMES[i]: int(counts[i]) for i in range(4)}
        results.append((
            "SNP / INS / DEL and Normal all represented",
            all(counts[i] > 0 for i in range(4)),
            str(present),
        ))

        median_depth = float(depth.median())
        results.append((
            "depth distribution plausible",
            0 < median_depth < 1000,
            f"min={depth.min():.0f} median={median_depth:.0f} max={depth.max():.0f}",
        ))

        base_fractions = features[:, 0:4].sum(dim=1)
        covered = depth > 0
        results.append((
            "base fractions + del fraction sum to ~1 where covered",
            bool((((base_fractions + features[:, 4])[covered] - 1.0).abs() < 1e-3).all()),
            f"mean sum={(base_fractions + features[:, 4])[covered].mean():.4f}",
        ))

        # --- the leak checks ---------------------------------------------
        reference_onehot = features[:, 9:14]
        # Reconstruct "which base does the reference say" and "which base does
        # the evidence most support", then confirm no single channel already
        # encodes their disagreement.
        ref_index = reference_onehot.argmax(dim=1)          # 0=N,1=A,2=C,3=G,4=T
        ref_base_fraction = torch.zeros_like(base_fractions)
        for slot in range(1, 5):                             # A,C,G,T -> cols 0..3
            mask = ref_index == slot
            ref_base_fraction[mask] = features[mask, slot - 1]

        differs = torch.zeros_like(labels, dtype=torch.bool)
        differs[covered] = ref_base_fraction[covered] < 0.85  # loose proxy for "non-reference"

        worst_corr = 0.0
        worst_name = ""
        target = differs.float()
        for index in range(FEATURE_DIM):
            channel = features[:, index]
            if channel.std() < 1e-8:
                continue
            corr = float(
                ((channel - channel.mean()) * (target - target.mean())).mean()
                / (channel.std() * target.std() + 1e-12)
            )
            if abs(corr) > abs(worst_corr):
                worst_corr, worst_name = corr, FEATURE_NAMES[index]
        results.append((
            "no single channel is a 'differs from reference' flag",
            abs(worst_corr) < 0.99,
            f"max |corr| = {abs(worst_corr):.4f} ({worst_name})",
        ))

        exact = [
            FEATURE_NAMES[i] for i in range(FEATURE_DIM)
            if set(torch.unique(features[:, i]).tolist()) <= {0.0, 1.0}
            and FEATURE_NAMES[i] not in
            ("ref_is_N", "ref_is_A", "ref_is_C", "ref_is_G", "ref_is_T")
        ]
        results.append((
            "no binary non-reference indicator channel",
            not exact,
            f"binary channels outside reference one-hot: {exact}",
        ))

        # --- example loci -------------------------------------------------
        print("\nExample loci (reference | raw evidence | label). Truth is printed"
              "\nfor inspection only; it plays no part in building the tensor.")
        for target_label in range(4):
            hits = (labels == target_label).nonzero().flatten()
            if not len(hits):
                continue
            index = int(hits[len(hits) // 2])
            row = features[index]
            ref_name = ["N", "A", "C", "G", "T"][int(row[9:14].argmax())]
            print(
                f"  {LABEL_NAMES[target_label]:<9} ref={ref_name} "
                f"A={row[0]:.2f} C={row[1]:.2f} G={row[2]:.2f} T={row[3]:.2f} "
                f"del={row[4]:.2f} ins={row[5]:.2f} depth={depth[index]:.0f} "
                f"bq={row[7]:.2f} mq={row[8]:.2f}"
            )

    print()
    ok = True
    for name, passed, detail in results:
        print(f"  [{'PASS' if passed else 'FAIL'}] {name:<58} {detail}")
        ok &= passed
    print(f"\n{'ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED'}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
