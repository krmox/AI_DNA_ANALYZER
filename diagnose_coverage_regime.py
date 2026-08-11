"""Is a coverage regime trivially separable by an allele-frequency threshold?

Run this BEFORE training. It answers the only question that decides whether a
training run is worth doing: on the held-out test region, do the
non-reference-fraction distributions of Normal and SNP loci overlap?

If they do not overlap, a one-parameter rule is a perfect variant caller and
no model result from that regime can be informative -- the same trap the 74x
experiment fell into, where max(Normal)=0.2500 sat below min(SNP)=0.3077 with
a clean 0.0577 gap.

Read-only: loads no checkpoint, writes nothing, trains nothing.
"""

from __future__ import annotations

import argparse

import torch

from config import LABEL_NAMES, LABEL_SNP
from baseline_raw_pileup import load_region, non_reference_fraction


def quantiles(values: torch.Tensor, points=(0.5, 0.9, 0.99, 0.999, 1.0)) -> str:
    """Format a few order statistics of a 1-D tensor."""
    if values.numel() == 0:
        return "(empty)"
    ordered = values.sort().values
    out = []
    for point in points:
        index = min(int(point * (ordered.numel() - 1)), ordered.numel() - 1)
        out.append(f"p{point*100:g}={ordered[index]:.4f}")
    return " ".join(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--region", type=int, nargs=2, required=True)
    parser.add_argument("--label", default="regime")
    parser.add_argument("--min-depth", type=int, default=5,
                        help="Provider's min_depth, for the uncallable-fraction report.")
    args = parser.parse_args()

    features, labels = load_region(
        args.fasta, args.bam, args.vcf, args.bed, tuple(args.region))
    fraction = non_reference_fraction(features)
    depth = torch.expm1(features[:, 6] * 5.0)

    print(f"\n================ {args.label} ================")
    print(f"loci={labels.numel():,} "
          f"labels={{{', '.join(f'{LABEL_NAMES[i]}:{int((labels==i).sum())}' for i in range(4))}}}")

    print("\n-- coverage --")
    print(f"  min={depth.min():.0f} median={depth.median():.0f} mean={depth.mean():.2f} "
          f"max={depth.max():.0f}")
    print(f"  {quantiles(depth, (0.01, 0.1, 0.5, 0.9, 0.99))}")
    below = int((depth < args.min_depth).sum())
    print(f"  loci below min_depth={args.min_depth}: {below} "
          f"({100*below/depth.numel():.4f}%)")
    hist = torch.histc(depth.clamp(max=150), bins=15, min=0, max=150)
    for i, count in enumerate(hist.tolist()):
        lo, hi = i * 10, (i + 1) * 10
        print(f"    depth {lo:>3}-{hi:<3}: {int(count):>7} {'#' * int(60 * count / hist.max())}")

    normal = fraction[labels == 0]
    snp = fraction[labels == LABEL_SNP]

    print("\n-- non-reference fraction distributions --")
    print(f"  Normal (n={normal.numel():,}): mean={normal.mean():.5f} {quantiles(normal)}")
    print(f"  SNP    (n={snp.numel():,}): mean={snp.mean():.5f} "
          f"min={snp.min():.4f} {quantiles(snp, (0.001, 0.01, 0.1, 0.5))}")

    max_normal = float(normal.max())
    min_snp = float(snp.min())
    gap = min_snp - max_normal
    overlap_normal = int((normal >= min_snp).sum())
    overlap_snp = int((snp <= max_normal).sum())
    print("\n-- overlap --")
    print(f"  max(Normal) = {max_normal:.4f}")
    print(f"  min(SNP)    = {min_snp:.4f}")
    print(f"  gap         = {gap:+.4f}  ({'SEPARABLE' if gap > 0 else 'OVERLAPPING'})")
    print(f"  Normal loci at or above min(SNP): {overlap_normal}")
    print(f"  SNP loci at or below max(Normal): {overlap_snp}")

    print("\n-- best achievable one-parameter rule (SNP vs Normal only, ORACLE) --")
    keep = (labels == 0) | (labels == LABEL_SNP)
    x, y = fraction[keep], labels[keep] == LABEL_SNP
    best = (0.0, -1.0, 0, 0, 0)
    for candidate in [i / 200 for i in range(1, 200)]:
        predicted = x >= candidate
        tp = int((predicted & y).sum())
        fp = int((predicted & ~y).sum())
        fn = int((~predicted & y).sum())
        precision = tp / max(tp + fp, 1)
        recall = tp / max(tp + fn, 1)
        f1 = 2 * precision * recall / max(precision + recall, 1e-12)
        if f1 > best[1]:
            best = (candidate, f1, tp, fp, fn)
    threshold, f1, tp, fp, fn = best
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    print(f"  thr={threshold:.3f} -> P={precision:.4f} R={recall:.4f} F1={f1:.4f} "
          f"TP={tp} FP={fp} FN={fn}")
    print(f"\n  VERDICT: {'STILL TRIVIAL - do not train' if f1 > 0.995 else 'non-trivial regime - training is informative'}")


if __name__ == "__main__":
    main()
