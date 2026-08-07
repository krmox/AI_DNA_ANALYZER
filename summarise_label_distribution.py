"""Pre-training diagnostic for real alignment data.

Run this **before** any training on a new BAM/VCF pair. It answers three
questions that otherwise only surface after a wasted training run:

1. *Are the coordinates right?* A near-total absence of one class almost
   always means a contig-naming or off-by-one bug rather than a quiet
   region. Real chr21 carries variants roughly every 1000 bases; a scan
   reporting zero SNPs over a megabase is a bug, not biology.
2. *How severe is the imbalance?* Synthetic density was 1:20. Real density
   is nearer 1:1000, and loss weights calibrated on the former do not
   transfer.
3. *Does the pileup carry signal?* Depth and VAF distributions split by
   class show directly whether the consensus is working: variant loci should
   sit near VAF 0.5 (heterozygous) or 1.0 (homozygous), background near 0.

Usage::

    python -m bimamba_variant_caller.summarise_label_distribution \\
        --fasta GRCh38.fa --bam HG002.chr21.bam --vcf HG002_benchmark.vcf.gz \\
        --bed HG002_highconf.bed --contig chr21 --start 5000000 --end 5200000
"""

from __future__ import annotations

import argparse
from collections import Counter
from typing import Dict, List, Tuple

import torch

from config import LABEL_NAMES
from loss import weights_from_counts
from providers import GiabAlignmentProvider


def scan_provider(provider: GiabAlignmentProvider) -> Tuple[Counter, Dict[str, List[float]], Dict[str, List[float]]]:
    """Walk every window and accumulate per-class statistics.

    Args:
        provider: An opened provider to scan.

    Returns:
        Tuple of ``(label_counts, depth_by_class, vaf_by_class)``.
    """
    label_counts: Counter[str] = Counter()
    depth_by_class: Dict[str, List[float]] = {name: [] for name in LABEL_NAMES}
    vaf_by_class: Dict[str, List[float]] = {name: [] for name in LABEL_NAMES}

    for window in provider:
        labels = window.labels.tolist()
        depths = window.depth.tolist()
        vafs = window.vaf.tolist() if window.vaf is not None else [0.0] * len(labels)

        for label, depth, vaf in zip(labels, depths, vafs):
            name = LABEL_NAMES[int(label)]
            label_counts[name] += 1
            depth_by_class[name].append(depth)
            vaf_by_class[name].append(vaf)

    return label_counts, depth_by_class, vaf_by_class


def _summarise(values: List[float]) -> str:
    """Render mean and median of a list, or a placeholder when empty.

    Args:
        values: Sample values.

    Returns:
        A short formatted summary.
    """
    if not values:
        return "     n/a       n/a"
    ordered = sorted(values)
    mean = sum(values) / len(values)
    median = ordered[len(ordered) // 2]
    return f"{mean:>9.2f} {median:>9.2f}"


def report(
    label_counts: Counter,
    depth_by_class: Dict[str, List[float]],
    vaf_by_class: Dict[str, List[float]],
    alpha: float,
    max_weight_cap: float,
) -> None:
    """Print the diagnostic tables.

    Args:
        label_counts: Per-class token counts.
        depth_by_class: Depth samples per class.
        vaf_by_class: VAF samples per class.
        alpha: Exponent for the recommended weights.
        max_weight_cap: Cap for the recommended weights.
    """
    total = sum(label_counts.values())

    print("=" * 74)
    print("CLASS DISTRIBUTION")
    print("=" * 74)
    print(f"{'class':<12}{'count':>12}{'fraction':>12}{'1 in':>12}")
    print("-" * 74)
    for name in LABEL_NAMES:
        count = label_counts[name]
        fraction = count / total if total else 0.0
        one_in = f"{1 / fraction:,.0f}" if fraction > 0 else "never"
        print(f"{name:<12}{count:>12,}{fraction:>12.6f}{one_in:>12}")
    print(f"{'TOTAL':<12}{total:>12,}")

    variants = total - label_counts["Normal"]
    if total:
        density = variants / total
        print(
            f"\nVariant density: {density:.6f}"
            + (f"  (1 in {1 / density:,.0f} columns)" if density > 0 else "")
        )

    print("\n" + "=" * 74)
    print("PILEUP SIGNAL BY CLASS  (depth and VAF should separate the classes)")
    print("=" * 74)
    print(f"{'class':<12}{'depth mean':>10}{'median':>10}{'VAF mean':>10}{'median':>10}")
    print("-" * 74)
    for name in LABEL_NAMES:
        print(f"{name:<12}{_summarise(depth_by_class[name])}{_summarise(vaf_by_class[name])}")

    counts = torch.tensor([float(label_counts[name]) for name in LABEL_NAMES])
    print("\n" + "=" * 74)
    print("RECOMMENDED LOSS WEIGHTS")
    print("=" * 74)
    normalised = weights_from_counts(counts, alpha=alpha, max_weight_cap=max_weight_cap)
    raw = weights_from_counts(counts, alpha=alpha, max_weight_cap=max_weight_cap, normalise=False)
    print(f"alpha={alpha}, cap={max_weight_cap}")
    for index, name in enumerate(LABEL_NAMES):
        print(f"  {name:<12} normalised={float(normalised[index]):>8.4f}  raw={float(raw[index]):>8.4f}")
    spread = float(normalised.max() / normalised.min()) if float(normalised.min()) > 0 else 0.0
    print(f"  weight ratio (max/min): {spread:.1f}")

    print("\n" + "=" * 74)
    print("SANITY CHECKS")
    print("=" * 74)
    problems: List[str] = []
    if label_counts["SNP"] == 0:
        problems.append("No SNPs found — suspect contig naming or VCF coordinates.")
    if variants == 0:
        problems.append("No variants at all — the VCF query is almost certainly wrong.")
    if total and label_counts["Normal"] / total > 0.9999:
        problems.append("Variant density below 1 in 10,000 — lower than real chr21.")
    mean_variant_depth = depth_by_class["SNP"]
    if mean_variant_depth and sum(mean_variant_depth) / len(mean_variant_depth) < 5:
        problems.append("SNP loci average under 5x depth — pileup may be misaligned.")

    if problems:
        for problem in problems:
            print(f"  [!] {problem}")
    else:
        print("  All checks passed.")


def main(argv: List[str] | None = None) -> int:
    """Scan a region and print the diagnostic report.

    Args:
        argv: Argument list, or ``None`` to read from ``sys.argv``.

    Returns:
        Process exit code.
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--fasta", required=True, help="Indexed reference FASTA.")
    parser.add_argument("--bam", required=True, help="Coordinate-sorted indexed BAM.")
    parser.add_argument("--vcf", required=True, help="Benchmark VCF.")
    parser.add_argument("--bed", default=None, help="Optional high-confidence BED.")
    parser.add_argument("--contig", default="chr21", help="Contig to scan.")
    parser.add_argument("--start", type=int, default=None, help="0-based region start.")
    parser.add_argument("--end", type=int, default=None, help="Exclusive region end.")
    parser.add_argument("--seq-len", type=int, default=64, help="Window width.")
    parser.add_argument("--stride", type=int, default=None, help="Window stride.")
    parser.add_argument("--max-windows", type=int, default=2000, help="Cap on windows scanned.")
    parser.add_argument("--min-mapq", type=int, default=20, help="MAPQ floor.")
    parser.add_argument("--min-depth", type=int, default=5, help="Depth floor for a call.")
    parser.add_argument("--vaf-threshold", type=float, default=0.15, help="Consensus VAF floor.")
    parser.add_argument("--alpha", type=float, default=0.5, help="Weight exponent.")
    parser.add_argument("--max-weight-cap", type=float, default=10.0, help="Weight cap.")
    args = parser.parse_args(argv)

    region = (args.start, args.end) if args.start is not None and args.end is not None else None

    provider = GiabAlignmentProvider(
        fasta_path=args.fasta,
        bam_path=args.bam,
        vcf_path=args.vcf,
        seq_len=args.seq_len,
        contig=args.contig,
        stride=args.stride,
        min_mapping_quality=args.min_mapq,
        min_depth=args.min_depth,
        vaf_threshold=args.vaf_threshold,
        high_confidence_bed=args.bed,
        max_windows=args.max_windows,
        region=region,
    )

    with provider:
        window_count = len(provider)
        print(f"Scanning {window_count} windows of {args.seq_len} columns on {args.contig}...\n")
        if window_count == 0:
            print("[!] No windows produced. Check the contig name, region bounds and BED.")
            return 1
        label_counts, depth_by_class, vaf_by_class = scan_provider(provider)

    report(label_counts, depth_by_class, vaf_by_class, args.alpha, args.max_weight_cap)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
