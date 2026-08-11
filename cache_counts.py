"""Materialize raw pileup counts once and cache them to .npz.

The pileup walk costs ~8 minutes for the 2 Mb training region and dominates
every experiment. Caching it means the real-LLR run, the shuffled-LLR control
and the evaluation all consume *byte-identical* inputs, which is a
reproducibility property as much as a speed one.

Only counts and labels are stored. Labels are kept alongside for evaluation
but play no part in feature construction downstream.
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

import numpy as np

from pileup_counts import COUNT_COLUMNS, load_counts

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--region", type=int, nargs=2, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    start = time.perf_counter()
    counts, labels = load_counts(args.fasta, args.bam, args.vcf, args.bed, tuple(args.region))
    elapsed = time.perf_counter() - start

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out, counts=counts, labels=labels,
        region=np.asarray(args.region), columns=np.asarray(COUNT_COLUMNS),
        bam=np.asarray([args.bam]), vcf=np.asarray([args.vcf]), bed=np.asarray([args.bed]),
    )
    logger.info("Cached %d loci (%d SNP) from %s in %.1fs -> %s",
                labels.size, int((labels == 1).sum()), args.bam, elapsed, args.out)


if __name__ == "__main__":
    main()
