"""Extract and cache the read-level tensor once, so every consumer is identical.

The pileup walk is the dominant cost and the read-level walk is heavier still,
so it is paid once per region. Caching is also a correctness property, not only
a speed one: the primary model, all six controls and the evaluation then read
*byte-identical* input, which is what makes an ablation a measurement of
information rather than of extraction jitter.

Stored as ``uint8`` compact columns (see ``read_level_pileup.READ_COLUMNS``);
a float32 expansion of the training region would be ~7 GB.

Labels are stored alongside for evaluation only and play no part in feature
construction.
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

import numpy as np

from read_level_pileup import MAX_READS, READ_COLUMNS, load_reads

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--region", type=int, nargs=2, required=True)
    parser.add_argument("--max-reads", type=int, default=MAX_READS)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    start = time.perf_counter()
    reads, labels = load_reads(args.fasta, args.bam, args.vcf, args.bed,
                               tuple(args.region), max_reads=args.max_reads)
    elapsed = time.perf_counter() - start

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out, reads=reads, labels=labels,
        region=np.asarray(args.region), columns=np.asarray(READ_COLUMNS),
        max_reads=np.asarray([args.max_reads]),
        bam=np.asarray([args.bam]), vcf=np.asarray([args.vcf]), bed=np.asarray([args.bed]),
    )
    logger.info("Cached %d loci x %d reads (%d SNP) from %s in %.1fs -> %s (%.0f MB raw)",
                labels.size, args.max_reads, int((labels == 1).sum()), args.bam, elapsed,
                args.out, reads.nbytes / 1e6)


if __name__ == "__main__":
    main()
