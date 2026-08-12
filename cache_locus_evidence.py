"""Build and cache the 44-channel locus evidence plus the Poisson-binomial prior.

Both are deterministic functions of the already-cached read tensor and count
matrix, so caching is a reproducibility property as much as a speed one: every
variant, control and evaluation then consumes byte-identical input and an
architecture comparison measures architecture rather than extraction jitter.

The PB LLR is the expensive part (an exact 48-step convolution per locus,
~4 minutes for the 1.88 M-locus training region), and it is needed identically
by every arm, so it is computed once here.

Labels are copied through for evaluation only and play no part in feature or
prior construction.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import time
from pathlib import Path

import numpy as np

from evaluate_quality_error import load_region
from locus_evidence import FEATURE_NAMES, build_locus_features
from quality_error_model import (
    candidate_alt,
    extract_quality_evidence,
    poisson_binomial_llr,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


#: Loci per chunk. The per-read intermediates are ``[chunk, 48]`` float64, so
#: the whole 1.88 M-locus training region at once needs tens of GB. Every
#: quantity here is a function of one locus alone -- no window, no neighbour --
#: so chunking is *exact*, not an approximation; ``test_locus_evidence`` asserts
#: chunk-invariance.
CHUNK = 200_000


def build(counts_npz: str, reads_npz: str, chunk: int = CHUNK) -> dict:
    """Materialize features, PB prior and labels for one region.

    Args:
        counts_npz: Cached count matrix.
        reads_npz: Cached read-level tensor.
        chunk: Loci processed per block.

    Returns:
        Dict of arrays ready for ``np.savez_compressed``.
    """
    counts, reads, labels = load_region(counts_npz, reads_npz)
    n_loci = labels.size

    feature_blocks, llr_blocks = [], []
    feature_seconds = prior_seconds = 0.0
    for start_index in range(0, n_loci, chunk):
        stop = min(start_index + chunk, n_loci)
        block_counts = counts[start_index:stop]
        block_reads = reads[start_index:stop]
        reference_index = block_counts[:, 9].astype(int)

        start = time.perf_counter()
        feature_blocks.append(build_locus_features(block_reads, block_counts))
        feature_seconds += time.perf_counter() - start

        start = time.perf_counter()
        evidence = extract_quality_evidence(block_reads, reference_index)
        _, k, _ = candidate_alt(block_counts[:, 0:4], reference_index)
        llr_blocks.append(poisson_binomial_llr(evidence, k.astype(np.int64))["llr"])
        prior_seconds += time.perf_counter() - start
        logger.info("  %d/%d loci", stop, n_loci)

    features = np.concatenate(feature_blocks)
    pb_llr = np.concatenate(llr_blocks)

    logger.info("features %.1fs (%d loci x %d ch), PB prior %.1fs",
                feature_seconds, features.shape[0], features.shape[1], prior_seconds)
    return {
        "features": features,
        "pb_llr": pb_llr.astype(np.float64),
        "labels": labels,
        "depth": counts[:, 6].astype(np.float32),
        "feature_names": np.asarray(FEATURE_NAMES),
        "counts_npz": np.asarray([counts_npz]),
        "reads_npz": np.asarray([reads_npz]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--reads", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    payload = build(args.counts, args.reads)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out, **payload)

    digest = hashlib.sha256(np.ascontiguousarray(payload["features"]).tobytes()).hexdigest()
    logger.info("Wrote %s | %d loci | feature sha256 %s", args.out,
                payload["labels"].size, digest[:16])


if __name__ == "__main__":
    main()
