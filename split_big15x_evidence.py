"""Materialize the train and validation evidence files for the powered benchmark.

``train_pb_residual`` loads whole ``.npz`` files, so the block split has to be
written out rather than applied in memory. Only train and validation are
written: the test blocks are never materialized as a file, which removes the
most obvious way to accidentally train on them.

The masks come from :mod:`genomic_split` at *window* granularity, so each output
file holds a whole number of 64-locus windows and no window mixes two roles.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np

from genomic_split import assign_window_roles, split_summary
from train_raw_pileup import SEQ_LEN

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", default="cache/big15x_evidence.npz")
    parser.add_argument("--positions", default="cache/big15x_positions.npy")
    parser.add_argument("--region-start", type=int, default=32000000)
    parser.add_argument("--out-prefix", default="cache/big15x")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    blob = np.load(args.evidence, allow_pickle=True)
    positions = np.load(args.positions)
    labels = blob["labels"]
    assert positions.size == labels.size, "position/locus count mismatch"

    masks = assign_window_roles(positions, args.region_start, seq_len=SEQ_LEN)
    summary = split_summary(positions, labels, args.region_start, seq_len=SEQ_LEN)
    for role, entry in summary.items():
        logger.info("%-11s %9d loci (%.2f Mb) SNP=%5d blocks=%s",
                    role, entry["loci"], entry["megabases"], entry["snp"], entry["blocks"])

    for role in ("train", "validation"):
        mask = masks[role]
        assert mask.sum() % SEQ_LEN == 0, f"{role} is not window-aligned"
        out = f"{args.out_prefix}_{role}_evidence.npz"
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            out,
            features=blob["features"][mask], pb_llr=blob["pb_llr"][mask],
            labels=labels[mask], depth=blob["depth"][mask],
            feature_names=blob["feature_names"], positions=positions[mask],
            role=np.asarray([role]), source=np.asarray([args.evidence]),
        )
        logger.info("Wrote %s | %d loci | %d SNP", out, int(mask.sum()),
                    int((labels[mask] == 1).sum()))

    logger.info("Test blocks (%d loci) deliberately not written to disk.",
                int(masks["test"].sum()))


if __name__ == "__main__":
    main()
