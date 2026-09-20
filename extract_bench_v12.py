"""Extract per-locus calling evidence for the bench_v12 benchmark region.

Why this exists (and why it is not ``extract_region_pb.py``)
------------------------------------------------------------
``extract_region_pb.py`` (devlog 14) retains the 44-channel locus features only
for windows the *frozen Stage-1 router* touches. That is exactly the constraint
that made devlog 16's gen14 Stage-2 pool "feature-constrained": a wider
PB-uncertain candidate pool could not be scored by the neural arm because its
features were never written to disk.

This script keeps everything ``extract_region_pb.py`` keeps and additionally
retains the 44-channel features for **every 64-locus window containing at least
one locus inside a pre-registered PB-uncertain candidate pool**
(``|pb_llr - FROZEN_PB_THRESHOLD| <= FEATURE_POOL_MARGIN``), so a genuinely
wide Stage-2 pool can be evaluated without re-touching the BAM.

Both retention rules depend only on caller scores, never on truth labels, so
nothing about the retained set leaks the answer.

Every statistical quantity delegates to the existing implementations
(``load_counts``, ``load_reads``, ``poisson_binomial_llr``,
``BinomialVariantCaller``, ``build_locus_features``) -- no formula is
duplicated here.
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

import numpy as np

from binomial_baseline import BinomialVariantCaller
from locus_evidence import build_locus_features
from pileup_counts import load_counts
from quality_error_model import candidate_alt, extract_quality_evidence, poisson_binomial_llr
from read_level_pileup import load_reads
from train_raw_pileup import SEQ_LEN

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)

CHUNK_BP = 500_000
PB_BLOCK = 200_000

#: Pre-registered PB-margin used *only* to decide which windows' features are
#: kept on disk. Chosen as devlog 16's ``CANDIDATE_MARGIN`` (13.5) so the
#: Stage-2 pool of this benchmark is a superset of that study's.
FEATURE_POOL_MARGIN = 13.5

FROZEN_BINOMIAL_THRESHOLD = 7.0
FROZEN_PB_THRESHOLD = 10.5
FROZEN_ROUTER_CUTOFF = 5.411872376933351


def poisson_binomial_for_chunk(reads: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """Poisson-binomial LLR for one chunk, blocked to bound peak memory."""
    return poisson_binomial_diagnostics(reads, counts)["pb_llr"]


def poisson_binomial_diagnostics(reads: np.ndarray, counts: np.ndarray,
                                 legacy: bool = False) -> dict[str, np.ndarray]:
    """PB LLR plus the per-locus quantities devlog 13's failure analysis needs.

    ``k_full`` is the ALT count over *all* reads (what ``candidate_alt`` returns
    from the uncapped count matrix); ``k_retained`` is the ALT count over the
    <=48 reads the read tensor kept, i.e. the count the Poisson-binomial pmf is
    actually a distribution over. Their difference is the depth-clipping
    condition, so both are recorded rather than inferred.

    Args:
        reads: ``[N, R, READ_DIM]`` read tensor for the chunk.
        counts: ``[N, 10]`` count matrix for the same loci.
        legacy: Additionally compute the pre-fix (defective) LLR, so the fix can
            be A/B-compared on identical reads in a single pass.

    Returns:
        Dict of ``[N]`` arrays: ``pb_llr``, ``k_full``, ``k_retained``,
        ``n_counted``, and (when ``legacy``) ``pb_llr_legacy``.
    """
    keys = ("pb_llr", "k_full", "k_retained", "n_counted") + (("pb_llr_legacy",) if legacy else ())
    blocks: dict[str, list[np.ndarray]] = {name: [] for name in keys}
    # The legacy arm is a benchmark artefact, so its cost is timed out of the
    # production path rather than folded into it.
    seconds = {"fixed": 0.0, "legacy": 0.0}
    for start in range(0, counts.shape[0], PB_BLOCK):
        stop = min(start + PB_BLOCK, counts.shape[0])
        block_counts = counts[start:stop]
        reference_index = block_counts[:, 9].astype(int)
        clock = time.perf_counter()
        evidence = extract_quality_evidence(reads[start:stop], reference_index)
        alt_index, k, _ = candidate_alt(block_counts[:, 0:4], reference_index)
        result = poisson_binomial_llr(evidence, k.astype(np.int64), alt_index=alt_index)
        seconds["fixed"] += time.perf_counter() - clock
        blocks["pb_llr"].append(result["llr"])
        blocks["k_full"].append(k.astype(np.int32))
        blocks["k_retained"].append(result["k_effective"].astype(np.int32))
        blocks["n_counted"].append(evidence.n_counted.astype(np.int32))
        if legacy:
            clock = time.perf_counter()
            blocks["pb_llr_legacy"].append(poisson_binomial_llr(
                evidence, k.astype(np.int64), legacy_truncation=True)["llr"])
            seconds["legacy"] += time.perf_counter() - clock
    out = {name: np.concatenate(values) for name, values in blocks.items()}
    out["pb_llr"] = out["pb_llr"].astype(np.float64)
    if legacy:
        out["pb_llr_legacy"] = out["pb_llr_legacy"].astype(np.float64)
    out["_seconds_fixed"] = seconds["fixed"]
    out["_seconds_legacy"] = seconds["legacy"]
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--region", type=int, nargs=2, required=True)
    parser.add_argument("--contig", default="chr21",
                        help=("Contig to extract. Defaults to chr21, the only contig used "
                              "before devlog 14, so every previous invocation of this "
                              "script is reproduced byte-for-byte by omitting the flag."))
    parser.add_argument("--chunk-bp", type=int, default=CHUNK_BP)
    parser.add_argument("--features", choices=("full", "none", "scores"), default="full",
                        help=("full: keep the 44-channel features for every pool window "
                              "(~150 MB/Mb of genome). none: Stage-1 evidence only. "
                              "scores: score the pool windows with the benchmark's "
                              "checkpoints in-process and keep only the per-locus "
                              "decision scores -- same numbers, ~30x less disk."))
    parser.add_argument("--legacy-pb", action="store_true",
                        help=("also store the pre-devlog-13 (defective) PB LLR under "
                              "``pb_llr_legacy``, computed on the same reads in the same "
                              "pass, so the depth-clipping fix can be A/B-compared exactly."))
    parser.add_argument("--workers", type=int, default=1,
                        help=("process chunks in this many worker processes (chunks are independent "
                              "and concatenated in order, so output is byte-identical to --workers 1)."))
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def score_windows(features: np.ndarray, llr: np.ndarray) -> dict:
    """Mamba decision scores for one chunk's retained windows.

    Kept identical to what ``bench_v12_stage23`` would compute from stored
    features: same ``model_scores`` call, same checkpoints, same window
    reshaping. Chunk boundaries are whole windows, so per-chunk scoring is
    exact rather than an approximation.
    """
    import torch

    from bench_v12_stage23 import CHECKPOINTS, LEAKED_CHECKPOINTS
    from evaluate_pb_residual import model_scores

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    out = {}
    for name, path in {**CHECKPOINTS, **LEAKED_CHECKPOINTS}.items():
        score, diagnostics, _ = model_scores(Path(path), features, llr, device)
        out[f"score__{name}"] = score.astype(np.float32)
        out[f"own_threshold__{name}"] = np.asarray(
            [diagnostics["validation_threshold"]], dtype=np.float64)
    return out


def process_chunk(args: argparse.Namespace, chunk_start: int, chunk_end: int) -> dict | None:
    """One chunk's evidence. A pure function of its arguments (no shared state), so chunks can run
    in any order or in separate processes and be concatenated in chunk order with a result
    identical to the serial loop."""
    caller = BinomialVariantCaller(error_rate=0.01)
    blocks: dict[str, np.ndarray] = {}
    timing = {"counts_seconds": 0.0, "reads_seconds": 0.0, "pb_seconds": 0.0,
              "feature_seconds": 0.0}

    start = time.perf_counter()
    counts, labels, positions = load_counts(args.fasta, args.bam, args.vcf, args.bed,
                                            (chunk_start, chunk_end), contig=args.contig,
                                            return_positions=True)
    timing["counts_seconds"] += time.perf_counter() - start
    if labels.size == 0:
        return None  # no 64 bp window of this chunk lies in the high-confidence BED

    start = time.perf_counter()
    reads, read_labels = load_reads(args.fasta, args.bam, args.vcf, args.bed,
                                    (chunk_start, chunk_end), contig=args.contig)
    timing["reads_seconds"] += time.perf_counter() - start
    assert np.array_equal(read_labels, labels), "count/read label mismatch"

    pb_diagnostics = poisson_binomial_diagnostics(reads, counts, legacy=args.legacy_pb)
    pb = pb_diagnostics["pb_llr"]
    timing["pb_seconds"] += pb_diagnostics["_seconds_fixed"]
    for name in ("k_full", "k_retained", "n_counted", "pb_llr_legacy"):
        if name in pb_diagnostics:
            blocks[name] = pb_diagnostics[name]
    binomial = caller.score_counts(counts[:, 0:4], counts[:, 9].astype(int)).llr

    pool = np.abs(pb - FROZEN_PB_THRESHOLD) <= FEATURE_POOL_MARGIN
    window_touched = pool.reshape(-1, SEQ_LEN).any(axis=1)
    if window_touched.any() and args.features != "none":
        selected = np.repeat(window_touched, SEQ_LEN)
        start = time.perf_counter()
        chunk_features = build_locus_features(reads[selected], counts[selected])
        timing["feature_seconds"] += time.perf_counter() - start
        blocks["pool_window_positions"] = positions[np.nonzero(selected)[0]].astype(np.int64)
        if args.features == "full":
            blocks["pool_window_features"] = chunk_features
        else:
            for key, value in score_windows(chunk_features, pb[selected]).items():
                blocks[key] = value
        del chunk_features

    router_touched = (np.abs(binomial - FROZEN_BINOMIAL_THRESHOLD)
                      <= FROZEN_ROUTER_CUTOFF).reshape(-1, SEQ_LEN).any(axis=1)
    if router_touched.any():
        blocks["router_window_positions"] = positions[
            np.nonzero(np.repeat(router_touched, SEQ_LEN))[0]].astype(np.int64)

    del reads

    if args.features != "scores":
        blocks["counts"] = counts.astype(np.float32)
    blocks["pb_llr"] = pb
    blocks["binomial_llr"] = binomial
    blocks["depth_"] = counts[:, 6].astype(np.float32)
    blocks["labels"] = labels
    blocks["positions"] = positions
    return {"blocks": blocks, "timing": timing, "n_loci": int(labels.size),
            "n_snp": int((labels == 1).sum()), "n_pool_windows": int(window_touched.sum())}


def main() -> None:
    args = parse_args()
    region_start, region_end = args.region

    blocks: dict[str, list[np.ndarray]] = {
        "counts": [], "pb_llr": [], "binomial_llr": [], "labels": [], "positions": [],
        "pool_window_features": [], "pool_window_positions": [],
        "router_window_positions": [],
    }
    depth_blocks: list[np.ndarray] = []
    timing = {"counts_seconds": 0.0, "reads_seconds": 0.0, "pb_seconds": 0.0,
              "feature_seconds": 0.0}
    started = time.perf_counter()

    chunks = [(c, min(c + args.chunk_bp, region_end))
              for c in range(region_start, region_end, args.chunk_bp)]
    if args.workers > 1:
        assert args.features == "none", "--workers > 1 is only supported with --features none"
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        executor = ProcessPoolExecutor(args.workers, mp_context=mp.get_context("fork"))
        results = executor.map(process_chunk, [args] * len(chunks),
                               [c[0] for c in chunks], [c[1] for c in chunks])
    else:
        executor = None
        results = (process_chunk(args, c0, c1) for c0, c1 in chunks)

    for (chunk_start, chunk_end), result in zip(chunks, results):
        if result is None:
            logger.info("%d-%d: no confident windows, skipped", chunk_start, chunk_end)
            continue
        for name, value in result["blocks"].items():
            if name == "depth_":
                depth_blocks.append(value)
            else:
                blocks.setdefault(name, []).append(value)
        for name, value in result["timing"].items():
            timing[name] += value
        logger.info("%d-%d: %d loci (%d SNP) | pool windows %d | cumulative %.1f min",
                    chunk_start, chunk_end, result["n_loci"], result["n_snp"],
                    result["n_pool_windows"], (time.perf_counter() - started) / 60)
    if executor is not None:
        executor.shutdown()
    wall_seconds = time.perf_counter() - started

    payload = {name: np.concatenate(values) for name, values in blocks.items() if values}
    if "counts" in payload:
        payload["depth"] = payload["counts"][:, 6].astype(np.float32)
    else:
        payload["depth"] = np.concatenate(depth_blocks).astype(np.float32)
    assert payload["labels"].size % SEQ_LEN == 0, "region is not window aligned"

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out, **payload,
        region=np.asarray(args.region), contig=np.asarray([args.contig]),
        bam=np.asarray([args.bam]),
        vcf=np.asarray([args.vcf]), bed=np.asarray([args.bed]),
        fasta=np.asarray([args.fasta]),
        feature_pool_margin=np.asarray([FEATURE_POOL_MARGIN]),
        timing=np.asarray([timing["counts_seconds"], timing["reads_seconds"],
                           timing["pb_seconds"], timing["feature_seconds"]]),
    )
    logger.info("Wrote %s | %d loci | %d SNP | counts %.0fs reads %.0fs pb %.0fs feat %.0fs "
                "(stage seconds summed over workers; wall %.0fs)",
                args.out, payload["labels"].size, int((payload["labels"] == 1).sum()),
                timing["counts_seconds"], timing["reads_seconds"], timing["pb_seconds"],
                timing["feature_seconds"], wall_seconds)


if __name__ == "__main__":
    main()
