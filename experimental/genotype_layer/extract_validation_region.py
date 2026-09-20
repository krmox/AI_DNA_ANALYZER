"""Extract per-locus evidence (counts, binomial LLR, PB LLR, cascade calls) for
the TRUE hold-out validation region chr21:32,000,000-44,000,000
(data/giab_hg002_chr21_12Mb), mirroring extract_dev_region.py exactly (same
frozen imports: pileup_counts.PileupCountsProvider, read_level_pileup.ReadLevelPileupProvider,
binomial_baseline.BinomialVariantCaller, cascade.py's frozen constants) but
chunked over the region (12Mb, vs dev's 0.47Mb) the way
experimental/unified_happy/recover_full_positions.py chunked its counts-only
recomputation, to keep memory bounded. Nothing here modifies any frozen file.

Output: experimental/genotype_layer/cache/validation_region.npz (NEW file;
does not touch cache/dev_region.npz).
"""
import sys
import time

import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
sys.path.insert(0, _ROOT)

import numpy as np

from binomial_baseline import BinomialVariantCaller
from cascade import FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD, FROZEN_ROUTER_CUTOFF
from pileup_counts import PileupCountsProvider
from quality_error_model import candidate_alt, extract_quality_evidence, poisson_binomial_llr
from read_level_pileup import ReadLevelPileupProvider

REPO = _ROOT
DATA = f"{REPO}/data/giab_hg002_chr21_12Mb"
FASTA = f"{REPO}/data/reference/chr21_full.fa"
BAM = f"{DATA}/hg002_chr21_32_44M_15x.bam"
VCF = f"{DATA}/hg002_chr21_32_44M.vcf.gz"
BED = f"{DATA}/hg002_chr21_32_44M_highconf.bed"

REGION_START, REGION_END = 32_000_000, 44_000_000
CHUNK_BP = 500_000
PB_BLOCK = 200_000

started = time.perf_counter()

counts_blocks, labels_blocks, pos_blocks, pbllr_blocks = [], [], [], []
binllr_blocks = []

for chunk_start in range(REGION_START, REGION_END, CHUNK_BP):
    chunk_end = min(chunk_start + CHUNK_BP, REGION_END)

    counts_provider = PileupCountsProvider(
        fasta_path=FASTA, bam_path=BAM, vcf_path=VCF, contig="chr21",
        region=(chunk_start, chunk_end), seq_len=64, high_confidence_bed=BED,
    )
    with counts_provider:
        windows = counts_provider._build_windows()
        c_list, l_list, p_list = [], [], []
        for w in windows:
            ref = counts_provider._fetch_reference(w)
            c_list.append(counts_provider.window_counts(w, ref))
            l_list.append(np.asarray(counts_provider._fetch_labels(w), dtype=np.int64))
            p_list.append(np.arange(w.start, w.end, dtype=np.int64))
        if not c_list:
            continue
        chunk_counts = np.concatenate(c_list).astype(np.float64)
        chunk_labels = np.concatenate(l_list)
        chunk_pos = np.concatenate(p_list)

    reads_provider = ReadLevelPileupProvider(
        fasta_path=FASTA, bam_path=BAM, vcf_path=VCF, contig="chr21",
        region=(chunk_start, chunk_end), seq_len=64, high_confidence_bed=BED,
    )
    with reads_provider:
        r_list = [reads_provider[i]["reads"] for i in range(len(reads_provider))]
    chunk_reads = np.concatenate(r_list)
    assert chunk_reads.shape[0] == chunk_counts.shape[0], (chunk_reads.shape, chunk_counts.shape)

    caller = BinomialVariantCaller(error_rate=0.01)
    chunk_binomial_llr = caller.score_counts(chunk_counts[:, 0:4], chunk_counts[:, 9]).llr

    pb_llr_parts = []
    for start in range(0, chunk_counts.shape[0], PB_BLOCK):
        stop = min(start + PB_BLOCK, chunk_counts.shape[0])
        block_counts = chunk_counts[start:stop]
        reference_index = block_counts[:, 9].astype(int)
        evidence = extract_quality_evidence(chunk_reads[start:stop], reference_index)
        alt_index, k, _ = candidate_alt(block_counts[:, 0:4], reference_index)
        result = poisson_binomial_llr(evidence, k.astype(np.int64), alt_index=alt_index)
        pb_llr_parts.append(result["llr"])
    chunk_pb_llr = np.concatenate(pb_llr_parts).astype(np.float64)

    counts_blocks.append(chunk_counts.astype(np.float32))
    labels_blocks.append(chunk_labels)
    pos_blocks.append(chunk_pos)
    binllr_blocks.append(chunk_binomial_llr)
    pbllr_blocks.append(chunk_pb_llr)

    elapsed = (time.perf_counter() - started) / 60
    print(f"{chunk_start}-{chunk_end}: {len(windows)} windows | cumulative {elapsed:.1f} min", flush=True)

counts = np.concatenate(counts_blocks)
labels = np.concatenate(labels_blocks)
positions = np.concatenate(pos_blocks)
binomial_llr = np.concatenate(binllr_blocks)
pb_llr = np.concatenate(pbllr_blocks)

assert labels.size == positions.size == counts.shape[0] == binomial_llr.size == pb_llr.size

binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
pb_calls_raw = pb_llr >= FROZEN_PB_THRESHOLD
routed = np.abs(binomial_llr - FROZEN_BINOMIAL_THRESHOLD) <= FROZEN_ROUTER_CUTOFF
cascade_calls_raw = np.where(routed, pb_calls_raw, binomial_calls)

np.savez_compressed(
    f"{REPO}/experimental/genotype_layer/cache/validation_region.npz",
    counts=counts.astype(np.float32), labels=labels, positions=positions,
    binomial_llr=binomial_llr, pb_llr=pb_llr, routed=routed,
    binomial_calls=binomial_calls, pb_calls=pb_calls_raw, cascade_calls=cascade_calls_raw,
)
print(f"Saved. loci={labels.size} binomial_calls={binomial_calls.sum()} pb_calls={pb_calls_raw.sum()} "
      f"cascade_calls={cascade_calls_raw.sum()} routed={routed.sum()} "
      f"total_time={(time.perf_counter() - started)/60:.1f}min")
