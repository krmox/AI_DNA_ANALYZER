"""Extract per-locus evidence (counts, binomial LLR, PB LLR, cascade calls) for
the dev region chr21:30,000,000-30,469,999 (data/giab_hg002_real_30M).

Mirrors experimental/unified_happy/recover_full_positions.py's approach
(import pileup_counts.PileupCountsProvider unmodified to get genomically
correct positions) plus extract_bench_v12.py's PB-LLR computation (also
imported unmodified). Nothing here modifies any frozen/production file; all
statistics are computed by calling into binomial_baseline.py /
quality_error_model.py / cascade.py exactly as those files already define
them. This script only exists to run those frozen functions against a new
region (the existing repo has no precomputed per-locus cache for this region,
unlike the 32-44M held-out region used by experimental/unified_happy).

Output: experimental/genotype_layer/cache/dev_region.npz with counts, labels,
positions, binomial_llr, pb_llr, pb_calls, cascade_calls, binomial_calls.
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
DATA = f"{REPO}/data/giab_hg002_real_30M"
FASTA = f"{REPO}/data/reference/chr21_full.fa"
BAM = f"{DATA}/hg002_chr21_30M.bam"
VCF = f"{DATA}/hg002_chr21_30M.vcf.gz"
BED = f"{DATA}/hg002_chr21_30M_highconf.bed"

REGION_START, REGION_END = 30_000_000, 30_469_999
PB_BLOCK = 200_000

started = time.perf_counter()

# --- counts + labels + positions (exactly recover_full_positions.py's method)
counts_provider = PileupCountsProvider(
    fasta_path=FASTA, bam_path=BAM, vcf_path=VCF, contig="chr21",
    region=(REGION_START, REGION_END), seq_len=64, high_confidence_bed=BED,
)
counts_list, labels_list, pos_list = [], [], []
with counts_provider:
    windows = counts_provider._build_windows()
    for w in windows:
        ref = counts_provider._fetch_reference(w)
        counts_list.append(counts_provider.window_counts(w, ref))
        labels_list.append(np.asarray(counts_provider._fetch_labels(w), dtype=np.int64))
        pos_list.append(np.arange(w.start, w.end, dtype=np.int64))

counts = np.concatenate(counts_list).astype(np.float64)
labels = np.concatenate(labels_list)
positions = np.concatenate(pos_list)
print(f"counts extracted: {labels.size} loci, {(labels == 1).sum()} SNP-labeled "
      f"({time.perf_counter() - started:.1f}s elapsed)")

# --- reads (needed for PB LLR), same region/tiling
reads_provider = ReadLevelPileupProvider(
    fasta_path=FASTA, bam_path=BAM, vcf_path=VCF, contig="chr21",
    region=(REGION_START, REGION_END), seq_len=64, high_confidence_bed=BED,
)
reads_list = []
with reads_provider:
    for index in range(len(reads_provider)):
        reads_list.append(reads_provider[index]["reads"])
reads = np.concatenate(reads_list)
assert reads.shape[0] == counts.shape[0], (reads.shape, counts.shape)
print(f"reads extracted ({time.perf_counter() - started:.1f}s elapsed)")

# --- binomial LLR (BinomialVariantCaller, error_rate=0.01, exactly extract_bench_v12.main)
caller = BinomialVariantCaller(error_rate=0.01)
binomial_score = caller.score_counts(counts[:, 0:4], counts[:, 9])
binomial_llr = binomial_score.llr

# --- PB LLR, blocked exactly like extract_bench_v12.poisson_binomial_diagnostics
pb_llr_blocks = []
for start in range(0, counts.shape[0], PB_BLOCK):
    stop = min(start + PB_BLOCK, counts.shape[0])
    block_counts = counts[start:stop]
    reference_index = block_counts[:, 9].astype(int)
    evidence = extract_quality_evidence(reads[start:stop], reference_index)
    alt_index, k, _ = candidate_alt(block_counts[:, 0:4], reference_index)
    result = poisson_binomial_llr(evidence, k.astype(np.int64), alt_index=alt_index)
    pb_llr_blocks.append(result["llr"])
pb_llr = np.concatenate(pb_llr_blocks).astype(np.float64)
print(f"PB LLR computed ({time.perf_counter() - started:.1f}s elapsed)")

# --- frozen two-stage cascade (binomial <-> PB), exactly as
# experimental/unified_happy/build_vcfs.py's docstring specifies:
# route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF) -> pb_calls else binomial_calls
binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
pb_calls_raw = pb_llr >= FROZEN_PB_THRESHOLD
routed = np.abs(binomial_llr - FROZEN_BINOMIAL_THRESHOLD) <= FROZEN_ROUTER_CUTOFF
cascade_calls_raw = np.where(routed, pb_calls_raw, binomial_calls)

np.savez_compressed(
    f"{REPO}/experimental/genotype_layer/cache/dev_region.npz",
    counts=counts.astype(np.float32), labels=labels, positions=positions,
    binomial_llr=binomial_llr, pb_llr=pb_llr, routed=routed,
    binomial_calls=binomial_calls, pb_calls=pb_calls_raw, cascade_calls=cascade_calls_raw,
)
print(f"Saved. binomial_calls={binomial_calls.sum()} pb_calls={pb_calls_raw.sum()} "
      f"cascade_calls={cascade_calls_raw.sum()} routed={routed.sum()} "
      f"total_time={time.perf_counter() - started:.1f}s")
