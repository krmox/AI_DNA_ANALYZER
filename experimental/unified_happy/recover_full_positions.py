"""Recompute counts/labels/TRUE-positions across chr21:32-44M using the same
chunking as extract_bench_v12.py, but replaying window starts explicitly so
positions are genomic-coordinate-correct (see SANITY_CHECK.md).

Imports pileup_counts.PileupCountsProvider unmodified (frozen code). Does not
call load_reads / poisson_binomial (not needed: pb_llr/binomial_llr/calls are
reused as-is from the existing experimental/head_to_head/ai_cascade_per_locus.npz,
we only need to recover correct positions + REF/ALT to align with them).
"""
import sys
import time
sys.path.insert(0, "/home/mark/Documents/Projects/AI_DNA_ANALYZER")

import numpy as np

from pileup_counts import PileupCountsProvider

FASTA = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/reference/chr21_full.fa"
BAM = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_15x.bam"
VCF = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz"
BED = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed"

REGION_START, REGION_END = 32_000_000, 44_000_000
CHUNK_BP = 500_000

counts_blocks, labels_blocks, pos_blocks = [], [], []
started = time.perf_counter()

for chunk_start in range(REGION_START, REGION_END, CHUNK_BP):
    chunk_end = min(chunk_start + CHUNK_BP, REGION_END)
    provider = PileupCountsProvider(
        fasta_path=FASTA, bam_path=BAM, vcf_path=VCF, contig="chr21",
        region=(chunk_start, chunk_end), seq_len=64, high_confidence_bed=BED,
    )
    with provider:
        windows = provider._build_windows()
        c_list, l_list, p_list = [], [], []
        for w in windows:
            ref = provider._fetch_reference(w)
            c_list.append(provider.window_counts(w, ref))
            l_list.append(np.asarray(provider._fetch_labels(w), dtype=np.int64))
            p_list.append(np.arange(w.start, w.end, dtype=np.int64))
        if c_list:
            counts_blocks.append(np.concatenate(c_list))
            labels_blocks.append(np.concatenate(l_list))
            pos_blocks.append(np.concatenate(p_list))
    elapsed = (time.perf_counter() - started) / 60
    print(f"{chunk_start}-{chunk_end}: {len(windows)} windows | cumulative {elapsed:.1f} min", flush=True)

counts = np.concatenate(counts_blocks).astype(np.float32)
labels = np.concatenate(labels_blocks)
positions = np.concatenate(pos_blocks)

assert labels.size == positions.size == counts.shape[0]
print(f"TOTAL loci: {labels.size}, SNP labels: {int((labels==1).sum())}")

np.savez_compressed(
    "/home/mark/Documents/Projects/AI_DNA_ANALYZER/experimental/unified_happy/recovered_full.npz",
    counts=counts, labels=labels, positions=positions,
)
print("Saved recovered_full.npz")
