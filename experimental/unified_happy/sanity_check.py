"""Step 1 sanity check: verify true-position recovery for load_counts windows.

Imports frozen code (providers.py, pileup_counts.py) unmodified. Does not edit
any production file. Writes findings to SANITY_CHECK.md.
"""
import sys
sys.path.insert(0, "/home/mark/Documents/Projects/AI_DNA_ANALYZER")

import numpy as np
import pysam

from pileup_counts import PileupCountsProvider

FASTA = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/reference/chr21_full.fa"
BAM = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_15x.bam"
VCF = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz"
BED = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed"

SUB_START, SUB_END = 32_000_000, 32_100_000  # 100kb sub-region, several 64bp windows

provider = PileupCountsProvider(
    fasta_path=FASTA, bam_path=BAM, vcf_path=VCF, contig="chr21",
    region=(SUB_START, SUB_END), seq_len=64, high_confidence_bed=BED,
)

with provider:
    windows = provider._build_windows()
    print(f"Sub-region {SUB_START}-{SUB_END}: {len(windows)} windows kept "
          f"out of {(SUB_END - SUB_START) // 64} possible (stride=seq_len=64)")
    kept_starts = [w.start for w in windows]
    all_possible_starts = list(range(SUB_START, SUB_END - 64 + 1, 64))
    dropped_starts = sorted(set(all_possible_starts) - set(kept_starts))
    print(f"Dropped window starts (first 10): {dropped_starts[:10]}")
    print(f"Total dropped: {len(dropped_starts)}")

    # Build counts/labels + TRUE positions by replaying window starts.
    counts_list, labels_list, true_pos_list = [], [], []
    for w in windows:
        ref = provider._fetch_reference(w)
        counts_list.append(provider.window_counts(w, ref))
        labels_list.append(np.asarray(provider._fetch_labels(w), dtype=np.int64))
        true_pos_list.append(np.arange(w.start, w.end, dtype=np.int64))
    counts = np.concatenate(counts_list)
    labels = np.concatenate(labels_list)
    true_pos = np.concatenate(true_pos_list)

    # NAIVE formula from extract_bench_v12.py: chunk_start + arange(labels.size)
    naive_pos = SUB_START + np.arange(labels.size, dtype=np.int64)

print(f"\nlabels.size = {labels.size}, naive last pos = {naive_pos[-1]}, "
      f"true last pos = {true_pos[-1]}")
print(f"naive == true everywhere? {np.array_equal(naive_pos, true_pos)}")
if not np.array_equal(naive_pos, true_pos):
    first_diff = np.argmax(naive_pos != true_pos)
    print(f"First divergence at array index {first_diff}: "
          f"naive={naive_pos[first_diff]} true={true_pos[first_diff]}")

# --- Cross-checks using true_pos ---
fasta = pysam.FastaFile(FASTA)
bed_regions = []
with open(BED) as fh:
    for line in fh:
        c, s, e = line.split()[:3]
        if c in ("chr21", "21"):
            bed_regions.append((int(s), int(e)))

def in_bed(pos0):
    return any(s <= pos0 < e for s, e in bed_regions)

vcf = pysam.VariantFile(VCF)

print("\n--- Spot checks (10 loci spread across sub-region) ---")
idxs = np.linspace(0, len(true_pos) - 1, 12).astype(int)
report_lines = []
for i in idxs:
    pos0 = int(true_pos[i])  # 0-based
    ref_base = fasta.fetch("chr21" if "chr21" in fasta.references else "21", pos0, pos0 + 1).upper()
    count_row = counts[i]
    ref_idx_from_counts = int(count_row[9])  # 0=N,1=A,2=C,3=G,4=T
    tok = ("N", "A", "C", "G", "T")[ref_idx_from_counts]
    bed_hit = in_bed(pos0)
    label = int(labels[i])
    line = (f"idx={i} true_pos_0based={pos0} vcf_pos_1based={pos0+1} "
            f"ref_fasta={ref_base} ref_from_counts_matrix={tok} match={ref_base==tok} "
            f"in_bed={bed_hit} label={label}")
    print(line)
    report_lines.append(line)

# Confirm a real truth SNP in this sub-region lands at the expected true_pos index
print("\n--- Truth VCF SNPs in sub-region (first 5) ---")
snp_checks = []
count = 0
for rec in vcf.fetch("chr21", SUB_START, SUB_END):
    if rec.alts is None:
        continue
    if len(rec.ref) != 1 or any(len(a) != 1 for a in rec.alts if a):
        continue
    pos0 = rec.start  # pysam 0-based
    if pos0 in true_pos:
        idx = int(np.searchsorted(true_pos, pos0))
        matches = true_pos[idx] == pos0
        line = (f"truth_snp pos1based={rec.pos} pos0based={pos0} ref={rec.ref} "
                f"alts={rec.alts} found_in_true_pos_array={matches} label_at_idx={int(labels[idx]) if matches else 'N/A'}")
    else:
        line = f"truth_snp pos1based={rec.pos} pos0based={pos0} NOT in kept true_pos array (window dropped by BED or filtered)"
    print(line)
    snp_checks.append(line)
    count += 1
    if count >= 5:
        break

with open("/home/mark/Documents/Projects/AI_DNA_ANALYZER/experimental/unified_happy/_sanity_raw_output.txt", "w") as f:
    f.write(f"windows kept: {len(windows)} / possible {len(all_possible_starts)}\n")
    f.write(f"dropped starts (first 20): {dropped_starts[:20]}\n")
    f.write(f"naive==true: {np.array_equal(naive_pos, true_pos)}\n")
    f.write("\n".join(report_lines) + "\n\n")
    f.write("\n".join(snp_checks) + "\n")
