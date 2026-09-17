# Step 1 — Coordinate sanity check

Script: `sanity_check.py` (imports `pileup_counts.PileupCountsProvider` unmodified;
does not edit any frozen file).

## Method

Sub-region tested: `chr21:32,000,000-32,100,000` (100kb, well inside the
benchmark region and inside `data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed`
for most of its span, with a real BED-driven window drop nearby — see below).

`PileupCountsProvider(..., region=(32_000_000, 32_100_000), seq_len=64,
high_confidence_bed=<highconf.bed>)` was constructed directly (same class
`load_counts` in `pileup_counts.py` uses internally). `provider._build_windows()`
was called to get the real list of kept 64bp windows (`GenomicWindow.start`
values are literal 0-based genomic coordinates on `chr21`).

Two position arrays were then built over the same window/label/count arrays:

1. **True position array**: `concat([arange(w.start, w.end) for w in windows])`
   — one true genomic 0-based coordinate per emitted locus.
2. **Naive formula** exactly as `extract_bench_v12.py` line ~234 computes it:
   `chunk_start + arange(labels.size)`.

## Finding: root cause CONFIRMED

- Of 1562 possible 64bp windows in this 100kb sub-region, only **1393 were kept**;
  **169 were dropped** by `_is_confident()` (whole-window BED masking in
  `providers.py::_build_windows`, not partial-row masking — confirmed by code
  reading of `providers.py` lines 550-599: `_is_confident` gates entire windows,
  and `window_counts`/`_fetch_labels` never drop rows inside a kept window since
  ``len(window)`` rows are always emitted per window whether or not there is
  pileup coverage).
- First dropped window observed at genomic start `32,005,568`.
- **`naive_pos == true_pos` is `False`.** First divergence is at array index
  `5568` (`5568 = 32_005_568 - 32_000_000`, exactly the offset of the first
  dropped window): `naive=32,005,568` but `true=32,005,760` (a 192bp = 3-window
  shift, since 3 consecutive windows were dropped there:
  `32005568, 32005632, 32005696`). This numerically confirms
  `ERROR_ANALYSIS.csv`'s claim and the task's stated root cause: the naive
  per-chunk formula desyncs from true genomic coordinates by the cumulative
  width of every window dropped so far in that chunk, silently (no
  discontinuity appears in the naive array itself — it keeps incrementing by
  1 — so the corruption is invisible without an independent check like this
  one).

## Finding: TRUE position reconstruction is correct (positive control)

Using `true_pos = concat([arange(w.start, w.end) ...])`, 12 loci spread across
the sub-region (indices 0, 8104, 16209, ..., 89151) were independently checked:

```
idx=0     true_pos_0based=32000000 vcf_pos_1based=32000001 ref_fasta=A ref_from_counts_matrix=A match=True in_bed=True label=0
idx=8104  true_pos_0based=32008488 vcf_pos_1based=32008489 ref_fasta=A ref_from_counts_matrix=A match=True in_bed=True label=0
idx=16209 true_pos_0based=32016593 vcf_pos_1based=32016594 ref_fasta=G ref_from_counts_matrix=G match=True in_bed=True label=0
idx=24313 true_pos_0based=32032697 vcf_pos_1based=32032698 ref_fasta=C ref_from_counts_matrix=C match=True in_bed=True label=0
idx=32418 true_pos_0based=32040866 vcf_pos_1based=32040867 ref_fasta=C ref_from_counts_matrix=C match=True in_bed=True label=0
idx=40523 true_pos_0based=32049099 vcf_pos_1based=32049100 ref_fasta=A ref_from_counts_matrix=A match=True in_bed=True label=0
idx=48627 true_pos_0based=32057587 vcf_pos_1based=32057588 ref_fasta=T ref_from_counts_matrix=T match=True in_bed=True label=0
idx=56732 true_pos_0based=32065692 vcf_pos_1based=32065693 ref_fasta=T ref_from_counts_matrix=T match=True in_bed=True label=0
idx=64837 true_pos_0based=32073797 vcf_pos_1based=32073798 ref_fasta=A ref_from_counts_matrix=A match=True in_bed=True label=0
idx=72941 true_pos_0based=32082669 vcf_pos_1based=32082670 ref_fasta=C ref_from_counts_matrix=C match=True in_bed=True label=0
idx=81046 true_pos_0based=32091350 vcf_pos_1based=32091351 ref_fasta=A ref_from_counts_matrix=A match=True in_bed=True label=0
idx=89151 true_pos_0based=32099967 vcf_pos_1based=32099968 ref_fasta=G ref_from_counts_matrix=G match=True in_bed=True label=0
```

For every checked locus: `samtools`/`pysam.FastaFile.fetch` reference base
(fetched 0-based half-open, i.e. `fetch(pos0, pos0+1)`) **matches exactly**
the `reference_index` column baked into the count matrix by
`PileupCountsProvider._count_locus` (independent cross-check: the count
matrix's own reference-base channel, computed from the same FASTA slice at
build time, agrees with a fresh, separate `pysam.FastaFile` read at the
candidate position). BED containment (`in_bed`) was independently
recomputed by parsing the raw BED file, not reusing the provider's cached
regions.

Five real truth SNPs in the sub-region were independently pulled directly
from the truth VCF (`pysam.VariantFile(...).fetch("chr21", 32_000_000,
32_100_000)`, filtered to single-base REF/ALT):

```
truth_snp pos1based=32000038 pos0based=32000037 ref=A alts=('G',) found_in_true_pos_array=True label_at_idx=1
truth_snp pos1based=32000436 pos0based=32000435 ref=T alts=('C',) found_in_true_pos_array=True label_at_idx=1
truth_snp pos1based=32000664 pos0based=32000663 ref=A alts=('T',) found_in_true_pos_array=True label_at_idx=1
truth_snp pos1based=32000810 pos0based=32000809 ref=C alts=('T',) found_in_true_pos_array=True label_at_idx=1
truth_snp pos1based=32000815 pos0based=32000814 ref=C alts=('G',) found_in_true_pos_array=True label_at_idx=1
```

Every one of these truth-VCF SNP coordinates was found at the exact matching
index in `true_pos`, AND the frozen `_fetch_labels()` output at that same
index is `1` (`LABEL_SNP`), independently confirming: (a) the true-position
array is genomically correct, (b) the frozen label-fetching logic agrees
with the truth VCF at that coordinate, and (c) 1-based/0-based conventions
are handled consistently (`pysam.VariantFile` record `.start`/`.pos` are
0-based/1-based respectively; `pysam.FastaFile.fetch(start, end)` is 0-based
half-open; all comparisons above were done in the 0-based frame and
converted to 1-based only for display).

## Verdict: PASS

The true-position reconstruction method — replaying
`provider._build_windows()` window starts and emitting `arange(w.start,
w.end)` per kept window, in the same chunk order `extract_bench_v12.py`
uses — is empirically verified correct on this sub-region, including across
a real BED-driven window drop. The naive `chunk_start + arange(labels.size)`
formula used by `extract_bench_v12.py`'s `positions` output is confirmed
broken (diverges silently after the first dropped window in a chunk), which
matches the root cause hypothesis in the task brief and in
`experimental/head_to_head/ERROR_ANALYSIS.csv`.

**Proceeding to Step 2** using the true-position reconstruction method
(implemented in `recover_full_positions.py`), not the naive formula.

Note on `(a)`/`(b)` from the task brief: the answer is **rows ARE dropped**
(whole 64bp windows, not individual loci) and true positions **are**
recoverable — from `GenomicWindow.start`, which `PileupCountsProvider`
already computes internally and exposes via the (private but unmodified)
`_build_windows()` method — without touching any frozen production file.
