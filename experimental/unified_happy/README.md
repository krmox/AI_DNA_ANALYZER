# Unified hap.py benchmark (AI PB-only / AI-Cascade vs DeepVariant / Clair3 / GATK)

This directory fixes the methodological flaw in `experimental/head_to_head/`:
AI PB-only and AI-Cascade were previously scored with the frozen internal
`accuracy_block()` evaluator against a different truth denominator (16,597
loci) than the three external callers (16,898 loci via hap.py), on a
non-genomic-coordinate scoring frame. Here, all five callers are scored by
the identical hap.py 0.3.15 (xcmp engine) run against the identical truth
VCF / high-confidence BED / reference / region.

## What changed vs `experimental/head_to_head/`

Nothing about the three external callers. Their VCFs and hap.py outputs
(`experimental/head_to_head/happy/{deepvariant,clair3,gatk}/`) are reused
as-is, not recomputed.

What's new: AI PB-only and AI-Cascade were converted into real per-locus VCFs
with true chr21 genomic coordinates (see `SANITY_CHECK.md` for why the
existing `extract_bench_v12.py` per-chunk `positions` array could not be
trusted for this, and how true positions were instead recovered from
`pileup_counts.PileupCountsProvider`'s own window boundaries), then run
through the same hap.py 0.3.15 invocation as the external callers.

## Root cause of the original flaw (see `SANITY_CHECK.md` for full derivation)

`extract_bench_v12.py` computes `positions = chunk_start + arange(labels.size)`
per 500kb chunk. This is only correct if no 64bp window inside that chunk was
dropped by BED confidence. `providers.py::_build_windows` drops whole windows
(never partial rows within a kept window) when `_is_confident()` is false or
the window is mostly `N`. Once any window in a chunk is dropped, every
downstream position in that chunk silently decouples from the true genomic
coordinate (no visible discontinuity — the naive array keeps incrementing by
1) by the cumulative width of all windows dropped so far. Confirmed
empirically in `SANITY_CHECK.md`.

## How AI PB-only / AI-Cascade VCFs were built

1. `recover_full_positions.py` replays the exact same chunk loop
   (`range(32_000_000, 44_000_000, 500_000)`) and provider construction
   `extract_bench_v12.py` uses (via `pileup_counts.PileupCountsProvider`,
   imported unmodified), but records `arange(window.start, window.end)` per
   *kept* window instead of trusting the chunk-relative naive formula. This
   reproduces the exact same flattened locus ordering `extract_bench_v12.py`
   itself produces (same provider, same chunking, same window tiling), which
   is verified by an exact label-array equality check against the existing
   frozen `experimental/head_to_head/ai_cascade_per_locus.npz` (see
   `build_vcfs.py`).
2. `build_vcfs.py` zips the recovered true positions with the **existing,
   unmodified** `pb_calls` / `cascade_calls` boolean arrays already computed
   in `experimental/head_to_head/ai_cascade_per_locus.npz` by
   `experimental/head_to_head/analyze_cascade.py` (which itself only imports
   `route_mask`/`accuracy_block` from `robustness_benchmark.py` and the frozen
   thresholds from `cascade.py` — nothing here recomputes or refits any
   threshold or routing decision). No new evidence, threshold, or routing
   logic was introduced; only the position/VCF-serialization layer is new.
3. For every locus flagged `True` in `pb_calls` (AI PB-only) or
   `cascade_calls` (AI-Cascade), a VCF record is emitted:
   - `REF`: the reference base baked into the count matrix's
     `reference_index` channel (independently cross-validated against a
     fresh `pysam.FastaFile` read in `SANITY_CHECK.md`).
   - `ALT`: the majority non-reference base among the count matrix's raw
     `count_A/C/G/T` channels at that locus (the same "which base is actually
     observed" signal the frozen `_pileup_locus` consensus logic in
     `providers.py` uses for its own VAF-threshold decision, but read off the
     provider's own count matrix here rather than reimplemented).
   - `GT`: **always `0/1`** (heterozygous). **Limitation, stated plainly**:
     the frozen PB-threshold / binomial-threshold decision logic is a binary
     variant/no-variant call over an LLR, not a genotype likelihood model —
     it has no zygosity signal to report. `0/1` is a documented default, not
     a derived genotype. This is a genuine methodological asymmetry vs.
     DeepVariant/Clair3/GATK, which do emit genotype-likelihood-based GT, and
     is called out in the final report's Problems section.
   - SNP-only: loci with an insertion/deletion truth label are excluded from
     the eval frame identically to how `analyze_cascade.py`'s `frame` mask
     already excluded them (`labels in {NORMAL, SNP}` only) — this mirrors,
     not extends, the frozen evaluator's own scope. The frozen evidence
     pipeline being used here (Binomial/PB LLR over substitution counts) does
     not natively emit indel calls at all, so no separate indel-filtering
     policy decision was needed.

## Files

- `sanity_check.py`, `SANITY_CHECK.md` — Step 1.
- `recover_full_positions.py`, `recovered_full.npz` — Step 2 position recovery.
- `build_vcfs.py`, `ai_pb_only.vcf(.gz)`, `ai_cascade.vcf(.gz)` — Step 2/3 VCF construction + normalization.
- `run_happy.sh`, `happy/ai_pb_only/`, `happy/ai_cascade/` — Step 4.
- `UNIFIED_HAPPY_BENCHMARK.csv` / `.md` — Step 6 final table.
- `DISAGREEMENTS.csv` — PB vs Cascade locus-level differences (Step 6).
- `runtime.csv` — Step 5.
- `COMMANDS.md`, `environment.txt`, `checksums.txt` — Step 7 reproducibility.
