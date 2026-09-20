> **SUPERSEDED IN PART (2026-09-20).** Historical methods note; `PAPER_MANUSCRIPT.md` §2 is authoritative (includes the native `load_reads` backend and the counts-backend limitation).

# Methods

## 1. Architecture

```
BAM/CRAM (Illumina, GRCh38-aligned)
        |
   Python/pysam pileup extraction (pileup_counts.py / providers.py)
        |  -> per-locus integer counts (A/C/G/T, depth, quality sums,
        |     reference index), tiled in 64bp windows
        v
   Binomial screening caller (binomial_baseline.py)
        |  LLR against a fixed per-base error rate; threshold = 7.0
        v
   Cheap router (cascade.py::route_mask)
        |  routes a locus to the Poisson-Binomial arm iff
        |  |binomial_LLR - 7.0| <= 5.411872376933351
        v
   Poisson-Binomial caller (quality_error_model.py)      [only for routed loci]
        |  per-read quality-aware error model; threshold = 10.5
        v
   Cascade decision = PB call if routed, else Binomial call
        |
   Method C genotype layer (experimental/genotype_layer/method_c_regression.py)
        |  binomial posterior argmax over {0/0, 0/1, 1/1}, eps=0.01,
        |  GQ = 10*log10(P_best/P_second) capped at 99
        v
   VCF (SNP-only, PASS)
```

All three detection thresholds and the genotype-layer formula are frozen
constants, defined once in `cascade.py` and
`experimental/genotype_layer/method_c_regression.py`, and were not changed at
any point in the work this dossier reports.

## 2. Data provider

`providers.py::GiabAlignmentProvider` tiles a requested genomic region into
non-overlapping 64bp windows, drops a window entirely if any base of it lies
outside an optional high-confidence BED or if it is mostly `N`, and exposes
per-locus read counts (`pileup_counts.py::PileupCountsProvider`) and
per-read tensors (`read_level_pileup.py::ReadLevelPileupProvider`) under
identical read-admission filters: `min_mapping_quality=20`,
`min_base_quality=13`, duplicate/QC-fail/unmapped/secondary/supplementary
reads excluded. Genomic coordinates are 0-based half-open throughout
(documented convention in `providers.py`'s module docstring), reconciling
FASTA/BAM 0-based coordinates against VCF's 1-based convention exactly once,
at the boundary (`record.start`, never `record.pos`).

## 3. Evaluation methodology

Two distinct evaluators are used across this project's history and are
**never treated as interchangeable** in this dossier:

- **Internal (`frozen_locus_evaluator`)**: `robustness_benchmark.py`'s
  `accuracy_block()`/`prf()`. Pure row-index boolean agreement between a
  call mask and a truth-label mask built from the same window tiling. No
  REF/ALT/GT semantics; fast, used for the pre-registered `bench_v13`-`v19`
  sweeps (many depths x regions x samples).
- **External (`hap.py_0.3.15_xcmp`)**: GA4GH-standard genotype-aware VCF
  comparison, used whenever a VCF was actually written and compared against
  the truth VCF (unified_happy HG002 benchmark, external-caller comparison,
  and the HG005 post-fix result this dossier centers on).

Every table in `BENCHMARK_MASTER.csv` carries an explicit `evaluator`
column; rows using different evaluators are marked `NOT DIRECTLY COMPARABLE`
in their notes.

## 4. Region and sample selection protocol

Regions for `bench_v14` onward were selected mechanically from public
annotation (GC content, mappability, GIAB difficulty strata) by a dedicated
`select_v1N_region*.py` script *before* the frozen cascade was ever run on
them, and the independence of this selection from any accuracy outcome is
verified programmatically per-run (`independence_audit.json` in each
`results/bench_v1N/` directory). The HG005 region
(chr1:1,000,001-4,000,000) was selected the same way in an earlier session
(documented in `experimental/stress_test/HG005_MAX_FEASIBLE_STRESS_TEST.md`
Sec. 2) and was **not** moved, expanded, or re-selected after this session's
extraction fix, per explicit instruction.

## 5. Statistical methodology

Where a paired comparison between Binomial-only, PB-only, and Cascade is
reported with a confidence interval (`bench_v13`-`v19`), it uses a paired
bootstrap over scored loci (10,000 resamples, fixed seed `20260812`),
reporting `delta_f1`, its 95% CI, and a two-sided p-value. No bootstrap was
re-run for the HG005 post-fix result in this session (sample size at the
per-experiment level is 1: one region, one sample); this is stated
explicitly as a limitation rather than glossed over.

## 6. Extraction correctness fix (this session)

See `experimental/stress_test/HG005_EXTRACTION_BUG_ROOT_CAUSE.md` for full
detail. Summary: `extract_bench_v12.py` previously derived each extracted
locus's genomic coordinate as `chunk_start + row_index`, which silently
desyncs the moment a 64bp window is dropped by the BED/N-content filter
described in Sec. 2. `pileup_counts.py::load_counts` was given an opt-in
`return_positions=True` path that returns each row's true coordinate
(`window.start + offset`, read directly from the surviving window object),
and `extract_bench_v12.py` was updated to consume it at its 3 coordinate
call sites. The fix was independently verified (full-population 0-mismatch
FASTA cross-check, byte-identical repeated extraction, 0 VCF REF mismatches)
and proven not to alter `counts`/`labels` arrays on real HG004 data
(byte-identical, pre- vs post-fix).

## 7. C/htslib extraction acceleration (this session, integrating a
prior worktree's validated experiment)

A separate, previously-completed worktree
(`.claude/worktrees/agent-a4cedfdeffb76e3dd/experimental/performance/`)
had implemented and rigorously validated a C extension
(`Python.h` + htslib's C API, `native/pileup_native.c`) that replaces
`PileupCountsProvider.window_counts`/`_count_locus`'s per-64bp-window
`pysam.pileup()` calls with a single streaming `bam_mplp_*` pass per
requested region, proven bit-exact against the pure-Python reference over a
full 12 Mb HG002 benchmark region (0/109,813,120 cell mismatches) and
byte-identical downstream VCF/hap.py output — but that work was never
integrated into the production codebase. This session located it, rebuilt
it from source (`native/build.sh`, portable, no hard-coded paths, no
external network dependency beyond the one-time header download when a
system `python3-devel` package is unavailable), independently re-verified
its equivalence via the actual `load_counts()` production entry point on
both HG002 and HG005 data, and wired it in as the default extraction path
(`pileup_counts.NATIVE_AVAILABLE`), with an automatic fallback to the
pure-Python path in any environment without a working C toolchain. Full
detail: `experimental/performance/CHTSLIB_INTEGRATION_REPORT.md`,
`PIPELINE_INTEGRITY_AUDIT.md`'s correction note, and
`test_native_pileup_equivalence.py` (12-category permanent regression
suite).

## 8. Performance methodology

All timing figures in `PERFORMANCE_RESULTS.csv` are wall-clock, single
invocation unless stated otherwise (`htslib_thread_scaling` and
`multiprocess_chunking` rows report `n=3` repetitions per configuration from
their source benchmark, which reported mean and stdev; this dossier reports
the mean). No performance figure is extrapolated to whole-genome scale
except the rows explicitly labeled `EXTRAPOLATION` in that CSV. A C/htslib
extraction extension (`native/pileup_native.c`) is integrated into
production `pileup_counts.load_counts` as of this session, proven bit-exact
against the pure-Python reference (see Sec. 6 and
`experimental/performance/CHTSLIB_INTEGRATION_REPORT.md`); its measured
speedup is reported for the extraction stage specifically and is not
conflated with full-pipeline wall time, which also depends on the
unaccelerated `read_level_pileup.load_reads` stage.
