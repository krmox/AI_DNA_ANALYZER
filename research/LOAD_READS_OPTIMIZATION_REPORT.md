# `load_reads` optimisation report

Scope: reduce wall-clock time of `read_level_pileup.load_reads` with zero scientific difference. No threshold, router cutoff, PB threshold, Method C parameter, truth set or protocol was touched. Nothing is committed.

Provenance note. The first version of this work (and the uncommitted HG005-fix worktree it lived in) was lost when the machine rebooted and cleared `/tmp`. Everything below was **re-implemented and re-measured** in `.claude/worktrees/opt-chr20` (on `/home`). The re-run reproduced the equivalence hashes of the lost run exactly (e.g. tensor sha prefixes `bb62a40c…`, `cd65bc97…`, `45edfa27…`, `92cf4512…`). Absolute times differ from the lost run because the desktop load on this i7-6700 (4 cores / 8 threads) differed; ratios agree.

## HYPOTHESIS

`load_reads` is CPU-bound in one Python thread, dominated by per-read Python work (pysam object/property access, per-read hashing, row arithmetic), not by I/O. If so: asyncio and Python threads cannot help it (GIL), multiprocessing helps in proportion to physical cores, and only moving the per-read work out of Python (C over htslib) gives an order-of-magnitude gain. Candidate ideas were not assumed to work; each was benchmarked.

## METHOD

- Inputs: HG005 30x chr1:1.000-1.400 Mb (8 x 50 kb chunks) as the main matrix; HG002 15x chr21:32.0-32.4 Mb as a second depth; identical BAM, reference, region, chunking and machine for every candidate; fresh process per repetition; 3 repetitions (min / median / max); CPU cores = (user+sys)/wall; peak RSS sampled over the whole process tree.
- Every run hashes the full tensor and labels (sha256) and compares with the baseline: a candidate that changes one byte is invalid.
- Baseline = the original code path: pure-Python `window_reads` and the original O(#BED intervals) `_is_confident` scan.
- Pipeline-level modes (async, threads, processes) on HG005 chr1:1.0-1.8 Mb, 8 x 100 kb chunks, stages: `load_counts` -> `load_reads` -> Binomial + PB.

## A1. Profile of the current production `load_reads` (baseline)

HG005 30x, chr1:1,000,000-1,060,000: 742 windows, 47,488 loci, 1,094,632 read-locus entries (57.4 M for the 3 Mb region; mean depth 23.9, max 87; 3,175 loci exceed the 48-read cap). Unprofiled wall 10.76 s, CPU 10.68 s => **1.0 core, no wait**.

| Component (cumulative decomposition on identical windows) | Seconds | Share of `window_reads` (10.98 s) |
|---|--:|--:|
| Window building (BED scan + reference) | 0.08 | (outside `window_reads`) |
| pysam/htslib column iteration only | 1.90 | 17.3 % |
| + `PileupRead` object creation | 0.44 | 4.0 % |
| + `.alignment` (AlignedSegment) and flag/MAPQ attributes | 0.77 | 7.0 % |
| + `indel`, `query_length`, `get_tag("NM")` | 0.33 | 3.0 % |
| + `query_position`, `query_qualities[..]`, `query_sequence[..]` | 1.24 | 11.3 % |
| read-name hash (`query_name` + BLAKE2b) | 1.25 | 11.4 % |
| everything else in Python (row arithmetic, `min`/`max`/`round`, list building, `np.asarray`, `argsort`, empty-locus allocation, `np.stack`) | 5.05 | 46.0 % |

(An earlier, heavier-loaded session measured the same region at 17.1 s with the same proportions; a `cProfile` there put `np.argsort` at 0.7 %, `np.asarray` at 3.4 %, zero-fill allocation at about 2 %.)

- **I/O**: evicting the BAM from the page cache (`posix_fadvise DONTNEED`) changed nothing: python 18.49 s cold vs 17.46 s warm; native 0.461 s vs 0.431 s. Not I/O-bound.
- **GIL / CPU**: about 85 % of the time is interpreter work; only the 17 % htslib pileup is C. The process never exceeds 1.05 cores. CPU-bound, GIL-bound.
- **Sorting, hashing, allocation, serialization**: sorting < 1 %, hashing 11 %, allocation about 2 %, serialization not present in a single process; in multi-process runs the worker returns one array set per chunk (about 30 MB per 500 kb chunk in the real driver).
- **Independence**: windows and chunks are independent (each 64-bp window is its own pileup pass); results do not depend on chunk order or thread count.
- **`_is_confident` scan** (found by profiling the native path): the original `any()` over every BED interval per window costs 0.39-0.5 ms per query on chr20 (10,192 intervals), growing with chromosome position; after the native change it was 21 % of the remaining `load_reads` time and it also affects `load_counts`.
- Verdict: **CPU-bound, Python/GIL-bound, single thread.** Not I/O-bound, not allocation-bound, not serialization-bound.

## RESULT

### A2 / A3 candidate matrix (`load_reads` only; HG005 30x 400 kb; medians of 3, spread = min-max)

| Config | Median s | Min-max s | Speedup | CPU cores | Peak MB | Output identical |
|---|--:|--:|--:|--:|--:|---|
| **Baseline (serial, pure Python)** | 79.95 | 79.6-80.3 | 1.00x | 1.05 | 642 | reference |
| Pure-Python micro-opt (one pileup per run of windows + name-hash memo) | 74.85 | 61.9-82.6 | 1.07x (noisy; best rep 1.29x) | 1.05 | 639 | yes on this data* |
| Python threads x4 over chunks (GIL test) | 109.14 | 106.2-112.2 | **0.73x** | 1.10 | 704 | yes |
| Python multiprocessing x4 | 27.95 | 27.6-28.3 | 2.86x | 3.70 | 2,037 | yes |
| Python multiprocessing x8 | 18.22 | 18.0-18.2 | 4.39x | 5.62 | 3,250 | yes |
| Native C only, original BED scan | 3.77 | 3.60-3.78 | 21.2x | 1.05 | 641 | yes |
| **Native C + bisect BED lookup (serial)** | **2.98** | 2.72-3.06 | **26.9x** | 1.05 | 642 | yes |
| Native + 4 threads inside the C call | 1.31 | 1.28-1.35 | 61.2x | 2.54 | 654 | yes |
| Native + Python threads x4 over chunks | 1.05 | 1.04-1.11 | 76.4x | 2.93 | 670 | yes |
| Native + multiprocessing x4 | 1.15 | 1.09-1.22 | 69.6x | 3.23 | 1,966 | yes |
| Native + multiprocessing x8 | 0.95 | 0.85-1.01 | 83.8x | 4.16 | 3,142 | yes |
| HG002 15x, baseline | 89.63 | 89.1-90.2 (n = 2) | 1.00x | 1.04 | 664 | reference |
| HG002 15x, native serial | 2.56 | 2.53-2.62 | **35.1x** | 1.06 | 666 | yes |
| HG002 15x, native + multiprocessing x4 | 0.94 | 0.85-1.12 | 95.9x | 3.39 | 1,936 | yes |

\* The micro-opt candidate runs one pileup pass over a run of windows instead of one per window. That is **not window-exact in general** (see below), so it was never a production candidate; it is included only to measure what pure-Python tuning can give: almost nothing.

Notes: at native speeds the region takes about 1-3 s, so the threaded/process rows are dominated by start-up and (for processes) pickling the returned 19 MB tensors per 50 kb chunk; they understate real scaling and should not be over-read. Multi-process peak RSS is the sum over workers.

### A2 async / prefetch, threaded prefetch, multiprocessing at pipeline level (`load_counts` + `load_reads` + Binomial + PB; HG005 chr1:1.0-1.8 Mb, 8 x 100 kb chunks)

| Mode | `load_reads` | Median wall s (n) | Speedup vs serial-Python | CPU cores | Peak MB |
|---|---|--:|--:|--:|--:|
| serial | pure Python (OLD) | 271.7 (2) | 1.00x | 1.05 | 1,117 |
| **async prefetch** (asyncio, next chunk loads while PB runs) | pure Python | 224.8 (2) | **1.21x** | 1.35 | 1,131 |
| threads x4 over chunks | pure Python | 246.4 (1) | 1.10x | 1.55 | 2,474 |
| processes x6 | pure Python | 78.4 (2) | 3.47x | 4.33 | 5,548 |
| serial | native | 99.9 (3) | 2.72x | 1.06 | 1,080 |
| async prefetch | native | 94.7 (3) | 2.87x (1.05x over native serial) | 1.13 | 1,152 |
| threads x4 | native | 42.9 (3) | 6.33x | 3.78 | 2,443 |
| processes x4 | native | 38.2 (3) | **7.11x** | 3.17 | 3,841 |
| processes x6 | native | 39.3 (3) | 6.91x | 3.16 | 5,406 |

All nine configurations produced the **same output hash** for `counts`, `labels`, `positions`, `pb_llr`, `binomial_llr`.

Reading it honestly:
- **Async on `load_reads` itself: expected 1.00x, not separately timed.** There is no I/O to hide (cold cache = warm cache) and the work is GIL-bound Python, so an `asyncio`/thread prefetch of `load_reads` alone cannot speed it up; the measured proxy is Python threads over chunks, which made `load_reads` alone 0.73x (contention). The only async measurement is the pipeline-level one below.
- **Async at pipeline level is a real but modest 1.21x** with Python reads, and it is *not* a speed-up of `load_reads`: it comes from overlapping `load_reads` (Python) with PB (numpy, which releases the GIL). Both stages become slower in isolation (`reads` 175-190 s -> 204-210 s; `pb` 83-86 s -> 191-198 s summed) but overlap wins. With native reads there is nothing left to overlap: 1.05x.
- **Multiprocessing works** (2.9x/4.4x on `load_reads` with 4/8 workers, 3.5x on the pipeline with 6) and scales with the 4 physical cores: 6 workers were no faster than 4 on this CPU.
- **Native C is required for a large `load_reads` gain**: 27-35x serial, an order of magnitude beyond anything Python-level.

## A4. Exact-equivalence gate (OLD vs NEW)

| # | Item | Evidence | Differences |
|---|---|---|---|
| 1 | positions | `_is_confident` bisect vs original scan: 253,707 real-BED queries on the HG005, HG002-chr21 and chr20 BEDs (incl. every interval boundary +/-1 window) + 320,000 queries over 200 randomised (overlapping, unsorted, adjacent) interval sets; chr20 chunk gate (arrays incl. `positions`) | 0 |
| 2 | read IDs / order | order is encoded by the BLAKE2b-8 sort (names are not stored); C hash equals `hashlib.blake2b(digest_size=8)` on 262 lengths x 5 random inputs; identical tensors imply identical order | 0 |
| 3 | read tensors | real data, byte for byte (pure-Python vs native `load_reads`): HG005 chr1:1.0-1.5 Mb (391,488 loci, 9,218,195 valid rows, includes loci at the 48-read cap), HG005 chr1:3.5-4.0 Mb (437,824 loci, 10,488,246 rows), HG002 15x chr21:32.0-32.5 and 40.0-40.5 Mb (470,464 and 455,296 loci), chr20:46.0-46.5 Mb: **0 mismatching cells** of 6.7e8+ cells; plus the 17-test suite | 0 |
| 4 | counts | chr20 chunks (`counts` array, 980,928 loci), all pipeline-mode runs share one hash | 0 |
| 5 | labels | every comparison above | 0 |
| 6, 7 | Binomial LLR, PB LLR | chr20 chunks `binomial_llr`, `pb_llr`, `k_full`, `k_retained`, `n_counted`, `depth` bytewise equal old-path vs production cache; 9 pipeline modes one hash | 0 |
| 8 | router mask | pure function of `binomial_llr` | 0 |
| 9-12 | final calls, genotype, GQ, VCF | pure functions of the arrays above; HG005 3 Mb VCF records identical between 1 and 4 workers; chr20 arrays identical for the original vs optimised path | 0 |

**A candidate that failed the gate and was rejected.** A first native design made one pileup pass over a whole run of windows. On a randomised BAM with repeated read names it differed from the original in 30 cells: htslib's mate-overlap quality adjustment pairs reads by name among those alive in the pileup, so a longer pass can pair reads differently from the original per-64-bp-window passes. The final design runs one fresh pass per window inside C (open BAM/index once), which is the same computation by construction. This is the reason the gate exists.

Also gated: thread count inside C (1, 2, 4) gives identical bytes; unsupported inputs (SEQ `*`, non-integer `NM` tag) make the C call raise so the wrapper falls back to the original Python path rather than diverge.

Caveat: my first session's stronger gate (all arrays and VCF records of a full 2.4 M-locus HG005 region compared with the previously validated pure-Python output) was lost with `/tmp` and was **not re-run** here (it would cost about 25 min of pure-Python extraction and the gold cache no longer exists). It is replaced by the per-chunk tensor gates on 5 real chunks and by the chr20 array gate on two 500 kb chunks against the original path.

## A5 / A6. New end-to-end performance and Amdahl

HG005 chr1:1.0-4.0 Mb (2,403,008 loci, 3,949 SNP), whole pipeline (BAM -> counts -> read tensor -> Binomial -> PB -> Method C -> VCF), new code:

| Config | Wall | Stage seconds (summed over workers) |
|---|--:|---|
| 1 worker | **309 s** (5 min 15 s incl. start-up) | counts 5, reads 15, PB 287, Method C + VCF < 0.1 |
| 4 workers | **146 s** | counts 8, reads 21, PB 423 |

Same rescue and routing as the HG005 record (1,600 routed, 142 rescuable, 140 rescued); 1 and 4 workers give identical arrays and VCF records. For reference, the earlier recorded pipeline times for this region were 1,392 s (pure Python) and 1,140 s (native counts only); those were not re-measured on today's machine load, so the same-day comparison is the 800 kb pipeline matrix above: **2.72x serial** from native reads and BED lookup, **7.1x** with 4 processes.

Amdahl (800 kb same-day stage sums, serial): OLD counts 3.7 s, reads about 183 s, PB about 85 s, Binomial 0.3 s (total 272 s). NEW counts 2.0, reads 4.8, PB 95, Binomial 0.3 (total 102 s). Reading `load_reads` as fraction p = 183/272 = 0.67 of the old runtime, the ceiling for making it free is 1/(1-p) = 3.06x; the measured 2.72x is 89 % of that ceiling, which is what Amdahl predicts. `load_reads` sped up 38x but the pipeline only 2.7x.

**New bottleneck: Poisson-binomial LLR = 93 % of runtime** (287 of 309 s on HG005; 95.9 % of stage compute on chr20). Ceiling from making PB free: about 14x. Native counts and reads together are now 6.5 % of the serial pipeline. The per-locus PB cost (0.12-0.15 ms) was not touched by this work.

chr20 (56.3 M loci, 6 workers): counts 194 s, reads 444 s, PB 15,036 s summed, wall 2,640 s. Native `load_reads` cost about 8 microseconds per locus there.

## EQUIVALENCE (summary)

Zero semantic differences were found in every comparison listed in A4. One earlier design was rejected for a real 30-cell divergence on adversarial input; the shipped design is exact by construction and tested.

## DECISION

- **Selected**: native per-window C tensor builder (`native/reads_native.c`, called from `load_reads`) + exactly-equivalent bisect BED lookup in `providers._is_confident`. Measured 26.9x (HG005 30x) and 35.1x (HG002 15x) on `load_reads`, 2.7x on the whole pipeline, byte-identical output. Falls back to the original path when the extension is absent, disabled (`AI_DNA_ANALYZER_DISABLE_NATIVE_READS=1`) or when an input cannot be reproduced.
- **Driver-level `--workers N`** in `extract_bench_v12.py` (chunks concatenated in order; byte-identical) is the simplest architecture that gives whole-pipeline scaling (7.1x at 4 workers); it helps PB as well, which C-level `load_reads` work cannot. Optional in-C threads (`AI_DNA_ANALYZER_READ_THREADS`, default 1) are available but redundant with `--workers`.
- **Not selected**: asyncio/prefetch (1.00x on `load_reads`; 1.21x at pipeline level only with Python reads, 1.05x with native), Python threads (0.73x / 1.10x), pure-Python micro-optimisation (1.07x median), and Python multiprocessing alone (2.9-4.4x, superseded by native).
- Not done: PB acceleration (out of scope, and it changes numerics). It is the next target.

## Files changed (all untracked/uncommitted, in `.claude/worktrees/opt-chr20`)

Modified: `read_level_pileup.py`, `providers.py`, `pileup_counts.py` (native counts path + `return_positions`, rebuilt from the previously validated design), `extract_bench_v12.py` (`process_chunk`, `--workers`, skip chunks with no confident windows, true positions). New: `native/reads_native.c`, `native/build.sh` (builds both extensions), `test_native_reads_equivalence.py`, `experimental/load_reads_opt/*` (benchmarks, equivalence scripts, results), `experimental/chr20_validation/*`, `bench_v20_chr20.py`.

## Tests

- `test_native_reads_equivalence.py`: 17 passed (8 fuzz BAMs, small max_reads, BED gaps, 2 thread counts, fallback, missing-SEQ behaviour, BLAKE2b, 2 real-data slices).
- Existing suites `test_bench_v12`, `test_cascade`, `test_cheap_router`, `test_robustness_benchmark`, `test_bench_v13`: 79 passed. (`test_providers.py` cannot be collected because of a missing `bimamba_variant_caller` module; unrelated and pre-existing.)
- Not restored: `test_native_pileup_equivalence.py` and `test_hg005_extraction_integrity.py` (lost with `/tmp`; I never had their source in this session).

## Limitations

- Desktop load on the benchmark machine (about 1-2 cores) adds noise; the micro-opt row shows it (61.9-82.6 s). HG002 baseline has n = 2.
- One dataset family (Illumina novoalign 2x250, no supplementary/secondary reads). The per-window design is exact by construction for repeated read names, but that path is only exercised by synthetic data.
- Multi-process speedups are limited to about 4x-5x by the CPU's 4 physical cores.
- The 48-read cap path is verified on real data (HG005 has loci at the cap) and synthetic data.
