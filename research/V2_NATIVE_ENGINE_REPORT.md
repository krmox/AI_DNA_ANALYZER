# AI_DNA_ANALYZER V2 — Native Engine: engineering report

Date: 2026-09-24. Scope: **same science, different execution engine.** Frozen V1 thresholds, Binomial /
router / Poisson-binomial / Method C semantics, MAX_READS = 48, SNP-only scope, reference, truth, BED and the
benchmark region are unchanged. Nothing in V1, the manuscript, Zenodo or the release was modified.

Companion files: `V2_ROUTING_ANALYSIS.md` (routing / rescue / loss analysis, written earlier),
`V2_CORRECTNESS_RESULTS.csv` (every correctness check, 139 rows), `V2_PERFORMANCE_RESULTS.csv` (every timed run).

## 0. Status

| Item | Result |
|---|---|
| Math equivalence (kernels) | **PASS** — bit-identical to frozen V1 NumPy/SciPy |
| Real-locus equivalence (3,491-locus gate) | **PASS** — 17/17 field checks, 0 mismatches, incl. 1,340,544 V1-tensor bytes |
| Full chr20 correctness | **PASS** — Binomial and PB LLR bit-identical on all 56,266,816 loci; cascade and PB-only VCFs byte-identical to V1 |
| Runtime, memory measured | **DONE** (Section 9) |
| Routing / rescue / losses analysed | **DONE** (`V2_ROUTING_ANALYSIS.md`) |
| hap.py 0.3.15 repeated | **BLOCKED** — Docker image store is an unmounted loopback (Section 10). Not run. |

Per the brief's final rule, V2 is **not yet fully "successful"**: criterion 8 (hap.py) is open. Byte-identical VCFs
make identical hap.py output a logical consequence (hap.py is deterministic in its inputs), but that is an
inference, not the requested re-run.

## 1. Architecture

```
BAM ──> ONE streaming htslib pileup per task (contiguous run of ≤2048 windows)
          per column:  count_column ──> Binomial LLR ──> router ──┬─ not routed ─> Binomial call
                                                                  └─ routed ─────> 48 retained reads (V1 order)
                                                                                    └> PB LLR (scalar, exact) ─> call
        ──> Method C genotype ──> VCF text  (streamed, in genomic order)
```

* No read tensor exists at any point. Per locus, PB receives ≤48 `(base, Phred)` pairs straight from the pileup.
* Python appears only in tests/oracle/analysis. The engine is a C++20 executable (`dnav2`), no Python boundary.
* Workers own their BAM handle, FASTA handle and buffers; tasks are claimed from an atomic counter; results are
  flushed strictly in genomic order, so output is deterministic for any thread count (verified: 1 vs 4 threads
  give byte-identical dumps and VCFs).
* Two extraction paths: **fused** (default, one pileup per task) and **`--exact-windows`** (V1-literal: separate
  counts pass + a *fresh* pileup per 64-bp window). Equivalence of the two was measured, not assumed (Section 7).
* Modes: `cascade` (PB only where the router sends it — 707 loci) and `pb_all` (PB on every locus; also writes the
  PB-only VCF; needed for the control arm and for the strongest equivalence test).

## 2. Changed / new files (all under `native_v2/` in worktree branch `native_v2`; nothing tracked was edited)

* `include/dnav2/numerics.hpp` — cephes `gammaln` and `log1p` ports (SciPy's, bit-exact), `logaddexp`, `xlogy`, `xlog1py`
* `include/dnav2/binomial.hpp` — Binomial score, frozen constants, router predicate
* `include/dnav2/pb.hpp` — `pb_reference` (transparent) and `pb_native_v2`
* `include/dnav2/methodc.hpp` — Method C genotype
* `include/dnav2/blake2b.hpp` — read-name ordering key (V1 algorithm)
* `include/dnav2/extract.hpp`, `src/extract.cpp` — counts/rows/selection helpers, fused pass, exact-window pass
* `include/dnav2/frame.hpp`, `src/frame.cpp` — window list (BED/N filter), truth-derived frame labels, FASTA access
* `include/dnav2/engine.hpp`, `src/engine.cpp`, `src/main.cpp` — task scheduler, workers, ordered output, CLI
* `src/capi.cpp` — C ABI used only by the tests
* `oracle/v1_oracle.py` — frozen V1 semantics as functions (V1 modules unmodified; Method C/VCF copied verbatim)
* `tests/*.py`, `scripts/*`, `tools/*`, `Makefile`
* `research/V2_*` — this report and the two CSVs

Not committed (standing rule: no commits unless asked). Backup: tag `v1-frozen-pre-v2-20260924` +
`/mnt/archive/AI_DNA_ANALYZER_v2/backup/` (uncommitted patch + tar of untracked files).
Everything V2 lives on the HDD (`/mnt/archive/AI_DNA_ANALYZER_v2/`), as requested.

## 3. Profiling evidence (no perf/valgrind on this machine; sudo unavailable)

Tools used: stage timers inside the engine, a profiling build (`-DDNAV2_PROFILE`), `/usr/bin/time -v`,
`bamscan` (BGZF+decode floor), `itrprobe` (iterator cost), GCC `-fopt-info-vec`. No hardware counters or memory-bandwidth
figures were available; **no claim about memory-bound vs compute-bound behaviour at the hardware-counter level is made.**

Slice benchmark (chr20:31,818,624–32,398,272, 8,192 contiguous windows, 524,288 loci, 144.2 M read-locus pairs, 1 thread, warm):

| Stage (fused path) | cascade | pb_all |
|---|---|---|
| htslib pileup incl. BAM decode | 3.29 s (52 %) | 3.38 s (33 %) |
| `count_column` | 1.34 s | 1.34 s |
| Binomial + router | 0.18 s | 0.21 s |
| BLAKE2b per read | 0.33 s → removed in cascade by lazy hashing | 0.33 s |
| row build + 48-read selection | ≈0 | ≈3.5 s |
| **PB arithmetic** | ≈0 | **0.33 s (3 %)** |
| wall | 6.29 s → 5.84 s after lazy hashing | 10.3 s |

Other measurements: BAM decode floor for this region 0.86 s warm / 1.82 s cold (0.28 s with 4 decompression threads);
a fresh 64-bp iterator query costs ≈2.1 ms and returns 362 reads (BAI linear-index granularity), which is why V1's
per-window passes are expensive (47.9 s vs 11.4 s fused for the same output).

**Consequence:** the Poisson-binomial *kernel* is not the bottleneck in V2. The time is spent in htslib's pileup, per-read
counting and (in `pb_all`) read selection.

## 4. PB implementation

`pb_reference`: 1:1 port of V1's loops (full 49-wide log-space DP over all 48 slots). `pb_native_v2`: identical arithmetic with
exact work reduction — (a) unused slots skipped, (b) DP truncated to columns 0..k (column *j* depends only on *j* and *j−1*, only column *k* is read),
(c) column 0 is a running sum (`logaddexp(x, −∞) ≡ x`), (d) `log q`, `log1p(−q)` from a 256-entry Phred table holding the very same doubles.
Contract: `poisson_binomial_llr(..., alt_index=…)` — **not** the `min(k, n_counted)` variant (2,400 of 3,491 loci differ; verified again here).

Kernel tests (`tests/test_kernels.py`, all bit-level): `gammaln` (210,001 arguments), cephes `log1p` (400,003), Binomial incl. router mask
at depth ≤ 5/60/400/1,500/9,000 (196,000 loci), PB reference and native on a synthetic grid of 8/16/24/32/48/random reads × zero/low/het/hom/any VAF ×
{Q25–39, Q0–93, extreme {0,1,2,12,13,60,93,255}, Q13–21} (18,000 loci, all bit-identical, decisions identical), Method C (100,000), BLAKE2b (20,008).

Two numerical facts found while testing (both would have silently broken exactness): SciPy's `binom.logpmf` uses **cephes `log1p`**, which differs from glibc `log1p`
for `p_het` in ~94 % of counts; NumPy's PB path uses glibc. Both are reproduced. Build must use `-ffp-contract=off`.

## 5. SIMD — measured, null result

* GCC auto-vectorization report (`scripts/vec_report.sh`, `logs/vec_report_*`): the hot loops (`count_column`, read selection, PB DP) are **not** auto-vectorized in either generic or AVX2 builds; only initialisation loops and small SLP fragments are.
* AVX2+FMA build (`-mavx2 -mfma`, FMA contraction still off): outputs byte-identical to the generic build, but **not faster** — cascade 6.06–6.27 s vs 5.83–5.99 s, pb_all 10.6–13.4 s vs 10.2–10.7 s (3 reps each).
* Hand-written intrinsics were **not** written: PB arithmetic is ≤3 % of `pb_all` and ≈0 % of `cascade`, so even an infinite PB speed-up is bounded by ≈3 %; and a SIMD `exp/log1p` cannot be bit-identical to glibc, so it would need a numerical-tolerance policy (`MODEL_CHANGE`-adjacent) for no measurable gain.
* No assembly. **The V2 engine currently ships without SIMD.**

## 6. Parallelization and memory architecture

Threads (not processes): workers share nothing mutable except the atomic task counter and the ordered-output mutex. Task = ≤2048 contiguous windows
(131 kb); per-task state ≈ a few MB (`CountRow` 36 B/locus + per-locus scalars). Output is streamed, so peak RSS is bounded by (workers × task) + the
truth-label array (64 MB, one byte per contig position) + window list (7 MB).

Scaling, full chr20 (1 run each; slice-level repeat spread ≈ ±2–5 %):

| mode | threads | cache | wall | speed-up vs 1 thread | CPU util | peak RSS |
|---|---|---|---|---|---|---|
| cascade | 1 | cold | 700.6 s | — | 78 % (I/O wait) | 110 MB |
| cascade | 2 | warm | 310.7 s | 2.26× *(vs cold t=1; see note)* | 198 % | 127 MB |
| cascade | 4 | warm | 165.0 s | 4.25× *(same note)* | 395 % | 161 MB |
| cascade | 8 | warm | 119.6 s | 5.86× *(same note)* | 669 % | 225 MB |
| cascade, 4 workers + 2 BGZF threads each | 4 | warm | 176.6 s | slower than plain 4 | 452 % | 133 MB |
| pb_all | 1 | warm | 1059.5 s | — | 99 % | 114 MB |
| pb_all | 4 | warm | 309.0 s | 3.43× | 392 % | 162 MB |
| pb_all | 8 | warm | 182.5 s | 5.81× | 721 % | 225 MB |

Note: the cascade 1-thread run was cold (BAM evicted from page cache) and includes HDD wait, so cascade speed-ups against it are overstated; the pb_all
column (1-thread run warm) is the clean scaling measurement: 3.43× on 4 workers (4 physical cores), 5.81× on 8 logical threads.
Hyper-threads add ≈1.4–1.7× over 4 workers here (measured, not assumed). Extra BGZF decompression threads made cascade slower. Cold 4-thread cascade: 281.4 s; cold 8-thread pb_all: 497.5 s at only 256 % CPU — the cold
runs are I/O-limited on the 5400-rpm HDD (this is a property of the test disk, not of the engine).

## 7. Correctness validation

Gate populations (frozen in `research/benchmark_chr20_300x/routing_analysis/loci_of_interest.tsv`): all 707 routed, all 892 rescuable (incl. the 191 rescues),
226 frame true-SNP losses, 1,906 PB≠Binomial and 1,593 cascade≠PB-only loci, 1,000 seeded-random loci = 3,491 loci.

* **Frame**: window list (879,169 windows) and truth-derived labels equal the V1 cache on all 56,266,816 loci; frame count 56,242,697 (`tests/test_frame.py`).
* **Oracle validity**: the reconstructed V1 adapter reproduces the archived V1 VCFs byte-for-byte from V1's own caches (72,872 / 72,661 records) — the original 300x adapter script no longer exists; its text was preserved and re-validated.
* **3,491-locus gate**, per locus, exact: counts cols 0–8, reference index, label, frame flag, Binomial LLR (bit), candidate ALT, *k*, *n*, router mask, Binomial call, PB LLR (bit, max |Δ| = 0), PB call, cascade call, V1 tensor rows (48×8 bytes, 1,340,544 bytes), n_counted, *k* in tensor; 37 loci with no admitted read hold exactly V1's padding row. 1-thread, 4-thread, fused and exact-window dumps are byte-identical to each other.
* **Fused vs V1's per-window pileup**: on 524,288 contiguous loci both paths are bit-identical to the V1 cache (Binomial + PB), and byte-identical VCFs; on the **whole chromosome** the fused path is bit-identical to the cached V1 PB on 56,266,816/56,266,816 loci (0 mismatches). This is an empirical result for this BAM, not a proof in general: V1's own source warns that per-window and long passes can differ when read names repeat. The `--exact-windows` path is kept for that reason.
* **Full chr20 `pb_all`** (8 threads): Binomial LLR, PB LLR bit-exact on all loci; router mask, Binomial/PB/cascade calls exact on 56,242,697 frame loci; exactly 707 routed loci (same positions as V1); cascade VCF (3,471,975 B) and PB-only VCF (3,481,518 B) **byte-identical** to the archived V1 VCFs. Every timed full run re-checked the cascade VCF against V1 (all identical).
* **Real V1 vs V2, same conditions (slice)**: the real V1 pipeline (native backends asserted active) and V2 produce byte-identical cascade and PB-only VCFs.
* No divergence was found, so nothing had to be explained away. During development three of my *own* mistakes were caught by these tests and fixed before any result was recorded (SciPy `log1p` ≠ glibc; a comparator that mis-handled the 37 empty loci; a test-script path that silently ran V1 in pure Python). Those intermediate failures are not in the CSV (it records runs from the final code).

## 8. Routing, rescue, losses, "564"

See `V2_ROUTING_ANALYSIS.md`. Unchanged facts used here: 707 routed = 191 rescue (184 FP-avoided, 7 FN-avoided) + 122 harmful (all FP-created) + 277 both-correct + 117 both-wrong; net +69 decisions. **"564 unnecessary loci" does not exist** (frozen: rate 0.564, 399 loci; non-rescue routes 516). The V2 routed table (`routed_loci_chr20_300x.tsv`) reproduces the same 707 positions and calls.

## 9. V1 vs V2 — runtime and memory

Three different quantities, kept separate:

**(a) Algorithmic work reduction** (compute PB only where the cascade uses it; no tensor). Same engine, 1 thread, warm:
`pb_all` 1059.5 s → `cascade` ≈ 540 CPU-s (user time; the only 1-thread cascade run was cold, wall 700.6 s). ≈ 1.5–1.9× fewer seconds; **this is not an implementation speed-up and must not be described as a faster PB.** The PB-only control arm still needs PB on all loci.

**(b) Implementation speed-up, identical work** (V1 computes counts + tensor + PB everywhere; so does V2 `pb_all`):
* Slice, same machine, same day, real V1 vs V2 1 thread: **V1 158.2 s → V2 pb_all 10.6 s ≈ 14.9×** (V2 wall includes ≈1.1 s of one-time window/label setup that is negligible on a full chromosome). V1 stages: counts 6.3, read-tensor 59.8, PB 91.8, rest 0.3 s.
* Full chromosome, V1 from the frozen benchmark log (not re-run; single-threaded stages summing to 18,578 s ≈ 5 h 10 min) vs V2 `pb_all`: 1 thread 1,059.5 s (**≈17.5×**), 4 threads 309.0 s (**≈60×**), 8 threads 182.5 s (**≈102×**). The V1 numbers come from another day/cache state; treat ratios as indicative to ±15 %.

**(c) End-to-end, the deliverable V1 actually produced (both VCFs) vs the V2 cascade product** — combines (a), (b) and parallelism:

| | wall | note |
|---|---|---|
| V1 (frozen record) | 18,578 s ≈ 5 h 10 min | 1 thread |
| V2 cascade, 4 threads, cold | 281.4 s | |
| V2 cascade, 4 threads, warm | 165.0 s | ≈113× |
| V2 cascade, 8 threads, warm | 119.6 s | ≈155× |
| V2 pb_all (all arms, both VCFs), 8 threads, warm | 182.5 s | ≈102× |

Stage table (V1 frozen vs V2): extraction 790.5 + 6,553.8 s → part of the single fused pass; PB 11,137.8 s → 0 s (cascade: 707 loci ≈ ms; pb_all: inside 182–1,059 s); router+genotype+output 96 s → 0.4 s of CPU.
V2 does not report separate extraction/read-preparation/PB times for the fused pass because they are interleaved per column; the profile in Section 3 gives the CPU split.

**Memory**: V1 read tensor 21.6 GB raw; V1 PB streaming ≈ 8.0–8.2 GB RSS (frozen record); V1 on the slice peaked at 1,291 MB for 0.93 % of the chromosome. V2 peak RSS on the full chromosome: 110 MB (1 thread) – 225 MB (8 threads), cascade and pb_all alike. Ratio vs V1 PB-stage RSS ≈ 36–75×. CPU utilisation, wall, user/sys time: table in Section 6 and the CSV. Memory bandwidth: not measured (no counters).

## 10. hap.py — not repeated

The Docker `data-root` is `/mnt/archive/AI_DNA_ANALYZER_benchmark/container_storage_ext4/docker`, the mount point of the ext4 loopback image
`container_storage.img`. After the reboot/remount the loopback is not mounted, so `docker images` is empty and `hap.py:0.3.15` is unavailable.
Re-mounting needs `sudo` and a Docker restart — a system-level change I did not make without your permission (and cannot: no passwordless sudo).
Once you approve, the V1 protocol can be replayed unchanged (xcmp engine, `-l chr20`, truth NISTv4.2.1, BED `chr20_highconf.bed`, 8 threads, same reference) on the already-produced V2 VCFs in `/mnt/archive/AI_DNA_ANALYZER_v2/work/full/pb_all_t8/`.
Until then: V1 cascade F1 0.972429 / PB-only 0.971104 stand as recorded, and V2's VCFs are byte-identical to the ones those numbers came from.

## 11. Remaining bottlenecks

1. htslib pileup + BAM decode: 52 % of single-thread `cascade`. Replacing it means re-implementing htslib's mate-overlap quality adjustment exactly — high exactness risk; the full-chromosome comparison against V1 would be the safety net. Not attempted.
2. `count_column`: ≈ 22 % (≈ 9 ns per read-locus, pointer chasing through pileup nodes).
3. `pb_all` only: read selection (≈ 3.5 s of 10.3 s per slice) — an incremental (sliding) selection would cut it; not done because `pb_all` is a control mode.
4. Cold-cache HDD I/O (CPU 78 % at 1 thread; 256 % at 8 threads).
5. One-time setup ≈ 1.1 s (window list + labels), serial.

## 12. Known limitations

* The engine is not truth-free: it reads the truth VCF **only** to reproduce V1's scored-locus frame (labels 0/1); it never affects a score, threshold or call. A production caller without truth would drop that filter — a different output (labels 2/3 loci would be emitted) → `MODEL_CHANGE_REQUIRED`-class decision, not made.
* M-1 (64-bp windows fully inside BED) is reproduced deliberately; 1,041 of 1,321 hap.py FN remain unreachable by design.
* MAX_READS = 48, SNP-only, fixed ε: unchanged. The observation that PB's outcome depends on which 48 reads survive at 300× is descriptive; testing it is a separate `MODEL_CHANGE` experiment.
* Exactness relies on glibc's libm matching NumPy's on this machine (verified by the tests) and on `-ffp-contract=off`; another libm/CPU may need the tests re-run.
* Like V1's native path, the engine refuses inputs V1 could not reproduce (`l_qseq == 0`, non-integer `NM`) instead of guessing.
* Fused-vs-per-window equivalence was verified on this BAM only.
* Built against pysam-bundled htslib headers and the system `libhts.so.3` (1.23.1), same as V1's native code.
* Timing caveats: single runs at full scale; V1 full-chromosome baseline is the frozen record, not a re-run; page-cache state differs between rows (column `cache`).
* Code is uncommitted in the `native_v2` worktree (tag + backup exist).
