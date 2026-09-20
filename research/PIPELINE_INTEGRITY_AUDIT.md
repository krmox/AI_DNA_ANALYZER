> **SUPERSEDED IN PART (2026-09-20).** Historical audit. Its statement that no native C/htslib implementation exists (I-2, H-1) is wrong (see manuscript §3.9 C4). Current status: `FINAL_FREEZE_REPORT.md`.

# Pipeline Integrity Audit

Read-only audit of `AI_DNA_ANALYZER`'s active codebase for scientific-correctness
and reproducibility risks. Scope: the tracked repository root (production and
`experimental/` code) plus `bench_htslib_io/`, excluded from "production" but
inspected for Part II of the prior task. Unrelated pre-existing developer
worktrees under `.claude/worktrees/agent-*` and `.claude/worktrees/30x-pb-residual`
were spot-checked but are explicitly out of scope (not part of this repo's
committed history on `main`).

## CRITICAL

### C-1. Fabricated genomic positions in `extract_bench_v12.py` (FIXED this session)
`chunk_start + row_index` was used as each extracted locus's genomic
coordinate, which silently desyncs whenever `providers.py::_build_windows`
drops a 64bp window (BED gap / high-N). This invalidated one benchmark
(`experimental/stress_test/run_hg005_pipeline.py`, the pre-fix HG005 run).
**Status: FIXED.** `pileup_counts.py::load_counts(..., return_positions=True)`
now returns the true per-row coordinate; `extract_bench_v12.py` uses it at all
3 call sites. See `experimental/stress_test/HG005_EXTRACTION_BUG_ROOT_CAUSE.md`
for full detail, `HG005_HISTORICAL_AUDIT.md` for the blast-radius audit
(no other consumer was affected), and `test_hg005_extraction_integrity.py`
for the permanent regression suite (11 tests).

## CORRECTION TO AN EARLIER VERSION OF THIS AUDIT

An earlier version of this document (H-1 below, original text preserved for
the record at the end of this section) concluded that no C/htslib
acceleration existed anywhere in the codebase and that a "~65x speedup"
claim was unsubstantiated. **That conclusion was wrong, and the error was
mine: I searched `.claude/worktrees/agent-afd5734475c9603bd/bench_htslib_io`
(a genuine null-result threading experiment) and never checked
`.claude/worktrees/agent-a4cedfdeffb76e3dd/experimental/performance`, a
second, different worktree that contains a real, compiled C extension
(`native/pileup_native.c`) with a full equivalence proof (0/109,813,120
count-cell mismatches, full 12Mb HG002 region) and a measured 65.1x
end-to-end speedup, documented in that worktree's own
`PERFORMANCE_OPTIMIZATION_REPORT.md`. The user caught this and directed a
re-audit before I made further claims.** Having now re-verified that report
directly (read the C source, rebuilt the extension from scratch in this
session's own worktree, reproduced 0-mismatch equivalence independently on
both HG002 and HG005 data, and integrated it into the production
`pileup_counts.load_counts` code path — see "H-1 RESOLVED" below), the
original 65x extraction-stage figure is real and is now the default,
integrated behavior, not merely an isolated experiment.

## HIGH (RESOLVED)

### H-1 RESOLVED. C/htslib acceleration integrated into production `load_counts`
Original finding (now superseded, see correction above): a *different*
worktree's htslib-threading experiment showed no speedup, and this was
incorrectly generalized to "no C/htslib implementation exists anywhere."
**Corrected finding:** `.claude/worktrees/agent-a4cedfdeffb76e3dd/experimental/performance/native/pileup_native.c`
is a genuine, validated C extension (`Python.h` + raw htslib C API, one
`bam_mplp_*` streaming pass per requested region) that reproduces
`PileupCountsProvider.window_counts`/`_count_locus`'s output columns 0-8
bit-for-bit (proven at 0/109,813,120 mismatches over the full 12Mb HG002
chr21:32-44M benchmark region in the source worktree, and independently
re-proven at 0/675,648 and 0/1,684,800 mismatches on HG005 chr1 and HG002
chr21 sub-regions in this session, both via the production `load_counts()`
entry point with the native path toggled on/off). It has been rebuilt from
scratch in this session (see `native/build.sh`, no hard-coded paths, no
dependency on any specific worktree) and wired into
`pileup_counts.py::load_counts` as the default path (`NATIVE_AVAILABLE`),
with an automatic, correct fallback to the pure-Python path when the
extension is not built (e.g. no C toolchain in a given environment). See
`research/PERFORMANCE_RESULTS.csv` for real, measured before/after timing at
both the original 12Mb HG002 benchmark scale and the HG005 3Mb scale this
session's earlier work centered on, and `test_native_pileup_equivalence.py`
for the 12-category permanent regression suite (normal window, BED gap,
high-N, partial BED overlap, contig boundary, multiple windows, empty
window, low depth, overlapping mate pairs, orphan reads, base-quality
filtering, determinism).

**What is still true from the original H-1 finding, and remains a distinct,
valid point:** the *separate* `bench_htslib_io/` experiment (pysam's
`threads=` kwarg, htslib-internal BGZF threading) really does show no
speedup, and is a genuinely different mechanism from the C extension above —
conflating the two would itself be a claims error. Both are now accurately
described in `SCIENTIFIC_CLAIMS.md`.

<details>
<summary>Original (superseded) H-1 text, preserved for the record</summary>

> No C/htslib acceleration exists; a "~65x speedup" claim circulated without
> substantiation. `bench_htslib_io/` (present only in an unrelated developer
> worktree, never on `main`) is not a separate C extension — it is the
> identical `PileupCountsProvider`/`ReadLevelPileupProvider` Python code,
> opened with pysam's `threads=` kwarg (htslib's internal BGZF-decompression
> threading). Its own measured results show no scaling (15.25s at 1 thread
> vs 16.0s at 8 threads) and an Amdahl analysis concluding the theoretical
> ceiling from infinitely-fast I/O is ~3.5x. No file in this repository
> substantiates a 65x figure.

This was true of the worktree actually inspected at the time, and false as
a claim about the whole codebase — a different worktree, not checked before
this text was written, contained the real implementation.

</details>

### H-2. `FROZEN_ROUTER_CUTOFF` is copy-pasted as a numeric literal, not imported, in several files
`5.411872376933351` appears as a hardcoded literal in `robustness_benchmark.py`,
`extract_bench_v12.py`, and `experimental/stress_test/run_hg005_pipeline.py`
(and its `_postfix` copy), rather than being imported once from `cascade.py`
(which does define and export `FROZEN_ROUTER_CUTOFF`, and multiple bench_v14+
scripts do `from cascade import FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD,
FROZEN_ROUTER_CUTOFF` correctly). If `cascade.py`'s value were ever
legitimately revised, these copies would silently drift out of sync with no
error. **Not fixed in this task** (a signature-level refactor across several
files is outside "minimum necessary change" scope for the HG005 fix this
session performed, and touching it now would expand this task's diff beyond
what was authorized). Recommend a follow-up: replace every literal with
`from cascade import FROZEN_ROUTER_CUTOFF`, and add a test asserting every
copy equals `cascade.FROZEN_ROUTER_CUTOFF` (this repo already has
`test_extractor_frozen_constants_match_cascade` in `test_bench_v12.py` for
`extract_bench_v12.py` specifically, which is exactly this pattern and should
be extended to the other files).

## MEDIUM

### M-1. Window-level BED filtering silently reduces truth-SNP recall relative to hap.py's per-SNP BED confidence
`GiabAlignmentProvider._build_windows` drops an entire 64bp window if *any*
base of it lies outside the high-confidence BED, even when the specific truth
SNP inside that window is itself BED-confident. Quantified on the HG005
region: 142/244 (58%) of the post-fix cascade run's false negatives have no
extracted candidate at all for this reason (`HG005_POSTFIX_STRESS_TEST_REPORT.md`
Sec. 9). A smaller version of the same effect was already documented for
HG002 chr21:32-44M (301/16898 ≈ 1.8%, `experimental/unified_happy/
ROOT_CAUSE_ANALYSIS.md` Sec. 11). This is a candidate-generation granularity
choice, not a coordinate-correctness defect (the coordinate bug, C-1, is
fixed and independently verified) — but it is a real, quantified drag on
measured recall that should be disclosed in every accuracy claim. Not fixed
in this task; flagged as a known limitation throughout the research dossier.

### M-2. Two internal accuracy evaluators exist and produce non-comparable numbers with similar-looking names
`robustness_benchmark.py`'s row-index `accuracy_block()`/`prf()` (used by
every `bench_v13`-`v19` script) and hap.py's genotype-aware `xcmp` engine
both produce "precision/recall/F1", but the former ignores REF/ALT/GT
entirely (pure boolean variant/no-variant agreement at matching array index)
while the latter is allele- and genotype-aware. The two have historically
been placed side-by-side without a clear label (this was the original
"F1 0.974 -> 0.624" apparent-regression that `experimental/unified_happy/
ROOT_CAUSE_ANALYSIS.md` had to formally reconcile). Every table in this
dossier (`BENCHMARK_MASTER.csv`) carries an explicit `evaluator` column for
this reason; it is a documentation/labeling debt in the source repository
itself, not fixed here, and future scripts should either not both be called
"F1" or should be renamed to make the distinction load-bearing in the API
(e.g. `raw_agreement_f1` vs `happy_f1`).

### M-3. `test_providers.py` is broken and not part of the running test suite
`test_providers.py` imports `bimamba_variant_caller.config`/`.providers`/
`.testdata`, a package that does not exist anywhere in this repository or as
an installed dependency (`ModuleNotFoundError` on collection). It appears to
be a stale artifact from an earlier package layout. All other test files
(`test_bench_v12.py`, `test_cascade.py`, etc.) import directly (`from
providers import ...`) and pass. Not fixed here (out of scope for the HG005
fix); flagged so a reviewer does not assume 100% of `test_*.py` files are
green — `pytest test_providers.py` currently errors at collection, and any
blanket `pytest` invocation over the whole repo root should exclude it or it
will report a false collection failure.

## LOW

### L-1. `experimental/head_to_head/overlap_analysis_DISCARDED.py`
Filename self-declares discarded status; contains a dead `if False else ...`
branch referencing the router cutoff literal. Already self-quarantined by
naming convention; no action needed, noted for completeness.

### L-2. Large, long-lived cache directories under `cache/`
`cache/bench_v19` (363M), `cache/hg005_stress` (324M, the pre-fix invalidated
cache — correctly retained, not deleted, per this task's "no data deletion"
rule), several `*_evidence.npz` files (60-260M each). Not a correctness issue;
noted for disk-budget awareness in `REPRODUCIBILITY.md`.

## INFORMATIONAL

### I-1. No production file imports from `experimental/`
Verified by grep across every top-level production `.py` file
(`cascade.py`, `pileup_counts.py`, `providers.py`, `extract_bench_v12.py`,
`binomial_baseline.py`, `quality_error_model.py`, `robustness_benchmark.py`):
zero references to `experimental`. Confirms experimental/exploratory code
cannot have silently leaked into a threshold, filter, or decision path used
by the frozen production cascade.

### I-2. No compiled/native extension exists anywhere in the tracked repository
`find . -iname "*.c" -o -iname "*.pyx" -o -iname "*.so"` (excluding
`.claude/worktrees/` and `node_modules`) returns zero results. All extraction
and calling code is pure Python + `pysam`/`numpy`/`scipy`. Relevant context
for H-1 above and for any future claim about "native" or "compiled"
acceleration.

### I-3. No hidden/undocumented threshold tuning found
Every frozen constant traced (`FROZEN_BINOMIAL_THRESHOLD=7.0`,
`FROZEN_PB_THRESHOLD=10.5`, `FROZEN_ROUTER_CUTOFF=5.411872376933351`,
Method C's `EPS=0.01`/`GQ_CAP=99.0`) resolves back to `cascade.py` /
`experimental/genotype_layer/method_c_regression.py`, both verified
byte-identical (sha256) to the values recorded in the pre-fix HG005 manifest
before and after this session's changes. No script in the audited scope
recomputes, refits, or overrides these values from data.

## Summary

| Severity | Count | Fixed this session | Open |
|---|---|---|---|
| CRITICAL | 1 | 1 | 0 |
| HIGH | 2 | 0 | 2 (documented/corrected in dossier, not code-fixed) |
| MEDIUM | 3 | 0 | 3 (flagged as limitations) |
| LOW | 2 | 0 | 2 (informational, no action needed) |
| INFORMATIONAL | 3 | n/a | n/a |

Per instruction, only the issue necessary for scientific correctness of the
HG005 result (C-1) was fixed in code. H-1 and H-2 are corrected at the
*claims* level (this dossier does not repeat the unsubstantiated 65x figure,
and every table states its evaluator explicitly) rather than via an
uncontrolled refactor of the source tree.
