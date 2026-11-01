# Release readiness - AI_DNA_ANALYZER / Varion (R2 SNP caller, nindels, R3.1)

Date: 2026-11-01. Evidence: `BUILD_AND_TEST_REPORT.md`, `PLATFORM_SUPPORT.md`, `PORTABILITY_AUDIT.md`.
Written against the working tree on top of `main` @ 94ae97d, before the portability commit that adds this file; the `r3.1` worktree statements refer to its own uncommitted files.

## 1. Exact test results (CachyOS x86-64, GCC 16.2.1, Python 3.14.8)

`bash verify.sh full` with htslib **1.23.1**, working tree and extracted source archive (same result both times):

| Group | Passed | Failed | Skipped |
|---|---|---|---|
| V1 (`115 passed, 2 skipped, 21 xfailed, 5 xpassed`) | 120 | 0 | 2 (documented) |
| V2/V2.x kernels bit-identical to V1 | 17 | 0 | 0 |
| R2 (`40186 checks, 0 failures`; frame suite passed) | 40,187 | 0 | 0 |
| nindels | 153 | 0 | 0 |
| **Total** | **40,477** | **0** | **2** |

Not run: `verify.sh --data` (samtools/bcftools missing). R3.1 branch (separate worktree, commit 31e14df): `nindels-tests` 156 passed, `tests_varion` 22 passed.
Old-vs-new build systems give byte-identical R2 and V2.x VCFs and routed-loci files on the shipped chr21 slice, under htslib 1.21 **and** 1.23.1 (and the
hashes are the same across the two htslib versions).

## 2. htslib actually used
**1.23.1**, built from the upstream tarball by `scripts/build_local_htslib.sh` (SHA-256 verified) into `third_party/htslib`. Headers and library come from that
prefix (`verify.sh` report; `hts_version()` = 1.23.1); all five native artefacts link `third_party/htslib/lib/libhts.so.3`; no reference to the 1.21 prefix; no
`LD_LIBRARY_PATH`. pysam 0.24.0 carries its own private htslib 1.23.1 for Python only.

## 3. Conda verification: **NOT verified**
No conda/mamba on the machine; downloading and running micromamba was denied by the sandbox (twice; not worked around). Unverified: the dependency solve for
Python 3.14, the installed htslib/samtools/bcftools versions, native compilation inside the env, `verify.sh full` inside the env. Only channel metadata was
checked: htslib/samtools/bcftools 1.23.1 exist for linux-64, linux-aarch64, osx-64, osx-arm64. A researcher who runs `conda env create -f environment.yml`
is the first tester of that route.

## 4. Source archive: **verified (Linux, one machine)**
`scripts/make_source_archive.sh --worktree` -> 429 files, 2.1 MB; contains sources, headers, tests, models, build scripts, README and docs; no research data
(largest file 1.3 MB). Extracted into an empty directory, built and fully verified (table above). The default mode archives committed `HEAD` and therefore does
**not** contain the portability changes until they are committed. The 857 MB GitHub-generated archive still contains all evidence by design.

## 5. Remaining limitations
- Tested on one Linux distribution only; Fedora 43, Ubuntu, macOS, WSL2 are documented-but-untested; native Windows unsupported.
- Docker image and conda environment never built; `--data` path never run after the changes.
- `research/stage4/R2_LOCK.md` per-file hash for `native_r2/Makefile` no longer matches (outputs verified identical); `verify_identity.sh` was already stale.
- Campaign/benchmark scripts keep the author's `/mnt/archive`, `/mnt/hdd` paths (by design, provenance); full-benchmark reproduction still needs editing them.
- htslib other than 1.23.x only warns; Ubuntu/Debian distro packages will take that path.
- Clone is ~0.5 GB (history) / 1.3 GB checked out; unchanged.
- Smoke-level equivalence check only; no published benchmark was re-run in this audit.
- Merging `r3.1` into `main` (dry run against committed `main`): conflict in `.gitignore`; the branch predates the portability changes, so README/`verify.sh` will need reconciling.

## 6. R3.1 (branch `r3.1`): final test completed 2026-11-01 - is it ready for public release? **Not as an improvement; only as an experimental feature.**
Full report: `docs/R3.1_FINAL_TEST_REPORT.md` in the `r3.1` worktree (raw tables `results/r31/test_report.*`, uncommitted there).
- **Pre-registered final test:** all required clean cells are scored - HG001 chr17, HG007 chr17, HG007 chr19 x depth ladder x {baseline, A, A+B, adaptive} = 68/68. The one missing cell (HG007 chr19 10x adaptive) was run; it chose A and equals A.
  Adaptive R3.1 equals A in all 12 cells below the frozen 67.37x switch and equals A+B in all 5 cells above it (record-identical calls in all 17).
- **Outcome:** A+B vs A at >= 67.37x was reliably better in 1/5 cells (+0.0009), indistinguishable in 3/5, reliably worse in 1/5 (HG001 chr17 native, -0.0021). No depth range shows a consistent advantage for A+B.
  The protocol defines **no numeric acceptance criterion**, none was added, so there is no pass/fail verdict - and **no evidence that R3.1 is better than configuration A**. Its benefit is avoiding A+B's loss at low depth.
- **Documented gaps:** repeat-context breakdown unavailable for `HG007_chr19 d10 R31` (two stratification BEDs were on the unmounted data disk and are not reproducible from GIAB; overall metrics verified unaffected);
  local-depth breakdown unavailable for the two HG007 native cells (BAMs on the unmounted disk); R3.1 not run on "previously seen" HG001 chr19 (truth/BED on the unmounted disk); optional HG003 chr8 not run.
  Three experiment-script edits needed for this were reverted and saved as `results/r31/session_2026-11-01_script_changes.patch`.
- Regression tests: `nindels-tests` 156 passed / 0 failed, `tests_varion` 22 passed. Tools in the experiment: htslib/samtools 1.21, hap.py 0.3.15 (not 1.23.1).
- Still true: the branch is unmerged and changes `nindels/src` (algorithm code), which this audit did not review; `docs/R3.1_THRESHOLD_SELECTION.md`, cited in `varion/thresholds.json`, is missing from the repository.

What **is** releasable after you commit and tag it yourself: the portability/build infrastructure and the documented R2 + nindels (Baseline/A/B/A+B) state on `main`,
described as "research software, Linux x86-64 verified, other platforms untested".

## 7. Release-blocking defects
1. **R3.1 branch unmerged, its new report/results uncommitted, `R3.1_THRESHOLD_SELECTION.md` missing.** The earlier "final test incomplete" blocker is resolved (section 6), but the result is a null/negative one:
   any wording that R3.1 or A+B *improves* accuracy over A is unsupported and must not be used.
2. **Portability work is committed together with this document** (14 modified, 11 new files). Before that commit nothing a researcher clones contained it, and `scripts/make_source_archive.sh` in committed mode lacked it (use `--worktree` to archive uncommitted state).
3. Not blocking a Linux source release, but blocking any "cross-platform" or "conda/Docker supported" wording: conda and Docker untested.
No defect was found in the verified scientific path (no failing test, no output difference).
