# Build and test report

Machine (the only one available): CachyOS (Arch), x86-64, Linux 7.3, GCC 16.2.1, CMake, Python 3.14.8.
No conda/mamba/micromamba, no Docker/Podman, no samtools/bcftools, no system htslib.
Dates: 2026-11-01. Statuses: **Tested** = executed here and passed; **Blocked** = attempted, prevented; **Not run**.

## Run 0 - baseline before any change (first audit, unchanged scratch copy)

| Step | Result |
|---|---|
| `pip install -r requirements.txt` | succeeded |
| `bash verify.sh` | FAIL at step 0: `system libhts.so.3 not found` |
| `bash native/build.sh` | FAIL: `cannot find -l:libhts.so.3` |
| V1 pytest subset, no native backends | 66 passed, 77 skipped |
| `make -C native_r2 test-tail` / `test-frame` (htslib 1.21 supplied by hand) | 40,186 checks, 0 failures / passed |
| nindels unit tests | 153 passed, 0 failed |

## Run 1 - htslib 1.21 (existing conda-style package prefix), scratch copy

`bash verify.sh full`: 40,477 passed, 0 failed, 2 skipped; WARN that htslib 1.21 is not the 1.23.x series. (First attempt failed in the V2 group
because `test_blake2b.py` wrote to the author's `/mnt/archive`; fixed, see `PORTABILITY_AUDIT.md` B1.)

## Run 2 - htslib 1.23.1 built by `scripts/build_local_htslib.sh` (working tree) - **Tested**

`bash scripts/build_local_htslib.sh` ran to completion (51 s; the pinned SHA-256 of the tarball matched, otherwise the script aborts).
Installed in `third_party/htslib` (git-ignored). Then `HTS_PREFIX=$PWD/third_party/htslib bash verify.sh full` (1 min 55 s):

| Group | Passed | Failed | Skipped | Detail |
|---|---|---|---|---|
| v1 (V1 C backends + pytest) | 120 | 0 | 2 | `115 passed, 2 skipped, 21 xfailed, 5 xpassed` - identical to the documented expectation; the 2 skips are the documented ones |
| v2 (V2/V2.x kernels bit-identical to V1 NumPy/SciPy) | 17 | 0 | 0 | 17 `[PASS]` |
| r2 | 40,187 | 0 | 0 | `40186 checks, 0 failures` + `FRAME TESTS PASSED` (frame suite prints no count; counted as 1) |
| nindels | 153 | 0 | 0 | `153 passed, 0 failed`; built with htslib |
| **Total** | **40,477** | **0** | **2** | `COMPLETE OFFLINE VERIFICATION PASSED (the optional real-data test was not run)` |

Header/library provenance: `verify.sh` reported `htslib 1.23.1 via HTS_PREFIX=...third_party/htslib` and PASS for "matches the 1.23.x series".
`ldd` of `native_r2/build/dnav2`, `native_v2x/build/dnav2`, `nindels/build/nindels-call`, `native/pileup_native*.so`, `native/reads_native*.so` shows
`libhts.so.3` resolved from `third_party/htslib/lib` in all five; `strings` finds no reference to the 1.21 prefix in any of them (no fallback).
The binaries start with `LD_LIBRARY_PATH` unset (embedded rpath). `libhts.so.3` reports version 1.23.1 through `hts_version()`.
All native backends compile: V1 C backends (2), V2.x engine, R1.1/R2 engine, nindels (core, CLI with `call`/`estimate-str`, tests).

## Run 3 - extracted source archive, clean directory - **Tested**

`bash scripts/make_source_archive.sh --worktree` -> 429 files, 2.1 MB tar.gz (5.0 MB extracted; largest file 1.3 MB, `ai_A_C_chr8_300x.vcf.gz`, which `--data` needs).
Extracted into an empty directory (no `.git`, `research/`, `data/`, `experimental/`, `checkpoints*/`), same venv, same `HTS_PREFIX`:
`bash verify.sh full` -> **identical result**: 40,477 passed, 0 failed, 2 skipped (1 min 44 s), `COMPLETE OFFLINE VERIFICATION PASSED`.
The archive contains the README, `verify.sh`, `environment.yml`, `Dockerfile`, `scripts/`, `mk/`, `docs/` (without `docs/history`), all `src/ include/ tail/ tests/`
directories of `native/`, `native_v2/`, `native_v2x/`, `native_r2/`, `nindels/`, the V2.x models and `research/stage4/frozen_R2.model`.
Note: the default (committed) mode archives `HEAD`; until the portability changes are committed it lacks them. Use `--worktree` before committing.

## Run 4 - R3.1 branch (`.worktrees/r3.1`, commit 31e14df) - **Partially tested**

Built with `-DHTSLIB_ROOT=third_party/htslib` (that branch has the older CMake, which supports `HTSLIB_ROOT`): `nindels-tests` 156 passed, 0 failed
(153 + 3 new); `tests_varion` (Python front end) 22 passed. Not run: R3.1 on real BAMs, the final test matrix (incomplete by the branch's own records:
`results/r31/test_report_partial.*` reaches only HG007 chr19 at 10x/20x), `verify.sh` from that branch (it does not contain the portability changes).

## Scientific outputs unchanged

Old build system (Makefiles from `HEAD`, htslib flags on the command line, `LD_LIBRARY_PATH`) vs new build system, same sources, input
`data/giab_chr21/chr21_slice.{bam,fa,vcf.gz}` shipped in the repo, `--mode cascade --threads 4`:

| Engine | `ai_cascade_t.vcf` (3,781 calls) | `routed_loci_t.tsv` | htslib 1.21 | htslib 1.23.1 |
|---|---|---|---|---|
| R2 (`frozen_R2.model --r2`) | sha256 0f5163fd... | sha256 2e6f2e6f... | old = new, byte-identical | old = new, byte-identical, same hash as 1.21 |
| V2.x A+C (`A_C.model`) | sha256 af3b576f... | sha256 2e6f2e6f... | old = new, byte-identical | old = new, byte-identical, same hash as 1.21 |

Smoke-scale (a synthetic 5 Mb slice), not a re-run of any published benchmark. `git diff HEAD` touches no file under `src/`, `include/`, `tail/`, `tools/`,
no `*.model`, nothing in `research/`, `results/`, `data/`, `experimental/`, `datasets/`. The only test files touched are the three `test_blake2b.py`
(scratch directory `/mnt/archive/...` -> `$V2_WORK_DIR` or the system temp dir; the assertions are unchanged; no test was deleted or skipped).
Build files changed (sha256 prefix `HEAD` -> now): `native_r2/Makefile` 7bc02d24->d08e7f77, `native_v2x/Makefile` 681d489c->98cd64e2,
`native_v2/Makefile` a1fe34ef->5e1f1436, `native/build.sh` cfd79d78->f1752cf8, `nindels/CMakeLists.txt` 8dec9435->e92ba353.
`research/stage4/R2_LOCK.md` locks a per-file hash of `native_r2/Makefile` in an external file; that one entry no longer matches
(`verify_identity.sh` already could not pass at the current `HEAD`: it expects commit `3cdfe5a` and files on `/mnt/hdd`).

## Other modes of `verify.sh` (htslib 1.21 run)

`python` -> 105 passed, 17 skipped, banner PARTIAL. `r2 nindels` -> 40,340 passed, 0 failed, banner PARTIAL listing `v1 v2 data` as not run.
`env` -> diagnostics, banner PARTIAL. `bogus` -> exit 2. `--data` without samtools/bcftools -> FAIL in step 0 before any download.
With no htslib configured, `verify.sh full` stops in step 0 with install instructions.

## Blocked or not run

| Item | Status | Exact reason |
|---|---|---|
| Conda/Mamba environment from `environment.yml` (solve, htslib version, samtools/bcftools, build, `verify.sh full`) | **Blocked, not tested** | no conda/mamba on the machine; downloading the official micromamba binary and running it was denied by the sandbox twice (not worked around). Package availability of htslib/samtools/bcftools 1.23.1 on bioconda (linux-64, linux-aarch64, osx-64, osx-arm64) was verified from channel metadata only; a Python 3.14 + `pip -r requirements.txt` solve is unverified. |
| `verify.sh --data` (real GIAB download) | Not run | samtools and bcftools not installed; installing them needs a package manager action not taken |
| Docker image | Not run | no container runtime |
| Fedora 43, Ubuntu, macOS, WSL2, native Windows | Not run | not available |
| `HTS_LIB` soname fallback `-l:libhts.so.3` (runtime-only installs) | Not exercised | the tested prefixes ship `libhts.so` |
| Any published benchmark re-run | Not run | out of scope |
