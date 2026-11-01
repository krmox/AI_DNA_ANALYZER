# Portability audit (written before the changes)

Scope: installation, dependency resolution, native compilation, tests. No algorithm, threshold,
statistic or result file is in scope. Method: `git grep` for `lib64|libhts|HTS_|rpath|/mnt/`,
reading every Makefile/CMake/shell build script, and following the README in a clean scratch copy
(CachyOS, Python 3.14, GCC 16.2). "Verified" = reproduced on that machine.

## A. Observed findings re-checked against the repository

| # | Reported finding | Verdict | Evidence |
|---|---|---|---|
| 1 | Only Fedora 43 tested | **Confirmed** | README "Fedora 43 (tested) / Debian (untested)"; no CI, no Dockerfile, no `environment.yml` existed. |
| 2 | `native/build.sh` ignores `HTS_INC`/`HTS_LIB` | **Confirmed** | `native/build.sh` (old lines 10-14): headers from `pysam.get_include()`, link line hard-coded `-L/usr/lib64 -l:libhts.so.3`. |
| 3 | Hard-coded `/usr/lib64` | **Confirmed, wider than reported** | `native_r2/Makefile` lines 18,22,28,33,47,56,72; `native_v2/Makefile` 18,22,28,33; `native_v2x/Makefile` 18,22,28,33,47; `native/build.sh` 14. Both a `-L` and an `-Wl,-rpath,/usr/lib64` in every link line. |
| 4 | Custom htslib needs runtime library path | **Confirmed** | rpath was fixed at `/usr/lib64`; overriding `HTS_LIB` linked fine but the binary then failed to load `libhts` without `LD_LIBRARY_PATH` (reproduced). |
| 5 | `verify.sh --data` checks dependencies too late | **Not confirmed as stated** | `verify.sh` old line 28 already checked `samtools/bcftools/curl` in step 0. Real weaknesses: presence only (no version or functionality), the htslib check aborted before reporting other gaps, no OS/compiler/cmake diagnostics. |
| 6 | Repository ~1.4 GB | **Confirmed, with nuance** | See section D. Clone transfer ≈ 490 MB (pack), working tree ≈ 1.3 GB, `git archive` tar.gz 857 MB; the local checkout is larger (≈ 4 GB) only because of two local worktrees. |

## B. Additional issues found

| ID | File and location | Root cause | OS affected | Correction | Changes science? |
|---|---|---|---|---|---|
| B1 | `native_v2x/tests/test_blake2b.py:9` (identical copies in `native_v2/`, `native_r2/`) | `tempfile.mkdtemp(dir="/mnt/archive/AI_DNA_ANALYZER_v2/work")` – the author's disk. `make -C native_v2x test` (verify.sh step 3) failed on every other machine even though all 16 kernel checks passed. | all except author's | temp dir from `$V2_WORK_DIR`, default system temp | No (test scratch location only) |
| B2 | `native/build.sh:13` | `-I /usr/include/python3.14` hard-coded Python version; `-march=native` not overridable | all / macOS arm64 | derive from `sys.version_info`; `NATIVE_ARCH_FLAG` override (default unchanged) | No |
| B3 | `native/build.sh:13-14` | `-l:libhts.so.3` is GNU-ld syntax; macOS needs `-undefined dynamic_lookup` for Python extensions and uses `.dylib` | macOS | portable `-lhts` when a dev symlink exists, soname fallback otherwise; Darwin flag | No |
| B4 | `native_*/scripts/vec_report.sh:5` | htslib headers from pysam only | all | uses `scripts/hts_env.sh` unless `HTS_INC` set | No |
| B5 | `native_*/Makefile:4` | `HTS_INC` from `pysam.get_include()` requires pysam importable in the *build* Python, and pairs pysam's htslib 1.23.1 headers with whatever system library exists (version skew → link errors such as `tbx_itr_next1`) | all | headers and library now always come from the same prefix | No |
| B6 | `nindels/CMakeLists.txt:36-47` | htslib only via `-DHTSLIB_ROOT` or pkg-config; silently builds a CLI without `call`/`estimate-str` otherwise | all | also `$HTS_PREFIX`, `$CONDA_PREFIX`, default paths; louder warning; `verify.sh` fails if htslib is missing | No |
| B7 | `verify.sh:33-44` (old) | requires pysam-bundled htslib minor == system libhts minor; this is a symptom of B5, not a scientific requirement | all | replaced by header/library consistency + a WARN when the series differs from 1.23.x | No |
| B8 | `verify.sh` (old) | single fixed sequence; R2 and nindels never verified; skip counts of V1 tests not checked (77 skips on a machine without native backends would still print PASS) | all | groups, exact counts, skip budget (see below) | No |
| B9 | `native_r2/Makefile` (no `.gitignore` entry) | `native_r2/build/` is not git-ignored (only `native_v2*/build/`) | all | added | No |
| B10 | `nindels/scripts/env.sh`, `native_*/scripts/env.sh`, `run_*.sh`, `native_r2/analysis/*.py`, `scripts/data/*.sh` | benchmark-campaign scripts hard-code `/mnt/archive`, `/mnt/hdd`, `$HOME/varion_mamba` | all | **Not changed**: these reproduce the author's recorded runs and the README already says to edit the variables; changing them risks altering recorded provenance | n/a |
| B11 | `native_r2/Makefile` per-file hash | `research/stage4/R2_LOCK.md` records a per-file sha256 of `native_r2/{src,include,tail,Makefile}` in an external file | all | Makefile edit is unavoidable (B3/B5); the 20 other locked files are untouched. Old/new hashes in `BUILD_AND_TEST_REPORT.md`; outputs verified byte-identical | No (build file only) |
| B12 | `requirements.txt` | pins `torch==2.12.1` etc. from a CPU wheel index; fine on Linux x86-64 and (Python 3.14) installs; macOS arm64 / aarch64 wheels not checked | macOS, aarch64 | documented, not changed | No |

## C. Root-cause summary

Native code was always compiled against **pysam's bundled htslib headers** but linked against **the
system `libhts.so.3`**, with the system path fixed to Fedora's `/usr/lib64`. pysam wheels ship their
own statically-linked htslib inside Python extension modules; they provide neither a `libhts.so`
nor a usable header set for the C++ engines. Any machine whose distro lacks `htslib` at exactly
Fedora's path, or has a different htslib series, fails or silently mixes versions.

## D. Repository size (measured)

| Quantity | Size | Notes |
|---|---|---|
| Git pack (what `git clone` transfers) | ≈ 490 MB | `git count-objects -vH`; loose objects add ≈ 400 MB locally |
| Tracked tree at HEAD | 1,287 MB | `research/` 712 MB, `checkpoints*/` ≈ 440 MB (abandoned neural runs), `experimental/` 83 MB, `data/` 20 MB |
| `git archive` of everything | 1,290 MB tar / 857 MB tar.gz | what GitHub's auto-generated "Source code" archive contains |
| Lightweight source archive (`scripts/make_source_archive.sh`) | ≈ 2 MB tar.gz | code, tests, models, the 3 files `verify.sh` reads; **no** evidence, checkpoints, data |
| Reference data for `verify.sh --data` | ≈ 60 MB download | HG003 chr8 250 kb slice, truth, reference region |
| Full benchmark data | 100s of GB | GIAB 300× BAMs etc.; never in git |

Largest history blobs: `experimental/head_to_head/shared_ref/chr21_chrname.fa` 47 MB, `data/giab_chr21/chr21_slice.bam` (3 versions, 13 MB each), HG00x `.bai` files (now removed from the tree, still in history), DeepVariant/GATK `happy` outputs. **Nothing was removed or rewritten**; shrinking the clone would need history rewriting (out of scope) or moving `checkpoints*/` to a release asset (author's decision). `.gitattributes` `export-ignore` was deliberately **not** used: the Zenodo release archive is generated from the same machinery and would silently lose the evidence.

## E. Planned changes (all infrastructure)

`scripts/hts_env.sh` (single htslib resolver), `mk/htslib.mk`, edits to the three engine Makefiles,
`native/build.sh`, `vec_report.sh` ×3, `nindels/CMakeLists.txt`, the `test_blake2b.py` temp dir,
`verify.sh` (groups/diagnostics/counts), `environment.yml`, `Dockerfile`, `.dockerignore`,
`scripts/build_local_htslib.sh` (opt-in), `scripts/make_source_archive.sh`, README install section,
`.gitignore`. No file under any `src/`, `include/`, `tail/`, `tools/`, model, result or truth path changes.

## F. Status after implementation and verification (2026-11-01)

| ID | Issue | Status |
|---|---|---|
| A1-A4 | hard-coded `/usr/lib64`, ignored `HTS_*`, rpath | **Fixed and tested** (CachyOS, htslib 1.21 and 1.23.1; five linked artefacts resolve the intended library, no LD_LIBRARY_PATH) |
| A5 | `--data` dependency checks | Existing step-0 check kept and extended; early failure **tested** (without samtools/bcftools); the download itself **not run** |
| A6 | repository size | **Measured, not changed**; lightweight archive **tested** (2.1 MB, clean build + `verify.sh full` identical result) |
| B1 | `/mnt/archive` scratch dir in `test_blake2b.py` | **Fixed and tested** |
| B2-B4 | Python version / `-march=native` / macOS flags / vec_report | Linux **tested**; macOS branch **documented but untested**; `vec_report.sh` edit not executed (profiling helper, needs vectorisation tooling) |
| B5-B7 | pysam-header/system-lib skew, nindels lookup, version gate | **Fixed and tested** |
| B8 | verify.sh groups, skip budget, counts | **Tested** (all groups and modes except `data`) |
| B9 | `.gitignore` | Done |
| B10 | campaign scripts with `/mnt/...` paths | **Deliberately unchanged** (recorded provenance); still not portable |
| B11 | `native_r2/Makefile` in the R2 lock | Hash of that one file changed; outputs byte-identical (1.21 and 1.23.1) |
| B12 | requirements on macOS/aarch64 | Unverified |
| new | `make_source_archive.sh` only archived committed `HEAD` | Added `--worktree` mode; both modes **tested** |
| new | conda environment / Docker | **Not verified** (blocked / no runtime) |
