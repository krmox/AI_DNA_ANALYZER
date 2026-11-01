# Platform support

"Tested" means the exact commands were executed on that platform by the person who wrote this
file. Everything else is labelled **UNTESTED** and may fail. Test details: `BUILD_AND_TEST_REPORT.md`.

Legend: **Tested** = commands executed and passed on that platform; **Partially tested** = some routes or versions only;
**Documented but untested** = instructions exist, never executed; **Unsupported**.

| Platform | Status | What was actually run |
|---|---|---|
| Arch Linux / CachyOS x86-64 (GCC 16.2, Python 3.14) | **Tested** (offline verification) | `bash verify.sh full` = 40,477 passed, 0 failed, 2 skipped with htslib **1.23.1** built by `scripts/build_local_htslib.sh`, from the working tree and from an extracted source archive; also with an htslib 1.21 prefix. **Not run here:** `--data` (no samtools/bcftools), conda, Docker. |
| Fedora 43 | **Documented but untested** after the changes | It was the original reference environment; the legacy layout (system `htslib` in `/usr/lib64`) is the last-resort fallback and was not re-run. |
| Ubuntu 22.04/24.04 LTS | **Documented but untested** | Distro `libhts-dev` is older than the 1.23 series (recalled, not checked); `Dockerfile` builds htslib 1.23.1 from source (never built). |
| macOS (Intel / Apple silicon) | **Documented but untested** | `native/build.sh` has a Darwin branch; bioconda has 1.23.1 for osx-64/osx-arm64 (metadata only). |
| Windows 10/11 + WSL2 (Ubuntu) | **Documented but untested** | Expected to follow the Linux path; keep data inside the WSL filesystem. |
| Conda/Mamba route (`environment.yml`) | **Documented but untested** | Could not be run (see `BUILD_AND_TEST_REPORT.md`, "Blocked"). |
| Project-local htslib (`scripts/build_local_htslib.sh`) | **Tested** on CachyOS only | Completed, SHA-256 verified, used for the full verification. Other OSes not tested. |
| Docker / Podman | **Documented but untested** | No container runtime available; does not claim bit-identical results. |
| Native Windows (MSVC/MinGW) | **Unsupported** | Not attempted: POSIX headers (`pthread.h`, `sys/resource.h`), bash/make drivers, no supported MSVC htslib build. Use WSL2 or Docker. |

## Dependency model (read this first)

| Component | What it is | Used by |
|---|---|---|
| Python packages (`requirements.txt`) | numpy, scipy, torch (CPU), scikit-learn, pytest, **pysam** | V1 Python code and tests |
| htslib **inside the pysam wheel** | a private static copy (1.23.1 for pysam 0.24.0) linked into pysam's extension modules | Python code only. **Not** a `libhts.so`, not usable by the C/C++ engines |
| **System/conda htslib** (`libhts` + `htslib/*.h`) | shared library and headers | `native/` (V1 C backends), `native_v2*/`, `native_r2/`, `nindels-call` |
| samtools, bcftools | CLI tools | only `verify.sh --data` (slicing remote data) |

The native build uses **one** htslib prefix for headers *and* library (`scripts/hts_env.sh`), searched in this order:
`HTS_PREFIX` → active conda env → `pkg-config htslib` → legacy (pysam headers + system lib).
Explicit overrides: `HTS_INC` (directory containing `htslib/`), `HTS_LIBDIR`, `HTS_LIB` (full linker flags), `HTS_RPATH`.
Binaries get an rpath to the library directory, so `LD_LIBRARY_PATH` is not needed. Check what will be used:
`bash scripts/hts_env.sh report` or `make -C native_r2 hts-info`.

## Install recipes

**A. Conda/Mamba (UNTESTED - could not be run on the audit machine; recommended on paper because it supplies htslib 1.23.1)**
```bash
conda env create -f environment.yml && conda activate varion
pip install -r requirements.txt
bash verify.sh full
```

**B. System packages (distro htslib)**
```bash
# Fedora 43 (original environment)  sudo dnf install gcc gcc-c++ make cmake python3-devel htslib htslib-devel samtools bcftools curl
# Debian/Ubuntu (UNTESTED)          sudo apt install build-essential cmake python3-dev python3-venv libhts-dev samtools bcftools curl
# Arch/CachyOS: htslib was not in the configured repositories on the audit machine -> use A or C
python3 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
bash verify.sh full
```

**C. Project-local htslib 1.23.1 (opt-in; TESTED on CachyOS only)** - downloads a pinned, checksum-verified source tarball and compiles it into `third_party/htslib` (needs zlib/bzip2/xz/libcurl dev headers):
```bash
bash scripts/build_local_htslib.sh && export HTS_PREFIX=$PWD/third_party/htslib && bash verify.sh full
```

**D. Docker (UNTESTED)** – see the header of `Dockerfile`.

**Tested recipe** (as run on the audit machine; works with any htslib prefix, 1.23.1 via recipe C was the final run):
```bash
python3 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
export HTS_PREFIX=/path/to/prefix   # contains include/htslib/sam.h and lib/libhts.so
bash verify.sh full
```

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `htslib headers not found` / `library not found` | no htslib on the search path → recipe A/C, or `export HTS_PREFIX=…` |
| `cannot find -lhts` / `-l:libhts.so.3` | only a runtime lib is installed (no `libhts.so` dev symlink); install the `-devel`/`-dev` package or use a prefix |
| `error while loading shared libraries: libhts.so.3` | binary built elsewhere; rebuild (the rpath is set at link time) or set `LD_LIBRARY_PATH` |
| `undefined reference to tbx_itr_next1` (old builds) | headers and library from different htslib versions; now prevented, rebuild with `make clean` |
| V1 tests: many skips | native backends not built; `verify.sh v1` fails when more than 2 tests skip |
| `WARN htslib X differs from 1.23.x` | allowed; published numbers are for 1.23.x. The equivalence tests (`v1`, `v2`) must still pass |
| `nindels-call` lacks `call` | built without htslib (CMake warning); set `HTS_PREFIX`/`CONDA_PREFIX`; `verify.sh nindels` fails in this case |
| `--data` fails early | install `samtools`, `bcftools`, `curl` (conda env includes them) |
