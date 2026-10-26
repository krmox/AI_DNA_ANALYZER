# Varion R1.1 — implementation lock for Stage 3 (2026-10-26)

Nothing below may change during Stage 3. Any deviation is a new version and must be reported as such.

## Source identity
| Item | Value |
|---|---|
| Git commit | `659df8298a82ccb6f3db542ff15a28d51215a3b4` ("Add Varion R1.1 call domain: BED-edge windows and truth-free frame") on `main` |
| Working tree | clean at freeze (`git status --short --untracked-files=all` empty); only new, untracked files under `research/stage3/` were added afterwards (scripts, plan, frozen model copy). No change to `native_v2x/`. |
| Source hashes | `native_v2x/{src,include,Makefile}` per-file sha256 in `/mnt/archive/AI_DNA_ANALYZER_stage3/gates/source_sha256.txt` (17 files); concatenated-content sha256 `a498f42d7fad198c0545a5f339952cfdbc6b647f2f7c5b16c2ff63b7704ee3c6` |
| Build | `make -C native_v2x BUILD=/mnt/archive/AI_DNA_ANALYZER_stage3/build engine kernels test-frame` |
| Compiler / flags | g++ 15.3.1 (Fedora 43); `-O3 -std=c++20 -Wall -Wextra -fPIC -ffp-contract=off -fno-fast-math -fno-associative-math -frounding-math`, generic x86-64 (no `-march`) |
| Dependencies | htslib 1.23.1 (`/usr/lib64/libhts.so.3`, dynamic), libstdc++, zlib, bz2, lzma; pysam 0.24.0 headers (compile only) |
| Executable | `/mnt/archive/AI_DNA_ANALYZER_stage3/build/dnav2`, sha256 `e83d87fd2f810eb7a9d096fb797c6d3451aaef86d14a0071dc2ad9b91a3a39df` |
| Correspondence | The freshly built binary is **bit-identical** to the Stage 2 R1.1 binary `AI_DNA_ANALYZER_stage2/build_r11/dnav2` (same sha256), i.e. the documented source reproduces the executable. |
| Model | `research/stage3/frozen/R1.model` (copy of `AI_DNA_ANALYZER_depth_robust/R1.model`), sha256 `63816707c4af7794138c3917ddcc22ddce67f800f724db8a004ad60ee4b711ed` |
| Frozen constants | Binomial threshold 7.0, router cutoff 5.411872376933351, PB threshold 10.5, Method C eps 0.01 (unchanged since Stage 1; compiled in) |

## Pipeline / flags (R1.1 defaults)
```
dnav2 run --bam BAM --ref REF --bed CONTIG_CONFIDENT_BED --out-dir D --contig C --contig-len L \
          --sample S --tag T --threads 8 --mode cascade --model R1.model
```
No `--truth` argument: `--truth-free` (default) never reads truth. `--edge-windows` is the default. The confident BED defines the call domain (input, not a label). Output: SNP-only VCF, converted with `bcftools view -Oz`, scored unchanged by hap.py.

## Machine (all Stage 3 measurements)
Intel Core i7-6700 (4 cores / 8 threads, max 4.0 GHz, `intel_pstate` governor `powersave`), 31 GB RAM, Fedora Linux 43, kernel 6.19.12, all data on an HDD (`/mnt/archive`, ext/NTFS-class mount, 21 GB free at start); root SSD full. Desktop session with browser and other processes running (load average 5–7 before the first run). `--threads 8` for Varion; bcftools baseline sharded over 8 parallel processes.

## Regression gates run on the freshly built binary (all before any evaluation run)
| Gate | Result |
|---|---|
| `make test-frame` (`FRAME TESTS PASSED`) | pass |
| `tests/test_kernels.py`, `tests/test_blake2b.py` (kernel library built to stage3 build dir) | pass (all `[PASS]`) |
| R1 compatibility: `--strict-windows --truth-mask` on HG002 chr20 odd 30x (`HG002_chr20_test2_f010_s7.bam`, HG002 truth) vs the archived Stage 1 run `depth_robust/runs/t2_chr20_30x_R1` | VCF body identical (69,891 records); the header differs only by the added `##varion_mode=edge_windows=0,truth_free=0` line |
| `scripts/test_truth_free.sh` on HG004 chr21 30x slice | PASS: default VCF byte-identical with real vs nonexistent truth file (3,208 records); R1 mode differs (3,206 records) |
| Default mode vs Stage 2 H1+H2 VCF (HG004 chr21 ~30x) | evaluated as part of the anchor run (see `RESULTS.md`) |
Outputs: `/mnt/archive/AI_DNA_ANALYZER_stage3/gates/`.

## Defect policy
A defect that invalidates an experiment is documented and the experiment stopped; the caller is not patched in this stage.
