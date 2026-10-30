# DeepVariant feasibility audit (2026-10-29)

## Findings
| Question | Finding |
|---|---|
| Was DeepVariant run before? | Yes: v1.10.0 (`google/deepvariant:1.10.0`) in Docker for HG002 chr20 300× (`research/benchmark_chr20_300x`, 1 h 28 m wall, CPU only) and a chr21 15× head-to-head. Only HG002 — no held-out sample. |
| Is a runtime available now? | **No.** docker, podman, singularity, apptainer are not installed. |
| Can the old container image be reused? | **No.** Its storage is `/mnt/hdd/docker-data`, `containerd-data`, `AI_DNA_ANALYZER_benchmark/container_storage*/docker` — root-owned (`drwx--x---`/`drwx------`), unreadable without sudo. The `.img` needs a root loop mount. Not attempted (no sudo, per standing rules). |
| conda/bioconda route? | **Not usable.** `bioconda::deepvariant` is a ~0 MB `noarch` stub (1.8–1.10) with no binaries; the only linux-64 build is 1.5.0 on Python 3.6 / TensorFlow-GPU 2.0 (different version, unmaintained stack). |
| Source build? | Bazel + TensorFlow build; impractical, not attempted. |
| **User-space OCI route** | **Feasible and tested.** The registry manifest for `google/deepvariant:1.10.0` (digest `sha256:962e5a83b1d76aae6990625d47102785f791603f2138aa1fa9aa4fb6a2eecbe6`, config `sha256:96e40015…`) has 12 layers, 2.1 GB compressed; all layers were fetched with curl and their sha256 verified (12/12 `LAYER_OK`), unpacked as an ordinary user to `~/varion_data/dv/rootfs` (4.7 GB), and executed with `bubblewrap` (`/usr/bin/bwrap`, user namespaces work). No root, no system change. |
| Same version as before? | Yes, 1.10.0 (`run_deepvariant --version` → "DeepVariant version 1.10.0"). Model: `/opt/models/wgs`. |
| Smoke test | HG007 chr19:20,000,000–21,000,000, native BAM, `--model_type=WGS --num_shards=8`: exit 0, 101 s wall (shared with a download), 3,077 records, hap.py 0.3.15 accepted the VCF: SNP PASS TP 1779 / FN 11 / FP 1, F1 0.99664 (a 1-Mb sanity check only; not a result). |
| Evaluator | The same hap.py 0.3.15 (conda env `hp`), same truth, BED, reference as the other callers; DV VCF read directly. |
| Resources | i7-6700, 4C/8T, 31 GB RAM, no GPU (CPU TensorFlow with AVX2/FMA). Earlier: 1 h 28 m for a 64-Mb contig at 300×. Estimate for 30× chr17 (83 Mb) / chr19 (59 Mb): roughly 0.5–1.5 h each (not measured on these contigs). 4 runs at 30× fit a working day; per-run abort limit in the preregistration: 4 h wall or > 24 GB RSS. Disk: image 4.7 GB (+2.1 GB blobs, deletable) and a few GB of intermediates per run. |

## Differences that cannot be removed (to state in any comparison)
- Container image, not a bioconda build: the exact image used before, executed through bubblewrap rather than Docker (no GPU either way). Behaviour is expected to be identical; this was not cross-checked against a Docker run.
- DV is given the whole contig (`--regions CONTIG`); Varion restricts calls to the confident BED; hap.py scores inside the BED only for all callers.
- DV is a SNP+INDEL caller; only SNP PASS rows are compared. Varion is SNP-only.
- Resource figures come from `wait4` rusage of the whole `bwrap` process tree; CPU and RSS include the TensorFlow stages; they are not directly comparable to Varion's single binary and no efficiency ratio will be claimed beyond this machine.

## Decision
DeepVariant is **feasible**; it will be run as preregistered (30×, seed 101, 4 cells, optional 10×/20×), in a separate table. No claim of superiority/non-inferiority vs DeepVariant is preregistered. Final status is recorded in `expanded_benchmark_report.md`.

## Outcome (2026-10-29, appended after the runs)
DeepVariant ran successfully on **HG007 chr19** (1,879 s wall, 10,885 s CPU, 8.9 GB peak) and **HG001 chr19** (1,688 s, 10,659 s, 9.2 GB) at 30× seed 101, scored with the same hap.py/truth/BED/reference. The two chr17 runs were **not completed: stopped by the user (too slow)** — recorded as UNASSESSED in `FAILED_RUNS_EXPANDED.md`. Measured cost: ~28–31 min per chr19 run (59 Mb), against the 0.5–1.5 h estimate. Results: `DEEPVARIANT_TABLE.md`, `expanded_benchmark_report.md` §5.
