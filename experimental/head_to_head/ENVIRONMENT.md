# Environment — Head-to-Head Benchmark

- Host: Linux fedora 6.19.12-200.fc43.x86_64, 8 CPUs, 31GiB RAM
- Disk: `/home` filesystem, 94GB total
- Git HEAD: `16b221ee82d64f66de25119aafa4e5944cfe63cc` (branch `main`), clean tree throughout
- Work performed in isolated worktree: `/tmp/head_to_head_work` (detached HEAD, `data/` symlinked
  from the main working tree, not copied)
- Python: 3.14.3 (host); Python 2.7.18 (built via `pyenv install 2.7.18` solely to run Strelka2's
  bundled python2 workflow scripts — this did not touch system Python)
- numpy 2.4.4, scipy 1.17.1, pysam 0.24.0
- samtools 1.23.1, bcftools 1.23.1
- Docker 29.4.1
- Java: OpenJDK 25.0.2 (used successfully with GATK 4.5.0.0 despite GATK's docs targeting
  older JDKs — no compatibility issue observed)
- No CUDA/GPU used for any caller (DeepVariant CPU image explicitly, no `-gpu` suffix)

## External caller versions/digests

| Caller | Version/tag | Image digest | Delivery |
|---|---|---|---|
| DeepVariant | 1.6.1 | `sha256:ccab95548e6c3ec28c75232987f31209ff1392027d67732435ce1ba3d0b55c68` | Docker (`google/deepvariant`) |
| Clair3 | v1.0.10 (retry; `:latest` was broken, see COMMANDS.md) | `sha256:57cf5d20f2ee39c1b91493ad1fb5c1b9fa838691efce818c3139caa5e6c6b974` | Docker (`hkubal/clair3`) |
| GATK | 4.5.0.0 | N/A (static zip release, not Docker) | `github.com/broadinstitute/gatk` release asset |
| Strelka2 | 2.9.10 | N/A (static tarball, not Docker) — BLOCKED before producing output | `github.com/Illumina/strelka` release asset |
| hap.py | 0.3.15 (`py27hcb73b3d_0`); 0.3.7 (`py27_1`) tried first and failed | `sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0` (0.3.15) | Docker (`quay.io/biocontainers/hap.py`) |

## Random seeds

`paired_bootstrap` (from `evaluate_quality_error.py`) used its module-default seed `20260812`,
unchanged — same seed already used by `bench_v12_stage1.py` for comparable prior results.
10,000 resamples, as in all prior bench_v* rounds.

## Disk-budget baseline at start of this task

13GB free / 94GB total (87% used) — already below the 20% (18.8GB) target margin, carried
forward from the prior FINAL_RESEARCH_REPORT.md pass. The user explicitly approved a sequential
pull→run→save→cleanup→verify protocol to work within this constraint rather than requiring the
full margin be restored first. See COMMANDS.md's "Sequential disk protocol log" table for the
full before/after trace.

## Reproducibility verification performed at the end of this task

```
git status --short   # only experimental/head_to_head/** untracked; nothing else
git diff              # empty (no tracked file modified)
```
Confirmed: no production/frozen file (`cascade.py`, `robustness_benchmark.py`,
`extract_bench_v12.py`, `evaluate_quality_error.py`, `residual_metrics.py`, `providers.py`, etc.)
was edited. All new content lives under `experimental/head_to_head/` (this task's own output) plus
the four files from the immediately prior FINAL_RESEARCH_REPORT.md task, both sets untracked and
uncommitted, as instructed.
