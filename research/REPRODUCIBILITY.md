> **SUPERSEDED IN PART (2026-09-20).** Update 2026-09-20: the `/tmp/claude-…/hg005-fix-worktree` path below no longer exists (wiped by reboot). Current locations: main working tree (uncommitted) and `.claude/worktrees/opt-chr20`; experiment scripts now resolve the project root from their own location. Lost raw HG005 post-fix artefacts are listed in `PAPER_MANUSCRIPT.md` §5.

# Reproducibility

## Code

- Repository: `AI_DNA_ANALYZER` (local; no public remote configured in this
  session).
- Main branch HEAD at session start: `a3d5761b9ad74a04d0344d94b4ec088c3440aff6`.
- This dossier and the HG005 extraction fix were produced in an isolated
  worktree: `fix/hg005-extraction-integrity`, path
  `/tmp/claude-1000/-home-mark-Documents-Projects-AI-DNA-ANALYZER/hg005-fix-worktree`.
  The main checkout was never modified (`git diff --stat` against
  `extract_bench_v12.py`/`pileup_counts.py` in the main checkout is empty).
- Files modified from `main` HEAD: `extract_bench_v12.py` (11 lines),
  `pileup_counts.py` (50 lines). Files added: `test_hg005_extraction_integrity.py`.
  No commits were made; `git status --short` in the worktree shows the
  expected 2 modified + untracked new files only.

## Hardware

- CPU: Intel Core i7-6700 @ 3.40GHz, 4 cores / 8 threads (`nproc` = 8).
- RAM: 31 GiB total.
- OS: Linux (Fedora), kernel `6.19.12-200.fc43.x86_64`.
- Disk: local NVMe, `/home` partition ~94G total, ~11G free at time of the
  HG005 rerun (safety margin maintained above the ≥4.5GB requirement
  throughout).

## Software versions

| Component | Version |
|---|---|
| Python | 3.14.3 |
| pysam | 0.24.0 |
| numpy | 2.4.4 |
| scipy | 1.17.1 |
| samtools / bcftools | 1.23.1 |
| gcc | 15.3.1 (20260722, Red Hat 15.3.1-1) |
| Docker | present; `quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0`, digest `sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0` |
| HTSlib | as vendored by pysam 0.24.0 (no standalone HTSlib build; no C extension in this repo — see `PIPELINE_INTEGRITY_AUDIT.md` I-2) |

## Datasets and checksums (HG005 chr1:1,000,001-4,000,000 experiment)

| File | sha256 |
|---|---|
| `data/giab_hg005_stress/chr1_24Mb_30x.bam` | `e3f2637b2659e10864bd2a43a065f225da29bd8fa612336b8c2752c049e34963` |
| `data/giab_hg005_stress/chr1_24Mb.vcf.gz` | `10e886563bdb9b5487a6ef9fac733524aa94e8d58a57c92231edaf009868dfab` |
| `data/giab_hg005_stress/chr1_24Mb_highconf.bed` | `a8031902d2c5045c45fc19df647e995217a7f2f353fcb8c5b2443a556f012642` |
| `experimental/stress_test/shared_ref/chr1_chrname.fa` | `749354128afbf131d512ed4daaee2a13aa60f3ddf983f90e8039931a974dae74` |

`cascade.py` sha256 `b0ee9f4b24fe06dc4d677ca80fdd9bedaa9885e52038679812f1f550c5286131`;
`experimental/genotype_layer/method_c_regression.py` sha256
`b5cea1f1660eb472190267c743f78ef4c330872b6c4fbede9f136f4cca5882d2` — both
verified byte-identical to the pre-fix session's recorded values (proves no
threshold or genotype-formula change occurred anywhere in this work).

For HG002/HG003/HG004 experiments (`bench_v13`-`bench_v19`), inputs are
documented per-script in `data/giab_hg002_v14/`, `data/giab_hg002_v15_test/`
(or equivalent), `data/giab_hg004_v19/` etc.; each script's own docstring and
`results/bench_v1*/region_selection.json` records the exact region-selection
protocol and timestamp.

## Exact commands (HG005 post-fix)

**Steps 2 and 3 are not currently runnable.** `run_hg005_pipeline_postfix.py` and
`run_happy_postfix_regionscoped.sh` were lost in the same 2026-09-20 `/tmp`-worktree wipe that took the
raw post-fix HG005 artefacts (see `FINAL_SCIENTIFIC_AUDIT.md` Finding 5, and
`LIMITATIONS_AND_OPEN_QUESTIONS.md` item A.13). Step 1 (extraction) and step 4 (regression tests) still
run against the current codebase. The HG005 numbers this pipeline once produced are preserved as a
record in `research/PAPER_MANUSCRIPT_V2.md` and flagged `RECORD_ONLY_RAW_UNAVAILABLE` in
`CLAIM_EVIDENCE_MAP.csv`; they cannot be independently regenerated from this repository.

```bash
# 1. Extraction (fixed extractor)
python3 extract_bench_v12.py \
  --fasta experimental/stress_test/shared_ref/chr1_chrname.fa \
  --bam data/giab_hg005_stress/chr1_24Mb_30x.bam \
  --vcf data/giab_hg005_stress/chr1_24Mb.vcf.gz \
  --bed data/giab_hg005_stress/chr1_24Mb_highconf.bed \
  --region 1000000 4000000 --contig chr1 --chunk-bp 500000 --features none \
  --out cache/hg005_stress_postfix/chr1_1000000_4000000.npz

# 2. Frozen cascade + Method C -> VCF
python3 experimental/stress_test/run_hg005_pipeline_postfix.py

# 3. hap.py (region-scoped, primary)
bash experimental/stress_test/run_happy_postfix_regionscoped.sh

# 4. Regression tests
python3 -m pytest test_hg005_extraction_integrity.py test_bench_v12.py \
  test_cascade.py test_cheap_router.py test_robustness_benchmark.py \
  test_bench_v13.py -q
```

No random seed is consumed by the extraction or cascade/Method-C path itself
(deterministic, proven in `HG005_COORDINATE_INTEGRITY_REPORT.md` Gate C).
Where a seed is used elsewhere in the repo (e.g. `robustness_benchmark.py`'s
paired bootstrap, `seed=20260812`; BAM downsampling, `samtools view -s 42.09`),
it is recorded per-experiment in that experiment's own report/devlog.

## Environment variables

None of the scripts used in this work read environment variables for
scientific parameters (thresholds, seeds, region, or file paths are all
CLI arguments or in-file constants). `LANG=C.UTF-8` was set inside the hap.py
Docker container only (standard biocontainers image default), not a
host-side dependency.

## What a third party needs to reconstruct this work

1. This repository at commit `a3d5761b9ad74a04d0344d94b4ec088c3440aff6` (or
   later, with the two-file diff from this session applied — see
   `git diff main fix/hg005-extraction-integrity -- extract_bench_v12.py
   pileup_counts.py`).
2. The four checksummed HG005 input files above (not redistributed in this
   dossier; obtainable from the GIAB FTP URLs recorded in
   `experimental/stress_test/HG005_MANIFEST.json`).
3. The software versions table above (pinned via a venv/conda environment;
   no `requirements.txt`/`environment.yml` currently exists in this repo —
   noted as a reproducibility gap, see `PEER_REVIEW_SELF_AUDIT.md`).
4. Docker, for the hap.py evaluation step.
5. ~12GB free disk for the regional BAM + hap.py image.
