# Varion (AI_DNA_ANALYZER)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22978775.svg)](https://doi.org/10.5281/zenodo.22978775) ![license](https://img.shields.io/badge/license-GPL--3.0--only-blue) ![platform](https://img.shields.io/badge/platform-Linux%20x86--64-lightgrey)

Varion is a research germline variant caller for Illumina short reads (GRCh38), built from statistical components and benchmarked against GIAB truth sets. It has two separate engines:

- a **SNP caller** (`dnav2`; versions V1 → V2 → V2.x → R1.1 → **R2**), and
- an **indel caller** (`nindels`, C++20; configurations Baseline / A / B / A+B, and the depth-adaptive **R3.1** layer).

**Intended use:** methods research and reproduction of the reported benchmarks. **It is not a production or clinical caller.** All results below are on GIAB high-confidence regions of a few samples and chromosomes; they are not guarantees for other data.

> **On the name.** "AI" is historical: the project began as a neural sequence model. On real reads every neural stage was null or negative against the statistical baselines, and no version described here contains a neural network. (`checkpoints*/` and `model_*.py` are from those abandoned experiments and are unused.)

> **Release status.** Tag `v2.0.0` (Zenodo DOI) covers V1, V2 and V2.x. **R1.1 and R2 are newer than that tag.** R2 is frozen (lock record `research/stage4/R2_LOCK.md`; the validation reports cite commit `3cdfe5a`). `nindels` Baseline/A/B/A+B is on `main`. **R3.1 lives only on branch `r3.1`** (not merged; final test run incomplete) — `python3 -m varion …` does not exist on `main`.

## Contents

1. [Versions](#versions) · 2. [Architecture](#architecture-and-algorithms) · 3. [Install](#installation-and-compilation) · 4. [Usage](#usage) · 5. [Formats](#inputoutput-formats) · 6. [Results](#benchmark-summary) · 7. [Limitations](#known-limitations) · 8. [Reproducibility](#reproducibility) · 9. [Layout](#repository-layout) · 10. [Citation](#citation) · 11. [License](#license)

## Versions

| Version | What it is | Changes calls vs. previous? |
|---|---|---|
| **V1** (`v1.0.0`) | Python + C/htslib. A fixed-error binomial screen calls every locus; a router sends loci near its decision boundary to a per-read Poisson-binomial (PB) model. Three frozen constants. | reference |
| **V2** | C++20 engine `dnav2`, same decision logic, multi-threaded. | No (byte-identical VCFs on HG002 chr20 300×) |
| **V2.x A+C** | V2 + read-level evidence + a 32-input logistic classifier fitted on HG002 chr20. | Yes |
| **R1.1** | V2.x with a corrected call domain: BED-edge windows (`--edge-windows`) and a truth-free frame (`--truth-free`), both default. | Yes |
| **R2** | R1.1 + per-sample self-calibration, an out-of-distribution (OOD) guard, a pattern tail (59 instead of 51 features) and a retrained model (`research/stage4/frozen_R2.model`). Enabled with `--r2`; without it output is byte-identical to R1.1. | Yes |

R1 only scanned 64-bp tiles fully inside one BED interval and masked positions where the truth VCF has an indel, which flattered precision. R1.1 fixes both, so **R1.1/R2 numbers are not directly comparable with V1/V2.x numbers.**

## Architecture and algorithms

**SNP caller (`dnav2`).** Per 64-bp window the engine builds base counts from the BAM (one fused streaming pileup) and runs a *cascade*: (1) a cheap binomial log-likelihood-ratio screen with a fixed error rate decides every locus; (2) a confidence router sends only loci near the decision boundary (a small fraction) to a per-read Poisson-binomial model that uses base qualities; (3) from V2.x on, a logistic classifier over read-level evidence features makes the final decision. R2 adds sample-level calibration of the error model, a guard that falls back to the plain binomial caller when the sample is out of distribution (reason recorded in the VCF header), and extra pattern features. Floating-point behaviour is fixed (`-ffp-contract=off`, no fast-math) so results are bit-reproducible against the V1 reference.

**Indel caller (`nindels`).** Candidates come from CIGAR I/D events in the BAM (≥ 2 supporting reads, alt fraction ≥ 0.03); each is scored with a banded PairHMM. It does not call SNPs.

| Configuration | Description | Flags |
|---|---|---|
| Baseline | Biallelic genotyper, one indel error rate everywhere. Frozen reference only; not recommended. | `call …` |
| **A** | Joint multi-allelic genotyping at repeat loci. | `--genotyper multi` |
| **B** | Repeat-specific PairHMM gap-open estimated from the BAM being called (no truth labels). Valid only with A. | `--str-model S.tsv` (from `estimate-str`) |
| **A+B** | Both. | both |

The pattern-memory and learned logistic components from the original spec are implemented but **not used** in any evaluated configuration. Rationale and error analysis: `nindels/REPORT.md`, `nindels/EXPERIMENTS.md`.

**R3.1 (branch `r3.1`).** Asks whether there is a depth above which A+B is reliably better than A. `--mode auto` samples contig depth and uses A below a frozen threshold (67.4×) and A+B above it; B is also gated per locus by read count. The protocol was written before any R3.1 run (`docs/R3.1_PROTOCOL.md` on that branch). Development finding: A+B is *worse* than A below ~30× (HG004 chr20 ΔF1 −0.049 at 10×, −0.0045 at 30×; 90 % intervals exclude 0) and only ≈ +0.002–0.003 at 40× to native depth, with intervals touching 0. The final test matrix is incomplete (the partial HG007 chr19 report reaches only 10× and 20×; at 20× R3.1 chose A, F1 0.9610). **No claim about the R3.1 test set is made.**

## Installation and compilation

Prerequisites (Fedora 43 tested; Debian/Ubuntu untested): `gcc gcc-c++ make cmake python3-devel htslib samtools bcftools curl`.

```bash
git clone https://github.com/krmox/AI_DNA_ANALYZER.git && cd AI_DNA_ANALYZER
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt        # CPU torch/scikit-learn only for V1 data layer and tests
```

The clone is ≈ 1.4 GB because frozen result files are in the repository.

**htslib.** The SNP engines compile against the htslib headers bundled with pysam and link the system `libhts.so.3`; versions must match (pysam 0.24.0 bundles htslib 1.23.1, matching Fedora 43). `HTS_INC` is taken from `pysam.get_include()`. `HTS_LIB` is a linker-flag string, default `-L/usr/lib64 -l:libhts.so.3`; override both for another htslib:

```bash
make -C native_r2 engine HTS_INC=/path/to/htslib/include HTS_LIB="-L/path/to/htslib/lib -l:libhts.so.3"
```

**Build targets**

```bash
make -C native_r2 engine         # -> native_r2/build/dnav2   (R1.1 and R2)
make -C native_r2 test-tail      # R2 tail unit tests
make -C native_r2 test-frame     # call-domain tests
make -C native_v2x engine        # -> native_v2x/build/dnav2  (V2 / V2.x)
cmake -S nindels -B nindels/build && cmake --build nindels/build -j   # nindels-call, nindels-tests
nindels/build/nindels-tests    # unit tests (needs no htslib)
bash verify.sh [--data]          # V1/V2/V2.x checks; does NOT cover R1.1/R2
```

`nindels-call call` and `estimate-str` need htslib (`-DHTSLIB_ROOT=…` or pkg-config); without it only the core library, tests and the table subcommands build.

## Usage

### SNP caller

```bash
mkdir -p out                                   # the output directory must exist
COMMON="--bam sample.bam --ref GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
  --bed benchmark.bed --truth benchmark.vcf.gz \
  --contig chr19 --contig-len 58617616 --sample SAMPLE --tag chr19 \
  --threads 8 --mode cascade"

native_r2/build/dnav2 run $COMMON                                                   # R1.1
native_r2/build/dnav2 run $COMMON --model research/stage4/frozen_R2.model --r2      # R2 (frozen)
native_v2x/build/dnav2 run $COMMON --model native_v2x/models/A_C.model              # V2.x A+C
```

`--bed` and `--truth` define the scoring frame. R1.1/R2 do not use the truth VCF for calls (`--truth-free` default), but `dnav2` still requires both, so it is a **benchmark tool, not a BAM-only caller**. Run `dnav2` without arguments for all options.

### Indel caller

```bash
# A: joint multi-allelic genotyping
nindels-call call --bam B --ref R --region chr20 --force-reassess --genotyper multi --out-vcf a.vcf

# A+B: also estimate the repeat error table from the same BAM
nindels-call estimate-str --bam B --ref R --region chr20 --out S.tsv
nindels-call call --bam B --ref R --region chr20 --force-reassess --genotyper multi --str-model S.tsv --out-vcf ab.vcf
```

On branch `r3.1` only: `python3 -m varion indel --bam B --ref R --out-dir out/ --region chr20 --mode auto` and `python3 -m varion snp …`.

### Scoring (hap.py)

```bash
docker run --rm -v "$PWD":/w quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0 \
  hap.py /w/truth.vcf.gz /w/out/ai_cascade_chr19.vcf -f /w/benchmark.bed \
  -r /w/ref.fna -o /w/out/happy --engine=vcfeval --pass-only
```

Read the `SNP` (or `INDEL`) / `PASS` row of `out/happy.summary.csv`. Published runs used hap.py 0.3.15 (Podman works too).

## Input/output formats

| | Format |
|---|---|
| Input | coordinate-sorted, indexed BAM (Illumina, GRCh38); indexed FASTA (`.fai`); GIAB high-confidence BED and truth VCF (`.vcf.gz` + `.tbi`) for `dnav2` |
| `dnav2` output | `out/ai_cascade_<tag>.vcf` (calls), `out/summary_<tag>.json`, `out/routed_loci_<tag>.tsv` (loci sent to PB); with `--mode pb_all` also a PB-only VCF |
| `nindels-call` output | VCF (`--out-vcf`; INFO `K`, `DP` declared in header), optional per-candidate table (`--out-table`) |
| Model files | `research/stage4/frozen_R2.model`, `native_v2x/models/*.model` (checksums in `research/stage4/R2_LOCK.md`, `research/V2X_MODEL_FILES_SHA256.csv`) |

## Benchmark summary

Sources are given for every number; full tables and intervals are in the cited files. SNP numbers: GRCh38, GIAB v4.2.1 high-confidence regions, hap.py 0.3.15, PASS.

### R2 (SNP) — source: `research/final_validation/evidence_package/`

Expanded benchmark on HG001 and HG007, chr17 and chr19, 10–50× (2-Mb block-bootstrap 95 % CIs; contig-level, not between-genome):

| SNP F1 at 30×, pooled over 4 cells | Δ [95 % CI] |
|---|---|
| R2 − R1.1 | +0.00134 [+0.00073, +0.00211] |
| R2 − bcftools 1.21 | +0.00249 [+0.00141, +0.00398] |

- R2 > bcftools at every tested depth; margin shrinks from +0.035 F1 at 10× to +0.002 at 30–50×. R2 F1 ≥ 0.990 at ≥ 30× in all 8 tested cells.
- R2's gain over R1.1 is real (CI excludes 0) but small: the **pre-registered target of +0.003 over R1.1 was not met**.
- DeepVariant 1.10.0 was **more accurate than R2** on the two cells where it ran (chr19 30×; HG001 F1 0.99559 vs R2 0.99315, HG007 0.99652 vs 0.99379) at ≈ 150× the compute (`DEEPVARIANT_TABLE.md`). No superiority claim over DeepVariant is possible.
- Cost (8 threads, i7-6700): R2 ≈ R1.1 (≈ 10 s per contig, ≈ 216 MB); bcftools ≈ 5.6× slower in wall time.

### V2.x A+C by evidence tier — sources: `research/V2X_FINAL_REPORT.md`, `research/final_independent_eval/FINAL_EVALUATION_REPORT.md`

| Tier | Data | Precision | Recall | F1 |
|---|---|---|---|---|
| Development (in-sample) | HG002 chr20, even 2-Mb blocks | 0.998662 | 0.980511 | 0.989503 |
| Holdout | HG002 chr20, odd 2-Mb blocks | 0.996362 | 0.981781 | 0.989017 |
| Final evaluation (one pre-registered run) | HG003 chr8 | 0.997567 | 0.980603 | 0.989012 |

HG003 is HG002's father on the same platform/aligner, so this separates sample and chromosome, not platform or ancestry. HG003 chr8 later became R2 training data, so it is **not independent for R2**.

### External callers, HG002 chr20 300× (descriptive, not a ranking) — `research/benchmark_chr20_300x/BENCHMARK_CHR20_300X_FINAL.md`

SNP F1: DeepVariant 1.10.0 0.997613, GATK 4.6.2.0 0.991712, V2.x A+C 0.989271 (includes fitting blocks), V1/V2 0.972429. DeepVariant and GATK also call indels.

### Speed, V2 vs. V1 (HG002 chr20 300×) — `research/V2_NATIVE_ENGINE_REPORT.md`

Against the V1 frozen record (18,578 s ≈ 5 h 10 min, not re-run): ≈ 17.5× on 1 thread, ≈ 102× on 8 threads (`pb_all`), ≈ 155× for the cascade on 8 threads (which also does less work). Treat ratios as indicative to ±15 %.

### `nindels` indel F1 (hap.py xcmp, PASS) — `nindels/REPORT.md`, `nindels/EXPERIMENTS.md`, `nindels/results/HG004_chr20_regression_tables.md`, `nindels/results_validation_hg001_chr19.md`

| Data | Baseline | A | A+B | bcftools 1.21 |
|---|---|---|---|---|
| HG004 chr20, 30× (development) | 0.8743 | **0.9438** | 0.9397 | 0.9330 |
| HG004 chr20, 67× native (development) | 0.8390 | 0.9516 | **0.9548** | 0.9548 |
| HG001 chr19, 30× (independent check, 1 run) | 0.8918 | 0.9609 | **0.9637** | 0.9552 |

The A vs A+B differences are small and **none of these comparisons has a confidence interval on `main`**; do not read them as "A+B is better". B without A is clearly worse than A (HG004 30×: 0.8871; 67×: 0.8697). bcftools is equal at 67× and is a general caller that also calls SNPs. `nindels` is faster in wall time (≈ 30 s per contig at 30×, 8 threads) but uses more memory (≈ 0.5–1 GB vs ≈ 0.1 GB).

## Known limitations

- **R2 is SNP-only**; MNPs and structural variants are neither called nor scored. The experimental `nindels` evidence is weaker: 3 comparisons on 2 samples, no confidence intervals on `main`, no pre-registered independent test.
- **No sample is untouched any more.** HG002/3/4 trained or tuned V2.x/R2, HG005 was used for debugging, HG006 chr16, HG001 and HG007 have been evaluated. The R2 expanded benchmark is a **replication, not a blind test**; an unseen germline GIAB genome is still needed.
- R2 evidence: two samples × two chromosomes, confident regions only, no whole-genome run, one lab/platform/aligner family (Illumina HiSeq, NovoAlign), GRCh38.
- Seeds 101–103 are subsamples of one BAM, not replicates. 50× has one seed; DeepVariant has one seed and 2 of 4 cells.
- bcftools baselines used 1.23.1 (HG006) or 1.21 (HG001/HG007); one resource-guardrail cell (HG001 chr19 50×, +18 % wall) failed and was not repeated.
- 1.4–3.5 % of truth SNPs fall in windows dropped by the V1/V2.x frame; R1.1/R2 recall is not comparable with older figures.
- V1 constants are specific to ≈ 15× depth; V2.x was fitted on one sample/chromosome/depth/platform. The V1 native counts backend is not exact on adversarial BAMs with repeated read names.
- HG005 raw artefacts were lost; those numbers are record-only (see `experimental/stress_test/INVALIDATED_prefix_coordinate_bug/`).
- R3.1 is not merged and its test run is incomplete (see above).
- Environment: Makefile defaults assume `/usr/lib64`; some scripts hard-code `/mnt/archive`; no lockfile for the native toolchain.

Full lists: `research/LIMITATIONS_AND_OPEN_QUESTIONS.md`, `research/final_validation/evidence_package/limitations.md`, `claims_audit.md`.

## Reproducibility

Input data are public and not in git: GIAB HG001–HG007 BAMs, NIST v4.2.1 GRCh38 truth VCFs and BEDs, and the `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna` reference (fetch scripts: `fetch_*.sh`, `research/final_validation/fetch_hg007.sh`, `research/final_validation/expanded/fetch_hg001_nvme.sh`). Checksums: `research/final_independent_eval/checksums.txt`, `research/stage4/R2_LOCK.md`, `research/final_validation/evidence_package/`.

| Level | Command | Expected |
|---|---|---|
| Build + unit tests | `bash verify.sh` (~2 min) | `ALL CHECKS PASSED`; V1 tests `115 passed, 2 skipped, 21 xfailed, 5 xpassed` |
| Re-run V2.x on real data | `bash verify.sh --data` (~3–40 min, ≈ 60 MB download) | calls identical to `research/final_independent_eval/run/ai_A_C_chr8_300x.vcf.gz` (calls, not accuracy) |
| R2 unit tests | `make -C native_r2 test-tail test-frame` | `40186 checks, 0 failures`; `FRAME TESTS PASSED` |
| `nindels` unit tests |  `nindels/build/nindels-tests` | `153 passed, 0 failed` |

| Result | Scripts | Recorded outputs |
|---|---|---|
| R2 frozen identity (hashes, build recipe) | `research/final_validation/expanded/verify_identity.sh` | `research/stage4/R2_LOCK.md` |
| R2 fresh-data validation (HG007 chr17/19) | `research/final_validation/` | `REPORT.md`, `hg007/SUMMARY.md` |
| R2 expanded benchmark | `research/final_validation/expanded/` | `research/final_validation/evidence_package/` (`reproduction_guide.md`) |
| R2 blind-new (HG006 chr16) | `research/stage4/` | `research/stage4/BLIND_HG006_chr16.md` |
| V2.x final evaluation (HG003 chr8) | `research/final_independent_eval/scripts/` | `research/final_independent_eval/` |
| V2.x dev/holdout (HG002 chr20) | `native_v2x/scripts/`, `native_v2x/analysis/` | `research/V2X_*` |
| V1 regional experiments v13–v20 | `fetch_*.sh`, `prepare_*.sh`, `run_*_extract.sh`, `bench_v13…v20*.py` | `results/bench_v13…v20/`, `logs_*/` |
| `nindels` regression / validation | `nindels/scripts/run_regression.sh`, `run_validation.sh` | `nindels/results/`, `nindels/results_*.md` |

Full-run scripts record the author's absolute paths (`/mnt/archive/...`); edit the variables at the top of each script. Development, validation and held-out data are kept apart by directory and by the pre-registrations (`research/stage4/PREREG_R2.md`, `research/final_independent_eval/`).

## Repository layout

The repository is a frozen research record: Python modules import each other from the root and result files are cited by path in the reports, so the layout is intentionally **not** reorganised.

```text
verify.sh                          one-command check for V1/V2/V2.x (not R1.1/R2)
*.py (root)                        V1 pipeline, benchmarks (bench_v12…v20), tests (test_*.py)
fetch_*.sh prepare_*.sh run_*.sh   data download / preparation per benchmark round
native/                            V1 C/htslib backends
native_v2/ native_v2x/ native_r2/  C++20 engines: V2, V2/V2.x/R1.1, R2 (src/, include/, tests/, scripts/)
nindels/                           indel caller (src/, include/, tests/, scripts/, results/, REPORT.md)
experimental/                      genotype layer, chr20/HG005 pipelines, head-to-head (Clair3 etc.)
results/  logs_*/                  frozen V1 outputs and run logs
History/                           development logs
research/                          manuscripts, audits, claim–evidence map, stage4 (R2), final_validation, final_independent_eval, latex
FINAL_*.md, VALIDATION_ROUND_2.md, METHOD_C_*.md, DEEPMIND_*.md   V1-era reports (cited by research/)
HG00*_benchmark.vcf.gz.tbi         truth-VCF indices (checksummed in research/final_independent_eval/logs/truth_sha256.txt)
```

## Citation

Cite via `CITATION.cff` (GitHub "Cite this repository").

- Software v2.0.0 (V1, V2, V2.x): DOI [10.5281/zenodo.22978775](https://doi.org/10.5281/zenodo.22978775)
- Historical V1 release v1.0.0: DOI [10.5281/zenodo.22895746](https://doi.org/10.5281/zenodo.22895746)
- Manuscript: `research/latex/paper.pdf` (source `research/PAPER_MANUSCRIPT_V2.md`); covers up to V2.x, **not** R1.1/R2/nindels
- Claim–evidence map: `research/CLAIM_EVIDENCE_MAP.csv`; AI-assistance disclosure: `research/AI_ASSISTANCE_DISCLOSURE.md`

## License

GPL-3.0-only. See `LICENSE` and `NOTICE`.
