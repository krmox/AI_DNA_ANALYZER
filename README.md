# AI_DNA_ANALYZER (Varion)

### A cascaded statistical SNP caller (R1.1 → R2) and an experimental indel caller (`nindels`), evaluated on GIAB samples

![scope](https://img.shields.io/badge/scope-SNP%20%2B%20indel%20(experimental)-blue) ![lang](https://img.shields.io/badge/Python%20%2B%20C%2B%2B20%2Fhtslib-3776ab) ![platform](https://img.shields.io/badge/platform-Linux%20x86--64-lightgrey) ![release](https://img.shields.io/badge/latest%20tag-v2.0.0-informational) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22978775.svg)](https://doi.org/10.5281/zenodo.22978775)

A germline variant caller for Illumina short reads (GRCh38), built from statistical components and benchmarked against GIAB truth sets. It has two separate engines: a **SNP caller** (`dnav2`, V1 → R2) and an **indel caller** (`nindels`, C++20). It is a research project with frozen protocols and a manuscript, **not a production or clinical caller**.

> **On the name.** "AI" is historical: the project began as a neural sequence model. On real reads every neural stage was null or negative against the statistical baselines. No version described here contains a neural network.

> **Status.** `v2.0.0` (tag + Zenodo DOI) covers V1, V2 and V2.x. **R1.1 and R2 are newer than that tag** and are not part of the DOI release. R2 is frozen (commit `3cdfe5a`) and its validation is in `research/final_validation/`. The indel caller `nindels` is on `main` (merged from the `nindels-test` branch) and is **experimental**. **R3.1** (depth-adaptive indel configuration + a `varion` front-end command) lives on the separate branch `r3.1`; it is **not merged into `main`** and its final test run is not complete, so it is described below but not in the results.

## Versions at a glance

| Version | What it is | Changes the calls vs. the one before? | Headline result |
|---|---|---|---|
| **V1** (v1.0.0) | Python + C/htslib. A fixed-error binomial screen calls every locus; a router sends the few loci near its decision boundary to a per-read Poisson-binomial (PB) model. Three constants, never changed. | Reference behaviour | Routing kept F1 within 0.0027 of PB-everywhere in eight experiments; ten true SNPs lost. |
| **V2** | C++20 engine `dnav2`, same decision logic, multi-threaded. | No (byte-identical VCFs on HG002 chr20 300×) | ~2 min on 8 threads vs 5 h 10 min for the V1 record. |
| **V2.x A+C** | V2 + read-level evidence + a 32-input logistic classifier fitted on HG002 chr20. | Yes | Pre-registered HG003 chr8 run: F1 0.989012. |
| **R1.1** | V2.x with a corrected call domain: BED-edge windows (`--edge-windows`) and a truth-free frame (`--truth-free`), both on by default. | Yes (fixes two frame defects, see below) | Baseline for R2. |
| **R2** | `native_r2/`: R1.1 plus per-sample self-calibration, an out-of-distribution guard, a pattern tail (59 features instead of 51) and a retrained model (`research/stage4/frozen_R2.model`). Enabled with `--r2`. | Yes | See "R2 results". Without `--r2`, output is byte-identical to R1.1. |

Why R1.1 exists: R1 only scanned 64-bp tiles fully inside one BED interval (68–73 % of false negatives on HG003 chr8 were in dropped edge tiles), and it masked positions where the truth VCF has an indel, which flattered precision. R1.1 fixes both, so numbers from R1.1/R2 are **not directly comparable** with V1/V2.x numbers.

---

## Indel caller: `nindels` (experimental) and R3.1

`nindels/` is an independent indel caller, separate from the SNP engine and with its own VCF. It finds candidates from CIGAR I/D events in the BAM (≥ 2 supporting reads, alt fraction ≥ 0.03) and scores every candidate with a banded PairHMM. It does **not** call SNPs. Build and unit tests: see `nindels/README.md` (`cmake -S nindels -B build && cmake --build build -j`).

| Configuration | What it is | Flags |
|---|---|---|
| **Baseline** | Original biallelic genotyper, one indel error rate everywhere. Kept only as the frozen reference; not recommended. | `nindels-call call …` |
| **A** | Joint multi-allelic genotyping: candidates at one repeat locus are genotyped together, so a read from a third allele no longer counts as reference. | `--genotyper multi` |
| **B** | Repeat-specific PairHMM gap-open, estimated from the BAM being called (no truth labels). Only valid together with A. | `--str-model <table>` (made by `nindels-call estimate-str`) |
| **A+B** | A plus B. | both |

Why A and B exist (measured in `nindels/REPORT.md`): the baseline had low precision, falling with depth (3.5× more false positives than bcftools at 30×, 7× at 67× on HG004 chr20). Most errors were multi-allelic sites genotyped one allele at a time, and homopolymers/tandem repeats treated with one error rate. The pattern memory and the learned logistic model from the original spec are implemented but **not used** in any evaluated configuration (they need labelled training data).

Run (R3.1-free, as on `main`):

```bash
# A  (joint multi-allelic genotyping)
nindels-call call --bam B --ref R --region chr20 --force-reassess --genotyper multi --out-vcf a.vcf

# A+B  (also estimate the repeat error table from the same BAM)
nindels-call estimate-str --bam B --ref R --region chr20 --out S.tsv
nindels-call call --bam B --ref R --region chr20 --force-reassess --genotyper multi --str-model S.tsv --out-vcf ab.vcf
```

### Indel results (hap.py 0.3.15 xcmp, GIAB v4.2.1, PASS calls, indel F1)

| Data | Baseline | A | A+B | bcftools 1.21 |
|---|---|---|---|---|
| HG004 chr20, 30× (development) | 0.8743 | **0.9438** | 0.9397 | 0.9330 |
| HG004 chr20, 67× (native, development) | 0.8390 | 0.9516 | **0.9548** | 0.9548 |
| HG001 chr19, 30× (independent check, 1 run) | 0.8918 | 0.9609 | **0.9637** | 0.9552 |

Sources: `nindels/REPORT.md`, `nindels/results/HG004_chr20_regression_tables.md`, `nindels/results_validation_hg001_chr19.md`. The A vs A+B differences are small and **none of these three comparisons has a confidence interval** in the `main` results; do not read them as "A+B is better". bcftools is ahead of or equal to nindels at 67× and is a general-purpose caller that also calls SNPs; nindels is much faster in wall time (about 30 s per contig at 30×, 8 threads) but uses more memory (about 0.5–1 GB vs 0.1 GB).

### R3.1 (branch `r3.1`, not merged, final test incomplete)

R3.1 asks one question: is there a sequencing depth above which A+B is reliably better than A? It adds a depth-adaptive switch and a front-end command:

```bash
python3 -m varion indel --bam B --ref R --out-dir out/ --region chr20 --mode auto   # or A / A+B
python3 -m varion snp   --bam B --ref R --bed BED --out-dir out/ --contig chr20 --sample S
```

- `--mode auto` samples the contig depth and uses A below `global_min_depth`, A+B above it; a per-locus gate applies B only where the locus has enough reads. The choice is logged and written to the VCF header.
- The protocol was written before any R3.1 run (`docs/R3.1_PROTOCOL.md` on that branch); thresholds are chosen on HG004 chr20/chr21 only, then frozen. Test data: HG001 chr17, HG007 chr17/chr19, HG003 chr8.
- Development result: A+B is **worse** than A below ~30× (HG004 chr20 ΔF1 −0.049 at 10×, −0.0045 at 30×, 90 % intervals exclude 0) and only about +0.002–0.003 at 40× to native, with intervals touching 0. The frozen global threshold is **67.4×**; no local-depth bin met the rule. So R3.1 mostly means "use A, switch to A+B only at very high depth".
- Final test: the partial HG007 chr19 report only reaches 10× and 20×. At 20× R3.1 chose A (F1 0.9610), as intended. The remaining test matrix has not been completed and no claim about the test set is made here.

Details and the old-command → new-command table are in `docs/ARCHITECTURE.md` and `docs/COMMAND_MIGRATION.md` on the `r3.1` branch.

---

## Quick start for researchers

Three levels, from "does it build" to "does it reproduce the paper". Each lists the time, what it needs and what you should see.

| # | Command | Time | Needs | Expected last line / output |
|---|---|---|---|---|
| 1 | `bash verify.sh` | ~2 min | system prerequisites (below) | `ALL CHECKS PASSED` |
| 2 | `bash verify.sh --data` | ~3 min to ~40 min (download speed) | + internet, ~60 MB | `ALL CHECKS PASSED`; V2.x calls identical to the archived HG003 chr8 VCF |
| 3 | R2 build + unit tests (below) | ~1–2 min | same toolchain | `40186` checks, 0 failures |

Level 1–2 cover V1, V2 and V2.x. **`verify.sh` does not cover R1.1/R2**; use level 3 for those.

### 0. Prerequisites (once)

```bash
# Fedora 43 (tested)
sudo dnf install gcc gcc-c++ make python3-devel htslib samtools bcftools curl
# Debian/Ubuntu (untested; see the htslib note)
sudo apt install build-essential python3-dev python3-venv libhts3 samtools bcftools curl
```

### 1. Clone and install

```bash
git clone https://github.com/krmox/AI_DNA_ANALYZER.git
cd AI_DNA_ANALYZER
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` pins the versions used for the results. It includes CPU PyTorch and scikit-learn only because the V1 data layer imports `torch` and one V1 test uses `sklearn`; no neural network is used. The clone is about 1.4 GB (frozen result files are in the repo).

> **htslib note.** Native code is compiled against the htslib headers bundled with pysam and linked to the system `libhts.so.3`; the versions must match. pysam 0.24.0 bundles htslib 1.23.1 (matches Fedora 43). With a mismatched pysam the engine fails to link (e.g. `tbx_itr_next1`). `verify.sh` checks this first. Otherwise install the pysam release matching your system htslib, or override the paths: `make -C native_r2 engine HTS_INC=/path/to/htslib/include HTS_LIB=/path/to/htslib/lib`.

### 2. Run the checks

```bash
bash verify.sh            # offline: build everything, run the tests
bash verify.sh --data     # also re-run V2.x on real GIAB data
```

| Step | What it shows |
|---|---|
| 0 | Prerequisites present, htslib versions match |
| 1 | V1 native backends and the V2/V2.x engine compile |
| 2 | V1 tests: expected `115 passed, 2 skipped, 21 xfailed, 5 xpassed` (skips need a local HG002 slice; xfails are documented adversarial cases) |
| 3 | V2 numeric kernels bit-identical to the frozen V1 NumPy/SciPy code |
| 4 (`--data`) | Downloads HG003 300× chr8:20,000,001–20,249,984 + truth + reference, runs V2.x A+C, checks calls are **identical** to `research/final_independent_eval/run/ai_A_C_chr8_300x.vcf.gz` |

Each step prints PASS or FAIL. Logs: `data/verify_logs/`; step-4 downloads: `data/verify_hg003_chr8/` (git-ignored). The download is resumable: re-run `bash verify.sh --data` if it stops. Step 4 reproduces calls, not accuracy (F1 on 250 kb means little).

### 3. Build and test R1.1 / R2

```bash
make -C native_r2 engine        # builds native_r2/build/dnav2
make -C native_r2 test-tail     # R2 tail unit tests: expect 40,186 checks, 0 failures
make -C native_r2 test-frame    # call-domain (frame) tests
```

Not covered by these commands: full `make test` (needs `torch`), and the ASan/TSan/valgrind gates (only ASan/UBSan were run by the author).

### 4. Run a caller on your own data

Pick a version by flags. Same command shape for all; the binary is `native_r2/build/dnav2` (R1.1 and R2) or `native_v2x/build/dnav2` (V2/V2.x).

```bash
mkdir -p out                                   # the output directory must exist

COMMON="--bam sample.bam --ref GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
  --bed benchmark.bed --truth benchmark.vcf.gz \
  --contig chr19 --contig-len 58617616 --sample SAMPLE --tag chr19 \
  --threads 8 --mode cascade"

# R1.1 (default model, no --r2): corrected call domain
native_r2/build/dnav2 run $COMMON

# R2 (frozen): self-calibrating, OOD-guarded, pattern tail
native_r2/build/dnav2 run $COMMON --model research/stage4/frozen_R2.model --r2

# V2.x A+C (older accuracy layer)
native_v2x/build/dnav2 run $COMMON --model native_v2x/models/A_C.model
```

Outputs: `out/ai_cascade_<tag>.vcf` (calls), `out/summary_<tag>.json`, `out/routed_loci_<tag>.tsv`. R1.1/R2 do not read the truth VCF for calls (`--truth-free` is the default), but `--truth` and `--bed` are still the scoring frame, so `dnav2` is a **benchmark tool, not a stand-alone caller that takes only a BAM**. R2 may fall back to the plain binomial caller when its OOD guard fires; the reason is written in the VCF header.

### 5. Score a VCF (SNP, GIAB)

```bash
docker run --rm -v "$PWD":/w quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0 \
  hap.py /w/truth.vcf.gz /w/out/ai_cascade_chr19.vcf -f /w/benchmark.bed \
  -r /w/ref.fna -o /w/out/happy --engine=vcfeval --pass-only
```

Read the `SNP` / `PASS` row of `out/happy.summary.csv`. (The published runs used hap.py 0.3.15; Podman works the same way.)

---

## Input data (public, not in git)

| Used for | Item | Source |
|---|---|---|
| V2, V2.x, final eval | HG002/HG003 NHGRI Illumina 300× NovoAlign BAMs (2×148 bp), GRCh38 | GIAB `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/<sample>/…/NHGRI_Illumina300X_AJtrio_novoalign_bams/` |
| R2 validation | HG007 (NA24695) HiSeq 2×148 100× BAM, GRCh38 | GIAB `…/ChineseTrio/HG007_NA24695_mother/…` (script: `research/final_validation/fetch_hg007.sh`) |
| R2 validation | HG001 natives | `research/final_validation/expanded/fetch_hg001_nvme.sh` |
| V1 experiments | HG002/3/4 Illumina 2×250 BAMs; Ensembl release 110 GRCh38 FASTA | GIAB; `ftp.ensembl.org/pub/release-110/` |
| all | Reference `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna` | NCBI GenBank `seqs_for_alignment_pipelines.ucsc_ids/` |
| all | Truth VCFs + confident-region BEDs, NIST v4.2.1 GRCh38 | GIAB `release/<trio>/<sample>/NISTv4.2.1/GRCh38/` |
| scoring | hap.py 0.3.15 | `quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0` |

Checksums: `research/final_independent_eval/checksums.txt`, `research/stage4/R2_LOCK.md`, `research/final_validation/evidence_package/`.

## Where each published result comes from

| Result | Scripts | Recorded outputs |
|---|---|---|
| **R2 frozen identity** (source/binary/model hashes, build recipe) | `research/final_validation/expanded/verify_identity.sh` | `research/stage4/R2_LOCK.md` |
| **R2 fresh-data validation** (HG007 chr17/chr19) | `research/final_validation/` (`fetch_hg007.sh`, `hg007/`) | `research/final_validation/REPORT.md`, `hg007/SUMMARY.md` |
| **R2 expanded benchmark** (HG001 + HG007, 10–50×, DeepVariant on 2 cells) | `research/final_validation/expanded/` (`run_expanded.py`, `matrix_expanded.sh`, `analyze_expanded.py`) | `research/final_validation/evidence_package/` (`RESULTS.md`, `expanded_benchmark_report.md`, `claims_audit.md`, `reproduction_guide.md`) |
| R2 blind-new (HG006 chr16) | `research/stage4/` | `research/stage4/BLIND_HG006_chr16.md` |
| Final evaluation HG003 chr8 (V2.x) | `research/final_independent_eval/scripts/` | `research/final_independent_eval/` |
| V2.x development/holdout (HG002 chr20) | `native_v2x/scripts/`, `native_v2x/analysis/` | `research/V2X_*`, `research/V2X_FINAL_REPORT.md` |
| V2 equivalence/timing | `native_v2x/scripts/run_gates.sh`, `bench_full.sh` | `research/V2_*` |
| DeepVariant/GATK chr20 300× | commands in `research/benchmark_chr20_300x/BENCHMARK_CHR20_300X_FINAL.md` | `research/benchmark_chr20_300x/happy/` |
| V1 regional experiments v13–v20 | `fetch_*.sh`, `prepare_*.sh`, `run_*_extract.sh`, `bench_v13…v20*.py` | `results/bench_v13…v20/` |

Full-run scripts record the author's absolute paths (`/mnt/archive/...`); edit the variables at the top of each script. The R2 reproduction steps, with hashes, are in `research/final_validation/evidence_package/reproduction_guide.md`.

A small V1 run (Python path, HG002 chr20 15×, 200 kb, ~35 s):

```bash
python3 extract_bench_v12.py --fasta data/reference/chr20_full.fa \
  --bam data/giab_hg002_chr20/chr20_15x.bam --vcf data/giab_hg002_chr20/chr20.vcf.gz \
  --bed data/giab_hg002_chr20/chr20_highconf.bed --contig chr20 \
  --region 33000000 33200000 --chunk-bp 500000 --features none --out cache/smoke.npz
mkdir -p out && PIPE_OUT=out PIPE_CACHE=cache/smoke.npz PIPE_TAG=smoke \
  python3 experimental/chr20_validation/run_chr20_pipeline.py
```

Expected: `Wrote cache/smoke.npz | 193792 loci | 203 SNP`, then `out/ai_cascade_smoke.vcf` and `out/ai_pb_only_smoke.vcf`.

---

## Results

All numbers are SNP-only, GRCh38, GIAB v4.2.1 high-confidence regions, hap.py 0.3.15 unless noted. Full tables and intervals are in the manuscript and the evidence packages.

### R2 results (frozen R2 vs. R1.1 vs. bcftools 1.21)

Expanded benchmark on **HG001 and HG007, chr17 and chr19** (2-Mb block-bootstrap 95 % CIs; one contig-level CI, not a between-genome interval):

| Comparison, SNP F1 at 30× | Pooled (4 cells) |
|---|---|
| R2 − R1.1 | +0.00134 [+0.00073, +0.00211] |
| R2 − bcftools | +0.00249 [+0.00141, +0.00398] |

- R2 beats bcftools at every tested depth, by a margin that shrinks from **+0.035 F1 at 10×** to **+0.002 at 30–50×**.
- R2 reaches F1 ≥ 0.990 at ≥ 30× in all 8 tested cells.
- R2's gain over R1.1 is real (CI excludes 0) but small: **the pre-registered target of +0.003 over R1.1 was not met** (anywhere it was tested).
- DeepVariant 1.10.0 was **more accurate than R2** on the two cells where it ran (HG001/HG007 chr19, 30×; F1 0.99559 / 0.99652 vs 0.99315 / 0.99379), at about 150× the compute. No superiority claim over DeepVariant is possible.
- Cost (8 threads, i7-6700): R2 ≈ R1.1 (~10 s per contig, ~216 MB); bcftools is ~5.6× slower in wall time.

### V2.x A+C by evidence tier

| Tier | Data | Precision | Recall | F1 |
|---|---|---|---|---|
| Development (in-sample) | HG002 chr20, even 2-Mb blocks | 0.998662 | 0.980511 | 0.989503 |
| Holdout | HG002 chr20, odd 2-Mb blocks | 0.996362 | 0.981781 | 0.989017 |
| Final evaluation (one pre-registered run) | HG003 chr8 | 0.997567 | 0.980603 | 0.989012 |

HG003 is HG002's father, same platform/aligner/read length, so this separates sample and chromosome but not platform or ancestry. HG003 chr8 later became R2 training data, so it is **not** independent for R2.

### External callers on HG002 chr20 at 300× (descriptive, not a ranking)

DeepVariant 1.10.0 F1 0.997613, GATK 4.6.2.0 0.991712, V2.x A+C 0.989271 (includes the fitting blocks), V1/V2 0.972429. DeepVariant and GATK also call indels.

### Speed (V2 vs. V1, HG002 chr20 300×)

About 17.5× on 1 thread and 102× on 8 threads for the same work; 155× for the cascade on 8 threads, which also does less work.

---

## Limitations

- **The SNP caller (R2) is SNP-only.** MNPs and structural variants are neither called nor scored. Indels are handled by the separate, experimental `nindels` caller, whose evidence is weaker (see its section): 3 comparisons on 2 samples, no confidence intervals on `main`, and no pre-registered independent test yet.
- **No sample is untouched any more.** HG002/3/4 trained or tuned V2.x/R2, HG005 was used for debugging, HG006 chr16 and now HG001/HG007 have been evaluated. The R2 expanded benchmark is a **replication, not a blind test**; an unseen germline GIAB genome is still needed.
- Scope of the R2 evidence: two samples × two chromosomes, confident regions only, no whole-genome run, one lab/platform/aligner family (Illumina HiSeq, NovoAlign), GRCh38.
- Seeds 101–103 are subsamples of one BAM, not replicates. 50× has one seed; DeepVariant has one seed and 2 of 4 cells (chr17 cells unassessed).
- bcftools baselines used 1.23.1 (Stage 3, HG006) or 1.21 (HG001/HG007). Same command line, different versions.
- Timings: one machine, mostly single runs; one resource-guardrail cell (HG001 chr19 50×, +18 % wall, equal CPU) failed and was not repeated.
- About 1.4–3.5 % of truth SNPs fall in windows dropped by the V1/V2.x frame; R1.1 recovers most of them, so R1.1/R2 recall is not comparable with older figures.
- V1 constants are specific to ~15× depth. V2.x was fitted on one sample/chromosome/depth/platform.
- The V1 native counts backend is not exact on adversarial BAMs with repeated read names (matches pysam on all validated real data).
- HG005 raw artefacts were lost; those numbers are record-only.
- Environment: the Makefile defaults assume `/usr/lib64` for libhts (use `HTS_INC`/`HTS_LIB`); some scripts hard-code `/mnt/archive`; no lockfile for the native toolchain.

Full lists: `research/LIMITATIONS_AND_OPEN_QUESTIONS.md`, `research/final_validation/evidence_package/limitations.md`, `claims_audit.md`, and the manuscript (Table 17).

---

## Repository layout

```text
AI_DNA_ANALYZER/
├── verify.sh                  one-command check for V1/V2/V2.x (not R1.1/R2)
├── cascade.py, binomial_baseline.py, quality_error_model.py   V1 frozen logic
├── pileup_counts.py, read_level_pileup.py, extract_bench_v12.py   V1 extraction
├── native/                    V1 C/htslib backends and build.sh
├── native_v2/                 V2 engine as frozen at the V2 milestone
├── native_v2x/                V2 / V2.x / R1.1 engine (dnav2), models/, tests/, scripts/
├── native_r2/                 R2 engine: calibration, OOD guard, pattern tail, tests/
├── nindels/                   indel caller (C++20): src/, scripts/, REPORT.md, results/
├── bench_v13 … bench_v20      V1 validation experiments
├── experimental/              genotype layer (Method C), chr20 and HG005 pipelines
├── results/                   frozen V1 outputs
└── research/
    ├── stage4/                R2 prereg, lock, frozen_R2.model, HG006 chr16 blind result
    ├── final_validation/      R2 fresh-data campaign: HG007, HG001, DeepVariant, claims audit
    ├── final_independent_eval/  V2.x HG003 chr8 pre-registered evaluation
    └── manuscript, LaTeX, figures, reports
```

`checkpoints*/` and `model_*.py` come from the earlier neural experiments and are not used.

## Publication and citation

- Manuscript: `research/latex/paper.pdf` (source `research/PAPER_MANUSCRIPT_V2.md`); covers up to V2.x
- Software v2.0.0 (V1, V2, V2.x): DOI [10.5281/zenodo.22978775](https://doi.org/10.5281/zenodo.22978775)
- Historical V1 release v1.0.0: DOI [10.5281/zenodo.22895746](https://doi.org/10.5281/zenodo.22895746)
- Claim–evidence map: `research/CLAIM_EVIDENCE_MAP.csv`; AI-assistance disclosure: `research/AI_ASSISTANCE_DISCLOSURE.md`
- Citation metadata: `CITATION.cff`

## License

GPL-3.0-only. See `LICENSE`.
