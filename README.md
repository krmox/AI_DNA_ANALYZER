# AI_DNA_ANALYZER

### A cascaded statistical SNP caller, a native execution engine and a lightweight accuracy layer, evaluated on GIAB samples

![scope](https://img.shields.io/badge/scope-SNP--only-blue) ![lang](https://img.shields.io/badge/Python%20%2B%20C%2B%2B20%2Fhtslib-3776ab) ![platform](https://img.shields.io/badge/platform-Linux%20x86--64-lightgrey) ![release](https://img.shields.io/badge/release-v2.0.0-informational) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22978775.svg)](https://doi.org/10.5281/zenodo.22978775)

A germline SNP caller for Illumina short reads (GRCh38), built from statistical components and benchmarked against GIAB truth sets. It is a research project with a frozen protocol and a manuscript, **not a production or clinical caller**.

> **On the name.** "AI" is historical: the project began as a neural sequence model. On real reads every neural stage was null or negative against the statistical baselines. No version described here contains a neural network.

## Three versions, kept separate

| Version | What it is | Changes the calls? | Main result |
|---|---|---|---|
| **V1** (frozen protocol, v1.0.0) | Python + C/htslib. A fixed-error binomial screen calls every locus; a router sends the few loci near its decision boundary (0.04–0.16 %) to a per-read Poisson-binomial (PB) model. Three constants, selected once and never changed. | Defines the reference behaviour | Routing kept F1 within 0.0027 of PB-everywhere in eight experiments (2 PRESERVED, 5 IMPROVED, 1 DEGRADED; v14 experiment-level verdict NEGATIVE). Ten true SNPs were lost. |
| **V2** (engine rewrite) | C++20 engine `dnav2`, same decision logic, multi-threaded, no read tensor. | No. VCFs byte-identical to V1 on HG002 chr20 at 300× | About 2 min on 8 threads vs 5 h 10 min for the V1 record (indicative; other day and cache state). |
| **V2.x A+C** (accuracy layer) | V2 plus cheap read-level evidence and a 32-input logistic classifier, fitted on the even 2-Mb blocks of HG002 chr20. Not part of the V1 protocol. | Yes | Holdout (odd blocks) F1 0.97195 → 0.98902. Pre-registered single run on HG003 chr8: F1 0.989012 (95 % CI 0.98713–0.99048). |

V2 and V2.x are the same binary (`native_v2x/build/dnav2`): without `--model` it is V2, with `--model native_v2x/models/A_C.model` it is V2.x A+C.

---

## Check that it works (one command)

Tested on a fresh clone with a clean Python venv on Fedora 43 (Python 3.14, system htslib 1.23.1). Other distributions should work if the htslib versions match (see below), but were not tested.

**1. Install the system prerequisites** (once):

```bash
# Fedora 43
sudo dnf install gcc gcc-c++ make python3-devel htslib samtools bcftools curl
# Debian/Ubuntu (untested; see the htslib note below)
sudo apt install build-essential python3-dev python3-venv libhts3 samtools bcftools curl
```

**2. Clone and install the Python packages:**

```bash
git clone https://github.com/krmox/AI_DNA_ANALYZER.git
cd AI_DNA_ANALYZER
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` pins the versions the results were produced with. It includes the CPU build of PyTorch and scikit-learn: no neural network is used, but the V1 data layer (`providers.py`) imports `torch` and a V1 test uses `sklearn`. The clone is about 1.4 GB, because frozen result files are kept in the repository.

> **htslib note.** The native code is compiled against the htslib headers bundled with pysam and linked to the system `libhts.so.3`. The two versions must match: pysam 0.24.0 bundles htslib 1.23.1, which matches Fedora 43. With a newer pysam (e.g. 0.24.1) the V2 engine fails to link. `verify.sh` checks this first and tells you if the versions differ. On another distribution, install the pysam release whose bundled htslib matches your system htslib, or install htslib 1.23 from source.

**3. Run the check:**

```bash
bash verify.sh           # offline: builds everything and runs the tests (~2 min)
bash verify.sh --data    # also re-runs V2.x on real GIAB data (~60 MB download; ~3 min in total on a good connection)
```

What each step checks:

| Step | What it shows | Needs |
|---|---|---|
| 0 | Prerequisites are present and the htslib versions match | — |
| 1 | V1 native backends and the V2/V2.x engine compile | — |
| 2 | V1 test suite: frozen cascade, native vs pysam equivalence. Expected on a fresh clone: `115 passed, 2 skipped, 21 xfailed, 5 xpassed` (the 2 skips need a local HG002 slice; the xfails are documented adversarial cases) | — |
| 3 | V2 numeric kernels are bit-identical to the frozen V1 NumPy/SciPy code (gammaln, log1p, binomial LLR, router, PB, Method C, BLAKE2b ordering) | — |
| 4 (`--data`) | Downloads HG003 300× reads for chr8:20,000,001–20,249,984, the GIAB v4.2.1 truth and the GRCh38 reference sequence for that region, runs V2.x A+C, and checks that the calls are **identical** to the archived final-evaluation VCF (`research/final_independent_eval/run/ai_A_C_chr8_300x.vcf.gz`) on that region | internet, samtools, bcftools, curl |

Each step prints PASS or FAIL, and the last line reads `ALL CHECKS PASSED` or `SOME CHECKS FAILED`. Logs are in `data/verify_logs/`, and step 4 keeps its downloads in `data/verify_hg003_chr8/` (both git-ignored).

The NCBI server speed varies: in testing, step 4 took under a minute on one day and about 40 minutes on another. The download is resumable: if it fails or you stop it, run `bash verify.sh --data` again and it continues where it stopped.

Step 4 reproduces calls, not accuracy: F1 on a 250-kb region would mean little. To reproduce the published accuracy numbers, see the next section.

---

## Reproducing the published results

The full runs need whole-chromosome 300× BAMs (the HG003 chr8 slice alone is 28 GB) and, for scoring, Docker or Podman with hap.py. The scripts record the exact commands of the published runs with the author's absolute paths (`/mnt/archive/...`); change the variables at the top of each script to your own paths.

### Input data (public, not in git)

| Used for | Item | Source |
|---|---|---|
| V2, V2.x, final evaluation | HG002 and HG003 NHGRI Illumina 300× NovoAlign BAMs (HiSeq, 2×148 bp), GRCh38 | GIAB, `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/<sample>/NIST_HiSeq_<sample>_Homogeneity-*/NHGRI_Illumina300X_AJtrio_novoalign_bams/` |
| V2, V2.x, final evaluation | Reference `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna` | NCBI GenBank `seqs_for_alignment_pipelines.ucsc_ids/` |
| V1 experiments | HG002/HG003/HG004 Illumina 2×250 NovoAlign BAMs, GRCh38 | GIAB, same tree |
| V1 experiments | Ensembl release 110 GRCh38 per-chromosome FASTA | `ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/` |
| all | Truth VCFs and confident-region BEDs, NIST v4.2.1 GRCh38 (`*_benchmark.vcf.gz`, `*_benchmark_noinconsistent.bed`) | GIAB `release/AshkenazimTrio/<sample>/NISTv4.2.1/GRCh38/` |
| scoring | hap.py 0.3.15 image `quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0` | biocontainers |

Checksums of the inputs used: `research/final_independent_eval/checksums.txt`, `research/final_independent_eval/logs/`, `experimental/chr20_validation/FROZEN_MANIFEST.json`, `experimental/stress_test/HG005_MANIFEST.json`.

### Running the engine on your own data

```bash
make -C native_v2x engine           # builds native_v2x/build/dnav2 (the plain `make` builds only the kernel library)
mkdir -p out                        # the output directory must exist
native_v2x/build/dnav2 run \
  --bam  HG003.GRCh38.300x_chr8.bam \
  --ref  GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
  --truth HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  --bed  HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  --contig chr8 --contig-len 145138636 --sample HG003 --tag chr8_300x \
  --threads 8 --mode cascade \
  --model native_v2x/models/A_C.model   # omit --model for V2 (= V1 output)
```

Outputs: `out/ai_cascade_<tag>.vcf` (calls), `out/summary_<tag>.json` (loci, routed loci, timings) and `out/routed_loci_<tag>.tsv`. `--mode pb_all` also writes the PB-only control VCF. The truth VCF is read only to reproduce V1's scored-locus frame; the decision logic never uses it. In other words, `dnav2` is a benchmark tool, not a stand-alone caller that takes only a BAM.

### Where each published result comes from

| Result | Scripts | Recorded outputs |
|---|---|---|
| Final evaluation, HG003 chr8 (V2.x) | `research/final_independent_eval/scripts/` (`download_bam.sh`, `run_final.sh`, `happy_final.sh`, `analyze_final.py`, `frame_audit.py`) | `research/final_independent_eval/` (pre-registration, exact command, VCF, hap.py output, report) |
| V2.x development and holdout (HG002 chr20) | `native_v2x/scripts/` (`env.sh` holds the paths; `run_all_finalists.sh`, `run_happy.sh`), `native_v2x/analysis/` | `research/V2X_*.csv`, `research/V2X_FINAL_REPORT.md` |
| V2 equivalence and timing | `native_v2x/scripts/run_gates.sh`, `bench_full.sh`; `native_v2x/tests/*_compare.py` | `research/V2_*` |
| DeepVariant / GATK on chr20 300× | no script; the exact commands are written out in `research/benchmark_chr20_300x/BENCHMARK_CHR20_300X_FINAL.md` | `research/benchmark_chr20_300x/happy/` |
| V1 chr20 15× (DEGRADED verdict) | `experimental/chr20_validation/run_extract_chr20.sh`, `run_chr20_pipeline.py`, `run_happy_chr20.sh`, `bench_v20_chr20.py` | `results/bench_v20/`, `experimental/chr20_validation/` |
| V1 regional experiments v13–v19 | `fetch_*.sh`, `prepare_*.sh`, `run_*_extract.sh`, `bench_v13…v19*.py` | `results/bench_v13…v19/` |

A small V1 run (Python path, HG002 chr20 15×, 200 kb, about 35 s):

```bash
python3 extract_bench_v12.py --fasta data/reference/chr20_full.fa \
  --bam data/giab_hg002_chr20/chr20_15x.bam --vcf data/giab_hg002_chr20/chr20.vcf.gz \
  --bed data/giab_hg002_chr20/chr20_highconf.bed --contig chr20 \
  --region 33000000 33200000 --chunk-bp 500000 --features none --out cache/smoke.npz
mkdir -p out && PIPE_OUT=out PIPE_CACHE=cache/smoke.npz PIPE_TAG=smoke \
  python3 experimental/chr20_validation/run_chr20_pipeline.py
```

Expected: `Wrote cache/smoke.npz | 193792 loci | 203 SNP`, then `out/ai_cascade_smoke.vcf` and `out/ai_pb_only_smoke.vcf`.

Not re-run for this README: the full-chromosome chains above and hap.py scoring. `verify.sh` covers the build, the tests and a call-level reproduction on one region.

---

## Results

All numbers are SNP-only, GRCh38, GIAB v4.2.1 high-confidence regions. Full tables with evaluators and intervals are in the manuscript.

**V2.x A+C by evidence tier** (hap.py 0.3.15):

| Tier | Data | Precision | Recall | F1 |
|---|---|---|---|---|
| Development (in-sample) | HG002 chr20, even 2-Mb blocks | 0.998662 | 0.980511 | 0.989503 |
| Holdout | HG002 chr20, odd 2-Mb blocks | 0.996362 | 0.981781 | 0.989017 |
| Final evaluation (single pre-registered run) | HG003 chr8, whole chromosome | 0.997567 | 0.980603 | 0.989012 |

HG003 is the father of HG002 and was sequenced on the same platform, aligner and read length. The final evaluation therefore separates sample and chromosome from the development data, but not platform, ancestry or technology. No V2 control was run on HG003, so it does not show that the gain over V2 transfers.

**External callers on HG002 chr20 at 300×** (descriptive, not a ranking; runtimes not matched): DeepVariant 1.10.0 F1 0.997613, GATK 4.6.2.0 0.991712, V2.x A+C 0.989271 (whole chromosome, includes the fitting blocks), V1/V2 0.972429. DeepVariant and GATK also call indels.

**V1 cascade vs PB-only** (paired ΔF1): see the manuscript, Table 4. The whole-chromosome chr20 run at 15× was DEGRADED by ΔF1 −5.6×10⁻⁵ (10 discordant loci in 56.2 M loci).

**Speed.** V2 vs the V1 record on HG002 chr20 300×: about 17.5× on 1 thread and 102× on 8 threads for the same work (PB at every locus); 155× for the cascade on 8 threads, which also does less work. V2.x A+C adds about 1 % CPU time and 2.4 MB memory over V2.

---

## Limitations

- SNP-only; indels, MNPs and structural variants are neither called nor scored.
- About 1.4–3.5 % of truth SNPs lie in 64-bp windows that are dropped because they are not wholly inside the confident BED. They can never be called, which caps recall (0.986 on chr20, 0.9846 on chr8).
- V2.x was fitted on one sample, chromosome, depth and platform; the final evaluation adds one related sample on the same platform. No cross-platform, cross-ancestry, genome-wide or clinical validation.
- The V1 constants are specific to about 15× depth; at 300× the router engages on 0.001 % of loci.
- V1 loses a few true SNPs at low VAF in segmental duplications; V2.x base-quality rules are specific to this data set's quality profile.
- The V1 native counts backend is not exact on adversarial BAMs with repeated read names and inconsistent mate fields (it matches pysam on all validated real data).
- HG005 raw artefacts were lost; those numbers are record-only.
- One machine, mostly single timing runs; no environment lockfile.

Full list: `research/LIMITATIONS_AND_OPEN_QUESTIONS.md` and the manuscript, Table 17.

---

## Repository layout

```text
AI_DNA_ANALYZER/
├── verify.sh                  one-command check (see above)
├── cascade.py                 V1 frozen constants and routing
├── binomial_baseline.py       V1 binomial screen
├── quality_error_model.py     V1 Poisson-binomial LLR
├── pileup_counts.py, read_level_pileup.py, extract_bench_v12.py   V1 extraction
├── native/                    V1 C/htslib backends and build.sh
├── native_v2/                 V2 engine as frozen at the V2 milestone
├── native_v2x/                V2 + V2.x engine (dnav2), models/, tests/, scripts/, analysis/
├── bench_v13 … bench_v20      V1 validation experiments
├── fetch_*.sh / prepare_*.sh  V1 data recipes
├── test_*.py                  V1 tests
├── experimental/              genotype layer (Method C), chr20 and HG005 pipelines
├── results/                   frozen V1 outputs
└── research/                  manuscript, LaTeX, figures, reports, final evaluation
```

`checkpoints*/` and `model_*.py` are from the earlier neural experiments and are not used.

## Publication and citation

- Manuscript: `research/latex/paper.pdf` (source `research/PAPER_MANUSCRIPT_V2.md`)
- Software, v2.0.0 (V1, V2, V2.x, final evaluation): DOI [10.5281/zenodo.22978775](https://doi.org/10.5281/zenodo.22978775)
- Historical V1 release, v1.0.0: DOI [10.5281/zenodo.22895746](https://doi.org/10.5281/zenodo.22895746)
- Claim–evidence map: `research/CLAIM_EVIDENCE_MAP.csv`; AI-assistance disclosure: `research/AI_ASSISTANCE_DISCLOSURE.md`
- Citation metadata: `CITATION.cff`

## License

GPL-3.0-only. See `LICENSE`.
