# Varion Stage 2 — recall recovery: results

Pre-specification: `STAGE2_PREREG.md` (frozen 2026-10-26 before any run on HG004). Raw data: `METRICS_STAGE2.csv`, `TIMING_RAW.tsv`, `BOOTSTRAP_STAGE2.txt`; commands: `scripts/`. Artifacts: `/mnt/archive/AI_DNA_ANALYZER_stage2/` (runs, happy/out, timing, bam, truth). Status: **implementation complete, locally tested, benchmarked on HG004 chr20/chr21; not independently validated; not production-ready** (flags off by default, uncommitted).

Labels: Varion R1 = Stage 1 `R1.model` (sha256 63816707…b711ed). "Stage-1 R2" (R1 + log depth) is not used. B0 = R1 as published (truth-masked call set). H1 = `--edge-windows`. H2 = `--truth-free`. Reference for the honest claim = H2 (R1 without truth-derived mask).

## 1. Observations

### 1.1 Code audit
* `frame.cpp` `build_windows` kept only 64-bp tiles fully inside one merged BED interval. BED positions in edge tiles were never scanned.
* `engine.cpp` `score_locus` set `frame[i] = truth_label <= SNP`: positions where the truth VCF has an indel were never called (truth-derived call mask).
* Truth VCF otherwise used only for the `label` column in feature/locus dumps.

### 1.2 Exploratory (already-exposed data; HG003 chr8 30x, seed 7)
FN breakdown of R1 (4,212 FN): 2,873 (68 %) in edge tiles, 0 TP there; 1,339 in scanned tiles. At 300x: 2,873 of 3,909 (73 %). Truth SNPs at indel-labelled positions: 0 FN.
With H1+H2: 1,543 FN remain. Of these (descriptive): 560 have zero ALT reads after filters (176+ with zero raw depth; of 216 with raw depth ≥10, 179 have no ALT read even unfiltered); 418 have 1–2 ALT reads (median depth 10–16); 422 have ≥3 ALT reads but are not called (338 of them would pass the binomial 7.0 test); 153 have a call with the correct ALT but the wrong genotype (het→1/1: 71, hom→0/1: 80), counted by hap.py as FN **and** FP (all 153 are also FP positions; 153 of 650 FP).

### 1.3 Confirmatory accuracy (HG004, SNP PASS, hap.py 0.3.15, 2-Mb block bootstrap, 4,000 resamples, seed 20261026)
| Set | Variant | TP | FN | FP | Precision | Recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| chr21 30x | B0 | 44538 | 1235 | 336 | 0.9925 | 0.9730 | 0.98267 |
| chr21 30x | H2 (honest baseline) | 44538 | 1235 | 360 | 0.9920 | 0.9730 | 0.98241 |
| chr21 30x | H1 | 45314 | 459 | 376 | 0.9918 | 0.9900 | 0.99087 |
| chr21 30x | H1+H2 | 45314 | 459 | 400 | 0.9913 | 0.9900 | 0.99061 |
| chr21 30x | bcftools QUAL≥30 | 45250 | 523 | 270 | 0.9941 | 0.9886 | 0.99131 |
| chr20 30x | B0 | 70012 | 1647 | 570 | 0.9919 | 0.9770 | 0.98441 |
| chr20 30x | H2 | 70012 | 1647 | 604 | 0.9914 | 0.9770 | 0.98418 |
| chr20 30x | H1 | 71105 | 554 | 594 | 0.9917 | 0.9923 | 0.99199 |
| chr20 30x | H1+H2 | 71105 | 554 | 630 | 0.9912 | 0.9923 | 0.99174 |
| chr20 30x | bcftools QUAL≥30 | 71068 | 591 | 566 | 0.9921 | 0.9918 | 0.99193 |

Other depths (full table in `METRICS_STAGE2.csv`): H1+H2 F1 at ~56x 0.99147 (chr21) / 0.99235 (chr20); at ~6x 0.85084 / 0.84099 (precision 0.947 / 0.941).

Paired bootstrap, H1+H2 minus H2, ~30x (`BOOTSTRAP_STAGE2.txt`):
| Set | dF1 [95 % CI] | dRecall [95 % CI] | dPrecision [95 % CI] |
|---|---|---|---|
| chr21 | +0.00820 [+0.00671, +0.00994] | +0.01695 [+0.01421, +0.02068] | −0.00073 [−0.00200, −0.00013] |
| chr20 | +0.00756 [+0.00681, +0.00838] | +0.01525 [+0.01373, +0.01690] | −0.00023 [−0.00052, −0.00004] |
Same direction at ~56x and ~6x on both chromosomes (dF1 +0.0054…+0.0085).

### 1.4 Performance (HG004 ~30x, `--threads 8`, 5 interleaved repetitions, median; whole run incl. VCF writing)
| Contig | Variant | Wall s (range) | CPU s (user+sys) | Peak RSS MB | dWall | dCPU | dRSS |
|---|---|---|---|---|---:|---:|---:|
| chr21 | B0 | 13.29 (12.63–14.37) | 63.7 | 181.0 | – | – | – |
| chr21 | H1 | 12.94 | 63.0 | 193.9 | −2.6 % | −1.0 % | +7.1 % |
| chr21 | H1+H2 | 12.73 | 62.8 | 185.0 | −4.2 % | −1.4 % | +2.2 % |
| chr20 | B0 | 22.71 (22.29–27.10) | 113.7 | 209.5 | – | – | – |
| chr20 | H1 | 22.52 | 112.3 | 225.8 | −0.8 % | −1.2 % | +7.8 % |
| chr20 | H1+H2 | 22.34 | 111.8 | 216.2 | −1.6 % | −1.7 % | +3.2 % |
bcftools (1 thread): chr21 405 s wall / 412 s CPU; chr20 658 s / 682 s. Output VCFs identical across the 5 repetitions for every variant.

### 1.5 Tests and gates
* Default-flag binary reproduces the archived R1 VCF byte-for-byte (HG002 chr20 odd 30x, `cmp`).
* `make test-frame` passed; `scripts/test_truth_free.sh` passed (VCF independent of the truth file with `--truth-free`; 3,926 records vs 3,922 truth-masked).
* `make test` (Python kernel tests) passed on the R1.1 build; not run: sanitizers, VCF validator beyond hap.py/bcftools parsing.

## 2. Pre-registered criteria
| Criterion | chr21 30x | chr20 30x |
|---|---|---|
| P1 recall ≥ +0.004 and F1 ≥ +0.004, CI lower bound > 0 | met (+0.0170 / +0.0082) | met (+0.0153 / +0.0076) |
| P2 precision loss ≤ 0.002 | met (−0.0007) | met (−0.0002) |
| Guardrails wall/CPU ≤ +10 %, RSS ≤ +15 % | met | met |
| S3 same direction at all depths | met | met |
| S4 F1 ≥ 0.990 at ~30x (research target) | 0.99061 (point estimate) | 0.99174 (point estimate) |
S1 (H1): truth SNPs in unscanned tiles — not re-measured on HG004 (FN fell from 1,235→459 and 1,647→554, consistent with H1); S2 (H2): ΔFP +24 / +34, ΔF1 −0.0003 / −0.0002 (30x).

## 3. Interpretation (supported by the data, not proven beyond it)
* The tiling rule, not the classifier, explained most of the recall deficit of R1 under the Stage 1 protocol: fixing it recovers ~1.5–1.7 points of recall at 30x on two HG004 contigs and ~2.9 k extra TP on HG003 chr8.
* The truth-derived mask inflated precision slightly; removing it costs ~0.0002–0.0007 precision at 30x.
* With both changes R1 is at parity with bcftools QUAL≥30 on F1 at 30x on HG004 (−0.0007 on chr21, −0.0002 on chr20), with higher recall and lower precision, at ~30x less wall time and ~6.5x less CPU time (bcftools single-threaded).
* The remaining ~1.5 k FN on HG003 chr8 are mostly data-limited (no or 1–2 ALT reads); the recoverable parts (genotype errors, model-rejected k≥3 candidates) are small.

## 4. Unverified / limitations
* F1 ≥ 0.990 is shown only as point estimates on two contigs of one sample (HG004), subsamples of one library; chr20 chromosome was used in R1 training (other sample). No CI for absolute F1.
* "Independent validation" is not achieved: HG004 is related to HG002/HG003; one platform and aligner.
* Timings taken on a desktop with background load (load average ~3); wall-time differences of a few percent are within the observed run-to-run spread (up to ~20 %). Single machine, 8 threads.
* ~6x depth is far out of distribution for R1 (trained on ≥30x): F1 ≈ 0.85, not a depth-robustness success.
* H1 clips calling to BED positions; behaviour without a BED or with a coarse BED is untested. SNPs only.
* The Stage 1 conclusion (primary criterion not met) is unchanged.

## 5. Decisions
| Change | Status |
|---|---|
| H1 `--edge-windows` | Supported on confirmatory data; kept as opt-in; making it the default is a pending owner decision (changes the evaluated call domain) |
| H2 `--truth-free` | Required for a valid claim; kept as opt-in; making it the default is a pending owner decision |
| Genotype layer, model threshold/retraining | Not implemented; need fresh data (HG004 chr20/chr21 now exposed); pending owner decision |
| Public renaming to Varion (README, CLI, targets) | Not done; no public identifiers changed; needs owner decision (compatibility) |
No commits, tags or releases were made. Historical documents (`bioadv_draft`, Stage 1 files) are unmodified.

## 6. Varion R1.1 defaults (owner decision B, 2026-10-26)
`--edge-windows` and `--truth-free` are now the defaults; `--strict-windows` and `--truth-mask` restore R1. Build `stage2/build_r11/dnav2`, sha256 below; the VCF header gains `##varion_mode=edge_windows=N,truth_free=N` (the only header change; `##source` unchanged). `--truth` is optional in the default mode; if given together with `--dump-features/--dump-loci/--setup-only` it only fills the `label` column and never reaches a call.
Gates (all run): (1) `--strict-windows --truth-mask` VCF body identical to the archived R1 run (HG002 chr20 odd 30x); (2) default mode, run **without** `--truth`, body identical to the H1+H2 run on HG004 chr21 30x; (3) dump label column populated (85 SNP / 217 normal in a 3,000-window slice); (4) `make test-frame`, `scripts/test_truth_free.sh` (default VCF independent of truth file; R1 mode differs: 3,922 vs 3,911 records on the slice), `make test` pass.
Not re-measured for R1.1: timing (the code path equals H1+H2 measured above) and accuracy beyond gate (2).
R1.1 binary sha256 e83d87fd2f810eb7a9d096fb797c6d3451aaef86d14a0071dc2ad9b91a3a39df; the Stage 2 measurements in sections 1.3–1.4 were taken with the flag build (`stage2/build/dnav2`, 88bf8ce3…758cd8), functionally identical to the R1.1 defaults.

### 6.1 R1.1 final timing (HG004 ~30x, `--threads 8`, 5 interleaved reps, median; raw: `TIMING_R11_RAW.tsv`)
| Contig | Variant | Wall s (range) | CPU s | Peak RSS MB | dWall | dCPU | dRSS (vs archived R1 binary) |
|---|---|---|---|---|---:|---:|---:|
| chr21 | archived R1 binary | 13.73 (12.77–14.25) | 63.1 | 181.3 | – | – | – |
| chr21 | R1.1 binary, R1 flags | 13.01 | 63.6 | 181.0 | −5.2 % | +0.7 % | −0.2 % |
| chr21 | R1.1 defaults | 13.33 (12.56–13.96) | 62.7 | 187.6 | −2.9 % | −0.7 % | +3.5 % |
| chr20 | archived R1 binary | 23.01 (21.99–24.16) | 113.4 | 210.1 | – | – | – |
| chr20 | R1.1 binary, R1 flags | 23.33 | 113.6 | 210.0 | +1.4 % | +0.1 % | −0.0 % |
| chr20 | R1.1 defaults | 22.80 (21.53–23.09) | 111.2 | 215.4 | −0.9 % | −1.9 % | +2.5 % |
All three guardrails met (wall ≤ +10 %, CPU ≤ +10 %, RSS ≤ +15 %). Body of the R1-flag VCF equals the archived-binary VCF in every repetition; defaults VCF identical across repetitions. Background load during this run was higher than in 1.4 (unrelated process ~1.4 cores, load average ≈ 6): differences of a few percent in wall time are within noise.
