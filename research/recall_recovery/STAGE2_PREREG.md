# Varion Stage 2 — recall recovery: pre-specification

Written 2026-10-26, BEFORE any Varion run on the HG004 confirmatory data. Separate analysis; it does not alter or restate any Stage 1 result.

Naming: **Varion V1** = original implementation; **Varion R1** = Stage 1 depth-robust classifier (`R1.model`); **Varion R2** = reserved for a future implementation. Stage 1 `RESULTS.md` uses the label "R2" for R1 + log1p(depth); in this stage that historical model is called **S1-R2** and is not used. Historical files keep their original names.

## Hypotheses
* **H1 (window tiling).** `build_windows` (`native_v2x/src/frame.cpp`) keeps only 64-bp windows lying fully inside one merged BED interval, so BED positions in edge windows can never be called. Prediction: on HG003 chr8 (exploratory), 2,873 of R1's 4,212 FN (30x) and of 3,909 FN (300x) are such positions, with 0 TP there (descriptive audit, session 2026-10-26, `fn_breakdown.py`).
* **H2 (truth-derived call mask).** `engine.cpp` (`score_locus`) sets `frame[i] = truth_label <= SNP`; truth-indel positions are never called. Prediction: removal costs no SNP recall and raises FP (precision loss). Not a recall hypothesis; a validity requirement.

## Variants (all use R1.model sha256 63816707…b711ed, frozen constants 7.0 / router / PB 10.5 unchanged; no retraining, no threshold change in this stage)
| Label | Flags added to the R1 command |
|---|---|
| B0 R1-as-published | none |
| H1 | `--edge-windows` |
| H2 | `--truth-free` |
| H1+H2 | `--edge-windows --truth-free` (run only after H1 and H2 are each understood; interaction reported) |

Reference for the honest claim is **R1 with `--truth-free`** (call set without truth masking); B0 is reported for continuity with Stage 1.

## Data and exposure ledger
| Set | Status |
|---|---|
| HG002 chr20 even blocks | train/dev in Stage 1 (exposed) |
| HG002 chr20 odd 30x (seed 7), HG003 chr8 30x/300x | exposed (Stage 1 test + H1 audit); **exploratory only** |
| HG004 chr21 (new sample, new chromosome) | untouched; **confirmatory C1** |
| HG004 chr20 (new sample, chromosome used in training for HG002) | untouched sample; **confirmatory C2**; chromosome exposed, so C2 is weaker than C1 |
Depths: HG004 2x250 is ~56x genome-wide. Full contig (~56x), `samtools view -s 21.53` (~30x), `-s 21.10` (~6x); seed 21 fixed now. Truth/BED: GIAB v4.2.1 HG004 `_benchmark.vcf.gz`, `_benchmark_noinconsistent.bed` (sha256 in `/mnt/archive/AI_DNA_ANALYZER_stage2/truth/SHA256SUMS.txt`). Scorer: hap.py 0.3.15, same image/options as Stage 1, SNP PASS row; same BED and `-l` for every method.
No threshold, feature, model or flag is changed after a C1/C2 number has been seen. Development work (gates, unit tests, timing) uses exploratory sets only.

## Endpoints (fixed now)
Primary (C1 and C2 each at ~30x): **H1+H2 vs R1 truth-free**.
* P1: recall gain ≥ +0.004 absolute AND F1 gain ≥ +0.004, paired 2-Mb-block bootstrap (4,000 resamples) lower 95 % bound > 0 for both.
* P2: precision loss ≤ 0.002 relative to R1 truth-free (no recall-for-precision trade-off disguised as F1 gain).
Secondary:
* S1 (H1 alone): fraction of truth SNP in BED lying in unscanned windows ≤ 0.1 % of truth SNP; precision change ≤ 0.001.
* S2 (H2 alone): report ΔFP, ΔP, ΔF1 vs B0 (expected negative; descriptive).
* S3 (consistency): C1, C2 at ~56x and ~6x reported; improvement direction must agree in all conditions to be called depth-robust.
* S4 (research target, no pass/fail claim until shown): F1 ≥ 0.990 on C1 and C2 at ~30x. bcftools QUAL≥30 (`research/depth30_baseline` command) run on the same BAMs as comparator.
Not claimed: any F1 ≥ 0.990 statement on exploratory sets counts as exploratory.

## Performance protocol and guardrails (vs B0 on the same BAM)
Same machine, `--threads 8`, no other heavy process, 5 repetitions per variant interleaved, median and range reported; wall, user+sys CPU, peak RSS from `/usr/bin/time -v`. Guardrails: wall ≤ +10 %, CPU ≤ +10 %, RSS ≤ +15 %. A breach stops acceptance and goes to the owner with the F1/recall benefit.

## Gates and tests
* Default-flag binary must produce a VCF byte-identical to the archived R1 run (done 2026-10-26 on HG002 chr20 odd 30x; recorded below).
* `make test-frame` (synthetic BED edge test) must pass.
* Truth-separation regression: with `--truth-free`, the call VCF must be byte-identical for the real truth file and a decoy truth file (to be added to `native_v2x/scripts/`).
* VCF validated with `bcftools view` and by hap.py parsing; deterministic across two repeats.

## Source identity
HEAD 091cd37 plus uncommitted `native_v2x` changes (diff sha256 of `git diff -- native_v2x` recorded per run); binary `stage2/build/dnav2` sha256 88bf8ce3cc3ce025670bf571c8f9cf0d86778758563732b6d1c86d7fa7758cd8 (built with the H1/H2 flags, defaults unchanged); archived R1 binary `v2x/build_new/dnav2` sha256 9852d184…. Gate result: stage2 binary with default flags and R1.model gave a VCF identical to the archived `t2_chr20_30x_R1` (`cmp`).

## Known limitations and confounders
* Subsamples of one library are not independent low-coverage libraries; one platform, aligner and benchmark programme; SNPs only; indels not called.
* HG004 is a relative of HG002/HG003 (Ashkenazi trio); population overlap can inflate performance relative to unrelated samples.
* Edge windows change the evaluated call domain only by calling BED positions that were previously unscanned; hap.py evaluation regions are identical for all methods.
* Confidence intervals from block bootstrap treat blocks as exchangeable; they are not independent biological replicates.
* Whole-genome use without a BED (or with a very coarse BED) is not tested here.
* If H1+H2 does not meet P1/P2 it is reported as failed; criteria are not adjusted.
