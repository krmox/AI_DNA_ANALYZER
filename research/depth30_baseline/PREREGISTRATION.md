# Supplementary extension: 30x downsample and simple-baseline comparison — pre-specification

Written 2026-10-25, BEFORE any output of these runs exists. Requested by the author after the closed final evaluation (research/final_independent_eval); it is a NEW, separately labelled analysis and does not alter any earlier result.

## Questions
1. Does the frozen V2.x A+C model keep its gain over the V2 control at a typical depth (~30x)?
2. Does a standard off-the-shelf simple caller (bcftools) with a fixed filter match V2.x on the same data?

## Data (fixed)
* HG002 chr20 GRCh38, the 300x BAM used throughout (`HG002.GRCh38.300x_chr20.bam`), truth GIAB v4.2.1, scorer hap.py 0.3.15 (xcmp, SNP row, PASS), same image/options as before.
* 30x set: `samtools view -s 42.1` (seed 42, fraction 0.1) of the 300x BAM, once. No second seed. Achieved depth is measured and reported.
* Scoring: (a) whole chr20 (mixed dev + holdout; engineering) and (b) holdout only = odd 2-Mb blocks (block = floor(pos0/2,000,000), odd = holdout), as in SPLIT_PREREGISTRATION.md, via hap.py `-T` on a BED of odd blocks.

## Methods compared (all frozen before running)
* V2 control: `dnav2 run ... --mode cascade` (no model). Constants unchanged.
* V2.x A+C: same command with `--model A_C.model` (sha256 e9942a22…2d51ce). Model, threshold and constants NOT refitted for 30x.
* Simple baseline (bcftools, SNPs only): `bcftools mpileup -f REF -r chr20 -q 20 -Q 13 -d 10000 -a FORMAT/AD` then `bcftools call -mv`, keep SNPs, keep QUAL >= 30. MAPQ/baseQ admission (20/13) matches V1. QUAL >= 30 is a conventional value chosen now and not tuned. The same baseline is also run once on the 300x BAM for reference.

## Rules
* One run per method per depth. No tuning, no repeat for a better number, no change of threshold after seeing results. A failure is reported as is.
* All outcomes are reported whether favourable or not (including if V2.x degrades at 30x, or the baseline is equal or better).
* Runtime is descriptive only (single runs, warm/cold cache not controlled).
* Limits stated in advance: one individual, one chromosome, the 30x set is a subsample of the same sequencing run (not an independent 30x library), the model was fitted on 300x data (even blocks), so dev blocks in whole-chr20 numbers are not independent.
