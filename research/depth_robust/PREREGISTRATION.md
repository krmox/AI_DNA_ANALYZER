# Stage 1 — depth-robust second-stage classifier: pre-specification

Written 2026-10-25, BEFORE any multi-depth feature table, model fit, or test-set run of this stage exists. Separate, new analysis; it does not alter any earlier result. Stage 2 (frame / truth-mask, M-1) gets its own pre-specification after Stage 1 is finished.

## Motivation (facts already in hand)
* Source audit: `A_C.model` standardises `log1p(depth, n, k, nref)` with 300x means; at ~28x the depth z-score is about -10 and the model calls 98.6 % of k>=3 candidates (V2.x F1 0.9662 vs V2 control 0.9774 at 28x, seed-42 downsample of HG002 chr20, odd blocks).
* Error taxonomy of the V2 control, DEV blocks (even 2-Mb blocks) of HG002 chr20 300x only: 1,326 FP among 37,914 candidate calls; 1,296 (98 %) are positions with no truth SNP; 92.5 % have VAF 0.1-0.3; median ALT base quality 19.1 (TP 36.2); 88.5 % have ALT soft-clip fraction > 0.2 (FP rate 97.8 % when > 0.5); median ALT NM 8.0 (TP 1.3); median strand_diff 0.78 (TP 0.08). All of these are means/fractions, i.e. depth-invariant in principle.
* Caveat stated in advance: the diagnosis of the 28x failure was made after seeing the chr20 odd-block 28x result, so chr20 odd blocks are NOT a clean test for this idea. Primary test data is therefore a different sample and chromosome.

## Data (fixed)
* TRAIN/DEV: HG002 chr20 GRCh38, truth v4.2.1, even 2-Mb blocks (block = floor(pos0/2,000,000)). Depths: `samtools view -s S.F` of the 300x BAM with (F, S) = (0.10, 101), (0.20, 102), (0.33, 103), and the full 300x. One set each, no re-draws.
* TEST-1 (primary, independent sample + chromosome): HG003 chr8 GRCh38, 300x BAM already on disk (`AI_DNA_ANALYZER_final_eval`), truth HG003 v4.2.1 and its BED; full 300x and a 0.10 downsample (`-s 7.10`). Whole chr8; the model never saw chr8 or HG003.
* TEST-2 (secondary, same sample): HG002 chr20 ODD blocks, 0.10 downsample with a NEW seed (`-s 7.10`), scored with hap.py `-T` on the odd-block BED as before.
* Scorer: hap.py 0.3.15 (same image/options as every earlier run), SNP PASS row. Depth is measured with `samtools coverage` and reported.

## Models (definitions frozen now)
Logistic regression, StandardScaler, on frame candidates with ALT count k>=3 (existing candidate rule); all other loci keep the V2 control decision. Frozen constants 7.0 / router / PB 10.5 are unchanged and apply wherever the model does not.
* **R1 (scale-free)**: vaf; bllr_pr = bllr / max(n,1); kcap = min(k, 10); lowbq_frac = alt_lowbq / (k + alt_lowbq); and the existing A features that are means/fractions: alt_bq_mean, ref_bq_mean, bq_diff, alt_bq_hi_frac, alt_mq_mean, ref_mq_mean, mq_diff, alt_mq_lt40_frac, alt_fwd_frac, ref_fwd_frac, strand_diff, alt_min_strand, alt_end_mean, ref_end_mean, alt_near_end_frac, alt_nm_mean, ref_nm_mean, alt_clip_frac, ref_clip_frac, alt_uniq_frac, gap_frac, ins_frac, del_next_frac. No raw depth, n, nref, k2, no `routed` flag, no context (D) features.
* **R2 (control for the mechanism)**: R1 + log1p(depth), trained on the same multi-depth data. Tests whether multi-depth training alone fixes the depth dependence.
* Fitting (DEV only, all four depths pooled): C in {0.03, 0.1, 0.3, 1, 3}; GroupKFold(8) grouped by block, so all depths of one block stay in the same fold; one probability threshold chosen on out-of-fold predictions to maximise the MEAN of the four per-depth dev F1 values; final model refit on all dev rows. No other search.
* Finalist = higher mean dev OOF F1 (tie -> R1). Both are reported on the test sets regardless; the finalist is named before any test run.

## Comparators (same BAMs, same hap.py)
V2 control; old `A_C.model` (sha256 e9942a22...2d51ce); bcftools `mpileup -q 20 -Q 13 -d 10000 -a FORMAT/AD | call -mv`, SNPs, QUAL>=30 (as in `research/depth30_baseline`).

## Endpoints and decision rules (fixed before any test run)
* Primary: TEST-1, 0.10 downsample (~30x): F1(finalist) - F1(V2 control), with a 95 % paired block-bootstrap CI (2-Mb blocks, 4,000 resamples, from hap.py per-record output). SUCCESS if the difference is >= +0.010 and the CI lower bound > 0.
* Secondary 1: TEST-1 at 300x: finalist must not lose more than 0.002 F1 against old A+C (0.989012 on this data).
* Secondary 2: TEST-2 same criterion as primary.
* Also reported without a criterion: finalist vs bcftools, precision/recall/FP/FN, wall and CPU seconds, peak RSS.
* A failure on any endpoint is reported as is. No threshold, feature, C or seed is changed after a test number has been seen. If Stage 1 fails the primary endpoint, the next step is a diagnosis, not a retune.

## Limits stated in advance
Two samples (HG002, HG003) from one data-generation programme, one platform, one aligner; downsamples are subsamples of the same libraries, not independent low-coverage libraries; SNPs only; the V2 frame truth-mask (audit D2) still applies in this stage and is treated in Stage 2.
