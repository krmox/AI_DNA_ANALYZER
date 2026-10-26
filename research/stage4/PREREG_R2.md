# Varion R2 — preregistration (written 2026-10-26 BEFORE any R2 model fit or R2 run on blind data)

R1.1 stays frozen (commit 659df82, dnav2 e83d87fd, R1.model 63816707). R2 lives in `native_r2/`; the new "tail" lives in `native_r2/tail/`. No commits.

## Data roles (exposure ledger)
| Role | Sets | Use |
|---|---|---|
| TRAIN | HG002 chr20 (30x/60x/100x), HG003 chr8 (30x), HG004 chr21 (6,10,15,20,30,56x), HG004 chr20 (6,30,56x) | fit models; blocked CV by 2-Mb block inside TRAIN |
| VALID (already inspected in Stage 3 diagnostics) | HG004 chr22 (10,15,20,30,68x) | choose arm, threshold, tail strength; may be looked at |
| BLIND | chr18 of HG002/HG003/HG004 (10–native, 3 seeds) and HG005 chr22/chr18 (10–100x) | run ONCE after R2 is frozen; no error inspection before freeze |
Stage 3 results for BLIND conditions already exist for R1.1 and bcftools (aggregate counts only; no per-site error analysis was done on them).

## Arms (offline ablation on TRAIN via blocked CV, then VALID)
* A0 R1 model as is (control = R1.1 behaviour).
* A0T NO retraining: frozen R1.model score z plus a tail adjustment with <= 3 parameters (depth_rel penalty, weak-VAF penalty, candidate-cluster penalty) fitted on TRAIN only. Added at the owner's request (2026-10-26) as the lightest alternative; compared with the retrained arms by the same selection rule. Model size/runtime: retraining changes only the coefficient file, not the runtime cost.
* A1 same 29 features, retrained on TRAIN without the truth-indel mask (truth-free frame).
* A2 A1 + tail block 1: depth_rel (= depth / contig mean depth from BAM index statistics), nc25, lmr25, nearest, depth_ratio25, neigh_max_vaf25, hp_len.
* A3 A2 + tail block 2 (pattern layer): candidate/low-VAF density at three scales (±100 bp, ±5 kb inside the task, contig mean) as ratios.
* A4 A3 + genotype layer (allele-balance aware).
* B  best feature set with a small tree ensemble (<= 200 trees, depth <= 4) instead of logistic.
Selection rule (fixed now): maximise the mean over VALID depths of F1; ties -> simpler arm. Decision threshold chosen by blocked CV on TRAIN only. No other search.

## Endpoints (BLIND, R2 vs R1.1 and vs bcftools baseline of Stage 3, same BAMs, hap.py SNP PASS, 2-Mb block bootstrap)
* P1 (accuracy gain): at 30x, mean over contigs/seeds F1(R2) − F1(R1.1) ≥ +0.003, CI lower bound > 0.
* P2 (non-inferiority to bcftools): at every depth 10, 15, 20, 30, 50x, paired ΔF1 (R2 − bcftools) lower 95 % bound ≥ −0.0005.
* P3 (target): F1(R2) ≥ 0.990 at ≥30x is reported per contig/sample; claimed only where the bcftools baseline also reaches it (ceiling check), otherwise reported as data-limited.
* P4 (no regression at low depth): F1(R2) ≥ F1(R1.1) − 0.002 at 10x and 15x.
* Guardrails vs R1.1 (same BAM, `--threads 8`, 5 interleaved reps): wall ≤ +10 % or ≤ +2 s (whichever larger), CPU ≤ +10 %, peak RSS ≤ +15 %.
A failure of any endpoint is reported as is; nothing is re-tuned after a BLIND number is seen.

## Engineering gates (before the BLIND run)
* With tail disabled and R1.model, R2 output is byte-identical to R1.1 (VCF body).
* Output independent of thread count and of truth file; deterministic across repeats.
* Builds clean with `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wformat=2`; unit tests for tail; ASan+UBSan and TSan runs on a slice; valgrind/memcheck on a small slice; fixed-seed fuzz on tail inputs (empty windows, depth 0, N).
* Optimisation: profile first; ASM/AVX2 only for a measured hot loop with a scalar reference and exact-equivalence test.

## Limits stated in advance
Same lab/platform/aligner family as training; related trio; HG005 shares platform; SNP-only; the 0.99 target may be unreachable on contigs where the data cannot support it (see Stage 3 diagnostics: ~0.989 best-of-both ceiling on HG004 chr22 30x).

## Addendum 2026-10-26 (before any R2 fit): random new genome and added R2 components
* Owner asked for a blind test on a "random new genome". Candidate pool: GIAB samples with GRCh38 v4.2.1 truth and a GRCh38 Illumina BAM that no project stage has used: HG006 and HG007 (Chinese parents, HiSeq 2x148 novoalign 100x; different read length from all training data). Draw with `random.Random(20261027)`: sample **HG006**; contig drawn from {chr16, chr17, chr19} -> **chr16**. (chr16 had small V1-era windows on HG002 only, different model; HG006 itself is untouched.)
* BLIND-NEW = HG006 chr16: native (~100x, `HG006.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam`), 50x, 30x, 20x, 15x, 10x with seeds 101-103 (50x/native: one). Not downloaded or inspected before R2 freeze.
* Added R2 components (owner's decision): (1) per-sample self-calibration pre-pass (depth, MAPQ and BQ distributions, error rate); (2) error rate estimated from data instead of fixed 0.01 (mode chosen on VALID; clamp [0.002, 0.03]); (3) OOD guard: if calibration statistics fall outside the training envelope, the VCF header records `##varion_ood` and calling falls back to the conservative Binomial caller with the estimated error rate; (4) the owner's "tail": difficult-locus detector + three-scale pattern layer with researcher parameters (`--tail-regions`, `--tail-strength`).
* Domain-shift test is therefore two-fold: read length (148 vs 250), run/library (different lab protocol), and an unrelated genome.

## Addendum 2 (2026-10-26, after offline ablation, BEFORE any real-binary VALID run and before any BLIND run)
* Offline ablation on TRAIN/VALID (candidate-level approximation; `research/stage4/results/arms_results.csv`): valid mean F1 over 5 depths: A0 0.95752, A0T 0.95889, A1 0.96006, A2 0.96051, A2R 0.96038, A3 0.96059, A3R 0.96052. At 30x / 68x: A0 0.98612 / 0.98816; A3R 0.99016 / 0.99261.
* Selection: the literal rule (max mean) names A3 (0.96059); A3R differs by 0.00007, far below the resolution of a one-contig candidate-level approximation. **R2 = A3R** (A3 with the absolute MAPQ/BQ features replaced by sample-relative ones) because the owner asked for universality and the robust variant costs nothing within noise. A3 is run on the blind HG006 30x condition as a sensitivity check only.
* OOD envelope fixed from the physics of the model, not from the data: eps_hat in [0.0005, 0.02], read MAPQ p90 in [40, 80], BQ median in [15, 60], median depth in [4, 800]. (A first draft used 0.5x/3x of the training range; rejected before use because it would route any other Illumina run to the fallback.)
* The pattern layer's marginal value over A2 is ~+0.0001 F1 at 30x; it stays in R2 because it is cheap and is the owner's design, but no claim of benefit is made for it.
* Clarification (2026-10-26): downloading the BLIND-NEW data before the R2 freeze is allowed (it cannot leak model information); running any caller on it or inspecting errors before the freeze is not. HG006 v4.2.1 truth/BED downloaded (md5 checked against the release md5.in).
* Reduced BLIND-NEW schedule (owner needs results quickly, 2026-10-26): HG006 chr16 conditions run first = 30x seeds 101-103, native (~100x), 10x/15x/20x/50x seed 101 only. Seeds 102-103 for 10x/15x/20x/50x are added afterwards if time permits; until then those depths are reported with one repeat. Endpoints are unchanged. Blind runs of chr18 (HG002/3/4) and HG005 follow HG006.
