# Genotype Layer -- Dev Region Results (chr21:30.0-30.47Mb)

Dev region only (data/giab_hg002_real_30M). Held-out region (data/giab_hg002_chr21_12Mb) NOT touched, per instructions.

## Table 1: Method x Genome[dev] x Precision/Recall/F1/GTacc/0/1acc/1/1acc/Runtime

| Method | Precision | Recall | F1 | GT Acc (overall) | 0/0 Acc | 0/1 Acc | 1/1 Acc | Genotype-layer runtime |
|---|---|---|---|---|---|---|---|---|
| A (baseline, always 0/1) | 0.2883 | 0.2893 | 0.2888 | 0.2914 | nan | 1.0000 | 0.0000 | n/a (no genotype layer) |
| B (VAF-threshold) | 0.9892 | 0.9928 | 0.9910 | 1.0000 | nan | 1.0000 | 1.0000 | 3.37 ms |
| C (binomial posterior) | 0.9892 | 0.9928 | 0.9910 | 1.0000 | nan | 1.0000 | 1.0000 | 172.84 ms |

## Table 2: Method x AlleleF1 x GA4GH-F1 x delta x extra runtime

| Method | Allele-level F1 (pos+ALT only) | GA4GH/hap.py F1 | Delta (GA4GH - Allele) | Extra genotype-layer runtime vs. A |
|---|---|---|---|---|
| A | 0.9037 | 0.2888 | -0.6149 | 0 ms (baseline) |
| B | 0.9037 | 0.9910 | +0.0873 | 3.37 ms |
| C | 0.9037 | 0.9910 | +0.0873 | 172.84 ms |

Allele-level TP=549 FP=6 FN=111 (identical across A/B/C by construction -- same 555-call set).

## GT confusion matrices (rows=truth, cols=predicted; order 0/0,0/1,1/1)

**Method A** (n=549 loci with known truth GT):

| truth\pred | 0/0 | 0/1 | 1/1 |
|---|---|---|---|
| 0/0 | 0 | 0 | 0 |
| 0/1 | 0 | 160 | 0 |
| 1/1 | 0 | 389 | 0 |

**Method B** (n=549 loci with known truth GT):

| truth\pred | 0/0 | 0/1 | 1/1 |
|---|---|---|---|
| 0/0 | 0 | 0 | 0 |
| 0/1 | 0 | 160 | 0 |
| 1/1 | 0 | 0 | 389 |

**Method C** (n=549 loci with known truth GT):

| truth\pred | 0/0 | 0/1 | 1/1 |
|---|---|---|---|
| 0/0 | 0 | 0 | 0 |
| 0/1 | 0 | 160 | 0 |
| 1/1 | 0 | 0 | 389 |

## Fitted parameters

- Method B: fitted tau = **0.69** (swept 0.05-0.95 step 0.01 on dev region; GT accuracy at fit = 1.0000)
- Method C: eps = **0.01** (fixed constant, not fit to labels)

## High-uncertainty / flagged sites

- Method C: theta=0 (0/0) won the posterior at **3/555** called sites. These are NOT suppressed -- they are kept as GT 0/1 with low GQ (<=5), per the safety requirement that a genotype layer must never drop a detection.

## Two-level invariant (allele call set unchanged across methods)

- Allele call set (CHROM:POS:REF:ALT) identical across A/B/C: **True**
- N called loci: 555 (0 skipped for ref=N/no-alt-support)
- Allele-level P/R/F1 identical by construction across A/B/C: P=0.9892 R=0.8318 F1=0.9037

