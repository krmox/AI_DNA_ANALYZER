# Genotype Layer -- TRUE Hold-Out Validation Results (chr21:32-44Mb)

This is a strict hold-out test: tau=0.69 was fit on the dev region (chr21:30.0-30.47Mb) and is used here UNCHANGED -- not refit, not tuned, regardless of the results below.

Data used: reference data/reference/chr21_full.fa, BAM data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_15x.bam, truth data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz(.tbi), high-confidence BED data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed (matching naming pattern to the dev region's *_highconf.bed -- found directly alongside the BAM/truth VCF, no substitute needed).

## Table 1: Method x Precision/Recall/F1/GTacc/0/1acc/1/1acc/Runtime

| Method | Precision | Recall | F1 | GT Acc (overall) | 0/0 Acc | 0/1 Acc | 1/1 Acc | Genotype-layer runtime |
|---|---|---|---|---|---|---|---|---|
| A (baseline, always 0/1) | 0.6382 | 0.6089 | 0.6232 | 0.6464 | nan | 1.0000 | 0.0000 | n/a (no genotype layer) |
| B (VAF-threshold, tau=0.69 frozen) | 0.9135 | 0.8716 | 0.8921 | 0.9250 | nan | 0.8840 | 0.9998 | 15.68 ms |

## Table 2: Method x AlleleF1 x GA4GH-F1 x delta x extra runtime

| Method | Allele-level F1 (pos+ALT only) | GA4GH/hap.py F1 | Delta (GA4GH - Allele) | Extra genotype-layer runtime vs. A |
|---|---|---|---|---|
| A | 0.8644 | 0.6232 | -0.2411 | 0 ms (baseline) |
| B | 0.8644 | 0.8921 | +0.0277 | 15.68 ms |

Allele-level TP=15921 FP=203 FN=4794 (identical across A/B by construction -- same 16124-call set).

## GT confusion matrices (rows=truth, cols=predicted; order 0/0,0/1,1/1)

**Method A** (n=15940 loci with known truth GT):

| truth\pred | 0/0 | 0/1 | 1/1 |
|---|---|---|---|
| 0/0 | 0 | 0 | 0 |
| 0/1 | 0 | 10303 | 0 |
| 1/1 | 0 | 5637 | 0 |

**Method B** (n=15940 loci with known truth GT):

| truth\pred | 0/0 | 0/1 | 1/1 |
|---|---|---|---|
| 0/0 | 0 | 0 | 0 |
| 0/1 | 0 | 9108 | 1195 |
| 1/1 | 0 | 1 | 5636 |

## Frozen parameters (NOT tuned on this region)

- Method B: tau = **0.69** (fit on dev region chr21:30.0-30.47Mb; used as-is here)

## Two-level invariant (allele call set unchanged between A and B)

- Allele call set (CHROM:POS:REF:ALT) identical between A and B: **True**
- Independent set comparison: |A|=16124 |B|=16124 A-B=0 B-A=0
- N called loci: 16124 (0 skipped for ref=N/no-alt-support)
- Allele-level P/R/F1 identical by construction between A/B: P=0.9874 R=0.7686 F1=0.8644

## Stratified checks on Method B's GT accuracy

**1. Low-VAF sites (VAF<0.3):** n=777, GT acc=1.0000, per-class={'0/0': nan, '0/1': np.float64(1.0), '1/1': nan}

**2. High-VAF het-truth sites (truth=0/1 but VAF>0.69):** 1195 of 10303 total het-truth sites. Method B wrongly called 1/1 at 1195 of these (1195 candidates), rate=1.0000 (by construction this rate should be ~1.0, since VAF>tau triggers 1/1 by definition -- this quantifies exactly how many true hets get mis-genotyped by the VAF rule).

**3. Homozygous-truth (1/1) accuracy:** n=5637, GT acc=0.9998, 1/1-class acc=0.9998

**4. Low-depth loci (depth<10):** n=2708, GT acc=0.8785, per-class={'0/0': nan, '0/1': np.float64(0.7940991839296924), '1/1': np.float64(0.9991031390134529)}

**5. Segdup/low-mappability overlap:** used data/strat_v31/lowmap_segdup.bed.gz (1906 intervals overlapping chr21:32-44Mb). n=417 called loci overlap segdup/low-mappability regions, GT acc=0.9353, per-class={'0/0': nan, '0/1': np.float64(0.9039145907473309), '1/1': np.float64(1.0)}

## Did tau=0.69 generalize?

- Method A GT accuracy (baseline, always 0/1): 0.6464
- Method B GT accuracy (VAF-threshold, tau=0.69 frozen): 0.9250
- Delta: +0.2786
- hap.py SNP F1: A=0.6232, B=0.8921, delta=+0.2689
- See generalization verdict and failure-mode description in the accompanying report.

