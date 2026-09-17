# Genotype Layer -- Method C (Binomial Genotype Posterior) on TRUE Hold-Out Validation (chr21:32-44Mb)

Method C reuses the EXACT SAME 16,124 cascade calls used for A/B (cache/validation_region.npz) and applies the identical Method C logic from build_vcfs.py verbatim: eps=0.01 fixed (not fit), theta in {0.0, 0.5, 1-eps} representing {0/0, 0/1, 1/1}, GT=argmax likelihood, GQ=10*log10(P_best/P_second) capped at 99. Sites where theta=0 (0/0) wins are NEVER dropped -- kept as GT 0/1 with low GQ and flagged (never-suppress-a-detection invariant).

A and B numbers below are pulled verbatim from RESULTS_VALIDATION.md (not recomputed here).

## Table 1: A/B/C x GT Acc x 0/1 Acc x 1/1 Acc x Allele F1 x hap.py F1 x Runtime

| Method | GT Acc (overall) | 0/1 Acc | 1/1 Acc | Allele-level F1 | hap.py SNP F1 | Genotype-layer runtime |
|---|---|---|---|---|---|---|
| A (baseline, always 0/1) | 0.6464 | 1.0000 | 0.0000 | 0.8644 | 0.6232 | n/a |
| B (VAF-threshold, tau=0.69 frozen) | 0.9250 | 0.8840 | 0.9998 | 0.8644 | 0.8921 | 15.68 ms |
| C (binomial posterior, eps=0.01 fixed) | 0.9888 | 0.9836 | 0.9984 | 0.8644 | 0.9538 | 7289.23 ms |

## Error-type breakdown (truth->pred mismatches, B vs C)

| Error type | B | C |
|---|---|---|
| 0/1 -> 1/1 | 1195 | 169 |
| 1/1 -> 0/1 | 1 | 9 |
| other | 0 | 0 |
| **total errors** | 1196 | 178 |

## GT confusion matrix, Method C (rows=truth, cols=predicted; order 0/0,0/1,1/1)

n=15940 loci with known truth GT:

| truth\pred | 0/0 | 0/1 | 1/1 |
|---|---|---|---|
| 0/0 | 0 | 0 | 0 |
| 0/1 | 0 | 10134 | 169 |
| 1/1 | 0 | 9 | 5628 |

## GT confusion matrix, Method B (for direct comparison; from records_validation.npz)

n=15940 loci with known truth GT:

| truth\pred | 0/0 | 0/1 | 1/1 |
|---|---|---|---|
| 0/0 | 0 | 0 | 0 |
| 0/1 | 0 | 9108 | 1195 |
| 1/1 | 0 | 1 | 5636 |

## Direct comparison: what changed between B and C

- GT accuracy: B=0.9250, C=0.9888, delta=+0.0639
- 0/1 accuracy: B=0.8840, C=0.9836, delta=+0.0996
- 1/1 accuracy: B=0.9998, C=0.9984, delta=-0.0014
- hap.py SNP F1: B=0.8921, C=0.9538, delta=+0.0617
- hap.py Precision: B=0.9135, C=0.9767, delta=+0.0632
- hap.py Recall: B=0.8716, C=0.9319, delta=+0.0603
- Error count: B=1196, C=178, delta=-1018
- theta=0(0/0)-favored sites (C only): 7 of 16124 (0.043%); flagged, kept as GT 0/1, low GQ, never dropped. Of these, truth=0/1 at 3, truth=0/0 at 0 (truth=0/0 not in call set otherwise by construction -- these are candidate false detections upstream of the genotype layer).

## Stratified GT accuracy, C vs B

| Stratum | n | B GT acc | C GT acc | Delta |
|---|---|---|---|---|
| Het truth, VAF<0.3 | 777 | 1.0000 | 1.0000 | +0.0000 |
| Het truth, VAF 0.3-0.5 | 3790 | 1.0000 | 1.0000 | +0.0000 |
| Het truth, VAF 0.5-0.69 | 4541 | 1.0000 | 1.0000 | +0.0000 |
| Het truth, VAF>0.69 | 1195 | 0.0000 | 0.8586 | +0.8586 |
| Hom truth (1/1), VAF>=0.95 (near-VAF-1) | 5570 | 1.0000 | 1.0000 | +0.0000 |
| Hom truth (1/1), depth<10 (low-coverage) | 1115 | 0.9991 | 0.9928 | -0.0063 |
| Depth<10 | 2708 | 0.8785 | 0.9642 | +0.0857 |
| Depth 10-19 | 11230 | 0.9293 | 0.9928 | +0.0635 |
| Depth 20-29 | 1969 | 0.9639 | 1.0000 | +0.0361 |
| Depth 30+ | 33 | 0.9394 | 1.0000 | +0.0606 |
| Segdup/low-mappability overlap | 417 | 0.9353 | 0.9856 | +0.0504 |
| Tandem repeat overlap | 897 | 0.9130 | 0.9766 | +0.0635 |

Difficult-region strat used: data/strat_v31/lowmap_segdup.bed.gz (1906 intervals overlapping chr21:32-44Mb) and data/strat_v31/tandemrepeats.bed.gz (18314 intervals overlapping chr21:32-44Mb). No local GIAB HG004 stratification BED applicable to this chr21/HG002 region was found under data/giab_hg004_v19/ (that directory holds chr2/chr3/chr5 HG004 test data unrelated to this analysis).

## Efficiency: C vs B

- Method B total genotype-layer runtime: 15.68 ms (0.9725 us/site)
- Method C total genotype-layer runtime: 7289.23 ms (452.0731 us/site)
- Runtime ratio C/B: 464.9x
- Both measured identically: wall-clock time.perf_counter() around the pure-numpy/scipy genotyping loop over already-extracted (k, n) counts, per-site, excluding BAM-scanning/detection (a cost already paid identically for A/B/C by the shared cascade extraction).

## Invariant confirmation

- Call-set size before (B): 16124
- Call-set size after (C): 16124
- Symmetric difference B vs C: |B-C|=0, |C-B|=0 (both zero -> identical call sets)
- n_records in C build: 16124, n_skipped_no_alt: 0

