# nindels — controlled experiments on precision and genotyping

Diagnosis and source locations: `DIAGNOSIS.md`. Baseline = tag `nindels-baseline-18a2010`
(the default flags still reproduce it). Untrained decision rule throughout: PASS iff the called
genotype is non-reference and QUAL ≥ 20. No threshold was tuned.

**Protocol**
- **Development (exploratory, all decisions made here):** HG004 chr20, Illumina 2×250, 30× and
  67×.
- **Rule for D, fixed before the combination results:** D = the components that each improved
  F1 over the baseline at *both* 30× and 67×.
- **Independent validation, run once after freezing (commit `1dc2ed4`):** HG001 chr19, Illumina
  2×148, 30×. That is a different individual, chromosome, read length and sequencing run.
- **Scoring:** hap.py 0.3.15 (xcmp), GIAB v4.2.1 confident regions, PASS calls.
- **GT concordance** = TP / (TP + truth indels whose allele was found but whose genotype
  differs).
- **Timing:** wall time = median of 3 runs at 30× and 1 run at 67×; 8 threads. bcftools was
  timed once and is single-threaded.

## Experiments

| Exp | Flag | Component |
|---|---|---|
| A | `--genotyper multi` | joint multi-allelic genotyping on a haplotype set, posterior QUAL |
| B | `--str-model FILE` (from `estimate-str`) | repeat slippage ε(period, units) estimated from the BAM with no labels; position-specific gap-open in the PairHMM |
| C | `--fragments` | collapse unmarked duplicates; one mate per fragment (higher BQ) |
| C2 | `--fragment-pick end` | as C, but keep the mate whose alignment puts the event farther from the read end |
| D | `--genotyper multi --str-model FILE` | A + B |

## Development results — HG004 chr20 (exploratory)

| Run | F1 30× | P 30× | R 30× | FP 30× | FN 30× | GT conc. 30× | Wall / CPU s 30× | RSS MB 30× | F1 67× | P 67× | R 67× | FP 67× | FN 67× |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Baseline | 0.8743 | 0.8334 | 0.9194 | 2023 | 887 | 0.9879 | 29.2 / 161 | 502 | 0.8390 | 0.7678 | 0.9248 | 3078 | 827 |
| A | **0.9438** | 0.9412 | **0.9464** | 674 | **590** | 0.9863 | 28.9 / 185 | 549 | 0.9516 | 0.9445 | 0.9589 | 644 | 452 |
| B | 0.8871 | 0.8646 | 0.9108 | 1573 | 981 | 0.9927 | 36.3 / 237 | 498 | 0.8697 | 0.8212 | 0.9244 | 2217 | 832 |
| C | 0.8438 | 0.8794 | 0.8110 | 1227 | 2079 | 0.9825 | 18.6 / 125 | 512 | 0.8575 | 0.8227 | 0.8954 | 2124 | 1151 |
| C2 | 0.8455 | 0.8727 | 0.8199 | 1319 | 1981 | 0.9829 | 21.9 / 130 | 543 | 0.8521 | 0.8087 | 0.9003 | 2344 | 1097 |
| **D = A+B** | 0.9397 | **0.9515** | 0.9282 | **539** | 790 | **0.9913** | 39.1 / 216 | 519 | **0.9548** | **0.9504** | **0.9594** | **572** | **447** |
| A+B+C2 | 0.8697 | 0.9566 | 0.7973 | 410 | 2230 | 0.9903 | 25.5 / 138 | 538 | 0.9394 | 0.9548 | 0.9244 | 498 | 832 |
| bcftools | 0.9330 | 0.9483 | 0.9181 | 578 | 901 | 0.9747 | 660.5 / 581 | 110 | 0.9548 | 0.9627 | 0.9470 | 426 | 583 |

STR estimation (`estimate-str`), part of B/D runtime:
- 30×: 49 s wall in a first run, before the timer existed.
- 67×: 10.8 s wall / 52.6 s CPU / 241 MB (the first run had warmed the page cache).

### Development results by repeat context, 30× (R / P / F1)

| Run | Homopolymer ≥ 7 (4908) | TR ≤ 50 bp (1449) | TR 51–200 bp (653) | Not in TR/HP (3876) |
|---|---|---|---|---|
| Baseline | 0.892 / 0.782 / 0.833 | 0.846 / 0.728 / 0.782 | 0.855 / 0.748 / 0.798 | 0.991 / 0.967 / 0.979 |
| A | 0.916 / 0.914 / 0.915 | 0.934 / 0.923 / 0.928 | 0.948 / 0.943 / 0.945 | 0.989 / 0.984 / 0.987 |
| B | 0.877 / 0.839 / 0.858 | 0.839 / 0.742 / 0.787 | 0.841 / 0.757 / 0.796 | 0.991 / 0.968 / 0.979 |
| C | 0.721 / 0.850 / 0.780 | 0.762 / 0.751 / 0.756 | 0.770 / 0.779 / 0.775 | 0.949 / 0.980 / 0.965 |
| D | 0.878 / 0.931 / 0.904 | 0.930 / 0.934 / 0.932 | 0.934 / 0.948 / 0.941 | 0.989 / 0.985 / 0.987 |
| bcftools | 0.857 / 0.914 / 0.885 | 0.946 / 0.949 / 0.947 | 0.928 / 0.937 / 0.933 | 0.986 / 0.994 / 0.990 |

At 67×, D in homopolymers gives 0.942 / 0.934 / 0.938, against A's 0.940 / 0.922 / 0.931.

### Selection of D (rule fixed beforehand)

| Component alone | ΔF1 30× | ΔF1 67× | Included |
|---|---|---|---|
| A | +0.0695 | +0.1126 | yes |
| B | +0.0128 | +0.0307 | yes |
| C | −0.0305 | +0.0185 | no |
| C2 | −0.0288 | +0.0131 | no |

**D vs A on development (honest):**
- At 30×, D has lower F1 than A (0.9397 vs 0.9438). B removes 135 FPs but loses 200 TPs, mostly
  low-depth hets in long homopolymers (homopolymer recall 0.916 → 0.878).
- At 67×, D is better (0.9548 vs 0.9516).
- D was validated as primary because the rule says so. A was validated as secondary.

## Independent validation — HG001 chr19 30× (run once)

| Run | Wall (s) | CPU (s) | Peak RSS (MB) | INDEL F1 | Precision | Recall | TP | FP | FN | GT conc. |
|---|---|---|---|---|---|---|---|---|---|---|
| Baseline | 20.4 | 109.3 | 484 | 0.8918 | 0.8560 | 0.9307 | 9848 | 1657 | 733 | 0.9918 (82 err) |
| A | 19.2 | 105.4 | 500 | 0.9609 | 0.9594 | 0.9623 | 10182 | 445 | 399 | 0.9933 (69 err) |
| **D = A+B** | 21.5 | 121.4 | 491 | **0.9637** | **0.9689** | 0.9586 | 10143 | **337** | 438 | **0.9967 (34 err)** |
| bcftools | 306.7 (n=1) | 325.7 | 102 | 0.9552 | 0.9619 | 0.9486 | 10037 | 417 | 544 | 0.9772 (234 err) |

D's STR model for HG001 was estimated from the HG001 BAM itself, with no labels. That is legitimate
at test time.

| Run, by repeat context (R / P / F1) | Homopolymer ≥ 7 (4948) | TR ≤ 50 bp (1172) | TR 51–200 bp (549) | Not in TR/HP (3722) |
|---|---|---|---|---|
| Baseline | 0.919 / 0.812 / 0.862 | 0.850 / 0.757 / 0.801 | 0.829 / 0.742 / 0.783 | 0.986 / 0.975 / 0.980 |
| A | 0.960 / 0.950 / 0.955 | 0.928 / 0.937 / 0.933 | 0.907 / 0.925 / 0.916 | 0.984 / 0.986 / 0.985 |
| **D** | 0.955 / **0.966** / **0.961** | 0.922 / **0.949** / **0.935** | 0.893 / 0.932 / 0.912 | 0.984 / 0.986 / 0.985 |
| bcftools | 0.940 / 0.957 / 0.948 | 0.916 / 0.930 / 0.923 | 0.856 / 0.894 / 0.875 | 0.984 / 0.992 / 0.988 |

## Conclusions (what the experiments support)

1. **Multi-allelic genotyping (A) is the main fix.** On validation it gives:
   - FP 1657 → 445 (−73%) and FN 733 → 399, while *raising* recall.
   - Development precision no longer falls with depth (0.941 at 30×, 0.945 at 67×).
   - Runtime and RAM are essentially unchanged.
2. **The STR error model (B) is a smaller, real precision and genotype gain on top of A.**
   - On validation: FP 445 → 337, genotype errors 69 → 34, at the cost of 39 extra FNs and +12%
     wall time.
   - On development at 30× it cost more recall than it gained in precision. Its benefit therefore
     depends on depth and data. The D vs A F1 difference on validation (0.0028) is small, and
     **no confidence interval was computed**, so it is not established that D beats A.
3. **Fragment collapse (C, C2) hurts as implemented.**
   - It halves the evidence per site, but the QUAL ≥ 20 rule was set for read-level evidence.
   - Mates also disagree near read ends: 2 of 35 fragments at 5 true het sites.
   - The correlated-evidence problem is real (85% of pairs overlap; ≈1% unmarked duplicates).
     It needs a fragment likelihood that combines both mates' independent sequencing errors
     while sharing PCR errors, not simple mate selection. Not done.
4. **Compared with bcftools on validation, D has:**
   - higher F1 (0.964 vs 0.955), higher precision (0.969 vs 0.962) and higher recall (0.959 vs
     0.949);
   - 7× fewer genotype errors;
   - about 14× less wall time on 8 threads, and 2.7× less CPU time.

   bcftools remains better outside repeats (F1 0.988 vs 0.985) and uses about 5× less RAM.

**Limits:**
- One validation sample at one depth and on one platform family (Illumina HiSeq).
- No confidence intervals.
- hap.py xcmp instead of vcfeval.
- Indels > 50 bp are not represented in these strata.

## Remaining errors after A (development 30×, from `scripts/analyze.py`)

- **FP (674):**
  - 463 wrong allele length next to a real variant (288 in HP ≥ 7);
  - 145 genotype errors, mostly HP ≥ 7 hom-alt called 0/1. B targets these, and genotype errors
    fall to 92 in D;
  - only 66 FPs with no true variant nearby (746 in the baseline).
- **FN (590):**
  - 165 different allele / representation;
  - 148 multi-allelic sites still missed (K ≤ 5 truncation or footprint split at chunk edges, not
    yet separated);
  - 130 genotype;
  - 116 never a candidate.
- **Clustered variants:** recall is still 0.80 vs 0.96 isolated. Nearby SNPs are not on the
  haplotypes.

## Next steps (ranked)

1. **Compound haplotypes.** Put nearby candidate indels and pileup SNPs on the haplotypes
   (clustered recall 0.80, plus the representation FPs).
2. **Fragment likelihood:** P(frag | H) = P(mate1 | H) · P(mate2 | H) with the PCR-error part of
   the indel error shared, instead of mate selection. Re-test C on that.
3. **Depth-aware or posterior-calibrated QUAL.** Check how calibrated the posterior is (QUAL vs
   empirical FP rate) before any threshold work.
4. **Bootstrap confidence intervals** for D vs A vs bcftools; a second validation sample (HG007
   chr19 is on the HDD) and a second depth.
