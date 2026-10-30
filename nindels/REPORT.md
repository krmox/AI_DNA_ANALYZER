# nindels — first benchmark, accuracy diagnosis and fixes

Date: 2026-10-30. Branch `nindels-test`. All numbers below are measured. Anything not
measured is marked **unknown**.

## 1. What was tested

- **nindels** = independent candidate discovery from BAM CIGARs (≥ 2 supporting reads, f ≥ 0.03)
  → every candidate reassessed with the banded PairHMM → diploid GL(0/0, 0/1, 1/1).
- **No training.** Decision = *untrained rule*: PASS if the genotype is not 0/0 and
  QUAL ≥ 20 (phred). The threshold was fixed before the runs and not tuned. The pattern memory,
  logistic model and aggressiveness A are implemented but **not used here**, because they need
  labelled training data.
- **Baseline:** `bcftools mpileup -a AD,DP | bcftools call -mv` (1.21).
- **Data:** HG004 (GIAB, Illumina HiSeq 2×250, NovoAlign), chr20 (64.4 Mb), native 66.7× and
  seeded subsets at 30× and 15× (`samtools view -s 21.x`).
- **Truth:** NIST v4.2.1, `benchmark_noinconsistent.bed`, 11,000 indels on chr20.
- **Scoring:** hap.py 0.3.15, xcmp engine (rtg/vcfeval not installed). Stratification uses the
  GIAB v3 GRCh38 homopolymer / tandem-repeat BEDs.
- **Machine:** 8 threads, NVMe for code, NTFS HDD for data (page cache warm after run 1).
- **Commands:** `nindels/scripts/{prepare_data.sh,run_callers.sh,happy.sh,analyze.py}`. Raw
  outputs: `/mnt/hdd/nindels/runs/<caller>_<depth>/`.

## 2. Benchmark table (PASS, chr20 confident regions)

| Caller | Wall (s) | CPU (s) | Peak RSS (MB) | SNP F1 | INDEL F1 | INDEL Precision | INDEL Recall | FP | FN |
|---|---|---|---|---|---|---|---|---|---|
| bcftools 15× | 371.1 (n=1) | 339.4 | 110 | 0.9692 | **0.8857** | 0.9242 | 0.8502 | 799 | 1648 |
| nindels 15× | 12.7 (median, n=3) | 71.0 | 329 | N/A | 0.8491 | 0.8274 | **0.8719** | 2001 | 1409 |
| bcftools 30× | 660.5 (n=1) | 581.0 | 110 | 0.9911 | **0.9330** | 0.9483 | 0.9181 | 578 | 901 |
| nindels 30× | 29.2 (median, n=3) | 160.7 | 502 | N/A | 0.8743 | 0.8334 | **0.9194** | 2023 | 887 |
| bcftools 67× | 1019.8 (n=1) | 1019.2 | 111 | 0.9926 | **0.9548** | 0.9627 | 0.9470 | 426 | 583 |
| nindels 67× | 64.9 (median, n=3) | 367.5 | 960 | N/A | 0.8390 | 0.7678 | 0.9248 | 3078 | 827 |

Caveats:
- bcftools was timed once and is single-threaded. It also calls SNPs, so wall time does not
  compare like with like. CPU seconds are the fairer comparison.
- Peak RSS for bcftools is the larger of the two pipeline processes, not their sum.
- nindels calls no SNPs.

**Verdict:**
- nindels' recall is equal to or better than bcftools at every depth (+2.2 points at 15×).
- Its precision is far worse: **3.5× more FPs at 30×, 7× at 67×**.
- Its **precision falls as depth rises** (0.833 → 0.768). A correct caller gains precision with
  depth, as bcftools does (0.948 → 0.963). This points to a broken error model, not to noise.

### By repeat context and length (30×, PASS)

| Stratum | Truth | nindels R / P / F1 | bcftools R / P / F1 |
|---|---|---|---|
| not in TR or HP | 3876 | 0.991 / 0.967 / **0.979** | 0.986 / 0.994 / **0.990** |
| homopolymer ≥ 7 | 4908 | 0.892 / **0.782** / 0.833 | 0.857 / 0.914 / 0.885 |
| tandem repeat ≤ 50 bp | 1449 | 0.846 / **0.728** / 0.782 | 0.946 / 0.949 / 0.947 |
| tandem repeat 51–200 bp | 653 | 0.855 / 0.748 / 0.798 | 0.928 / 0.937 / 0.933 |
| DEL 1–5 | 5066 | 0.922 / 0.807 / 0.861 | 0.907 / 0.948 / 0.927 |
| DEL 6–15 | 610 | 0.818 / 0.746 / 0.781 | 0.934 / 0.920 / 0.927 |
| DEL ≥ 16 | 136 | 0.853 / 0.859 / 0.856 | 0.934 / 0.985 / 0.959 |
| INS 1–5 | 4941 | 0.917 / 0.890 / 0.903 | 0.918 / 0.955 / 0.937 |
| INS 6–15 | 557 | 0.776 / 0.720 / 0.747 | 0.908 / 0.920 / 0.914 |
| INS ≥ 16 | 104 | 0.885 / 0.868 / 0.876 | 0.904 / 0.904 / 0.904 |
| clustered (another truth variant ≤ 10 bp) | 648 | recall 0.829 | unknown (not computed) |
| isolated | 10352 | recall 0.928 | unknown |

Outside repeats, nindels is close to bcftools (F1 0.979 vs 0.990). **Nearly the whole gap is
in homopolymers and tandem repeats.**

- Genotype strata (het / hom-alt / het-alt): hap.py xcmp did not produce them in
  `extended.csv`, so they are **unknown**. The genotype taxonomy below covers part of this.
- Indels > 50 bp: there are none in this truth subset with these bins (hap.py's ≥ 16 bins are
  shown instead).

## 3. Why the accuracy is what it is (nindels 30×, every error classified)

**False negatives (887):**

| Cause | Count | Share |
|---|---|---|
| multi-allelic site (truth 1/2) — alleles genotyped separately | 484 | 56% |
| different allele / representation called nearby | 153 | 18% |
| allele found, genotype wrong | 114 | 13% |
| never a candidate (discovery miss) | 105 | 12% |
| candidate found but filtered | 2 | 0% |

**False positives (2023):**

| Cause | Count | Of which HP ≥ 7 | Called genotype |
|---|---|---|---|
| wrong allele / representation next to a real variant | 1153 | 725 | 0/1: 631, 1/1: 422, 2/1: 100 |
| no truth variant nearby (pure FP) | 746 | 400 | 1/1: 384, 0/1: 359 |
| right allele, wrong genotype | 124 | 100 | 0/1: 92, 1/1: 32 |

Pure FPs vs TPs (medians): NALT 5 vs 15, **NREF 0 vs 9, f 1.00 vs 0.58**. Of the pure FPs,
67% have NREF = 0 and 84% have reads carrying a *different* indel at the same locus (NOTH > 0).
At 20 sampled pure-FP sites, raw BAM depth is 27 and 23 reads are counted. The reads are not
lost; they are assigned to "other alleles".

### Accuracy bottleneck 1 — no multi-allelic genotyping (the largest)
Each candidate allele is genotyped alone against the reference. At a 1/2 site
(e.g. `chr20:232033 A>AAAATAAAT,AAAATAAATAAAT`, truth 1/2), each allele gets
NALT ≈ 6–9, NREF = 0, and the other allele's reads are set aside. Each allele is called 1/1,
and hap.py scores the site as FN plus FPs. This alone explains **56% of FNs**, plus the
100 "2/1" FPs and a large part of the 1/1 representation FPs.

### Accuracy bottleneck 2 — f excludes competing alleles (spec §3 as written)
f = n_alt / (n_alt + n_ref) ignores reads with a different indel at the same locus. In
homopolymers and STRs, slippage spreads reads over ±1-unit alleles. Every allele then shows
f ≈ 1 with NREF = 0, and the diploid model calls it **1/1**. 384 pure FPs and 422
representation FPs are called 1/1.

- Example: `chr20:736736`, two candidate alleles, NALT 14 / 21, NREF 0, NOTH 24 / 17.

*Correction to my earlier claim:* I first blamed the rule that a reference read must span the
whole repeat. The depth check shows that rule drops only about 4 of 27 reads. Excluding
competing alleles is the main cause.

### Accuracy bottleneck 3 — one indel error rate for every context
The PairHMM uses gap-open Q45 (3×10⁻⁵) everywhere. Real slippage in a 10–15 bp
homopolymer is about 10⁻² per read. Under 0/0, two or more slipped reads look almost
impossible, so they get called heterozygous.

- At 30×, 1,225 of 2,023 FPs (61%) are in HP ≥ 7, where precision is 0.782 vs 0.967 outside
  repeats.
- More depth gives more slipped reads, which gives more FPs. That explains why precision falls
  from 0.833 to 0.768 at 67×.
- 2-Mb probe example: `chr20:31270770 C>CA` with NALT 2, NREF 13 (f = 0.13) gets HMM QUAL
  well above threshold (HQ 22), while the count-based genotype gives −0.7.

### Smaller causes
- **Clustered variants:** recall 0.829 vs 0.928 isolated. Each candidate is tested on a
  haplotype that holds only that one indel, so a nearby SNP or indel counts as mismatches
  against both haplotypes.
- **Discovery misses:** 105 at 30×, 31 at 67×. Most are long insertions in tandem repeats that
  the aligner soft-clips (e.g. `chr20:527001 C>CAGACAGAGAGAGAGAG`), or have fewer than 2
  supporting reads at low depth.
- **Band approximation:** the ±(20 + L) DP band changed QUAL at 5 of 2,203 probe candidates,
  all in tandem repeats. No PASS/FAIL change was measured on chr20; that check is **unknown**.

## 4. Performance bottlenecks (measured with stage timers; `perf` is not installed)

1. **PairHMM, before the fix: 99.6% of time** (848 of 851 thread-seconds on a 2-Mb probe).
   The cause was `pow()` in the inner loop plus a full DP matrix.
   - Fixed with a quality lookup table and a ±(20 + L) band.
   - Result: 7 min 38 s → 5.0 s wall (90×). Banded = full DP to 1e-6 on test fixtures.
2. **PairHMM, after the fix: still 87% of time.** At 30×, the engine takes 176.7 thread-s and
   scan + evidence takes 25.9. Every candidate is reassessed, because A is not used without a
   model.
3. **Peak RSS grows with depth:** 329 → 502 → 960 MB. All reads of 8 chunks (8 × 250 kb) are
   held at once, plus a copy of every read for each candidate that overlaps it.
4. **Thread utilisation was 138% before the fix** (per-chunk work units). Now 530–560%, using
   per-candidate scheduling.
5. Cache misses: **unknown** (needs `perf stat -e cache-misses`; `perf` is not installed).
6. I/O: **unknown in isolation.** Scan + BAM decode is about 13% of thread time on a warm page
   cache.

## 5. Fixes, ranked by expected impact / effort

| # | Fix | Targets | Effort | How to measure |
|---|---|---|---|---|
| 1 | **Count competing alleles:** f = n_alt / (n_alt + n_ref + n_other), and count n_other as non-alt in the count genotype | ~800 FPs called 1/1 (pure + representation) | very low (2 lines) | FP called 1/1 should fall; HP and TR precision |
| 2 | **Joint multi-allelic genotyping:** group candidates at one locus or repeat; build one haplotype per allele; compute the PairHMM per read per haplotype; take the GL over all genotype pairs (0/0 … 1/2 … 2/2); emit one record | 484 FN (56%), 100 "2/1" FPs, part of 1,153 representation FPs | medium | hap.py het-alt recall; FP.al |
| 3 | **Context-dependent indel error model:** gap-open and per-read slippage set by homopolymer or STR length and period (e.g. Q45 → about Q20 for a 12-bp homopolymer). The table can be estimated *from the same BAM without labels*: the rate of ±1-unit reads at sites that are clearly hom-ref | 1,225 HP ≥ 7 FPs; the fall in precision with depth | medium | HP ≥ 7 precision; precision at 15 / 30 / 67× should rise with depth |
| 4 | **Local haplotypes for clustered variants:** add nearby candidates (and SNPs from read pileup) to the haplotypes, or do small local assembly in active windows | clustered recall 0.829 → isolated level; representation FPs | high | clustered vs isolated recall |
| 5 | **Discovery from soft clips:** realign clipped tails against the reference near the clip to recover long insertions in STRs; at low depth, allow n_alt = 1 when BQ and MAPQ are high | 105 discovery misses at 30× | medium | "never a candidate" count |
| 6 | **Speed: reassess only uncertain candidates** (the spec's A, or a simple count-GL gate), SIMD/float PairHMM, cache haplotypes per locus | 87% HMM time | low–medium | wall and CPU vs F1 for each gate setting |
| 7 | **Memory:** smaller batches, and store read indices per candidate instead of copies | 960 MB at 67× | low | peak RSS |

Fix 1 is a two-line change and should be done first. Fixes 2 and 3 are where most of the
remaining gap to bcftools lies.

## 6. What the original spec got wrong (measured, not assumed)

- **Δℓ (§8) with one alternative haplotype** fails for heterozygotes. The test fixture shows
  Δℓ ≈ 0 for a clean 0/1 deletion; on real data `DLO` is negative for some true hets
  (e.g. chr21 mock `7000 A>AT`, DLO = −9.9). Replaced by diploid GL.
- **f (§3) as a two-allele fraction** is the direct cause of accuracy bottleneck 2. With
  several indel alleles at one locus, the competing reads must count against each allele.
- **The spec has no multi-allelic sites at all** (one candidate = one REF/ALT pair). Real
  repeats need joint genotyping (bottleneck 1).
- **The pattern memory (§6–7) does not target the main error sources.** The top three errors
  come from the genotyping model and allele representation, not from "this kind of pattern is
  often false". A memory could learn that HP ≥ 7 is risky, but that patches bottleneck 3
  instead of fixing it. **Not tested:** that would need labelled training data, which was out
  of scope.
- **N_eff (§7.1)** does not count independent loci (shown by unit test: 40 records of one
  locus → Kish 40, distinct-locus 1). Untested on real data.
- **Aggressiveness A (§9) and the evidence features (§2.1)** were not evaluated, because there
  is no fitted model. Unknown whether they help.

## 7. Next experiments (in order)

1. Fix 1, rerun 15 / 30 / 67×, compare FP counts by genotype.
2. Fix 3 with the error table estimated from the BAM; check precision vs depth.
3. Fix 2 (joint multi-allelic genotyping); check het-alt recall and FP.al.
4. Gate the HMM to uncertain candidates; measure the speed vs F1 curve.
5. Repeat on an independent sample and platform: HG001 chr19 300× (NVMe), HG007 chr19 (HDD).
   No change has been shown to generalise until it is tested on one of these.
