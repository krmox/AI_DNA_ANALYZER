# Varion Stage 3 — evaluation plan (written 2026-10-26, before the full evaluation)

Evaluation only. The caller, model, thresholds, filters and flags are frozen (`IMPLEMENTATION_LOCK.md`). Nothing in this stage may change a prediction.

## 1. Questions
1. Does Varion R1.1 keep SNP F1 ≥ 0.990 outside the conditions used during development?
2. How do precision/recall move with depth (10–~100x); smooth degradation or a failure regime?
3. How does it vary by sample, contig and GIAB difficulty stratum?
4. How does it compare with bcftools (pre-specified baseline) in accuracy, wall/CPU time and memory?

## 2. Data (all GRCh38, reference `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna`, contig names `chrN`)
| Set | Samples | Contigs | Source reads | Independence class |
|---|---|---|---|---|
| **N1: new contigs, trio** | HG002, HG003, HG004 | chr22, chr18 | GIAB NIST Illumina HiSeq 2x250 novoalign, ~56x genome-wide (FTP `data/AshkenazimTrio/*/NIST_Illumina_2x250bps/novoalign_bams`) | **Region-level independent only.** Samples were exposed (HG002 chr20 = R1 training; HG003 chr8 = R1 test + H1 audit; HG004 chr20/21 = Stage 2 confirmatory). chr22 and chr18 appear in no earlier experiment of any model generation (repository grep). The trio is related; the three samples are not independent of each other. |
| **N2: unrelated sample** | HG005 | chr22 (full ladder), chr18 (reduced) | GIAB HiSeq 2500 2x250 novoalign 300x (`data/ChineseTrio/HG005_NA24631_son/.../NHGRI_Illumina300X_Chinesetrio_novoalign_bams`, `...300x.bam`) | **Sample-level new for R1/R1.1** (never used to train, select thresholds or diagnose R1/R1.1). Project-level touch: a V1-era stress test on chr1:1.0–4.0 Mb (different model, invalidated pipeline). Same lab/platform/aligner family as training data, so not technology-independent. |
| **E1: exposed anchors** | HG004 | chr20, chr21 | existing Stage 2 BAMs (native ~56x, `-s 21.53` ≈30x, `-s 21.10` ≈6x) | Exposed (Stage 2 confirmatory set; contig chr20 also R1 training contig on HG002). Re-run through the Stage 3 pipeline for continuity; **not** counted as independent evidence. |
Not available / not done: an independent sequencing technology (e.g. NovaSeq/PCR-free 2x150 for the trio at depth is available but unused here), WGS-scale runs, chrX/chrY, other population samples with GRCh38 v4.2.1 truth other than HG001/HG006/HG007 (no suitable Illumina depth fetched in this stage).

Truth/BED: GIAB v4.2.1 GRCh38 `*_benchmark.vcf.gz` and, for HG002–4, `*_benchmark_noinconsistent.bed` (HG005 release has only `*_benchmark.bed`). Stratifications: GIAB genome-stratifications v3.6 `GRCh38@all` (46 pre-selected strata, `strat/selected_strata.txt`, md5 verified against the release md5 list). Strata were selected by class (coding, GC, homopolymer, tandem repeat, satellite, mappability, segdup, MHC/KIR/VDJ/L1H, difficult unions) before any result existed.

## 3. Coverage design
Observed depth = mean `samtools bedcov` depth inside the sample/contig confident BED. Native BED depth is measured first; a target T uses `samtools view -s SEED.FRACTION` with FRACTION = T / native_depth (deterministic hash on read name).
* Trio (N1) and HG005 chr22: native; 50x (seed 101, one repeat: ≥75 % of reads shared with native); 30x, 20x, 15x, 10x with seeds 101, 102, 103.
* HG005 additionally 100x (seeds 101–103) derived from a ~110x base subsample of the 300x library (`-s 100.NNN`; lower depths are derived from that base).
* HG005 chr18 (reduced to bound cost): 100x s101, 30x ×3, 10x s101.
* E1: the three existing Stage 2 BAMs per contig.
Seeds are technical replicates of one library, not independent samples.
Target 15x/10x are below the depth range where R1 was trained (≥28x); they probe out-of-distribution behaviour and are reported in the primary tables.

## 4. Callers (same BAM, reference, BED, hap.py settings)
* **Varion R1.1**: `dnav2 run ... --mode cascade --model R1.model --threads 8` (defaults: `--edge-windows --truth-free`, no `--truth` argument). SNPs only. Needs the confident BED to define the call domain (an input, not a label).
* **bcftools** (pre-specified baseline, identical to `research/depth30_baseline`): `bcftools mpileup -f REF -r REGION -q 20 -Q 13 -d 10000 -a FORMAT/AD -Ou | bcftools call -mv -Ou | bcftools view -v snps -i 'QUAL>=30'`, version 1.23.1. Primary comparison runs it as 32 contiguous shards on 8 parallel processes (equal 8-thread budget) and concatenates; sharded output is verified against the unsharded run on one condition. No other bcftools configuration is run in the primary tables. Single-process timings are measured separately.
Both outputs are SNP-only and unfiltered by hap.py beyond the callers' own filters (Varion's VCF is `PASS`; bcftools QUAL≥30 as above). Neither output is normalised differently; hap.py performs its standard normalisation for both.

## 5. Evaluation
hap.py 0.3.15 (image `quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0`, engine `xcmp`, default), `-f` confident BED, `-l CONTIG`, `--stratification` TSV of the selected strata, `--threads 4`. Headline row: `Type=SNP, Filter=PASS, Subset=*`. TP := `TRUTH.TP`, FN := `TRUTH.FN`, FP := `QUERY.FP`; `QUERY.TP` also stored. F1 := 2·TP/(2·TP+FP+FN) from integers (hap.py's own `METRIC.F1_Score` stored for cross-check). Genotype-level: hap.py's `FP.gt` and genotype-mismatch counts are stored. INDEL rows are stored for completeness; Varion does not call indels (INDEL recall 0 by construction), bcftools indels are filtered to SNPs by the baseline command.
Uncertainty: Wilson 95 % intervals for recall/precision (locus-level, optimistic because loci are correlated); for primary conditions a 2-Mb block bootstrap (4,000 resamples, seed 20261027) of F1 and of the paired F1 difference Varion − bcftools using the hap.py annotated VCFs (the approach used in Stage 2). Seeds are summarised as mean and range, never as independent replicates. Strata with < 100 truth SNPs are flagged "sparse".
Aggregates: no pooled headline. Per sample×contig×depth rows are primary. An optional pooled count-based F1 (sum of TP/FP/FN over a stated set of conditions) is shown next to them with its weighting.

## 6. Pre-specified interpretation rules
* "Retains F1 ≥ 0.990" is judged per sample×contig at ~30x (mean over three seeds, with the range), on SNP PASS in the confident region. A single sample×contig below 0.990 is reported as a failure of the claim for that condition, not averaged away.
* Depth regimes: the depth at which F1 falls more than 0.005 below its 30x value is reported as a degradation point; no smoothing/fitting beyond reporting.
* Varion vs bcftools: difference is called "not distinguishable" if the paired block-bootstrap 95 % CI includes 0.
* Verdicts on strata use TP/FP/FN counts in the stratum; sparse strata are labelled, not interpreted.

## 7. Performance protocol
Same machine (Intel i7-6700, 4C/8T, 31 GB, governor `powersave`/intel_pstate, desktop session with background load recorded via `/proc/loadavg` before every timed run). Sequential execution, no concurrent job during timing, `/usr/bin/time -v` (wall, user+sys CPU, peak RSS; peak RSS for the sharded bcftools is the largest single process, not the sum). Calling time excludes hap.py and downsampling; bcftools time includes the concat step. Accuracy runs (Phase 1) give one timing observation per condition, taken while a background download may be active; the formal timing campaign (Phase 3) runs five interleaved repetitions on a quiet-as-possible machine: HG002 chr22 {10x, 20x, 30x, native}, HG002 chr18 {30x, native}, HG005 chr22 100x; bcftools single-process (1 thread) five repetitions at HG002 chr22 30x.

## 8. Phases
1. Freeze (`IMPLEMENTATION_LOCK.md`), regression gates.
2. Fetch data, derive depths, run callers (accuracy).
3. hap.py for every caller×condition.
4. Timing campaign.
5. Aggregation (`collect.py`), bootstrap, report.

## 9. Known limitations (declared up front)
SNP-only caller; call domain requires the confident BED; one sequencing platform/aligner family; subsamples of single libraries; trio relatedness; HG005 shares lab/platform with training data; hap.py/truth bias toward easier genome (truth excludes the hardest regions), so strata fully outside the confident BED cannot be evaluated; chrX/Y and WGS-scale behaviour not tested; background load on a desktop machine.
