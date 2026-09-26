# Final independent evaluation of V2.x A+C — HG003 chr8, 300×

Date of run: 2026-09-26. Status: **single run, completed, no failure.** No parameter, model coefficient, threshold, feature or preprocessing step was changed, before or after the run. No V1, V2 control, PB-only, ablation or other-model run was made (by instruction). Nothing was committed.

The design was written before any V2.x output existed for these data: `PREREGISTRATION.md` (sections 1–3; the addendum in section 4 was also written before the run and only adds descriptive metrics).

## 1. Fixed before the run

| Item | Value |
|---|---|
| Sample | HG003 / NA24149 (GIAB, father of HG002). V2.x used no HG003 data in any way. |
| Region | Whole chr8 (chr8:1–145,138,636, GRCh38). Not used by V2.x development (HG002 chr20 only) and not among the chromosomes of the V1 stress cells. |
| Reference | `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna` (same file as the chr20 300× benchmark) |
| Sequencing depth | ~300× nominal (NHGRI Illumina 300× NovoAlign V3.02.07 BAM, NIST HiSeq Homogeneity dataset). **Measured after the run: mean depth 302.693× over the contig** (`samtools coverage`: 297,503,843 reads, 99.16 % of positions covered, mean BQ 35.5, mean MAPQ 69.5); reads are 148 bp. |
| Truth set | GIAB v4.2.1 GRCh38 `HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`, confident regions `..._benchmark_noinconsistent.bed` (SHA-256 in `logs/truth_sha256.txt`) |
| Model | V2.x A+C: engine `dnav2` SHA-256 `adb2f4a00f60a476fca1ec43c5af2d8c6b86e5876e8c95207c297f943fbb7bdb`, model `A_C.model` SHA-256 `e9942a22e6fee7980b3881750316291bb0d62b4c8be000377f72a2280c5d51ce` (equals `research/V2X_MODEL_FILES_SHA256.csv`); worktree base commit `ca60e12e19b1c135fd12ee790ebee9b6628e47ba` + uncommitted `native_v2x/` (source digest in `hardware_and_model.txt`) |
| Parameters | frozen: Binomial 7.0, router cutoff 5.411872376933351, PB 10.5, Method C ε = 0.01, MAX_READS = 48, SNP-only; V2.x candidate rule k ≥ 3; 8 threads; mode `cascade` |
| Exact command | `run/exact_command.txt` (only `--contig chr8 --contig-len 145138636 --sample HG003 --tag chr8_300x` differ from the chr20 runs) |
| Scorer | hap.py 0.3.15, xcmp, `-f` confident BED, `-l chr8`, `--threads 8` (rootless Podman; image `quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0`, manifest digest `sha256:d63b963a…f96d0`, identical to the digest recorded for the earlier hap.py runs) |
| Why independent of V2.x development | Different sample, different chromosome. All V2.x features, coefficients, thresholds and the choice of A+C were fitted or selected on HG002 chr20 (even 2-Mb blocks) and confirmed on the odd blocks. Limits: same platform family (NovoAlign, HiSeq, 148 bp), same depth class, HG003 is the father of HG002, autosome, high-confidence regions only. |

## 2. Data integrity

BAM slice: 16 disjoint start-position slices from the official NCBI file, concatenated (`scripts/download_bam.sh`). `samtools quickcheck` OK. Read counts of the local slice equal the official index of the parent file exactly: 297,503,843 mapped, 2,622,496 unmapped-placed (`logs/remote_idxstats_chr8.txt`). SHA-256 `248bb9aad482a342ed5bcbb9ced9de45eb828edefadf3ee90ffa1548489e4626` (BAM), `941f5519…8f14` (BAI). There is no official checksum for this slice.

## 3. Result (hap.py 0.3.15, SNP, PASS = ALL)

| Metric | Value |
|---|---|
| TRUTH.TOTAL (SNP) | 186,522 |
| TP / FP / FN | **182,904 / 446 / 3,618** |
| QUERY.TOTAL / UNK | 183,350 / 0 |
| Precision | **0.997567** |
| Recall | **0.980603** |
| F1 | **0.989012** |
| 95 % block-bootstrap CI (73 whole 2-Mb blocks, 4,000 resamples, seed 20260926) | F1 [0.98713, 0.99048]; precision [0.99655, 0.99828]; recall [0.97770, 0.98299] |
| FP.gt / FP.al | 91 / 65 |
| Ti/Tv (query / truth) | 1.9322 / 1.9302 |
| het/hom (query / truth) | 1.7008 / 1.6968 |
| INDEL rows | 0 calls by construction (SNP-only caller) |

Recomputation of P/R/F1 from the per-record annotation of the hap.py VCF gives identical values (`results_chr8_analysis.json`).

## 4. Processing and cost (one run, 8 threads; caches not dropped)

| Metric | Value |
|---|---|
| Loci extracted / scored | 133,987,200 (2,093,550 windows) / 133,935,319 |
| Candidate loci (k ≥ 3) | 338,184 (0.25 % of scored) |
| Calls (`n_calls`) | 183,350 (183,334 from the classifier path); 219 genotypes forced 0/0 → 0/1 |
| Routed loci / PB evaluations | 1,018 / 102 |
| Wall time | 412.485 s (6:52.51), 324,829 loci/s |
| CPU user + system | 1,903.599 + 26.569 = 1,930.168 s (467 % of one core) |
| Peak RSS | 333,900 kB (326 MB) |
| Output VCF | `run/ai_cascade_chr8_300x.vcf.gz` (bgzip copy; decompressed SHA-256 `383761cbbfe2dd16cca9a516ae81485634a81a2422c4caeaee3fb7d62fed32e5`, 8,648,703 bytes) |
| Full stdout / stderr | `run/stdout.json`, `run/stderr.txt`, `run/time.txt` (`/usr/bin/time -v`) |

The BAM was read from the 5400-rpm data disk with part of it in the page cache; CPU use was 467 %, against 641 % in the warm chr20 runs. The timing is a single observation and is not comparable with the warm-cache medians of the chr20 benchmark.

## 5. Frame audit (M-1), descriptive, after the run

`scripts/frame_audit.py` rebuilds the 64-bp window list from the confident BED (global tiling, whole windows inside the BED) and reproduces the engine's window count exactly (2,093,550). 2,873 of the 186,522 truth SNPs (1.54 %) lie outside all retained windows; all 2,873 are false negatives (79.4 % of the 3,618). Recall is bounded by 0.98460 by construction. The remaining 745 false negatives are in retained windows and were not decomposed.

## 6. What this run does and does not show

It shows that a model fitted on HG002 chr20 produced F1 0.989012, precision 0.997567 and recall 0.980603 on HG003 chr8, a different sample and chromosome of the same platform family and depth class. The holdout F1 on HG002 chr20 (0.989017) lies inside the chr8 interval.

It does not show the gain over V2 on HG003 (no V2 control was run), any comparison with other callers on this data, variance between samples (one sample, one run), or transfer to other platforms, depths, ancestries or sample types. It is a sample-level and chromosome-level separation from the development data, not an independent biological replication and not a population-level validation.

## 7. Problems and deviations (none affected the result)

* The engine warns `[W::hts_idx_load3] The index file is older than the data file` for the truth VCF (download timestamps). The index loaded and worked. The warning is in `run/stderr.txt`.
* The first version of `frame_audit.py` assumed the window list was tiled per 500-kb chunk and reproduced 2,093,446 windows against the engine's 2,093,550. It was corrected after reading `native_v2x/src/frame.cpp` (global 64-bp tiling), which then matched exactly. The script only audits the frame after the run. It does not touch the engine, the model or the VCF.
* A large `samtools coverage` job was run in parallel with hap.py, after the engine run had finished.

## 8. Files

`PREREGISTRATION.md`, `hardware_and_model.txt`, `run/` (command, JSON, time, VCF, checksums), `happy/` (hap.py summary, extended, annotated VCF, log), `results_chr8_analysis.json`, `results_chr8_frame_audit.json`, `scripts/`, `logs/`, `checksums.txt` (SHA-256 of every file in this directory).
