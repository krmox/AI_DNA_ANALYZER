# Final independent evaluation of V2.x A+C — pre-registration

Written 2026-09-26, BEFORE any V2.x output was produced on this data. No result from this sample/chromosome had been seen.
This is the single, final benchmark. No tuning, no re-selection of region, no repeated attempts to obtain a better number.
If the run fails, the failure is recorded and reported; the method is not modified.

## 1. Data (fixed before the run)
| Item | Value |
|---|---|
| Sample | HG003 / NA24149 (GIAB Ashkenazim Trio, father). NOT HG002 (the only sample used for V2.x development). |
| Region | **entire chromosome chr8** (chr8:1-145,138,636; GRCh38), one run, no sub-selection. |
| Sequencing | Illumina HiSeq 2500, 2x148 bp PCR-free, NovoAlign V3.02.07, nominal ~300x (genome-wide GIAB dataset "NHGRI_Illumina300X_AJtrio_novoalign_bams"). Achieved depth is measured and reported after the run (descriptive only). |
| BAM | `HG003.GRCh38.300x.bam` from NCBI GIAB FTP (md5 of the 569 GB parent listed in the GIAB alignment index), chr8 slice obtained by ranged reads: 16 disjoint slices by read start position, concatenated (`download_bam.sh`). |
| Reference | GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna (same file as chr20 benchmark; contig `chr8`). |
| Truth | GIAB v4.2.1 GRCh38 `HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`; confident regions `..._benchmark_noinconsistent.bed`. SHA-256 in `truth/SHA256SUMS.txt`. |
| Scorer | hap.py 0.3.15 (xcmp, `--threads 8`, `-f` confident BED, `-l chr8`) — the same tool/settings as the chr20 benchmarks; SNP-restricted results are also reported as hap.py's SNP row (the engine is SNP-only, as in V1/V2). |

## 2. Model (fixed)
V2.x **A+C** only: engine `build_x/dnav2` (sha256 `adb2f4a0…bb7bdb`), model `native_v2x/models/A_C.model` (sha256 `e9942a22…2d51ce`, identical to `research/V2X_MODEL_FILES_SHA256.csv`), worktree `/mnt/archive/AI_DNA_ANALYZER_v2x/repo` branch `native_v2x` on base commit `ca60e12e19b1c135fd12ee790ebee9b6628e47ba` plus UNCOMMITTED `native_v2x/` (source tree digest recorded in `hardware_and_model.txt`).
Frozen constants unchanged: Binomial 7.0, router cutoff 5.411872376933351, PB 10.5, Method C eps 0.01, MAX_READS 48, SNP-only, 64-bp window frame.
Command (only `--contig/--contig-len/--sample/--tag` differ from the chr20 runs; no parameter or model coefficient is changed):

    /usr/bin/time -v dnav2 run --bam HG003.GRCh38.300x_chr8.bam --ref <GRCh38 analysis set fna> --truth HG003_..._benchmark.vcf.gz \
        --bed HG003_..._benchmark_noinconsistent.bed --out-dir runs/A_C_chr8 --contig chr8 --contig-len 145138636 \
        --sample HG003 --tag chr8_300x --threads 8 --mode cascade --model native_v2x/models/A_C.model

Not run (by instruction): V1, V2 control, PB-only, other V2.x variants, ablations. Single run, warm/cold page-cache state recorded as observed (caches are not dropped).

## 3. Why this is independent of V2.x development
* V2.x features, logistic coefficients, threshold, and the choice of A+C were fitted / selected only on even 2-Mb blocks of **HG002 chr20** (dev) and confirmed on odd blocks (holdout). No HG003 data and no chr8 data was read by any V2.x analysis script (`native_v2x/analysis/*` load only chr20 HG002 dumps).
* **Sample-level separation**: HG003 is a different individual, different library/flowcell than HG002. Caveat: HG003 is the father of HG002 (GIAB trio); the two share the sequencing platform, aligner and reference, so this is not a cross-platform test.
* **Chromosome-level separation**: chr8 was never used in the V2.x programme. It was also not among the chromosomes of the earlier V1 stress cells (chr1,2,3,4,5,7,13-17,19,20,21 in `BENCHMARK_MASTER.csv`); the frozen V1 constants were not fitted on it.
* Limits of independence: same sequencing platform and aligner as the development data; germline diploid autosome; GIAB high-confidence regions only.
* Selection of the region was made before the run: chr8 is a large autosome (145 Mb, ~2.3x chr20), untouched by any previous experiment; whole chromosome, so there is no region selection inside it.

## 4. Addendum (written before the run; adds descriptive metrics only, changes nothing above)
* Expected chr8 read count from the remote BAM index (`samtools idxstats`, logs/remote_idxstats_chr8.txt): 297,503,843 mapped + 2,622,496 unmapped-placed reads. The local slice must reproduce these two counts exactly; a mismatch is reported, and the run then does not start.
* Depth is reported descriptively after the run (index read count x read length / contig length, and `samtools coverage` on the local slice). It is not used for any decision.
* Additional descriptive metrics, computed after the single run with a method fixed here: hap.py SNP precision/recall/F1 with TP/FP/FN, FP.gt / FP.al, Ti/Tv and het/hom of the query, QUERY.TOTAL / UNK; a 95 % paired-free block bootstrap of F1 / precision / recall (whole 2-Mb blocks of chr8 as resampling units, 4,000 resamples, seed 20260926, computed from the per-record hap.py TP/FP/FN annotation), the same estimator family as used for the V2.x holdout; number of loci processed, candidates, calls, PB evaluations from the engine's JSON; wall time, CPU time, peak RSS from `/usr/bin/time -v`.
* One engine run, one hap.py run. No repeat is scheduled. If either fails, the failure is documented as is.
