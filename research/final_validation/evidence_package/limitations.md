# Limitations (expanded benchmark)
1. **Not independent.** HG001 and HG007 were evaluated before this preregistration; 10/20/30× numbers were known when endpoints were written. Only 15×, 50×, HG007 CIs, DeepVariant and the pooled analysis are new.
2. **Two samples × two contigs** (chr17, chr19; confident regions only). No whole-genome run; no hard-region (outside GIAB) claim; short-read Illumina HiSeq novoalign only; GRCh38 only.
3. **SNPs only.** Varion calls no INDELs; INDEL columns in the CSV are for the baselines' completeness only.
4. **CIs** are 2-Mb block bootstraps within a contig. Seeds are subsamples of one BAM, not replicates. The pooled interval covers 4 cells of 2 samples; it is not a population interval.
5. **50×** has one seed. **DeepVariant** has one seed and 2 of 4 cells (user stopped the rest); no paired test; chr17 UNASSESSED.
6. **DeepVariant configuration:** official 1.10.0 image executed through bubblewrap, not Docker (not cross-checked against Docker); whole-contig regions vs Varion's BED call domain; CPU-only; default sharding. Its resource numbers include TensorFlow stages.
7. **Timing:** single machine, single run per condition (3 seeds per condition provide the only repeat), EXISTING 10/20/30× timings carry no load record; HG001 chr19 EXISTING timings overlapped a download. One NEW guardrail cell (HG001 chr19 50×, +18.2% wall) fails and was not repeated.
8. **Data handling:** HG001 natives re-fetched; whole-file sha256 differs (header `@PG` path); identity by BED-mean depth (Addendum A). bcftools 1.21 differs from the 1.23.1 used in the Stage 3 / HG006 runs.
9. **Harness history:** earlier driver bugs (archived in `independent_genome/` and `hg007/`) did not affect these runs; the expanded driver is a parametrised copy of the fixed one.
10. `make test` (full) remains blocked by missing `torch`; ASan/TSan/valgrind were not re-run.
