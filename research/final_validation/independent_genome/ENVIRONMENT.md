# ENVIRONMENT (HG001 independent-genome campaign)

Date: 2026-10-28  Host: Linux 7.3.0-4-cachyos  CPU: Intel(R) Core(TM) i7-6700 CPU @ 3.40GHz (8 threads)  RAM: 31 GB

## Frozen artefacts (unchanged; original files used directly)
- Repo HEAD: 3cdfe5abef503ce5fb779dc471cf571a5a27420f (3cdfe5a), tracked tree clean (only untracked research/final_validation outputs)
- native_r2 sources: 21/21 files match /mnt/hdd/AI_DNA_ANALYZER_r2/gates/source_sha256.txt
- R2 binary: /mnt/hdd/AI_DNA_ANALYZER_r2/build/dnav2 sha256 6c675174854917ca3bbbfc1192c9205ac0043a696bd37d6ff2b2c08c3a3304be (R2_LOCK: 6c675174…3304be)
- R2 model: research/stage4/frozen_R2.model sha256 649215d52e57cf00d9f69b8631d7c0d9a1848db2550fc593b7d464c1b4cab6ac
- R1.1 binary: /mnt/hdd/AI_DNA_ANALYZER_stage3/build/dnav2 sha256 e83d87fd2f810eb7a9d096fb797c6d3451aaef86d14a0071dc2ad9b91a3a39df (PREREG: e83d87fd)
- R1.1 model: research/stage3/frozen/R1.model sha256 63816707c4af7794138c3917ddcc22ddce67f800f724db8a004ad60ee4b711ed (PREREG: 63816707)
- R2 command: dnav2 run --bam B --ref REF --bed CONTIG.bed --out-dir D --contig C --contig-len L --sample HG001 --tag COND --threads 8 --mode cascade --model frozen_R2.model --r2   (identical to research/stage4/scripts/run_r2.py; no --truth)
- R1.1 command: same without --r2 and with R1.model
- bcftools baseline: research/stage3/scripts/bcftools_sharded.sh (32 shards, 8 parallel): mpileup -q 20 -Q 13 -d 10000 -a FORMAT/AD | call -mv | view -v snps -i 'QUAL>=30'
- Original binary needs libhts.so.3: run with LD_LIBRARY_PATH=/home/kramkost/varion_mamba/envs/vt/lib (htslib 1.21); binary itself is not modified

## Tools (user-space micromamba, no root)
- samtools 1.21, bcftools 1.21, htslib bgzip (htslib) 1.21
- hap.py: Hap.py  (bioconda hap.py 0.3.15 py27, xcmp engine; same version as the container used in earlier stages)
- Python: Python 3.14.8 (driver), Python 3.14.10 (analysis, numpy 2.5.4)
- hap.py command: hap.py TRUTH.vcf.gz QUERY.vcf.gz -f CONTIG.bed -r REF -o PREFIX -l CONTIG --threads 4 (no stratification, no vcfeval; default left-normalisation; SNP Filter=PASS rows of h.summary.csv)

## Data identifiers
- Reference: /mnt/hdd/AI_DNA_ANALYZER_benchmark/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna (GRCh38 no_alt + hs38d1; sha256 592502666381c995…); same file for alignments (BAM header uses GRCh38_full_plus_hs38d1_analysis_set_minus_alts; contigs chr1..; identical primary contig names/lengths)
- Reads: GIAB NA12878/HG001 Illumina HiSeq 2x148 300x novoalign (v3.02.07) GRCh38: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam
- Per-contig native BAMs (streamed with samtools view -X), sha256 in /mnt/hdd/AI_DNA_ANALYZER_hg001/logs/fetch.log:
    chr19 10293591289 17c6e14e69f9800366abe5e2993be46d056747e6599d2a5f47636f256ef05368
    FETCH_DONE
- Truth: GIAB release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.{vcf.gz,bed} (no *_noinconsistent.bed provided for HG001; same .bed convention as HG005–HG007); sha256:
    dc3485e60447a3e863c34dfc063c3710ea5ed5efae50855c5197ba62538263b9  HG001_GRCh38_1_22_v4.2.1_benchmark.bed
    93bc4c2c696eaf13515ab058caecc064bfed704f85bac7482330ca91bc730daa  HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
    68a91b78042e66ed9a53ae4069ac7686abf02ea9a1f88da00e0054a7db06dccc  HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

## Downsampling
samtools view -b -s SEED.FFFF (FFFF = round(10000 × target / native BED-mean-depth)), native BED-mean depth by samtools bedcov over the contig-restricted confident BED. Same recipe as research/stage3/scripts/stage3.py cmd_derive. Recorded per run: flag, fraction, BAM sha256, achieved BED-mean depth (RESULTS.jsonl).
