# Paths shared by the nindels benchmark scripts. Data and outputs live on the HDD.
export NINDELS_ROOT=/mnt/hdd/nindels
export DATA=$NINDELS_ROOT/data
export RUNS=$NINDELS_ROOT/runs
export TOOLS=$HOME/varion_mamba/envs/vt/bin        # samtools / bcftools / tabix 1.21
export HAPPY_ENV=$HOME/varion_mamba/envs/hp        # hap.py (python 2.7)
export NINDELS_BIN=/home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER/.worktrees/nindels-test/nindels/build/nindels-call

# Inputs already on the HDD from earlier AI_DNA_ANALYZER work (GIAB, public).
export SRC_BAM=/mnt/hdd/AI_DNA_ANALYZER_stage2/bam/HG004_chr20_full.bam   # HG004 Illumina 2x250 NovoAlign, 66.7x
export REF=/mnt/hdd/AI_DNA_ANALYZER_benchmark/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
export TRUTH=/mnt/hdd/AI_DNA_ANALYZER_stage2/truth/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
export CONF_BED=/mnt/hdd/AI_DNA_ANALYZER_stage2/truth/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
export STRAT_DIR=/mnt/hdd/AI_DNA_ANALYZER_stage3/strat
export CHROM=chr20
