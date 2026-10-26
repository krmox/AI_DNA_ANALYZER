# Data locations for V2 runs (all on the HDD).  Nothing here is written to the main disk.
export A=/mnt/archive/AI_DNA_ANALYZER_benchmark
export V2W=/mnt/archive/AI_DNA_ANALYZER_v2x/work
export DNAV2=/mnt/archive/AI_DNA_ANALYZER_v2x/build/dnav2
export BAM=$A/HG002_GRCh38_chr20/HG002.GRCh38.300x_chr20.bam
export REF=$A/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
export TRUTH=$A/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
export BED=$A/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
export V1RUN=$A/ai_dna_analyzer_run
