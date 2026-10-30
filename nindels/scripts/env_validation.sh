# Independent validation set: HG001 chr19, NHGRI HiSeq 2x148 NovoAlign (~254x),
# downsampled to 30x. Never used for development decisions.
source "$(dirname "${BASH_SOURCE[0]:-$0}")/env.sh" 2>/dev/null || source /home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER/.worktrees/nindels-test/nindels/scripts/env.sh
export SRC_BAM=/home/kramkost/varion_data/hg001_bam/HG001_chr19_native.bam
export TRUTH=/mnt/hdd/AI_DNA_ANALYZER_hg001/truth/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
export CONF_BED=/mnt/hdd/AI_DNA_ANALYZER_hg001/truth/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
export CHROM=chr19
export SAMPLE=HG001
export DATA=$NINDELS_ROOT/data_validation
export RUNS=$NINDELS_ROOT/runs_validation
