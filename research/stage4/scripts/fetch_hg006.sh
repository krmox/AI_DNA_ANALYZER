#!/usr/bin/env bash
# BLIND-NEW: HG006 chr16, HiSeq 2x148 novoalign ~100x GRCh38 (stream one contig). No caller is run on it here.
set -euo pipefail
D=/mnt/archive/AI_DNA_ANALYZER_r2/blind; mkdir -p $D
U=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG006_NA24694-huCA017E_father/NA24694_Father_HiSeq100x/NHGRI_Illumina100X_Chinesetrio_novoalign_bams
BAM=HG006.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam
[ -s $D/HG006_chr16_native.bam ] && exit 0
curl -sf -o $D/$BAM.bai $U/$BAM.bai
for i in 1 2 3; do samtools view -X -b -o $D/HG006_chr16_native.bam.part $U/$BAM $D/$BAM.bai chr16 && break; echo retry $i; done
mv $D/HG006_chr16_native.bam.part $D/HG006_chr16_native.bam && samtools index $D/HG006_chr16_native.bam
echo "HG006 chr16 $U/$BAM $(stat -c %s $D/HG006_chr16_native.bam) $(sha256sum $D/HG006_chr16_native.bam | cut -d' ' -f1)" >> /mnt/archive/AI_DNA_ANALYZER_r2/logs/fetch_blind.log
