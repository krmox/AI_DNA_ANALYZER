#!/usr/bin/env bash
# Stage 2 confirmatory data (HG004 2x250, ~56x genome-wide). Usage: fetch_hg004.sh CONTIG
# Full contig is streamed once; ~30x and ~6x are derived locally with fixed seed 21 (selection is a hash of read name + seed).
set -euo pipefail
D=/mnt/archive/AI_DNA_ANALYZER_stage2
B=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
BAM=$B/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.GRCh38.2x250.bam
IDX=/home/mark/Documents/Projects/AI_DNA_ANALYZER/HG004.GRCh38.2x250.bam.bai
c=$1
samtools view -X -b -o $D/bam/HG004_${c}_full.bam $BAM $IDX $c
samtools index $D/bam/HG004_${c}_full.bam
for spec in f053:21.53 f010:21.10; do
  n=${spec%%:*}; s=${spec#*:}
  samtools view -b -s $s -o $D/bam/HG004_${c}_${n}_s21.bam $D/bam/HG004_${c}_full.bam
  samtools index $D/bam/HG004_${c}_${n}_s21.bam
done
for f in $D/bam/HG004_${c}_*.bam; do samtools coverage $f > ${f%.bam}.coverage.txt; done
echo "$c done" >> $D/logs/fetch.log
