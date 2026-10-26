#!/usr/bin/env bash
# Stream one contig of a GIAB Illumina 2x250 novoalign BAM (GRCh38) from NCBI FTP using the local .bai.
# usage: fetch_contig.sh HG002|HG003|HG004 CONTIG
set -euo pipefail
S3=/mnt/archive/AI_DNA_ANALYZER_stage3
REPO=/home/mark/Documents/Projects/AI_DNA_ANALYZER
sample=$1; contig=$2
case $sample in
  HG002) dir=HG002_NA24385_son;;
  HG003) dir=HG003_NA24149_father;;
  HG004) dir=HG004_NA24143_mother;;
esac
url=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/$dir/NIST_Illumina_2x250bps/novoalign_bams/$sample.GRCh38.2x250.bam
out=$S3/bam/${sample}_${contig}_native.bam
[ -f "$out" ] && { echo "exists: $out"; exit 0; }
samtools view -X -b -o "$out.part" "$url" "$REPO/$sample.GRCh38.2x250.bam.bai" "$contig"
mv "$out.part" "$out"
samtools index "$out"
echo "$sample $contig $url $(stat -c %s "$out") $(sha256sum "$out" | cut -d' ' -f1)" >> $S3/logs/fetch.log
