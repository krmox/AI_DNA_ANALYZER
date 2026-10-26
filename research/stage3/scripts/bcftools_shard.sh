#!/usr/bin/env bash
# One shard of the pre-specified bcftools baseline (identical flags to research/depth30_baseline).
# usage: bcftools_shard.sh BAM REF CONTIG START END OUT.vcf.gz
set -euo pipefail
bam=$1; ref=$2; contig=$3; start=$4; end=$5; out=$6
bcftools mpileup -f "$ref" -r "$contig:$start-$end" -q 20 -Q 13 -d 10000 -a FORMAT/AD -Ou "$bam" \
  | bcftools call -mv -Ou | bcftools view -v snps -i 'QUAL>=30' -Oz -o "$out"
bcftools index -t "$out"
