#!/usr/bin/env bash
# bcftools baseline split into SHARDS contiguous regions run with P parallel processes, then concatenated.
# usage: bcftools_sharded.sh BAM REF CONTIG LEN OUTDIR SHARDS P
set -euo pipefail
bam=$1; ref=$2; contig=$3; len=$4; out=$5; shards=$6; par=$7
here=$(dirname "$(readlink -f "$0")")
mkdir -p "$out/shards"
for ((i = 0; i < shards; i++)); do
  printf '%s %s %s %s %s %s\n' "$bam" "$ref" "$contig" $((i * len / shards + 1)) $(((i + 1) * len / shards)) "$out/shards/$(printf %03d "$i").vcf.gz"
done | xargs -P "$par" -L 1 "$here/bcftools_shard.sh"
bcftools concat -Oz -o "$out/bcftools.vcf.gz" "$out"/shards/*.vcf.gz
bcftools index -t "$out/bcftools.vcf.gz"
rm -rf "$out/shards"
