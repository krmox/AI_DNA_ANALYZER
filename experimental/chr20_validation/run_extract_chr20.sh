#!/bin/bash
# chr20 holdout extraction: BAM -> native pileup counts -> native read tensor -> Binomial LLR -> PB LLR, whole chromosome.
set -eu
cd "$(dirname "$0")/../.."
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../.." && pwd)"  # project root (path-independent)
M=$ROOT
D=$M/data/giab_hg002_chr20
/usr/bin/time -v python3 extract_bench_v12.py \
  --fasta $M/data/reference/chr20_full.fa --bam $D/chr20_15x.bam --vcf $D/chr20.vcf.gz --bed $D/chr20_highconf.bed \
  --region 0 64444167 --contig chr20 --chunk-bp 500000 --features none --workers ${WORKERS:-6} \
  --out cache/chr20_validation/chr20_15x.npz > logs_chr20/extract_chr20.log 2> logs_chr20/extract_chr20.time
echo EXTRACT_DONE > logs_chr20/extract.done
