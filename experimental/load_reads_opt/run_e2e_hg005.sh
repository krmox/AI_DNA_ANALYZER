#!/bin/bash
# BAM -> extraction -> Binomial -> PB -> Method C -> VCF on HG005 chr1:1.0-4.0 Mb, new code; per worker count.
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../.." && pwd)"  # project root (path-independent)
W=$ROOT; cd $W
OUT=experimental/load_reads_opt/results/e2e; mkdir -p $OUT cache/e2e_hg005
for workers in "$@"; do
  tag=w$workers
  /usr/bin/time -v python3 extract_bench_v12.py --fasta experimental/stress_test/shared_ref/chr1_chrname.fa \
     --bam data/giab_hg005_stress/chr1_24Mb_30x.bam --vcf data/giab_hg005_stress/chr1_24Mb.vcf.gz \
     --bed data/giab_hg005_stress/chr1_24Mb_highconf.bed --region 1000000 4000000 --contig chr1 \
     --chunk-bp 500000 --features none --workers $workers --out cache/e2e_hg005/$tag.npz > $OUT/extract_$tag.log 2> $OUT/extract_$tag.time
  mkdir -p $OUT/pipe_$tag
  PIPE_CACHE=cache/e2e_hg005/$tag.npz PIPE_OUT=$W/$OUT/pipe_$tag PIPE_CONTIG=chr1 PIPE_CONTIG_LEN=248956422 PIPE_SAMPLE=HG005 PIPE_TAG=hg005_$tag \
     python3 experimental/chr20_validation/run_chr20_pipeline.py > $OUT/pipe_$tag.log 2>&1
done
