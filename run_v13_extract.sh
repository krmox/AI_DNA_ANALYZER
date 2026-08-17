#!/bin/bash
# Group A extractions, pre-registered in History/13_DEVLOG.md §4.2.
# Each region is extracted exactly once, with --legacy-pb so the depth-clipping
# fix can be A/B compared on identical reads.
set -u
FA=data/reference/chr21_full.fa
V33=data/giab_hg002_real_train/hg002_chr21_31_33M.vcf.gz
B33=data/giab_hg002_real_train/hg002_chr21_31_33M_highconf.bed
V30=data/giab_hg002_real_30M/hg002_chr21_30M.vcf.gz
B30=data/giab_hg002_real_30M/hg002_chr21_30M_highconf.bed

run () { # name bam vcf bed start end
  echo "START $1 $(date -Is)"
  python3 extract_bench_v12.py --fasta $FA --bam "$2" --vcf "$3" --bed "$4" \
    --region "$5" "$6" --features none --legacy-pb \
    --out "cache/bench_v13/$1.npz" > "logs_v13/$1.log" 2>&1
  echo "DONE $1 rc=$? $(date -Is)"
}

run new32_33M_full  data/giab_hg002_real_train/hg002_chr21_31_33M.bam     $V33 $B33 32000000 33000000 &
run re31_32M_full   data/giab_hg002_real_train/hg002_chr21_31_33M.bam     $V33 $B33 31000000 32000000 &
run r30M_full       data/giab_hg002_real_30M/hg002_chr21_30M.bam          $V30 $B30 29999999 30469999 &
run new32_33M_30x   data/giab_hg002_real_30x/hg002_chr21_31_33M_30x.bam   $V33 $B33 32000000 33000000 &
wait
run new32_33M_15x   data/giab_hg002_real_15x/hg002_chr21_31_33M_15x.bam   $V33 $B33 32000000 33000000 &
run re31_32M_30x    data/giab_hg002_real_30x/hg002_chr21_31_33M_30x.bam   $V33 $B33 31000000 32000000 &
wait
echo "ALL DONE $(date -Is)"
