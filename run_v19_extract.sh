#!/bin/bash
# bench_v19 extraction (Validation Round 2, HG004): 3 regions x 3 depths, each
# run exactly once. Same extractor, same flags, same atomic-scratch resume
# guard as run_v14/16/17/18_extract.sh.
set -u
D=data/giab_hg004_v19
mkdir -p cache/bench_v19 logs_v19

run () {
  local name=$1 contig=$2 start=$3 stop=$4 tag=$5
  local key="${name}_${tag}"
  if [ -s "cache/bench_v19/${key}.npz" ]; then
    echo "SKIP $key (already extracted) $(date -Is)"
    return 0
  fi
  echo "START $key $(date -Is)"
  local scratch="cache/bench_v19/${key}.partial.npz"
  rm -f "$scratch"
  if python3 extract_bench_v12.py \
      --fasta "data/reference/${contig}_full.fa" \
      --bam "$D/${name}_${tag}.bam" \
      --vcf "$D/${name}.vcf.gz" \
      --bed "$D/${name}_highconf.bed" \
      --contig "$contig" --region "$start" "$stop" \
      --features none \
      --out "$scratch" > "logs_v19/${key}.log" 2>&1; then
    mv "$scratch" "cache/bench_v19/${key}.npz"
    echo "DONE $key rc=0 $(date -Is)"
  else
    rm -f "$scratch"
    echo "FAILED $key $(date -Is)"
    return 1
  fi
}

for tag in full 30x 15x; do
  run chr2_ordinary chr2 210000000 213000000 $tag &
  run chr3_difficult chr3 75000000 78000000 $tag &
  run chr5_stress chr5 27000000 30000000 $tag &
done
wait
echo "ALL EXTRACTED $(date -Is)"
ls -l cache/bench_v19/
