#!/bin/bash
# bench_v15 extraction: 4 regions x 3 depths, each run exactly once
# (devlog 15 §7.7). Coordinates come from results/bench_v15/region_selection.json.
# Same extractor, same flags and the same resume guard as run_v14_extract.sh, so
# the test cells are produced by byte-identical code to the validation cells.
set -u
D=data/giab_hg002_v15
mkdir -p cache/bench_v15 logs_v15

run () {
  local name=$1 contig=$2 start=$3 stop=$4 tag=$5
  local key="${name}_${tag}"
  if [ -s "cache/bench_v15/${key}.npz" ]; then
    echo "SKIP $key (already extracted) $(date -Is)"
    return 0
  fi
  echo "START $key $(date -Is)"
  python3 extract_bench_v12.py \
    --fasta "data/reference/${contig}_full.fa" \
    --bam "$D/${name}_${tag}.bam" \
    --vcf "$D/${name}.vcf.gz" \
    --bed "$D/${name}_highconf.bed" \
    --contig "$contig" --region "$start" "$stop" \
    --features none \
    --out "cache/bench_v15/${key}.npz" > "logs_v15/${key}.log" 2>&1
  echo "DONE $key rc=$? $(date -Is)"
}

for tag in full 30x 15x; do
  run chr16_segdup  chr16 29000000 30000000 $tag &
  run chr15_segdup  chr15 30000000 31000000 $tag &
  run chr7_segdup   chr7  57000000 58000000 $tag &
  run chr14_neutral chr14 53000000 54000000 $tag &
  wait
  echo "BATCH $tag COMPLETE $(date -Is)"
done
echo "ALL EXTRACTED $(date -Is)"
ls -l cache/bench_v15/
