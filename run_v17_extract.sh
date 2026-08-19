#!/bin/bash
# bench_v17 extraction: 1 region x 3 depths, each run exactly once.
# Coordinates come from results/bench_v17/region_selection.json. Same extractor,
# same flags and the same resume guard as run_v14_extract.sh / run_v15_extract.sh,
# so the v17 test cells are produced by byte-identical code to every earlier cell.
#
# extract_bench_v12.py writes its .npz in place with np.savez_compressed, which
# is NOT atomic: an interrupted run can leave a truncated file that a plain
# `[ -s ... ]` resume guard would silently reuse. The extractor is frozen, so
# instead of touching it this script extracts to a scratch path and moves the
# result into place only after the extractor exits 0. A partial run therefore
# leaves nothing in cache/bench_v17/ to resume from, and the guard is safe.
set -u
D=data/giab_hg002_v17
mkdir -p cache/bench_v17 logs_v17

run () {
  local name=$1 contig=$2 start=$3 stop=$4 tag=$5
  local key="${name}_${tag}"
  if [ -s "cache/bench_v17/${key}.npz" ]; then
    echo "SKIP $key (already extracted) $(date -Is)"
    return 0
  fi
  echo "START $key $(date -Is)"
  local scratch="cache/bench_v17/${key}.partial.npz"
  rm -f "$scratch"
  if python3 extract_bench_v12.py \
      --fasta "data/reference/${contig}_full.fa" \
      --bam "$D/${name}_${tag}.bam" \
      --vcf "$D/${name}.vcf.gz" \
      --bed "$D/${name}_highconf.bed" \
      --contig "$contig" --region "$start" "$stop" \
      --features none \
      --out "$scratch" > "logs_v17/${key}.log" 2>&1; then
    mv "$scratch" "cache/bench_v17/${key}.npz"
    echo "DONE $key rc=0 $(date -Is)"
  else
    rm -f "$scratch"
    echo "FAILED $key $(date -Is)"
    return 1
  fi
}

for tag in full 30x 15x; do
  run chr17_segdup chr17 18000000 21000000 $tag &
done
wait
echo "ALL EXTRACTED $(date -Is)"
ls -l cache/bench_v17/
