#!/bin/bash
# bench_v14 extraction: 4 regions x 3 depths, each run exactly once
# (devlog 14 §3.5). Coordinates come from results/bench_v14/region_selection.json.
#
# --features none  : Stage-1 evidence only; this benchmark has no neural arm.
# no --legacy-pb   : the PB depth-clipping A/B is settled in devlog 13; every
#                    cell here uses the repaired PB and only the repaired PB.
set -u
D=data/giab_hg002_v14
mkdir -p cache/bench_v14 logs_v14

# name contig start stop depthtag
run () {
  local name=$1 contig=$2 start=$3 stop=$4 tag=$5
  local key="${name}_${tag}"
  # Resume guard. extract_bench_v12 writes its .npz once, atomically, after the
  # final chunk, so an existing file is necessarily a complete file and can be
  # kept. An interrupted run therefore resumes instead of redoing finished
  # cells, and no completed artefact is ever overwritten.
  if [ -s "cache/bench_v14/${key}.npz" ]; then
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
    --out "cache/bench_v14/${key}.npz" > "logs_v14/${key}.log" 2>&1
  echo "DONE $key rc=$? $(date -Is)"
}

for tag in full 30x 15x; do
  run chr20_neutral chr20  33000000  34000000 $tag &
  run chr19_gcrich  chr19   1000000   2000000 $tag &
  run chr4_atrich   chr4  101000000 102000000 $tag &
  run chr1_segdup   chr1  120000000 121000000 $tag &
  wait
  echo "BATCH $tag COMPLETE $(date -Is)"
done
echo "ALL EXTRACTED $(date -Is)"
ls -l cache/bench_v14/
