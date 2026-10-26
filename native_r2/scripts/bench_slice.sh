#!/bin/bash
# bench_slice.sh BINARY LABEL MODE THREADS REPS [extra args...]  -> appends rows to $V2_PERF_CSV
# Region: 8192 contiguous windows starting at window 400000 (chr20:31,818,624-32,398,272), warm page cache
# unless COLD=1 (then the BAM is evicted from the page cache before every rep).
set -euo pipefail
source "$(dirname "$0")/env.sh"
BIN=$1; LABEL=$2; MODE=$3; THREADS=$4; REPS=$5; shift 5
CSV=${V2_PERF_CSV:-/mnt/archive/AI_DNA_ANALYZER_v2/repo/research/V2_PERFORMANCE_RESULTS.csv}
[ -f "$CSV" ] || echo "experiment,date,label,mode,threads,rep,cache,region,wall_s,user_s,sys_s,maxrss_kb,n_loci,n_routed,n_cascade_calls,notes" > "$CSV"
OUT=$V2W/bench/$LABEL-$MODE-t$THREADS; mkdir -p "$OUT"
for r in $(seq 1 "$REPS"); do
  cache=warm
  if [ "${COLD:-0}" = 1 ]; then python3 "$(dirname "$0")/drop_cache.py" "$BAM" >/dev/null; cache=cold; fi
  /usr/bin/time -f "%e %U %S %M" -o "$OUT/time_$r.txt" "$BIN" run --bam "$BAM" --ref "$REF" --truth "$TRUTH" --bed "$BED" \
      --out-dir "$OUT" --tag chr20_300x --mode "$MODE" --skip-windows 400000 --limit-windows 8192 --threads "$THREADS" "$@" > "$OUT/stdout_$r.json"
  read w u s m < "$OUT/time_$r.txt"
  nl=$(python3 -c "import json;d=json.load(open('$OUT/stdout_$r.json'));print(d['n_loci'],d['n_routed_pb'],d['cascade']['n_calls'])")
  echo "slice_bench,$(date +%F),$LABEL,$MODE,$THREADS,$r,$cache,chr20:31818624-32398272,$w,$u,$s,$m,${nl// /,},\"$*\"" >> "$CSV"
  echo "$LABEL $MODE t=$THREADS rep=$r $cache wall=${w}s user=${u}s rss=${m}KB"
done
