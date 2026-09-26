#!/bin/bash
# bench_full.sh LABEL MODE THREADS [COLD=0|1] [extra dnav2 args...]  -> full chr20, appends one row to $V2_PERF_CSV
set -euo pipefail
source "$(dirname "$0")/env.sh"
LABEL=$1; MODE=$2; THREADS=$3; COLD=${4:-0}; shift 4 || true
CSV=${V2_PERF_CSV:-/mnt/archive/AI_DNA_ANALYZER_v2/repo/research/V2_PERFORMANCE_RESULTS.csv}
[ -f "$CSV" ] || echo "experiment,date,label,mode,threads,rep,cache,region,wall_s,user_s,sys_s,maxrss_kb,n_loci,n_routed,n_cascade_calls,notes" > "$CSV"
OUT=$V2W/bench_full/$LABEL-$MODE-t$THREADS; mkdir -p "$OUT"
cache=warm
if [ "$COLD" = 1 ]; then python3 "$(dirname "$0")/drop_cache.py" "$BAM" >/dev/null; cache=cold; fi
/usr/bin/time -f "%e %U %S %M" -o "$OUT/time.txt" "$DNAV2" run --bam "$BAM" --ref "$REF" --truth "$TRUTH" --bed "$BED" \
    --out-dir "$OUT" --tag chr20_300x --mode "$MODE" --threads "$THREADS" "$@" > "$OUT/stdout.json"
read w u s m < "$OUT/time.txt"
nl=$(python3 -c "import json;d=json.load(open('$OUT/stdout.json'));print(d['n_loci'],d['n_routed_pb'],d['cascade']['n_calls'])")
# byte-identity to V1 is re-checked on every timed run (a fast wrong answer must not enter the table)
same=$(cmp -s "$OUT/ai_cascade_chr20_300x.vcf" "$V1RUN/ai_cascade_chr20_300x.vcf" && echo cascade_vcf_identical_to_V1 || echo CASCADE_VCF_DIFFERS)
echo "full_bench,$(date +%F),$LABEL,$MODE,$THREADS,1,$cache,chr20:full,$w,$u,$s,$m,${nl// /,},\"$* ; $same\"" >> "$CSV"
echo "$(date +%T) $LABEL $MODE t=$THREADS $cache wall=${w}s user=${u}s sys=${s}s rss=${m}KB $same"
