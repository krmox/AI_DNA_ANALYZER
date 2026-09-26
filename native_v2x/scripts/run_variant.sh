#!/bin/bash
# run_variant.sh LABEL THREADS COLD(0|1) [dnav2 args...]  -> $V2W/runs/LABEL_tTHREADS_{cold|warm}/ ; prints wall/user/sys/rss
set -euo pipefail
source "$(dirname "$0")/env.sh"
LABEL=$1; THREADS=$2; COLD=$3; shift 3
BIN=${BIN:-$DNAV2}
OUT=$V2W/runs/${LABEL}_t${THREADS}_$([ "$COLD" = 1 ] && echo cold || echo warm); mkdir -p "$OUT"
[ "$COLD" = 1 ] && python3 "$(dirname "$0")/drop_cache.py" "$BAM" >/dev/null
/usr/bin/time -v -o "$OUT/time.txt" "$BIN" run --bam "$BAM" --ref "$REF" --truth "$TRUTH" --bed "$BED" \
   --out-dir "$OUT" --tag chr20_300x --threads "$THREADS" "$@" > "$OUT/stdout.json" 2> "$OUT/stderr.txt"
w=$(grep Elapsed "$OUT/time.txt" | awk '{print $NF}'); m=$(grep Maximum "$OUT/time.txt" | awk '{print $NF}'); c=$(grep "Percent of CPU" "$OUT/time.txt" | awk '{print $NF}')
echo "$LABEL threads=$THREADS cold=$COLD wall=$w rss_kb=$m cpu=$c out=$OUT"
