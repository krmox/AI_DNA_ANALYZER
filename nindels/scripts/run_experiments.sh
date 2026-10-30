#!/usr/bin/env bash
# Controlled experiments on the development set (HG004 chr20). Each experiment
# changes exactly one component against the baseline. 30x: 3 timed runs, 67x: 1.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/env.sh
mkdir -p "$RUNS/str"
for d in 30 67; do
    s=$RUNS/str/HG004_${d}x.tsv
    [ -s "$s" ] || python3 -I scripts/timed.py "$RUNS/str/time_${d}x.json" -- "$NINDELS_BIN" estimate-str \
        --bam "$DATA/HG004_${CHROM}_${d}x.bam" --ref "$DATA/ref_${CHROM}.fa" --region "$CHROM" --threads 8 --out "$s"
done
run() {  # tag depth reps args...
    local tag=$1 d=$2 reps=$3; shift 3
    TAG=$tag bash scripts/run_callers.sh nindels "$d" "$reps" "$@"
    bash scripts/happy.sh "$tag"
}
for d in 30 67; do
    reps=$([ "$d" = 30 ] && echo 3 || echo 1)
    run expA_multi_${d}x "$d" "$reps" --genotyper multi
    run expB_str_${d}x "$d" "$reps" --str-model "$RUNS/str/HG004_${d}x.tsv"
    run expC_frag_${d}x "$d" "$reps" --fragments
done
echo EXPERIMENTS_DONE
