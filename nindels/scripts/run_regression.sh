#!/usr/bin/env bash
# Controlled regression on HG004 chr20: Baseline, A, B, A+B with ONE binary,
# identical BAM / reference / truth / BED / hap.py settings; 30x and 67x, 3 timed
# runs each. STR models are re-estimated here from each BAM (no labels).
set -euo pipefail
cd "$(dirname "$0")/.."
export ENV_FILE=$PWD/scripts/env_regression.sh
source "$ENV_FILE"
mkdir -p "$RUNS/str"
{
  echo "date: $(date -Is)"
  echo "git_commit: $(git rev-parse HEAD)"
  echo "binary_sha256: $(sha256sum "$NINDELS_BIN" | cut -d' ' -f1)"
  for f in "$DATA"/HG004_chr20_30x.bam "$SRC_BAM" "$DATA/ref_chr20.fa" "$TRUTH" "$DATA/conf_chr20.bed"; do
    echo "input: $f size=$(stat -L -c %s "$f")"
  done
  echo "hap.py: $("$HAPPY_ENV/bin/python" "$HAPPY_ENV/bin/hap.py" --version 2>/dev/null | tail -1) engine=xcmp"
} > "$RUNS/manifest.txt"
for d in 30 67; do
    python3 -I scripts/timed.py "$RUNS/str/time_${d}x.json" -- "$NINDELS_BIN" estimate-str \
        --bam "$DATA/HG004_chr20_${d}x.bam" --ref "$DATA/ref_chr20.fa" --region chr20 --threads 8 \
        --out "$RUNS/str/HG004_${d}x.tsv" 2> "$RUNS/str/log_${d}x.txt"
done
run() { local tag=$1 d=$2; shift 2; TAG=$tag bash scripts/run_callers.sh nindels "$d" 3 "$@"; bash scripts/happy.sh "$tag"; }
for d in 30 67; do
    run baseline_${d}x "$d"
    run A_${d}x "$d" --genotyper multi
    run B_${d}x "$d" --str-model "$RUNS/str/HG004_${d}x.tsv"
    run AB_${d}x "$d" --genotyper multi --str-model "$RUNS/str/HG004_${d}x.tsv"
done
echo REGRESSION_DONE
