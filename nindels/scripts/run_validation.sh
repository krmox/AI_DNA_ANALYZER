#!/usr/bin/env bash
# Independent validation, run ONCE after the configuration was frozen on HG004 chr20:
# HG001 chr19 30x. Primary: D = A+B (rule fixed before validation); secondary: A alone.
set -euo pipefail
cd "$(dirname "$0")/.."
export ENV_FILE=$PWD/scripts/env_validation.sh
source "$ENV_FILE"
mkdir -p "$RUNS/str"
s=$RUNS/str/${SAMPLE}_30x.tsv
python3 -I scripts/timed.py "$RUNS/str/time_30x.json" -- "$NINDELS_BIN" estimate-str \
    --bam "$DATA/${SAMPLE}_${CHROM}_30x.bam" --ref "$DATA/ref_${CHROM}.fa" --region "$CHROM" --threads 8 --out "$s"
run() { local tag=$1 c=$2 reps=$3; shift 3; TAG=$tag bash scripts/run_callers.sh "$c" 30 "$reps" "$@"; bash scripts/happy.sh "$tag"; }
run val_baseline_30x nindels 3
run val_A_30x nindels 3 --genotyper multi
run val_D_30x nindels 3 --genotyper multi --str-model "$s"
run val_bcftools_30x bcftools 1
echo VALIDATION_DONE
