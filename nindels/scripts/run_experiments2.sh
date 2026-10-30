#!/usr/bin/env bash
# C2 alone, and combinations, on the development set. D is chosen by the rule
# fixed beforehand: components that each improved F1 vs baseline at 30x and 67x.
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/env.sh
run() { local tag=$1 d=$2 reps=$3; shift 3; TAG=$tag bash scripts/run_callers.sh nindels "$d" "$reps" "$@"; bash scripts/happy.sh "$tag"; }
for d in 30 67; do
    reps=$([ "$d" = 30 ] && echo 3 || echo 1)
    run expC2_fragend_${d}x "$d" "$reps" --fragment-pick end
    run expAB_${d}x "$d" "$reps" --genotyper multi --str-model "$RUNS/str/HG004_${d}x.tsv"
    run expABC2_${d}x "$d" "$reps" --genotyper multi --str-model "$RUNS/str/HG004_${d}x.tsv" --fragment-pick end
done
echo EXPERIMENTS2_DONE
