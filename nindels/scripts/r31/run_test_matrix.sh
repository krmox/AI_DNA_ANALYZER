#!/usr/bin/env bash
# Final test matrix: {baseline, A, AB, R31} x depth ladder on held-out samples/chromosomes.
# Run ONCE, after varion/thresholds.json has been frozen and committed. Guards below refuse otherwise.
#   usage: run_test_matrix.sh PHASE DATASET [DATASET ...]      (PHASE = test | seen)
set -u
cd "$(dirname "$0")"
REPO=$(git rev-parse --show-toplevel)
if grep -q unselected "$REPO/varion/thresholds.json"; then echo "thresholds not frozen"; exit 3; fi
if [ -n "$(git status --porcelain -- "$REPO/varion/thresholds.json")" ]; then echo "thresholds.json has uncommitted changes"; exit 3; fi
phase=$1; shift
for ds in "$@"; do
  for dk in $(python3 -c "import json,sys;print(' '.join(json.load(open('/home/kramkost/r31/$ds/ladder.json'))['ladder']))"); do
    reps=1; [ "$dk" = d30 ] && reps=3
    for cfg in baseline A AB R31; do
      python3 run_one.py "$phase" "$ds" "$dk" "$cfg" --reps $reps || echo "FAILED $ds $dk $cfg"
    done
  done
done
echo TEST_MATRIX_DONE "$phase" "$@"
