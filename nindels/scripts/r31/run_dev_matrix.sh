#!/usr/bin/env bash
# Development matrix (threshold selection): HG004 chr20/chr21 x depth ladder x {A, AB, baseline}.
# Sequential on purpose (8 threads each); dev timing is not used for any claim.
set -u
cd "$(dirname "$0")"
for ds in "$@"; do
  for dk in d10 d15 d20 d25 d30 d40 d50 native; do
    for cfg in A AB baseline; do
      python3 run_one.py dev "$ds" "$dk" "$cfg" || echo "FAILED $ds $dk $cfg"
    done
  done
done
echo DEV_MATRIX_DONE "$@"
