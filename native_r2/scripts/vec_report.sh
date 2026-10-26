#!/bin/bash
# GCC auto-vectorization report for every V2 translation unit, generic x86-64 (SSE2) vs AVX2+FMA.
set -euo pipefail
cd "$(dirname "$0")/.."
HTS_INC=$(python3 -c 'import pysam; print(pysam.get_include()[1])')
BASE=(g++ -O3 -std=c++20 -ffp-contract=off -fno-fast-math -frounding-math -Iinclude -isystem "$HTS_INC" -fopt-info-vec-optimized -c)
mkdir -p logs
for f in extract engine frame capi; do
  "${BASE[@]}" src/$f.cpp -o /dev/null 2> logs/vec_report_${f}_generic.txt || true
  "${BASE[@]}" -mavx2 -mfma src/$f.cpp -o /dev/null 2> logs/vec_report_${f}_avx2.txt || true
done
for t in generic avx2; do
  echo "== $t"
  cat logs/vec_report_*_${t}.txt | grep -E "optimized: (loop|basic block)" | sed -E 's/^([^:]*):([0-9]+):[0-9]+: optimized: (.*)$/\1:\2 \3/' | sort | uniq -c | sort -rn | head -40
done
