#!/bin/bash
cd "$(dirname "$0")"
./bench_full.sh full cascade 1 1
./bench_full.sh full cascade 2 0
./bench_full.sh full cascade 4 0
./bench_full.sh full cascade 8 0
./bench_full.sh full_d2 cascade 4 0 --decomp-threads 2
./bench_full.sh full pb_all 4 0
./bench_full.sh full pb_all 8 0
./bench_full.sh full pb_all 1 0
echo CHAIN_DONE
