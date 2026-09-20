#!/bin/bash
cd "$(dirname "$0")"; R=results
python3 equiv_windows.py > $R/equiv_windows.jsonl 2>&1
: > $R/equiv_reads_chunks.jsonl
for a in "hg005 1000000 1500000" "hg005 3500000 4000000" "hg002_15x 32000000 32500000" "hg002_15x 40000000 40500000" "chr20 46000000 46500000"; do python3 equiv_reads.py $a >> $R/equiv_reads_chunks.jsonl 2>> $R/equiv_reads_chunks.err; done
python3 profile_stages.py > $R/profile_stages_baseline.txt 2>&1
python3 cold_cache.py > $R/cold_cache.json 2>&1
./run_matrix.sh
./run_pipeline_matrix.sh
./run_e2e_hg005.sh 1 4
touch $R/partA.done
