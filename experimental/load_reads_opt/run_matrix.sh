#!/bin/bash
# load_reads candidate matrix: identical BAM/region/chunking/hardware; 3 reps each, fresh process per rep.
cd "$(dirname "$0")"
OUT=results/matrix_load_reads.jsonl
run () {
  local cand=$1 mode=$2 workers=$3 ds=$4 s=$5 e=$6 chunk=$7; shift 7
  for rep in 1 2 3; do
    env BENCH_MODE=$mode "$@" python3 bench_candidates.py $cand $ds $s $e $chunk $workers | python3 -c "
import sys,json; d=json.loads(sys.stdin.read()); d['rep']=$rep; d['label']='$cand/$mode/w$workers/$*'; print(json.dumps(d))" >> $OUT
  done
}
: > $OUT
H5="hg005 1000000 1400000 50000"
run baseline serial 1 $H5
run microopt serial 1 $H5
run baseline threads 4 $H5
run baseline procs 4 $H5
run baseline procs 8 $H5
run native_origscan serial 1 $H5
run native serial 1 $H5
run native serial 1 $H5 AI_DNA_ANALYZER_READ_THREADS=4
run native threads 4 $H5
run native procs 4 $H5
run native procs 8 $H5
H2="hg002_15x 32000000 32400000 50000"
run baseline serial 1 $H2
run native serial 1 $H2
run native procs 4 $H2
touch results/matrix.done
