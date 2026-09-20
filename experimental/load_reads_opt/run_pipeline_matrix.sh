#!/bin/bash
# E2E chunk-pipeline modes, HG005 chr1:1.0-1.8 Mb, 8 x 100 kb chunks. Same day, same machine.
cd "$(dirname "$0")"
OUT=results/matrix_pipeline.jsonl
run () { local reps=$1; shift; for rep in $(seq $reps); do python3 bench_pipeline.py "$@" | python3 -c "
import sys,json; d=json.loads(sys.stdin.read()); d['rep']=$rep; print(json.dumps(d))" >> $OUT; done; }
: > $OUT
S=1000000; E=1800000
run 2 serial python $S $E 1
run 2 async  python $S $E 1
run 1 threads python $S $E 4
run 3 serial native $S $E 1
run 3 async  native $S $E 1
run 3 threads native $S $E 4
run 3 procs  native $S $E 4
run 3 procs  native $S $E 6
run 2 procs  python $S $E 6
touch results/matrix_pipeline.done
