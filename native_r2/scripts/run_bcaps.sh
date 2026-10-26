#!/bin/bash
source "$(dirname "$0")/env.sh"; export BIN=/mnt/archive/AI_DNA_ANALYZER_v2x/build_x/dnav2; cd "$(dirname "$0")/.."
for m in Bcap24 Bcap48 Bcap64 Bcap96 Bcap128 Bcapall Bcapblkmed Bcapblkmax; do
  L=M_$m; scripts/run_variant.sh $L 8 0 --mode cascade --model models/$m.model
  scripts/run_happy.sh $L
done
echo BCAPS_DONE
