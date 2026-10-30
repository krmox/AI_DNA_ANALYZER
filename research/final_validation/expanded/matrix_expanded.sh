#!/usr/bin/env bash
# usage: matrix_expanded.sh SAMPLE CONTIG [dv]   (resume-safe, strictly sequential: timed runs must not overlap anything else)
#  default : NEW conditions 15x (seeds 101-103) and 50x (seed 101) x {r2,r11,bcftools}
#  dv      : DeepVariant at 30x seed 101 (BAM regenerated deterministically; sha256 compared with the first campaign afterwards)
export SAMPLE=$1; c=$2
[ "$SAMPLE" = HG001 ] && export BAMDIR=$HOME/varion_data/hg001_bam
cd /home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER/research/final_validation/expanded
if [ "${3:-}" = dv ]; then
  python3 run_expanded.py $c 30 101 --callers deepvariant
else
  for s in 101 102 103; do python3 run_expanded.py $c 15 $s; done
  python3 run_expanded.py $c 50 101
fi
echo MATRIX_DONE $SAMPLE $c ${3:-}
