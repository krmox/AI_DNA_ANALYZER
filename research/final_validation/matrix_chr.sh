#!/usr/bin/env bash
# usage: matrix_chr.sh CONTIG   (runs sequentially; resume-safe)
c=$1; cd /home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER/research/final_validation
for t in 10 20 30; do for s in 101 102 103; do python3 run_hg007.py $c $t $s; done; done
python3 run_hg007.py $c native 0
echo MATRIX_DONE $c
