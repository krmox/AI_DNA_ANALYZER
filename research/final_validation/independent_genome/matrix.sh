#!/usr/bin/env bash
# usage: matrix.sh CONTIG   (resume-safe, sequential; 3 depths x 3 seeds x {r2,r11,bcftools})
c=$1; cd /home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER/research/final_validation/independent_genome
for t in 10 20 30; do for s in 101 102 103; do python3 run_campaign.py $c $t $s; done; done
echo MATRIX_DONE $c
