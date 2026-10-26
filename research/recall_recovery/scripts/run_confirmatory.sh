#!/usr/bin/env bash
# Confirmatory accuracy runs per STAGE2_PREREG.md: 2 contigs x 3 depths x 4 variants, then hap.py. Sequential, no timing claims.
set -uo pipefail
S=/mnt/archive/AI_DNA_ANALYZER_stage2; R=/mnt/archive/AI_DNA_ANALYZER_depth_robust
T=$S/truth/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz; BD=$S/truth/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
declare -A LEN=([chr21]=46709983 [chr20]=64444167)
for c in chr21 chr20; do for d in f053_s21 full f010_s21; do
  bam=$S/bam/HG004_${c}_${d}.bam
  [ $d = full ] && bam=$S/bam/HG004_${c}_full.bam
  for v in B0: H1:--edge-windows H2:--truth-free H1H2:"--edge-windows --truth-free"; do
    n=${v%%:*}; f=${v#*:}
    $S/run_engine.sh c_${c}_${d}_$n $bam $c ${LEN[$c]} $T $BD HG004 $R/R1.model $f
  done
  $S/happy.sh hg004_$c c_${c}_${d}_B0 c_${c}_${d}_H1 c_${c}_${d}_H2 c_${c}_${d}_H1H2
done; done
echo CONFIRM_DONE >> $S/logs/happy.log
