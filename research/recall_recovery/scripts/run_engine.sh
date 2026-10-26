#!/usr/bin/env bash
# usage: run_engine.sh TAG BAM CONTIG LEN TRUTH BED SAMPLE MODEL [extra flags...]   (outputs in stage2/runs/TAG)
set -uo pipefail
S=/mnt/archive/AI_DNA_ANALYZER_stage2; B=/mnt/archive/AI_DNA_ANALYZER_benchmark
REF=$B/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
tag=$1; bam=$2; ctg=$3; len=$4; truth=$5; bed=$6; sample=$7; model=$8; shift 8
out=$S/runs/$tag; mkdir -p $out
/usr/bin/time -v -o $out/time.txt $S/build/dnav2 run --bam $bam --ref $REF --truth $truth --bed $bed --out-dir $out \
  --contig $ctg --contig-len $len --sample $sample --tag $tag --threads 8 --mode cascade --model $model "$@" > $out/stdout.json 2> $out/stderr.txt
echo "$tag rc=$? flags: $*" >> $S/logs/engine.log
