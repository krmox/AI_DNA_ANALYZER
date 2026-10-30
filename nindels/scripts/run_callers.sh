#!/usr/bin/env bash
# Runs nindels (untrained rule, every candidate reassessed) and bcftools on
# HG004 chr20 at each depth. REPS runs each for timing; outputs and timings go
# to $RUNS/<tag>/.
#   usage: run_callers.sh <caller: nindels|bcftools> <depth: 15|30|67> [reps] [extra nindels args...]
set -euo pipefail
source "$(dirname "$0")/env.sh"
caller=$1; depth=$2; reps=${3:-3}; shift 3 || true
tag=${TAG:-${caller}_${depth}x}
out=$RUNS/$tag
mkdir -p "$out"
bam=$DATA/HG004_${CHROM}_${depth}x.bam
ref=$DATA/ref_${CHROM}.fa
timed="python3 -I $(dirname "$0")/timed.py"

for r in $(seq 1 "$reps"); do
    sync
    case $caller in
        nindels)
            $timed "$out/time_$r.json" -- "$NINDELS_BIN" call --bam "$bam" --ref "$ref" --region "$CHROM" \
                --threads 8 --force-reassess --out-vcf "$out/calls.vcf" "$@" 2> "$out/log_$r.txt"
            ;;
        bcftools)
            $timed "$out/time_$r.json" -- "$TOOLS/bcftools" mpileup -f "$ref" -r "$CHROM" -a AD,DP -Ou "$bam" \
                "|" "$TOOLS/bcftools" call -mv -Ov -o "$out/calls.vcf" 2> "$out/log_$r.txt"
            ;;
    esac
done
"$TOOLS/bgzip" -f "$out/calls.vcf"
"$TOOLS/tabix" -f -p vcf "$out/calls.vcf.gz"
echo "DONE $tag"
