#!/usr/bin/env bash
# Scores $RUNS/<tag>/calls.vcf.gz against the GIAB v4.2.1 truth of $ENV_FILE (default HG004 chr20) with
# hap.py (xcmp engine; rtg/vcfeval is not installed), stratified by repeat context.
#   usage: happy.sh <tag>
set -euo pipefail
source "${ENV_FILE:-$(dirname "$0")/env.sh}"
tag=$1
out=$RUNS/$tag
export PATH=$HAPPY_ENV/bin:$PATH
export HGREF=$DATA/ref_${CHROM}.fa
"$HAPPY_ENV/bin/python" "$HAPPY_ENV/bin/hap.py" "$TRUTH" "$out/calls.vcf.gz" \
    -r "$DATA/ref_${CHROM}.fa" -f "$DATA/conf_${CHROM}.bed" -l "$CHROM" \
    --stratification "$DATA/strata.tsv" --threads 8 \
    -o "$out/happy" > "$out/happy.log" 2>&1
grep -E '^(INDEL|SNP),(ALL|PASS),' "$out/happy.summary.csv" | cut -d, -f1-4,6,7,8,11,12,13,14
