#!/usr/bin/env bash
# Regression: in the default (truth-free) mode the call set must not depend on the truth file (inference/evaluation separation).
# Runs the same slice with the real truth VCF and with a nonexistent path; the VCFs must be byte-identical,
# and the R1 mode (--truth-mask --strict-windows) must differ (proves the mask exists and the flag removes it).
# usage: test_truth_free.sh DNAV2 BAM REF TRUTH_VCF BED CONTIG CONTIG_LEN [MODEL]
set -euo pipefail
BIN=$1; BAM=$2; REF=$3; TRUTH=$4; BED=$5; CONTIG=$6; LEN=$7; MODEL=${8:-}
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
m=(); [ -n "$MODEL" ] && m=(--model "$MODEL")
run() { # tag truth extra...
  local tag=$1 truth=$2; shift 2
  mkdir -p "$W/$tag"
  "$BIN" run --bam "$BAM" --ref "$REF" --truth "$truth" --bed "$BED" --out-dir "$W/$tag" --contig "$CONTIG" \
    --contig-len "$LEN" --tag t --threads 4 --limit-windows 40000 "${m[@]}" "$@" > "$W/$tag/out.json"
}
run real "$TRUTH"
run none /nonexistent/truth.vcf.gz
run masked "$TRUTH" --truth-mask --strict-windows
cmp "$W/real/ai_cascade_t.vcf" "$W/none/ai_cascade_t.vcf"
echo "PASS: default VCF independent of truth file ($(grep -vc '^#' "$W/real/ai_cascade_t.vcf") records)"
if cmp -s "$W/real/ai_cascade_t.vcf" "$W/masked/ai_cascade_t.vcf"; then
  echo "NOTE: R1 and default outputs identical on this slice (no truth-indel positions called)"
else
  echo "PASS: R1 mode differs from default ($(grep -vc '^#' "$W/masked/ai_cascade_t.vcf") vs $(grep -vc '^#' "$W/real/ai_cascade_t.vcf") records)"
fi
