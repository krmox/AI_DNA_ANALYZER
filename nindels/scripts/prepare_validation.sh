#!/usr/bin/env bash
# HG001 chr19 at 30x (seed 21), chr19 reference, confident regions, strata.
set -euo pipefail
source /home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER/.worktrees/nindels-test/nindels/scripts/env_validation.sh
mkdir -p "$DATA" "$RUNS"
out=$DATA/${SAMPLE}_${CHROM}_30x.bam
if [ ! -s "$out.bai" ]; then
    # 100,572,395 reads x 148 bp / 58,617,616 bp = 254x -> 30/254 = 0.118
    "$TOOLS/samtools" view -@4 -b -s 21.118 -o "$out" "$SRC_BAM"
    "$TOOLS/samtools" index "$out"
fi
[ -s "$DATA/ref_${CHROM}.fa.fai" ] || {
    "$TOOLS/samtools" faidx "$REF" "$CHROM" > "$DATA/ref_${CHROM}.fa"
    "$TOOLS/samtools" faidx "$DATA/ref_${CHROM}.fa"
}
awk -v c="$CHROM" '$1 == c' "$CONF_BED" > "$DATA/conf_${CHROM}.bed"
cp /mnt/hdd/nindels/data/strata.tsv "$DATA/strata.tsv"
echo "$(basename "$out") $("$TOOLS/samtools" idxstats "$out" | awk -v c="$CHROM" '$1==c{print $3}') reads"
echo PREPARE_VALIDATION_DONE
