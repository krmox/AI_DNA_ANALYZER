#!/usr/bin/env bash
# Builds the benchmark inputs under /mnt/hdd/nindels/data from files already on
# the HDD: depth subsets of HG004 chr20 (fixed seed), the chr20 reference, the
# chr20 confident regions and the stratification BEDs used for the breakdown.
set -euo pipefail
source "$(dirname "$0")/env.sh"
mkdir -p "$DATA"

# Depth subsets. samtools -s SEED.FRACTION keeps read pairs together.
src_depth=66.7
for d in 30 15; do
    out=$DATA/HG004_${CHROM}_${d}x.bam
    [ -s "$out.bai" ] && continue
    frac=$(python3 -c "print(f'{$d/$src_depth:.4f}'.split('.')[1])")
    "$TOOLS/samtools" view -@4 -b -s "21.$frac" -o "$out" "$SRC_BAM"
    "$TOOLS/samtools" index "$out"
done
ln -sf "$SRC_BAM" "$DATA/HG004_${CHROM}_67x.bam"
ln -sf "$SRC_BAM.bai" "$DATA/HG004_${CHROM}_67x.bam.bai"

# Reference restricted to the contig (smaller page-cache footprint, same coordinates).
[ -s "$DATA/ref_${CHROM}.fa.fai" ] || {
    "$TOOLS/samtools" faidx "$REF" "$CHROM" > "$DATA/ref_${CHROM}.fa"
    "$TOOLS/samtools" faidx "$DATA/ref_${CHROM}.fa"
}

awk -v c="$CHROM" '$1 == c' "$CONF_BED" > "$DATA/conf_${CHROM}.bed"

# Stratifications for hap.py (name <TAB> bed), GIAB v3.x GRCh38 set.
L=$STRAT_DIR/LowComplexity
cat > "$DATA/strata.tsv" <<EOF
homopolymer_ge7	$L/GRCh38_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz
tandem_repeat_le50	$L/GRCh38_AllTandemRepeats_le50bp_slop5.bed.gz
tandem_repeat_51to200	$L/GRCh38_AllTandemRepeats_51to200bp_slop5.bed.gz
not_in_TR_or_HP	$L/GRCh38_notinAllTandemRepeatsandHomopolymers_slop5.bed.gz
EOF
for f in $(cut -f2 "$DATA/strata.tsv"); do [ -s "$f" ] || { echo "missing $f"; exit 1; }; done

for f in "$DATA"/*.bam; do echo "$(basename "$f") $("$TOOLS/samtools" idxstats "$f" | awk -v c="$CHROM" '$1==c{print $3}') reads"; done
echo PREPARE_DONE
