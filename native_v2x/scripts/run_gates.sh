#!/bin/bash
# Full correctness gate sequence (needs the archive).  Appends every check to research/V2_CORRECTNESS_RESULTS.csv.
set -euo pipefail
cd "$(dirname "$0")/.."; source scripts/env.sh
export V2_CORRECTNESS_CSV=${V2_CORRECTNESS_CSV:-$PWD/../research/V2_CORRECTNESS_RESULTS.csv}
python3 tests/test_kernels.py
python3 tests/test_blake2b.py
python3 tests/test_reconstructed_v1_vcf.py
mkdir -p $V2W/setup && $DNAV2 run --bam $BAM --ref $REF --truth $TRUTH --bed $BED --out-dir $V2W/setup --setup-only && python3 tests/test_frame.py
[ -f $V2W/gate/v1_gate.npz ] || python3 scripts/make_gate_v1.py
G=$V2W/gate/run_gate; mkdir -p $G
$DNAV2 run --bam $BAM --ref $REF --truth $TRUTH --bed $BED --out-dir $G --tag chr20_300x --loci-file $V2W/gate/gate_positions.txt --dump-loci $G/dump.tsv --threads 4 >/dev/null
python3 tests/gate_compare.py $G/dump.tsv gate_3491_loci
V2_LABEL=gate_vcf python3 tests/vcf_compare.py $G/ai_cascade_chr20_300x.vcf $V1RUN/ai_cascade_chr20_300x.vcf $V2W/gate/gate_positions.txt
echo "gates done (full-chromosome comparison: python3 tests/full_compare.py <pb_all run dir>)"
