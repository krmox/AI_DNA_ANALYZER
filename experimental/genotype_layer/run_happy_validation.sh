#!/bin/bash
# hap.py evaluation for methods A/B on the TRUE hold-out validation region
# (chr21:32-44Mb), mirroring run_happy.sh exactly (same image/digest, same
# flag structure), only truth/BED/query differ (validation region instead of
# the dev region).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../.." && pwd)"  # project root (path-independent)
cd $ROOT/experimental/genotype_layer

TRUTH=$ROOT/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz
BED=$ROOT/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed
REF=$ROOT/experimental/head_to_head/shared_ref/chr21_chrname.fa
IMAGE=quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0

mkdir -p happy/validation_A happy/validation_B

TRUTH_DIR=$(dirname "$TRUTH")
REF_DIR=$(dirname "$REF")

for name in A B; do
  echo "== running hap.py on validation_${name} =="
  time docker run --rm \
    -v ${TRUTH_DIR}:/truth \
    -v $(pwd)/vcfs:/query \
    -v ${REF_DIR}:/ref \
    -v $(pwd)/happy/validation_${name}:/out \
    ${IMAGE} \
    hap.py \
    /truth/hg002_chr21_32_44M.vcf.gz \
    /query/validation_${name}.vcf.gz \
    -f /truth/hg002_chr21_32_44M_highconf.bed \
    -r /ref/chr21_chrname.fa \
    -o /out/happy
done

echo "== checksum outputs =="
sha256sum happy/validation_A/*.csv happy/validation_B/*.csv > checksums_happy_outputs_validation.txt
cat checksums_happy_outputs_validation.txt
