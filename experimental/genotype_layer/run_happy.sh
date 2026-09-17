#!/bin/bash
# hap.py evaluation for methods A/B/C on the dev region (chr21:30.0-30.47Mb),
# mirroring experimental/unified_happy/run_happy.sh exactly (same image/digest,
# same flag structure), only the truth/BED/query differ (dev region instead of
# the held-out 32-44M region).
set -euo pipefail
cd /home/mark/Documents/Projects/AI_DNA_ANALYZER/experimental/genotype_layer

TRUTH=/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/giab_hg002_real_30M/hg002_chr21_30M.vcf.gz
BED=/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/giab_hg002_real_30M/hg002_chr21_30M_highconf.bed
REF=/home/mark/Documents/Projects/AI_DNA_ANALYZER/experimental/head_to_head/shared_ref/chr21_chrname.fa
IMAGE=quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0

mkdir -p happy/A happy/B happy/C

TRUTH_DIR=$(dirname "$TRUTH")
REF_DIR=$(dirname "$REF")

for name in A B C; do
  echo "== running hap.py on ${name} =="
  time docker run --rm \
    -v ${TRUTH_DIR}:/truth \
    -v $(pwd)/vcfs:/query \
    -v ${REF_DIR}:/ref \
    -v $(pwd)/happy/${name}:/out \
    ${IMAGE} \
    hap.py \
    /truth/hg002_chr21_30M.vcf.gz \
    /query/${name}.vcf.gz \
    -f /truth/hg002_chr21_30M_highconf.bed \
    -r /ref/chr21_chrname.fa \
    -o /out/happy
done

echo "== checksum outputs =="
sha256sum happy/A/*.csv happy/B/*.csv happy/C/*.csv > checksums_happy_outputs.txt
cat checksums_happy_outputs.txt
