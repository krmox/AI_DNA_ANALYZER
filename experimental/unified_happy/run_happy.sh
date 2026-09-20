#!/bin/bash
# Step 3 (normalize) + Step 4 (hap.py) for ai_pb_only.vcf and ai_cascade.vcf.
# Mirrors experimental/head_to_head/COMMANDS.md section 6 exactly: same image,
# same digest, same flag structure (-f BED, -r reference), only query VCF differs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../.." && pwd)"  # project root (path-independent)
cd $ROOT/experimental/unified_happy

TRUTH=$ROOT/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz
BED=$ROOT/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed
REF=$ROOT/experimental/head_to_head/shared_ref/chr21_chrname.fa
IMAGE=quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0
DIGEST=sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0

echo "== normalize VCFs (sort, bgzip, tabix -- no bcftools norm, matching how DeepVariant/Clair3/GATK VCFs were consumed by hap.py per COMMANDS.md section 6) =="
for name in ai_pb_only ai_cascade; do
  bcftools sort -Oz -o ${name}.vcf.gz ${name}.vcf
  bcftools index -t ${name}.vcf.gz
done

echo "== disk before pull =="
df -h /home | tail -1

docker pull ${IMAGE}

echo "== disk after pull =="
df -h /home | tail -1

mkdir -p happy/ai_pb_only happy/ai_cascade

TRUTH_DIR=$(dirname "$TRUTH")
REF_DIR=$(dirname "$REF")

for name in ai_pb_only ai_cascade; do
  echo "== running hap.py on ${name} =="
  mkdir -p happy/${name}
  time docker run --rm \
    -v ${TRUTH_DIR}:/truth \
    -v $(pwd):/query \
    -v ${REF_DIR}:/ref \
    -v $(pwd)/happy/${name}:/out \
    ${IMAGE} \
    hap.py \
    /truth/hg002_chr21_32_44M.vcf.gz \
    /query/${name}.vcf.gz \
    -f /truth/hg002_chr21_32_44M_highconf.bed \
    -r /ref/chr21_chrname.fa \
    -o /out/happy
done

echo "== checksum outputs =="
sha256sum happy/ai_pb_only/*.csv happy/ai_cascade/*.csv > checksums_happy_outputs.txt
cat checksums_happy_outputs.txt

echo "== cleanup docker image per sequential disk protocol =="
docker rmi ${IMAGE}

echo "== disk after cleanup =="
df -h /home | tail -1
