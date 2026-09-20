#!/bin/bash
# hap.py (same image/digest as every prior project run) on the chr20 holdout: PB-only+MethodC and cascade+MethodC.
# Same truth VCF, same confident-region BED, same reference, same default settings for both arms.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../.." && pwd)"  # project root (path-independent)
HERE=$ROOT/experimental/chr20_validation
D=$ROOT/data/giab_hg002_chr20
IMAGE=quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0
cd "$HERE"
for name in ai_pb_only_chr20 ai_cascade_chr20; do
  bcftools view -Oz -o ${name}.vcf.gz ${name}.vcf && bcftools index -f -t ${name}.vcf.gz
  mkdir -p happy/${name}
  echo "== hap.py ${name} =="
  /usr/bin/time -v docker run --rm \
    -v ${D}:/truth -v ${D}/shared_ref:/ref -v ${HERE}:/query -v ${HERE}/happy/${name}:/out \
    $IMAGE hap.py /truth/chr20.vcf.gz /query/${name}.vcf.gz \
    -f /truth/chr20_highconf.bed -r /ref/chr20_chrname.fa -o /out/happy 2> happy/${name}.time
done
sha256sum happy/ai_pb_only_chr20/*.csv happy/ai_cascade_chr20/*.csv > checksums_happy_outputs_chr20.txt
