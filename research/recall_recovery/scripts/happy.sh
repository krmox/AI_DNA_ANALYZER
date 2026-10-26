#!/usr/bin/env bash
# usage: happy.sh TRUTHSET TAG...   TRUTHSET in {hg003_chr8, hg002_chr20_odd, hg004}; same hap.py image/options as Stage 1
set -uo pipefail
S=/mnt/archive/AI_DNA_ANALYZER_stage2; B=/mnt/archive/AI_DNA_ANALYZER_benchmark; E=/mnt/archive/AI_DNA_ANALYZER_final_eval
IMG=quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0
set_=$1; shift
case $set_ in
  hg003_chr8) TR=$E/truth; TF=HG003_GRCh38_1_22_v4.2.1_benchmark; L=chr8; extra=();;
  hg002_chr20_odd) TR=$B/truth_HG002_GRCh38_v4.2.1; TF=HG002_GRCh38_1_22_v4.2.1_benchmark; L=chr20; extra=(-T /w/holdout_odd_blocks.bed);;
  hg004_*) TR=$S/truth; TF=HG004_GRCh38_1_22_v4.2.1_benchmark; L=${set_#hg004_}; extra=();;
esac
mkdir -p $S/happy/query $S/happy/out
for tag in "$@"; do
  if [[ $tag == bcftools_* ]]; then cp $S/runs/$tag/bc.vcf.gz $S/happy/query/$tag.vcf.gz; cp $S/runs/$tag/bc.vcf.gz.tbi $S/happy/query/$tag.vcf.gz.tbi
  else bcftools view -Oz -o $S/happy/query/$tag.vcf.gz $S/runs/$tag/ai_cascade_$tag.vcf && bcftools index -f -t $S/happy/query/$tag.vcf.gz; fi
  podman --root /dev/shm/podman_root --storage-driver vfs run --rm --security-opt label=disable -v $TR:/truth -v $B/reference:/ref -v $S/happy/query:/query -v $S/happy/out:/out -v /mnt/archive/AI_DNA_ANALYZER_depth30:/w $IMG \
    /usr/local/bin/hap.py /truth/$TF.vcf.gz /query/$tag.vcf.gz -f /truth/${TF}_noinconsistent.bed -r /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna -o /out/$tag -l $L --threads 4 "${extra[@]}" > $S/happy/out/$tag.log 2>&1
  echo "$tag happy rc=$?" >> $S/logs/happy.log
done
