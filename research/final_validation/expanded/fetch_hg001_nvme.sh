#!/usr/bin/env bash
# Re-fetch HG001 native chr BAMs (same remote slices as the first campaign) onto the NVMe; verify against recorded sha256.
set -uo pipefail
export PATH=$HOME/varion_mamba/envs/vt/bin:$PATH
D=/mnt/hdd/AI_DNA_ANALYZER_hg001; O=$HOME/varion_data/hg001_bam; U=$(cat $D/logs/source_url.txt)
L=/home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER/research/final_validation/independent_genome/logs/fetch_sha256.log
for c in "$@"; do
  [ -s $O/HG001_${c}_native.bam ] || {
  for i in 1 2 3; do samtools view -X -b -o $O/HG001_${c}_native.bam.part $U $D/bam/HG001.300x.bam.bai $c && break; echo "retry $i $c"; done
  mv $O/HG001_${c}_native.bam.part $O/HG001_${c}_native.bam && samtools index $O/HG001_${c}_native.bam; }
  got=$(sha256sum $O/HG001_${c}_native.bam | cut -d' ' -f1); want=$(grep "^$c " $L | cut -d' ' -f3)
  [ "$got" = "$want" ] && echo "$c SHA_MATCH $got" || echo "$c SHA_MISMATCH got=$got want=$want"
done
echo FETCH_DONE
