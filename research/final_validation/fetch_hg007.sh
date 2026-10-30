#!/usr/bin/env bash
# Fresh independent data: HG007 (NA24695) HiSeq 2x148 novoalign ~100x GRCh38, streamed per contig (same recipe as research/stage4/scripts/fetch_hg006.sh).
set -uo pipefail
export PATH=$HOME/varion_mamba/envs/vt/bin:$PATH
D=/mnt/hdd/AI_DNA_ANALYZER_hg007; U=$(cat $D/logs/source_url.txt)
for c in "$@"; do
  [ -s $D/bam/HG007_${c}_native.bam ] && continue
  for i in 1 2 3; do samtools view -X -b -o $D/bam/HG007_${c}_native.bam.part $U $D/bam/HG007.100x.bam.bai $c && break; echo "retry $i $c"; done
  mv $D/bam/HG007_${c}_native.bam.part $D/bam/HG007_${c}_native.bam && samtools index $D/bam/HG007_${c}_native.bam \
    && echo "$c $(stat -c %s $D/bam/HG007_${c}_native.bam) $(sha256sum $D/bam/HG007_${c}_native.bam | cut -d' ' -f1)" >> $D/logs/fetch.log
done
echo FETCH_DONE >> $D/logs/fetch.log
