#!/usr/bin/env bash
# Download HG003 GRCh38 300x NovoAlign BAM, chr8 only, in 16 disjoint coordinate slices (reads assigned by start pos), then concatenate.
set -euo pipefail
cd /mnt/archive/AI_DNA_ANALYZER_final_eval/bam
URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.300x.bam
LEN=145138636; N=16; STEP=$(( (LEN + N - 1) / N ))
mkdir -p parts
one() { i=$1; s=$(( i*STEP )); e=$(( (i+1)*STEP )); [ $e -gt $LEN ] && e=$LEN
  [ -f parts/p$(printf %02d $i).done ] && return 0
  for try in 1 2 3 4 5; do
    samtools view --no-PG -b -e "pos>=$((s+1)) && pos<=$e" -o parts/p$(printf %02d $i).bam -X $URL HG003.GRCh38.300x.bam.bai "chr8:$((s+1))-$e" && touch parts/p$(printf %02d $i).done && return 0
    sleep 10
  done; return 1; }
export -f one; export URL LEN STEP
seq 0 $((N-1)) | xargs -P 8 -I{} bash -c 'one {}'
samtools view --no-PG -H $URL > /dev/null 2>&1 || true
samtools cat -o HG003.GRCh38.300x_chr8.bam parts/p*.bam
samtools index HG003.GRCh38.300x_chr8.bam
samtools quickcheck -v HG003.GRCh38.300x_chr8.bam && echo QUICKCHECK_OK
sha256sum HG003.GRCh38.300x_chr8.bam HG003.GRCh38.300x_chr8.bam.bai > SHA256SUMS.local.txt
echo DOWNLOAD_DONE
