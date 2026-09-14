#!/bin/bash
# bench_v18 data acquisition — CROSS-SAMPLE generalization test. Same region
# (chr13:70,000,000-73,000,000) as devlog 16's "representative" HG002 cell,
# same platform (Illumina 2x250) and same aligner (novoalign), but sample
# swapped to HG003 (trio father). This isolates the sample variable: nothing
# about the router, the thresholds or the region-selection rule changes.
# Every byte here comes from a public GIAB URL, same recipe as fetch_v17_data.sh.
set -u
D=data/giab_hg002_v18
BAM_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.GRCh38.2x250.bam
VCF_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
HC_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
HC=$D/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

mkdir -p "$D" logs_v18
if [ ! -s "$HC" ]; then
  curl -sS "$HC_URL" -o "$HC.part" && mv "$HC.part" "$HC"
fi

# name contig start stop   (start/stop are 0-based half-open, as selected)
fetch () {
  local name=$1 contig=$2 start=$3 stop=$4
  local region="${contig}:$((start + 1))-${stop}"
  echo "START $name $region $(date -Is)"
  if [ ! -s "$D/${name}_full.bam" ]; then
    samtools view -b "$BAM_URL" "$region" -o "$D/${name}_full.bam.part" \
      && mv "$D/${name}_full.bam.part" "$D/${name}_full.bam" \
      && samtools index "$D/${name}_full.bam"
  fi
  if [ ! -s "$D/${name}.vcf.gz" ]; then
    bcftools view -r "$region" "$VCF_URL" -O z -o "$D/${name}.vcf.gz.part" \
      && mv "$D/${name}.vcf.gz.part" "$D/${name}.vcf.gz" \
      && bcftools index -t "$D/${name}.vcf.gz"
  fi
  if [ ! -s "$D/${name}_highconf.bed" ]; then
    awk -v c="$contig" -v s="$start" -v e="$stop" \
      'BEGIN{OFS="\t"} $1==c && $3>s && $2<e {print $1,$2,$3}' "$HC" > "$D/${name}_highconf.bed"
  fi
  echo "DONE $name rc=$? $(date -Is)"
}

fetch chr13_crosssample chr13 70000000 73000000 > logs_v18/fetch_chr13.log 2>&1
echo "ALL FETCHED $(date -Is)"
ls -l $D
