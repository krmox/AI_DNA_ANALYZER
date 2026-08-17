#!/bin/bash
# bench_v14 data acquisition. Every byte here comes from a public GIAB URL and
# the recipe is the same one recorded in data/giab_hg002_real_30M/PROVENANCE.txt,
# only with the contig and coordinates changed. Coordinates come from
# select_v14_regions.py (results/bench_v14/region_selection.json) and are not
# edited by hand.
set -u
D=data/giab_hg002_v14
BAM_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam
VCF_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
HC=$D/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

mkdir -p "$D" logs_v14

# name contig start stop
fetch () {
  local name=$1 contig=$2 start=$3 stop=$4
  local region="${contig}:$((start + 1))-${stop}"
  echo "START $name $region $(date -Is)"
  if [ ! -s "$D/${name}_full.bam" ]; then
    samtools view -b "$BAM_URL" "$region" -o "$D/${name}_full.bam" \
      && samtools index "$D/${name}_full.bam"
  fi
  if [ ! -s "$D/${name}.vcf.gz" ]; then
    bcftools view -r "$region" "$VCF_URL" -O z -o "$D/${name}.vcf.gz" \
      && bcftools index -t "$D/${name}.vcf.gz"
  fi
  if [ ! -s "$D/${name}_highconf.bed" ]; then
    awk -v c="$contig" -v s="$start" -v e="$stop" \
      'BEGIN{OFS="\t"} $1==c && $3>s && $2<e {print $1,$2,$3}' "$HC" > "$D/${name}_highconf.bed"
  fi
  echo "DONE $name rc=$? $(date -Is)"
}

fetch chr20_neutral chr20  33000000  34000000 > logs_v14/fetch_chr20.log 2>&1 &
fetch chr19_gcrich  chr19   1000000   2000000 > logs_v14/fetch_chr19.log 2>&1 &
fetch chr4_atrich   chr4  101000000 102000000 > logs_v14/fetch_chr4.log  2>&1 &
fetch chr1_segdup   chr1  120000000 121000000 > logs_v14/fetch_chr1.log  2>&1 &
wait
echo "ALL FETCHED $(date -Is)"
ls -l $D
