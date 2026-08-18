#!/bin/bash
# bench_v16 reference acquisition. Same source and release as every existing
# data/reference/chr*_full.fa: Ensembl release-110 GRCh38 per-chromosome FASTA.
# The contig comes from select_v16_region.py Rule C16
# (results/bench_v16/region_selection.json) and is not edited by hand.
set -u
mkdir -p data/reference logs_v16
BASE=https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna

for c in 13; do
  out=data/reference/chr${c}_full.fa
  if [ -s "$out" ]; then echo "have $out"; continue; fi
  echo "START chr${c} $(date -Is)"
  curl -A "Mozilla/5.0" --retry 5 --retry-connrefused -sSL \
    "$BASE/Homo_sapiens.GRCh38.dna.chromosome.${c}.fa.gz" -o "${out}.gz" \
    && gunzip -f "${out}.gz" && samtools faidx "$out"
  echo "DONE chr${c} rc=$? $(date -Is)"
done
echo "ALL REFERENCES $(date -Is)"
ls -l data/reference/chr13_full.fa*
