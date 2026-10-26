#!/usr/bin/env bash
# R2 development data: (1) extra downsamples of existing HG004 BAMs, (2) feature dumps from the FROZEN R1.1 binary
# (truth-free frame, so candidates at indel sites are included; labels from truth only fill the label column).
# Blind sets (chr18 of HG002/3/4, HG005) are NOT touched here.
set -uo pipefail
R2=/mnt/archive/AI_DNA_ANALYZER_r2; S3=/mnt/archive/AI_DNA_ANALYZER_stage3; ST2=/mnt/archive/AI_DNA_ANALYZER_stage2; DR=/mnt/archive/AI_DNA_ANALYZER_depth_robust
REF=/mnt/archive/AI_DNA_ANALYZER_benchmark/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
M=/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage3/frozen/R1.model
BIN=/mnt/archive/AI_DNA_ANALYZER_r2/build/dnav2
truth() { case $1 in
  HG002) echo /mnt/archive/AI_DNA_ANALYZER_benchmark/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz;;
  HG003) echo /mnt/archive/AI_DNA_ANALYZER_final_eval/truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz;;
  HG004) echo $ST2/truth/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz;; esac; }
bedof() { case $1 in
  HG002) echo /mnt/archive/AI_DNA_ANALYZER_benchmark/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed;;
  HG003) echo /mnt/archive/AI_DNA_ANALYZER_final_eval/truth/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed;;
  HG004) echo $ST2/truth/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed;; esac; }
declare -A LEN=([chr20]=64444167 [chr21]=46709983 [chr22]=50818468 [chr8]=145138636)
dump() { # name sample contig bam
  local name=$1 s=$2 c=$3 bam=$4 out=$R2/dump2/$1
  [ -s $out/features.tsv ] && return 0
  mkdir -p $out; grep -P "^$c\t" $(bedof $s) > $out/conf.bed
  nice -n 10 $BIN run --bam $bam --ref $REF --truth $(truth $s) --bed $out/conf.bed --out-dir $out --contig $c --contig-len ${LEN[$c]} \
    --sample $s --tag $name --threads 4 --mode cascade --model $M --no-feat-pb --r2 --ood off --dump-features $out/features.tsv > $out/stdout.json 2> $out/stderr.txt
  echo "$name rc=$? $(wc -l < $out/features.tsv) rows" >> $R2/logs/dump.log
  rm -f $out/ai_cascade_*.vcf $out/routed_loci_*.tsv
}
down() { # src dst target_depth seed contig sample bed -> downsample by fraction target/native
  local src=$1 dst=$2 tgt=$3 seed=$4 c=$5 s=$6
  [ -s $dst ] && return 0
  local nat; nat=$(grep -P "^$c\t" $(bedof $s) > $R2/logs/$c.$s.bed; samtools bedcov $R2/logs/$c.$s.bed $src | awk '{b+=$3-$2; d+=$4} END{printf "%.3f", d/b}')
  local frac; frac=$(python3 -c "print(f'{$tgt/$nat:.4f}'[1:])")
  samtools view -b -s $seed$frac -@ 2 -o $dst $src && samtools index $dst
  echo "down $dst native=$nat seed=$seed$frac" >> $R2/logs/dump.log
}
# --- training/dev BAMs
for d in 10 15 20; do down $ST2/bam/HG004_chr21_full.bam $R2/bam/HG004_chr21_${d}x_s301.bam $d 301 chr21 HG004; done
for d in 10 15 20; do down $S3/bam/HG004_chr22_native.bam $R2/bam/HG004_chr22_${d}x_s301.bam $d 301 chr22 HG004; done
# --- dumps (role in name: tr = training, va = validation)
dump tr_hg002_chr20_30x HG002 chr20 $DR/bam/HG002_chr20_f010.bam
dump tr_hg002_chr20_60x HG002 chr20 $DR/bam/HG002_chr20_f020.bam
dump tr_hg002_chr20_100x HG002 chr20 $DR/bam/HG002_chr20_f033.bam
dump tr_hg003_chr8_30x HG003 chr8 $DR/bam/HG003_chr8_test1_f010_s7.bam
dump tr_hg004_chr21_6x HG004 chr21 $ST2/bam/HG004_chr21_f010_s21.bam
dump tr_hg004_chr21_30x HG004 chr21 $ST2/bam/HG004_chr21_f053_s21.bam
dump tr_hg004_chr21_56x HG004 chr21 $ST2/bam/HG004_chr21_full.bam
for d in 10 15 20; do dump tr_hg004_chr21_${d}x HG004 chr21 $R2/bam/HG004_chr21_${d}x_s301.bam; done
dump tr_hg004_chr20_6x HG004 chr20 $ST2/bam/HG004_chr20_f010_s21.bam
dump tr_hg004_chr20_30x HG004 chr20 $ST2/bam/HG004_chr20_f053_s21.bam
dump tr_hg004_chr20_56x HG004 chr20 $ST2/bam/HG004_chr20_full.bam
dump va_hg004_chr22_30x HG004 chr22 $S3/bam/HG004_chr22_30x_s101.bam
dump va_hg004_chr22_68x HG004 chr22 $S3/bam/HG004_chr22_native.bam
for d in 10 15 20; do dump va_hg004_chr22_${d}x HG004 chr22 $R2/bam/HG004_chr22_${d}x_s301.bam; done
echo ALLDONE >> $R2/logs/dump.log
