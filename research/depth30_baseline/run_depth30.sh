#!/usr/bin/env bash
# See PREREGISTRATION.md. Single pass, no retries.
set -uo pipefail
W=/mnt/archive/AI_DNA_ANALYZER_depth30
B=/mnt/archive/AI_DNA_ANALYZER_benchmark
BIN=/mnt/archive/AI_DNA_ANALYZER_v2x/build_x/dnav2
MODEL=/mnt/archive/AI_DNA_ANALYZER_v2x/repo/native_v2x/models/A_C.model
REF=$B/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
TRUTH=$B/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
BED=$B/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
BAM300=$B/HG002_GRCh38_chr20/HG002.GRCh38.300x_chr20.bam
BAM30=$W/bam/HG002.GRCh38.30x_chr20.bam
cd $W
date -Is > logs/start.txt
sha256sum $BIN $MODEL > logs/binary_model_sha256.txt
samtools view -@4 -b -s 42.1 -o $BAM30 $BAM300 chr20 && samtools index -@4 $BAM30
samtools coverage -r chr20 $BAM30 > logs/coverage_30x.txt
samtools flagstat $BAM30 > logs/flagstat_30x.txt
python3 - <<PY
with open("$W/holdout_odd_blocks.bed","w") as f:
    L=64444167; b=1
    while b*2000000<L:
        f.write(f"chr20\t{b*2000000}\t{min((b+1)*2000000,L)}\n"); b+=2
PY
run_engine() { # tag bam extra...
  tag=$1; bam=$2; shift 2; out=$W/runs/$tag; mkdir -p $out
  /usr/bin/time -v -o $out/time.txt $BIN run --bam $bam --ref $REF --truth $TRUTH --bed $BED --out-dir $out \
    --contig chr20 --contig-len 64444167 --sample HG002 --tag $tag --threads 8 --mode cascade "$@" > $out/stdout.json 2> $out/stderr.txt
  echo "rc=$?" > $out/rc.txt
}
run_engine v2control_30x $BAM30
run_engine v2x_AC_30x $BAM30 --model $MODEL
run_baseline() { # tag bam
  tag=$1; bam=$2; out=$W/baseline/$tag; mkdir -p $out
  /usr/bin/time -v -o $out/time.txt bash -c "bcftools mpileup -f $REF -r chr20 -q 20 -Q 13 -d 10000 -a FORMAT/AD -Ou $bam | bcftools call -mv -Ou | bcftools view -v snps -i 'QUAL>=30' -Oz -o $out/$tag.vcf.gz" 2> $out/stderr.txt
  echo "rc=$?" > $out/rc.txt; bcftools index -t $out/$tag.vcf.gz
}
run_baseline bcftools_30x $BAM30
date -Is > logs/engine_done.txt
