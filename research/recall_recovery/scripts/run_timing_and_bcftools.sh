#!/usr/bin/env bash
# Stage 2 timing (STAGE2_PREREG.md): 5 interleaved reps of B0/H1/H2/H1H2, --threads 8, 30x HG004 chr21+chr20; then bcftools comparator.
set -uo pipefail
S=/mnt/archive/AI_DNA_ANALYZER_stage2; R=/mnt/archive/AI_DNA_ANALYZER_depth_robust; B=/mnt/archive/AI_DNA_ANALYZER_benchmark
REF=$B/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
T=$S/truth/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz; BD=$S/truth/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
declare -A LEN=([chr21]=46709983 [chr20]=64444167)
mkdir -p $S/timing; : > $S/timing/timing.tsv
echo -e "contig\trep\tvariant\twall_s\tuser_s\tsys_s\tmaxrss_kb" >> $S/timing/timing.tsv
for rep in 1 2 3 4 5; do for c in chr21 chr20; do for v in B0: H1:--edge-windows H2:--truth-free H1H2:"--edge-windows --truth-free"; do
  n=${v%%:*}; f=${v#*:}; o=$S/timing/${c}_${n}_r$rep; mkdir -p $o
  /usr/bin/time -v -o $o/time.txt $S/build/dnav2 run --bam $S/bam/HG004_${c}_f053_s21.bam --ref $REF --truth $T --bed $BD --out-dir $o --contig $c --contig-len ${LEN[$c]} --sample HG004 --tag t --threads 8 --mode cascade --model $R/R1.model $f > $o/out.json 2> $o/err.txt
  w=$(grep "Elapsed" $o/time.txt | sed 's/.*: //'); u=$(grep "User time" $o/time.txt | awk '{print $NF}'); s=$(grep "System time" $o/time.txt | awk '{print $NF}'); m=$(grep "Maximum resident" $o/time.txt | awk '{print $NF}')
  echo -e "$c\t$rep\t$n\t$w\t$u\t$s\t$m" >> $S/timing/timing.tsv
  cmp -s $o/ai_cascade_t.vcf $S/timing/${c}_${n}_r1/ai_cascade_t.vcf || echo "NONDETERMINISTIC $c $n rep$rep" >> $S/timing/timing.tsv
done; done; done
echo TIMING_DONE >> $S/logs/timing.log
for c in chr21 chr20; do
  out=$S/runs/bcftools_${c}_30x; mkdir -p $out
  /usr/bin/time -v -o $out/time.txt bash -c "bcftools mpileup -f $REF -r $c -q 20 -Q 13 -d 10000 -a FORMAT/AD -Ou $S/bam/HG004_${c}_f053_s21.bam | bcftools call -mv -Ou | bcftools view -v snps -i 'QUAL>=30' -Oz -o $out/bc.vcf.gz" 2> $out/stderr.txt
  bcftools index -t $out/bc.vcf.gz
done
echo BCFTOOLS_DONE >> $S/logs/timing.log
