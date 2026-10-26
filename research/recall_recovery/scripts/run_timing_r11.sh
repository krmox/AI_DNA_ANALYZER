#!/usr/bin/env bash
# R1.1 final timing: 5 interleaved reps, --threads 8, HG004 ~30x. Variants: R11 (defaults), R1flags (same binary, --strict-windows --truth-mask), R1archived (archived V2.x binary, truth required).
set -uo pipefail
S=/mnt/archive/AI_DNA_ANALYZER_stage2; R=/mnt/archive/AI_DNA_ANALYZER_depth_robust; B=/mnt/archive/AI_DNA_ANALYZER_benchmark
REF=$B/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
T=$S/truth/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz; BD=$S/truth/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
NEW=$S/build_r11/dnav2; OLD=/mnt/archive/AI_DNA_ANALYZER_v2x/build_new/dnav2
declare -A LEN=([chr21]=46709983 [chr20]=64444167)
mkdir -p $S/timing_r11; printf "contig\trep\tvariant\twall_s\tuser_s\tsys_s\tmaxrss_kb\n" > $S/timing_r11/timing.tsv
for rep in 1 2 3 4 5; do for c in chr21 chr20; do
  for v in R1archived R1flags R11; do
    o=$S/timing_r11/${c}_${v}_r$rep; mkdir -p $o
    case $v in R1archived) bin=$OLD; x="--truth $T";; R1flags) bin=$NEW; x="--truth $T --strict-windows --truth-mask";; R11) bin=$NEW; x="";; esac
    /usr/bin/time -v -o $o/time.txt $bin run --bam $S/bam/HG004_${c}_f053_s21.bam --ref $REF --bed $BD $x --out-dir $o --contig $c --contig-len ${LEN[$c]} --sample HG004 --tag t --threads 8 --mode cascade --model $R/R1.model > $o/out.json 2> $o/err.txt
    w=$(grep Elapsed $o/time.txt | sed 's/.*: //'); u=$(grep "User time" $o/time.txt | awk '{print $NF}'); s=$(grep "System time" $o/time.txt | awk '{print $NF}'); m=$(grep "Maximum resident" $o/time.txt | awk '{print $NF}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $c $rep $v $w $u $s $m >> $S/timing_r11/timing.tsv
  done
  # determinism + equivalence
  b=$S/timing_r11/${c}
  cmp -s <(grep -v '^#' $b\_R11_r$rep/ai_cascade_t.vcf) <(grep -v '^#' $b\_R11_r1/ai_cascade_t.vcf) || echo "NONDET R11 $c rep$rep" >> $S/timing_r11/flags.txt
  cmp -s <(grep -v '^#' $b\_R1flags_r$rep/ai_cascade_t.vcf) <(grep -v '^#' $b\_R1archived_r$rep/ai_cascade_t.vcf) || echo "R1flags!=R1archived $c rep$rep" >> $S/timing_r11/flags.txt
done; done
echo R11TIMING_DONE >> $S/logs/timing.log
