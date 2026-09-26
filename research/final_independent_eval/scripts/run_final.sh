#!/usr/bin/env bash
# Single final run of V2.x A+C on HG003 chr8 (see PREREGISTRATION.md). No retries, no parameter changes.
set -uo pipefail
E=/mnt/archive/AI_DNA_ANALYZER_final_eval
BIN=/mnt/archive/AI_DNA_ANALYZER_v2x/build_x/dnav2
MODEL=/mnt/archive/AI_DNA_ANALYZER_v2x/repo/native_v2x/models/A_C.model
REF=/mnt/archive/AI_DNA_ANALYZER_benchmark/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
BAM=$E/bam/HG003.GRCh38.300x_chr8.bam
TRUTH=$E/truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
BED=$E/truth/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
OUT=$E/runs/A_C_chr8; mkdir -p $OUT
sha256sum $BIN $MODEL > $OUT/binary_model_sha256.txt
free -m > $OUT/free_before.txt; date -Is > $OUT/start.txt
CMD=(/usr/bin/time -v -o $OUT/time.txt $BIN run --bam $BAM --ref $REF --truth $TRUTH --bed $BED --out-dir $OUT --contig chr8 --contig-len 145138636 --sample HG003 --tag chr8_300x --threads 8 --mode cascade --model $MODEL)
printf '%q ' "${CMD[@]}" > $OUT/exact_command.txt; echo >> $OUT/exact_command.txt
"${CMD[@]}" > $OUT/stdout.json 2> $OUT/stderr.txt
echo "engine_rc=$?" > $OUT/rc.txt; date -Is > $OUT/end.txt
md5sum $OUT/ai_cascade_chr8_300x.vcf; sha256sum $OUT/ai_cascade_chr8_300x.vcf | tee $OUT/vcf_sha256.txt
