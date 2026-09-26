#!/bin/bash
# engine (8 thr, warm) -> engine-vs-offline gate -> hap.py (full/dev/holdout) once per distinct VCF
source "$(dirname "$0")/env.sh"
export BIN=/mnt/archive/AI_DNA_ANALYZER_v2x/build_x/dnav2
cd "$(dirname "$0")/.."
declare -A SEEN
NAMES=(A_C A B C D A_B A_D B_C B_D C_D A_B_C A_B_D A_C_D B_C_D A_B_C_D)
for m in "${NAMES[@]}"; do
  L=M_$m
  [ -f $V2W/runs/${L}_t8_warm/stdout.json ] || scripts/run_variant.sh $L 8 0 --mode cascade --model models/$m.model
  python3 analysis/s07_check_engine.py $L "${m//_/+}" | tee -a $V2W/gate_engine_vs_offline.log
  h=$(grep -v '^##' $V2W/runs/${L}_t8_warm/ai_cascade_chr20_300x.vcf | md5sum | cut -d' ' -f1)
  if [ -n "${SEEN[$h]:-}" ]; then echo "$L identical VCF to ${SEEN[$h]} -> hap.py reused" | tee -a $V2W/happy_reuse.log; echo "${SEEN[$h]}" > $V2W/runs/${L}_t8_warm/SAME_VCF_AS
  else SEEN[$h]=$L; [ -f $V2W/happy/$L/holdout/ai_$L.summary.csv ] || scripts/run_happy.sh $L; fi
done
echo ALL_FINALISTS_DONE
