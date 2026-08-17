#!/bin/bash
# bench_v15 depth realisations. Identical recipe to prepare_v14_depths.sh: each
# region's native BAM is downsampled to nominal 30x and 15x with
# `samtools view -s SEED.FRACTION`, where FRACTION comes from the region's
# *measured* native mean depth (samtools depth -a), not guessed. SEED is fixed
# at 42 for every region and depth.
#
# These are downsamples of ONE library. They are three regimes of one region,
# never three independent regions (devlog 15 §5.3).
set -u
D=data/giab_hg002_v15
SEED=42
mkdir -p logs_v15

prepare () {
  local name=$1 contig=$2 start=$3 stop=$4
  local region="${contig}:$((start + 1))-${stop}"
  local native
  native=$(samtools depth -a -r "$region" "$D/${name}_full.bam" \
           | awk '{s+=$3; n++} END{printf "%.4f", (n ? s/n : 0)}')
  echo "$name native_mean_depth=$native"
  for target in 30 15; do
    local out="$D/${name}_${target}x.bam"
    local frac
    frac=$(python3 -c "print(f'{min(0.999, $target/$native):.4f}'[1:])")
    if [ ! -s "$out" ]; then
      samtools view -b -s "${SEED}${frac}" "$D/${name}_full.bam" -o "$out" \
        && samtools index "$out"
    fi
    local got
    got=$(samtools depth -a -r "$region" "$out" | awk '{s+=$3; n++} END{printf "%.2f", (n ? s/n : 0)}')
    echo "  ${target}x: fraction=${SEED}${frac} achieved_mean_depth=$got"
  done
}

{
  prepare chr16_segdup  chr16 29000000 30000000
  prepare chr15_segdup  chr15 30000000 31000000
  prepare chr7_segdup   chr7  57000000 58000000
  prepare chr14_neutral chr14 53000000 54000000
} 2>&1 | tee logs_v15/depths.log
echo "DEPTHS READY $(date -Is)"
