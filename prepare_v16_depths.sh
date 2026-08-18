#!/bin/bash
# bench_v16 depth realisations. Identical recipe to prepare_v14_depths.sh /
# prepare_v15_depths.sh: the native BAM is downsampled to nominal 30x and 15x
# with `samtools view -s SEED.FRACTION`, where FRACTION comes from the region's
# *measured* native mean depth (samtools depth -a), not guessed. SEED is fixed
# at 42.
#
# These are downsamples of ONE library: three depth regimes of ONE region, never
# three independent regions.
set -u
D=data/giab_hg002_v16
SEED=42
mkdir -p logs_v16

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
  prepare chr13_representative chr13 70000000 73000000
} 2>&1 | tee logs_v16/depths.log
echo "DEPTHS READY $(date -Is)"
