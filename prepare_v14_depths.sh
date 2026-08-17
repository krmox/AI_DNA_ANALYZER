#!/bin/bash
# bench_v14 depth realisations. Each region's native BAM is downsampled to
# nominal 30x and 15x with `samtools view -s SEED.FRACTION`, where FRACTION is
# derived from the region's *measured* native mean depth (samtools depth -a over
# the region), not guessed. SEED is fixed at 42 for every region and depth, so
# the subsampling is reproducible and the 15x set is a subset of nothing in
# particular -- each realisation is an independent draw from the same library,
# exactly as in devlogs 6-13.
#
# These are downsamples of ONE library. They are three regimes of one region,
# never three independent regions (devlog 14 §5.4).
set -u
D=data/giab_hg002_v14
SEED=42
mkdir -p logs_v14

# name contig start stop
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
  prepare chr20_neutral chr20  33000000  34000000
  prepare chr19_gcrich  chr19   1000000   2000000
  prepare chr4_atrich   chr4  101000000 102000000
  prepare chr1_segdup   chr1  120000000 121000000
} 2>&1 | tee logs_v14/depths.log
echo "DEPTHS READY $(date -Is)"
