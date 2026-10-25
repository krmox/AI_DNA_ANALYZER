#!/usr/bin/env bash
# hap.py 0.3.15 (same image/options as earlier benchmarks), SNP row. Whole chr20 and holdout (odd 2-Mb blocks via -T).
set -uo pipefail
W=/mnt/archive/AI_DNA_ANALYZER_depth30; B=/mnt/archive/AI_DNA_ANALYZER_benchmark
T=$B/truth_HG002_GRCh38_v4.2.1; IMG=quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0
mkdir -p $W/happy/query $W/happy/out
declare -A SRC=( [v2control_30x]=$W/runs/v2control_30x/ai_cascade_v2control_30x.vcf [v2x_AC_30x]=$W/runs/v2x_AC_30x/ai_cascade_v2x_AC_30x.vcf )
for k in "${!SRC[@]}"; do bcftools view -Oz -o $W/happy/query/$k.vcf.gz ${SRC[$k]}; bcftools index -f -t $W/happy/query/$k.vcf.gz; done
for f in $W/baseline/*/*.vcf.gz; do n=$(basename $f .vcf.gz); cp $f $W/happy/query/$n.vcf.gz; cp $f.tbi $W/happy/query/$n.vcf.gz.tbi; done
for q in $W/happy/query/*.vcf.gz; do n=$(basename $q .vcf.gz)
 for scope in whole holdout; do
  extra=(); [ $scope = holdout ] && extra=(-T /w/holdout_odd_blocks.bed)
  podman --root /dev/shm/podman_root --storage-driver vfs run --rm --security-opt label=disable \
   -v $T:/truth -v $B/reference:/ref -v $W/happy/query:/query -v $W/happy/out:/out -v $W:/w $IMG \
   /usr/local/bin/hap.py /truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz /query/$n.vcf.gz \
   -f /truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -r /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
   -o /out/${n}_$scope -l chr20 --threads 8 "${extra[@]}" > $W/happy/out/${n}_$scope.log 2>&1
  echo "$n $scope rc=$?" >> $W/happy/out/rc.txt
 done
done
echo done > $W/happy/out/ALLDONE
