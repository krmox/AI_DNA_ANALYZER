#!/usr/bin/env bash
# Pre-specified reference run (PREREGISTRATION.md): same bcftools baseline on the 300x BAM, then hap.py.
set -uo pipefail
W=/mnt/archive/AI_DNA_ANALYZER_depth30; B=/mnt/archive/AI_DNA_ANALYZER_benchmark
REF=$B/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
BAM=$B/HG002_GRCh38_chr20/HG002.GRCh38.300x_chr20.bam
T=$B/truth_HG002_GRCh38_v4.2.1; IMG=quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0
out=$W/baseline/bcftools_300x; mkdir -p $out
/usr/bin/time -v -o $out/time.txt bash -c "bcftools mpileup -f $REF -r chr20 -q 20 -Q 13 -d 10000 -a FORMAT/AD -Ou $BAM | bcftools call -mv -Ou | bcftools view -v snps -i 'QUAL>=30' -Oz -o $out/bcftools_300x.vcf.gz" 2> $out/stderr.txt
echo "rc=$?" > $out/rc.txt; bcftools index -t $out/bcftools_300x.vcf.gz
cp $out/bcftools_300x.vcf.gz* $W/happy/query/
for scope in whole holdout; do
  extra=(); [ $scope = holdout ] && extra=(-T /w/holdout_odd_blocks.bed)
  podman --root /dev/shm/podman_root --storage-driver vfs run --rm --security-opt label=disable \
   -v $T:/truth -v $B/reference:/ref -v $W/happy/query:/query -v $W/happy/out:/out -v $W:/w $IMG \
   /usr/local/bin/hap.py /truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz /query/bcftools_300x.vcf.gz \
   -f /truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -r /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
   -o /out/bcftools_300x_$scope -l chr20 --threads 8 "${extra[@]}" > $W/happy/out/bcftools_300x_$scope.log 2>&1
done
echo done > $W/happy/out/BC300_DONE
