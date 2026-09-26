#!/usr/bin/env bash
# hap.py 0.3.15 (same image and options as the chr20 benchmarks) on the single chr8 VCF. rootless podman, vfs storage in RAM.
set -uo pipefail
E=/mnt/archive/AI_DNA_ANALYZER_final_eval; RUN=$E/runs/A_C_chr8; H=$E/happy; mkdir -p $H/query $H/out
IMG=quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0
python3 - <<PY
import pysam,shutil
dst="$H/query/ai_A_C_chr8_300x.vcf"; shutil.copy("$RUN/ai_cascade_chr8_300x.vcf",dst); pysam.tabix_index(dst,preset="vcf",force=True)
PY
podman --root /dev/shm/podman_root --storage-driver vfs run --rm --security-opt label=disable \
  -v $E/truth:/truth -v /mnt/archive/AI_DNA_ANALYZER_benchmark/reference:/ref -v $H/query:/query -v $H/out:/out $IMG \
  /usr/local/bin/hap.py /truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz /query/ai_A_C_chr8_300x.vcf.gz \
  -f /truth/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -r /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
  -o /out/ai_A_C_chr8 -l chr8 --threads 8 > $H/out/hap.log 2>&1
echo "happy_rc=$?" > $H/out/rc.txt
