#!/bin/bash
# run_happy.sh LABEL   : real hap.py 0.3.15 (same image, truth, reference, xcmp, --threads 8) on RUNDIR/LABEL's VCF, for
# NOTE: docker data-root (loopback ext4) is not mounted after the reboot and needs sudo; the SAME image
# (quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0, digest of the pulled manifest) is run with rootless podman instead.
# region = full | dev | holdout.  Only the -f BED differs between regions (dev/holdout = the confident BED cut at the
# pre-registered 2-Mb blocks).  Output: $V2W/happy/LABEL/{full,dev,holdout}/
set -uo pipefail
source "$(dirname "$0")/env.sh"
LABEL=$1; RUN=${2:-$V2W/runs/${LABEL}_t8_warm}
IMG=quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0
H=$V2W/happy/$LABEL; mkdir -p $H/query
python3 - <<PY
import pysam,shutil
src="$RUN/ai_cascade_chr20_300x.vcf"; dst="$H/query/ai_${LABEL}_chr20_300x.vcf"
shutil.copy(src,dst); pysam.tabix_index(dst,preset="vcf",force=True)
PY
for R in full dev holdout; do
  mkdir -p $H/$R
  /usr/bin/time -v -o $H/$R/time.txt podman --root /dev/shm/podman_root --runroot /dev/shm/podman_run --storage-driver vfs run --rm --security-opt label=disable \
    -v $A/truth_HG002_GRCh38_v4.2.1:/truth -v $A/reference:/ref -v $V2W/beds:/bed -v $H/query:/query -v $H/$R:/out \
    $IMG /usr/local/bin/hap.py /truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz /query/ai_${LABEL}_chr20_300x.vcf.gz \
    -f /bed/chr20_$R.bed -r /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna -o /out/ai_$LABEL -l chr20 --threads 8 \
    > $H/$R/log.txt 2>&1
  echo "$LABEL $R rc=$? $(grep Elapsed $H/$R/time.txt | awk '{print $NF}')"
done
