#!/usr/bin/env bash
# One-command check of a checkout of AI_DNA_ANALYZER.
#
#   bash verify.sh          build everything and run the offline tests (no downloads, ~2 min)
#   bash verify.sh --data   also re-run V2.x A+C on 250 kb of real GIAB data and compare the calls
#                           with the archived final-evaluation VCF (~60 MB download)
#
# Every step prints PASS or FAIL; the script exits non-zero if anything failed.
# Downloads are resumable: if step 4 is interrupted or a download fails, just run it again.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$PWD
WITH_DATA=0; [ "${1:-}" = "--data" ] && WITH_DATA=1
LOG=$ROOT/data/verify_logs; mkdir -p "$LOG"
FAILED=0
step() { printf '\n=== %s ===\n' "$1"; }
ok()   { printf 'PASS  %s\n' "$1"; }
bad()  { printf 'FAIL  %s\n' "$1"; FAILED=1; }
retry() { local i; for i in 1 2 3 4 5; do "$@" && return 0; echo "  retry $i: $*" | cut -c1-120; sleep 5; done; return 1; }

# ---------------------------------------------------------------------------------------------
step "0. Prerequisites"
need() { command -v "$1" >/dev/null 2>&1 && ok "$1 found" || bad "$1 not found ($2)"; }
need python3 "install Python 3"
need gcc     "install gcc"
need g++     "install g++ (C++20)"
need make    "install make"
if [ $WITH_DATA = 1 ]; then need samtools "needed for --data"; need bcftools "needed for --data"; need curl "needed for --data"; fi
python3 -c 'import numpy, scipy, pysam, pytest, torch' 2>/dev/null \
  && ok "python packages (numpy scipy pysam pytest torch)" \
  || bad "python packages missing: pip install -r requirements.txt"
# The native code is compiled against pysam's bundled htslib headers and linked to the system
# libhts.so.3. The two htslib versions must match, or the build fails or misbehaves.
HTS=$(python3 - <<'EOF' 2>/dev/null
import ctypes, pysam.version as v
lib = ctypes.CDLL("libhts.so.3"); lib.hts_version.restype = ctypes.c_char_p
print(v.__htslib_version__, lib.hts_version().decode())
EOF
)
if [ -z "$HTS" ]; then
  bad "system libhts.so.3 not found (Fedora: dnf install htslib; Debian/Ubuntu: apt install libhts3)"
else
  set -- $HTS
  if [ "${1%.*}" = "${2%.*}" ]; then ok "htslib versions match (pysam $1, system $2)"
  else bad "htslib mismatch: pysam bundles $1, system has $2 (use pysam==0.24.0 with system htslib 1.23, see README)"; fi
fi
[ $FAILED = 1 ] && { echo; echo "Fix the prerequisites above and re-run."; exit 1; }

# ---------------------------------------------------------------------------------------------
step "1. Build"
bash native/build.sh > "$LOG/build_v1.log" 2>&1 \
  && ok "V1 native backends (native/build.sh)" || bad "V1 native build, see $LOG/build_v1.log"
make -C native_v2x engine kernels > "$LOG/build_v2.log" 2>&1 \
  && ok "V2/V2.x engine native_v2x/build/dnav2" || bad "V2 engine build, see $LOG/build_v2.log"

# ---------------------------------------------------------------------------------------------
step "2. V1 test suite (frozen cascade, native backends)"
python3 -m pytest -q -p no:cacheprovider \
  test_cascade.py test_cheap_router.py test_native_pileup_equivalence.py \
  test_native_reads_equivalence.py test_hg005_extraction_integrity.py \
  test_counts_backend_fuzz.py test_bench_v12.py test_robustness_benchmark.py > "$LOG/test_v1.log" 2>&1
rc=$?; tail -1 "$LOG/test_v1.log"
[ $rc = 0 ] && ok "V1 tests (xfailed = documented adversarial cases)" || bad "V1 tests, see $LOG/test_v1.log"

# ---------------------------------------------------------------------------------------------
step "3. V2 kernels bit-identical to the frozen V1 NumPy/SciPy code"
make -s -C native_v2x test > "$LOG/test_v2.log" 2>&1; rc=$?
grep -E '^\[(PASS|FAIL)\]' "$LOG/test_v2.log" | cut -c1-90
[ $rc = 0 ] && ! grep -q '^\[FAIL\]' "$LOG/test_v2.log" && ok "V2 kernel tests" || bad "V2 kernel tests, see $LOG/test_v2.log"

# ---------------------------------------------------------------------------------------------
if [ $WITH_DATA = 1 ]; then
  step "4. Real data: V2.x A+C on HG003 chr8:20,000,001-20,249,984 (300x) vs archived VCF"
  D=$ROOT/data/verify_hg003_chr8; mkdir -p "$D"; cd "$D"
  # Region bounds are multiples of 64 so that no 64-bp window straddles them.
  CONTIG=chr8; LEN=145138636; S=20000000; E=20249984; REG="$CONTIG:$((S+1))-$E"
  GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
  BAM_URL=$GIAB/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.300x.bam
  TRUTH_URL=$GIAB/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark
  # Same GRCh38 analysis-set sequence as the paper's reference (identical on this region).
  REF_URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
  CURL=(curl -fsSL --retry 5 --retry-delay 5 --connect-timeout 30)

  # Reads. Piece i holds the reads that START in its sub-range; piece 0 also keeps the reads that
  # start before the region and overlap it. Concatenated in order, this is the full sorted slice.
  if [ ! -s slice.bam.bai ]; then
    echo "downloading reads for $REG (~45 MB in 8 parallel pieces; NCBI can be slow) ..."
    [ -s remote.bai ] || { retry "${CURL[@]}" -o remote.bai.tmp "$BAM_URL.bai" && mv remote.bai.tmp remote.bai; }
    piece() {
      local i=$1 a=$2 b=$3 filt=()
      [ -f part$i.done ] && return 0
      [ "$i" -gt 0 ] && filt=(-e "pos>=$((a+1))")
      samtools view --no-PG -b "${filt[@]}" -o part$i.bam -X "$BAM_URL" remote.bai "$CONTIG:$((a+1))-$b" 2>> download.log \
        && samtools quickcheck part$i.bam && touch part$i.done
    }
    N=8; STEP=$(( (E - S) / N ))
    for i in $(seq 0 $((N-1))); do
      a=$(( S + i*STEP )); b=$(( i == N-1 ? E : a + STEP ))
      retry piece $i $a $b > /dev/null &
    done
    wait
    if [ "$(ls part?.done 2>/dev/null | wc -l)" = $N ]; then
      samtools cat -o slice.bam part0.bam part1.bam part2.bam part3.bam part4.bam part5.bam part6.bam part7.bam \
        && samtools index slice.bam && rm -f part?.bam part?.done
    fi
  fi
  [ -s slice.bam.bai ] && ok "reads downloaded ($(du -h slice.bam | cut -f1))" || bad "read download incomplete; run the same command again to resume"

  # Truth VCF (region only) and confident-region BED clipped to the region.
  if [ ! -s truth.vcf.gz.tbi ]; then
    retry bcftools view -r "$REG" -Oz -o truth.vcf.gz "$TRUTH_URL.vcf.gz" && bcftools index -t truth.vcf.gz
  fi
  if [ ! -s region.bed ]; then
    retry "${CURL[@]}" -o full.bed.tmp "${TRUTH_URL}_noinconsistent.bed" \
      && awk -v s=$S -v e=$E '$1=="chr8"{a=$2<s?s:$2; b=$3>e?e:$3; if(a<b) print $1"\t"a"\t"b}' full.bed.tmp > region.bed \
      && rm -f full.bed.tmp
  fi
  [ -s truth.vcf.gz.tbi ] && [ -s region.bed ] && ok "truth VCF and BED downloaded" || bad "truth download; run again to resume"

  # Reference: full-length chr8 with the real sequence on the region and N elsewhere (only the region is scored).
  if [ ! -s ref.fa.fai ]; then
    retry samtools faidx -o region.fa "$REF_URL" "$REG" && python3 - "$LEN" "$S" <<'EOF' && samtools faidx ref.fa
import sys
n, s = int(sys.argv[1]), int(sys.argv[2])
seq = "".join(l.strip() for l in open("region.fa") if not l.startswith(">"))
assert len(seq) > 0, "empty reference download"
full = "N" * s + seq + "N" * (n - s - len(seq))
with open("ref.fa", "w") as f:
    f.write(">chr8\n")
    for i in range(0, n, 60):
        f.write(full[i:i + 60] + "\n")
EOF
  fi
  [ -s ref.fa.fai ] && ok "reference sequence downloaded" || bad "reference download; run again to resume"

  if [ $FAILED = 0 ]; then
    rm -rf out; mkdir -p out
    "$ROOT/native_v2x/build/dnav2" run --bam slice.bam --ref ref.fa --truth truth.vcf.gz --bed region.bed \
        --out-dir out --contig $CONTIG --contig-len $LEN --sample HG003 --tag chr8_300x --threads 8 \
        --mode cascade --model "$ROOT/native_v2x/models/A_C.model" > out/summary_stdout.json 2> out/stderr.txt \
      && ok "dnav2 run" || bad "dnav2 run, see $D/out/stderr.txt"
    ARCH=$ROOT/research/final_independent_eval/run/ai_A_C_chr8_300x.vcf.gz
    grep -v '^#' out/ai_cascade_chr8_300x.vcf > new.txt 2>/dev/null
    bcftools view -H -r "$REG" "$ARCH" > archived.txt
    echo "SNP calls on $REG: this run $(wc -l < new.txt), archived v2.0.0 result $(wc -l < archived.txt)"
    if [ -s archived.txt ] && cmp -s new.txt archived.txt; then
      ok "V2.x calls identical to the archived HG003 chr8 final-evaluation VCF"
    else
      bad "V2.x calls differ from the archive (diff $D/new.txt $D/archived.txt)"
    fi
  fi
  cd "$ROOT"
fi

# ---------------------------------------------------------------------------------------------
echo
if [ $FAILED = 0 ]; then echo "ALL CHECKS PASSED"; else echo "SOME CHECKS FAILED"; fi
exit $FAILED
