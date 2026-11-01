#!/usr/bin/env bash
# Verification of a checkout of AI_DNA_ANALYZER, in separately selectable groups.
#
#   bash verify.sh                 default = groups "v1 v2" (same checks as earlier releases; PARTIAL, see below)
#   bash verify.sh full            env + v1 + v2 + r2 + nindels   (complete offline verification)
#   bash verify.sh full --data     ... plus the real-data integration test (~60 MB download)
#   bash verify.sh r2 nindels      any combination of groups
#   bash verify.sh --data          same as "v1 v2 data"
#
# Groups:
#   env      machine/toolchain/htslib diagnostics only (always run first; changes nothing)
#   python   Python-only tests: V1 suite WITHOUT native backends; native-dependent tests are skipped and
#            counted as skipped, so this group can never complete a verification
#   native   compile everything (V1 C backends, V2.x and R2 engines, nindels); no tests
#   v1       build V1 C backends + full V1 pytest suite (expects exactly the documented skips)
#   v2       build V2/V2.x engine + V2 kernels bit-identical to the frozen V1 NumPy/SciPy code
#   r2       build R1.1/R2 engine + R2 tail tests + call-domain (frame) tests
#   nindels  build nindels (CMake) + its unit tests
#   data     real GIAB data: V2.x A+C on HG003 chr8:20,000,001-20,249,984 vs the archived VCF
#   full     env v1 v2 r2 nindels   (+ data with --data)
#
# Every step prints PASS, FAIL or SKIP and the script ends with exact passed/failed/skipped counts.
# It says "COMPLETE VERIFICATION PASSED" only when `full` ran without failures; any other
# selection ends with "PARTIAL: ..." and lists what was NOT run.  Nothing is installed or downloaded
# except by the `data` group (resumable; re-run to continue).  Run `bash scripts/hts_env.sh report`
# to see which htslib is used.
set -uo pipefail
cd "$(dirname "$0")"
ROOT=$PWD
LOG=$ROOT/data/verify_logs; mkdir -p "$LOG"
FAILED=0
REF_HTS_MINOR="1.23"          # htslib series used for the published results
EXPECTED_V1_SKIPS=2           # documented: need a local HG002 slice (README "Reproducibility")

step() { printf '\n=== %s ===\n' "$1"; }
ok()   { printf 'PASS  %s\n' "$1"; }
bad()  { printf 'FAIL  %s\n' "$1"; FAILED=1; }
skip() { printf 'SKIP  %s\n' "$1"; }
warn() { printf 'WARN  %s\n' "$1"; }
retry() { local i; for i in 1 2 3 4 5; do "$@" && return 0; echo "  retry $i: $*" | cut -c1-120; sleep 5; done; return 1; }

# ---- argument parsing -------------------------------------------------------------------------
WITH_DATA=0; GROUPS_REQ=()
for a in "$@"; do
  case "$a" in
    --data) WITH_DATA=1 ;;
    -h|--help) sed -n '2,/^set -uo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit 0 ;;
    env|python|native|v1|v2|r2|nindels|data|full) GROUPS_REQ+=("$a") ;;
    *) echo "unknown argument '$a' (try --help)"; exit 2 ;;
  esac
done
if [ ${#GROUPS_REQ[@]} = 0 ]; then GROUPS_REQ=(v1 v2); fi
want() { local g; for g in "${GROUPS_REQ[@]}"; do [ "$g" = "$1" ] && return 0; done; return 1; }
if want full; then GROUPS_REQ+=(v1 v2 r2 nindels); fi
if [ $WITH_DATA = 1 ]; then GROUPS_REQ+=(data); fi
PY_NEEDED=0; for g in python v1 v2; do want $g && PY_NEEDED=1; done
IS_FULL=0; want full && IS_FULL=1

# ---- result bookkeeping -----------------------------------------------------------------------
SUMMARY=()   # "group|passed|failed|skipped|note"
record() { SUMMARY+=("$1|$2|$3|$4|${5:-}"); }

# ---- install hints ----------------------------------------------------------------------------
OS=$(uname -s); ARCH=$(uname -m); DISTRO=""
if [ "$OS" = Linux ] && [ -r /etc/os-release ]; then DISTRO=$(. /etc/os-release; echo "${ID:-linux} ${ID_LIKE:-}"); fi
hint() {
  echo "  Recommended on every platform: conda env create -f environment.yml && conda activate varion"
  case "$OS:$DISTRO" in
    Darwin:*) echo "  macOS: xcode-select --install; brew install htslib samtools bcftools cmake  (untested)" ;;
    Linux:*fedora*) echo "  Fedora: sudo dnf install gcc gcc-c++ make cmake python3-devel htslib htslib-devel samtools bcftools curl" ;;
    Linux:*ubuntu*|Linux:*debian*) echo "  Ubuntu/Debian: sudo apt install build-essential cmake python3-dev python3-venv libhts-dev samtools bcftools curl"
                                   echo "  (distro htslib is older than $REF_HTS_MINOR on current LTS releases; the conda env is preferred)" ;;
    Linux:*arch*|Linux:*cachyos*) echo "  Arch/CachyOS: sudo pacman -S base-devel cmake samtools bcftools curl; htslib may need the AUR/conda" ;;
  esac
}
need() { command -v "$1" >/dev/null 2>&1 && ok "$1 found" || { bad "$1 not found ($2)"; MISSING=1; }; }
opt()  { command -v "$1" >/dev/null 2>&1 && ok "$1 found (optional)" || skip "$1 not found (optional: $2)"; }

# ---- 0. environment diagnostics ---------------------------------------------------------------
HTS_OK=0
env_diag() {
  step "0. Environment"
  MISSING=0
  echo "platform : $OS $ARCH ${DISTRO:+($DISTRO)}"
  [ "$OS" = Linux ] && grep -m1 PRETTY_NAME /etc/os-release 2>/dev/null | sed 's/^/os       : /'
  [ -n "${WSL_DISTRO_NAME:-}" ] && echo "WSL      : $WSL_DISTRO_NAME"
  [ -n "${CONDA_PREFIX:-}" ] && echo "conda env: $CONDA_PREFIX"
  case "$OS" in Linux|Darwin) ;; *) bad "unsupported OS '$OS' (use WSL2 Ubuntu or Docker on Windows)"; MISSING=1 ;; esac
  need python3 "install Python 3"
  command -v python3 >/dev/null 2>&1 && echo "python   : $(python3 --version 2>&1)"
  need gcc "install gcc / Xcode command line tools"; need g++ "install g++ (C++20)"; need make "install make"
  command -v g++ >/dev/null 2>&1 && echo "g++      : $(g++ --version | head -1)"
  if command -v g++ >/dev/null 2>&1; then
    echo 'int main(){return 0;}' | g++ -std=c++20 -x c++ - -o /dev/null 2>/dev/null \
      && ok "g++ accepts -std=c++20" || { bad "g++ cannot compile C++20 (need GCC >= 10 / Clang >= 13)"; MISSING=1; }
  fi
  if want nindels; then need cmake "needed to build nindels"; fi
  if want data; then
    need samtools "needed for the data group"; need bcftools "needed for the data group"; need curl "needed for the data group"
    for t in samtools bcftools; do command -v $t >/dev/null 2>&1 && echo "$t : $($t --version | head -1)"; done
  else
    opt samtools "only needed for the data group"; opt bcftools "only needed for the data group"
  fi
  if [ $PY_NEEDED = 1 ]; then
    miss=$(python3 - <<'PYEOF' 2>/dev/null
import importlib
bad = [m for m in ("numpy", "scipy", "pysam", "pytest", "torch") if importlib.util.find_spec(m) is None]
print(" ".join(bad))
PYEOF
)
    if [ -z "$miss" ]; then ok "python packages (numpy scipy pysam pytest torch)"
    else bad "python packages missing: $miss (pip install -r requirements.txt)"; MISSING=1; fi
  fi
  # htslib: needed by every native group.  Headers and library must come from the SAME prefix.
  NEED_HTS=0; for g in env native v1 v2 r2 nindels data; do want $g && NEED_HTS=1; done
  if [ $NEED_HTS = 1 ]; then
    HTS_ERR=$(mktemp); source "$ROOT/scripts/hts_env.sh"; hts_resolve
    echo "htslib   : ${HTS_VERSION:-unknown} via $HTS_SOURCE"
    echo "           include=$HTS_INC  link=$HTS_LIB"
    if hts_check 2>"$HTS_ERR"; then
      ok "htslib headers + library found"; HTS_OK=1
      case "${HTS_VERSION:-}" in
        $REF_HTS_MINOR*) ok "htslib $HTS_VERSION matches the series used for the published results ($REF_HTS_MINOR.x)" ;;
        *) warn "htslib ${HTS_VERSION:-unknown} differs from the $REF_HTS_MINOR.x used for the published results; the equivalence tests below decide whether outputs are unchanged" ;;
      esac
      PYHTS=$(python3 -c 'import pysam.version as v; print(v.__htslib_version__)' 2>/dev/null)
      [ -n "$PYHTS" ] && echo "pysam    : bundles its own htslib $PYHTS (used by Python code only; independent of the native libhts above)"
    else
      bad "htslib not usable"; sed "s/^/      /" "$HTS_ERR"; MISSING=1
    fi
    rm -f "$HTS_ERR"
  fi
  if [ $MISSING = 1 ]; then echo; echo "Missing prerequisites:"; hint; echo; echo "Fix the items above and re-run."; print_summary; exit 1; fi
}

# ---- builders (memoised) ----------------------------------------------------------------------
B_V1=""; B_V2X=""; B_R2=""; B_NIN=""
build_v1() { [ -n "$B_V1" ] && return $B_V1; bash native/build.sh > "$LOG/build_v1.log" 2>&1; B_V1=$?
  [ $B_V1 = 0 ] && ok "V1 native backends (native/build.sh)" || bad "V1 native build, see $LOG/build_v1.log"; return $B_V1; }
build_v2x() { [ -n "$B_V2X" ] && return $B_V2X; make -C native_v2x engine kernels > "$LOG/build_v2.log" 2>&1; B_V2X=$?
  [ $B_V2X = 0 ] && ok "V2/V2.x engine native_v2x/build/dnav2" || bad "V2 engine build, see $LOG/build_v2.log"
  [ $B_V2X = 0 ] && runtime_check native_v2x/build/dnav2; return $B_V2X; }
build_r2() { [ -n "$B_R2" ] && return $B_R2; make -C native_r2 engine > "$LOG/build_r2.log" 2>&1; B_R2=$?
  [ $B_R2 = 0 ] && ok "R1.1/R2 engine native_r2/build/dnav2" || bad "R2 engine build, see $LOG/build_r2.log"
  [ $B_R2 = 0 ] && runtime_check native_r2/build/dnav2; return $B_R2; }
build_nin() { [ -n "$B_NIN" ] && return $B_NIN
  { cmake -S nindels -B nindels/build -DCMAKE_BUILD_TYPE=RelWithDebInfo ${HTS_PREFIX:+-DHTSLIB_ROOT="$HTS_PREFIX"} && cmake --build nindels/build -j; } > "$LOG/build_nindels.log" 2>&1; B_NIN=$?
  [ $B_NIN = 0 ] && ok "nindels (nindels/build/nindels-call, nindels-tests)" || bad "nindels build, see $LOG/build_nindels.log"
  if [ $B_NIN = 0 ] && grep -q 'built without BAM/VCF' "$LOG/build_nindels.log"; then
    bad "nindels was built WITHOUT htslib: 'call' and 'estimate-str' are missing (see $LOG/build_nindels.log)"; B_NIN=1; fi
  return $B_NIN; }
# The binary must start: catches an unresolved libhts at run time (rpath / LD_LIBRARY_PATH problems).
runtime_check() { local out; out=$("$1" 2>&1 </dev/null | head -3)
  if printf '%s' "$out" | grep -qiE 'error while loading|cannot open shared|image not found|dyld'; then
    bad "$1 cannot load libhts at run time: $out"; else ok "$1 starts (libhts resolved at run time)"; fi; }

# ---- result helpers ---------------------------------------------------------------------------
print_summary() {
  echo; echo "=== Summary ==="
  printf '%-10s %8s %8s %8s  %s\n' group passed failed skipped note
  local row g p f s n tp=0 tf=0 ts=0
  for row in "${SUMMARY[@]+"${SUMMARY[@]}"}"; do
    IFS='|' read -r g p f s n <<<"$row"; printf '%-10s %8s %8s %8s  %s\n' "$g" "$p" "$f" "$s" "$n"
    tp=$((tp+p)); tf=$((tf+f)); ts=$((ts+s))
  done
  printf '%-10s %8s %8s %8s\n' TOTAL $tp $tf $ts
}

# ---- groups -----------------------------------------------------------------------------------
V1_TESTS="test_cascade.py test_cheap_router.py test_native_pileup_equivalence.py test_native_reads_equivalence.py test_hg005_extraction_integrity.py test_counts_backend_fuzz.py test_bench_v12.py test_robustness_benchmark.py"
pytest_counts() { # echo "passed failed skipped" from a pytest log
  python3 - "$1" <<'PYEOF'
import re, sys
t = open(sys.argv[1]).read().strip().splitlines()[-1]
n = lambda k: int(next(iter(re.findall(r"(\d+) " + k, t)), 0))
print(n("passed") + n("xpassed"), n("failed") + n("error"), n("skipped"))
PYEOF
}
group_python() {
  step "Python-only tests (V1 suite, native backends disabled via AI_DNA_ANALYZER_DISABLE_NATIVE_PILEUP=1)"
  AI_DNA_ANALYZER_DISABLE_NATIVE_PILEUP=1 python3 -m pytest -q -p no:cacheprovider $V1_TESTS > "$LOG/test_python.log" 2>&1; local rc=$?
  tail -1 "$LOG/test_python.log"; read -r p f s <<<"$(pytest_counts "$LOG/test_python.log")"
  [ $rc = 0 ] && ok "Python-only tests ($p passed, $s skipped because they need native backends)" || bad "Python-only tests, see $LOG/test_python.log"
  record python "$p" "$f" "$s" "skips need native backends: run group v1"
}
group_native() {
  step "Native build (no tests)"
  build_v1; build_v2x; build_r2; build_nin
  record native 4 0 0 "compile only"
}
group_v1() {
  step "1. V1 native backends + V1 test suite (frozen cascade)"
  build_v1 || { record v1 0 1 0 "build failed"; return; }
  python3 -m pytest -q -p no:cacheprovider $V1_TESTS > "$LOG/test_v1.log" 2>&1; local rc=$?
  tail -1 "$LOG/test_v1.log"; read -r p f s <<<"$(pytest_counts "$LOG/test_v1.log")"
  if [ $rc != 0 ]; then bad "V1 tests, see $LOG/test_v1.log"
  elif [ "$s" -gt $EXPECTED_V1_SKIPS ]; then bad "V1 tests: $s skipped, only $EXPECTED_V1_SKIPS documented skips are allowed (native backend not loaded?), see $LOG/test_v1.log"
  else ok "V1 tests ($p passed, $s skipped = documented; xfailed = documented adversarial cases)"; fi
  record v1 "$p" "$f" "$s" "$(tail -1 "$LOG/test_v1.log" | sed 's/ in [0-9.]*s.*//')"
}
group_v2() {
  step "2. V2/V2.x engine + kernels bit-identical to the frozen V1 NumPy/SciPy code"
  build_v2x || { record v2 0 1 0 "build failed"; return; }
  make -s -C native_v2x test > "$LOG/test_v2.log" 2>&1; local rc=$?
  grep -E '^\[(PASS|FAIL)\]' "$LOG/test_v2.log" | cut -c1-90
  local p f; p=$(grep -c '^\[PASS\]' "$LOG/test_v2.log"); f=$(grep -c '^\[FAIL\]' "$LOG/test_v2.log")
  if [ $rc = 0 ] && [ "$f" = 0 ] && [ "$p" -gt 0 ]; then ok "V2 kernel tests ($p passed)"; else bad "V2 kernel tests, see $LOG/test_v2.log"; [ "$f" = 0 ] && f=1; fi
  record v2 "$p" "$f" 0
}
group_r2() {
  step "3. R1.1/R2 engine: tail tests and call-domain (frame) tests"
  build_r2 || { record r2 0 1 0 "build failed"; return; }
  make -s -C native_r2 test-tail > "$LOG/test_r2_tail.log" 2>&1; local rc=$?
  local line; line=$(grep -E '[0-9]+ checks, [0-9]+ failures' "$LOG/test_r2_tail.log" | tail -1); echo "$line"
  local checks fails; checks=$(sed -n 's/^\([0-9]*\) checks.*/\1/p' <<<"$line"); fails=$(sed -n 's/.* \([0-9]*\) failures.*/\1/p' <<<"$line")
  local pass=0 fail=0
  if [ $rc = 0 ] && [ -n "$checks" ] && [ "${fails:-1}" = 0 ]; then ok "R2 tail tests ($checks checks, 0 failures)"; pass=$((pass+checks)); else bad "R2 tail tests, see $LOG/test_r2_tail.log"; fail=$((fail+${fails:-1})); fi
  make -s -C native_r2 test-frame > "$LOG/test_r2_frame.log" 2>&1; rc=$?
  if [ $rc = 0 ] && grep -q 'FRAME TESTS PASSED' "$LOG/test_r2_frame.log"; then ok "R2 frame tests (FRAME TESTS PASSED)"; pass=$((pass+1)); else bad "R2 frame tests, see $LOG/test_r2_frame.log"; fail=$((fail+1)); fi
  record r2 "$pass" "$fail" 0 "tail = $checks checks; frame = 1 suite (no per-check count printed)"
}
group_nindels() {
  step "4. nindels unit tests"
  build_nin || { record nindels 0 1 0 "build failed"; return; }
  nindels/build/nindels-tests > "$LOG/test_nindels.log" 2>&1; local rc=$?
  local line; line=$(grep -E '^[0-9]+ passed, [0-9]+ failed' "$LOG/test_nindels.log" | tail -1); echo "$line"
  local p f; p=$(sed -n 's/^\([0-9]*\) passed.*/\1/p' <<<"$line"); f=$(sed -n 's/.* \([0-9]*\) failed.*/\1/p' <<<"$line")
  if [ $rc = 0 ] && [ -n "$p" ] && [ "${f:-1}" = 0 ]; then ok "nindels tests ($p passed, 0 failed)"; else bad "nindels tests, see $LOG/test_nindels.log"; fi
  record nindels "${p:-0}" "${f:-1}" 0
}
group_data() {
  step "5. Real data: V2.x A+C on HG003 chr8:20,000,001-20,249,984 (300x) vs archived VCF"
  build_v2x || { record data 0 1 0 "needs V2.x engine"; return; }
  local before=$FAILED; FAILED=0
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
  local after=$FAILED; FAILED=$((before | after))
  if [ $after = 0 ]; then record data 1 0 0 "calls identical to archive (calls, not accuracy)"; else record data 0 1 0 "see messages above"; fi
}

# ---- run --------------------------------------------------------------------------------------
env_diag
want python  && group_python
want native  && group_native
want v1      && group_v1
want v2      && group_v2
want r2      && group_r2
want nindels && group_nindels
want data    && group_data

print_summary
NOTRUN=""
for g in v1 v2 r2 nindels; do want $g || NOTRUN="$NOTRUN $g"; done
want data || NOTRUN="$NOTRUN data(optional: --data)"
echo
if [ $FAILED != 0 ]; then
  echo "SOME CHECKS FAILED"
elif [ $IS_FULL = 1 ]; then
  if want data; then echo "COMPLETE VERIFICATION PASSED (all offline groups and the real-data test)"
  else echo "COMPLETE OFFLINE VERIFICATION PASSED (the optional real-data test was not run)"; fi
else
  echo "ALL SELECTED CHECKS PASSED (groups: ${GROUPS_REQ[*]})"
  echo "PARTIAL: this is NOT a complete verification. Not run:${NOTRUN:- none of the offline groups; use 'full'}. Run: bash verify.sh full"
fi
exit $FAILED
