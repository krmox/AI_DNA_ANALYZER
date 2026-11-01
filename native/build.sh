#!/bin/bash
# Builds pileup_native (counts) and reads_native (read tensor), the optional C/htslib accelerations.
# Needs a C compiler, Python.h and htslib (headers + shared library).  htslib is located by
# scripts/hts_env.sh: HTS_PREFIX, then the active conda env, then pkg-config, then the legacy
# Fedora layout.  Override with HTS_INC (dir containing htslib/), HTS_LIB (linker flags) and
# HTS_RPATH.  The compiled modules embed an rpath to the htslib directory, so no
# LD_LIBRARY_PATH is needed.  Run `bash scripts/hts_env.sh report` to see what is used.
set -euo pipefail
cd "$(dirname "$0")"
source ../scripts/hts_env.sh
hts_resolve
hts_check || exit 1
echo "htslib ${HTS_VERSION:-?} via $HTS_SOURCE ($HTS_INC ; $HTS_LIB)"

SOEXT="$(python3 -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')"
PYINC="$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["include"])')"
PYVER="$(python3 -c 'import sys; print("%d.%d" % sys.version_info[:2])')"
[ -f "$PYINC/Python.h" ] || PYINC="$HOME/.cache/pydevel/usr/include/python$PYVER"
[ -f "$PYINC/Python.h" ] || { echo "Python.h not found in $PYINC (install python3-devel / python3-dev, or use the conda env)"; exit 1; }
CC="${CC:-gcc}"
# -march=native is the original setting; set NATIVE_ARCH_FLAG= (empty) to build for generic CPUs.
ARCHFLAG="${NATIVE_ARCH_FLAG--march=native}"
EXTRA=()
[ "$(uname -s)" = Darwin ] && EXTRA=(-undefined dynamic_lookup)
for m in pileup_native reads_native; do
  echo "Building $m$SOEXT ..."
  # shellcheck disable=SC2086
  "$CC" -O3 $ARCHFLAG -fPIC -shared -I "$PYINC" -I "/usr/include/python$PYVER" -I "$HTS_INC" \
      -o "$m$SOEXT" "$m.c" $HTS_LIB $HTS_RPATH -lm -lpthread ${EXTRA[@]+"${EXTRA[@]}"}
done
python3 -c "
import sys; sys.path.insert(0,'.')
import pileup_native, reads_native; print('imports OK')"
