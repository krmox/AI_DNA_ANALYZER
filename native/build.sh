#!/bin/bash
# Builds pileup_native (counts) and reads_native (read tensor), the optional C/htslib accelerations.
# Uses pysam's bundled htslib headers and the system libhts.so.3; needs gcc and Python.h.
set -euo pipefail
cd "$(dirname "$0")"
SOEXT="$(python3 -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')"
PYINC="$(python3 -c 'import sysconfig; print(sysconfig.get_paths()["include"])')"
[ -f "$PYINC/Python.h" ] || PYINC="$HOME/.cache/pydevel/usr/include/python3.14"
[ -f "$PYINC/Python.h" ] || { echo "Python.h not found in $PYINC (install python3-devel)"; exit 1; }
HTS_INCLUDE_DIR="$(python3 -c 'import pysam; print(pysam.get_include()[1])')"
for m in pileup_native reads_native; do
  echo "Building $m$SOEXT ..."
  gcc -O3 -march=native -fPIC -shared -I "$PYINC" -I /usr/include/python3.14 -I "$HTS_INCLUDE_DIR" -o "$m$SOEXT" "$m.c" \
      -L/usr/lib64 -l:libhts.so.3 -lm -lpthread
done
python3 -c "
import sys; sys.path.insert(0,'.')
import pileup_native, reads_native; print('imports OK')"
