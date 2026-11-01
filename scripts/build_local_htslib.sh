#!/usr/bin/env bash
# OPTIONAL (Option B), UNTESTED: build htslib 1.23.1 from the upstream source tarball into a
# project-local prefix.  Nothing is downloaded or compiled unless YOU run this script.
#
#   bash scripts/build_local_htslib.sh [PREFIX]        default PREFIX = ./third_party/htslib
#   export HTS_PREFIX=$PWD/third_party/htslib          then: bash verify.sh full
#
# Source  : https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2
# SHA-256 : f8a3f36effeec38f043c53ab1f2d9ed45064f14205c5ef8e3c815763b90803c4 (checked before unpacking;
#           this value was computed from a download of that URL, the build itself has not been run)
# License : htslib is MIT/Expat (bundled CRAM code BSD-style); see COPYING in the tarball.
# Needs   : C compiler, make, curl, and the zlib, bzip2, xz (liblzma) and libcurl development headers
#           (Debian/Ubuntu: zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev;
#            Fedora: zlib-devel bzip2-devel xz-devel libcurl-devel).
# Linking : the Makefiles/CMake add an rpath to $PREFIX/lib, so no LD_LIBRARY_PATH is needed.
set -euo pipefail
cd "$(dirname "$0")/.."
VER=1.23.1
SHA=f8a3f36effeec38f043c53ab1f2d9ed45064f14205c5ef8e3c815763b90803c4
URL=https://github.com/samtools/htslib/releases/download/$VER/htslib-$VER.tar.bz2
PREFIX=$(mkdir -p "${1:-third_party/htslib}" && cd "${1:-third_party/htslib}" && pwd)
WORK=$(mktemp -d); trap 'rm -rf "$WORK"' EXIT
echo "Downloading $URL"
curl -fsSL --retry 3 -o "$WORK/htslib.tar.bz2" "$URL"
if command -v sha256sum >/dev/null; then got=$(sha256sum "$WORK/htslib.tar.bz2" | cut -d' ' -f1); else got=$(shasum -a 256 "$WORK/htslib.tar.bz2" | cut -d' ' -f1); fi
[ "$got" = "$SHA" ] || { echo "SHA-256 mismatch: got $got, expected $SHA" >&2; exit 1; }
tar -xjf "$WORK/htslib.tar.bz2" -C "$WORK"
cd "$WORK/htslib-$VER"
JOBS=$( (nproc || sysctl -n hw.ncpu) 2>/dev/null || echo 2)
./configure --prefix="$PREFIX" > "$WORK/configure.log" 2>&1 || { tail -20 "$WORK/configure.log" >&2; echo "configure failed (missing dev headers? see the list at the top of this script)" >&2; exit 1; }
make -j"$JOBS" > "$WORK/make.log" 2>&1 || { tail -20 "$WORK/make.log" >&2; exit 1; }
make install > "$WORK/install.log" 2>&1
echo "htslib $VER installed in $PREFIX"
echo "Next: export HTS_PREFIX=$PREFIX && bash scripts/hts_env.sh report"
