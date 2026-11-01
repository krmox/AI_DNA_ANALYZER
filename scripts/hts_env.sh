#!/usr/bin/env bash
# Single source of truth for locating htslib (headers + shared library).
# Used by native/build.sh, verify.sh and (through mk/htslib.mk) every Makefile.
#
# Resolution order (first match wins):
#   1. HTS_PREFIX            explicit install prefix (include/, lib/)
#   2. CONDA_PREFIX          active conda/mamba env that contains include/htslib/sam.h
#   3. pkg-config htslib
#   4. legacy fallback       pysam's bundled headers + system libhts in a standard lib dir
#                            (the original Fedora layout)
# Individual pieces can still be overridden: HTS_INC (dir containing htslib/),
# HTS_LIBDIR (dir containing libhts), HTS_LIB (complete linker flags), HTS_RPATH.
#
# Usage:  scripts/hts_env.sh {inc|libdir|lib|rpath|version|source|report|check}
# Pure bash + grep/sed; no network access, nothing is installed or modified.

_hts_has_hdr()  { [ -f "$1/htslib/sam.h" ]; }
_hts_has_lib()  { local d=$1 f; for f in libhts.so libhts.so.3 libhts.dylib libhts.a; do [ -e "$d/$f" ] && return 0; done; ls "$d"/libhts.*.dylib >/dev/null 2>&1; }

_hts_resolve() {
  [ -n "${_HTS_DONE:-}" ] && return 0
  local inc="" libdir="" src="" d
  if [ -n "${HTS_PREFIX:-}" ]; then
    inc="$HTS_PREFIX/include"; src="HTS_PREFIX=$HTS_PREFIX"
    for d in "$HTS_PREFIX/lib" "$HTS_PREFIX/lib64"; do _hts_has_lib "$d" && { libdir=$d; break; }; done
    [ -z "$libdir" ] && libdir="$HTS_PREFIX/lib"
  elif [ -n "${CONDA_PREFIX:-}" ] && _hts_has_hdr "$CONDA_PREFIX/include"; then
    inc="$CONDA_PREFIX/include"; libdir="$CONDA_PREFIX/lib"; src="conda env $CONDA_PREFIX"
  elif command -v pkg-config >/dev/null 2>&1 && pkg-config --exists htslib 2>/dev/null; then
    inc=$(pkg-config --variable=includedir htslib); libdir=$(pkg-config --variable=libdir htslib); src="pkg-config"
  else
    src="legacy fallback (pysam headers + system libhts)"
    inc=$(python3 -c 'import pysam; print(pysam.get_include()[1])' 2>/dev/null || true)
    for d in /usr/lib64 /usr/lib/x86_64-linux-gnu /usr/lib/aarch64-linux-gnu /usr/lib /usr/local/lib /opt/homebrew/lib; do
      _hts_has_lib "$d" && { libdir=$d; break; }
    done
    [ -z "$libdir" ] && libdir=/usr/lib64
  fi
  HTS_INC=${HTS_INC:-$inc}; HTS_LIBDIR=${HTS_LIBDIR:-$libdir}; HTS_SOURCE=$src
  if [ -z "${HTS_LIB:-}" ]; then
    # Prefer the portable "-lhts" when the development symlink exists; fall back to the
    # versioned soname (GNU ld only) for runtime-only installs such as Fedora's htslib.
    if [ -e "$HTS_LIBDIR/libhts.so" ] || [ -e "$HTS_LIBDIR/libhts.dylib" ] || [ -e "$HTS_LIBDIR/libhts.a" ]; then
      HTS_LIB="-L$HTS_LIBDIR -lhts"
    elif [ -e "$HTS_LIBDIR/libhts.so.3" ]; then
      HTS_LIB="-L$HTS_LIBDIR -l:libhts.so.3"
    else
      HTS_LIB="-L$HTS_LIBDIR -lhts"
    fi
  fi
  HTS_RPATH=${HTS_RPATH--Wl,-rpath,$HTS_LIBDIR}
  HTS_VERSION=
  if [ -f "$HTS_INC/htslib/hts.h" ]; then
    # hts.h encodes the version as MMmmpp (e.g. 102301 = 1.23.1)
    local n; n=$(sed -n 's/^#define HTS_VERSION \([0-9][0-9]*\).*/\1/p' "$HTS_INC/htslib/hts.h" | head -1)
    if [ -n "$n" ]; then
      HTS_VERSION="$((n / 100000)).$(((n / 100) % 1000))"
      [ $((n % 100)) -ne 0 ] && HTS_VERSION="$HTS_VERSION.$((n % 100))"
    fi
  fi
  _HTS_DONE=1
}

hts_resolve() { _hts_resolve; export HTS_INC HTS_LIBDIR HTS_LIB HTS_RPATH HTS_VERSION HTS_SOURCE; }

# 0 = headers and library found; otherwise prints what is missing and returns 1.
hts_check() {
  _hts_resolve; local rc=0
  if ! _hts_has_hdr "$HTS_INC"; then
    echo "htslib headers not found: '${HTS_INC:-<empty>}/htslib/sam.h' does not exist." >&2; rc=1
  fi
  if ! _hts_has_lib "$HTS_LIBDIR"; then
    echo "htslib library not found: no libhts.{so,so.3,dylib,a} in '$HTS_LIBDIR'." >&2; rc=1
  fi
  if [ $rc = 1 ]; then
    cat >&2 <<MSG
Resolved via: $HTS_SOURCE
Fix one of:
  conda env create -f environment.yml && conda activate varion       (recommended, all platforms)
  export HTS_PREFIX=/path/to/htslib/install                          (prefix with include/ and lib/)
  make ... HTS_INC=/path/include HTS_LIB="-L/path/lib -lhts"         (explicit)
  Fedora: dnf install htslib htslib-devel | Debian/Ubuntu: apt install libhts-dev | macOS: brew install htslib
See docs/PLATFORM_SUPPORT.md.
MSG
  fi
  return $rc
}

# Run as a command (not when sourced).
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  _hts_resolve
  case "${1:-report}" in
    inc) echo "$HTS_INC" ;;
    libdir) echo "$HTS_LIBDIR" ;;
    lib) echo "$HTS_LIB" ;;
    rpath) echo "$HTS_RPATH" ;;
    version) echo "$HTS_VERSION" ;;
    source) echo "$HTS_SOURCE" ;;
    check) hts_check ;;
    report) printf 'source : %s\ninclude: %s\nlibdir : %s\nlink   : %s\nrpath  : %s\nversion: %s (headers)\n' \
              "$HTS_SOURCE" "$HTS_INC" "$HTS_LIBDIR" "$HTS_LIB" "$HTS_RPATH" "${HTS_VERSION:-unknown}" ;;
    *) echo "usage: $0 {inc|libdir|lib|rpath|version|source|report|check}" >&2; exit 2 ;;
  esac
fi
