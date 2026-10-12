#!/usr/bin/env bash
# toolchain_fingerprint.sh - one-line fingerprint of the toolchain that bitwise results depend on (task P10.4, D-044).
#
# Usage: toolchain_fingerprint.sh [--binary FILE] [--field glibc|libm_md5|gcc|gfortran|root|built_by]
#   (no option)    print "glibc 2.43 (libm md5 <md5>); gcc <ver>; gfortran <ver>; ROOT <ver> (<ROOTSYS>)"
#                  ROOT = $ROOTSYS/bin/root-config, default /home/ayerbe/Builds/ROOT6_install (P10.5; it only affects tree I/O)
#   --binary FILE  append "; <FILE basename> built by GCC <ver>" (from the ELF .comment section)
#   --field NAME   print only that field (built_by needs --binary)
# Bitwise references (Fortran tables/captures in reference/, C++ regression baselines) are only valid for the
# toolchain they were made with: libm results can change by ~1 ulp between glibc versions and the compiler can
# change code generation (gcc 16 no longer folds pow(x,2.) into x*x). After an OS compiler/libc update,
# regenerate them (D-044).
set -uo pipefail
BIN=""; FIELD=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary) BIN="$2"; shift ;;
    --field) FIELD="$2"; shift ;;
    -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
  esac; shift
done
glibc="$(getconf GNU_LIBC_VERSION 2>/dev/null | awk '{print $2}')"
libm=""
for p in /usr/lib64/libm.so.6 /lib64/libm.so.6 /lib/x86_64-linux-gnu/libm.so.6 /usr/lib/libm.so.6; do
  [[ -f "$p" ]] && { libm="$(md5sum "$(readlink -f "$p")" | awk '{print $1}')"; break; }
done
ver() { command -v "$1" >/dev/null 2>&1 && "$1" --version 2>/dev/null | head -1 | sed -E 's/^.*\) *([0-9][0-9.]* [0-9]{8}).*$/\1/;s/^.*\) *([0-9][0-9.]*).*$/\1/'; }
gcc_v="$(ver gcc)"; gf_v="$(ver gfortran)"
rootsys="${ROOTSYS:-/home/ayerbe/Builds/ROOT6_install}"
root_v="$("$rootsys/bin/root-config" --version 2>/dev/null)"
built=""
if [[ -n "$BIN" && -f "$BIN" ]]; then
  built="$(readelf -p .comment "$BIN" 2>/dev/null | grep -m1 -oE 'GCC: \([^)]*\) [0-9][0-9.]*( [0-9]{8})?' | sed -E 's/^GCC: \([^)]*\) //')"
  [[ -n "$built" ]] || built="unknown"
fi
case "$FIELD" in
  glibc) echo "$glibc" ;;
  libm_md5) echo "$libm" ;;
  gcc) echo "$gcc_v" ;;
  gfortran) echo "$gf_v" ;;
  root) echo "$root_v" ;;
  built_by) echo "$built" ;;
  "") line="glibc ${glibc:-?} (libm md5 ${libm:-?}); gcc ${gcc_v:-none}; gfortran ${gf_v:-none}; ROOT ${root_v:-none} ($rootsys)"
      [[ -n "$BIN" ]] && line="$line; $(basename "$BIN") built by GCC ${built}"
      echo "$line" ;;
  *) echo "unknown field $FIELD" >&2; exit 2 ;;
esac
