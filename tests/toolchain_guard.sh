#!/usr/bin/env bash
# tests/toolchain_guard.sh - toolchain hint for tests that compare bitwise against toolchain-tied data (task P10.4, D-044).
#
# Usage:
#   toolchain_guard.sh --fortran-ref TOOLCHAIN.txt -- COMMAND [ARGS...]
#       CTest TEST_LAUNCHER of the L1 tests: runs COMMAND; if it fails (exit code not 0 and not 77 = skip), prints
#       the toolchain of the Fortran references (reference/TOOLCHAIN.txt: glibc, libm md5) next to this machine's
#       (and the compiler that built COMMAND), and exits with COMMAND's code.
#   toolchain_guard.sh --baseline-hint README.txt [--binary SIMCPP]
#       prints the same hint for a C++ regression baseline (its README.txt "toolchain:" line) and exits 0.
# The hint says DIFFERENT when the glibc version, the libm md5 or (baselines) the compiler that built simc++ differ:
# then the failure is expected and the references must be regenerated (D-044); SAME means the toolchain is not the
# explanation.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FP="$HERE/../tools/toolchain_fingerprint.sh"
MODE=""; REF=""; BIN=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fortran-ref) MODE=fortran; REF="$2"; shift 2 ;;
    --baseline-hint) MODE=baseline; REF="$2"; shift 2 ;;
    --binary) BIN="$2"; shift 2 ;;
    --) shift; break ;;
    *) break ;;
  esac
done

hint() {   # $1 = binary under test (may be empty)
  local here_glibc here_libm here_built ref_glibc="?" ref_libm="?" ref_built="" ref_desc verdict="SAME"
  here_glibc="$("$FP" --field glibc)"; here_libm="$("$FP" --field libm_md5)"
  [[ -n "$1" ]] && here_built="$("$FP" --binary "$1" --field built_by)" || here_built=""
  if [[ ! -f "$REF" ]]; then
    echo "TOOLCHAIN HINT: no toolchain record $REF; this machine has $("$FP" ${1:+--binary "$1"}). If the failure is a"
    echo "  ~1 ulp difference, regenerate the bitwise references on this toolchain per D-044."
    return
  fi
  if [[ $MODE == fortran ]]; then
    ref_glibc="$(grep -m1 -oE '^glibc/libm: *glibc-[0-9.]+' "$REF" | grep -oE '[0-9]+\.[0-9]+$')"
    ref_libm="$(grep -m1 -E '^libm\.so\.6:' "$REF" | grep -oE 'md5 [0-9a-f]{32}' | awk '{print $2}')"
    ref_desc="the Fortran references ($REF)"
  else
    local line; line="$(grep -m1 -E '^toolchain:' "$REF")"
    if [[ -z "$line" ]]; then
      echo "TOOLCHAIN HINT (DIFFERENT?): $REF has no 'toolchain:' line (baseline made before P10.4, i.e. not on a recorded"
      echo "  toolchain); this machine has $("$FP" ${1:+--binary "$1"}). Regenerate the baseline per D-044."
      return
    fi
    ref_glibc="$(grep -oE 'glibc [0-9.]+' <<<"$line" | awk '{print $2}')"
    ref_libm="$(grep -oE 'libm md5 [0-9a-f]{32}' <<<"$line" | awk '{print $3}')"
    ref_built="$(grep -oE 'built by GCC [0-9.]+( [0-9]{8})?' <<<"$line" | sed 's/built by GCC //')"
    ref_desc="the baseline ($(dirname "$REF"))"
  fi
  [[ "$ref_glibc" == "$here_glibc" && "$ref_libm" == "$here_libm" ]] || verdict="DIFFERENT"
  [[ -n "$ref_built" && -n "$here_built" && "$ref_built" != "$here_built" ]] && verdict="DIFFERENT"
  echo "TOOLCHAIN HINT ($verdict): references made on glibc ${ref_glibc:-?} (libm md5 ${ref_libm:-?})${ref_built:+, simc++ built by GCC $ref_built}"
  echo "  [$ref_desc]; this machine has glibc $here_glibc (libm md5 $here_libm)${here_built:+, $(basename "$1") built by GCC $here_built}."
  if [[ $verdict == DIFFERENT ]]; then
    echo "  The toolchain differs: ~1 ulp differences are expected; regenerate the references per D-044"
    echo "  (Fortran: reference/TOOLCHAIN.txt, task P10.F2; C++ baselines: <commit>_glibc<ver>_* folders, task P10.4)."
  else
    echo "  Same glibc/libm${ref_built:+ and compiler} as the references: the failure is NOT explained by a toolchain change (D-044)."
  fi
}

if [[ $MODE == baseline ]]; then hint "$BIN"; exit 0; fi
if [[ $MODE != fortran || $# -eq 0 ]]; then sed -n '2,14p' "$0"; exit 2; fi
"$@"; rc=$?
if [[ $rc -ne 0 && $rc -ne 77 ]]; then
  exe="$1"; [[ -f "$exe" ]] || exe="$(command -v "$1" 2>/dev/null || true)"
  hint "$exe"
fi
exit $rc
