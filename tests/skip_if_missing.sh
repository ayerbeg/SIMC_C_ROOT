#!/usr/bin/env bash
# tests/skip_if_missing.sh — run a test command only if its input files exist (task P15.D1, D-073).
# usage: skip_if_missing.sh FILE [FILE...] -- COMMAND [ARGS...]
# Exit 77 (CTest SKIP, with SKIP_RETURN_CODE 77) if a FILE is missing (e.g. SIMC_REFERENCE_DIR / SIMC_FORTRAN_DIR not
# configured on this machine), else the exit status of COMMAND.
files=()
while [[ $# -gt 0 && $1 != "--" ]]; do files+=("$1"); shift; done
[[ ${1:-} == "--" ]] && shift
for f in "${files[@]}"; do
  [[ -e $f ]] || { echo "SKIP: $f not found"; exit 77; }
done
exec "$@"
