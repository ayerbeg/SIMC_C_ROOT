#!/usr/bin/env bash
# tests/fortran_diff/cli_decay_non_pik.sh — D-037 F9 (task P6.4): doing_decay = 1 with a hadron that is neither a pion
# nor a kaon (here the proton of H(e,e'p), reference/runs/eep_h_elastic_norad_A with doing_decay = 1, ctau = 1 cm).
#   default (corrected): rejected at start-up, exit status 2, message "... is neither a pion nor a kaon ...";
#   --fortran-compat   : the Fortran behaviour, the run ends at the first proton decay (project.f:87-91 `stop`;
#                        SIMC++ simc::FortranStop, exit status 3, "error, cannot decay particle with mass=  938.27").
#
# usage: cli_decay_non_pik.sh <simc++> <reference dir> <workdir>
set -uo pipefail
if [[ $# -ne 3 ]]; then sed -n '2,9p' "$0"; exit 2; fi
SIMCPP=$1; REF=$2; WORK=$3
SRC=$REF/runs/eep_h_elastic_norad_A/eep_h_elastic_norad_A.inp
[[ -f $SRC ]] || { echo "SKIP: $SRC missing"; exit 77; }
rm -rf "$WORK"; mkdir -p "$WORK"
sed -e 's/^\( *doing_decay = \)0/\11/' -e 's/^\( *ctau = \)[0-9.]*/\11.0/' -e 's/^\( *ngen = \)[0-9-]*/\1-200/' \
    "$SRC" > "$WORK/eep_decay.inp"
grep -q "doing_decay = 1" "$WORK/eep_decay.inp" || { echo "FAIL: could not set doing_decay"; exit 1; }

fail=0
expect() {   # label, expected exit, expected message, command...
    local label=$1 want=$2 msg=$3; shift 3
    "$@" > "$WORK/$label.out" 2>&1; local rc=$?
    if [[ $rc -ne $want ]]; then echo "FAIL $label: exit $rc (expected $want)"; tail -5 "$WORK/$label.out"; fail=1; return; fi
    if ! grep -qF "$msg" "$WORK/$label.out"; then echo "FAIL $label: message not found"; tail -5 "$WORK/$label.out"; fail=1; return; fi
    echo "ok   $label: exit $rc, '$(grep -F "$msg" "$WORK/$label.out" | head -1)'"
}
expect corrected 2 "is neither a pion nor a kaon" "$SIMCPP" "$WORK/eep_decay.inp" -o "$WORK/out_c" --quiet --seed 4242
expect compat 3 "error, cannot decay particle with mass=  938.27" \
    "$SIMCPP" "$WORK/eep_decay.inp" -o "$WORK/out_f" --quiet --seed 4242 --fortran-compat
[[ $fail -eq 0 ]] && echo "DECAY OF A NON-PI/K HADRON (D-037 F9): PASS" || echo "DECAY OF A NON-PI/K HADRON (D-037 F9): FAIL"
exit $fail
