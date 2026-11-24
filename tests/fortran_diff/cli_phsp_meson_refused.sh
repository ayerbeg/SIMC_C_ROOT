#!/usr/bin/env bash
# tests/fortran_diff/cli_phsp_meson_refused.sh — FORTRAN_DIFFERENCES "PHYS.2C phsp mesons" (task FIX.1): doing_phsp = 1
# together with doing_pion, doing_kaon, doing_delta or doing_rho is a configuration error in both modes (exit status 2,
# "doing_phsp cannot be combined with doing_<x>"), with and without the tree (Nntu = 1 / 0): the refusal comes from the
# phase_space reaction's init(), before any sink chooses the ntuple layout. The input is the user-manual example
# h_elastic_rad.inp (hydrogen) with the two flags set.
#
# usage: cli_phsp_meson_refused.sh <simc++> <example .inp> <workdir> <flag: doing_pion|doing_kaon|doing_delta|doing_rho>
set -uo pipefail
if [[ $# -ne 4 ]]; then sed -n '2,9p' "$0"; exit 2; fi
SIMCPP=$1; SRC=$2; WORK=$3; FLAG=$4
[[ -f $SRC ]] || { echo "FAIL: $SRC missing"; exit 1; }
rm -rf "$WORK"; mkdir -p "$WORK"
fail=0
expect() {   # label, expected exit, expected message, command...
    local label=$1 want=$2 msg=$3; shift 3
    "$@" > "$WORK/$label.out" 2>&1; local rc=$?
    if [[ $rc -ne $want ]]; then echo "FAIL $label: exit $rc (expected $want)"; tail -5 "$WORK/$label.out"; fail=1; return; fi
    if ! grep -qF "$msg" "$WORK/$label.out"; then echo "FAIL $label: message not found"; tail -5 "$WORK/$label.out"; fail=1; return; fi
    echo "ok   $label: exit $rc, '$(grep -F "$msg" "$WORK/$label.out" | head -1)'"
}
for nntu in 1 0; do
    inp=$WORK/phsp_${FLAG}_nntu${nntu}.inp
    sed -e 's/^\( *doing_phsp = \)0/\11/' -e "s/^\( *${FLAG} = \)0/\11/" -e "s/^\( *Nntu = \)[0-9]*/\1${nntu}/" \
        -e 's/^\( *ngen = \)[0-9-]*/\1-100/' "$SRC" > "$inp"
    grep -q "^ *doing_phsp = 1" "$inp" && grep -q "^ *${FLAG} = 1" "$inp" && grep -q "^ *Nntu = ${nntu}" "$inp" \
        || { echo "FAIL: could not set the flags in $inp"; exit 1; }
    msg="doing_phsp cannot be combined with ${FLAG}"
    expect "default_nntu${nntu}" 2 "$msg" "$SIMCPP" "$inp" -o "$WORK/out_d${nntu}" --quiet --seed 4242
    expect "compat_nntu${nntu}" 2 "$msg" "$SIMCPP" "$inp" -o "$WORK/out_c${nntu}" --quiet --seed 4242 --fortran-compat
done
[[ $fail -eq 0 ]] && echo "PHSP + ${FLAG} REFUSED: PASS" || echo "PHSP + ${FLAG} REFUSED: FAIL"
exit $fail
