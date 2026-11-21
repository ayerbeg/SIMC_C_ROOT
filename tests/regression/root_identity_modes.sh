#!/usr/bin/env bash
# tests/regression/root_identity_modes.sh — the U.3 modes of simc_root_identity on Fortran reference files only (fast):
#   --bin   a make_root_tree tree against the .bin it was made from must be IDENTICAL in every column (the .bin reader and
#           the real*8 -> real*4 rounding); the A5 tree (1 entry of 200, D-083) against its .bin must be DIFFERENT (event count);
#   --first the 182-entry partial tree of u59_sidis_d_k1_pip_fermi_rad_A against itself: --first 182 IDENTICAL, --first 100
#           DIFFERENT (A must hold exactly N entries).
# usage: root_identity_modes.sh <simc_root_identity> <reference runs dir> <workdir>; exit 0 pass, 1 fail, 77 reference missing.
set -u
id=$1; ref=$2; work=$3
fe=$ref/u59_shipped_test_eep_fe/test_eep_fe; pz=$ref/u59_p15f2_poltar_pi0/p15f2_poltar_pi0; sf=$ref/u59_sidis_d_k1_pip_fermi_rad_A/sidis_d_k1_pip_fermi_rad_A
for f in $fe.root $fe.bin $pz.root $pz.bin $sf.root; do [[ -f $f ]] || { echo "reference $f missing: SKIP"; exit 77; }; done
mkdir -p "$work"; fail=0
expect() {   # want_exit label cmd...
  local want=$1 lab=$2; shift 2
  "$@" > "$work/$lab.txt" 2>&1; local rc=$?
  if [[ $rc == "$want" ]]; then echo "ok   $lab (exit $rc): $(tail -1 "$work/$lab.txt")"; else echo "FAIL $lab (exit $rc, want $want)"; cat "$work/$lab.txt"; fail=1; fi
}
expect 0 bin_eep_fe "$id" --bin $fe.root $fe.bin 1 999
expect 1 bin_a5_truncated "$id" --bin $pz.root $pz.bin 1 80
expect 0 first_182 "$id" --first 182 $sf.root $sf.root
expect 1 first_100 "$id" --first 100 $sf.root $sf.root
[[ $fail == 0 ]] && echo PASS || echo FAIL
exit $fail
