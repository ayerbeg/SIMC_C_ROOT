#!/bin/bash
# tests/l1/sidis_central_log.sh — P15.1 (D-051): the SIDIS central-kinematics lines of calculate_central (simc.f:1289-1303)
# in the simc++ log must equal the Fortran run log (reference/runs/<run>/<run>.out), fortran_compat, a 5-event run.
# usage: sidis_central_log.sh <simc++> <ref_run_dir> <run> <workdir>; exit 0 pass, 1 fail, 77 reference missing
set -u
simc=$1; ref=$2; run=$3; work=$4
out=$(ls "$ref/$run.out" "$ref/$run.out.xz" 2>/dev/null | head -1)
[ -f "$ref/$run.inp" ] && [ -n "$out" ] || { echo "reference missing: SKIP"; exit 77; }
rm -rf "$work"; mkdir -p "$work"
"$simc" "$ref/$run.inp" -o "$work" --ngen 5 --fortran-compat > "$work/log" 2>&1 || { tail -5 "$work/log"; exit 1; }
sed -n '/SIDIS: central/,/lab cross section/p' "$work/log" > "$work/cpp.txt"
xzcat -f "$out" | sed -n '/SIDIS: central/,/lab cross section/p' > "$work/fortran.txt"
if [ -s "$work/fortran.txt" ] && diff "$work/cpp.txt" "$work/fortran.txt"; then echo "PASS ($(wc -l < "$work/cpp.txt") lines identical)"; exit 0; fi
echo FAIL; exit 1
