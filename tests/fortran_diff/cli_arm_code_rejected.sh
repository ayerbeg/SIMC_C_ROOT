#!/usr/bin/env bash
# tests/fortran_diff/cli_arm_code_rejected.sh — D-032 (task P7.5): `simc++` and `inp2json` refuse the retired arm code
# (SHMS LSA tune) from a Fortran .inp and from a JSON input, with exit status 2 and the message
#   "<arm>_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS"
# in both modes (default and --fortran-compat: the Fortran itself STOPs, dbase.f:895-897/915-917).
#
# usage: cli_arm_code_rejected.sh <simc++> <inp2json> <workdir>
set -uo pipefail
if [[ $# -ne 3 ]]; then sed -n '2,8p' "$0"; exit 2; fi
SIMCPP=$1; INP2JSON=$2; WORK=$3
rm -rf "$WORK"; mkdir -p "$WORK"
MSG_E="electron_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS"
MSG_H="hadron_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS"

cat > "$WORK/arm6_electron.inp" <<'EOF'
begin parm experiment
  ngen = 10
  Ebeam = 10600.
end parm experiment
begin parm simulate
  electron_arm = 6
  hadron_arm = 1
end parm simulate
EOF
sed -e 's/electron_arm = 6/electron_arm = 5/' -e 's/hadron_arm = 1/hadron_arm = 6/' "$WORK/arm6_electron.inp" \
    > "$WORK/arm6_hadron.inp"
echo '{"parameters": {"ngen": 10, "Ebeam": 10600.0, "electron_arm": 6, "hadron_arm": 1}}' > "$WORK/arm6_electron.json"
echo '{"parameters": {"ngen": 10, "Ebeam": 10600.0, "electron_arm": 5, "hadron_arm": 6}}' > "$WORK/arm6_hadron.json"

fail=0
expect_reject() {   # label, expected message, command...
    local label=$1 msg=$2; shift 2
    "$@" > "$WORK/$label.out" 2>&1; local rc=$?
    if [[ $rc -ne 2 ]]; then echo "FAIL $label: exit $rc (expected 2)"; cat "$WORK/$label.out"; fail=1; return; fi
    if ! grep -qF "$msg" "$WORK/$label.out"; then echo "FAIL $label: message not found"; cat "$WORK/$label.out"; fail=1; return; fi
    echo "ok   $label: exit 2, '$(grep -F "$msg" "$WORK/$label.out" | head -1)'"
}
for f in electron hadron; do
    msg=$MSG_E; [[ $f == hadron ]] && msg=$MSG_H
    expect_reject "simcpp_inp_$f" "$msg" "$SIMCPP" "$WORK/arm6_$f.inp" -o "$WORK/out" --quiet
    expect_reject "simcpp_inp_${f}_compat" "$msg" "$SIMCPP" "$WORK/arm6_$f.inp" -o "$WORK/out" --quiet --fortran-compat
    expect_reject "simcpp_json_$f" "$msg" "$SIMCPP" "$WORK/arm6_$f.json" -o "$WORK/out" --quiet
    expect_reject "simcpp_json_${f}_compat" "$msg" "$SIMCPP" "$WORK/arm6_$f.json" -o "$WORK/out" --quiet --fortran-compat
    expect_reject "inp2json_$f" "$msg" "$INP2JSON" -q -o "$WORK/x.json" "$WORK/arm6_$f.inp"
done
[[ $fail -eq 0 ]] && echo "ARM CODE 6 REJECTED: PASS" || echo "ARM CODE 6 REJECTED: FAIL"
exit $fail
