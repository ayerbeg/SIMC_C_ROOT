#!/usr/bin/env bash
# tests/fortran_diff/cli_shipped_inputs.sh — task P15.2 (D-065, reference/P15.F1_report.md §3): the shipped inputs of
# simc_gfortran/infiles that the Fortran does not run as shipped, and the spect_mode warning (D-065 S2).
#   1. The 18 legacy dotted-name main inputs: the Fortran ignores the 46 '.' names and STOPs with `bad targ.can value`
#      (dbase.f:461). simc++ (both modes) gives the same STOP (simc::FortranStop, exit status 3), and additionally lists
#      the 46 ignored names and points to convert_inputfile.sh.
#   2. test_eep_fe, test_eep_fe_bh, test_hms_shms_validation (legacy names only in the extra file): the Fortran runs with
#      20 ignored names; simc++ runs too (exit 0) and warns loudly with the 20 names (short runs here; the full-run
#      identity is in the regression_shipped tests).
#   3. fpi12_deltapp (H(e,e'pi-)Delta++): --fortran-compat STOPs as the Fortran (dbase.f:823-826, `Can't do Delta+
#      production from neutron for pi0`); the corrected mode runs it (D-038 F1, Category A).
#   4. spect_mode = 2 (not 0/1/-1/-2): the Fortran runs with both arms on, no message (dbase.f:531-532). simc++ runs
#      (exit 0); the corrected mode prints a warning, --fortran-compat does not; and the outputs equal those of
#      spect_mode = 0 (same seed).
#   5./6. rho on A = 4 still refused (D-057 F16); pi0 in the SHMS with doing_decay = 1 (D-065 P1) in both modes.
# usage: cli_shipped_inputs.sh <simc++> <simc_root_identity> <simc_gfortran dir> <reference dir> <workdir>
# Exit 0 pass, 1 fail, 77 inputs missing.
set -uo pipefail
if [[ $# -ne 5 ]]; then sed -n '2,20p' "$0"; exit 2; fi
SIMCPP=$1; IDENT=$2; FDIR=$3; REF=$4; WORK=$5
INF=$FDIR/infiles
[[ -f $INF/test_eep_h.inp && -f $INF/fpi12_deltapp.inp ]] || { echo "SKIP: $INF missing"; exit 77; }
rm -rf "$WORK"; mkdir -p "$WORK"
fail=0
LEGACY=(inputfile_protoype mc_neg_d_z85_zbin mc_neg_h_z85_zbin rhotest test_delta test_delta2 test_delta_h test_delta_h2
        test_eep_d test_eep_h test_kaon_c12 test_kaon_d test_kaon_h test_kaon_he3 test_kaon_he4 test_pion_d test_pion_he
        test_pion_h)
check() {   # label, expected exit, log, fixed strings that must all appear...
    local label=$1 want=$2 rc=$3 log=$4; shift 4
    if [[ $rc -ne $want ]]; then echo "FAIL $label: exit $rc (expected $want)"; tail -4 "$log"; fail=1; return 1; fi
    local s; for s in "$@"; do
        grep -qF -- "$s" "$log" || { echo "FAIL $label: '$s' not in the output"; tail -4 "$log"; fail=1; return 1; }
    done
    return 0
}
# 1. legacy main inputs
n=0
for f in "${LEGACY[@]}"; do
    for mode in corrected compat; do
        extra=(); [[ $mode == compat ]] && extra=(--fortran-compat)
        log=$WORK/legacy_${f}_$mode.log
        "$SIMCPP" "$INF/$f.inp" -o "$WORK/legacy_${f}_$mode" --quiet "${extra[@]}" > "$log" 2>&1; rc=$?
        check "legacy $f ($mode)" 3 $rc "$log" "STOP bad targ.can value" \
            "WARNING: 46 input variable(s) are not SIMC variables and are IGNORED" "targ.can" "spec.e.P" \
            "convert_inputfile.sh" "ignored variable(s) (listed above)" && n=$((n + 1))
    done
done
echo "ok   legacy dotted-name inputs: $n/$(( ${#LEGACY[@]} * 2 )) runs STOP 'bad targ.can value' with the 46 ignored names listed"
# 2. legacy extra files only
for f in test_eep_fe test_eep_fe_bh test_hms_shms_validation; do
    log=$WORK/extra_$f.log
    "$SIMCPP" "$INF/$f.inp" -o "$WORK/extra_$f" --quiet --ngen 20 > "$log" 2>&1; rc=$?
    check "extra-file legacy $f" 0 $rc "$log" "WARNING: 20 input variable(s) are not SIMC variables and are IGNORED" \
        "convert_inputfile.sh" && echo "ok   $f: runs, 20 ignored names listed"
done
# 3. fpi12_deltapp
# U.R: D-038 F1 fixed upstream in 20d46fc: --fortran-compat (target 59ccdf6) runs; target 294f662 keeps the STOP
log=$WORK/deltapp_compat294.log
"$SIMCPP" "$INF/fpi12_deltapp.inp" -o "$WORK/deltapp_compat294" --quiet --fortran-compat --compat-target 294f662 > "$log" 2>&1; rc=$?
check "fpi12_deltapp (compat 294f662)" 3 $rc "$log" "STOP Can't do Delta+ production from neutron for pi0" \
    && echo "ok   fpi12_deltapp --fortran-compat --compat-target 294f662: Fortran STOP (294f662 dbase.f:823-826)"
log=$WORK/deltapp_compat.log
"$SIMCPP" "$INF/fpi12_deltapp.inp" -o "$WORK/deltapp_compat" --quiet --ngen 20 --seed 4242 --fortran-compat > "$log" 2>&1; rc=$?
check "fpi12_deltapp (compat 59ccdf6)" 0 $rc "$log" \
    && echo "ok   fpi12_deltapp --fortran-compat (59ccdf6): runs (20d46fc dbase.f:832)"
log=$WORK/deltapp_corrected.log
"$SIMCPP" "$INF/fpi12_deltapp.inp" -o "$WORK/deltapp_corrected" --quiet --ngen 20 --seed 4242 > "$log" 2>&1; rc=$?
check "fpi12_deltapp (corrected)" 0 $rc "$log" && echo "ok   fpi12_deltapp corrected: runs (D-038 F1)"
# 4. spect_mode = 2
SRC=$REF/runs/p15f1_eep_sm1_A/p15f1_eep_sm1_A.inp
if [[ -f $SRC ]]; then
    for sm in 0 2; do
        sed -e "s/^\( *spect_mode *= *\)-\{0,1\}[0-9]*/\1$sm/" "$SRC" > "$WORK/sm$sm.inp"
        grep -qE "^ *spect_mode *= *$sm\b" "$WORK/sm$sm.inp" || { echo "FAIL: could not set spect_mode $sm"; fail=1; }
        for mode in corrected compat; do
            extra=(); [[ $mode == compat ]] && extra=(--fortran-compat)
            "$SIMCPP" "$WORK/sm$sm.inp" -o "$WORK/sm${sm}_$mode" --quiet --ngen 500 --seed 4242 "${extra[@]}" \
                > "$WORK/sm${sm}_$mode.log" 2>&1; rc=$?
            check "spect_mode $sm ($mode)" 0 $rc "$WORK/sm${sm}_$mode.log" || continue
            if [[ $sm == 2 && $mode == corrected ]]; then
                check "spect_mode 2 warning" 0 0 "$WORK/sm2_corrected.log" "spect_mode = 2 is not one of 0, 1, -1, -2" \
                    && echo "ok   spect_mode 2 (corrected): warning printed, run continues"
            elif grep -q "is not one of 0, 1, -1, -2" "$WORK/sm${sm}_$mode.log"; then
                echo "FAIL spect_mode $sm ($mode): unexpected warning"; fail=1
            fi
        done
    done
    for mode in corrected compat; do
        if "$IDENT" "$WORK/sm0_$mode/sm0.root" "$WORK/sm2_$mode/sm2.root" > "$WORK/sm_ident_$mode.log" 2>&1 &&
           diff <(sed 2,3d "$WORK/sm0_$mode/sm0.hist" | grep -v "spect_mode\|sm0") \
                <(sed 2,3d "$WORK/sm2_$mode/sm2.hist" | grep -v "spect_mode\|sm2") > /dev/null; then
            echo "ok   spect_mode 2 = spect_mode 0 ($mode): tree and .hist identical"
        else echo "FAIL spect_mode 2 != spect_mode 0 ($mode)"; tail -3 "$WORK/sm_ident_$mode.log"; fail=1; fi
    done
else echo "note: $SRC missing, spect_mode part skipped"; fi
# 5. D-057 F16 kept: rho on A = 4 with a valid target can is still refused at start-up in both modes (exit 2); the refusal
#    now comes after the dbase.f:461/465 target STOPs (rhotest above gets the Fortran 'bad targ.can value').
SRC=$REF/runs/rho_h_k1_rad_A/rho_h_k1_rad_A.inp
if [[ -f $SRC ]]; then
    sed -e 's/^\( *targ%A *= *\)[0-9.]*/\14.0/' -e 's/^\( *targ%mass_amu *= *\)[0-9.]*/\14.002602/' "$SRC" > "$WORK/rho_a4.inp"
    for mode in corrected compat; do
        extra=(); [[ $mode == compat ]] && extra=(--fortran-compat)
        "$SIMCPP" "$WORK/rho_a4.inp" -o "$WORK/rho_a4_$mode" --quiet "${extra[@]}" > "$WORK/rho_a4_$mode.log" 2>&1; rc=$?
        check "rho A=4 ($mode)" 2 $rc "$WORK/rho_a4_$mode.log" "exclusive rho production is only set up for H, D and 3He" \
            && echo "ok   rho A=4 ($mode): refused at start-up (D-057 F16)"
    done
else echo "note: $SRC missing, rho part skipped"; fi
# 6. pi0 with the SHMS as hadron arm (D-065 P1, Category B) and doing_decay = 1: the Fortran transports the pi0 itself with
#    decay_flag = doing_decay and STOPs at the first decay test (project.f:87-91, as D-062 G1); --fortran-compat does the
#    same; the corrected mode runs (a pi0 is never decayed in flight, D-037 F9 decaysInFlight; the photons of pizero_decay
#    are not tracked by HMS/SHMS).
SRC=$REF/runs/p15f1_pi0_shms_A/p15f1_pi0_shms_A.inp
if [[ -f $SRC ]]; then
    sed -e 's/^\( *doing_decay = \)0/\11/' "$SRC" > "$WORK/pi0dec.inp"
    "$SIMCPP" "$WORK/pi0dec.inp" -o "$WORK/pi0dec_compat" --quiet --ngen 100 --seed 4242 --fortran-compat \
        > "$WORK/pi0dec_compat.log" 2>&1; rc=$?
    check "pi0 SHMS doing_decay (compat)" 3 $rc "$WORK/pi0dec_compat.log" "error, cannot decay particle with mass=  134.98" \
        && echo "ok   pi0 in the SHMS, doing_decay 1 (compat): Fortran STOP"
    "$SIMCPP" "$WORK/pi0dec.inp" -o "$WORK/pi0dec_corrected" --quiet --ngen 100 --seed 4242 > "$WORK/pi0dec_corrected.log" 2>&1
    check "pi0 SHMS doing_decay (corrected)" 0 $? "$WORK/pi0dec_corrected.log" \
        && echo "ok   pi0 in the SHMS, doing_decay 1 (corrected): runs, pi0 not decayed"
else echo "note: $SRC missing, pi0 part skipped"; fi
[[ $fail -eq 0 ]] && echo "SHIPPED INPUTS / spect_mode warning (P15.2): PASS" || echo "SHIPPED INPUTS / spect_mode warning (P15.2): FAIL"
exit $fail
