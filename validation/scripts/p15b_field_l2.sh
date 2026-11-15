#!/usr/bin/env bash
# validation/scripts/p15b_field_l2.sh — task P15.2b: L2 of the polarized-target field (using_tgt_field) runs, seed 4242
# (mt64), in both modes, against the Fortran archive of the shipped poltar_exclusivepion_example (RANLUX default seed,
# 50000 events): simc_compare --mode-b, and simc_mode_diff of the corrected run vs the compat run.
# Folders simc_runs/<timestamp>_P15.2b_L2_<mode>_poltar/ (run.log, git.txt, timing.txt, outputs, compare_vs_A,
# mode_diff_vs_compat.txt).
# Usage: p15b_field_l2.sh <build dir (bin/simc++)> <simc_compare> <simc_mode_diff>   (env SEED 4242, NGEN 50000)
# Expected corrected-mode differences: D-067 T1/T2 (not active for the example: targ_Bphi = 270, targ_Bangle > theta),
# T14 (no pi0 here), and the corrections of earlier phases (pion D-038 F5 CENTRAL.sigcc, D-006 pass_cuts).
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/simc_env.sh"   # SIMC_PROJECT_DIR, SIMC_REFERENCE, ROOTSYS -> LD_LIBRARY_PATH
simc_need SIMC_REFERENCE
PROJ=$SIMC_PROJECT_DIR
REF=$SIMC_REFERENCE/runs/shipped_poltar_exclusivepion_example
SRC="$PROJ/SIMC++_Claude"
BUILD="$1"; CMP="$2"; MD="$3"
CFG="$SRC/validation/config/compare_default.json"
SEED=${SEED:-4242}; NGEN=${NGEN:-50000}
STAMP=$(date +%Y%m%d_%H%M%S)
run=poltar_exclusivepion_example
for mode in compat corrected; do
  d="$PROJ/simc_runs/${STAMP}_P15.2b_L2_${mode}_poltar"; mkdir -p "$d"
  { echo "commit $(git -C "$SRC" rev-parse HEAD)"; git -C "$SRC" status --short; } > "$d/git.txt"
  args="--seed $SEED --ngen $NGEN --quiet"; [[ $mode == compat ]] && args="$args --fortran-compat"
  ( t0=$(date +%s); "$BUILD/bin/simc++" "$REF/$run.inp" $args -o "$d" > "$d/run.log" 2>&1; rc=$?
    echo "exit $rc wall_s $(( $(date +%s) - t0 ))" > "$d/timing.txt" ) &
done
wait
for mode in compat corrected; do
  d="$PROJ/simc_runs/${STAMP}_P15.2b_L2_${mode}_poltar"
  "$CMP" --a "$REF" --b "$d/$run.summary.json" --mode-b $mode --out "$d/compare_vs_A" --config "$CFG" > "$d/compare_vs_A.log" 2>&1
  [[ $mode == corrected ]] && "$MD" "${d/_corrected_/_compat_}/$run.root" "$d/$run.root" > "$d/mode_diff_vs_compat.txt" 2>&1
  echo "$mode $d $(cat "$d/timing.txt") $(grep -m1 VERDICT "$d/compare_vs_A/summary.txt" 2>/dev/null)"
done
