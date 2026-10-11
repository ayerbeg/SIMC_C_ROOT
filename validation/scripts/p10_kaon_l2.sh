#!/usr/bin/env bash
# validation/scripts/p10_kaon_l2.sh — task P10.2: the kaon L2 runs (seed 4242, 200000 successes) in both modes against
# the P10.F1 Fortran references, through the shared tools/run_simc.sh (serialized, archived in simc_runs/), then
# simc_compare (vs A, and vs B for h_k1_rad_dec0) and, for the corrected runs, simc_mode_diff vs the compat run.
# Usage: p10_kaon_l2.sh <simc_compare> <simc_mode_diff> [tag ...]   (D-040: the compare tools from the same commit)
set -uo pipefail
PROJ=/home/ayerbe/NewDisk/SIMC++_Project
REF=/home/ayerbe/NewDisk/reference/runs
CMP="$1"; MD="$2"; shift 2
CFG="$PROJ/SIMC++_Claude/validation/config/compare_default.json"
TAGS=("$@")
[[ ${#TAGS[@]} -eq 0 ]] && TAGS=(h_k1_rad_dec0 h_k1_rad_dec1 h_k2_rad_dec0 h_k2_rad_dec1 h_k1_sig0_rad_dec0 d_k1_sigm_rad_dec1
                              he3_k1_rad_dec1 he3coh_k2c_rad_dec1 h_k1_norad_dec0)
first=1
for tag in "${TAGS[@]}"; do
  run=kaon_${tag}_A
  for mode in compat corrected; do
    args="--seed 4242 --ngen 200000 --quiet"; [[ $mode == compat ]] && args="$args --fortran-compat"
    nb="--no-build"; [[ $first == 1 ]] && nb=""; first=0
    "$PROJ/tools/run_simc.sh" "$REF/$run/$run.inp" $nb --label "P10.2_L2_${mode}_${tag}" --args "$args" > /dev/null 2>&1
    d="$PROJ/simc_runs/$(readlink "$PROJ/simc_runs/latest")"
    "$CMP" --a "$REF/$run" --b "$d/$run.summary.json" --out "$d/compare_vs_A" --config "$CFG" > "$d/compare_vs_A.log" 2>&1
    if [[ -d "$REF/kaon_${tag}_B" ]]; then
      "$CMP" --a "$REF/kaon_${tag}_B" --b "$d/$run.summary.json" --out "$d/compare_vs_B" --config "$CFG" > "$d/compare_vs_B.log" 2>&1
    fi
    if [[ $mode == compat ]]; then cdir="$d"; else "$MD" "$cdir/$run.root" "$d/$run.root" > "$d/mode_diff_vs_compat.txt" 2>&1; fi
    echo "$tag $mode $d $(grep -m1 VERDICT "$d/compare_vs_A/summary.txt")"
  done
done
