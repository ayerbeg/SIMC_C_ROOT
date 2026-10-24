#!/usr/bin/env bash
# validation/scripts/p12_sidis_l2.sh — task P12.2: the SIDIS L2 runs (seed 4242, the reference ngen = 100000 successes)
# in both modes against the P12.F1 Fortran references (reference/runs/sidis_*), through the shared tools/run_simc.sh
# (serialized, archived in simc_runs/), then simc_compare (vs A, and vs B for h_k1_pip_rad) and, for the corrected runs,
# simc_mode_diff vs the compat run. sidis_d_k1_pip_fermi_rad_A has no Fortran .hist (the run STOPs, D-051 F2): only a
# corrected run is made (completes; no comparison).
# Usage: p12_sidis_l2.sh <simc_compare> <simc_mode_diff> [tag ...]   (D-040: the compare tools from the same commit)
set -uo pipefail
PROJ=/home/ayerbe/NewDisk/SIMC++_Project
REF=/home/ayerbe/NewDisk/reference/runs
CMP="$1"; MD="$2"; shift 2
CFG="$PROJ/SIMC++_Claude/validation/config/compare_default.json"
NGEN=${NGEN:-100000}
TAGS=("$@")
[[ ${#TAGS[@]} -eq 0 ]] && TAGS=(h_k1_pip_rad h_k1_pip_norad h_k1_pim_rad d_k1_pip_rad h_k1_kp_rad h_k2_pip_rad c12_k1_pip_rad d_k1_pip_fermi_rad)
first=1
for tag in "${TAGS[@]}"; do
  run=sidis_${tag}_A
  modes="compat corrected"; [[ ! -f "$REF/$run/$run.hist" ]] && modes="corrected"
  for mode in $modes; do
    args="--seed 4242 --ngen $NGEN --quiet"; [[ $mode == compat ]] && args="$args --fortran-compat"
    nb="--no-build"; [[ $first == 1 ]] && nb=""; first=0
    "$PROJ/tools/run_simc.sh" "$REF/$run/$run.inp" $nb --label "P12.2_L2_${mode}_${tag}" --args "$args" > /dev/null 2>&1
    d="$PROJ/simc_runs/$(readlink "$PROJ/simc_runs/latest")"
    if [[ -f "$REF/$run/$run.hist" ]]; then
      "$CMP" --a "$REF/$run" --b "$d/$run.summary.json" --mode-b $mode --out "$d/compare_vs_A" --config "$CFG" > "$d/compare_vs_A.log" 2>&1
      if [[ -d "$REF/sidis_${tag}_B" ]]; then
        "$CMP" --a "$REF/sidis_${tag}_B" --b "$d/$run.summary.json" --mode-b $mode --out "$d/compare_vs_B" --config "$CFG" > "$d/compare_vs_B.log" 2>&1
      fi
    fi
    if [[ $mode == compat ]]; then cdir="$d"; elif [[ -n "${cdir:-}" && $modes != corrected ]]; then
      "$MD" "$cdir/$run.root" "$d/$run.root" > "$d/mode_diff_vs_compat.txt" 2>&1; fi
    v="(no reference)"; [[ -f "$d/compare_vs_A/summary.txt" ]] && v="$(grep -m1 VERDICT "$d/compare_vs_A/summary.txt")"
    echo "$tag $mode $d $v"
  done
  cdir=""
done
