#!/usr/bin/env bash
# validation/scripts/p14_calo_l2.sh — task P14.3: the calorimeter (arms 7/8, pi0 -> 2 gamma) L2 runs (seed 4242, ngen = 100000 successes as
# the references) in both modes against the P14.F2 Fortran references reference/runs/p14f2_<tag>_A (and _B where it
# exists), then simc_compare with --mode-b, and simc_mode_diff of each corrected run vs its compat run.
# Each run gets its own folder simc_runs/<timestamp>_P14.3_L2_<mode>_<tag>/ (run.log, git.txt, timing.txt, outputs,
# compare_vs_A[/B], mode_diff_vs_compat.txt). Runs are made in parallel (JOBS, default 4) with the given build.
# Usage: p14_calo_l2.sh <build dir (bin/simc++)> <simc_compare> <simc_mode_diff> [tag ...]
#   (D-040: the compare tools must come from the same commit as the build)
#   tags (default): excl_h_ng2 (also vs B) excl_h_ng1 excl_d_p excl_d_n delta_h sidis_h eep_caloe
# Env: NGEN (100000), SEED (4242), JOBS (4).
# Note: in the corrected mode the pi0 runs differ from the Fortran by construction in FP%p%path (D-062 G3: one photon
# instead of the sum of both) and in the reaction corrections of the earlier phases (pion D-038, delta D-048, SIDIS D-051).
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/simc_env.sh"   # SIMC_PROJECT_DIR, SIMC_REFERENCE, ROOTSYS -> LD_LIBRARY_PATH
simc_need SIMC_REFERENCE
PROJ=$SIMC_PROJECT_DIR
REF=$SIMC_REFERENCE/runs
SRC="$PROJ/SIMC++_Claude"
BUILD="$1"; CMP="$2"; MD="$3"; shift 3
CFG="$SRC/validation/config/compare_default.json"
NGEN=${NGEN:-100000}; SEED=${SEED:-4242}; JOBS=${JOBS:-4}
TAGS=("$@")
[[ ${#TAGS[@]} -eq 0 ]] && TAGS=(excl_h_ng2 excl_h_ng1 excl_d_p excl_d_n delta_h sidis_h eep_caloe)
STAMP=$(date +%Y%m%d_%H%M%S)
LIST="$PROJ/simc_runs/p14_calo_l2_${STAMP}.list"
: > "$LIST"
for tag in "${TAGS[@]}"; do for mode in compat corrected; do
  echo "$tag $mode $PROJ/simc_runs/${STAMP}_P14.3_L2_${mode}_${tag}" >> "$LIST"; done; done

run_one() {   # tag mode dir
  local tag=$1 mode=$2 d=$3 run=p14f2_${1}_A
  mkdir -p "$d"
  { echo "commit $(git -C "$SRC" rev-parse HEAD)"; git -C "$SRC" status --short; } > "$d/git.txt"
  local args="--seed $SEED --ngen $NGEN --quiet"; [[ $mode == compat ]] && args="$args --fortran-compat"
  local t0=$(date +%s)
  "$BUILD/bin/simc++" "$REF/$run/$run.inp" $args -o "$d" > "$d/run.log" 2>&1
  local rc=$? t1=$(date +%s)
  echo "exit $rc wall_s $((t1 - t0)) start $(date -d @$t0 +%FT%T) end $(date -d @$t1 +%FT%T)" > "$d/timing.txt"
}
export -f run_one; export SRC REF SEED NGEN BUILD
xargs -P "$JOBS" -L 1 bash -c 'run_one "$0" "$1" "$2"' < "$LIST"

while read -r tag mode d; do
  run=p14f2_${tag}_A
  "$CMP" --a "$REF/$run" --b "$d/$run.summary.json" --mode-b $mode --out "$d/compare_vs_A" --config "$CFG" > "$d/compare_vs_A.log" 2>&1
  if [[ -d "$REF/p14f2_${tag}_B" ]]; then
    "$CMP" --a "$REF/p14f2_${tag}_B" --b "$d/$run.summary.json" --mode-b $mode --out "$d/compare_vs_B" --config "$CFG" > "$d/compare_vs_B.log" 2>&1
  fi
  if [[ $mode == corrected ]]; then
    c="${d/_corrected_/_compat_}"
    "$MD" "$c/$run.root" "$d/$run.root" > "$d/mode_diff_vs_compat.txt" 2>&1
  fi
  v="$(grep -m1 VERDICT "$d/compare_vs_A/summary.txt" 2>/dev/null)"
  echo "$tag $mode $d $(cat "$d/timing.txt") $v"
done < "$LIST"
