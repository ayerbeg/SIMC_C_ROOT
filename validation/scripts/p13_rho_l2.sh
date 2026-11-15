#!/usr/bin/env bash
# validation/scripts/p13_rho_l2.sh — task P13.2: the exclusive-rho L2 runs (seed 4242, ngen = 100000 successes as the
# references) in both modes against the P13.F1 Fortran references (reference/runs/rho_*), then simc_compare (vs A, and vs
# B for h_k1_rad) with --mode-b, and simc_mode_diff of each corrected run vs its compat run.
# Rho runs are slow (~2.9e8 tries per 1e5 events, about 2 h each for SIMC++ and 3 h for the Fortran), so the simc++ runs
# are made in parallel (JOBS, default 4) directly with the given build (not through tools/run_simc.sh, which serializes);
# each run gets its own folder simc_runs/<timestamp>_P13.2_L2_<mode>_<tag>/ (run.log, git.txt, timing.txt, outputs,
# compare_vs_A[/B], mode_diff_vs_compat.txt).
# Usage: p13_rho_l2.sh <build dir (bin/simc++)> <simc_compare> <simc_mode_diff> [tag ...]
#   (D-040: the compare tools must come from the same commit as the build)
# Env: NGEN (100000), SEED (4242), JOBS (4).
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
[[ ${#TAGS[@]} -eq 0 ]] && TAGS=(h_k1_rad h_k1_norad d_k1_rad he3_k1_rad h_k1_dec1_rad h_k1_pim_rad h_k2_rad)
STAMP=$(date +%Y%m%d_%H%M%S)
LIST="$PROJ/simc_runs/p13_rho_l2_${STAMP}.list"
: > "$LIST"
for tag in "${TAGS[@]}"; do for mode in compat corrected; do
  echo "$tag $mode $PROJ/simc_runs/${STAMP}_P13.2_L2_${mode}_${tag}" >> "$LIST"; done; done

run_one() {   # tag mode dir
  local tag=$1 mode=$2 d=$3 run=rho_${1}_A
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

# comparisons (serial, after all runs)
while read -r tag mode d; do
  run=rho_${tag}_A
  "$CMP" --a "$REF/$run" --b "$d/$run.summary.json" --mode-b $mode --out "$d/compare_vs_A" --config "$CFG" > "$d/compare_vs_A.log" 2>&1
  if [[ -d "$REF/rho_${tag}_B" ]]; then
    "$CMP" --a "$REF/rho_${tag}_B" --b "$d/$run.summary.json" --mode-b $mode --out "$d/compare_vs_B" --config "$CFG" > "$d/compare_vs_B.log" 2>&1
  fi
  if [[ $mode == corrected ]]; then
    c="${d/_corrected_/_compat_}"
    "$MD" "$c/$run.root" "$d/$run.root" > "$d/mode_diff_vs_compat.txt" 2>&1
  fi
  v="$(grep -m1 VERDICT "$d/compare_vs_A/summary.txt" 2>/dev/null)"
  echo "$tag $mode $d $(cat "$d/timing.txt") $v"
done < "$LIST"
