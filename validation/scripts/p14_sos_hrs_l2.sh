#!/usr/bin/env bash
# validation/scripts/p14_sos_hrs_l2.sh — tasks P14.1/P14.2: the SOS and HRS L2 runs (seed 4242, ngen = 100000 successes as
# the references) in both modes against the P14.F1 Fortran references reference/runs/p14f1_<tag>_A (and _B where it
# exists), then simc_compare with --mode-b, and simc_mode_diff of each corrected run vs its compat run.
# Each run gets its own folder simc_runs/<timestamp>_P14_L2_<mode>_<tag>/ (run.log, git.txt, timing.txt, outputs,
# compare_vs_A[/B], mode_diff_vs_compat.txt). Runs are made in parallel (JOBS, default 4) with the given build.
# Usage: p14_sos_hrs_l2.sh <build dir (bin/simc++)> <simc_compare> <simc_mode_diff> [tag ...]
#   (D-040: the compare tools must come from the same commit as the build)
#   tags (default): eep_sos pion_sos delta_sos kaon_sos sidis_sos hrs_LeRp_rad hrs_LeRp_norad hrs_ReLp_rad
#   (rho_sos is ~2.9e8 tries per 1e5 events: give it explicitly)
# Env: NGEN (100000), SEED (4242), JOBS (4).
# Note: in the corrected mode the HRS runs differ from the Fortran by construction (D-061 F1: FP path -3.134e+232 in the
# Fortran, physical in the corrected mode; F1e: 45 degree VDC rotation), so the corrected HRS comparison must show them.
set -uo pipefail
PROJ=/home/ayerbe/NewDisk/SIMC++_Project
REF=/home/ayerbe/NewDisk/reference/runs
SRC="$PROJ/SIMC++_Claude"
BUILD="$1"; CMP="$2"; MD="$3"; shift 3
CFG="$SRC/validation/config/compare_default.json"
NGEN=${NGEN:-100000}; SEED=${SEED:-4242}; JOBS=${JOBS:-4}
TAGS=("$@")
[[ ${#TAGS[@]} -eq 0 ]] && TAGS=(eep_sos pion_sos delta_sos kaon_sos sidis_sos hrs_LeRp_rad hrs_LeRp_norad hrs_ReLp_rad)
export LD_LIBRARY_PATH="${ROOTSYS:-/home/ayerbe/Builds/ROOT6_install}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
STAMP=$(date +%Y%m%d_%H%M%S)
LIST="$PROJ/simc_runs/p14_sos_hrs_l2_${STAMP}.list"
: > "$LIST"
for tag in "${TAGS[@]}"; do for mode in compat corrected; do
  echo "$tag $mode $PROJ/simc_runs/${STAMP}_P14_L2_${mode}_${tag}" >> "$LIST"; done; done

run_one() {   # tag mode dir
  local tag=$1 mode=$2 d=$3 run=p14f1_${1}_A
  mkdir -p "$d"
  { echo "commit $(git -C "$SRC" rev-parse HEAD)"; git -C "$SRC" status --short; } > "$d/git.txt"
  local args="--seed $SEED --ngen $NGEN --quiet"; [[ $mode == compat ]] && args="$args --fortran-compat"
  local t0=$(date +%s)
  # (compat HRS runs print the transp.f length warning at every call, as the Fortran: keep only the tail of the log)
  "$BUILD/bin/simc++" "$REF/$run/$run.inp" $args -o "$d" 2>&1 | grep -v "PROBLEM WITH TRANSFORMATION\|Central path length passed\|But the comments say" > "$d/run.log"
  local rc=${PIPESTATUS[0]} t1=$(date +%s)
  echo "exit $rc wall_s $((t1 - t0)) start $(date -d @$t0 +%FT%T) end $(date -d @$t1 +%FT%T)" > "$d/timing.txt"
}
export -f run_one; export SRC REF SEED NGEN BUILD
xargs -P "$JOBS" -L 1 bash -c 'run_one "$0" "$1" "$2"' < "$LIST"

while read -r tag mode d; do
  run=p14f1_${tag}_A
  "$CMP" --a "$REF/$run" --b "$d/$run.summary.json" --mode-b $mode --out "$d/compare_vs_A" --config "$CFG" > "$d/compare_vs_A.log" 2>&1
  if [[ -d "$REF/p14f1_${tag}_B" ]]; then
    "$CMP" --a "$REF/p14f1_${tag}_B" --b "$d/$run.summary.json" --mode-b $mode --out "$d/compare_vs_B" --config "$CFG" > "$d/compare_vs_B.log" 2>&1
  fi
  if [[ $mode == corrected ]]; then
    c="${d/_corrected_/_compat_}"
    "$MD" "$c/$run.root" "$d/$run.root" > "$d/mode_diff_vs_compat.txt" 2>&1
  fi
  v="$(grep -m1 VERDICT "$d/compare_vs_A/summary.txt" 2>/dev/null)"
  echo "$tag $mode $d $(cat "$d/timing.txt") $v"
done < "$LIST"
