#!/usr/bin/env bash
# validation/scripts/p154_l2.sh — task P15.4: L2 of the parallel event loop. For every case (reference/runs/<case>, its ngen,
# seed 4242, corrected mode, mt64): a --threads 1 run and a --threads N run (default 8), each in its own folder
# simc_runs/<stamp>_P15.4_L2_<case>_t<n>/ (run.log, git.txt, timing.txt, outputs), then simc_compare
#   t<N> vs the Fortran reference A   (compare_vs_A/)
#   t<N> vs t1 (different random streams: statistical)   (compare_vs_t1/)
#   t1 vs the Fortran reference A     (compare_vs_A/ of the t1 folder)
# Usage: p154_l2.sh <simc++> <simc_compare> case[:ngen] ...      (default ngen: the reference's)
# Env: THREADS (8), JOBS_T1 (number of --threads 1 runs made at the same time, default 1).
set -uo pipefail
PROJ=/home/ayerbe/NewDisk/SIMC++_Project
REF=/home/ayerbe/NewDisk/reference/runs
SRC="$PROJ/SIMC++_Claude"
SIMC="$1"; CMP="$2"; shift 2
CFG="$SRC/validation/config/compare_default.json"
N=${THREADS:-8}
export LD_LIBRARY_PATH="${ROOTSYS:-/home/ayerbe/Builds/ROOT6_install}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
STAMP=$(date +%Y%m%d_%H%M%S)

run_one() {   # case ngen threads dir
  local c=$1 n=$2 t=$3 d=$4 inp=""
  for i in "$REF/$c"/*.inp; do b=$(basename "$i" .inp); [[ -f $REF/$c/$b.hist ]] && inp=$i; done
  mkdir -p "$d"
  { echo "commit $(git -C "$SRC" rev-parse HEAD)"; git -C "$SRC" status --short; echo "binary $SIMC"; } > "$d/git.txt"
  local t0=$(date +%s.%N)
  "$SIMC" "$inp" --seed 4242 --ngen "$n" --threads "$t" --quiet -o "$d" > "$d/run.log" 2>&1
  local rc=$? t1=$(date +%s.%N)
  echo "exit $rc wall_s $(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f", b-a}') threads $t ngen $n" > "$d/timing.txt"
}

for e in "$@"; do
  c=${e%%:*}; n=${e#*:}; [[ $n == "$e" ]] && n=""
  h=""; for i in "$REF/$c"/*.hist; do h=$i; done
  [[ -z $n ]] && n=$(grep -m1 'Ngen (request)' "$h" | awk '{print $NF}')
  d1="$PROJ/simc_runs/${STAMP}_P15.4_L2_${c}_t1"; dN="$PROJ/simc_runs/${STAMP}_P15.4_L2_${c}_t$N"
  run_one "$c" "$n" "$N" "$dN"
  run_one "$c" "$n" 1 "$d1"
  stem=$(basename "$(ls "$dN"/*.summary.json | head -1)" .summary.json)
  "$CMP" --a "$REF/$c" --b "$dN/$stem.summary.json" --mode-b corrected --out "$dN/compare_vs_A" --config "$CFG" \
      --label-a "Fortran A" --label-b "SIMC++ t$N" > "$dN/compare_vs_A.log" 2>&1
  "$CMP" --a "$d1/$stem.summary.json" --b "$dN/$stem.summary.json" --mode-a corrected --mode-b corrected \
      --out "$dN/compare_vs_t1" --config "$CFG" --label-a "SIMC++ t1" --label-b "SIMC++ t$N" > "$dN/compare_vs_t1.log" 2>&1
  "$CMP" --a "$REF/$c" --b "$d1/$stem.summary.json" --mode-b corrected --out "$d1/compare_vs_A" --config "$CFG" \
      --label-a "Fortran A" --label-b "SIMC++ t1" > "$d1/compare_vs_A.log" 2>&1
  for x in "$dN/compare_vs_A" "$dN/compare_vs_t1" "$d1/compare_vs_A"; do
    echo "$c $(basename "$(dirname "$x")") $(basename "$x"): $(grep -m1 VERDICT "$x/summary.txt" 2>/dev/null)"
  done
  echo "$c t$N: $(cat "$dN/timing.txt")   t1: $(cat "$d1/timing.txt")"
done
