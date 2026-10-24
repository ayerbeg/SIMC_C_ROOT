#!/usr/bin/env bash
# validation/scripts/p12_wtcontr_recompare.sh — task P12.3 (D-053): re-run simc_compare (not the simulations) on the
# existing L2 run folders of P9.2, P10.2, P11.1, P12.2 and the P12 director seed runs, with the wtcontr sigma from the
# pass_cuts subset of each run (the run's mode taken from the folder label: compat / corrected). The old run folders
# are not modified: every comparison goes to a NEW folder simc_runs/<timestamp>_P12.3_wtcontr_recompare/<run folder>/.
# Prints a table: run folder, mode, reference, old wtcontr pull (the folder's compare_vs_A, if any), new pull, old and
# new sigma of A and B, the new verdict.
# Usage: p12_wtcontr_recompare.sh <simc_compare> [glob ...]   (D-040: simc_compare from the same commit)
set -uo pipefail
PROJ=/home/ayerbe/NewDisk/SIMC++_Project
RUNS=$PROJ/simc_runs
REF=/home/ayerbe/NewDisk/reference/runs
CMP="$(realpath "$1")"; shift
CFG="$PROJ/SIMC++_Claude/validation/config/compare_default.json"
GLOBS=("$@")
[[ ${#GLOBS[@]} -eq 0 ]] && GLOBS=("*_P9.2_L2_*" "*_P10.2_L2_*" "*_P11.1_L2_*" "*_P12.2_L2_*" "*_P12dir_*")
OUT="$RUNS/$(date +%Y%m%d_%H%M%S)_P12.3_wtcontr_recompare"
mkdir -p "$OUT"
printf "run\tmode\treference\told_pull\tnew_pull\told_sigA\tnew_sigA\told_sigB\tnew_sigB\tnew_verdict\n" | tee "$OUT/table.tsv"
cd "$RUNS"
for g in "${GLOBS[@]}"; do
  for d in $g; do
    [[ -d "$d" && "$d" != *ABORTED* ]] || continue
    js=( "$d"/*.summary.json ); [[ -f "${js[0]}" ]] || continue
    tag=$(basename "${js[0]}" .summary.json)
    ref="$REF/$tag/$tag.hist"; [[ -f "$ref" ]] || continue
    mode=corrected; [[ "$d" == *compat* ]] && mode=compat
    o="$OUT/$d"
    "$CMP" --a "$ref" --b "${js[0]}" --mode-b $mode --out "$o" --config "$CFG" --no-pdf > "$o.log" 2>&1
    row() { grep -m1 "^wtcontr " "$1" | awk '{print $4"\t"$6"\t"$7}'; }   # errA errB pull
    old="-\t-\t-"; [[ -f "$d/compare_vs_A/summary.txt" ]] && old=$(row "$d/compare_vs_A/summary.txt")
    new=$(row "$o/summary.txt")
    IFS=$'\t' read -r oA oB oP <<< "$(echo -e "$old")"
    IFS=$'\t' read -r nA nB nP <<< "$new"
    v=$(grep -m1 VERDICT "$o/summary.txt" | awk '{print $2}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$d" "$mode" "$tag" "$oP" "$nP" "$oA" "$nA" "$oB" "$nB" "$v" | tee -a "$OUT/table.tsv"
  done
done
echo "outputs: $OUT"
