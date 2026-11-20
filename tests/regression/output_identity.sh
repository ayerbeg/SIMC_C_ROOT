#!/usr/bin/env bash
# tests/regression/output_identity.sh — before/after output identity of a complete simc++ run (task P7.4).
#
# Usage: output_identity.sh SIMCPP ROOT_IDENTITY INPUT SEED RUN_NAME BASELINE_DIR WORK_DIR [SIMCPP_ARGS...]
#
# Runs   SIMCPP INPUT --seed SEED --run-name RUN_NAME -o WORK_DIR --quiet [SIMCPP_ARGS...]
# (P7.5: SIMCPP_ARGS, e.g. --fortran-compat for the exact-Fortran baselines; U.2: --tree-layout legacy is added unless
# SIMCPP_ARGS has --tree-layout, the baselines being in the legacy layout)
# and requires every output file to be identical to the file of the same name in BASELINE_DIR:
#   *.root                     bitwise identical trees (all branches, all entries) and histograms (ROOT_IDENTITY);
#                              the TFile/TKey date stamps are ignored
#   .hist .geni .gen .summary.txt .summary.json
#                              text-identical apart from the time stamps (BEGIN/END Time, begin_time/end_time)
# The set of output files must also be the same. BASELINE_DIR holds the files of a reference commit (README.txt
# and MD5SUMS there are ignored). Exit: 0 identical, 1 different, 77 (CTest SKIP) if the baseline or the input
# is missing, 2 if simc++ fails.
set -uo pipefail
if [[ $# -lt 7 ]]; then sed -n '2,17p' "$0"; exit 2; fi
SIMCPP=$1; ROOTID=$2; INPUT=$3; SEED=$4; RUN=$5; BASE=$6; WORK=$7
shift 7; EXTRA=("$@")
# U.2 (D-082): the stored baselines hold the legacy (294f662) tree layout; unless the caller names a layout, the run
# uses --tree-layout legacy (the upstream default layout is validated against the Fortran itself, tests/l1/eep.cmake).
[[ " ${EXTRA[*]-} " == *" --tree-layout "* ]] || EXTRA+=(--tree-layout legacy)
# U.R: the compat baselines were written with the 294f662 behaviour (D-006 Mk kaon decay); --fortran-compat now targets
# 59ccdf6, so unless the caller names a target the run uses --compat-target 294f662 (no effect without --fortran-compat).
[[ " ${EXTRA[*]-} " == *" --compat-target "* ]] || EXTRA+=(--compat-target 294f662)

if [[ ! -d "$BASE" || ! -f "$INPUT" ]]; then echo "SKIP: baseline $BASE or input $INPUT missing"; exit 77; fi
rm -rf "$WORK"; mkdir -p "$WORK"
if ! "$SIMCPP" "$INPUT" --seed "$SEED" --run-name "$RUN" -o "$WORK" --quiet "${EXTRA[@]}" > "$WORK/run.log" 2>&1; then
    echo "simc++ failed:"; tail -n 20 "$WORK/run.log"; exit 2
fi

strip_times() {   # the only lines that legitimately differ between two runs
    grep -v -E '^ (BEGIN|END) Time:|^ *"?(begin|end)_time"? *[=:]' "$1"
}

status=0
base_files=$(cd "$BASE" && ls | grep -v -E '^(README\.txt|MD5SUMS)$' | sort)
work_files=$(cd "$WORK" && ls | grep -v -E '^run\.log$' | sort)
if [[ "$base_files" != "$work_files" ]]; then
    echo "DIFFERENT file sets:"; diff <(echo "$base_files") <(echo "$work_files"); status=1
fi
for f in $base_files; do
    [[ -f "$WORK/$f" ]] || continue
    case "$f" in
        *.root)
            if "$ROOTID" "$BASE/$f" "$WORK/$f" > "$WORK/$f.identity.txt" 2>&1; then
                echo "identical  $f  ($(grep -c 'identical' "$WORK/$f.identity.txt") trees, $(grep -o '[0-9]* histograms' "$WORK/$f.identity.txt" || echo '0 histograms'))"
            else
                echo "DIFFERENT  $f"; cat "$WORK/$f.identity.txt"; status=1
            fi ;;
        *)
            if diff <(strip_times "$BASE/$f") <(strip_times "$WORK/$f") > "$WORK/$f.diff" 2>&1; then
                echo "identical  $f  ($(wc -l < "$BASE/$f") lines, time stamps ignored)"
            else
                echo "DIFFERENT  $f"; head -n 20 "$WORK/$f.diff"; status=1
            fi ;;
    esac
done
if [[ $status -eq 0 ]]; then
    echo "OUTPUT IDENTITY: PASS ($RUN vs $BASE)"
else
    echo "OUTPUT IDENTITY: FAIL ($RUN vs $BASE)"
    # P10.4 (D-044): the baseline is only bitwise valid for the toolchain recorded in its README.txt
    "$(dirname "${BASH_SOURCE[0]}")/../toolchain_guard.sh" --baseline-hint "$BASE/README.txt" --binary "$SIMCPP"
fi
exit $status
