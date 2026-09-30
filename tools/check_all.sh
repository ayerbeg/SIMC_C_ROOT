#!/usr/bin/env bash
# check_all.sh - Clean build of everything in SIMC++_Claude plus the full test suite, logged.
#
# Usage: tools/check_all.sh [--label NAME]
#
# Builds from scratch (fresh build dirs):
#   1. the main CMake project (library, tools, unit + L1 tests)  -> ctest
#   2. the comparison tool in validation/ (own CMake project)     -> ctest
# Writes simc_runs/<timestamp>_check[_label]/ {summary.md, git.txt, build.log, tests.log,
# val_build.log, val_tests.log} and appends one line to simc_runs/index.tsv
# (columns: run_id overall build run tests seconds config; run = n/a, config = check_all@<commit>).
# Exit code: 0 only if both builds and both test suites pass.

set -uo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC_DIR="$PROJECT_DIR/SIMC++_Claude"
RUNS_DIR="$PROJECT_DIR/simc_runs"
ROOT_DIR="${ROOTSYS:-/home/ayerbe/Builds/ROOT6_install}"
export LD_LIBRARY_PATH="$ROOT_DIR/lib:/usr/local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
LABEL=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --label) LABEL="$2"; shift ;;
    -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
  esac; shift
done

mkdir -p "$RUNS_DIR"
exec 9>"$RUNS_DIR/.lock"; flock 9          # same lock as run_simc.sh

COMMIT="$(git -C "$SRC_DIR" rev-parse --short HEAD)"
DIRTY="$(git -C "$SRC_DIR" status --porcelain | wc -l)"
RUN_ID="$(date +%Y%m%d_%H%M%S)_check${LABEL:+_$LABEL}"
OUT="$RUNS_DIR/$RUN_ID"; mkdir -p "$OUT"
BUILD="$OUT/build"; VBUILD="$OUT/val_build"
T0=$(date +%s)
{ echo "commit $COMMIT (branch $(git -C "$SRC_DIR" branch --show-current)), uncommitted files: $DIRTY";
  git -C "$SRC_DIR" status --short; } > "$OUT/git.txt"

jobs=$(nproc)
( cmake -S "$SRC_DIR" -B "$BUILD" -DCMAKE_PREFIX_PATH="$ROOT_DIR" -DCMAKE_BUILD_TYPE=Release &&
  cmake --build "$BUILD" -j"$jobs" ) > "$OUT/build.log" 2>&1; B=$?
T=skipped; if [[ $B -eq 0 ]]; then (cd "$BUILD" && ctest --output-on-failure) > "$OUT/tests.log" 2>&1; T=$?; fi
( cmake -S "$SRC_DIR/validation" -B "$VBUILD" -DCMAKE_PREFIX_PATH="$ROOT_DIR" -DCMAKE_BUILD_TYPE=Release &&
  cmake --build "$VBUILD" -j"$jobs" ) > "$OUT/val_build.log" 2>&1; VB=$?
VT=skipped; if [[ $VB -eq 0 ]]; then (cd "$VBUILD" && ctest --output-on-failure) > "$OUT/val_tests.log" 2>&1; VT=$?; fi
SECS=$(( $(date +%s) - T0 ))

ok() { [[ "$1" == 0 ]]; }
if ok $B && ok $T && ok $VB && ok $VT; then OVERALL=OK; else OVERALL=FAIL; fi
tline() { [[ -f "$1" ]] && grep -E "tests passed" "$1" | tail -1 || echo "not run"; }
WARN=$(grep -c "warning:" "$OUT/build.log" 2>/dev/null)
cat > "$OUT/summary.md" <<EOF
# check_all $RUN_ID : $OVERALL
- commit: $COMMIT (uncommitted files: $DIRTY)
- main build: exit $B, compiler warnings: $WARN
- main tests: $(tline "$OUT/tests.log")
- validation build: exit $VB
- validation tests: $(tline "$OUT/val_tests.log")
- time: ${SECS}s
EOF
# keep the logs, drop the build trees (they are large)
rm -rf "$BUILD" "$VBUILD"
[[ -f "$RUNS_DIR/index.tsv" ]] || printf "run_id\toverall\tbuild\trun\ttests\tseconds\tconfig\n" > "$RUNS_DIR/index.tsv"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$RUN_ID" "$OVERALL" "$B/$VB" "n/a" "$T/$VT" "$SECS" "check_all@$COMMIT" >> "$RUNS_DIR/index.tsv"
ln -sfn "$RUN_ID" "$RUNS_DIR/latest"
cat "$OUT/summary.md"
[[ $OVERALL == OK ]]
