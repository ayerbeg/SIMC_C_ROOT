#!/usr/bin/env bash
# tests/regression/threads_test.sh — the parallel event loop (task P15.4, docs/PERFORMANCE.md).
#   repro  <simc++> <simc_root_identity> <inp> <workdir> <N> <ngen>
#          two runs with the same (seed 4242, --threads N) must give identical outputs (.hist/.geni/.gen/summary without
#          the time stamps, the .root trees and histograms through simc_root_identity); the log must report N threads and
#          ngen successes; the summary JSON must carry "threads": N; no worker temporary file may be left behind
#   one    <simc++> <simc_root_identity> <inp> <workdir> <ngen>
#          --threads 1 must be bitwise the run without the option (the sequential loop)
#   refuse <simc++> <inp> <workdir>
#          --threads 2 with --fortran-compat, and with the pipeline rng "ranlux", must stop with exit 2 and a message
set -uo pipefail
mode=$1; shift
# the input (.inp) is argument 3 (repro, one) or 2 (refuse): skip (77) when it is absent (no SIMC_FORTRAN_DIR /
# SIMC_REFERENCE_DIR configured)
inp=$([[ $mode == refuse ]] && echo "${2:-}" || echo "${3:-}")
[[ -f $inp ]] || { echo "SKIP: input $inp not found"; exit 77; }

flt() { grep -v -E '^ (BEGIN|END) Time:|^ *"?(begin|end)_time"? *[=:]' "$1" | sed "s|$2|<DIR>|g"; }

same_dir() {   # same_dir A B IDENT: every file of A equal in B
  local A=$1 B=$2 IDENT=$3 bad=0 n=0 f
  while IFS= read -r f; do
    n=$((n+1))
    [[ -f $B/$f ]] || { echo "missing in $B: $f"; bad=1; continue; }
    case $f in
      *.root) "$IDENT" "$A/$f" "$B/$f" > /dev/null 2>&1 || { echo "DIFF $f"; bad=1; } ;;
      *) diff -q <(flt "$A/$f" "$A") <(flt "$B/$f" "$B") > /dev/null || { echo "DIFF $f"; bad=1; } ;;
    esac
  done < <(cd "$A" && find . -type f | sort)
  [[ $n -gt 0 ]] || { echo "no files in $A"; bad=1; }
  echo "$n files compared: $([[ $bad == 0 ]] && echo IDENTICAL || echo DIFFERENT)"
  return $bad
}

if [[ $mode == repro ]]; then
  SIMC=$1; IDENT=$2; INP=$3; W=$4; N=$5; NGEN=$6
  rm -rf "$W"; mkdir -p "$W"
  for r in a b; do
    "$SIMC" "$INP" -o "$W/$r" --seed 4242 --ngen "$NGEN" --threads "$N" --quiet > "$W/$r.log" 2>&1 ||
      { echo "run $r failed (exit $?)"; tail -5 "$W/$r.log"; exit 1; }
  done
  grep -q "event loop on $N threads" "$W/a.log" || { echo "the log does not report $N threads"; exit 1; }
  succ=$(grep -E -- '---> Last Event' "$W/a.log" | awk '{print $6}')
  [[ $succ == "$NGEN" ]] || { echo "successes $succ != ngen $NGEN"; exit 1; }
  grep -q "\"threads\": $N" "$W"/a/*.summary.json || { echo "summary.json without \"threads\": $N"; exit 1; }
  if ls "$W"/a/*.tmp > /dev/null 2>&1; then echo "worker temporary files left: $(ls "$W"/a/*.tmp)"; exit 1; fi
  same_dir "$W/a" "$W/b" "$IDENT"
elif [[ $mode == one ]]; then
  SIMC=$1; IDENT=$2; INP=$3; W=$4; NGEN=$5
  rm -rf "$W"; mkdir -p "$W"
  "$SIMC" "$INP" -o "$W/default" --seed 4242 --ngen "$NGEN" --quiet > "$W/default.log" 2>&1 || { echo "default run failed"; exit 1; }
  "$SIMC" "$INP" -o "$W/one" --seed 4242 --ngen "$NGEN" --threads 1 --quiet > "$W/one.log" 2>&1 || { echo "--threads 1 failed"; exit 1; }
  if grep -q '"threads"' "$W"/one/*.summary.json; then echo "--threads 1 wrote a threads key"; exit 1; fi
  same_dir "$W/default" "$W/one" "$IDENT"
elif [[ $mode == refuse ]]; then
  SIMC=$1; INP=$2; W=$3
  rm -rf "$W"; mkdir -p "$W"
  "$SIMC" "$INP" -o "$W/c" --ngen 100 --threads 2 --fortran-compat --quiet > "$W/compat.log" 2>&1; rc=$?
  [[ $rc == 2 ]] || { echo "--fortran-compat --threads 2: exit $rc (expected 2)"; exit 1; }
  grep -q "cannot be used with --fortran-compat" "$W/compat.log" || { echo "no compat message"; cat "$W/compat.log"; exit 1; }
  [[ -e $W/c ]] && { echo "--fortran-compat --threads 2 created outputs"; exit 1; }
  echo '{"rng": "ranlux"}' > "$W/ranlux.json"
  "$SIMC" "$INP" -o "$W/r" --ngen 100 --threads 2 --pipeline "$W/ranlux.json" --quiet > "$W/ranlux.log" 2>&1; rc=$?
  [[ $rc == 2 ]] || { echo "rng ranlux --threads 2: exit $rc (expected 2)"; exit 1; }
  grep -q 'needs the rng "mt64"' "$W/ranlux.log" || { echo "no ranlux message"; cat "$W/ranlux.log"; exit 1; }
  echo '{"rng": "ranlux", "threads": 4}' > "$W/ranlux4.json"
  "$SIMC" "$INP" -o "$W/r4" --ngen 100 --pipeline "$W/ranlux4.json" --quiet > "$W/ranlux4.log" 2>&1; rc=$?
  [[ $rc == 2 ]] || { echo "pipeline threads 4 + ranlux: exit $rc (expected 2)"; exit 1; }
  "$SIMC" "$INP" -o "$W/z" --threads 0 > "$W/zero.log" 2>&1; rc=$?
  [[ $rc == 1 ]] || { echo "--threads 0: exit $rc (expected 1)"; exit 1; }
  echo "refusals OK"
else
  echo "usage: $0 repro|one|refuse ..."; exit 2
fi
