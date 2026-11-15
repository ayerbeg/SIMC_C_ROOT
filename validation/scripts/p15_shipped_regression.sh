#!/usr/bin/env bash
# validation/scripts/p15_shipped_regression.sh — task P15.2: full-run identity of simc++ against Fortran reference runs.
# For every reference run folder (default: all reference/runs/p15f1_* and shipped_*) simc++ runs the archived .inp with the
# RANLUX port and --fortran-compat (the seed of the .inp, as the Fortran) and the outputs are compared with the Fortran:
#   .hist without the BEGIN/END Time lines (diff), .geni and .gen (cmp), the h10 tree (simc_root_identity, bitwise).
# The run name of a folder is the basename of its .inp that has a .hist next to it (shipped_<x> folders keep the
# shipped file name). Folders whose run the C++ refuses (NotPorted) are reported as REFUSED; since P15.2b (using_tgt_field
# ported) nothing is refused, so REFUSED counts as a failure unless ALLOW_REFUSED=1.
# usage: p15_shipped_regression.sh <simc++> <simc_root_identity> <outdir> [run folder name ...]
# env: JOBS (default 4), REF (default $SIMC_REFERENCE/runs; SIMC_REFERENCE from the environment or check_all.env, see simc_env.sh)
# writes <outdir>/<folder>/ (outputs, run.log, identity.txt) and <outdir>/summary.tsv; exit 0 if every run is IDENTICAL or
# (with ALLOW_REFUSED=1) REFUSED, 1 otherwise (77 when a single requested folder is absent: CTest skip).
# The CTest label regression_shipped (tests/regression/CMakeLists.txt) runs this script once per folder.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/simc_env.sh"   # SIMC_PROJECT_DIR, SIMC_REFERENCE, ROOTSYS -> LD_LIBRARY_PATH
SIMC=$1; IDENT=$2; OUT=$3; shift 3
[[ -n ${REF-} ]] || simc_need SIMC_REFERENCE
REF=${REF:-$SIMC_REFERENCE/runs}; JOBS=${JOBS:-4}
LIST=("$@")
# CTest regression_shipped_<folder> (one folder per test): skip (77) when that reference folder is absent
if [[ ${#LIST[@]} -eq 1 && ! -d "$REF/${LIST[0]}" ]]; then echo "reference $REF/${LIST[0]} missing: SKIP"; exit 77; fi
mkdir -p "$OUT"
if [[ ${#LIST[@]} -eq 0 ]]; then
  for d in "$REF"/p15f1_* "$REF"/shipped_*; do [[ -d $d ]] && LIST+=("$(basename "$d")"); done
fi
echo '{"pipeline": {"rng": {"type": "ranlux"}}}' > "$OUT/ranlux.json"
one() {   # folder
  local f=$1 d="$REF/$1" w="$OUT/$1" run="" inp t0 rc v
  for inp in "$d"/*.inp; do local b; b=$(basename "$inp" .inp); [[ -f "$d/$b.hist" ]] && run=$b; done
  rm -rf "$w"; mkdir -p "$w"
  if [[ -z "$run" ]]; then printf '%s\t-\tNO_REFERENCE\t0\n' "$f" > "$w/row.tsv"; return; fi
  t0=$(date +%s)
  "$SIMC" "$d/$run.inp" -o "$w" --fortran-compat --pipeline "$OUT/ranlux.json" --quiet > "$w/run.log" 2>&1
  rc=$?
  { echo "== $f ($run): simc++ RANLUX --fortran-compat, exit $rc, $(( $(date +%s) - t0 )) s"
    if [[ $rc -ne 0 ]]; then tail -3 "$w/run.log"; fi
    if [[ $rc -eq 0 ]]; then
      diff <(sed 2,3d "$w/$run.hist") <(sed 2,3d "$d/$run.hist") > "$w/hist.diff" && echo "hist: identical" || echo "hist: DIFFERENT"
      cmp "$w/$run.geni" "$d/$run.geni" > /dev/null && echo "geni: identical" || echo "geni: DIFFERENT"
      cmp "$w/$run.gen" "$d/$run.gen" > /dev/null && echo "gen: identical" || echo "gen: DIFFERENT"
      "$IDENT" "$w/$run.root" "$d/$run.root" > "$w/tree.txt" 2>&1 && echo "tree: identical ($(grep -o '[0-9]* entries x [0-9]* leaves' "$w/tree.txt" | head -1))" || { echo "tree: DIFFERENT"; head -5 "$w/tree.txt"; }
    fi; } > "$w/identity.txt"
  if [[ $rc -ne 0 ]] && grep -q "not ported\|NotPorted" "$w/run.log"; then v=REFUSED
  elif [[ $rc -ne 0 ]]; then v=ERROR
  elif grep -q DIFFERENT "$w/identity.txt"; then v=DIFFERENT
  else v=IDENTICAL; fi
  printf '%s\t%s\t%s\t%s\n' "$f" "$run" "$v" "$(( $(date +%s) - t0 ))" > "$w/row.tsv"
  rm -f "$w"/*.root.tmp
}
export -f one; export SIMC IDENT OUT REF
printf '%s\n' "${LIST[@]}" | xargs -P "$JOBS" -I{} bash -c 'one {}'
{ printf 'folder\trun\tverdict\tseconds\n'; for f in "${LIST[@]}"; do cat "$OUT/$f/row.tsv"; done; } > "$OUT/summary.tsv"
column -t -s$'\t' "$OUT/summary.tsv"
[[ ${#LIST[@]} -eq 1 ]] && { cat "$OUT/${LIST[0]}/identity.txt"; grep -h "not ported" "$OUT/${LIST[0]}/run.log"; }
bad='DIFFERENT|ERROR|NO_REFERENCE'; [[ ${ALLOW_REFUSED:-0} == 1 ]] || bad="$bad|REFUSED"
! grep -q -P "\t($bad)\t" "$OUT/summary.tsv"
