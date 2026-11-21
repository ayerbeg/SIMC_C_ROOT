#!/usr/bin/env bash
# validation/scripts/u3_l2_corrected.sh — task U.3: L2 of the default (corrected) mode against the upstream SIMC 59ccdf6
# references reference/runs/u59_* (the same reaction/spectrometer sets as the earlier corrected L2 entries: P9-P15
# scripts p1x_*_l2.sh, p15_spect_mode_l2.sh, p15b_field_l2.sh, which are hard-wired to the 294f662 folders
# <reference>/runs/<run>/<run>.inp and run both modes; this driver takes u59 folder names, runs only the corrected mode
# — the compat mode is bitwise identical to these references, u3_upstream_identity.sh — and uses the same comparison
# tool, simc_compare --mode-b corrected, with validation/config/compare_default.json).
# For every folder: simc++ <inp> --seed SEED --ngen <the reference "Ngen (request)", or NGEN> --quiet (default mt64 RNG,
# upstream tree layout), JOBS runs at a time, each in <outdir>/<folder>/ (run.log, timing.txt, outputs, compare_vs_A/),
# then simc_compare vs the reference folder. A _nontu reference (Nntu = 0: .hist only) is compared at the run level only.
# Writes <outdir>/summary.tsv: folder, verdict, N tests, F failed, gross, wall s, gross failures (name: A -> B, p).
# usage: u3_l2_corrected.sh <simc++> <simc_compare> <outdir> u59_folder ...
# env: SEED (4242), JOBS (4), NGEN (default: the reference's), REF (default $SIMC_REFERENCE/runs)
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/simc_env.sh"   # SIMC_REFERENCE, ROOTSYS -> LD_LIBRARY_PATH
SIMC=$1; CMP=$2; OUT=$3; shift 3
[[ -n ${REF-} ]] || simc_need SIMC_REFERENCE
REF=${REF:-$SIMC_REFERENCE/runs}; JOBS=${JOBS:-4}; SEED=${SEED:-4242}
CFG="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/compare_default.json"
mkdir -p "$OUT"
LIST=("$@")
run_one() {   # folder
  local f=$1 d="$REF/$1" w="$OUT/$1" run="" b n
  for i in "$d"/*.inp; do b=$(basename "$i" .inp); [[ -f "$d/$b.hist" ]] && run=$b; done
  [[ -z $run ]] && { echo "no .hist in $d" > "$w.missing"; return; }
  n=${NGEN:-$(grep -m1 'Ngen (request)' "$d/$run.hist" | awk '{print $NF}')}
  mkdir -p "$w"
  local t0; t0=$(date +%s)
  # (corrected HRS runs: drop the per-call transp.f warning lines from the log, as p14_sos_hrs_l2.sh)
  "$SIMC" "$d/$run.inp" --seed "$SEED" --ngen "$n" --quiet -o "$w" 2>&1 | grep -v "PROBLEM WITH TRANSFORMATION\|Central path length passed\|But the comments say" > "$w/run.log"
  local rc=${PIPESTATUS[0]}
  echo "exit $rc wall_s $(( $(date +%s) - t0 )) seed $SEED ngen $n" > "$w/timing.txt"
  "$CMP" --a "$d" --b "$w/$run.summary.json" --mode-b corrected --out "$w/compare_vs_A" --config "$CFG" \
      --label-a "Fortran 59ccdf6" --label-b "SIMC++ corrected" > "$w/compare_vs_A.log" 2>&1
}
export -f run_one; export SIMC CMP OUT REF SEED CFG NGEN
printf '%s\n' "${LIST[@]}" | xargs -P "$JOBS" -I{} bash -c 'run_one {}'
{ printf 'folder\tverdict\tN\tF\tgross\twall_s\tgross failures\n'
  for f in "${LIST[@]}"; do
    s="$OUT/$f/compare_vs_A/summary.json"
    if [[ ! -f $s ]]; then printf '%s\tNO_COMPARE\t-\t-\t-\t-\t%s\n' "$f" "$(tail -1 "$OUT/$f/compare_vs_A.log" 2>/dev/null)"; continue; fi
    python3 -I - "$s" "$f" "$(awk '{print $4}' "$OUT/$f/timing.txt")" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1])); c = d["counts"]; g = d["config"]["gross_p"]
gross = []
for t in d["tests"]:
    if not t.get("counted", True):
        continue
    if (t.get("method") == "det" and not t.get("pass", True)) or (t.get("p") is not None and t["p"] < g):
        a, b = t.get("a"), t.get("b")
        val = f"{a:.6g} -> {b:.6g}" if isinstance(a, (int, float)) and isinstance(b, (int, float)) else ""
        gross.append(f"{t['name']} ({val}, p={t['p']:.2g})" if t.get("p") is not None else f"{t['name']} ({val})")
print("\t".join([sys.argv[2], d["verdict"], str(c["n_tests"]), str(c["n_fail"]), str(c["n_gross"]), sys.argv[3], "; ".join(gross) or "-"]))
EOF
  done; } > "$OUT/summary.tsv"
column -t -s$'\t' "$OUT/summary.tsv"
