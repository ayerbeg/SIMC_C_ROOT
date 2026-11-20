#!/bin/bash
# tests/l1/ranlux_fullrun.sh — a complete simc++ run of a Fortran reference run with the RANLUX port (IRandom "ranlux",
# seeded from the .inp random_seed like the Fortran) and --fortran-compat, compared with the Fortran outputs (task P14.1):
#   <run>.hist (except the BEGIN/END Time lines 2-3), <run>.geni, <run>.gen text-identical, and the h10 tree of <run>.root
#   bitwise identical (simc_root_identity). Then the same run in the default (corrected) mode must succeed and its .geni
#   must print the stop counters of its own summary JSON (for the SOS: slit_hor before slit_vert, FORTRAN-BUG D-011).
# U.2: the tree layout follows the reference: upstream (59ccdf6) for a reference whose .hist has the upstream line
# `doing_nuc_elast=` (reference/runs/up_*), legacy (294f662) otherwise; env LAYOUT=upstream|legacy overrides (the 20d46fc
# nuclear-elastic references have the line but the 46-column legacy tree).
# U.R: --fortran-compat reproduces upstream 59ccdf6 by default; a 294f662 reference (no `doing_nuc_elast=` line) is run
# with --compat-target 294f662 (D-006 kaon decay mass Mk, D-038 F1 STOP kept as in the frozen tree). A reference with
# Nntu = 0 (no .root, only .hist/.geni/.gen) is accepted: the tree comparison is replaced by "no .root written".
# usage: ranlux_fullrun.sh <simc++> <simc_root_identity> <ref_run_dir> <run> <workdir>
# Exit: 0 pass, 1 fail, 77 reference missing.
set -u
simc=$1; ident=$2; ref=$3; run=$4; work=$5
[ -f "$ref/$run.inp" ] && { [ -f "$ref/$run.root" ] || [ -f "$ref/$run.hist" ]; } || { echo "reference $ref/$run.* missing: SKIP"; exit 77; }
rm -rf "$work"; mkdir -p "$work/compat" "$work/corrected"
echo '{"pipeline": {"rng": {"type": "ranlux"}}}' > "$work/ranlux.json"
fail=0
layout=legacy; grep -q 'doing_nuc_elast=' "$ref/$run.hist" && layout=upstream
[ -n "${LAYOUT:-}" ] && layout=$LAYOUT
target=59ccdf6; grep -q 'doing_nuc_elast=' "$ref/$run.hist" || target=294f662
echo "  tree layout: $layout, compat target: $target"
"$simc" "$ref/$run.inp" -o "$work/compat" --fortran-compat --compat-target "$target" --pipeline "$work/ranlux.json" --tree-layout "$layout" --quiet > "$work/compat.log" 2>&1 \
    || { echo "simc++ (compat) failed:"; tail -5 "$work/compat.log"; exit 1; }
# U.1: the .hist of SIMC++ has the line `doing_nuc_elast= F` of upstream simc.f:894 (4df7351); a reference of the frozen
# 294f662 Fortran does not have it, so the line is dropped from the C++ .hist only when the reference lacks it.
nucel='^ *doing_nuc_elast= [TF]$'
grep -q 'doing_nuc_elast=' "$ref/$run.hist" && nucel='^$NEVER'
if diff <(sed '2,3d' "$work/compat/$run.hist" | grep -v -E "$nucel") <(sed '2,3d' "$ref/$run.hist") > "$work/hist.diff"; then
    echo "  $run.hist (except BEGIN/END Time): identical"
else echo "  $run.hist: DIFFERENT"; head -20 "$work/hist.diff"; fail=1; fi
for x in geni gen; do
    if cmp -s "$work/compat/$run.$x" "$ref/$run.$x"; then echo "  $run.$x: identical"
    else echo "  $run.$x: DIFFERENT"; diff "$work/compat/$run.$x" "$ref/$run.$x" | head -20; fail=1; fi
done
# U.1: COMMON_LEAVES=1 compares only the h10 leaves present in both trees (a reference with the 59ccdf6 ntuple layout)
identopt=(); [ "${COMMON_LEAVES:-0}" = 1 ] && identopt=(--common-leaves)
# U.1: a 294f662 reference (no doing_nuc_elast line): accept the upstream z/pt2 and pfermi changes (tests/l1/Upstream294.h)
grep -q 'doing_nuc_elast=' "$ref/$run.hist" || identopt+=(--u1-upstream)
if [ ! -f "$ref/$run.root" ]; then   # U.R: Nntu = 0 reference
    if [ -f "$work/compat/$run.root" ]; then echo "  $run.root written but the reference has none"; fail=1
    else echo "  $run.root: none (Nntu = 0), as the reference"; fi
elif "$ident" "${identopt[@]}" "$work/compat/$run.root" "$ref/$run.root" > "$work/identity.log" 2>&1; then
    grep -E "common leaves|U.1 upstream" "$work/identity.log"
    echo "  $run.root: $(tail -1 "$work/identity.log")"
else echo "  $run.root: DIFFERENT"; tail -10 "$work/identity.log"; fail=1; fi
# default (corrected) mode
"$simc" "$ref/$run.inp" -o "$work/corrected" --pipeline "$work/ranlux.json" --tree-layout "$layout" --quiet > "$work/corrected.log" 2>&1 \
    || { echo "simc++ (corrected) failed:"; tail -5 "$work/corrected.log"; exit 1; }
sj="$work/corrected/$run.summary.json"
if grep -q "SOS Trials" "$work/corrected/$run.geni"; then
    want=$(python3 -I -c "import json,sys; s=json.load(open(sys.argv[1]))['stops']['sos']; print(s['slit_hor'], s['slit_vert'], s['slit_oct'])" "$sj")
    got=$(awk '/SOS Trials/{f=1} f && /Slit hor\/vert\/corners/{print $3, $4, $5; exit}' "$work/corrected/$run.geni")
    if [ "$want" = "$got" ]; then echo "  corrected .geni SOS slit hor/vert/corners = summary ($got)"
    else echo "  corrected .geni SOS slit line '$got' != summary '$want'"; fail=1; fi
fi
for blk in "HRSr Trials:hrsr" "HRSl Trials:hrsl"; do   # P14.2: HRS "Slit hor/vert" (FORTRAN-BUG D-011)
    lab=${blk%%:*}; key=${blk##*:}
    if grep -q "$lab" "$work/corrected/$run.geni"; then
        want=$(python3 -I -c "import json,sys; s=json.load(open(sys.argv[1]))['stops'][sys.argv[2]]; print(s['slit_hor'], s['slit_vert'])" "$sj" "$key")
        got=$(awk -v l="$lab" 'index($0,l){f=1} f && /Slit hor\/vert/{print $3, $4; exit}' "$work/corrected/$run.geni")
        if [ "$want" = "$got" ]; then echo "  corrected .geni $key slit hor/vert = summary ($got)"
        else echo "  corrected .geni $key slit line '$got' != summary '$want'"; fail=1; fi
    fi
done
nh=$(grep -c "HRS[rl] Trials" "$work/corrected/$run.geni"); ns=$(grep -c "Successes" "$work/corrected/$run.geni")
nt=$(grep -c "Trials" "$work/corrected/$run.geni")
if [ "$nh" -gt 0 ]; then   # D-061 F8: every block of the corrected .geni has its Successes line
    if [ "$ns" = "$nt" ]; then echo "  corrected .geni: a Successes line in each of the $nt blocks"
    else echo "  corrected .geni: $ns Successes lines for $nt blocks"; fail=1; fi
fi
[ $fail = 0 ] && echo "PASS" || echo "FAIL"
exit $fail
