#!/usr/bin/env bash
# examples/ideal_box/cli_check.sh SIMCPP INPUT PIPELINE OUTDIR — the ideal_box plug-in chosen purely by the JSON:
# simc++ INPUT --pipeline PIPELINE must run, print the resolved pipeline with ideal_box and write a .geni with the
# plug-in's counters. Exit 77 if INPUT is missing.
set -uo pipefail
SIMCPP=$1; INPUT=$2; PIPE=$3; OUT=$4
[[ -f "$INPUT" ]] || { echo "SKIP: $INPUT missing"; exit 77; }
rm -rf "$OUT"; mkdir -p "$OUT"
"$SIMCPP" "$INPUT" --pipeline "$PIPE" --ngen 500 --seed 7 --run-name cli -o "$OUT" --quiet --print-pipeline \
    > "$OUT/run.log" 2>&1 || { echo "simc++ failed"; tail -n 20 "$OUT/run.log"; exit 1; }
grep -q '"ideal_box"' "$OUT/run.log" || { echo "resolved pipeline does not show ideal_box"; exit 1; }
grep -q 'ideal_box stop counters' "$OUT/cli.geni" || { echo ".geni lacks the ideal_box block"; exit 1; }
echo "PASS: ideal_box selected through $PIPE"; cat "$OUT/cli.geni"
