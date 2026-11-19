#!/usr/bin/env bash
# validation/scripts/compare_baseline_sets.sh — task P15.1 (D-054): prove that a regenerated baseline set differs from the
# previous one ONLY by the new summary keys. For every NEW/<prefix>_<suffix> folder the OLD folder of the same suffix is
# compared file by file: .root / .gen.root with simc_root_identity (bitwise), the text files (.hist .geni .gen
# .summary.txt .summary.json) after removing the time stamps; any text difference must consist only of ADDED lines that
# carry the allowed keys (default: fortran_compat, rng). README.txt / MD5SUMS are skipped (they name the commit).
# usage: compare_baseline_sets.sh <simc_root_identity> <baselines dir> <new prefix> <old-prefix-for-suffix map file>
#   map file lines: "<suffix> <old prefix>"
# U.1: ALLOWED_KEYS (env, regex, default "fortran_compat|rng") = the keys of the allowed added text lines;
# ALLOWED_LEAVES (env, regex, default empty) = tree leaves that may differ in a .root (every other leaf must be bitwise
# identical; the differing leaves are listed).
set -uo pipefail
IDENT=$1; DIR=$2; NEW=$3; MAP=$4
ALLOWED=${ALLOWED_KEYS:-fortran_compat|rng}
ALLOWED_LEAVES=${ALLOWED_LEAVES:-}
strip_times() { grep -v -E '^ (BEGIN|END) Time:|^ *"?(begin|end)_time"? *[=:]' "$1"; }
status=0; npairs=0
while read -r suffix old; do
  [[ -z "$suffix" || "$suffix" == \#* ]] && continue
  N="$DIR/${NEW}_$suffix"; O="$DIR/${old}_$suffix"
  if [[ ! -d "$N" || ! -d "$O" ]]; then echo "MISSING pair for $suffix ($N / $O)"; status=1; continue; fi
  npairs=$((npairs+1)); bad=0
  for f in $(cd "$O" && ls | grep -v -E '^(README.txt|MD5SUMS)$'); do
    if [[ ! -f "$N/$f" ]]; then echo "  $suffix: $f missing in the new set"; bad=1; continue; fi
    case "$f" in
      *.root)
         if ! out=$("$IDENT" "$O/$f" "$N/$f" 2>&1); then
           leaves=$(echo "$out" | sed -n -E 's/^ +[^:]+: leaf ([^ ]+) differs in ([0-9]+) entries$/\1(\2)/p' | tr '\n' ' ')
           other=$(echo "$out" | grep -E '^  ' | grep -v -E ': leaf [^ ]+ differs (in|first)|TTree .*DIFFERENT|histograms compared')
           names=$(echo "$out" | sed -n -E 's/^ +[^:]+: leaf ([^ ]+) differs in [0-9]+ entries$/\1/p')
           if [[ -n "$ALLOWED_LEAVES" && -n "$names" && -z "$other" ]] && ! echo "$names" | grep -v -q -x -E "$ALLOWED_LEAVES"; then
             echo "  $suffix: $f differs only in the allowed leaves: $leaves"
           else
             echo "  $suffix: $f DIFFERENT $leaves"; echo "$other" | head -5; bad=1
           fi
         fi ;;
      *) d=$(diff <(strip_times "$O/$f") <(strip_times "$N/$f") | grep -E '^[<>]')
         if [[ -n "$d" ]]; then
           if echo "$d" | grep -q '^<' || echo "$d" | grep -v -q -E "^> *\"?($ALLOWED)\"? *[=:]"; then
             echo "  $suffix: $f differs beyond the new keys:"; echo "$d" | head -5; bad=1
           fi
         fi ;;
    esac
  done
  extra=$(comm -13 <(cd "$O" && ls | sort) <(cd "$N" && ls | sort))
  [[ -n "$extra" ]] && { echo "  $suffix: extra files in the new set: $extra"; bad=1; }
  [[ $bad == 0 ]] && echo "identical except the allowed keys/leaves: $suffix ($old -> $NEW)" || status=1
done < "$MAP"
echo "$npairs pairs, $([[ $status == 0 ]] && echo PASS || echo FAIL)"
exit $status
