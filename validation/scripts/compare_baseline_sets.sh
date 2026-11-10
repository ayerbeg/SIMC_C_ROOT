#!/usr/bin/env bash
# validation/scripts/compare_baseline_sets.sh — task P15.1 (D-054): prove that a regenerated baseline set differs from the
# previous one ONLY by the new summary keys. For every NEW/<prefix>_<suffix> folder the OLD folder of the same suffix is
# compared file by file: .root / .gen.root with simc_root_identity (bitwise), the text files (.hist .geni .gen
# .summary.txt .summary.json) after removing the time stamps; any text difference must consist only of ADDED lines that
# carry the allowed keys (default: fortran_compat, rng). README.txt / MD5SUMS are skipped (they name the commit).
# usage: compare_baseline_sets.sh <simc_root_identity> <baselines dir> <new prefix> <old-prefix-for-suffix map file>
#   map file lines: "<suffix> <old prefix>"
set -uo pipefail
IDENT=$1; DIR=$2; NEW=$3; MAP=$4
ALLOWED='fortran_compat|rng'
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
      *.root) "$IDENT" "$O/$f" "$N/$f" > /dev/null 2>&1 || { echo "  $suffix: $f DIFFERENT"; bad=1; } ;;
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
  [[ $bad == 0 ]] && echo "identical except the new keys: $suffix ($old -> $NEW)" || status=1
done < "$MAP"
echo "$npairs pairs, $([[ $status == 0 ]] && echo PASS || echo FAIL)"
exit $status
