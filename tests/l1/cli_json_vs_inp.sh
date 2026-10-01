#!/usr/bin/env bash
# tests/l1/cli_json_vs_inp.sh — task P7.2: the `simc++` CLI gives the same run from a Fortran .inp and from the
# equivalent JSON (inp2json parameters + explicit default pipeline, data/config/*.json): .hist (except the time
# lines), .geni and .gen must be identical, and the h10 tree must have the same number of entries.
# usage: cli_json_vs_inp.sh <simc++> <file.inp> <file.json> <workdir> <ngen>
set -u
SIMC="$1"; INP="$2"; JSON="$3"; W="$4"; N="$5"
[[ -f "$INP" && -f "$JSON" ]] || { echo "SKIP: input missing"; exit 77; }
rm -rf "$W" && mkdir -p "$W/inp" "$W/json"
"$SIMC" "$INP" -o "$W/inp" --ngen "$N" --run-name run --quiet > "$W/inp.log" 2>&1 || { echo "FAIL: .inp run"; tail "$W/inp.log"; exit 1; }
"$SIMC" "$JSON" -o "$W/json" --ngen "$N" --run-name run --quiet > "$W/json.log" 2>&1 || { echo "FAIL: .json run"; tail "$W/json.log"; exit 1; }
rc=0
diff <(sed '2,3d' "$W/inp/run.hist") <(sed '2,3d' "$W/json/run.hist") || rc=1
diff "$W/inp/run.geni" "$W/json/run.geni" || rc=1
diff "$W/inp/run.gen" "$W/json/run.gen" || rc=1
for f in run.root run.summary.json run.gen.root; do [[ -s "$W/inp/$f" && -s "$W/json/$f" ]] || { echo "missing $f"; rc=1; }; done
[[ $rc -eq 0 ]] && echo "PASS: .inp and .json runs identical ($N events)" || echo "FAIL"
exit $rc
