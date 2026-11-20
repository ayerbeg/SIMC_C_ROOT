#!/usr/bin/env bash
# tests/unit/cli_tree_layout.sh — task U.2 (D-082): the tree-layout switch of simc++.
#   default                         -> upstream 59ccdf6 layout, H(e,e'p): 53 h10 leaves
#   --tree-layout legacy            -> 294f662 layout: 46 leaves; the .hist/.geni/.gen equal the default run's
#   JSON "output": {"layout": "legacy"} -> 46 leaves; --tree-layout upstream wins over it (53)
#   --tree-layout bogus -> exit 1; JSON "output": {"layout": "bogus"} -> exit 2
# usage: cli_tree_layout.sh <simc++> <simc_root_identity> <inp2json> <input.inp> <workdir>
set -u
SIMC=$1; IDENT=$2; INP2JSON=$3; INP=$4; W=$5
[[ -f $INP ]] || { echo "input $INP missing: SKIP"; exit 77; }
rm -rf "$W"; mkdir -p "$W"; fail=0
leaves() { "$IDENT" "$1" "$1" | grep -o 'x [0-9]* leaves' | head -1 | tr -dc 0-9; }
run() { local d=$1; shift; "$SIMC" "$@" -o "$W/$d" --run-name r --ngen 200 --seed 4242 --quiet > "$W/$d.log" 2>&1; }
expect() { if [[ "$2" == "$3" ]]; then echo "ok   $1 ($2)"; else echo "FAIL $1: got '$2', expected '$3'"; fail=1; fi; }
run def "$INP"; expect "default layout: leaves" "$(leaves "$W/def/r.root")" 53
run leg "$INP" --tree-layout legacy; expect "--tree-layout legacy: leaves" "$(leaves "$W/leg/r.root")" 46
for x in hist geni gen; do
  expect "legacy .$x = default .$x" "$(diff <(sed 2,3d "$W/def/r.$x") <(sed 2,3d "$W/leg/r.$x") > /dev/null && echo same)" same
done
"$INP2JSON" "$INP" > "$W/p.json" 2> /dev/null || { echo "inp2json failed"; exit 1; }
python3 -I -c 'import json,sys; j=json.load(open(sys.argv[1])); j=j if "parameters" in j else {"parameters": j}; j["output"]={"layout": sys.argv[3]}; json.dump(j, open(sys.argv[2], "w"))' "$W/p.json" "$W/leg.json" legacy
python3 -I -c 'import json,sys; j=json.load(open(sys.argv[1])); j=j if "parameters" in j else {"parameters": j}; j["output"]={"layout": sys.argv[3]}; json.dump(j, open(sys.argv[2], "w"))' "$W/p.json" "$W/bad.json" bogus
run jleg "$W/leg.json"; expect "JSON output.layout legacy: leaves" "$(leaves "$W/jleg/r.root")" 46
run jup "$W/leg.json" --tree-layout upstream; expect "CLI upstream over JSON legacy: leaves" "$(leaves "$W/jup/r.root")" 53
run bad1 "$INP" --tree-layout bogus; expect "--tree-layout bogus: exit" "$?" 1
run bad2 "$W/bad.json"; expect "JSON layout bogus: exit" "$?" 2
exit $fail
