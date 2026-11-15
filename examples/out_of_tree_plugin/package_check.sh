#!/usr/bin/env bash
# examples/out_of_tree_plugin/package_check.sh — CTest package_out_of_tree_plugin (task P15.3):
#   1. cmake --install <simc build dir> --prefix <work>/prefix   (the SIMC CMake package)
#   2. configure + build this example against it with find_package(SIMC) (only the installed headers/targets)
#   3. run my_simc on <input> with the toy_window plug-in on the hadron arm (selected by the JSON pipeline) and check
#      that the plug-in was used (resolved pipeline, its .geni counter block) and that the run produced its outputs.
#   4. P15.4 (D-069): without --data-dir both the installed simc++ and my_simc take the INSTALLED data
#      (<prefix>/share/simc/data: --print-data-dir), SIMC_DATA_DIR overrides it (a run with a wrong one fails), and the
#      installed simc++ runs the input without --data-dir.
# usage: package_check.sh <simc build dir> <example source dir> <work dir> <input.inp> <ROOT_DIR (ROOTConfig.cmake dir)>
# Exit 0 pass, 1 fail, 77 input missing.
set -uo pipefail
BUILD=$1; SRC=$2; WORK=$3; INPUT=$4; ROOTP=$5
[[ -f $INPUT ]] || { echo "SKIP: $INPUT missing"; exit 77; }
rm -rf "$WORK"; mkdir -p "$WORK"
cmake --install "$BUILD" --prefix "$WORK/prefix" > "$WORK/install.log" 2>&1 || { echo "install failed"; tail "$WORK/install.log"; exit 1; }
for f in lib/cmake/SIMC/SIMCConfig.cmake lib/cmake/SIMC/SIMCConfigVersion.cmake include/simc/driver/CommandLine.h bin/simc++ \
         share/simc/data/field/trg_field_map.dat; do
    [[ -e $WORK/prefix/$f ]] || { echo "not installed: $f"; exit 1; }
done
echo "installed: $(find "$WORK/prefix" -type f | wc -l) files"
cmake -S "$SRC" -B "$WORK/build" -DCMAKE_PREFIX_PATH="$WORK/prefix" -DROOT_DIR="$ROOTP" -DCMAKE_BUILD_TYPE=Release \
    > "$WORK/configure.log" 2>&1 || { echo "configure failed"; tail -20 "$WORK/configure.log"; exit 1; }
cmake --build "$WORK/build" -j4 > "$WORK/build.log" 2>&1 || { echo "build failed"; tail -20 "$WORK/build.log"; exit 1; }
echo "built my_simc against the installed package ($(grep -c 'warning:' "$WORK/build.log") compiler warnings)"
"$WORK/build/my_simc" "$INPUT" --pipeline "$SRC/toy_window.json" --ngen 300 --seed 4242 --run-name toy -o "$WORK/out" \
    --quiet --print-pipeline > "$WORK/run.log" 2>&1 || { echo "my_simc failed"; tail -20 "$WORK/run.log"; exit 1; }
grep -q '"toy_window"' "$WORK/run.log" || { echo "resolved pipeline does not show toy_window"; exit 1; }
grep -q 'toy_window' "$WORK/out/toy.geni" || { echo ".geni lacks the toy_window counters"; cat "$WORK/out/toy.geni"; exit 1; }
for x in hist geni root summary.json; do [[ -s $WORK/out/toy.$x ]] || { echo "missing output toy.$x"; exit 1; }; done
# 4. the data directory of installed programs (P15.4, D-069)
unset SIMC_DATA_DIR
want=$(cd "$WORK/prefix/share/simc/data" && pwd)
for exe in "$WORK/prefix/bin/simc++" "$WORK/build/my_simc"; do
    got=$("$exe" --print-data-dir)
    [[ $(cd "$got" 2>/dev/null && pwd) == "$want" ]] || { echo "$exe: data directory '$got', expected the installed '$want'"; exit 1; }
    [[ $(SIMC_DATA_DIR=/tmp/simc_elsewhere "$exe" --print-data-dir) == /tmp/simc_elsewhere ]] ||
        { echo "$exe: SIMC_DATA_DIR does not override the data directory"; exit 1; }
done
SIMC_DATA_DIR="$WORK/no_such_data" "$WORK/build/my_simc" "$INPUT" --ngen 10 --seed 4242 --run-name bad -o "$WORK/out_bad" \
    --quiet > "$WORK/run_bad.log" 2>&1 && { echo "my_simc ran with SIMC_DATA_DIR pointing nowhere: data not from the data root"; exit 1; }
"$WORK/prefix/bin/simc++" "$INPUT" --ngen 300 --seed 4242 --run-name inst -o "$WORK/out_inst" --quiet > "$WORK/run_inst.log" 2>&1 ||
    { echo "installed simc++ failed without --data-dir"; tail -20 "$WORK/run_inst.log"; exit 1; }
[[ -s $WORK/out_inst/inst.hist ]] || { echo "installed simc++: no inst.hist"; exit 1; }
echo "data directory of the installed programs: $want (SIMC_DATA_DIR overrides it)"
echo "PASS: out-of-tree plug-in toy_window built with find_package(SIMC) and used by my_simc"
grep -A4 toy_window "$WORK/out/toy.geni"
