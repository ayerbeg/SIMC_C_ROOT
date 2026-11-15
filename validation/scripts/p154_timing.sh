#!/usr/bin/env bash
# validation/scripts/p154_timing.sh — task P15.4: wall time per try and per event of simc++ (corrected mode, mt64, seed 4242)
# on a set of reference inputs, for a list of thread counts, next to the Fortran reference run of the same input (its
# .hist BEGIN/END time and Ntried: the Fortran wall time per try; the reference ngen is larger).
# Usage: p154_timing.sh <simc++> <outdir> "<threads list>" case:ngen ...
#   e.g. p154_timing.sh build/bin/simc++ /tmp/t "1 2 4 8" eep_h_elastic_rad_A:50000 rho_h_k1_rad_A:1000
#   (thread count "seq": run without the --threads option, for a build before P15.4)
# Output: <outdir>/timing.tsv (case ngen threads wall_s ntried us_per_try ms_per_event speedup fortran_us_per_try) and
# the table on stdout. Run it on an otherwise idle machine.
set -uo pipefail
SIMC=$1; OUT=$2; THREADS=($3); shift 3
REF=/home/ayerbe/NewDisk/reference/runs
mkdir -p "$OUT"
TSV="$OUT/timing.tsv"
printf "case\tngen\tthreads\twall_s\tntried\tus_per_try\tms_per_event\tspeedup\tfortran_us_per_try\n" > "$TSV"
for e in "$@"; do
  f=${e%%:*}; n=${e##*:}; d=$REF/$f
  inp=""; for i in "$d"/*.inp; do b=$(basename "$i" .inp); [[ -f $d/$b.hist ]] && inp=$i; done
  [[ -n $inp ]] || { echo "no input with a .hist in $d"; continue; }
  hist=$d/$(basename "$inp" .inp).hist
  fb=$(grep -m1 'BEGIN Time:' "$hist" | sed 's/.*Time: *//'); fe=$(grep -m1 'END Time:' "$hist" | sed 's/.*Time: *//')
  fn=$(grep -m1 'Ntried' "$hist" | awk '{print $NF}')
  fus=$(awk -v s="$(( $(date -d "$fe" +%s) - $(date -d "$fb" +%s) ))" -v n="$fn" 'BEGIN{printf "%.2f", 1e6*s/n}')
  t1=""
  for t in "${THREADS[@]}"; do
    o="$OUT/${f}_t$t"; rm -rf "$o"
    s0=$(date +%s.%N)
    targs=(--threads "$t"); [[ $t == seq ]] && targs=()   # "seq": no --threads option (a build before P15.4)
    "$SIMC" "$inp" -o "$o" --seed 4242 --ngen "$n" "${targs[@]}" --quiet > "$o.log" 2>&1 || echo "run $f t=$t: exit $?"
    s1=$(date +%s.%N)
    w=$(awk -v a="$s0" -v b="$s1" 'BEGIN{printf "%.2f", b-a}')
    nt=$(grep -h -m1 'Ntried' "$o"/*.hist | awk '{print $NF}')
    [[ -z $t1 ]] && t1=$w
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$f" "$n" "$t" "$w" "$nt" \
      "$(awk -v w="$w" -v n="$nt" 'BEGIN{printf "%.2f", 1e6*w/n}')" "$(awk -v w="$w" -v n="$n" 'BEGIN{printf "%.3f", 1e3*w/n}')" \
      "$(awk -v a="$t1" -v w="$w" 'BEGIN{printf "%.2f", a/w}')" "$fus" >> "$TSV"
  done
done
column -t -s $'\t' "$TSV"
