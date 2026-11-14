#!/usr/bin/env bash
# validation/scripts/p153_golden.sh — task P15.3: before/after proof that a refactor changes no output.
#   make  <simc++> <simc_compare> <outdir> : runs a fixed set of inputs (every reaction family and arm type, the target
#          field) in both modes (seed 4242, mt64, short ngen) and one simc_compare, keeping every output
#   check <dirA> <dirB> <simc_root_identity> : cmp of every file; .hist/.log/.summary.* without the begin/end time lines, .root
#          trees with simc_root_identity (the ROOT file header holds a date), the compare PDF by size only; exit 0 if identical
set -uo pipefail
REF=/home/ayerbe/NewDisk/reference/runs
RUNS=(eep_h_elastic_rad_A:4000 acc_phsp_hms_shms_A:4000 eep_d_rad_A:2000 eep_fe56_rad_A:2000 eep_fe56_bh_rad_A:1000
      pion_h_k1_rad_dec1_A:2000 pion_d_k1_pim_rad_dec1_A:1000 kaon_h_k1_rad_dec1_A:1000 kaon_d_k1_sigm_rad_dec1_A:500
      delta_h_d1_rad_A:2000 sidis_h_k1_pip_rad_A:2000 rho_h_k1_rad_A:100 p14f2_excl_h_ng2_A:300 p14f2_sidis_h_A:200
      p14f1_eep_sos_A:2000 p14f1_hrs_LeRp_rad_A:2000 p15f1_eep_smm1_A:2000 p15f1_pi0_shms_A:300 p15f2_poltar_pi0:0
      shipped_poltar_exclusivepion_example:300)
mode=$1; shift
if [[ $mode == make ]]; then
  SIMC=$1; CMP=$2; OUT=$3; mkdir -p "$OUT"
  for e in "${RUNS[@]}"; do
    f=${e%%:*}; n=${e##*:}; d=$REF/$f; [[ -d $d ]] || continue
    inp=""; for i in "$d"/*.inp; do b=$(basename "$i" .inp); [[ -f $d/$b.hist ]] && inp=$i; done
    [[ -n $inp ]] || continue
    for m in compat corrected; do
      x=(--seed 4242 --quiet); [[ $n != 0 ]] && x+=(--ngen "$n"); [[ $m == compat ]] && x+=(--fortran-compat)
      "$SIMC" "$inp" -o "$OUT/$f/$m" "${x[@]}" > "$OUT/$f.$m.log" 2>&1 || echo "run $f $m: exit $?" >> "$OUT/errors.txt"
    done
  done
  "$CMP" --a "$REF/eep_h_elastic_rad_A" --b "$OUT/eep_h_elastic_rad_A/corrected/eep_h_elastic_rad_A.summary.json" \
     --mode-b corrected --out "$OUT/compare" > "$OUT/compare.log" 2>&1
  ls -R "$OUT" | wc -l
elif [[ $mode == check ]]; then
  A=$1; B=$2; IDENT=$3; bad=0; n=0
  while IFS= read -r f; do
    n=$((n+1)); a=$A/$f; b=$B/$f
    [[ -f $b ]] || { echo "missing in B: $f"; bad=1; continue; }
    case $f in
      *.hist|*.log|*.out|*.txt|*.json)   # (time stamps and the output directory names excluded)
        flt() { grep -v -E '^ (BEGIN|END) Time:|^ *"?(begin|end)_time"? *[=:]' "$1" | sed "s|$2|<DIR>|g"; }
        diff -q <(flt "$a" "$A") <(flt "$b" "$B") > /dev/null || { echo "DIFF $f"; bad=1; } ;;
      *.root) "$IDENT" "$a" "$b" > /dev/null 2>&1 || { echo "DIFF $f"; bad=1; } ;;
      *.pdf) [[ $(stat -c %s "$a") == $(stat -c %s "$b") ]] || { echo "PDF size differs $f"; } ;;
      *) cmp -s "$a" "$b" || { echo "DIFF $f"; bad=1; } ;;
    esac
  done < <(cd "$A" && find . -type f | sort)
  echo "$n files compared: $([[ $bad == 0 ]] && echo IDENTICAL || echo DIFFERENT)"
  exit $bad
fi
