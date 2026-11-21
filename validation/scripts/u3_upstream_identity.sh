#!/usr/bin/env bash
# validation/scripts/u3_upstream_identity.sh — task U.3: full-run identity of simc++ against the upstream SIMC 59ccdf6
# reference set reference/runs/u59_* (U.F2: the 294f662 reference inputs re-run with simc_gfortran 59ccdf6).
# Every run: simc++ <inp> --fortran-compat (compat target 59ccdf6, the default) with the RANLUX port seeded from the .inp,
# --tree-layout upstream. What is compared depends on what the Fortran produced (U.F2_report.md §1-§3):
#   FULL     .hist + .root     .hist without BEGIN/END Time (diff), .geni and .gen (cmp), every h10 column bitwise
#                              (simc_root_identity, upstream layout: same columns in the same order)
#   NONTU    .hist, no .root   (Nntu = 0 copies) .hist/.geni/.gen as above; simc++ must not write a .root either
#   A5       .bin > 80 columns (p15f2_poltar_pi0, U.F1 A4 / D-083 A5: ntu(80) written out of bounds, make_root_tree also
#                              overflows and keeps 1 entry) .hist/.geni/.gen as above; the simc++ tree is compared with
#                              the .bin itself (simc_root_identity --bin), columns 1-80 (inside ntu(80)) decisive,
#                              81-84 reported
#   STOP_A2  kaon + Nntu = 1   upstream STOPs at event 1 (results_write.f:317-319, 66 values for 65 tags, U.F1 A2); SIMC++
#                              does NOT reproduce that STOP (U-A2, corrected in both modes: the upstream layout declares
#                              the written sigcm column, 66 columns). Checked: simc++ runs (exit 0); its .hist/.geni/.gen
#                              (run name <run>_nontu) are identical to reference/runs/u59_<run>_nontu (Nntu = 0, same
#                              seed); its column list is the 65 tags of the Fortran .bin header plus sigcm after siglab;
#                              with decay off (no 1fae2d7 change) its tree equals the archived 294f662 tree
#                              reference/runs/<run>/<run>.root in the 54 common columns (--common-leaves --u1-upstream).
#   STOP     other STOP        simc++ must stop (exit != 0) with the Fortran STOP text (compare.txt "NEW STOP/crash
#                              text"); if the Fortran left a partial ntuple (sidis_d_k1_pip_fermi_rad_A: 182 events),
#                              simc++ is re-run with --ngen <events> and that tree must equal the partial Fortran tree.
#   SWEEP    u59_shipped_sweep the 18 legacy '.'-name inputs: simc++ must STOP with the Fortran STOP text and print the same
#                              "not registered" lines as the Fortran .out (in order); temp.inp must run and its tree equal
#                              u59_shipped_mc_neg_h_z37_ptbin (byte-identical .bin in the Fortran, P15.F1 W2).
#   PREFIX   FULL with env NGEN_PREFIX=N (CTest regression_upstream, for the 2-hour rho runs): simc++ --ngen N and the
#            N h10 entries must equal the first N of the reference tree (simc_root_identity --first N); the .hist/.geni/.gen
#            (whole-run sums) are not compared then.
# usage: u3_upstream_identity.sh <simc++> <simc_root_identity> <outdir> [u59 folder ...]   (default: every u59_* folder)
# env: JOBS (default 4, longest Fortran wall time first), REF (default $SIMC_REFERENCE/runs)
# writes <outdir>/<folder>/ (outputs, run.log, identity.txt, row.tsv) and <outdir>/summary.tsv
# (folder, run, class, verdict, seconds, note); exit 0 if every verdict is IDENTICAL / STOP_SAME / A2_NO_STOP /
# A5_BIN_IDENTICAL / SWEEP_SAME / PREFIX_IDENTICAL, 1 otherwise (77 when a single requested folder is absent: CTest skip).
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/simc_env.sh"   # SIMC_REFERENCE, ROOTSYS -> LD_LIBRARY_PATH
SIMC=$1; IDENT=$2; OUT=$3; shift 3
[[ -n ${REF-} ]] || simc_need SIMC_REFERENCE
REF=${REF:-$SIMC_REFERENCE/runs}; JOBS=${JOBS:-4}
LIST=("$@")
if [[ ${#LIST[@]} -eq 1 && ! -d "$REF/${LIST[0]}" ]]; then echo "reference $REF/${LIST[0]} missing: SKIP"; exit 77; fi
mkdir -p "$OUT"
if [[ ${#LIST[@]} -eq 0 ]]; then   # every u59_ folder, longest Fortran run first (.timing simc_wall)
  mapfile -t LIST < <(for d in "$REF"/u59_*; do
      w=$(grep -ho 'simc_wall=[0-9.]*' "$d"/*.timing 2>/dev/null | head -1 | cut -d= -f2); echo "${w:-0} $(basename "$d")"
    done | sort -k1,1gr | awk '{print $2}')
fi
echo '{"pipeline": {"rng": {"type": "ranlux"}}}' > "$OUT/ranlux.json"

one() {   # folder
  local f=$1 d="$REF/$1" w="$OUT/$1" run="" inp b t0 rc v cls note=""
  rm -rf "$w"; mkdir -p "$w"; t0=$(date +%s)
  local base=(--fortran-compat --pipeline "$OUT/ranlux.json" --tree-layout upstream --quiet)
  hist_same() {   # simc++ dir, ref dir, run name in the simc++ dir, ref run name -> prints 3 lines
    local x
    diff <(sed 2,3d "$1/$3.hist") <(sed 2,3d "$2/$4.hist") > "$w/hist.diff" 2>&1 && echo "hist: identical" || echo "hist: DIFFERENT"
    for x in geni gen; do cmp -s "$1/$3.$x" "$2/$4.$x" && echo "$x: identical" || echo "$x: DIFFERENT"; done
  }
  if [[ $f == u59_shipped_sweep ]]; then
    cls=SWEEP; local n=0 ok=0 l
    { for inp in "$d"/*.inp; do
        b=$(basename "$inp" .inp); [[ $b == extra_* || $b == fpi12_deltapp ]] && continue   # fpi12_deltapp: own folder
        mkdir -p "$w/$b"; "$SIMC" "$inp" -o "$w/$b" "${base[@]}" > "$w/$b/run.log" 2>&1; rc=$?
        if [[ $b == temp ]]; then
          if [[ $rc -eq 0 ]] && "$IDENT" "$w/$b/temp.root" "$REF/u59_shipped_mc_neg_h_z37_ptbin/mc_neg_h_z37_ptbin.root" > "$w/$b/tree.txt" 2>&1
          then echo "temp: runs, tree identical to u59_shipped_mc_neg_h_z37_ptbin ($(grep -o '[0-9]* entries x [0-9]* leaves' "$w/$b/tree.txt"))"
          else echo "temp: DIFFERENT (exit $rc)"; fi
          continue
        fi
        n=$((n + 1))
        l=$(grep -m1 '^ *STOP' "$d/$b.out" | sed 's/^ *//')
        if [[ $rc -ne 0 ]] && grep -qF "$l" "$w/$b/run.log" && \
           diff <(grep 'not registered' "$d/$b.out") <(grep 'not registered' "$w/$b/run.log") > /dev/null; then
          ok=$((ok + 1)); echo "$b: exit $rc, '$l', $(grep -c 'not registered' "$w/$b/run.log") not-registered lines as the Fortran"
        else echo "$b: DIFFERENT (exit $rc; Fortran '$l', $(grep -c 'not registered' "$d/$b.out") not-registered lines; simc++ $(grep -c 'not registered' "$w/$b/run.log"))"; fi
      done; } > "$w/identity.txt"
    run="(19 inputs)"; note="$ok/$n legacy inputs same STOP + not-registered list"
    grep -q DIFFERENT "$w/identity.txt" && v=DIFFERENT || v=SWEEP_SAME
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$f" "$run" "$cls" "$v" "$(( $(date +%s) - t0 ))" "$note" > "$w/row.tsv"; return
  fi
  for inp in "$d"/*.inp; do b=$(basename "$inp" .inp); [[ -f "$d/$b.hist" ]] && run=$b; done
  if [[ -z $run ]]; then for inp in "$d"/*.inp; do b=$(basename "$inp" .inp); [[ $b == extra_* ]] || run=$b; done; fi
  local stoptext; stoptext=$(sed -n 's/^NEW STOP\/crash text: //p' "$d/compare.txt" 2>/dev/null | head -1)
  [[ $stoptext == None ]] && stoptext=""
  if [[ -f "$d/$run.hist" && -f "$d/$run.root" ]]; then
    cls=FULL
    local nbin=0; [[ -f "$d/$run.bin" ]] && nbin=$(python3 -I -c 'import struct,sys; f=open(sys.argv[1],"rb"); f.read(4); print(struct.unpack("<i",f.read(4))[0])' "$d/$run.bin")
    [[ $nbin -gt 80 ]] && cls=A5
  elif [[ -f "$d/$run.hist" ]]; then cls=NONTU
  elif [[ $stoptext == *"number of output quantities"* ]]; then cls=STOP_A2
  elif [[ -n $stoptext ]]; then cls=STOP
  else cls=NO_REFERENCE; fi
  {
  case $cls in
  FULL|A5|NONTU)
    if [[ $cls == FULL && -n ${NGEN_PREFIX-} ]]; then
      cls=PREFIX
      "$SIMC" "$d/$run.inp" -o "$w" "${base[@]}" --ngen "$NGEN_PREFIX" > "$w/run.log" 2>&1; rc=$?
      echo "== $f ($run, $cls): simc++ RANLUX --fortran-compat (59ccdf6) --tree-layout upstream --ngen $NGEN_PREFIX, exit $rc"
      [[ $rc -ne 0 ]] && { echo "run: ERROR"; tail -3 "$w/run.log"; }
      [[ $rc -eq 0 ]] && { "$IDENT" --first "$NGEN_PREFIX" "$w/$run.root" "$d/$run.root" > "$w/tree.txt" 2>&1 \
        && echo "tree: first $NGEN_PREFIX entries identical ($(grep -o '[0-9]* entries x [0-9]* leaves' "$w/tree.txt" | head -1))" \
        || { echo "tree: DIFFERENT"; head -8 "$w/tree.txt"; }; }
    else
    "$SIMC" "$d/$run.inp" -o "$w" "${base[@]}" > "$w/run.log" 2>&1; rc=$?
    echo "== $f ($run, $cls): simc++ RANLUX --fortran-compat (59ccdf6) --tree-layout upstream, exit $rc"
    [[ $rc -ne 0 ]] && { echo "run: ERROR"; tail -3 "$w/run.log"; }
    if [[ $rc -eq 0 ]]; then
      hist_same "$w" "$d" "$run" "$run"
      if [[ $cls == NONTU ]]; then
        [[ -f "$w/$run.root" ]] && echo "tree: DIFFERENT (simc++ wrote a .root, the Nntu = 0 reference has none)" || echo "tree: none written (Nntu = 0), as the Fortran"
      elif [[ $cls == A5 ]]; then
        "$IDENT" --bin "$w/$run.root" "$d/$run.bin" 1 80 > "$w/tree.txt" 2>&1 && echo "tree vs .bin columns 1-80: identical" || echo "tree vs .bin: DIFFERENT"
        sed 's/^/  /' "$w/tree.txt"
        "$IDENT" "$w/$run.root" "$d/$run.root" > "$w/tree_root.txt" 2>&1
        echo "  (Fortran .root: $(grep -m1 -o 'entries [0-9]* vs [0-9]*' "$w/tree_root.txt" || echo 'see tree_root.txt'); make_root_tree overflow, D-083 A5)"
      else
        "$IDENT" "$w/$run.root" "$d/$run.root" > "$w/tree.txt" 2>&1 && echo "tree: identical ($(grep -o '[0-9]* entries x [0-9]* leaves' "$w/tree.txt" | head -1))" || { echo "tree: DIFFERENT"; head -8 "$w/tree.txt"; }
      fi
    fi
    fi ;;
  STOP_A2)
    local nt="$REF/${f}_nontu"
    "$SIMC" "$d/$run.inp" -o "$w" "${base[@]}" --run-name "${run}_nontu" > "$w/run.log" 2>&1; rc=$?
    echo "== $f ($run, $cls): Fortran: '$stoptext' at event 1 (U.F1 A2); simc++ exit $rc (U-A2: no STOP, 66-column tree)"
    [[ $rc -ne 0 ]] && { echo "run: ERROR"; tail -3 "$w/run.log"; }
    if [[ $rc -eq 0 ]]; then
      if [[ -d $nt ]]; then echo "vs ${f}_nontu (Nntu = 0, same seed):"; hist_same "$w" "$nt" "${run}_nontu" "${run}_nontu"
      else echo "hist: DIFFERENT (no ${f}_nontu reference)"; fi
      local tags cols
      tags=$(python3 -I -c '
import struct,sys
d=open(sys.argv[1],"rb").read();p=0;r=[]
while p<len(d):
    n=struct.unpack("<i",d[p:p+4])[0];r.append(d[p+4:p+4+n]);p+=8+n
n=struct.unpack("<i",r[0][:4])[0];t=[x.decode(errors="replace").strip() for x in r[1:1+n]]
i=t.index("siglab");print(" ".join(t[:i+1]+["sigcm"]+t[i+1:]))' "$d/$run.bin")
      cols=$("$ROOTSYS/bin/root" -l -b -q -e "TFile f(\"$w/${run}_nontu.root\"); auto t=(TTree*)f.Get(\"h10\"); for(auto l:*t->GetListOfLeaves()) printf(\"%s \",l->GetName());" 2>/dev/null | tail -1 | sed 's/ *$//')
      [[ "$tags" == "$cols" ]] && echo "columns: Fortran .bin header (65 tags) + sigcm after siglab = simc++ h10 ($(wc -w <<< "$cols") leaves)" || echo "columns: DIFFERENT (Fortran+sigcm: $tags | simc++: $cols)"
      if grep -qiE '^ *doing_decay *= *0' "$d/$run.inp" && [[ -f "$REF/$run/$run.root" ]]; then
        "$IDENT" --common-leaves --u1-upstream "$w/${run}_nontu.root" "$REF/$run/$run.root" > "$w/tree_old.txt" 2>&1 \
          && echo "tree vs archived 294f662 tree (decay off): identical in the common columns ($(grep -o '[0-9]* common leaves' "$w/tree_old.txt"), $(grep -o '[0-9]* entries' "$w/tree_old.txt" | head -1))" \
          || { echo "tree vs archived 294f662 tree: DIFFERENT"; head -5 "$w/tree_old.txt"; }
      else echo "tree: decay on, the 294f662 tree differs by 1fae2d7 and 59ccdf6 writes none; not compared"; fi
    fi ;;
  STOP)
    "$SIMC" "$d/$run.inp" -o "$w" "${base[@]}" > "$w/run.log" 2>&1; rc=$?
    echo "== $f ($run, $cls): Fortran STOP '$stoptext'; simc++ exit $rc"
    local key=${stoptext%%[0-9-]*}; key=$(sed 's/ *$//' <<< "$key")   # the text before a printed number
    if [[ $rc -ne 0 ]] && grep -qF "$key" "$w/run.log"; then
      echo "stop: same STOP ($(grep -m1 -F "$key" "$w/run.log" | sed 's/^ *//'))"
      if [[ -f "$d/$run.root" ]]; then
        local nev; nev=$("$ROOTSYS/bin/root" -l -b -q -e "TFile f(\"$d/$run.root\"); printf(\"%lld\n\", ((TTree*)f.Get(\"h10\"))->GetEntries());" 2>/dev/null | tail -1)
        mkdir -p "$w/ngen$nev"
        "$SIMC" "$d/$run.inp" -o "$w/ngen$nev" "${base[@]}" --ngen "$nev" > "$w/ngen$nev/run.log" 2>&1
        "$IDENT" "$w/ngen$nev/$run.root" "$d/$run.root" > "$w/tree.txt" 2>&1 \
          && echo "partial tree: the Fortran wrote $nev events before the STOP; simc++ --ngen $nev tree identical ($(grep -o '[0-9]* entries x [0-9]* leaves' "$w/tree.txt"))" \
          || { echo "partial tree: DIFFERENT"; head -5 "$w/tree.txt"; }
      fi
    else echo "stop: DIFFERENT"; tail -3 "$w/run.log"; fi ;;
  *) echo "== $f: no reference outputs"; rc=1 ;;
  esac; } > "$w/identity.txt"
  if grep -q 'DIFFERENT\|ERROR' "$w/identity.txt" || [[ $cls == NO_REFERENCE ]]; then v=DIFFERENT
  else case $cls in FULL|NONTU) v=IDENTICAL ;; PREFIX) v=PREFIX_IDENTICAL ;; A5) v=A5_BIN_IDENTICAL ;; STOP_A2) v=A2_NO_STOP ;; STOP) v=STOP_SAME ;; esac; fi
  grep -q '^run: ERROR' "$w/identity.txt" && v=ERROR
  case $cls in A5) note="cols 81-84: $(grep -o 'columns outside 1-80:.*' "$w/tree.txt" | cut -d: -f2-)" ;;
               STOP_A2) note="U-A2: simc++ runs" ;; STOP) note="$stoptext" ;; esac
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$f" "$run" "$cls" "$v" "$(( $(date +%s) - t0 ))" "$note" > "$w/row.tsv"
  rm -f "$w"/*.root.tmp
}
export -f one; export SIMC IDENT OUT REF NGEN_PREFIX
printf '%s\n' "${LIST[@]}" | xargs -P "$JOBS" -I{} bash -c 'one {}'
{ printf 'folder\trun\tclass\tverdict\tseconds\tnote\n'; for f in "${LIST[@]}"; do cat "$OUT/$f/row.tsv"; done; } > "$OUT/summary.tsv"
column -t -s$'\t' "$OUT/summary.tsv"
echo; cut -f3,4 "$OUT/summary.tsv" | tail -n +2 | sort | uniq -c
[[ ${#LIST[@]} -eq 1 ]] && cat "$OUT/${LIST[0]}/identity.txt"
! grep -q -P "\t(DIFFERENT|ERROR)\t" "$OUT/summary.tsv"
