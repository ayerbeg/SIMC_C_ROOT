#!/usr/bin/env bash
# validation/scripts/regen_baselines.sh - regenerate the C++ regression baselines on the current toolchain (task P10.4, D-044).
#
# Usage: regen_baselines.sh [--commit REV] [--tag TAG] [--jobs N] [--only SUBSTR]
#   --commit REV   commit to build (clean detached worktree, Release, like tools/check_all.sh); default HEAD
#   --tag TAG      toolchain tag in the folder names; default glibc<major><minor> of this machine (e.g. glibc243)
#   --jobs N       simc++ runs in parallel (default 4)
#   --only SUBSTR  only the baselines whose folder suffix contains SUBSTR
# Writes NEW folders SIMC++_Project/baselines/<rev>_<tag>_<suffix>/ (refuses to overwrite an existing one) with the
# outputs of exactly the command the regression test runs (tests/regression/CMakeLists.txt, output_identity.sh), plus
# README.txt (commit, command, toolchain line read by tests/toolchain_guard.sh) and MD5SUMS. Old folders are untouched.
set -uo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
SRC_DIR="$PROJECT_DIR/SIMC++_Claude"
BASEDIR="$PROJECT_DIR/baselines"
REF=/home/ayerbe/NewDisk/reference
FORTRAN=/home/ayerbe/NewDisk/simc_gfortran
ROOT_DIR="${ROOTSYS:-/home/ayerbe/Builds/ROOT6_install}"
export LD_LIBRARY_PATH="$ROOT_DIR/lib:/usr/local/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
FP="$SRC_DIR/tools/toolchain_fingerprint.sh"
REV=HEAD; TAG=""; JOBS=4; ONLY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --commit) REV="$2"; shift ;;
    --tag) TAG="$2"; shift ;;
    --jobs) JOBS="$2"; shift ;;
    --only) ONLY="$2"; shift ;;
    -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
  esac; shift
done
[[ -n "$TAG" ]] || TAG="glibc$("$FP" --field glibc | tr -d .)"
SHORT="$(git -C "$SRC_DIR" rev-parse --short "$REV")" || exit 2

# suffix | input | extra simc++ args | run name | task/test note
CONFIGS=(
  "rad_s4242|$FORTRAN/infiles/eep_h_elastic_rad_A.inp|--fortran-compat|dir_rad_s4242|regression_identity_rad (was f9cea51_rad_s4242)"
  "norad_s4242|$FORTRAN/infiles/eep_h_elastic_norad_A.inp|--fortran-compat|dir_norad_s4242|regression_identity_norad (was f9cea51_norad_s4242)"
  "phsp_s4242|$FORTRAN/infiles/acc_phsp_hms_shms_A.inp|--fortran-compat|dir_phsp_s4242|regression_identity_phsp (was f9cea51_phsp_s4242)"
  "corrected_rad_s4242|$FORTRAN/infiles/eep_h_elastic_rad_A.inp||dir_rad_s4242|regression_corrected_rad (was 08812e6_corrected_rad_s4242)"
  "corrected_norad_s4242|$FORTRAN/infiles/eep_h_elastic_norad_A.inp||dir_norad_s4242|regression_corrected_norad (was 08812e6_corrected_norad_s4242)"
  "corrected_phsp_s4242|$FORTRAN/infiles/acc_phsp_hms_shms_A.inp||dir_phsp_s4242|regression_corrected_phsp (was 08812e6_corrected_phsp_s4242)"
)
for t in d_rad fe56_rad; do
  CONFIGS+=("eep_${t}_n20k_s4242|$REF/runs/eep_${t}_A/eep_${t}_A.inp|--fortran-compat --ngen 20000|dir_eep_${t}_s4242|regression_identity_eep_${t} (was 508633a_eep_${t}_n20k_s4242)"
            "corrected_eep_${t}_n20k_s4242|$REF/runs/eep_${t}_A/eep_${t}_A.inp|--ngen 20000|dir_eep_${t}_s4242|regression_corrected_eep_${t} (was 508633a_corrected_eep_${t}_n20k_s4242)")
done
for t in h_k1_rad_dec1 d_k1_pim_rad_dec1; do
  CONFIGS+=("pion_${t}_n20k_s4242|$REF/runs/pion_${t}_A/pion_${t}_A.inp|--fortran-compat --ngen 20000|dir_pion_${t}_s4242|regression_identity_pion_${t} (was 3b27d5c_pion_${t}_n20k_s4242)"
            "corrected_pion_${t}_n20k_s4242|$REF/runs/pion_${t}_A/pion_${t}_A.inp|--ngen 20000|dir_pion_${t}_s4242|regression_corrected_pion_${t} (was 3b27d5c_corrected_pion_${t}_n20k_s4242)")
done
for t in h_k1_rad_dec1 d_k1_sigm_rad_dec1; do
  CONFIGS+=("kaon_${t}_n20k_s4242|$REF/runs/kaon_${t}_A/kaon_${t}_A.inp|--fortran-compat --ngen 20000|dir_kaon_${t}_s4242|regression_identity_kaon_${t} (was 1275770_kaon_${t}_n20k_s4242)"
            "corrected_kaon_${t}_n20k_s4242|$REF/runs/kaon_${t}_A/kaon_${t}_A.inp|--ngen 20000|dir_kaon_${t}_s4242|regression_corrected_kaon_${t} (was 1275770_corrected_kaon_${t}_n20k_s4242)")
done

WT="$(mktemp -d "${TMPDIR:-/tmp}/simc_regen_wt.XXXXXX")"
cleanup() { git -C "$SRC_DIR" worktree remove --force "$WT" >/dev/null 2>&1; git -C "$SRC_DIR" worktree prune; rm -rf "$WT.build"; }
trap cleanup EXIT
git -C "$SRC_DIR" worktree add -q --detach "$WT" "$REV" || exit 2
echo "building $SHORT in $WT.build ..."
if ! { cmake -S "$WT" -B "$WT.build" -DCMAKE_PREFIX_PATH="$ROOT_DIR" -DCMAKE_BUILD_TYPE=Release \
         -DSIMC_REFERENCE_DIR=$REF -DSIMC_FORTRAN_DIR=$FORTRAN > "$WT.build.log" 2>&1 &&
       cmake --build "$WT.build" -j"$(nproc)" --target simcpp >> "$WT.build.log" 2>&1; }; then
  echo "build failed, see $WT.build.log"; exit 2
fi
SIMCPP="$WT.build/bin/simc++"
TOOLCHAIN="$("$FP" --binary "$SIMCPP")"
echo "toolchain: $TOOLCHAIN"

run_one() {
  local suffix input extra run note out
  IFS='|' read -r suffix input extra run note <<<"$1"
  out="$BASEDIR/${SHORT}_${TAG}_${suffix}"
  if [[ -e "$out" ]]; then echo "exists, not touched: $out"; return 1; fi
  mkdir -p "$out"
  # shellcheck disable=SC2086
  if ! "$SIMCPP" "$input" --seed 4242 --run-name "$run" -o "$out" --quiet $extra > "$out/.run.log" 2>&1; then
    echo "FAILED: $suffix (log $out/.run.log)"; return 1
  fi
  rm -f "$out/.run.log"
  local cmd="simc++ ${input#/home/ayerbe/NewDisk/} --seed 4242 --run-name $run -o <dir> --quiet${extra:+ $extra}"
  {
    echo "simc++ (commit $SHORT, task P10.4, Release build of the committed tree in a clean worktree by cpp-translator,"
    echo "validation/scripts/regen_baselines.sh) on: $cmd"
    echo "Produced $(date -Iseconds) for the output-identity test $note."
    echo "Regenerated because the bitwise C++ baselines are tied to the toolchain (D-044, glibc 2.43 / gcc 16.2.1)."
    echo "toolchain: $TOOLCHAIN"
  } > "$out/README.txt"
  (cd "$out" && md5sum $(ls | grep -v -E '^(README\.txt|MD5SUMS)$') > MD5SUMS)
  echo "done: $out"
}
# f9cea51_K2K4_s4242 (P7.3 L2 runs of the K2-K4 kinematics; validation test compare_both_trees_empty reads
# cpp_acc_phsp_K2_A.summary.json): 7 runs cpp_<tag>.* plus their stdout cpp_<tag>.log, Fortran-compatible mode.
run_k2k4() {
  local out="$BASEDIR/${SHORT}_${TAG}_K2K4_s4242" t
  if [[ -e "$out" ]]; then echo "exists, not touched: $out"; return 1; fi
  mkdir -p "$out"
  for t in acc_phsp_K2_A eep_h_K4_norad_A eep_h_K2_norad_A eep_h_K3_norad_A eep_h_K2_rad_A eep_h_K4_rad_A eep_h_K3_rad_A; do
    if ! "$SIMCPP" "$REF/runs/$t/$t.inp" --seed 4242 --run-name "cpp_$t" -o "$out" --fortran-compat > "$out/cpp_$t.log" 2>&1; then
      echo "FAILED: K2K4 $t"; return 1
    fi
  done
  {
    echo "simc++ (commit $SHORT, task P10.4, Release build of the committed tree in a clean worktree by cpp-translator,"
    echo "validation/scripts/regen_baselines.sh) on: simc++ reference/runs/<tag>/<tag>.inp --seed 4242 --run-name cpp_<tag> -o <dir>"
    echo "--fortran-compat > cpp_<tag>.log, for <tag> = acc_phsp_K2_A eep_h_K{2,3,4}_{norad,rad}_A."
    echo "Produced $(date -Iseconds); replaces f9cea51_K2K4_s4242 (P7.3 L2 runs) for the validation test compare_both_trees_empty."
    echo "Regenerated because the C++ baselines are tied to the toolchain (D-044, glibc 2.43 / gcc 16.2.1)."
    echo "toolchain: $TOOLCHAIN"
  } > "$out/README.txt"
  (cd "$out" && md5sum $(ls | grep -v -E '^(README\.txt|MD5SUMS)$') > MD5SUMS)
  echo "done: $out"
}
export -f run_one run_k2k4
export SIMCPP SHORT TAG BASEDIR TOOLCHAIN REF
status=0
k2pid=""
if [[ -z "$ONLY" || "K2K4_s4242" == *"$ONLY"* ]]; then run_k2k4 & k2pid=$!; fi
for c in "${CONFIGS[@]}"; do
  [[ -z "$ONLY" || "${c%%|*}" == *"$ONLY"* ]] || continue
  printf '%s\0' "$c"
done | xargs -0 -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {} || status=1
if [[ -n "$k2pid" ]]; then wait "$k2pid" || status=1; fi
exit $status
