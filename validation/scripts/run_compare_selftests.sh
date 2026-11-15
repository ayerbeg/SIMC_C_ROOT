#!/usr/bin/env bash
# Run the simc_compare self-tests and store the outputs next to the references:
#   reference/comparisons/F_vs_F/<tag>/            Fortran A vs Fortran B (must PASS)
#   reference/comparisons/sensitivity/<name>/      rad vs norad, hsdelta+0.1 (must FAIL)
# Usage: validation/scripts/run_compare_selftests.sh [BUILD_DIR]
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:-$HERE/build}"
source "$(dirname "${BASH_SOURCE[0]}")/simc_env.sh"   # SIMC_PROJECT_DIR, SIMC_REFERENCE, ROOTSYS -> LD_LIBRARY_PATH
simc_need SIMC_REFERENCE
R=$SIMC_REFERENCE/runs
C=$SIMC_REFERENCE/comparisons
CFG="$HERE/config/compare_default.json"
EXE="$BUILD/simc_compare"
rc=0
for tag in eep_h_elastic_norad eep_h_elastic_rad acc_phsp_hms_shms; do
  "$EXE" --a "$R/${tag}_A/${tag}_A.hist" --b "$R/${tag}_B/${tag}_B.hist" \
         --label-a "${tag}_A" --label-b "${tag}_B" --out "$C/F_vs_F/$tag" --config "$CFG" --expect PASS || rc=1
done
"$EXE" --a "$R/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist" --b "$R/eep_h_elastic_norad_A/eep_h_elastic_norad_A.hist" \
       --label-a rad_A --label-b norad_A --out "$C/sensitivity/rad_A_vs_norad_A" --config "$CFG" --expect FAIL || rc=1
"$EXE" --a "$R/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist" --b "$R/eep_h_elastic_rad_A/eep_h_elastic_rad_A.hist" \
       --perturb hsdelta:0.1 --label-a rad_A --label-b "rad_A(hsdelta+0.1)" \
       --out "$C/sensitivity/rad_A_vs_rad_A_hsdelta+0.1" --config "$CFG" --expect FAIL --expect-flagged hsdelta || rc=1
exit $rc
