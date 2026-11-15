# validation/scripts/simc_env.sh — sourced (not executed) by the validation scripts in this directory: the locations
# of the validation data and of the project folders, without any machine-specific path in the repository.
#
# Sets (variables already in the environment take precedence over everything below):
#   SIMC_SRC_DIR      the SIMC++ repository (two levels above this file)
#   SIMC_PROJECT_DIR  the project folder that holds the repository, tools/run_simc.sh, simc_runs/, baselines/ and the
#                     untracked check_all.env (default: the parent of the repository)
#   SIMC_REFERENCE    Fortran reference data (runs/, functions/, comparisons/)        } from the environment, else from
#   SIMC_FORTRAN      the Fortran SIMC tree (infiles/, data files)                     } $SIMC_PROJECT_DIR/check_all.env
#   SIMC_BASELINES    the C++ regression baselines                                     } (shell syntax, as tools/check_all.sh)
#   ROOTSYS           the ROOT 6 install; its lib/ is prepended to LD_LIBRARY_PATH
# and defines simc_need VAR... : exit 2 with a message when one of the variables is empty (call it for the data a script
# really uses, so that a script that only needs the reference data does not fail on a missing Fortran tree).
_simc_env_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIMC_SRC_DIR="${SIMC_SRC_DIR:-$(cd "$_simc_env_here/../.." && pwd)}"
SIMC_PROJECT_DIR="${SIMC_PROJECT_DIR:-$(cd "$SIMC_SRC_DIR/.." && pwd)}"
if [[ -f $SIMC_PROJECT_DIR/check_all.env ]]; then
  _simc_r=${SIMC_REFERENCE-} _simc_f=${SIMC_FORTRAN-} _simc_b=${SIMC_BASELINES-} _simc_rs=${ROOTSYS-}
  # shellcheck disable=SC1091
  source "$SIMC_PROJECT_DIR/check_all.env"
  SIMC_REFERENCE=${_simc_r:-${SIMC_REFERENCE-}}; SIMC_FORTRAN=${_simc_f:-${SIMC_FORTRAN-}}
  SIMC_BASELINES=${_simc_b:-${SIMC_BASELINES-}}; ROOTSYS=${_simc_rs:-${ROOTSYS-}}
fi
export SIMC_SRC_DIR SIMC_PROJECT_DIR SIMC_REFERENCE="${SIMC_REFERENCE-}" SIMC_FORTRAN="${SIMC_FORTRAN-}" \
       SIMC_BASELINES="${SIMC_BASELINES-}"
if [[ -n ${ROOTSYS-} ]]; then
  export ROOTSYS
  export LD_LIBRARY_PATH="$ROOTSYS/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi
simc_need() {
  local v
  for v in "$@"; do
    [[ -n ${!v-} ]] || { echo "$(basename "$0"): set $v (environment or $SIMC_PROJECT_DIR/check_all.env)" >&2; exit 2; }
  done
}
