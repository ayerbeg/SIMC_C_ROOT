/**
 * @file
 * @ingroup core
 * @brief simc::RunSummary: the end-of-run numbers handed to the output sinks (counts, normalisation, stop counters).
 * @note Origin: Fortran SIMC simc.f l.360-400 (normalisation), simc.f:report (l.644-1142) and the `.geni` writing
 *       (l.462-532).
 */
#ifndef SIMC_CORE_RUNSUMMARY_H
#define SIMC_CORE_RUNSUMMARY_H

#include <map>
#include <string>
#include <vector>

#include "simc/core/Json.h"

namespace simc {

/** @addtogroup core
 *  @{ */

/// Per-aperture stop counters of one spectrometer. Keys are documented in docs/OUTPUT_FORMAT.md.
struct StopCounters {
    std::string spectrometer;            ///< "hms", "shms", "sos", "hrsr", "hrsl", "calo"
    std::string role;                    ///< "electron", "hadron" or "both"
    std::map<std::string, long> counts;  ///< e.g. {"trials": 3628935, "slit_hor": 2789053, ...}
};

/**
 * @brief End-of-run numbers: event counts, normalisation, luminosity and the aperture stop counters.
 *
 * Filled by the driver after the event loop (simc::driver::SimcRun::finish) and passed to every
 * simc::IOutputSink::endRun. The normalisation follows the original program: the event weights are multiplied by
 * #normfac, and #wtcontribute holds the normalised weight sum. Plain value type.
 */
struct RunSummary {
    std::string run_name;                ///< Run name (file stem of the outputs).
    std::string reaction_title;          ///< Reaction title line (simc::ReactionTitle).
    long ngen{0};                        ///< Ngen (request)
    long ntried{0};                      ///< Ntried
    long ncontribute{0};                 ///< Ncontribute
    long ncontribute_no_rad_proton{0};   ///< Nco_no_rad_prot
    long nevent{0};                      ///< Fortran nevent (successes, or tries if ngen<0)
    double wtcontribute{0};              ///< after *normfac (simc.f:398); .hist prints wtcontribute/nevent
    double sum_sigcc{0};                 ///< .hist prints sum_sigcc/nevent as AVERAGE.sigcc
    double central_sigcc{0};             ///< Cross section at the central kinematics (units of the physics model).
    double charge{0};                    ///< mC
    double targetfac{0};                 ///< Target factor: luminosity = charge / targetfac (simc.f:94-95).
    double luminosity{0};                ///< ub^-1
    double genvol{0};                    ///< Generation volume (product of the generated ranges).
    double normfac{0};                   ///< Normalisation factor: luminosity x genvol / ntried.
    int random_seed{0};                  ///< Seed of the run (as given in the input).
    std::vector<StopCounters> stops;     ///< Aperture stop counters, one entry per spectrometer arm.
    /// Optional summary sections already in the summary.json layout (kinematics, target, central,
    /// radiative, recon_summary, limits, ...). Merged verbatim into the output.
    Json extra = Json::object();
};

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_RUNSUMMARY_H
