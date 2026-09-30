// include/simc/core/RunSummary.h
//
// End-of-run numbers handed to the output sinks. Fortran: simc.f l.360-400
// (normalisation) and simc.f:report (l.644-1142) / .geni writing (l.462-532).
#ifndef SIMC_CORE_RUNSUMMARY_H
#define SIMC_CORE_RUNSUMMARY_H

#include <map>
#include <string>
#include <vector>

#include "simc/core/Json.h"

namespace simc {

/// Per-aperture stop counters of one spectrometer. Keys are documented in docs/OUTPUT_FORMAT.md.
struct StopCounters {
    std::string spectrometer;            ///< "hms", "shms", "sos", "hrsr", "hrsl", "calo"
    std::string role;                    ///< "electron", "hadron" or "both"
    std::map<std::string, long> counts;  ///< e.g. {"trials": 3628935, "slit_hor": 2789053, ...}
};

struct RunSummary {
    std::string run_name;
    std::string reaction_title;
    long ngen{0};                        ///< Ngen (request)
    long ntried{0};                      ///< Ntried
    long ncontribute{0};                 ///< Ncontribute
    long ncontribute_no_rad_proton{0};   ///< Nco_no_rad_prot
    long nevent{0};                      ///< Fortran nevent (successes, or tries if ngen<0)
    double wtcontribute{0};              ///< after *normfac (simc.f:398); .hist prints wtcontribute/nevent
    double sum_sigcc{0};                 ///< .hist prints sum_sigcc/nevent as AVERAGE.sigcc
    double central_sigcc{0};
    double charge{0};                    ///< mC
    double targetfac{0};
    double luminosity{0};                ///< ub^-1
    double genvol{0};
    double normfac{0};
    int random_seed{0};
    std::vector<StopCounters> stops;
    /// Optional summary sections already in the summary.json layout (kinematics, target, central,
    /// radiative, recon_summary, limits, ...). Merged verbatim into the output.
    Json extra = Json::object();
};

}  // namespace simc

#endif  // SIMC_CORE_RUNSUMMARY_H
