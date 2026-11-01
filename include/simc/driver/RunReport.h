// include/simc/driver/RunReport.h
//
// End-of-run bookkeeping of the SIMC main program and its text reports (task P7.2).
//   RunTotals      simc.f main-loop accumulators (l.160-352): ntried, nevent, ncontribute, npasscuts, wtcontribute,
//                  sum_sigcc, sumerr/sumerr2 (sums_twoarm), contrib (contribtype) and the slop ranges
//   Normalisation  simc.f l.360-432: normfac, genvol, the resolutions aveerr/resol; targetfac/luminosity (l.94-101)
//   FormatHist()   simc.f:report (l.644-1142): the .hist text, line for line in the gfortran format
//   FormatGeni()   simc.f l.462-532: the .geni stop-counter text
//   MakeRunSummary()  the RunSummary (and its "extra" sections) handed to the output sinks (summary JSON)
#ifndef SIMC_DRIVER_RUNREPORT_H
#define SIMC_DRIVER_RUNREPORT_H

#include <map>
#include <string>
#include <vector>

#include "simc/core/RunSetup.h"
#include "simc/core/RunSummary.h"
#include "simc/driver/Central.h"
#include "simc/reactions/EventChain.h"

namespace simc {
namespace driver {

/// Fortran: type sums_twoarm (structures.inc) — per arm delta, xptar, yptar, ytar.
struct TwoArmSums {
    struct Arm { double delta{0}, xptar{0}, yptar{0}, ytar{0}; };
    Arm e, p;
};

/// Accumulators of the main loop (simc.f l.160-352 and COMMON /simulate/ counters).
struct RunTotals {
    long ntried{0};                      ///< tries
    long nevent{0};                      ///< successes (ngen > 0) or tries (ngen < 0)
    long ncontribute{0};
    long ncontribute_no_rad_proton{0};
    long npasscuts{0};
    double wtcontribute{0};              ///< sum of main%weight over pass_cuts successes (before *normfac)
    double sum_sigcc{0};                 ///< sum of main%sigcc over successes
    TwoArmSums sumerr, sumerr2;          ///< pass_cuts successes
    reactions::Contrib contrib;          ///< "found" ranges (limits_update)
    Slops slop;                          ///< slop%MC / slop%total ranges (limits_update; .used from the setup)
};

/// simc.f l.94-101 and l.360-432.
struct Normalisation {
    double targetfac{0}, luminosity{0};  ///< from the setup (simc.f:94-101)
    double charge{0};                    ///< EXPER%charge [mC]
    double genvol{0}, normfac{0};
    double wtcontribute{0};              ///< wtcontribute*normfac (simc.f:398)
    TwoArmSums aveerr, resol;            ///< zero unless npasscuts > 1 (simc.f:405)
};

/// Everything simc.f:report prints that is not in the RunSetup.
struct RunReport {
    std::string run_name;
    std::string begin_time, end_time;    ///< ctime() strings (simc.f:162, 357)
    long ngen{0};
    int random_seed{0};
    std::string random_state_file;
    CentralEvent central;
    RunTotals totals;
    Normalisation norm;
    std::vector<StopCounters> stops;     ///< one entry per spectrometer type in use
};

/// Fortran: simc.f l.360-432 (normalisation and resolutions). `genvol` from IReaction::genVolume().
Normalisation Normalise(const RunSetup& setup, const RunTotals& t, double genvol);

/// Fortran: simc.f:report(iun, ...) (l.644-1142). Text identical to the .hist the Fortran writes.
/// Throws std::runtime_error where the Fortran report STOPs (unknown reaction).
std::string FormatHist(const RunSetup& setup, const RunReport& r);

/// Fortran: simc.f l.462-532 (.geni). The HMS, SOS, HRS-R/L and SHMS blocks are printed in the Fortran format (calo arms are
/// not supported by the driver yet); any other ISpectrometer plug-in gets a generic "key value" block.
std::string FormatGeni(const RunSetup& setup, const RunReport& r);

/// The RunSummary for the output sinks, with the optional summary sections (central, target, radiative,
/// recon_summary, limits, slop, kinematics) in the summary.json layout (docs/OUTPUT_FORMAT.md §2).
RunSummary MakeRunSummary(const RunSetup& setup, const RunReport& r);

/// ctime()-style time string as simc.f writes it (gfortran ctime: "Thu Oct  1 01:23:05 2026").
std::string FortranCtime(long epoch_seconds);

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_RUNREPORT_H
