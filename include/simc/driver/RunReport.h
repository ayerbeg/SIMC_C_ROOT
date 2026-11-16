/**
 * @file
 * @ingroup driver
 * @brief simc::driver::RunReport: run totals, normalisation and the text reports (.hist, .geni) of a run.
 *
 * End-of-run bookkeeping of the SIMC main program and its text reports (task P7.2).
 *   RunTotals      simc.f main-loop accumulators (l.160-352): ntried, nevent, ncontribute, npasscuts, wtcontribute,
 *                  sum_sigcc, sumerr/sumerr2 (sums_twoarm), contrib (contribtype) and the slop ranges
 *   Normalisation  simc.f l.360-432: normfac, genvol, the resolutions aveerr/resol; targetfac/luminosity (l.94-101)
 *   FormatHist()   simc.f:report (l.644-1142): the .hist text, line for line in the gfortran format
 *   FormatGeni()   simc.f l.462-532: the .geni stop-counter text
 *   MakeRunSummary()  the RunSummary (and its "extra" sections) handed to the output sinks (summary JSON)
 */
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

/// @addtogroup driver
/// @{

/// Fortran: type sums_twoarm (structures.inc) — per arm delta, xptar, yptar, ytar.
struct TwoArmSums {
    /// Sums for one arm: delta [%], xptar, yptar [rad], ytar [cm].
    struct Arm { double delta{0} /**< @copydoc ytar */, xptar{0} /**< @copydoc ytar */, yptar{0} /**< @copydoc ytar */, ytar{0}; /**< sums of delta, xptar, yptar, ytar */ };
    Arm e /**< @copydoc p */, p;  ///< e = electron arm, p = hadron arm.
};

/// Accumulators of the main loop (simc.f l.160-352 and COMMON /simulate/ counters).
struct RunTotals {
    long ntried{0};                      ///< tries
    long nevent{0};                      ///< successes (ngen > 0) or tries (ngen < 0)
    long ncontribute{0};  ///< successes with non-zero weight.
    long ncontribute_no_rad_proton{0};  ///< successes without a radiated proton.
    long npasscuts{0};  ///< successes that pass the SPedge cuts (PassCuts).
    double wtcontribute{0};              ///< sum of main%weight over pass_cuts successes (before *normfac)
    double sum_sigcc{0};                 ///< sum of main%sigcc over successes
    TwoArmSums sumerr /**< @copydoc sumerr2 */, sumerr2;          ///< pass_cuts successes
    reactions::Contrib contrib;          ///< "found" ranges (limits_update)
    Slops slop;                          ///< slop%MC / slop%total ranges (limits_update; .used from the setup)
};

/// P15.4 (parallel event loop): add the accumulators of one worker to the run's (counters and sums added, the
/// contrib/slop ranges combined with mergeContrib/mergeSlops). The workers are merged in the fixed order 0..N-1.
/// @param into the run's accumulators (updated).
/// @param worker the accumulators of one worker.
void MergeTotals(RunTotals& into, const RunTotals& worker);

/// simc.f l.94-101 and l.360-432.
struct Normalisation {
    double targetfac{0} /**< @copydoc luminosity */, luminosity{0};  ///< from the setup (simc.f:94-101)
    double charge{0};                    ///< EXPER%charge [mC]
    double genvol{0} /**< @copydoc normfac */, normfac{0};  ///< generation volume and normalisation factor (simc.f:360-398).
    double wtcontribute{0};              ///< wtcontribute*normfac (simc.f:398)
    TwoArmSums aveerr /**< @copydoc resol */, resol;            ///< zero unless npasscuts > 1 (simc.f:405)
};

/// Everything simc.f:report prints that is not in the RunSetup.
struct RunReport {
    std::string run_name;  ///< run name (output file stem).
    std::string begin_time /**< @copydoc end_time */, end_time;    ///< ctime() strings (simc.f:162, 357)
    long ngen{0};  ///< requested number of events (Fortran ngen).
    int random_seed{0};  ///< seed of the run.
    std::string rng_type;   ///< IRandom::name() of the run (D-054, P15.1)
    std::string random_state_file;  ///< file holding the random-generator state ("" if none).
    CentralEvent central;  ///< the central event.
    RunTotals totals;  ///< the accumulators of the main loop.
    Normalisation norm;  ///< normalisation and resolutions.
    std::vector<StopCounters> stops;     ///< one entry per spectrometer type in use
    int threads{1};                      ///< P15.4: worker threads of the event loop (1 = the sequential loop)
};

/// Fortran: simc.f l.360-432 (normalisation and resolutions). `genvol` from IReaction::genVolume().
/// @param setup the run set-up.
/// @param t the run totals.
/// @param genvol the generation volume.
/// @return the normalisation and resolutions.
Normalisation Normalise(const RunSetup& setup, const RunTotals& t, double genvol);

/// Fortran: simc.f:report(iun, ...) (l.644-1142). Text identical to the .hist the Fortran writes.
/// Throws std::runtime_error where the Fortran report STOPs (unknown reaction).
/// @param setup the run set-up.
/// @param r the run report.
/// @return the .hist file text.
/// @throws std::runtime_error for an unknown reaction.
std::string FormatHist(const RunSetup& setup, const RunReport& r);

/// Fortran: simc.f l.462-532 (.geni). The HMS, SOS, HRS-R/L and SHMS blocks are printed in the Fortran format (calo arms are
/// not supported by the driver yet); any other ISpectrometer plug-in gets a generic "key value" block.
/// @param setup the run set-up.
/// @param r the run report.
/// @return the .geni file text.
std::string FormatGeni(const RunSetup& setup, const RunReport& r);

/// The RunSummary for the output sinks, with the optional summary sections (central, target, radiative,
/// recon_summary, limits, slop, kinematics) in the summary.json layout (docs/OUTPUT_FORMAT.md §2).
/// @param setup the run set-up.
/// @param r the run report.
/// @return the run summary.
RunSummary MakeRunSummary(const RunSetup& setup, const RunReport& r);

/// ctime()-style time string as simc.f writes it (gfortran ctime: "Thu Oct  1 01:23:05 2026").
/// @param epoch_seconds time in seconds since the Unix epoch (local time zone).
/// @return the 24-character ctime string.
std::string FortranCtime(long epoch_seconds);

/// @}

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_RUNREPORT_H
