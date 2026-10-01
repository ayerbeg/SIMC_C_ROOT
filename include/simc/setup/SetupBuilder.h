// include/simc/setup/SetupBuilder.h
//
// BuildRunSetup: SimcParameters -> RunSetup, reproducing the Fortran start-up sequence of simc.f
// (l.53-101) minus calculate_central:
//   init.f:min_max_init (l.921-1193)        -> the ±1e10 defaults of the RunSetup structs
//   dbase.f:dbase_read post-read (l.119-966) -> flags, masses, unit conversions, cuts, switches
//     init.f:target_init (l.1-88)             (through ITargetLossProvider, default HallCCryoTarget)
//     init.f:theory_init (l.828-920)          (D(e,e'p), A(e,e'p) without Benhar SF)
//     dbase.f:571-596 momentum distributions  (deut/he pion, kaon, semi)
//     init.f:limits_init (l.91-575)           (extreme_trip_thru_target through the provider)
//     dbase.f:688-958 announcement STOPs
//   init.f:radc_init (l.576-654)
//   simc.f:94-101 targetfac, luminosity
// Not ported here: sf_lookup_init (Benhar spectral function, needed for A>2 pion/kaon limits),
// the kaon model tables (they belong to the kaon cross-section model) and calculate_central (P2.3).
#ifndef SIMC_SETUP_SETUPBUILDER_H
#define SIMC_SETUP_SETUPBUILDER_H

#include <string>

#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/stages/ITargetLossProvider.h"

namespace simc {

struct SetupOptions {
    std::string run_name{"simc"};
    /// Directory with the Fortran data files read at start-up (h2.theory, c12.theory, deut.dat,
    /// he3.dat, ...). The Fortran opens them in its working directory. Empty = current directory.
    std::string data_dir;
};

/// Build the run setup. `losses` = the target model's start-up losses (any ITarget; driver::SimcRun passes a dedicated
/// instance of the pipeline target); null uses material::HallCCryoTarget, see stages/ITargetLossProvider.h. Fortran STOP conditions are
/// returned in RunSetup::stops; the builder returns at the point where the Fortran would stop.
RunSetup BuildRunSetup(const SimcParameters& p, const SetupOptions& opt = {},
                       ITargetLossProvider* losses = nullptr);

/// Fortran: init.f:theory_init (l.828-920) reading `path`. `doing_heavy` enables the Em_int part.
/// Returns false on a read error or unequal Pm bins (the Fortran then STOPs 'THEORY_INIT failed!').
bool ReadTheoryFile(const std::string& path, bool doing_heavy, TheoryData& out, std::vector<std::string>* msgs = nullptr);

/// Fortran: dbase.f:584-593 (read up to 2000 (pval, mprob) pairs, normalise to the last mprob).
bool ReadMomentumDistribution(const std::string& path, MomentumDistribution& out);

}  // namespace simc

#endif  // SIMC_SETUP_SETUPBUILDER_H
