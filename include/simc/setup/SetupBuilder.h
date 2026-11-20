/**
 * @file
 * @ingroup input
 * @brief simc::BuildRunSetup: SimcParameters -> RunSetup, the Fortran start-up sequence (simc.f l.53-101).
 *
 * BuildRunSetup: SimcParameters -> RunSetup, reproducing the Fortran start-up sequence of simc.f
 * (l.53-101) minus calculate_central:
 *   init.f:min_max_init (l.921-1193)        -> the ±1e10 defaults of the RunSetup structs
 *   dbase.f:dbase_read post-read (l.119-966) -> flags, masses, unit conversions, cuts, switches
 *     init.f:target_init (l.1-88)             (through ITargetLossProvider, default HallCCryoTarget)
 *     init.f:theory_init (l.828-920)          (D(e,e'p), A(e,e'p) without Benhar SF)
 *     dbase.f:571-596 momentum distributions  (deut/he pion, kaon, semi)
 *     init.f:limits_init (l.91-575)           (extreme_trip_thru_target through the provider)
 *     dbase.f:688-958 announcement STOPs
 *   init.f:radc_init (l.576-654)
 *   simc.f:94-101 targetfac, luminosity
 *     dbase.f:592-621 spectral-function choice (the Benhar table itself is read by the "benhar" ISpectralFunction)
 * Not ported here: sf_lookup_init for A>2 pion/kaon (Emval(numEm) limit),
 * the kaon model tables (they belong to the kaon cross-section model) and calculate_central (P2.3).
 */
#ifndef SIMC_SETUP_SETUPBUILDER_H
#define SIMC_SETUP_SETUPBUILDER_H

#include <string>

#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/stages/ITargetLossProvider.h"

namespace simc {

/// @addtogroup input
/// @{

/// Options of BuildRunSetup().
struct SetupOptions {
    std::string run_name{"simc"};  ///< run name (output file stem).
    /// Directory with the Fortran data files read at start-up (h2.theory, c12.theory, deut.dat,
    /// he3.dat, ...). The Fortran opens them in its working directory. Empty = current directory.
    std::string data_dir;
    /// D-031: false = correct the clear Fortran bugs (default), true = exact Fortran behaviour (bitwise tests).
    /// Copied into RunSetup::fortran_compat and ExtremeTripArgs::fortran_compat.
    bool fortran_compat{false};
    /// Task U.2: ntuple layout (upstream 59ccdf6 = default, legacy 294f662). Copied into RunSetup::tree_layout.
    TreeLayout tree_layout{TreeLayout::Upstream};
    /// Task U.R: upstream commit reproduced by fortran_compat. Copied into RunSetup::compat_target.
    CompatTarget compat_target{CompatTarget::U59ccdf6};
    /// D(e,e'p)/A(e,e'p): the spectral-function file named explicitly (JSON pipeline.reaction.spectral_function.file;
    /// task P8.1, D-033 F1). Empty = the Fortran choice from nint(A) (init.f:839-851, dbase.f:598-616), with a warning
    /// when it falls back to the C-12 file.
    std::string sf_file;
};

/// Build the run setup. `losses` = the target model's start-up losses (any ITarget; driver::SimcRun passes a dedicated
/// instance of the pipeline target); null uses material::HallCCryoTarget, see stages/ITargetLossProvider.h. Fortran STOP conditions are
/// returned in RunSetup::stops; the builder returns at the point where the Fortran would stop.
/// @param p the run parameters.
/// @param opt build options.
/// @param losses target start-up losses (not owned), or nullptr for the Hall C cryotarget.
/// @return the run set-up (check RunSetup::stops).
RunSetup BuildRunSetup(const SimcParameters& p, const SetupOptions& opt = {},
                       ITargetLossProvider* losses = nullptr);

/// Fortran: init.f:theory_init (l.828-920) reading `path`. `doing_heavy` enables the Em_int part.
/// Returns false on a read error or unequal Pm bins (the Fortran then STOPs 'THEORY_INIT failed!').
/// @param path the theory file.
/// @param doing_heavy true for A(e,e'p) (reads the Em_int part).
/// @param[out] out the theory tables.
/// @param[out] msgs messages (optional, may be nullptr).
/// @return false on a read error.
bool ReadTheoryFile(const std::string& path, bool doing_heavy, TheoryData& out, std::vector<std::string>* msgs = nullptr);

/// Fortran: dbase.f:584-593 (read up to 2000 (pval, mprob) pairs, normalise to the last mprob).
/// @param path the momentum distribution file.
/// @param[out] out the distribution.
/// @return false on a read error.
bool ReadMomentumDistribution(const std::string& path, MomentumDistribution& out);

/// @}

}  // namespace simc

#endif  // SIMC_SETUP_SETUPBUILDER_H
