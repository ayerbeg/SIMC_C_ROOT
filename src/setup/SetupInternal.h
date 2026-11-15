// src/setup/SetupInternal.h — private to src/setup: the phases of BuildRunSetup (SetupBuilder.cpp), split by
// responsibility without changing the operation order of the Fortran start-up. Each phase returns false where the original single function
// returned the RunSetup early (a Fortran STOP recorded in RunSetup::stops).
//   SetupReaction.cpp : dbase.f:124-547   reaction flags, Mh, switches, kinematics_main, target masses, struck/recoil
//                                          masses, target geometry, offsets, SPedge, radiation tails, gen angles, flags
//   SetupTheory.cpp   : init.f:target_init, theory_init, the momentum distribution and spectral-function choice
//                       (dbase.f:551-681) and the data-file readers (ReadTheoryFile, ReadBenharEmMax, ReadMomentumDistribution)
//   SetupLimits.cpp   : init.f:limits_init (l.91-575) incl. extreme_trip_thru_target and the histogram axes
//   SetupBuilder.cpp  : BuildRunSetup (the sequence), dbase.f:689-958 reaction STOP checks, init.f:radc_init, simc.f:94-101
#ifndef SIMC_SRC_SETUP_SETUPINTERNAL_H
#define SIMC_SRC_SETUP_SETUPINTERNAL_H

#include <cmath>
#include <string>

#include "simc/setup/SetupBuilder.h"

namespace simc {
namespace setup_detail {

inline int nint(double x) { return static_cast<int>(std::lround(x)); }   // Fortran nint: half away from zero

inline std::string joinPath(const std::string& dir, const std::string& file) {
    if (dir.empty()) return file;
    return dir.back() == '/' ? dir + file : dir + "/" + file;
}

/// dbase.f:124-547 (SetupReaction.cpp).
bool BuildReactionAndKinematics(RunSetup& s, const SimcParameters& p);
/// init.f:target_init through `losses`, theory_init, momentum distribution, spectral function (SetupTheory.cpp).
/// benharEmMax = Emval(numEm) of the Benhar table when it is read (A > 2 pion/kaon), else 0.
bool BuildTargetAndTheory(RunSetup& s, const SimcParameters& p, const SetupOptions& opt, ITargetLossProvider* losses,
                          double& benharEmMax);
/// init.f:limits_init and the diagnostic histogram axes (SetupLimits.cpp).
bool BuildLimits(RunSetup& s, const SimcParameters& p, ITargetLossProvider* losses, double benharEmMax);

}  // namespace setup_detail
}  // namespace simc

#endif  // SIMC_SRC_SETUP_SETUPINTERNAL_H
