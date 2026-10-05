// include/simc/driver/Central.h
//
// The "central" event (task P7.2, from P2.3). Fortran: simc.f:calculate_central(central, vertex0) (l.1143-1309),
// type(event_central) in structures.inc.
//
// calculateCentral() builds vertex0 from the central spectrometer settings exactly like the Fortran (zero
// deltas/angles, Ein = Ebeam_vertex_ave, /pfermi/ = 0), runs complete_recon_ev on it (IReaction::completeRecon),
// then radc_init_ev (IRadiator::initEvent) with main0's teff(1:3) = the target averages, and complete_main with
// force_sigcc = .true. The radiation "Central Values" printed in the .hist (central%rad) are the radiator's
// currentValues() right after that radc_init_ev call (only when using_rad, simc.f:1268).
//
// FORTRAN-QUIRKs kept: main0 is a SAVEd local, so the fields calculate_central does not set are 0 (in particular
// main0%jacobian = 0, so main0%weight = 0; only main0%sigcc is used); complete_recon_ev is used instead of
// complete_ev ("JRA 2002"), so for H(e,e'p) the central Em/Pm/W/MM are those of the *reconstructed* central event
// (non-zero Em and Pm, negative MM in the reference .hist).
#ifndef SIMC_DRIVER_CENTRAL_H
#define SIMC_DRIVER_CENTRAL_H

#include <array>

#include "simc/core/Event.h"
#include "simc/core/RunSetup.h"
#include "simc/reactions/EventChain.h"
#include "simc/stages/IPhysicsModel.h"
#include "simc/stages/IRadiator.h"
#include "simc/stages/IReaction.h"

namespace simc {
namespace driver {

/// Fortran: central%rad (radc.inc values copied by calculate_central when using_rad): the radiator's
/// RadiationValues (tails 1..3 at [0..2]).
using CentralRad = RadiationValues;

/// Fortran: type event_central (structures.inc).
struct CentralEvent {
    double Q2{0}, q{0}, nu{0}, Em{0}, Pm{0}, W{0}, MM{0};
    struct Arm { double delta{0}, xptar{0}, yptar{0}; } e, p;
    double sigcc{0};
    CentralRad rad;
};

/// Result of calculate_central: the central event and vertex0 (passed to complete_main for every event).
struct CentralResult {
    CentralEvent central;
    Event vertex0;
    MainEvent main0;   ///< the SAVEd local main0 after the call (diagnostics)
};

/// Fortran: simc.f:calculate_central (l.1143-1309).
/// Throws std::runtime_error with the Fortran STOP text if complete_recon_ev fails.
/// complete_main is the reaction's (IReaction::completeMain, with the run's physics model); `model` only receives
/// beginEvent(0) (nevent = 0 before the loop).
CentralResult calculateCentral(const RunSetup& setup, IReaction& reaction, IRadiator& radiator, IPhysicsModel& model);

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_CENTRAL_H
