/**
 * @file
 * @ingroup driver
 * @brief The "central" event: the event at the central spectrometer settings (simc.f:calculate_central).
 *
 * The "central" event (task P7.2, from P2.3). Fortran: simc.f:calculate_central(central, vertex0) (l.1143-1309),
 * type(event_central) in structures.inc.
 *
 * calculateCentral() builds vertex0 from the central spectrometer settings exactly like the Fortran (zero
 * deltas/angles, Ein = Ebeam_vertex_ave, /pfermi/ = 0), runs complete_recon_ev on it (IReaction::completeRecon),
 * then radc_init_ev (IRadiator::initEvent) with main0's teff(1:3) = the target averages, and complete_main with
 * force_sigcc = .true. The radiation "Central Values" printed in the .hist (central%rad) are the radiator's
 * currentValues() right after that radc_init_ev call (only when using_rad, simc.f:1268).
 *
 * FORTRAN-QUIRKs kept: main0 is a SAVEd local, so the fields calculate_central does not set are 0 (in particular
 * main0%jacobian = 0, so main0%weight = 0; only main0%sigcc is used); complete_recon_ev is used instead of
 * complete_ev ("JRA 2002"), so for H(e,e'p) the central Em/Pm/W/MM are those of the *reconstructed* central event
 * (non-zero Em and Pm, negative MM in the reference .hist).
 */
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

/// @addtogroup driver
/// @{

/// Fortran: central%rad (radc.inc values copied by calculate_central when using_rad): the radiator's
/// RadiationValues (tails 1..3 at [0..2]).
using CentralRad = RadiationValues;

/// Fortran: type event_central (structures.inc).
struct CentralEvent {
    /// Central kinematics: Q2 [MeV^2], q, nu, Em, Pm, W, MM [MeV] (same units as Event).
    double Q2{0}, q{0} /**< @copydoc Q2 */, nu{0} /**< @copydoc Q2 */, Em{0} /**< @copydoc Q2 */, Pm{0} /**< @copydoc Q2 */, W{0} /**< @copydoc Q2 */, MM{0} /**< @copydoc Q2 */;
    /// Central reconstructed quantities of one arm (delta [%], xptar/yptar [rad]); `e` = electron arm, `p` = hadron arm.
    struct Arm { double delta{0} /**< @copydoc yptar */, xptar{0} /**< @copydoc yptar */, yptar{0}; /**< delta [%], xptar, yptar [rad] */ } e /**< @copydoc p */, p;  ///< e = electron arm, p = hadron arm
    double sigcc{0};  ///< central cross section (sigcc) from complete_main with force_sigcc = .true.
    CentralRad rad;  ///< radiation central values (filled only when radiation is on).
};

/// Result of calculate_central: the central event and vertex0 (passed to complete_main for every event).
struct CentralResult {
    CentralEvent central;  ///< the central event (printed in the .hist).
    Event vertex0;  ///< the central vertex, passed to complete_main for every event.
    MainEvent main0;   ///< the SAVEd local main0 after the call (diagnostics)
};

/// Fortran: simc.f:calculate_central (l.1143-1309).
/// Throws std::runtime_error with the Fortran STOP text if complete_recon_ev fails.
/// complete_main is the reaction's (IReaction::completeMain, with the run's physics model); `model` only receives
/// beginEvent(0) (nevent = 0 before the loop).
/// @param setup the run set-up (central spectrometer settings, target averages).
/// @param reaction the run's reaction (completeRecon / completeMain).
/// @param radiator the run's radiator (initEvent, currentValues).
/// @param model the run's physics model (receives beginEvent(0) only).
/// @return the central event, vertex0 and main0.
/// @throws std::runtime_error if complete_recon_ev fails.
CentralResult calculateCentral(const RunSetup& setup, IReaction& reaction, IRadiator& radiator, IPhysicsModel& model);

/// @}
}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_CENTRAL_H
