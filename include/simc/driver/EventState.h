/**
 * @file
 * @ingroup driver
 * @brief simc::driver::EventState: the per-event state of the SIMC main loop, shared by all tries of a run.
 *
 * The per-event state of the SIMC main loop (task P7.2), passed explicitly through the stages.
 *
 * Fortran: simc.f main program locals `main`, `vertex`, `orig`, `recon` (type(event_main) / type(event)), which are
 * static under -fno-automatic and are NOT cleared between tries. Several Fortran routines only write parts of
 * them, so a try can see values left by an EARLIER try. SIMC++ keeps exactly this behaviour by holding ONE
 * EventState for the whole run and handing the same object to every stage on every try. The persistence is
 * part of the physics output (D-024), not an accident of the implementation:
 *   * orig        is not written when generate_rad fails (radc.f) -> keeps the previous try's values;
 *   * main.SP / main.RECON / main.FP and recon arm fields: montecarlo returns early on a failed arm
 *     (simc.f:1310-1853), so the remaining fields keep earlier values (they feed pass_cuts and the geni
 *     histograms of failed tries, and limits_update of later successes only through fields that ARE rewritten);
 *   * vertex.p.xptar/yptar are read by generate before complete_ev rewrites them (event.f:325, harmless);
 *   * main.ntup: the loop resets radphot, radarm, resfac, sigcm, krel, mm, mmA, t, decdist, Mh2_final at the
 *     start of every try (simc.f:174-183) — done by SimcRun::resetForTry().
 * Stage-internal Fortran state (COMMON blocks / SAVEd locals of the radiation code, the spectrometer counters,
 * enerloss_new statics, the montecarlo tmpfact ...) lives inside the stage objects, not here.
 */
#ifndef SIMC_DRIVER_EVENTSTATE_H
#define SIMC_DRIVER_EVENTSTATE_H

#include "simc/core/Event.h"

namespace simc {
namespace driver {

/// @addtogroup driver
/// @{

/// The per-event state of the main loop (see the file description: values persist between tries on purpose).
struct EventState {
    MainEvent main;     ///< Fortran: main   (type(event_main))
    Event vertex;       ///< Fortran: vertex (generated physics vertex, after generate_rad)
    Event orig;         ///< Fortran: orig   (after radiation, before the spectrometers)
    Event recon;        ///< Fortran: recon  (reconstructed event)
    bool success{false};     ///< this try contributes (after montecarlo, complete_recon_ev, complete_main, hard_cuts)
    bool pass_cuts{false};   ///< simc.f:225-237 SPedge test on recon (computed for every try)
    long ntried{0};          ///< 1-based number of this try
};

/// @}

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_EVENTSTATE_H
