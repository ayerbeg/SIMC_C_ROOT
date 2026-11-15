/**
 * @file
 * @ingroup core
 * @brief Forward declarations of the framework data model.
 *
 * The stage interfaces use these types only by reference, so they include this header instead of the full
 * definitions: simc::SimcParameters (core/SimcParameters.h), simc::RunSetup (core/RunSetup.h), simc::Event and
 * simc::MainEvent (core/Event.h), simc::RunSummary (core/RunSummary.h) and simc::IRandom (stages/IRandom.h).
 */
#ifndef SIMC_CORE_FWD_H
#define SIMC_CORE_FWD_H

namespace simc {

struct SimcParameters;  // core/SimcParameters.h — every input variable of a .inp file
struct RunSetup;        // core/RunSetup.h — derived run constants
struct Event;           // core/Event.h — one set of event kinematics (vertex, orig, recon)
struct MainEvent;       // core/Event.h — weights, target, spectrometer quantities of the try
struct RunSummary;      // core/RunSummary.h — end-of-run numbers (ntried, normfac, stop counters ...)
class IRandom;

}  // namespace simc

#endif  // SIMC_CORE_FWD_H
