// include/simc/core/Fwd.h
// Forward declarations of the framework data model (docs/ARCHITECTURE.md §2).
// The full definitions (mirroring structures.inc) arrive with the driver tasks;
// stage interfaces only use them by reference.
#ifndef SIMC_CORE_FWD_H
#define SIMC_CORE_FWD_H

namespace simc {

struct SimcParameters;  // core/SimcParameters.h — every dbase.f:regallvars variable
struct RunSetup;        // derived setup (init.f / dbase.f equivalents)
struct Event;           // Fortran: type(event)       (vertex, orig, recon)
struct MainEvent;       // Fortran: type(event_main)
struct RunSummary;      // end-of-run numbers (ntried, normfac, stop counters ...)
class IRandom;

}  // namespace simc

#endif  // SIMC_CORE_FWD_H
