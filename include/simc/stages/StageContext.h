// include/simc/stages/StageContext.h
//
// What every stage receives at init() (task P7.4, explicit dependency injection; docs/ARCHITECTURE.md §4).
//
//   * setup     the run set-up (dbase.f / init.f quantities), read-only;
//   * target    the run's ITarget  (null when the stage being initialised IS the target);
//   * radiator  the run's IRadiator (null for the target and the radiator themselves).
//
// The driver (driver::SimcRun) creates the stages through the Registry, initialises the target and the radiator first
// and then hands them, through their INTERFACES, to the stages that need them (beam, reaction, ...). A stage that needs
// a collaborator keeps the pointer it gets here; it must not cast it to a concrete class.
//
// The constructor from a RunSetup is implicit on purpose: `stage.init(setup)` initialises a stage that has no
// collaborators (stand-alone use, unit tests). requireTarget()/requireRadiator() throw std::logic_error with the name
// of the stage when a needed collaborator is missing.
#ifndef SIMC_STAGES_STAGECONTEXT_H
#define SIMC_STAGES_STAGECONTEXT_H

#include <stdexcept>
#include <string>

#include "simc/core/Fwd.h"

namespace simc {

class ITarget;
class IRadiator;

struct StageContext {
    StageContext(const RunSetup& s) : setup(s) {}   // NOLINT(google-explicit-constructor): see the header comment

    const RunSetup& setup;
    ITarget* target{nullptr};
    IRadiator* radiator{nullptr};

    ITarget& requireTarget(const std::string& who) const {
        if (!target) throw std::logic_error(who + ": no target in the StageContext");
        return *target;
    }
    IRadiator& requireRadiator(const std::string& who) const {
        if (!radiator) throw std::logic_error(who + ": no radiator in the StageContext");
        return *radiator;
    }
};

}  // namespace simc

#endif  // SIMC_STAGES_STAGECONTEXT_H
