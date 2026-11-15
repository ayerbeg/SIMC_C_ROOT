/**
 * @file
 * @ingroup stages
 * @brief simc::StageContext: what every stage receives at init() (explicit dependency injection).
 *
 * The driver (simc::driver::SimcRun) creates the stages through the registry, initialises the target and the radiator
 * first and then hands them, through their INTERFACES, to the stages that need them (beam, reaction, ...). A stage
 * that needs a collaborator keeps the pointer it gets here; it must not cast it to a concrete class. Task P7.4.
 */
#ifndef SIMC_STAGES_STAGECONTEXT_H
#define SIMC_STAGES_STAGECONTEXT_H

#include <stdexcept>
#include <string>

#include "simc/core/Fwd.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

class ITarget;
class IRadiator;
class IPhysicsModel;

/**
 * @brief The run set-up and the collaborators handed to simc::IStage::init().
 *
 * A lightweight view: it holds a reference and non-owning pointers, all owned by the driver and valid for the whole
 * run. The collaborators are null when they do not exist yet (the target while the target itself is initialised, the
 * radiator for the target and the radiator, the physics model for the stages created before it).
 */
struct StageContext {
    /**
     * @brief A context with the set-up only (no collaborators).
     *
     * Implicit on purpose: `stage.init(setup)` initialises a stage that has no collaborators (stand-alone use, unit
     * tests).
     * @param s The run set-up; must outlive every stage initialised with this context.
     */
    StageContext(const RunSetup& s) : setup(s) {}   // NOLINT(google-explicit-constructor): see the header comment

    const RunSetup& setup;               ///< The run set-up (read-only).
    ITarget* target{nullptr};            ///< The run's target, or null (non-owning).
    IRadiator* radiator{nullptr};        ///< The run's radiator, or null (non-owning).
    IPhysicsModel* physics{nullptr};     ///< The run's physics model, or null (non-owning; the reaction weights the event with it).

    /**
     * @brief The target, which the calling stage requires.
     * @param who Name of the calling stage (for the error message).
     * @return The target.
     * @throws std::logic_error if #target is null.
     */
    ITarget& requireTarget(const std::string& who) const {
        if (!target) throw std::logic_error(who + ": no target in the StageContext");
        return *target;
    }
    /**
     * @brief The radiator, which the calling stage requires.
     * @param who Name of the calling stage (for the error message).
     * @return The radiator.
     * @throws std::logic_error if #radiator is null.
     */
    IRadiator& requireRadiator(const std::string& who) const {
        if (!radiator) throw std::logic_error(who + ": no radiator in the StageContext");
        return *radiator;
    }
    /**
     * @brief The physics model, which the calling stage requires.
     * @param who Name of the calling stage (for the error message).
     * @return The physics model.
     * @throws std::logic_error if #physics is null.
     */
    IPhysicsModel& requirePhysics(const std::string& who) const {
        if (!physics) throw std::logic_error(who + ": no physics model in the StageContext");
        return *physics;
    }
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_STAGECONTEXT_H
