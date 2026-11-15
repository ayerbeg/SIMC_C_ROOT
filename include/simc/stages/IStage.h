/**
 * @file
 * @ingroup stages
 * @brief simc::IStage: the common base of every pluggable stage, and its life cycle.
 *
 * Life cycle of a stage: create it by name (simc::Registry) -> configure(json, params) -> init(StageContext) -> the
 * per-event calls of its interface. See the developer manual, chapter 3, for the order of the calls in one try.
 */
#ifndef SIMC_STAGES_ISTAGE_H
#define SIMC_STAGES_ISTAGE_H

#include <string>

#include "simc/core/Fwd.h"
#include "simc/core/Json.h"
#include "simc/stages/StageContext.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief Common base of every pluggable stage of a run (random source, target, beam, reaction, physics model,
 *        radiator, spectrometer, output sink, ...).
 *
 * A stage is created by the driver through simc::Registry (default constructor), then configure() is called once with
 * its JSON configuration, then init() once with its simc::StageContext, then the per-event methods of the derived
 * interface any number of times. A stage object belongs to one run (one worker in the parallel loop) and is used by
 * one thread at a time; stages are owned by the driver (std::unique_ptr) and destroyed at the end of the run.
 */
class IStage {
public:
    /// Virtual destructor: stages are owned and deleted through interface pointers.
    virtual ~IStage() = default;
    /**
     * @brief The registry name of the implementation.
     * @return The name the stage is registered under (e.g. "hms", "hydrogen_elastic").
     */
    virtual std::string name() const = 0;
    /**
     * @brief Reads the stage's own JSON configuration and the global input parameters. Called once, before init().
     * @param config The stage's subtree of the JSON `pipeline` section (an empty object if none was given).
     * @param params The input parameters of the run (`.inp` units).
     * @throws std::invalid_argument for an invalid configuration (implementation-specific).
     */
    virtual void configure(const Json& config, const SimcParameters& params) = 0;
    /**
     * @brief Builds the derived set-up data (tables, matrices, limits) and takes the collaborators the stage needs.
     *
     * Called once after configure(). Collaborators (target, radiator, physics model) are taken from @p ctx and kept as
     * interface pointers; they outlive the stage. `stage.init(setup)` works for stages without collaborators (the
     * simc::StageContext constructor from a simc::RunSetup is implicit).
     * @param ctx The run set-up and the collaborators; the referenced simc::RunSetup must outlive the stage.
     * @throws std::logic_error if a needed collaborator is missing; simc::Error for a set-up the stage cannot handle.
     */
    virtual void init(const StageContext& ctx) = 0;
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_ISTAGE_H
