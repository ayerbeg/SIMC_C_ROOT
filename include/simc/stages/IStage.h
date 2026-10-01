// include/simc/stages/IStage.h
// Common base of every pluggable stage (docs/ARCHITECTURE.md §3).
// Life cycle: create (Registry) -> configure(json, params) -> init(StageContext) -> per-event calls.
#ifndef SIMC_STAGES_ISTAGE_H
#define SIMC_STAGES_ISTAGE_H

#include <string>

#include "simc/core/Fwd.h"
#include "simc/core/Json.h"
#include "simc/stages/StageContext.h"

namespace simc {

class IStage {
public:
    virtual ~IStage() = default;
    /// Registry name of the implementation (e.g. "hms", "hydrogen_elastic").
    virtual std::string name() const = 0;
    /// Read the stage's own JSON subtree and the global parameters. Called once.
    virtual void configure(const Json& config, const SimcParameters& params) = 0;
    /// Build derived setup data (tables, matrices, limits) and take the collaborators the stage needs from `ctx`
    /// (StageContext.h). Called once after configure(). `stage.init(setup)` works for stages without collaborators.
    virtual void init(const StageContext& ctx) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_ISTAGE_H
