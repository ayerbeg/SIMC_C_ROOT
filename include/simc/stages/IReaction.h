// include/simc/stages/IReaction.h
// Reaction: generates the physics vertex and completes the kinematics.
// Fortran: event.f:generate (l.126-429), complete_ev (l.432-1054),
//          complete_recon_ev (l.1056-1360); simc.f genvol (l.383-396).
#ifndef SIMC_STAGES_IREACTION_H
#define SIMC_STAGES_IREACTION_H

#include "simc/stages/IStage.h"

namespace simc {

class IReaction : public IStage {
public:
    /// Generate the variables of this reaction and complete the vertex (generate + complete_ev).
    virtual bool generate(MainEvent& main, Event& vertex, IRandom& rng) = 0;
    /// Complete the reconstructed event (complete_recon_ev).
    virtual bool completeRecon(Event& recon) = 0;
    /// Generated phase-space volume factor entering normfac (simc.f:383-396).
    virtual double genVolume() const = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_IREACTION_H
