// include/simc/stages/IPhysicsModel.h
// Cross-section model for a reaction.
// Fortran: physics_proton.f (sigep, deForest), physics_pion.f, physics_kaon.f, ...
//          called from event.f:complete_main (l.1363-1570).
#ifndef SIMC_STAGES_IPHYSICSMODEL_H
#define SIMC_STAGES_IPHYSICSMODEL_H

#include "simc/stages/IStage.h"

namespace simc {

class IPhysicsModel : public IStage {
public:
    /// main%sigcc for the vertex event.
    virtual double sigcc(const Event& vertex, const MainEvent& main) = 0;
    /// main%sigcc_recon for the reconstructed event (1.0 where the Fortran does not compute it).
    virtual double sigccRecon(const Event& recon, const MainEvent& main) = 0;
    /// Called by the driver before each event with the Fortran `nevent` counter (simc.f:165-349).
    /// Models only use it in WARNING messages (e.g. deForest, physics_proton.f:75). Default: ignore.
    virtual void beginEvent(long nevent) { (void)nevent; }
};

}  // namespace simc

#endif  // SIMC_STAGES_IPHYSICSMODEL_H
