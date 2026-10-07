// include/simc/stages/IReaction.h
// Reaction: generates the physics vertex, completes the kinematics and assembles the event weight.
// Fortran: event.f:generate (l.126-429), complete_ev (l.432-1054), complete_recon_ev (l.1056-1360),
//          complete_main (l.1363-1569); simc.f genvol (l.383-396).
#ifndef SIMC_STAGES_IREACTION_H
#define SIMC_STAGES_IREACTION_H

#include "simc/stages/IStage.h"

namespace simc {

class IReaction : public IStage {
public:
    /// Generate the variables of this reaction, complete the vertex and radiate it (event.f:generate after the
    /// beam part, l.209-428: generated variables, complete_ev, generate_rad or orig = vertex). `orig` is written
    /// only where the Fortran writes it (it keeps the previous try's values when generate_rad fails, D-024).
    /// The beam part (IBeam::generate) must have been called first for this try.
    virtual bool generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) = 0;
    /// Complete the vertex kinematics from the generated variables (complete_ev, event.f:432-1054).
    /// Also the CompleteEvFn handed to IRadiator::radiate.
    virtual bool completeEv(MainEvent& main, Event& vertex, IRandom& rng) = 0;
    /// Complete the reconstructed event (complete_recon_ev).
    virtual bool completeRecon(Event& recon) = 0;
    /// The event weight (event.f:complete_main, task P8.2, D-030): main%SF_weight (the reaction's spectral function),
    /// main%sigcc / sigcc_recon (the run's IPhysicsModel, StageContext::physics), the Coulomb focusing factor and
    /// main%weight = SF_weight*jacobian*gen_weight*sigcc*(nucleon count). Returns false (main%weight not written) when
    /// the event cannot contribute (SF_weight <= 0) unless `force_sigcc` (the central event).
    virtual bool completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0,
                              const Event& recon) = 0;
    /// Generated phase-space volume factor entering normfac (simc.f:383-396).
    virtual double genVolume() const = 0;
    /// Called by calculate_central (simc.f:1143-1300) after complete_recon_ev(vertex0) and before
    /// complete_main(.true., main0, ...): lets a reaction set the central main0 quantities its cross section reads
    /// (task P9.2: the pion reaction fills main0%epsilon, t, phi_pq in corrected mode, D-038 F5). Default: nothing.
    virtual void prepareCentral(MainEvent& main0, const Event& vertex0) {
        (void)main0;
        (void)vertex0;
    }
    /// The Fortran try counter `ntried` of the coming try (simc.f:190), before generate(). Only used in messages
    /// (event.f:591, 'two solutions' warning for ntried <= 5000). Default: ignored.
    virtual void beginTry(long ntried) { (void)ntried; }
};

}  // namespace simc

#endif  // SIMC_STAGES_IREACTION_H
