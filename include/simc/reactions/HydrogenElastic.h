// include/simc/reactions/HydrogenElastic.h
//
// IReaction "hydrogen_elastic" and "phase_space" (task P4.1; on EepReactionBase since P8.2): event.f for H(e,e'p) and
// doing_phsp.
//   generate()        event.f:generate l.214-428 (after the beam part): electron angles, physics_angles,
//                     /pfermi/ reset, complete_ev, main%Trec, generate_rad (IRadiator) or orig = vertex
//   spectralWeight()  event.f:complete_main l.1397-1398: SF_weight = 1 (no spectral function)
//   completeEv(), completeRecon(), completeMain(), genVolume(): EepReactionBase (the shared event.f branches)
//
// Both registry names create this class: like the Fortran, the code path is chosen by the RunSetup flags
// (doing_hyd_elast, doing_phsp, ...), not by the name. D-009: hydrogen + doing_phsp keeps doing_hyd_elast, so
// complete_ev takes the elastic branch (event.f:546 before 709), while complete_recon_ev overwrites the hadron
// momentum with spec%p%P (event.f:1132-1136) and complete_main sets sigcc = 1. D(e,e'p)/A(e,e'p) are the
// "deuterium_eep"/"nucleus_eep" plug-ins (NuclearEep.h); pion, kaon, delta, rho, semi are refused in init().
//
// FORTRAN-QUIRKs (D-024): physics_angles for the hadron arm uses the previous try's p%xptar/yptar (result
// overwritten); orig is not written when generate_rad fails (see also EepReactionBase.h).
#ifndef SIMC_REACTIONS_HYDROGENELASTIC_H
#define SIMC_REACTIONS_HYDROGENELASTIC_H

#include <string>

#include "simc/reactions/EepReactionBase.h"

namespace simc {
namespace reactions {

class HydrogenElastic : public EepReactionBase {
public:
    explicit HydrogenElastic(std::string name = "hydrogen_elastic");

    void init(const StageContext& ctx) override;
    bool generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) override;
    double spectralWeight(const Event& vertex) override;
};

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_HYDROGENELASTIC_H
