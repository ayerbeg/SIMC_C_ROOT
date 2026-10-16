// include/simc/reactions/DeltaProduction.h
//
// IReaction "delta_production" (task P11.1): (e,e'p)pi in the Delta region, the doing_delta branches of event.f. The
// PROTON is detected (Mh = Mp, dbase.f:182-183) and the pion is the missing 'recoil' (targ%Mrec_struck = Mpi = 139.57018,
// the charged-pion mass, dbase.f:297-300; FORTRAN-DIFF-CANDIDATE D-048 F7: for gamma* p -> p X the missing pion is a pi0).
// The Fortran shares generate/complete_ev/complete_recon_ev with the pion ("doing_delta: as doing_pion", event.f:230), so
// this class is PionProduction with the delta flags:
//   generate()      event.f:126-428 (electron E and angles, proton angles; Fermi momentum for D/A: corrected mode only).
//   completeEv()    event.f:432-1052: the pion quadratic with the proton as the hadron; for doing_delta ONE grnd picks the
//                   second root when grnd() > 0.5 (l.682-684, FORTRAN-DIFF-CANDIDATE D-048 F5: no factor 2); with
//                   which_pion 2/3 the Breit-Wigner Delta mass replaces targ%Mrec_struck first (l.612-623, D-048 F8: below
//                   threshold in the Delta region, every such call fails).
//   completeMain()  event.f:1363-1569, delta branch (l.1511-1513): SF_weight = 1, sigcc = peedelta(vertex, main) (the run's
//                   IPhysicsModel, "delta"), sigcc_recon = 1, tgtweight = 1 (not set on this branch), Coulomb focusing
//                   factor, weight. FORTRAN-DIFF-CANDIDATE(D-048 F4): peedelta contains NO cross section (placeholder 1.0).
//   prepareCentral() FORTRAN-BUG(D-048 F9 = D-038 F5): the central peedelta sees main0%epsilon = t = phi_pq = 0; corrected:
//                   main0 gets the central kinematics (inherited from PionProduction).
// D and A targets (doing_deutdelta, doing_hedelta) STOP in the Fortran (dbase.f:779-788, D-018 / D-048 F1); in corrected
// mode they run through the D/A branches of event.f, which no one has ever run: validated by self-consistency only.
// Random numbers: as the pion (9 grnd in generate for H) + 1 grnd for the root choice (when complete_ev reaches it).
#ifndef SIMC_REACTIONS_DELTAPRODUCTION_H
#define SIMC_REACTIONS_DELTAPRODUCTION_H

#include <string>

#include "simc/reactions/PionProduction.h"

namespace simc {
namespace reactions {

class DeltaProduction : public PionProduction {
public:
    explicit DeltaProduction(std::string name = "delta_production");
    bool completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0,
                      const Event& recon) override;
};

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_DELTAPRODUCTION_H
