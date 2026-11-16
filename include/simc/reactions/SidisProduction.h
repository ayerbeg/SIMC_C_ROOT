/**
 * @file
 * @ingroup reactions
 * @brief simc::reactions::SidisProduction: IReaction "sidis_production", semi-inclusive (e,e'h)X.
 *
 * IReaction "sidis_production" (task P12.2): semi-inclusive (e,e'h)X, the doing_semi branches of event.f (doing_semipi:
 * pi+/pi-, doing_semika: K+/K-; H, D (doing_deutsemi) and A > 2 (doing_hesemi)). The cross section is the run's
 * IPhysicsModel ("sidis" = semi_physics.f:peepiX, task P12.1).
 *   generate()       event.f:126-428, semi branches: electron and hadron angles (4 grnd), the HADRON energy first in the
 *                    sumEgen window (l.289-298, 1 grnd; gen_weight), the electron energy in its own generation limits
 *                    (l.301-317, no sumEgen adjustment for semi; 1 grnd), physics_angles, the deuteron Fermi momentum
 *                    (l.336-361, 4 grnd, drawn for doing_deutsemi even when do_fermi = 0; no Fermi momentum for A > 2),
 *                    complete_ev, generate_rad or orig = vertex.
 *   completeEv()     event.f:432-1052, semi branch: no hadron-momentum solution (the generated hadron is kept), W,
 *                    epsilon, theta_pq, t, tmin, phi_pq (l.723-786; pizero_decay for doing_pizero), the Pm vector, Pm = Pmiss
 *                    and Em = Emiss (l.952-954), the two-pion threshold (l.979-988), zhad, pt2 and the z > 1 test
 *                    (l.990-996), both 1/cos**3 jacobians, trip_thru_target 2/3 (8 gauss1), radc_init_ev.
 *   completeRecon()  EepReactionBase (complete_recon_ev: W from Mp, Em/mm/mmA/t, zhad, pt2).
 *   completeMain()   event.f:1363-1569, semi branch: SF_weight = 1, sigcc = peepiX (sigccMain), sigcc_recon = 1, ntup%dilu
 *                    = 1, tgtweight = 1, Coulomb factor, weight; * survivalprob for doing_semika without decay simulation
 *                    (l.1560-1561).
 *   prepareCentral() nothing: calculate_central's peepiX call uses vertex0 (z, pt2 of complete_recon_ev), which is
 *                    meaningful as it is (D-051 F13); main0 epsilon/t are not read by peepiX.
 * Random numbers (H): 10 grnd + 4 gauss1 in generate/beam, + 8 gauss1 when complete_ev passes the semi tests; D: + 4 grnd.
 *
 * fortran_compat (D-031; docs/FORTRAN_DIFFERENCES.md):
 *   FORTRAN-BUG(D-051 F9, as D-033 F5) event.f:922: PmPer from the previous try's vertex%Pm; corrected: from this Pm.
 *   FORTRAN-DIFF-CANDIDATE(D-051 F8) event.f:979-996: the two-pion threshold (Mp+Mpi0)^2 for every hadron, Mtar_struck
 *                    and no Fermi motion; the z > 1 test is unreachable. Kept.
 *   FORTRAN-DIFF-CANDIDATE(D-051 F16) init.f:504-506: hadron generation up to sumEgen%max - edge%e%E%min (delta up to
 *                    +116 %). Kept (set-up).
 *   FORTRAN-QUIRK main%sigcent is never assigned (ntuple sigcent, centjac = 0).
 */
#ifndef SIMC_REACTIONS_SIDISPRODUCTION_H
#define SIMC_REACTIONS_SIDISPRODUCTION_H

#include <string>

#include "simc/reactions/PionProduction.h"

namespace simc {
namespace reactions {

/// @addtogroup reactions
/// @{

/// IReaction for the doing_semi branches of event.f (pions or kaons; H, D or A). Not thread-safe.
class SidisProduction : public PionProduction {
public:
    /// @param name registry name.
    explicit SidisProduction(std::string name = "sidis_production");
    bool generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) override;
    bool completeEv(MainEvent& main, Event& vertex, IRandom& rng) override;
    bool completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0,
                      const Event& recon) override;
    void prepareCentral(MainEvent&, const Event&) override {}
    void finalTargetState(TargetInfo&) const override {}

    /// The SAVEd local survivalprob of complete_main (event.f:1372) as last used (test hook).
    /// @return the hadron decay survival probability.
    double survivalProbability() const { return survivalprob_; }

private:
    double survivalprob_{0};
};

/// @}

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_SIDISPRODUCTION_H
