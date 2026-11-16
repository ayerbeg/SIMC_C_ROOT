/**
 * @file
 * @ingroup reactions
 * @brief simc::reactions::KaonProduction: IReaction "kaon_production", kaon electroproduction e + A -> e' K+ Y.
 *
 * IReaction "kaon_production" (task P10.2): kaon electroproduction e + A -> e' K+ Y, the doing_kaon branches of event.f.
 * The Fortran shares generate/complete_ev/complete_recon_ev between pions and kaons ("doing_kaon: as doing_pion",
 * event.f:229), so this class reuses PionProduction with the kaon flags:
 *   which_kaon 0 / 1 / 2     K+ Lambda / K+ Sigma0 (from p) / K+ Sigma- (from n), quasifree: H (doing_hydkaon), D
 *                            (doing_deutkaon, Fermi momentum from deut.dat), A > 2 (doing_hekaon, Fermi momentum + Benhar Em)
 *   which_kaon 10 / 11 / 12  coherent production on A >= 3 (bound hypernucleus; treated as hydrogen with Mtar_struck = M_A)
 *   completeMain()  event.f:complete_main l.1363-1569, kaon branch: SF_weight = 1, sigcc = peeK(vertex, main,
 *                   survivalprob) (the run's IPhysicsModel, "kaon", through sigccMain), sigcc_recon = 1, tgtweight =
 *                   targ%N (which_kaon 2/12) or targ%Z (l.1502-1509; FORTRAN-DIFF-CANDIDATE D-043 F9 "OK for coherent???"),
 *                   Coulomb focusing factor, weight; when the decay is not simulated (doing_decay = 0) the weight is
 *                   multiplied by the survival probability (l.1560-1561), which peeK leaves in the COMMON decdist
 *                   (main.ntup.decdist, physics_kaon.f:164).
 *   prepareCentral() FORTRAN-BUG(D-043 F13) = D-038 F5: calculate_central evaluates the central cross section with
 *                   main0%epsilon = t = theta_pq = phi_pq = 0; corrected: the central kinematics (as for the pion).
 * Decay in flight (doing_decay = 1) happens in the hadron spectrometer (optics::CosyEngine, K -> mu nu / pi pi).
 * Random numbers: as PionProduction (no extra number on the kaon branch).
 */
#ifndef SIMC_REACTIONS_KAONPRODUCTION_H
#define SIMC_REACTIONS_KAONPRODUCTION_H

#include <string>

#include "simc/reactions/PionProduction.h"

namespace simc {
namespace reactions {

/// @addtogroup reactions
/// @{

/// IReaction for the doing_kaon branches of event.f (shares the pion event chain). Not thread-safe.
class KaonProduction : public PionProduction {
public:
    /// @param name registry name.
    explicit KaonProduction(std::string name = "kaon_production");
    bool completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0,
                      const Event& recon) override;
    /// The SAVEd local survivalprob of complete_main (event.f:1372) after the last call.
    /// @return the kaon decay survival probability.
    double survivalProbability() const { return survivalprob_; }

private:
    double survivalprob_{0};
};

/// @}

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_KAONPRODUCTION_H
