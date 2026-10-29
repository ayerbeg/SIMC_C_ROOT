// include/simc/reactions/RhoProduction.h
//
// IReaction "rho_production" (task P13.2): exclusive rho0 electroproduction, the doing_rho branches of event.f, with the
// rho kernels of task P13.1 (physics/rho: generate_rho, rho_decay) and the "rho" model (rho_physics.f:peerho).
// H (doing_hydrho), D (doing_deutrho) and 3He (doing_herho); A >= 4 is refused at input (D-057 F16).
//
//   generate()     event.f:126-428 after the beam part: electron angles (2 grnd), the electron energy in the sumEgen
//                  window (1 grnd, gen_weight), NO hadron angles/energy (l.280-298), physics_angles for both arms (the
//                  hadron on the stale vertex%p%xptar/yptar = the previous rho's lab theta/phi, D-057 F1), the Fermi block
//                  for D/3He (4 grnd; pfer = 0 because no momentum distribution is loaded, D-057 F4; 3He: generate_em with
//                  the empty table, 2 grnd, Em = 0), complete_ev, generate_rad or orig = vertex, and rho_decay on orig
//                  (event.f:420-422, D-057 F2).
//   completeEv()   event.f:432-1052, rho branch: ue, the q vector, generate_rho (BW mass, rho in 4pi in the CM, 3 grnd),
//                  W, epsilon, theta_pq, t and tmin with the BW Mh2, phi_pq, the Pm vector with the never-assigned
//                  vertex%Pm, Trec, krel, the ELECTRON jacobian only (l.1016-1022), trip_thru_target 2 and 3 (8 gauss1),
//                  radc_init_ev.
//   completeRecon() EepReactionBase (ntup%t with the COMMON Mh2 = Mpi2 after the decay).
//   completeMain() event.f:1363-1569, rho branch: SF_weight = 1, sigcc = peerho(vertex = the rho, main), sigcc_recon = 1,
//                  tgtweight = targ%Z + targ%N, the Coulomb factor, the weight.
//   prepareCentral() nothing: calculate_central's peerho sees main0%phi_pq = 0 (D-057 F17, documented, as the Fortran).
//   genVolume()    EepReactionBase (dOmega_e * 4 pi * dE_e, simc.f:375-391).
// COMMON /gnrl/ (Mh, Mh2, Mh2_final) changes per try: MainEvent::gnrl with gnrl_per_event = true; the radiation and the
// montecarlo use it (the BW rho mass in generate_rad, the decay pion in the spectrometer, recon E and Eloss).
//
// fortran_compat (D-031; D-057, docs/FORTRAN_DIFFERENCES.md):
//   FORTRAN-BUG(D-057 F1)  event.f:325, 1034: trip_thru_target(3) uses the stale hadron angle; corrected: the rho's lab
//                          angles (generate_rho.f:123-124) are used.
//   FORTRAN-BUG(D-057 F2)  event.f:420-422: rho_decay runs after a failed generate_rad on the stale orig (random numbers
//                          drawn, COMMON changed); corrected: only after a successful generation (no weight changes; the
//                          random stream differs after the first such try).
//   + the kernel corrections F3, F6, F9, F20, F21 (physics/rho/RhoPhysics.h) through the same switch.
#ifndef SIMC_REACTIONS_RHOPRODUCTION_H
#define SIMC_REACTIONS_RHOPRODUCTION_H

#include <string>

#include "simc/physics/rho/RhoPhysics.h"
#include "simc/reactions/PionProduction.h"

namespace simc {
namespace reactions {

class RhoProduction : public PionProduction {
public:
    explicit RhoProduction(std::string name = "rho_production");
    void init(const StageContext& ctx) override;
    bool generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) override;
    bool completeEv(MainEvent& main, Event& vertex, IRandom& rng) override;
    bool completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0,
                      const Event& recon) override;
    void prepareCentral(MainEvent&, const Event&) override {}
    void finalTargetState(TargetInfo&) const override {}

    /// event.f:420-422: rho_decay(orig, spec%p%P, main%epsilon, success) as generate calls it (test hook).
    bool decay(MainEvent& main, Event& orig, bool success, IRandom& rng);
    /// Test hooks: the COMMON Mh2 seen by complete_recon_ev; the kernel configurations.
    void setReconMh2(double m) { reconMh2_ = m; }
    rho::GenerateRhoConfig& generateConfig() { return gcfg_; }
    rho::RhoDecayConfig& decayConfig() { return dcfg_; }
    /// What the last generate() did (L1 tables): complete_ev succeeded, generate_rad succeeded (false if not called),
    /// rho_decay called / succeeded / uniforms drawn.
    bool lastCompleteEvOk() const { return lastCompleteEvOk_; }
    bool lastRadSuccess() const { return lastRadSuccess_; }
    bool lastDecayCalled() const { return lastDecayCalled_; }
    bool lastDecaySuccess() const { return lastDecaySuccess_; }
    int lastDecayNUnif() const { return lastDecayNUnif_; }
    /// Make main.gnrl the per-try COMMON /gnrl/ (initialised from the set-up the first time, dbase.f:210-222).
    void initGnrl(MainEvent& main) const;

protected:
    double commonMh2() const override { return reconMh2_; }

private:
    rho::GenerateRhoConfig gcfg_;
    rho::RhoDecayConfig dcfg_;
    double reconMh2_{0};
    bool fixStaleAngle_{true};   ///< D-057 F1
    bool decayOnlyOnSuccess_{true};   ///< D-057 F2
    bool lastCompleteEvOk_{false}, lastRadSuccess_{false}, lastDecayCalled_{false}, lastDecaySuccess_{false};
    int lastDecayNUnif_{0};
};

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_RHOPRODUCTION_H
