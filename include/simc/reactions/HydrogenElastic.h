// include/simc/reactions/HydrogenElastic.h
//
// IReaction "hydrogen_elastic" and "phase_space" (task P4.1): event.f for H(e,e'p) and doing_phsp.
//   generate()        event.f:generate l.214-428 (after the beam part): electron angles, physics_angles,
//                     /pfermi/ reset, complete_ev, main%Trec, generate_rad (IRadiator) or orig = vertex
//   completeEv()      event.f:complete_ev l.432-1054 (branches doing_hyd_elast and doing_phsp + the common part)
//   completeRecon()   event.f:complete_recon_ev l.1056-1359 (all branches; it is reaction-generic)
//   genVolume()       simc.f:383-396 (domega_e for H(e,e'p) / phsp)
//
// Both registry names create this class: like the Fortran, the code path is chosen by the RunSetup flags
// (doing_hyd_elast, doing_phsp, ...), not by the name. D-009: hydrogen + doing_phsp keeps doing_hyd_elast, so
// complete_ev takes the elastic branch (event.f:546 before 709), while complete_recon_ev overwrites the hadron
// momentum with spec%p%P (event.f:1132-1136) and complete_main sets sigcc = 1. Reactions whose generate/complete_ev
// branches are not ported (deuterium, A>1, pion, kaon, delta, rho, semi) are refused in init().
//
// Collaborators (set by the driver before init): the target (trip_thru_target for arms 2/3) and the radiator
// (radc_init_ev at the end of complete_ev, generate_rad when using_rad). complete_ev is handed to the radiator as
// the CompleteEvFn (D-022).
//
// FORTRAN-QUIRKs (D-024): x**2 -> x*x; complete_ev never sets W, epsilon, theta_pq, phi_pq for H; the only failure
// of the elastic branch is e%E > Ein (theta = 0 rounding); at theta = 0 the q-normalisation check sees NaN and does
// not STOP; physics_angles for the hadron arm uses the previous try's p%xptar/yptar (result overwritten);
// ntup%krel = the SAVEd local krel (never assigned for H: 0); orig is not written when generate_rad fails.
#ifndef SIMC_REACTIONS_HYDROGENELASTIC_H
#define SIMC_REACTIONS_HYDROGENELASTIC_H

#include <iosfwd>
#include <string>

#include "simc/core/RunSetup.h"
#include "simc/stages/IRadiator.h"
#include "simc/stages/IReaction.h"
#include "simc/stages/ITarget.h"

namespace simc {
namespace reactions {

/// Fortran COMMON /pfermi/ (simulate.inc): Fermi momentum of the struck nucleon (0 for H) and its energy.
struct FermiMomentum {
    double pfer{0}, pferx{0}, pfery{0}, pferz{0}, efer{0};
};

/// ntup%mm, mmA, t written by complete_recon_ev for pion/kaon/delta/rho/semi (event.f:1257-1266).
struct ReconNtup {
    double mm{0}, mmA{0}, t{0};
};

class HydrogenElastic : public IReaction {
public:
    explicit HydrogenElastic(std::string name = "hydrogen_elastic");

    std::string name() const override { return name_; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(RunSetup& setup) override;

    bool generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) override;
    bool completeEv(MainEvent& main, Event& vertex, IRandom& rng) override;
    bool completeRecon(Event& recon) override;
    double genVolume() const override;

    void setTarget(ITarget* t) { target_ = t; }
    void setRadiator(IRadiator* r) { radiator_ = r; }
    void setLog(std::ostream* os) { log_ = os; }

    const FermiMomentum& fermi() const { return fermi_; }
    double krel() const { return krel_; }
    const ReconNtup& reconNtup() const { return reconNtup_; }
    /// The CompleteEvFn for IRadiator::radiate (this->completeEv).
    CompleteEvFn completeEvFn();

private:
    std::string name_;
    const RunSetup* s_{nullptr};
    ITarget* target_{nullptr};
    IRadiator* radiator_{nullptr};
    std::ostream* log_;
    FermiMomentum fermi_;      ///< COMMON /pfermi/
    double krel_{0};           ///< complete_ev SAVEd local krel (-fno-automatic)
    ReconNtup reconNtup_;
    double targ_pol_{0};
};

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_HYDROGENELASTIC_H
