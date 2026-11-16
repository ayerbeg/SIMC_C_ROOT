/**
 * @file
 * @ingroup reactions
 * @brief simc::reactions::EepReactionBase: common base of the reactions (complete_ev, complete_recon_ev, complete_main).
 *
 * The code that the (e,e'p) reaction plug-ins share (task P8.2): "hydrogen_elastic"/"phase_space" (HydrogenElastic.h)
 * and "deuterium_eep"/"nucleus_eep" (NuclearEep.h). In the Fortran these are the flag branches of ONE set of routines
 * (event.f), so the base class holds the routines whose branches are selected by the RunSetup flags, exactly as there:
 *
 *   completeEv()     event.f:complete_ev (l.432-1054): branches doing_hyd_elast, doing_deuterium, doing_heavy, doing_phsp
 *                    and the common part (q vector, Pm components, jacobians, trip_thru_target 2/3, radc_init_ev)
 *   completeRecon()  event.f:complete_recon_ev (l.1056-1359), all branches (it is reaction-generic)
 *   completeMain()   event.f:complete_main (l.1363-1569): main%SF_weight = spectralWeight(vertex) (the reaction's own
 *                    part, l.1397-1430), then reactions::assembleWeight (cross section of the StageContext's
 *                    IPhysicsModel, Coulomb factor, weight) — D-030: each reaction owns its weight
 *   genVolume()      simc.f:383-396
 * The derived plug-ins implement generate() (event.f:generate after the beam part) and spectralWeight().
 *
 * Collaborators (StageContext at init(), through their interfaces): the target (trip_thru_target arms 2/3), the
 * radiator (radc_init_ev at the end of complete_ev, generate_rad), the physics model (complete_main).
 *
 * fortran_compat (D-031, docs/FORTRAN_DIFFERENCES.md):
 *   FORTRAN-BUG(D-033 F5)  event.f:926 vertex%PmPer from the previous try's vertex%Pm (D, A): corrected after l.940/943.
 *   FORTRAN-BUG(D-033 F7)  event.f:944 heavy Mrec = sqrt(Emiss**2-Pmiss**2) without a guard (NaN, success = .true.):
 *                          corrected: the try fails (success = false) when Emiss**2 - Pmiss**2 < 0 or NaN.
 *   FORTRAN-DIFF-CANDIDATE(D-033 F4) event.f:1020-1024: no hadron 1/cos**3 jacobian for D(e,e'p) (kept in both modes).
 * FORTRAN-QUIRKs (D-024, D-033 F8): x**2 -> x*x; complete_ev never sets W, epsilon, theta_pq, phi_pq for (e,e'p);
 * ntup%krel = the SAVEd local krel (never assigned on the (e,e'p) branches: 0).
 */
#ifndef SIMC_REACTIONS_EEPREACTIONBASE_H
#define SIMC_REACTIONS_EEPREACTIONBASE_H

#include <iosfwd>
#include <string>

#include "simc/core/RunSetup.h"
#include "simc/reactions/EventChain.h"
#include "simc/stages/IPhysicsModel.h"
#include "simc/stages/IRadiator.h"
#include "simc/stages/IReaction.h"
#include "simc/stages/ITarget.h"

namespace simc {
namespace reactions {

/// @addtogroup reactions
/// @{

/// Fortran COMMON /pfermi/ (simulate.inc): Fermi momentum of the struck nucleon (0 for (e,e'p)) and its energy.
struct FermiMomentum {
    double pfer{0} /**< @copydoc efer */, pferx{0} /**< @copydoc efer */, pfery{0} /**< @copydoc efer */, pferz{0} /**< @copydoc efer */, efer{0};  ///< |p|, px, py, pz [MeV] and energy [MeV] of the struck nucleon.
};

/// ntup%mm, mmA, t written by complete_recon_ev for pion/kaon/delta/rho/semi (event.f:1257-1266).
struct ReconNtup {
    double mm{0} /**< @copydoc t */, mmA{0} /**< @copydoc t */, t{0};  ///< missing mass, missing mass of the nucleus, Mandelstam t (ntuple units).
};

/// @brief Base class of the (e,e'p) and meson reactions: event completion, reconstruction and weight.
///
/// Holds the Fortran per-run state of the reaction (COMMON /pfermi/, SAVEd locals); collaborators are not owned.
/// Not thread-safe; one object per worker.
class EepReactionBase : public IReaction {
public:
    std::string name() const override { return name_; }
    void configure(const Json& config, const SimcParameters& params) override;
    /// Takes the set-up, the target, the radiator and the physics model (StageContext). Derived classes call it.
    /// @param ctx the stage context.
    void init(const StageContext& ctx) override;

    bool completeEv(MainEvent& main, Event& vertex, IRandom& rng) override;
    bool completeRecon(Event& recon) override;
    bool completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0,
                      const Event& recon) override;
    double genVolume() const override;
    void beginTry(long ntried) override { ntried_ = ntried; }

    /// main%SF_weight of the event (event.f:complete_main l.1397-1430).
    /// @param vertex the generated vertex.
    /// @return the spectral-function weight.
    virtual double spectralWeight(const Event& vertex) = 0;

    /// @param os stream for the warnings (not owned), or nullptr.
    void setLog(std::ostream* os) { log_ = os; }
    /// @return COMMON /pfermi/ after the last completeEv().
    const FermiMomentum& fermi() const { return fermi_; }
    /// @return the relative momentum krel of the last completeEv() [MeV].
    double krel() const { return krel_; }
    /// @return the ntuple values written by the last completeRecon().
    const ReconNtup& reconNtup() const { return reconNtup_; }
    /// The CompleteEvFn for IRadiator::radiate (this->completeEv).
    /// @return a callable bound to this object.
    CompleteEvFn completeEvFn();

    // ---- test hooks (L1 tables with harness overrides of the COMMON state) ----
    /// The run constants of complete_main (using_Coulomb, Coulomb%ave, Ebeam).
    /// @return the context (modifiable).
    MainContext& mainContext() { return mainCtx_; }
    /// Replace the cross-section model (e.g. a fixed deForest_flag).
    /// @param m the model (not owned).
    void setPhysics(IPhysicsModel* m) { physics_ = m; }
    /// @return the try number set by beginTry().
    long ntried() const { return ntried_; }

protected:
    /// COMMON /gnrl/ Mh2 as complete_recon_ev sees it (event.f:1329 ntup%t): the set-up value; the rho reaction returns
    /// the per-try value (the decay pion, task P13.2).
    /// @return Mh2 [MeV^2].
    virtual double commonMh2() const { return setup().Mh2; }

    /// @param name registry name of the reaction.
    explicit EepReactionBase(std::string name);
    /// @return the run set-up (after init()).
    const RunSetup& setup() const;

    std::string name_;  ///< registry name.
    const RunSetup* s_{nullptr};  ///< the run set-up (not owned).
    ITarget* target_{nullptr};          ///< not owned (StageContext::target)
    IRadiator* radiator_{nullptr};      ///< not owned (StageContext::radiator)
    IPhysicsModel* physics_{nullptr};   ///< not owned (StageContext::physics)
    std::ostream* log_;  ///< stream for the warnings (not owned).
    FermiMomentum fermi_;      ///< COMMON /pfermi/
    double krel_{0};           ///< complete_ev SAVEd local krel (-fno-automatic)
    ReconNtup reconNtup_;  ///< ntuple values of the last completeRecon().
    double targ_pol_{0};  ///< target polarization.
    long ntried_{0};           ///< COMMON ntried (beginTry)
    MainContext mainCtx_;  ///< the run constants of complete_main.
};

/// @}

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_EEPREACTIONBASE_H
