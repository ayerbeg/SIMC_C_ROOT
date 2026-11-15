/**
 * @file
 * @ingroup stages
 * @brief simc::IReaction: the reaction stage (event generation, kinematics completion, event weight).
 * @note Origin: Fortran SIMC event.f:generate (l.126-429), complete_ev (l.432-1054), complete_recon_ev
 *       (l.1056-1360), complete_main (l.1363-1569); simc.f genvol (l.383-396).
 */
#ifndef SIMC_STAGES_IREACTION_H
#define SIMC_STAGES_IREACTION_H

#include "simc/core/RunSetup.h"
#include "simc/stages/IStage.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief A reaction: generates the physics variables of a try, completes the kinematics and assembles the event
 *        weight.
 *
 * Per try the driver calls beginTry(), the beam (simc::IBeam::generate), generate() (which completes the vertex and
 * radiates it), the spectrometers, completeRecon() on the reconstructed event and completeMain() for the weight.
 * A reaction takes the radiator, the target and the physics model from its simc::StageContext. Implementations live
 * in simc::reactions (hydrogen elastic, nuclear (e,e'p), pion, kaon, Delta, SIDIS, rho, phase space). A reaction keeps
 * state between tries (as the original program does); not thread-safe.
 */
class IReaction : public IStage {
public:
    /**
     * @brief Generates the variables of this reaction, completes the vertex and radiates it.
     *
     * The beam part (IBeam::generate) must have been called first for this try. @p orig is written only where the
     * original program writes it (it keeps the previous try's values when the radiation fails, D-024).
     * @param main The quantities of the try (weights, target, ntuple fields).
     * @param vertex Out: the vertex kinematics.
     * @param orig Out: the kinematics after radiation (orig = vertex without radiation).
     * @param rng The random source.
     * @return False if the try fails (outside the generation limits, unphysical kinematics, failed radiation).
     * @note Origin: Fortran SIMC event.f:generate after the beam part (l.209-428).
     */
    virtual bool generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) = 0;
    /**
     * @brief Completes the vertex kinematics from the generated variables.
     *
     * Also the simc::CompleteEvFn handed to IRadiator::radiate.
     * @param main The quantities of the try (target energy losses are written).
     * @param vertex In/out: the vertex kinematics.
     * @param rng The random source (energy-loss draws).
     * @return False if the kinematics is unphysical (the try fails).
     * @note Origin: Fortran SIMC event.f:complete_ev (l.432-1054).
     */
    virtual bool completeEv(MainEvent& main, Event& vertex, IRandom& rng) = 0;
    /**
     * @brief Completes the reconstructed event from the reconstructed arm quantities.
     * @param recon In/out: the reconstructed kinematics (arm quantities set by the driver).
     * @return False if the reconstructed kinematics is unphysical.
     * @note Origin: Fortran SIMC event.f:complete_recon_ev.
     */
    virtual bool completeRecon(Event& recon) = 0;
    /**
     * @brief Assembles the event weight.
     *
     * main.SF_weight (the reaction's spectral function), main.sigcc / sigcc_recon (the run's simc::IPhysicsModel), the
     * Coulomb focusing factor and main.weight = SF_weight * jacobian * gen_weight * sigcc * (nucleon count). Task
     * P8.2, decision D-030.
     * @param force_sigcc True for the central event: the cross section is computed even if the event cannot
     *        contribute.
     * @param main In/out: the quantities of the try; weight, SF_weight, sigcc, sigcc_recon are written.
     * @param vertex The vertex kinematics.
     * @param vertex0 The central vertex kinematics (used by some models for reference quantities).
     * @param recon The reconstructed kinematics.
     * @return False (main.weight not written) when the event cannot contribute (SF_weight <= 0), unless
     *         @p force_sigcc.
     * @note Origin: Fortran SIMC event.f:complete_main (l.1363-1569).
     */
    virtual bool completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0,
                              const Event& recon) = 0;
    /**
     * @brief The generated phase-space volume entering the normalisation factor.
     * @return The product of the generated ranges of this reaction's variables (units of the jacobians).
     * @note Origin: Fortran SIMC simc.f:383-396 (genvol).
     */
    virtual double genVolume() const = 0;
    /**
     * @brief Lets a reaction set the central-event quantities its cross section reads.
     *
     * Called while the central event is computed, after completeRecon(vertex0) and before completeMain(true, main0,
     * ...). The pion reaction fills main0.epsilon, t, phi_pq in the corrected mode (D-038 F5). Default: nothing.
     * @param main0 In/out: the central-event quantities.
     * @param vertex0 The central vertex kinematics.
     * @note Origin: Fortran SIMC simc.f:calculate_central (l.1143-1300).
     */
    virtual void prepareCentral(MainEvent& main0, const Event& vertex0) {
        (void)main0;
        (void)vertex0;
    }
    /**
     * @brief Writes the target values the reaction changed during the run, as the report prints them.
     *
     * Called once at the end of the run, before the reports (e.g. Mrec_struck = the last Breit-Wigner Delta mass of the
     * pion Delta final states, simc.f:844). Default: nothing.
     * @param targ In/out: the target description printed in the report.
     */
    virtual void finalTargetState(TargetInfo& targ) const { (void)targ; }
    /**
     * @brief Announces the try counter of the coming try (before generate()).
     *
     * Only used in messages (event.f:591, the 'two solutions' warning for ntried <= 5000). Default: ignored.
     * @param ntried The number of the coming try (`ntried`, simc.f:190).
     */
    virtual void beginTry(long ntried) { (void)ntried; }
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_IREACTION_H
