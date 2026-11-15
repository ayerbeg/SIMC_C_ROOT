/**
 * @file
 * @ingroup stages
 * @brief simc::IPhysicsModel: the cross-section model of a reaction.
 * @note Origin: Fortran SIMC physics_proton.f (sigep, deForest), physics_pion.f, physics_kaon.f, ..., called from
 *       event.f:complete_main (l.1363-1570).
 */
#ifndef SIMC_STAGES_IPHYSICSMODEL_H
#define SIMC_STAGES_IPHYSICSMODEL_H

#include "simc/stages/IStage.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief Cross-section model of a reaction: the differential cross section of a try, at the vertex and at the
 *        reconstructed kinematics.
 *
 * The reaction stage calls it from its completeMain() (simc::IReaction) to assemble the event weight; the model is
 * handed to the reaction through simc::StageContext::physics. Implementations: the proton models ("deforest",
 * "deforest_cc1", "deforest_cc1_onshell", "deforest_cc2", "sigep"), simc::pion::PionModel, simc::kaon::KaonModel, simc::delta::DeltaModel, simc::sidis::SidisModel,
 * simc::rho::RhoModel. The units of the cross section are those of the model (documented per model); they only have
 * to be consistent with the reaction's jacobians. A model may keep state between calls; not thread-safe.
 */
class IPhysicsModel : public IStage {
public:
    /**
     * @brief The cross section of a try at the vertex kinematics (main.sigcc).
     * @param vertex The vertex kinematics.
     * @param main The weights and quantities of the try computed so far.
     * @return The cross section (units of the model).
     */
    virtual double sigcc(const Event& vertex, const MainEvent& main) = 0;
    /**
     * @brief The cross section at the reconstructed kinematics (main.sigcc_recon).
     * @param recon The reconstructed kinematics.
     * @param main The weights and quantities of the try.
     * @return The cross section, or 1.0 where the original program does not compute it.
     */
    virtual double sigccRecon(const Event& recon, const MainEvent& main) = 0;
    /**
     * @brief The vertex cross section, for the models that also WRITE fields of the try (what the reaction calls).
     *
     * The meson models set main.thetacm, phicm, pcm, davejac, johnjac, wcm and ntup.sigcm1/sigcm2/sigcm as their
     * original functions do (task P9.1).
     * @param vertex The vertex kinematics.
     * @param main The quantities of the try; the model may write its centre-of-mass and ntuple fields.
     * @return The cross section (units of the model). Default: sigcc(vertex, main), nothing written.
     * @note Origin: Fortran SIMC physics_pion.f:peepi (main%thetacm, ...).
     */
    virtual double sigccMain(const Event& vertex, MainEvent& main) { return sigcc(vertex, main); }
    /**
     * @brief Called by the driver before each event with the event counter of the original program.
     *
     * Models use it only in warning messages (e.g. deForest, physics_proton.f:75). Default: ignored.
     * @param nevent The number of the coming event (`nevent`, simc.f:165-349).
     */
    virtual void beginEvent(long nevent) { (void)nevent; }
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_IPHYSICSMODEL_H
