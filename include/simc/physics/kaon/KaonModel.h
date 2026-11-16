/**
 * @file
 * @ingroup physics
 * @brief simc::kaon::KaonModel: IPhysicsModel "kaon", the kaon electroproduction cross section (physics_kaon.f:peeK).
 *
 * IPhysicsModel plug-in "kaon" (task P10.1): the kaon electroproduction cross section of the Fortran SIMC,
 * physics_kaon.f:peeK, called by event.f:complete_main (l.1502-1509) for doing_kaon:
 *     main%sigcc       = peeK(vertex, main, survivalprob)   (sigccMain / sigcc)
 *     main%sigcc_recon = 1.0                                (sigccRecon)
 * The weight uses the factorized model (sig_factorized, ntup%sigcm2 = ntup%sigcm); the Saghai model (eekeek/eekeeks
 * with the /sigkaon/ tables) only fills ntup%sigcm1 (the "saghai" ntuple column). When the decay is not simulated
 * (doing_decay = 0), peeK also returns the survival probability of the kaon; it is written to the COMMON decdist
 * (main.ntup.decdist, physics_kaon.f:164), from where the reaction applies it to the weight (event.f:1560-1561).
 * The model reads the run state from the set-up: targ%Mtar_struck, targ%Mrec_struck, doing_decay, hadron_arm, ctau,
 * spec%p%P, Mh2.
 *
 * JSON config (all optional):
 *     { "data_dir": "<dir with weightc.dat weightd.dat saghai_proton.dat saghai_sigma0.dat>" (default data/kaon of the
 *       source tree, SIMC_KAON_DATA_DIR), "verbose": true }
 * fortran_compat (D-031, setup.fortran_compat) restores the Fortran as built; the default corrects D-043 F1 (fint
 * value), F2 (saghai_sigma0.dat records), F3 (eekeek s nodes), F5 (survival probability from the kaon momentum) and F6
 * (zaero for the SHMS, sqrt in the path). docs/FORTRAN_DIFFERENCES.md.
 */
#ifndef SIMC_PHYSICS_KAON_KAONMODEL_H
#define SIMC_PHYSICS_KAON_KAONMODEL_H

#include <string>

#include "simc/physics/kaon/KaonCrossSections.h"
#include "simc/stages/IPhysicsModel.h"

namespace simc {
namespace kaon {

/// @addtogroup physics
/// @{

/// @brief IPhysicsModel "kaon": peeK with the Saghai tables. Owns the tables and the SAVEd state; not thread-safe.
class KaonModel : public IPhysicsModel {
public:
    std::string name() const override { return "kaon"; }
    void configure(const Json& config, const SimcParameters& params) override;
    /// Reads the /sigkaon/ tables (dbase.f:624-679). Throws std::invalid_argument if they cannot be read or if
    /// targ%Mtar_struck is not set (not a kaon set-up).
    /// @param ctx the stage context.
    void init(const StageContext& ctx) override;

    /// peeK(vertex, main) on copies (main and the SAVEd state are not changed).
    /// @param vertex the vertex.
    /// @param main the main event record (copied, not written).
    /// @return the cross section.
    double sigcc(const Event& vertex, const MainEvent& main) override;
    /// main%sigcc = peeK(vertex, main, survivalprob): also writes main%thetacm ... wcm, ntup%sigcm* and (doing_decay
    /// = 0) the survival probability to ntup%decdist; updates the SAVEd locals (zaero, survivalprob).
    /// @param vertex the vertex.
    /// @param main the main event record (written).
    /// @return the cross section.
    double sigccMain(const Event& vertex, MainEvent& main) override;
    /// event.f:1504 main%sigcc_recon = 1.0
    /// @return 1.
    double sigccRecon(const Event&, const MainEvent&) override { return 1.0; }

    /// @return the peeK run constants.
    const PeekConfig& config() const { return cfg_; }
    /// @return the peeK run constants (modifiable).
    PeekConfig& config() { return cfg_; }   ///< test hook (L1 rows override doing_decay/hadron_arm/ctau)
    /// @return the Saghai tables.
    const SaghaiTables& tables() const { return tab_; }
    /// @return the SAVEd state of peeK.
    PeekState& state() { return st_; }
    /// The SAVEd local survivalprob of complete_main (event.f:1372), updated by peeK when doing_decay = 0.
    /// @return the survival probability.
    double survivalProbability() const { return survivalprob_; }
    /// @param p the survival probability.
    void setSurvivalProbability(double p) { survivalprob_ = p; }

private:
    std::string data_dir_;
    bool verbose_{true};
    PeekConfig cfg_;
    SaghaiTables tab_;
    PeekState st_;
    double survivalprob_{0};
};

/// @}

}  // namespace kaon
}  // namespace simc

#endif  // SIMC_PHYSICS_KAON_KAONMODEL_H
