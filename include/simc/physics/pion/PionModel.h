// include/simc/physics/pion/PionModel.h
//
// IPhysicsModel plug-in "pion" (task P9.1): the exclusive pion cross section of the Fortran SIMC,
// physics_pion.f:peepi, called by event.f:complete_main (l.1458-1505) for doing_pion:
//     main%sigcc       = peepi(vertex, main)       (sigccMain / sigcc)
//     main%sigcc_recon = 1.0                       (sigccRecon)
// The model choice is the Fortran's, from the run flags: P. Bosted's 2021 fit (sig_param_2021) for every final state,
// blended with MAID2007 (sigmaid) for W < 2 GeV; pi+/pi-/pi0 and the struck nucleon from which_pion and doing_pizero
// (setup.params.which_pion, setup.flags.doing_pizero), the nucleon mass targ%Mtar_struck (setup.targ.Mtar_struck).
// The Delta final-state factors (which_pion 2/3) and tgtweight (Z or N) belong to complete_main (the reaction).
//
// JSON config (all optional):
//     { "data_dir": "<dir with maid*.dat>" (default data/pion of the source tree, SIMC_PION_DATA_DIR),
//       "verbose": true }
// Fortran-compat (D-031): setup.fortran_compat; only D-038 F10 (stale ntup%sigcm2 for W >= 2 GeV) depends on it.
#ifndef SIMC_PHYSICS_PION_PIONMODEL_H
#define SIMC_PHYSICS_PION_PIONMODEL_H

#include <ostream>
#include <string>

#include "simc/physics/pion/PionCrossSections.h"
#include "simc/stages/IPhysicsModel.h"

namespace simc {
namespace pion {

class PionModel : public IPhysicsModel {
public:
    std::string name() const override { return "pion"; }
    void configure(const Json& config, const SimcParameters& params) override;
    /// Reads the MAID tables (the Fortran reads them at the first sigmaid call). Throws std::invalid_argument if they
    /// cannot be read or if targ%Mtar_struck is not set (not a pion set-up).
    void init(const StageContext& ctx) override;

    /// peepi(vertex, main) without writing main (evaluated on a copy).
    double sigcc(const Event& vertex, const MainEvent& main) override;
    /// main%sigcc = peepi(vertex, main): also writes main%thetacm, phicm, pcm, davejac, johnjac, wcm, ntup%sigcm*.
    double sigccMain(const Event& vertex, MainEvent& main) override;
    /// event.f:1505 main%sigcc_recon = 1.0
    double sigccRecon(const Event&, const MainEvent&) override { return 1.0; }

    const PeepiConfig& config() const { return cfg_; }
    const MaidTables& maid() const { return maid_; }

private:
    std::string data_dir_;
    bool verbose_{true};
    PeepiConfig cfg_;
    MaidTables maid_;
};

}  // namespace pion
}  // namespace simc

#endif  // SIMC_PHYSICS_PION_PIONMODEL_H
