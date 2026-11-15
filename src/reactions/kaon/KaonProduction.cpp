// src/reactions/kaon/KaonProduction.cpp — see include/simc/reactions/KaonProduction.h.
// The "kaon_production" reaction (K+ Lambda/Sigma0, Sigma-, coherent hypernuclei): PionProduction with the kaon flags;
// complete_main differs (the kaon cross section, the target weight and the survival weight).
// Fortran: simc_gfortran/event.f:complete_main (l.1363-1569, doing_kaon branch).
#include "simc/reactions/KaonProduction.h"

#include <stdexcept>

#include "simc/core/Registry.h"

namespace simc {
namespace reactions {

KaonProduction::KaonProduction(std::string name) : PionProduction(std::move(name), true) {}

// Fortran: event.f:complete_main (l.1363-1569), doing_kaon branch.
bool KaonProduction::completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& /*vertex0*/,
                                  const Event& recon) {
    if (!physics_) throw std::logic_error("reaction '" + name_ + "': no physics model (StageContext::physics)");
    const RunSetup& s = setup();
    // ntup%mm, mmA, t of complete_recon_ev (event.f:1323-1335, COMMON /ntuple_variables/)
    main.ntup.mm = reconNtup_.mm;
    main.ntup.mmA = reconNtup_.mmA;
    main.ntup.t = reconNtup_.t;

    main.SF_weight = 1.0;   // l.1397-1398
    if (main.SF_weight <= 0 && !force_sigcc) return false;   // l.1436

    double tgtweight = 1.0;
    main.sigcc = physics_->sigccMain(vertex, main);   // l.1503 peeK(vertex, main, survivalprob)
    // peeK sets survivalprob (and decdist = survivalprob) only when the decay is not simulated (physics_kaon.f:148-165)
    if (!s.params.doing_decay) survivalprob_ = main.ntup.decdist;
    main.sigcc_recon = physics_->sigccRecon(recon, main);   // l.1504: 1.0
    // FORTRAN-DIFF-CANDIDATE(D-043 F9): "OK for coherent???" (event.f:1505): Z or N also for coherent production.
    const int w = s.params.which_kaon;
    if (w == 2 || w == 12)
        tgtweight = s.targ.N;
    else
        tgtweight = s.targ.Z;

    // Coulomb focusing factor (l.1533-1535)
    if (mainCtx_.using_Coulomb) {
        const double c = 1.0 + mainCtx_.Coulomb_ave / mainCtx_.Ebeam;
        main.sigcc = main.sigcc * (c * c);
    }

    // l.1557-1561
    main.weight = main.SF_weight * main.jacobian * main.gen_weight * main.sigcc;
    main.weight = main.weight * tgtweight;   // correct for #/nucleons involved
    if (!s.params.doing_decay) main.weight = main.weight * survivalprob_;   // survival probability (no decay simulated)
    return true;
}

SIMC_REGISTER(IReaction, "kaon_production", KaonProduction)

}  // namespace reactions
}  // namespace simc
