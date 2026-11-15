// src/reactions/delta/DeltaProduction.cpp — see include/simc/reactions/DeltaProduction.h.
// The "delta_production" reaction, H/D/A(e,e'p)pi0 through the Delta(1232): the proton is detected, the pion is missing.
// Fortran: simc_gfortran/event.f:complete_main (l.1363-1569, doing_delta branch l.1511-1513). generate, complete_ev and
// complete_recon_ev are PionProduction / EepReactionBase with the delta flags.
#include "simc/reactions/DeltaProduction.h"

#include <stdexcept>

#include "simc/core/Registry.h"

namespace simc {
namespace reactions {

DeltaProduction::DeltaProduction(std::string name) : PionProduction(std::move(name), false) { doingDelta_ = true; }

// Fortran: event.f:complete_main (l.1363-1569), doing_delta branch.
bool DeltaProduction::completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& /*vertex0*/,
                                   const Event& recon) {
    if (!physics_) throw std::logic_error("reaction '" + name_ + "': no physics model (StageContext::physics)");
    // ntup%mm, mmA, t of complete_recon_ev (event.f:1323-1335, COMMON /ntuple_variables/)
    main.ntup.mm = reconNtup_.mm;
    main.ntup.mmA = reconNtup_.mmA;
    main.ntup.t = reconNtup_.t;

    main.SF_weight = 1.0;   // l.1397-1398
    if (main.SF_weight <= 0 && !force_sigcc) return false;   // l.1436

    const double tgtweight = 1.0;   // l.1442; not set on the doing_delta branch
    // l.1512: `main%sigcc = peedelta(vertex,main)	!Need new xsec model.`
    // FORTRAN-DIFF-CANDIDATE(D-048 F4): peedelta returns 1.0*jacobian*(gtpr*fac), no cross section (see DeltaModel.h).
    main.sigcc = physics_->sigccMain(vertex, main);
    main.sigcc_recon = physics_->sigccRecon(recon, main);   // l.1513: 1.0

    // Coulomb focusing factor (l.1533-1535)
    if (mainCtx_.using_Coulomb) {
        const double c = 1.0 + mainCtx_.Coulomb_ave / mainCtx_.Ebeam;
        main.sigcc = main.sigcc * (c * c);
    }

    // l.1557-1559
    main.weight = main.SF_weight * main.jacobian * main.gen_weight * main.sigcc;
    main.weight = main.weight * tgtweight;   // correct for #/nucleons involved
    return true;
}

SIMC_REGISTER(IReaction, "delta_production", DeltaProduction)

}  // namespace reactions
}  // namespace simc
