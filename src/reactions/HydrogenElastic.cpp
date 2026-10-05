// src/reactions/HydrogenElastic.cpp — see include/simc/reactions/HydrogenElastic.h.
// Fortran: simc_gfortran/event.f. Operation order follows the Fortran source; x**2 is x*x (D-024).
#include "simc/reactions/HydrogenElastic.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "simc/core/Errors.h"
#include "simc/core/Event.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace reactions {

HydrogenElastic::HydrogenElastic(std::string name) : EepReactionBase(std::move(name)) {}

void HydrogenElastic::init(const StageContext& ctx) {
    const ReactionFlags& f = ctx.setup.flags;
    if (f.doing_deuterium || f.doing_heavy)
        throw std::invalid_argument("reaction '" + name_ + "': D(e,e'p)/A(e,e'p) is the reaction 'deuterium_eep' / "
                                    "'nucleus_eep'");
    if (f.doing_pion || f.doing_kaon || f.doing_delta || f.doing_rho || f.doing_semi)
        throw NotPorted("reaction '" + name_ +
                                    "': only H(e,e'p) (doing_hyd_elast) and doing_phsp are ported (event.f, P4.1)");
    if (!f.doing_hyd_elast && !f.doing_phsp)
        throw NotPorted("reaction '" + name_ + "': the setup is neither doing_hyd_elast nor doing_phsp");
    EepReactionBase::init(ctx);
}

// Fortran: event.f:complete_main l.1397-1398 (doing_hyd_elast, doing_phsp: no spectral function)
double HydrogenElastic::spectralWeight(const Event&) { return 1.0; }

// Fortran: event.f:generate (l.214-428), after the beam part (StandardBeam).
bool HydrogenElastic::generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) {
    const RunSetup& s = setup();

    // l.216-217
    bool success = false;
    main.gen_weight = 1.0;

    // Generate Electron Angles (all cases) (l.276-277)
    vertex.e.yptar = s.gen.e.yptar.min + rng.uniform() * (s.gen.e.yptar.max - s.gen.e.yptar.min);
    vertex.e.xptar = s.gen.e.xptar.min + rng.uniform() * (s.gen.e.xptar.max - s.gen.e.xptar.min);
    // (hadron angles, hadron and electron energies: not generated for H(e,e'p) / phsp, l.280-313)

    // l.320-325. FORTRAN-QUIRK(D-024): the hadron call uses the previous try's p%xptar/yptar for H elastic;
    // complete_ev overwrites the result.
    PhysicsAngles(s.spec.e.theta, s.spec.e.phi, vertex.e.xptar, vertex.e.yptar, vertex.e.theta, vertex.e.phi, log_);
    PhysicsAngles(s.spec.p.theta, s.spec.p.phi, vertex.p.xptar, vertex.p.yptar, vertex.p.theta, vertex.p.phi, log_);

    // l.329-335 (/pfermi/ and Em; the Fermi-motion generation l.336-372 is for A(e,e'pi/K) only)
    fermi_.pfer = 0.0;
    fermi_.pferx = 0.0;
    fermi_.pfery = 0.0;
    fermi_.pferz = 0.0;
    vertex.Em = 0.0;
    fermi_.efer = s.targ.Mtar_struck;

    // l.380-393
    success = completeEv(main, vertex, rng);
    main.sigcc = 1.0;
    if (!success) return false;
    main.Trec = vertex.Trec;

    // l.400-415
    // FORTRAN-DIFF-CANDIDATE(D-024): when generate_rad fails, orig is not written and keeps the previous try's values
    // (event.f:400-415); harmless for the results (the try is rejected), visible only in per-try outputs.
    if (s.sw.using_rad) {
        success = radiator_->radiate(main, vertex, orig, rng, completeEvFn());
    } else {
        success = true;   // (the doing_heavy VERTEXedge test is not reachable here)
        orig = vertex;
    }
    return success;   // (rho decay l.421-423 not reachable)
}

namespace {
struct PhaseSpace : HydrogenElastic {
    PhaseSpace() : HydrogenElastic("phase_space") {}
};
}  // namespace

SIMC_REGISTER(IReaction, "hydrogen_elastic", HydrogenElastic)
SIMC_REGISTER(IReaction, "phase_space", PhaseSpace)

}  // namespace reactions
}  // namespace simc
