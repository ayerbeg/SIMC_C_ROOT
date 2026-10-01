// include/simc/reactions/EventChain.h
//
// The event-level routines of event.f that are not a single stage (task P4.1):
//   generateEvent   event.f:generate (l.126-428) = IBeam::generate (beam/vertex part, l.150-212) followed by
//                   IReaction::generate (generated variables, complete_ev, generate_rad, l.214-428)
//   completeMain    event.f:complete_main (l.1363-1569): spectral-function weight, sigcc/sigcc_recon from the
//                   IPhysicsModel, the Coulomb focusing factor and the event weight. It is the weight assembly of
//                   the driver's chain (docs/ARCHITECTURE.md §4 completeMain), kept here as a free function so the
//                   driver (P7) calls it with its own model; only the H(e,e'p)/phsp spectral-function branch is ported.
//   limitsUpdate    event.f:limits_update / update_range / update_slop_range (l.1-125), with Contrib = contribtype
//   transformToCm   jacobians.f:transform_to_cm (l.1-282) (pion/kaon/delta cross sections; not on the H path)
#ifndef SIMC_REACTIONS_EVENTCHAIN_H
#define SIMC_REACTIONS_EVENTCHAIN_H

#include <array>

#include "simc/core/Event.h"
#include "simc/core/RunSetup.h"
#include "simc/stages/IBeam.h"
#include "simc/stages/IPhysicsModel.h"
#include "simc/stages/IReaction.h"

namespace simc {
namespace reactions {

/// Fortran: event.f:generate(main, vertex, orig, success).
bool generateEvent(IBeam& beam, IReaction& reaction, MainEvent& main, Event& vertex, Event& orig, IRandom& rng);

/// Run constants read by complete_main.
struct MainContext {
    ReactionFlags flags;
    bool using_Coulomb{false};   ///< after dbase.f:547 (forced off for Z = 1)
    double Coulomb_ave{0};       ///< targ%Coulomb%ave [MeV]
    double Ebeam{0};             ///< [MeV]
    static MainContext FromSetup(const RunSetup& s);
};

/// Fortran: event.f:complete_main(force_sigcc, main, vertex, vertex0, recon, success).
/// Throws std::invalid_argument for reactions whose spectral-function branch is not ported (D(e,e'p), A(e,e'p)).
bool completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0, const Event& recon,
                  IPhysicsModel& model, const MainContext& ctx);

/// Fortran: structures.inc rangetype (lo, hi); min_max_init starts lo = +1e10, hi = -1e10.
struct Range {
    double lo{1.0e10}, hi{-1.0e10};
};
struct ArmRange {
    Range delta, yptar, xptar, z;
};
/// Fortran: structures.inc contribtype (init.f:min_max_init l.1080-1190).
struct Contrib {
    struct { ArmRange e, p; Range Trec, sumEgen; } gen;
    struct { struct { Range E, yptar, xptar; } e, p; Range Em, Pm, Trec; } tru;
    struct { ArmRange e, p; } SP;
    struct { Range Trec, Em, Pm; } vertex;
    struct { std::array<Range, 3> Egamma; Range Egamma_total; } rad;   ///< Egamma(1:3) -> [0..2]
};

/// Fortran: event.f:update_range (l.95-108).
void updateRange(double val, Range& range);
/// Fortran: event.f:update_slop_range (l.112-125).
void updateSlopRange(double val, SlopItem& slop);
/// Fortran: event.f:limits_update (l.1-91). Egamma_used(1:3) = the radiator's per-event values (index 1..3 used).
void limitsUpdate(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                  const ReactionFlags& f, const std::array<double, 4>& Egamma_used, Contrib& contrib, Slops& slop);

/// Outputs of jacobians.f:transform_to_cm.
struct CmTransform {
    double gstar{0}, bstar{0}, bstarx{0}, bstary{0}, bstarz{0};
    double nustar{0}, qstar{0}, qstarx{0}, qstary{0}, qstarz{0};
    double ehadcm{0}, phadcm{0}, phadcmx{0}, phadcmy{0}, phadcmz{0};
    double ebeamcm{0}, pbeamcm{0}, pbeamcmx{0}, pbeamcmy{0}, pbeamcmz{0};
    double etarcm{0}, ptarcm{0}, ptarcmx{0}, ptarcmy{0}, ptarcmz{0};
    double thetacm{0}, phicm{0}, phiqn{0}, jacobian{0}, jac_old{0};
};
/// Fortran: jacobians.f:transform_to_cm (l.1-282). Reads vertex Ein, nu, q, p%P, p%E, up, uq, main%phi_pq and
/// COMMON /pfermi/ (pfer, pferx..z, efer). The |beta|>1 warning goes to `log` (if non-null).
CmTransform transformToCm(const Event& vertex, const MainEvent& main, double pfer, double pferx, double pfery,
                          double pferz, double efer, std::ostream* log = nullptr);

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_EVENTCHAIN_H
