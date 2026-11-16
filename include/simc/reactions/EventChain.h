/**
 * @file
 * @ingroup reactions
 * @brief The event chain shared by the reactions: generate, complete_main weight, limits_update, transform_to_cm.
 *
 * The event-level routines of event.f that are not a single stage (task P4.1):
 *   generateEvent   event.f:generate (l.126-428) = IBeam::generate (beam/vertex part, l.150-212) followed by
 *                   IReaction::generate (generated variables, complete_ev, generate_rad, l.214-428)
 *   assembleWeight  event.f:complete_main (l.1432-1563), the reaction-independent part after the spectral-function
 *                   weight: the early return for SF_weight <= 0, sigcc/sigcc_recon from the IPhysicsModel, the Coulomb
 *                   focusing factor and the event weight. Since P8.2 (D-030) every reaction owns its weight
 *                   (IReaction::completeMain): it sets main%SF_weight and calls this kernel.
 *   limitsUpdate    event.f:limits_update / update_range / update_slop_range (l.1-125), with Contrib = contribtype
 *   transformToCm   jacobians.f:transform_to_cm (l.1-282) (pion/kaon/delta cross sections; not on the H path)
 */
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

/// @addtogroup reactions
/// @{

/// Fortran: event.f:generate(main, vertex, orig, success).
/// @param beam the beam stage.
/// @param reaction the reaction stage.
/// @param main the main event record (written).
/// @param vertex the generated vertex (written).
/// @param orig the event after radiation (written; not written when generate_rad fails).
/// @param rng the random generator.
/// @return true if the event was generated (Fortran success).
bool generateEvent(IBeam& beam, IReaction& reaction, MainEvent& main, Event& vertex, Event& orig, IRandom& rng);

/// Run constants read by complete_main.
struct MainContext {
    ReactionFlags flags;  ///< reaction flags of the run.
    bool using_Coulomb{false};   ///< after dbase.f:547 (forced off for Z = 1)
    double Coulomb_ave{0};       ///< targ%Coulomb%ave [MeV]
    double Ebeam{0};             ///< [MeV]
    /// @param s the run set-up.
    /// @return the context.
    static MainContext FromSetup(const RunSetup& s);
};

/// Fortran: event.f:complete_main (l.1432-1563) after main%SF_weight has been set by the reaction (l.1397-1430):
/// returns false if SF_weight <= 0 and not force_sigcc (l.1436); main%sigcc/sigcc_recon (phsp: 1; H, D, A(e,e'p):
/// the model), the Coulomb factor (l.1533-1535), main%weight (l.1557-1560, tgtweight = `tgtweight`).
/// Not used by the meson reactions (pion, kaon, Delta, rho, SIDIS), which assemble their weights themselves.
/// @param force_sigcc compute sigcc even when SF_weight <= 0 (central event).
/// @param main the main event record (sigcc, sigcc_recon, coulomb, weight written).
/// @param vertex the generated vertex.
/// @param vertex0 the central vertex.
/// @param recon the reconstructed event.
/// @param model the cross-section model.
/// @param ctx the run constants.
/// @param tgtweight target weight factor.
/// @return false if SF_weight <= 0 and not @p force_sigcc.
bool assembleWeight(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0, const Event& recon,
                    IPhysicsModel& model, const MainContext& ctx, double tgtweight = 1.0);

/// Fortran: structures.inc rangetype (lo, hi); min_max_init starts lo = +1e10, hi = -1e10.
struct Range {
    double lo{1.0e10}, hi{-1.0e10};  ///< lowest and highest value seen.
};
/// Ranges of one arm.
struct ArmRange {
    Range delta, yptar, xptar, z;  ///< delta [%], yptar, xptar [rad], z [cm].
};
/// Fortran: structures.inc contribtype (init.f:min_max_init l.1080-1190).
struct Contrib {
    struct { ArmRange e, p; Range Trec, sumEgen; } gen;  ///< generated values.
    struct { struct { Range E, yptar, xptar; } e, p; Range Em, Pm, Trec; } tru;  ///< true values (after radiation).
    struct { ArmRange e, p; } SP;  ///< spectrometer (SP) values.
    struct { Range Trec, Em, Pm; } vertex;  ///< vertex values.
    struct { std::array<Range, 3> Egamma; Range Egamma_total; } rad;   ///< Egamma(1:3) -> [0..2]
};

/// Fortran: event.f:update_range (l.95-108).
/// @param val the value.
/// @param[in,out] range the range to extend.
void updateRange(double val, Range& range);
/// Fortran: event.f:update_slop_range (l.112-125).
/// @param val the value.
/// @param[in,out] slop the slop range to extend.
void updateSlopRange(double val, SlopItem& slop);
/// Fortran: event.f:limits_update (l.1-91). Egamma_used = IRadiator::photonEnergies() (Fortran Egamma_used(1:3) at
/// [0..2]).
/// rho_gen_from_orig: D-057 F18 (doing_rho, corrected mode): the GENERATION P-arm values from orig%p (the decay pion).
/// @param main the main event record.
/// @param vertex the generated vertex.
/// @param orig the event after radiation.
/// @param recon the reconstructed event.
/// @param f the reaction flags.
/// @param Egamma_used the radiated photon energies [MeV].
/// @param[in,out] contrib the ranges to update.
/// @param[in,out] slop the slop ranges to update.
/// @param rho_gen_from_orig D-057 F18 (see above).
void limitsUpdate(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                  const ReactionFlags& f, const std::array<double, 3>& Egamma_used, Contrib& contrib, Slops& slop,
                  bool rho_gen_from_orig = false);

/// P15.4 (parallel event loop): combine the ranges found by one worker into the run's: lo = min, hi = max, with the
/// update_range comparisons (a worker that never updated a range, lo = +1e10 / hi = -1e10, changes nothing).
/// @param[in,out] into the run's range.
/// @param worker the worker's range.
void mergeRange(Range& into, const Range& worker);
/// @param[in,out] into the run's slop range.
/// @param worker the worker's slop range.
void mergeSlopRange(SlopItem& into, const SlopItem& worker);   ///< lo/hi as mergeRange; `used` is the set-up value
/// Merge all ranges of @p worker into @p into (mergeRange on every member).
/// @param[in,out] into the run's ranges.
/// @param worker the worker's ranges.
void mergeContrib(Contrib& into, const Contrib& worker);
/// Merge all slop ranges of @p worker into @p into (mergeSlopRange on every member).
/// @param[in,out] into the run's slop ranges.
/// @param worker the worker's slop ranges.
void mergeSlops(Slops& into, const Slops& worker);

/// Outputs of jacobians.f:transform_to_cm.
struct CmTransform {
    double gstar{0}, bstar{0}, bstarx{0}, bstary{0}, bstarz{0};  ///< gamma, beta and beta components of the boost to the CM frame.
    double nustar{0}, qstar{0}, qstarx{0}, qstary{0}, qstarz{0};  ///< energy, momentum and momentum components of the virtual photon in the CM frame [MeV].
    double ehadcm{0}, phadcm{0}, phadcmx{0}, phadcmy{0}, phadcmz{0};  ///< energy, momentum and momentum components of the hadron in the CM frame [MeV].
    double ebeamcm{0}, pbeamcm{0}, pbeamcmx{0}, pbeamcmy{0}, pbeamcmz{0};  ///< energy, momentum and momentum components of the beam in the CM frame [MeV].
    double etarcm{0}, ptarcm{0}, ptarcmx{0}, ptarcmy{0}, ptarcmz{0};  ///< energy, momentum and momentum components of the target nucleon in the CM frame [MeV].
    double thetacm{0}, phicm{0}, phiqn{0}, jacobian{0}, jac_old{0};  ///< CM angles [rad], and the CM-to-lab Jacobian (new and old formula).
};
/// Fortran: jacobians.f:transform_to_cm (l.1-282). Reads vertex Ein, nu, q, p%P, p%E, up, uq, main%phi_pq and
/// COMMON /pfermi/ (pfer, pferx..z, efer). The |beta|>1 warning goes to `log` (if non-null).
/// @param vertex the generated vertex.
/// @param main the main event record.
/// @param pfer,pferx,pfery,pferz Fermi momentum and its components [MeV].
/// @param efer energy of the struck nucleon [MeV].
/// @param log stream for the warning (not owned), or nullptr.
/// @return the CM quantities.
CmTransform transformToCm(const Event& vertex, const MainEvent& main, double pfer, double pferx, double pfery,
                          double pferz, double efer, std::ostream* log = nullptr);

/// @}

}  // namespace reactions
}  // namespace simc

#endif  // SIMC_REACTIONS_EVENTCHAIN_H
