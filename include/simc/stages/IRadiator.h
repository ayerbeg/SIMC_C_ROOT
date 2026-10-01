// include/simc/stages/IRadiator.h
// Radiative corrections.
// Fortran: init.f:radc_init (l.576), radc_init_ev (l.655), radc.f:generate_rad (l.120),
//          peaked_rad_weight (l.523).
#ifndef SIMC_STAGES_IRADIATOR_H
#define SIMC_STAGES_IRADIATOR_H

#include <array>
#include <functional>

#include "simc/stages/IStage.h"

namespace simc {

/// Fortran: event.f:complete_ev(main, vertex, success), as generate_rad calls it after the incoming-electron
/// photon has reduced vertex%Ein (radc.f:322-328). The callee recomputes the vertex kinematics, the target
/// quantities teff(2:3)/Eloss(2:3) (trip_thru_target: gauss1 draws from `rng`) and must call
/// IRadiator::initEvent(main, vertex) on the SAME radiator (event.f:1045 radc_init_ev), exactly like the
/// Fortran. Returns success. (P5.1: the reaction stage that owns complete_ev arrives in Phase 4.)
using CompleteEvFn = std::function<bool(MainEvent& main, Event& vertex, IRandom& rng)>;

/// The radiation quantities of the current event that the driver reports (Fortran radc.inc /radccom/, copied into
/// central%rad by calculate_central, simc.f:1268-1281, and printed in the .hist "Central Values"). The three tails
/// (1 incoming electron, 2 scattered electron, 3 hadron) are indexed 0..2; c_int/c_ext/c/g keep the Fortran (0:3).
struct RadiationValues {
    double hardcorfac{0}, etatzai{0}, g_int{0}, g_ext{0};
    std::array<double, 3> frac{};     ///< frac(1:3)
    std::array<double, 3> lambda{};   ///< lambda(1:3)
    std::array<double, 2> bt{};       ///< bt(1:2)
    std::array<double, 4> c_int{}, c_ext{}, c{}, g{};   ///< (0:3)
};

class IRadiator : public IStage {
public:
    /// Per-event initialisation from the vertex kinematics (radc_init_ev).
    virtual void initEvent(MainEvent& main, const Event& vertex) = 0;
    /// Radiate (generate_rad): vertex -> orig, multiplying main%gen_weight. False if no valid radiation.
    /// vertex is modified like in the Fortran (vertex%Ein reduced by the tail-1 photon, recompleted through
    /// complete_ev); the caller keeps using it as the vertex event.
    virtual bool radiate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng,
                         const CompleteEvFn& complete_ev) = 0;

    // ---- what the driver reads back (D-022 revisited in P7.4; defaults suit a radiator without these notions) ----
    /// Values after the last initEvent() (the driver takes them right after the central event's initEvent).
    virtual RadiationValues currentValues() const { return {}; }
    /// Photon energy (MeV) radiated in each tail by the last successful radiate() (Fortran Egamma_used(1:3), read by
    /// limits_update); tails 0..2 = Fortran 1..3.
    virtual std::array<double, 3> photonEnergies() const { return {}; }
    /// The hadron radiates in the current event (Fortran rad_proton_this_ev, set by radc_init_ev): such successes
    /// are not counted in ncontribute_no_rad_proton (simc.f:296).
    virtual bool hadronRadiates() const { return false; }
};

}  // namespace simc

#endif  // SIMC_STAGES_IRADIATOR_H
