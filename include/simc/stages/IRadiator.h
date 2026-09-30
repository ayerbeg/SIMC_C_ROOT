// include/simc/stages/IRadiator.h
// Radiative corrections.
// Fortran: init.f:radc_init (l.576), radc_init_ev (l.655), radc.f:generate_rad (l.120),
//          peaked_rad_weight (l.523).
#ifndef SIMC_STAGES_IRADIATOR_H
#define SIMC_STAGES_IRADIATOR_H

#include <functional>

#include "simc/stages/IStage.h"

namespace simc {

/// Fortran: event.f:complete_ev(main, vertex, success), as generate_rad calls it after the incoming-electron
/// photon has reduced vertex%Ein (radc.f:322-328). The callee recomputes the vertex kinematics, the target
/// quantities teff(2:3)/Eloss(2:3) (trip_thru_target: gauss1 draws from `rng`) and must call
/// IRadiator::initEvent(main, vertex) on the SAME radiator (event.f:1045 radc_init_ev), exactly like the
/// Fortran. Returns success. (P5.1: the reaction stage that owns complete_ev arrives in Phase 4.)
using CompleteEvFn = std::function<bool(MainEvent& main, Event& vertex, IRandom& rng)>;

class IRadiator : public IStage {
public:
    /// Per-event initialisation from the vertex kinematics (radc_init_ev).
    virtual void initEvent(MainEvent& main, const Event& vertex) = 0;
    /// Radiate (generate_rad): vertex -> orig, multiplying main%gen_weight. False if no valid radiation.
    /// vertex is modified like in the Fortran (vertex%Ein reduced by the tail-1 photon, recompleted through
    /// complete_ev); the caller keeps using it as the vertex event.
    virtual bool radiate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng,
                         const CompleteEvFn& complete_ev) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_IRADIATOR_H
