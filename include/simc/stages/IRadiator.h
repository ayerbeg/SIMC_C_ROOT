// include/simc/stages/IRadiator.h
// Radiative corrections.
// Fortran: init.f:radc_init (l.576), radc_init_ev (l.655), radc.f:generate_rad (l.120),
//          peaked_rad_weight (l.523).
#ifndef SIMC_STAGES_IRADIATOR_H
#define SIMC_STAGES_IRADIATOR_H

#include "simc/stages/IStage.h"

namespace simc {

class IRadiator : public IStage {
public:
    /// Per-event initialisation from the vertex kinematics (radc_init_ev).
    virtual void initEvent(MainEvent& main, const Event& vertex) = 0;
    /// Radiate: vertex -> orig, multiplying main%gen_weight (generate_rad). False if no valid radiation.
    virtual bool radiate(MainEvent& main, const Event& vertex, Event& orig, IRandom& rng) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_IRADIATOR_H
