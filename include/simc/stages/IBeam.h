// include/simc/stages/IBeam.h
// Beam: incident energy (spread, Coulomb/Eloss shifts applied by the driver) and the
// transverse vertex (beam width, raster, offsets).
// Fortran: event.f:generate, beam/vertex part (l.150-210).
#ifndef SIMC_STAGES_IBEAM_H
#define SIMC_STAGES_IBEAM_H

#include "simc/stages/IStage.h"

namespace simc {

class IBeam : public IStage {
public:
    /// Fill main%target%x/y/rasterx/rastery and the nominal incident energy of `vertex`.
    virtual void generate(MainEvent& main, Event& vertex, IRandom& rng) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_IBEAM_H
