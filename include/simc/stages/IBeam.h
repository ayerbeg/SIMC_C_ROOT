/**
 * @file
 * @ingroup stages
 * @brief simc::IBeam: the beam stage (incident energy and transverse vertex position).
 * @note Origin: Fortran SIMC event.f:generate, beam/vertex part (l.150-210).
 */
#ifndef SIMC_STAGES_IBEAM_H
#define SIMC_STAGES_IBEAM_H

#include "simc/stages/IStage.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief Generates the beam part of a try: the incident energy (with its spread) and the transverse vertex (beam
 *        width, raster, offsets).
 *
 * The energy shifts from Coulomb attraction and energy loss in the target are applied afterwards by the event chain.
 * Implementation: simc::StandardBeam ("standard").
 */
class IBeam : public IStage {
public:
    /**
     * @brief Generates the beam quantities of one try.
     * @param main Receives target.x, target.y, target.rasterx, target.rastery [cm] (and the beam energy shifts).
     * @param vertex Receives the nominal incident energy Ein [MeV].
     * @param rng The random source of the run (its call order is part of the defined behaviour).
     */
    virtual void generate(MainEvent& main, Event& vertex, IRandom& rng) = 0;
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_IBEAM_H
