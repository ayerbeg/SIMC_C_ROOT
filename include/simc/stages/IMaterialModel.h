/**
 * @file
 * @ingroup stages
 * @brief simc::IMaterialModel: energy loss and multiple scattering in a slab of material.
 * @note Origin: Fortran SIMC enerloss_new.f, shared/musc.f, shared/musc_ext.f.
 */
#ifndef SIMC_STAGES_IMATERIALMODEL_H
#define SIMC_STAGES_IMATERIALMODEL_H

#include "simc/stages/IStage.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief Ionisation energy loss (most probable loss with Landau straggling) and multiple scattering (Gaussian,
 *        Highland formula) of a particle crossing a slab of material.
 *
 * Implementation: simc::material::SimcMaterial ("simc_material"), the formulas of the original program. An implementation may keep
 * state between calls (the original keeps the Landau variable of the previous call); not thread-safe.
 */
class IMaterialModel : public IStage {
public:
    /**
     * @brief Energy lost by a particle crossing a slab.
     * @param len_cm Path length in the material [cm].
     * @param density Density of the material [g/cm^3].
     * @param zeff Effective atomic number.
     * @param aeff Effective mass number.
     * @param epart Total energy of the particle [MeV].
     * @param mpart Mass of the particle [MeV].
     * @param typeflag 1: sampled (one Gaussian draw as the Landau variable), 2: minimum, 3: maximum, 4: most probable
     *        energy loss. Any other value reuses the Landau variable of the previous call (FORTRAN-QUIRK D-012 F3).
     * @param rng Random source (used only for typeflag 1).
     * @return The energy loss [MeV] (>= 0).
     * @note Origin: Fortran SIMC enerloss_new(len, dens, zeff, aeff, epart, mpart, typeflag, Eloss).
     */
    virtual double energyLoss(double len_cm, double density, double zeff, double aeff, double epart,
                              double mpart, int typeflag, IRandom& rng) = 0;
    /**
     * @brief Multiple scattering in a thin slab: adds random scattering angles to the track direction.
     * @param m2 Mass squared of the particle [MeV^2].
     * @param p Momentum of the particle [MeV/c].
     * @param rad_len Thickness of the slab [radiation lengths].
     * @param dth In/out: first transverse angle [rad]; a Gaussian scattering angle is added.
     * @param dph In/out: second transverse angle [rad]; a Gaussian scattering angle is added.
     * @param rng Random source (two Gaussian draws; none if @p rad_len is 0).
     * @note Origin: Fortran SIMC shared/musc.f:musc(m2, p, rad_len, dth, dph).
     */
    virtual void musc(double m2, double p, double rad_len, double& dth, double& dph, IRandom& rng) = 0;
    /**
     * @brief Multiple scattering in a thick slab: adds the angles and the correlated transverse displacements.
     * @param m2 Mass squared of the particle [MeV^2].
     * @param p Momentum of the particle [MeV/c].
     * @param rad_len Thickness of the slab [radiation lengths].
     * @param x_len Length of the slab [cm].
     * @param dph In/out: horizontal angle (slope) [rad].
     * @param dth In/out: vertical angle (slope) [rad].
     * @param y In/out: horizontal position [cm].
     * @param x In/out: vertical position [cm].
     * @param rng Random source.
     * @note Origin: Fortran SIMC shared/musc_ext.f:musc_ext(m2, p, rad_len, x_len, dph, dth, y, x).
     */
    virtual void muscExt(double m2, double p, double rad_len, double x_len, double& dph, double& dth,
                         double& y, double& x, IRandom& rng) = 0;
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_IMATERIALMODEL_H
