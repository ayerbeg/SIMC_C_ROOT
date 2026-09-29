// include/simc/stages/IMaterialModel.h
// Energy loss and multiple scattering in a material slab.
// Fortran: enerloss_new.f, shared/musc.f, shared/musc_ext.f.
#ifndef SIMC_STAGES_IMATERIALMODEL_H
#define SIMC_STAGES_IMATERIALMODEL_H

#include "simc/stages/IStage.h"

namespace simc {

class IMaterialModel : public IStage {
public:
    /// enerloss_new(len, dens, zeff, aeff, epart, mpart, typeflag) -> Eloss (MeV).
    virtual double energyLoss(double len_cm, double density, double zeff, double aeff, double epart,
                              double mpart, int typeflag, IRandom& rng) = 0;
    /// musc(m2, p, rad_len, dth, dph): adds scattering angles (rad) to dth/dph.
    virtual void musc(double m2, double p, double rad_len, double& dth, double& dph, IRandom& rng) = 0;
    /// musc_ext(m2, p, rad_len, x_len, dph, dth, y, x): extended-slab version (also shifts x, y).
    virtual void muscExt(double m2, double p, double rad_len, double x_len, double& dph, double& dth,
                         double& y, double& x, IRandom& rng) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_IMATERIALMODEL_H
