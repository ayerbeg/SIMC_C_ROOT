// include/simc/stages/ITarget.h
// Target: geometry, z-vertex, material traversed by each particle.
// Fortran: target.f:trip_thru_target (l.1-306), extreme_trip_thru_target (l.310-544),
//          target_musc (l.548-589); event.f:generate z vertex (l.185).
#ifndef SIMC_STAGES_ITARGET_H
#define SIMC_STAGES_ITARGET_H

#include "simc/stages/IStage.h"

namespace simc {

class ITarget : public IStage {
public:
    /// Arm codes as in trip_thru_target: 1 = incoming beam, 2 = scattered electron, 3 = hadron.
    enum Arm { kBeam = 1, kElectron = 2, kHadron = 3 };
    /// Energy-loss mode as in enerloss_new typeflag: 1 sample, 2 min, 3 max, 4 most probable.
    enum LossMode { kSample = 1, kMin = 2, kMax = 3, kMostProbable = 4 };

    /// Generate the z position of the interaction (cm, including zoffset).
    virtual double generateZ(IRandom& rng) = 0;
    /// Energy loss (MeV) and thickness in radiation lengths along the path of `arm`.
    virtual void tripThruTarget(int arm, double zpos, double energy, double theta, double mass,
                                LossMode mode, IRandom& rng, double& eloss, double& radlen) = 0;
    /// Multiple-scattering angle offsets dangles[0..1] (rad) for momentum p (MeV/c).
    virtual void multipleScattering(double p, double beta, double teff, double dangles[2],
                                    IRandom& rng) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_ITARGET_H
