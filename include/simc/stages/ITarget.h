// include/simc/stages/ITarget.h
// Target: geometry, z-vertex, material traversed by each particle.
// Fortran: target.f:trip_thru_target (l.1-306), extreme_trip_thru_target (l.310-544),
//          target_musc (l.548-589); event.f:generate z vertex (l.185).
#ifndef SIMC_STAGES_ITARGET_H
#define SIMC_STAGES_ITARGET_H

#include "simc/stages/IStage.h"
#include "simc/stages/ITargetLossProvider.h"

namespace simc {

/// A target model. Besides the per-event calls below it provides the start-up losses of the set-up
/// (ITargetLossProvider: init.f:target_init, target.f:extreme_trip_thru_target), so that BuildRunSetup uses the same
/// model as the event loop.
class ITarget : public IStage, public ITargetLossProvider {
public:
    /// Arm codes as in trip_thru_target: 1 = incoming beam, 2 = scattered electron, 3 = hadron.
    enum Arm { kBeam = 1, kElectron = 2, kHadron = 3 };
    /// Energy-loss mode as in enerloss_new typeflag: 1 sample, 2 min, 3 max, 4 most probable.
    enum LossMode { kSample = 1, kMin = 2, kMax = 3, kMostProbable = 4 };

    /// Material traversed along one path: energy loss (MeV) and thickness (radiation lengths).
    struct Traversal {
        double eloss{0};
        double radlen{0};
    };
    /// Multiple-scattering angle kicks (rad) added to the slopes (Fortran target_musc dangles(1), dangles(2)).
    struct ScatteringAngles {
        double dyptar{0};
        double dxptar{0};
    };

    /// Generate the z position of the interaction (cm, including zoffset).
    virtual double generateZ(IRandom& rng) = 0;
    /// Path of `arm` from the vertex at `zpos` (cm, relative to the target centre) with total energy `energy` (MeV),
    /// angle `theta` (rad) and mass `mass` (MeV). Fortran: trip_thru_target(narm, zpos, energy, theta, Eloss, radlen,
    /// mass, typeflag).
    virtual Traversal tripThruTarget(Arm arm, double zpos, double energy, double theta, double mass, LossMode mode,
                                     IRandom& rng) = 0;
    /// Multiple-scattering kicks for momentum p (MeV/c), velocity beta and thickness teff (radiation lengths).
    /// Fortran: target_musc(p, beta, teff, dangles).
    virtual ScatteringAngles multipleScattering(double p, double beta, double teff, IRandom& rng) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_ITARGET_H
