/**
 * @file
 * @ingroup stages
 * @brief simc::ITarget: the target stage (geometry, vertex z, material traversed by each particle).
 * @note Origin: Fortran SIMC target.f:trip_thru_target (l.1-306), extreme_trip_thru_target (l.310-544), target_musc
 *       (l.548-589); event.f:generate z vertex (l.185).
 */
#ifndef SIMC_STAGES_ITARGET_H
#define SIMC_STAGES_ITARGET_H

#include "simc/stages/IStage.h"
#include "simc/stages/ITargetLossProvider.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief A target model: generates the interaction point and computes the material each particle crosses.
 *
 * Besides the per-event calls it provides the start-up losses of the set-up (simc::ITargetLossProvider), so that
 * simc::BuildRunSetup uses the same model as the event loop. Implementation: simc::material::HallCCryoTarget
 * ("hallc_cryo"). Not thread-safe (the energy loss may keep state).
 */
class ITarget : public IStage, public ITargetLossProvider {
public:
    /// The particle whose path is computed (codes of the original trip_thru_target).
    enum Arm {
        kBeam = 1,       ///< The incoming beam electron (vertex upstream).
        kElectron = 2,   ///< The scattered electron.
        kHadron = 3      ///< The hadron.
    };
    /// Energy-loss mode (the typeflag of the energy-loss formula).
    enum LossMode {
        kSample = 1,       ///< Sampled (random straggling).
        kMin = 2,          ///< Minimum loss.
        kMax = 3,          ///< Maximum loss.
        kMostProbable = 4  ///< Most probable loss.
    };

    /// Material traversed along one path.
    struct Traversal {
        double eloss{0};    ///< Energy loss [MeV].
        double radlen{0};   ///< Thickness [radiation lengths].
    };
    /// Multiple-scattering angle kicks added to the target slopes (original target_musc dangles(1), dangles(2)).
    struct ScatteringAngles {
        double dyptar{0};   ///< Kick of yptar [rad].
        double dxptar{0};   ///< Kick of xptar [rad].
    };

    /**
     * @brief Generates the z position of the interaction.
     * @param rng The random source.
     * @return z along the beam [cm], including the target zoffset.
     */
    virtual double generateZ(IRandom& rng) = 0;
    /**
     * @brief The material crossed by one particle from the vertex out of (or, for the beam, into) the target.
     * @param arm Which particle.
     * @param zpos Vertex z relative to the target centre [cm].
     * @param energy Total energy of the particle [MeV].
     * @param theta Lab angle of the particle [rad].
     * @param mass Mass of the particle [MeV].
     * @param mode Energy-loss mode.
     * @param rng The random source (used for kSample).
     * @return Energy loss and thickness along the path.
     * @note Origin: Fortran SIMC target.f:trip_thru_target(narm, zpos, energy, theta, Eloss, radlen, mass, typeflag).
     */
    virtual Traversal tripThruTarget(Arm arm, double zpos, double energy, double theta, double mass, LossMode mode,
                                     IRandom& rng) = 0;
    /**
     * @brief Multiple-scattering kicks of a particle leaving the target.
     * @param p Momentum [MeV/c].
     * @param beta Velocity v/c.
     * @param teff Thickness crossed [radiation lengths].
     * @param rng The random source.
     * @return The angle kicks.
     * @note Origin: Fortran SIMC target.f:target_musc(p, beta, teff, dangles).
     */
    virtual ScatteringAngles multipleScattering(double p, double beta, double teff, IRandom& rng) = 0;
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_ITARGET_H
