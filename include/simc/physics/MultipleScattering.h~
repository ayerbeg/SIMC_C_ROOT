// include/simc/MultipleScattering.h
// Multiple scattering calculations for charged particles through matter
// Ported from musc.f and musc_ext.f

#ifndef SIMC_MULTIPLE_SCATTERING_H
#define SIMC_MULTIPLE_SCATTERING_H

#include "SimcTypes.h"
#include "EnergyLoss.h"
#include "RandomGenerator.h"
#include <cmath>

namespace simc {

/**
 * @class MultipleScattering
 * @brief Calculate angular deflection due to Coulomb scattering
 * 
 * Implements multiple scattering using:
 * - Highland formula for RMS angle (small angles)
 * - Molière theory for full distribution (includes non-Gaussian tails)
 * - Lynch-Dahl corrections
 * 
 * Based on musc.f and musc_ext.f from original SIMC.
 * 
 * Physics References:
 * - W.T. Scott, Rev. Mod. Phys. 35, 231 (1963) - Molière theory
 * - V.L. Highland, Nucl. Instr. Meth. 129, 497 (1975) - Highland formula
 * - G.R. Lynch and O.I. Dahl, Nucl. Instr. Meth. B58, 6 (1991) - Corrections
 * - Particle Data Group, Section 33.3 "Passage of particles through matter"
 * 
 * Units:
 * - All angles in radians
 * - Momentum in MeV/c
 * - Length in cm
 * - Radiation length X0 in g/cm²
 */
class MultipleScattering {
public:
    
    /**
     * @struct ScatterAngles
     * @brief Result of multiple scattering calculation
     */
    struct ScatterAngles {
        double theta_x{0.0};    ///< Scattering angle in x (horizontal) [rad]
        double theta_y{0.0};    ///< Scattering angle in y (vertical) [rad]
        double theta_rms{0.0};  ///< RMS scattering angle [rad]
    };
    
    /**
     * @brief Calculate multiple scattering angles
     * 
     * @param length Path length through material (cm)
     * @param material Material properties
     * @param momentum Particle momentum (MeV/c)
     * @param charge Particle charge (in units of e)
     * @param mass Particle mass (MeV/c²)
     * @param rng Random number generator
     * @param use_moliere Use full Molière theory (true) or Gaussian approx (false)
     * @return ScatterAngles structure with theta_x, theta_y, theta_rms
     * 
     * Example:
     * @code
     * Material lh2 = Material::LiquidHydrogen();
     * RandomGenerator rng;
     * auto angles = MultipleScattering::Calculate(
     *     10.0,      // 10 cm path
     *     lh2,       // liquid hydrogen
     *     8000.0,    // 8 GeV/c
     *     1.0,       // electron charge
     *     0.511,     // electron mass
     *     rng,
     *     false      // Gaussian approximation
     * );
     * @endcode
     */
    static ScatterAngles Calculate(double length, 
                                   const Material& material,
                                   double momentum, 
                                   double charge, 
                                   double mass,
                                   RandomGenerator& rng,
                                   bool use_moliere = false);
    
    /**
     * @brief Calculate RMS scattering angle using Highland formula
     * 
     * Highland formula with Lynch-Dahl corrections:
     * theta_rms = (13.6 MeV/βcp) * z * sqrt(x/X0) * [1 + 0.038*ln(x/X0)]
     * 
     * where:
     * - βc = velocity
     * - p = momentum
     * - z = charge
     * - x/X0 = thickness in radiation lengths
     * 
     * @param thickness Path length in radiation lengths (x/X0)
     * @param momentum Particle momentum (MeV/c)
     * @param charge Particle charge (units of e)
     * @param mass Particle mass (MeV/c²)
     * @return RMS scattering angle in radians
     * 
     * Reference: Highland (1975), Lynch & Dahl (1991)
     */
    static double CalculateRMS(double thickness,
                              double momentum,
                              double charge,
                              double mass);
    
    /**
     * @brief Sample scattering angles from Gaussian distribution
     * 
     * Simple Gaussian approximation valid for:
     * - Small angles (< few degrees)
     * - Thick absorbers (x/X0 > 0.01)
     * 
     * @param theta_rms RMS scattering angle (rad)
     * @param rng Random number generator
     * @return ScatterAngles with Gaussian-sampled theta_x, theta_y
     */
    static ScatterAngles SampleGaussian(double theta_rms, RandomGenerator& rng);
    
    /**
     * @brief Sample scattering angles from Molière distribution
     * 
     * Full Molière theory including non-Gaussian tails.
     * More accurate for:
     * - Thin absorbers (x/X0 < 0.01)
     * - Large angles
     * - High precision requirements
     * 
     * @param thickness Path length in radiation lengths (x/X0)
     * @param momentum Particle momentum (MeV/c)
     * @param charge Particle charge (units of e)
     * @param mass Particle mass (MeV/c²)
     * @param rng Random number generator
     * @return ScatterAngles with Molière-sampled theta_x, theta_y
     * 
     * Reference: Scott (1963), PDG Section 33.3
     */
    static ScatterAngles SampleMoliere(double thickness,
                                      double momentum,
                                      double charge,
                                      double mass,
                                      RandomGenerator& rng);
    
    /**
     * @brief Calculate screening parameter chi_alpha
     * 
     * Used in Molière theory:
     * chi_alpha = 2.007e-5 * Z^(2/3) / (β² * (1 + 3.34*(Z*α/β)²))
     * 
     * @param Z Atomic number
     * @param beta Particle velocity (v/c)
     * @return Screening parameter chi_alpha
     */
    static double CalculateChiAlpha(double Z, double beta);
    
    /**
     * @brief Calculate Molière parameter B
     * 
     * B = ln(chi_c / chi_alpha) where chi_c ~ 1.13 + 3.76*alpha²
     * 
     * @param chi_alpha Screening parameter
     * @param Z Atomic number
     * @return Molière parameter B
     */
    static double CalculateMoliereB(double chi_alpha, double Z);
    
    /**
     * @brief Get radiation length in cm for a material
     * 
     * @param material Material properties
     * @return Radiation length in cm
     */
    static double GetRadiationLength(const Material& material);
    
    /**
     * @brief Calculate beta (v/c) from momentum and mass
     * 
     * @param momentum Particle momentum (MeV/c)
     * @param mass Particle mass (MeV/c²)
     * @return beta = v/c
     */
    static double CalculateBeta(double momentum, double mass);
    
    /**
     * @brief Calculate beta*gamma from momentum and mass
     * 
     * @param momentum Particle momentum (MeV/c)
     * @param mass Particle mass (MeV/c²)
     * @return beta*gamma
     */
    static double CalculateBetaGamma(double momentum, double mass);

private:
    /**
     * @brief Sample from Molière distribution (internal)
     * 
     * Uses inverse transform sampling from Molière's f(θ) distribution
     * 
     * @param B Molière parameter
     * @param chi_alpha Screening parameter
     * @param rng Random number generator
     * @return Scattering angle from Molière distribution (radians)
     */
    static double SampleMoliereAngle(double B, double chi_alpha, RandomGenerator& rng);
    
    /**
     * @brief Evaluate Molière cumulative distribution
     * 
     * @param theta Angle (radians)
     * @param B Molière parameter
     * @param chi_alpha Screening parameter
     * @return CDF value
     */
    static double MoliereCDF(double theta, double B, double chi_alpha);
    
    /**
     * @brief Calculate f0, f1, f2 functions for Molière theory
     * 
     * These are the Molière angular distribution functions
     * 
     * @param u Reduced angle variable
     * @param B Molière parameter
     * @return {f0, f1, f2} array
     */
    static std::array<double, 3> MoliereFunctions(double u, double B);
};

/**
 * @class MultipleScatteringTable
 * @brief Pre-computed multiple scattering tables for speed
 * 
 * For applications requiring many scattering calculations,
 * pre-compute and store tables of theta_rms vs. thickness
 */
class MultipleScatteringTable {
public:
    /**
     * @brief Constructor
     * @param material Material to tabulate
     * @param momentum_min Minimum momentum (MeV/c)
     * @param momentum_max Maximum momentum (MeV/c)
     * @param n_points Number of table points
     */
    MultipleScatteringTable(const Material& material,
                           double momentum_min,
                           double momentum_max,
                           int n_points = 100);
    
    /**
     * @brief Get RMS angle from table (interpolated)
     * @param momentum Particle momentum (MeV/c)
     * @param length Path length (cm)
     * @param charge Particle charge (units of e)
     * @param mass Particle mass (MeV/c²)
     * @return Interpolated theta_rms
     */
    double GetRMS(double momentum, double length, double charge, double mass) const;
    
private:
    Material material_;
    double momentum_min_;
    double momentum_max_;
    std::vector<double> momentum_points_;
    std::vector<double> rms_points_;
    
    double Interpolate(double momentum) const;
};

} // namespace simc

#endif // SIMC_MULTIPLE_SCATTERING_H
