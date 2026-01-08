// include/simc/EnergyLoss.h
// Energy loss calculations for charged particles through matter
// Ported from enerloss_new.f

#ifndef SIMC_ENERGY_LOSS_H
#define SIMC_ENERGY_LOSS_H

#include "SimcTypes.h"
#include "RandomGenerator.h"

namespace simc {

/**
 * @class EnergyLoss
 * @brief Calculate ionization energy loss through matter
 * 
 * This class implements the Bethe-Bloch formula with density effect
 * corrections and Landau fluctuations for calculating energy loss
 * of charged particles through matter.
 * 
 * Based on enerloss_new.f from Fortran SIMC.
 * 
 * Units:
 * - Energy: MeV
 * - Length: cm
 * - Density: g/cm³
 */
class EnergyLoss {
public:
  /**
     * @enum LossType
     * @brief Type of energy loss calculation
     * 
     * These correspond to different positions in the Landau distribution:
     * - SAMPLED: Random sampling from full Landau distribution
     * - MINIMUM: "Minimum ionizing" → lambda=3.0 (lower tail) → HIGHER energy loss
     * - MAXIMUM: lambda=0.0067 (upper tail) → LOWER energy loss
     * - MOST_PROBABLE: lambda=1.0 (peak of distribution) → INTERMEDIATE energy loss
     * 
     * Note: The naming is historical. "Minimum ionizing" refers to particles
     * that ionize infrequently but lose more energy per collision, resulting
     * in HIGHER total energy loss. The relationship is:
     * 
     *   MINIMUM > MOST_PROBABLE > MAXIMUM
     * 
     * This follows the Landau distribution physics where energy loss = lambda*xi + Eloss_mp
     */
  
    enum class LossType {
        SAMPLED = 1,      ///< Sample from Landau distribution
        MINIMUM = 2,      ///< Minimum energy loss
        MAXIMUM = 3,      ///< Maximum energy loss  
        MOST_PROBABLE = 4 ///< Most probable energy loss
    };
    
  /**
     * @brief Calculate energy loss through material
     * 
     * @param length Path length through material (cm)
     * @param density Material density (g/cm³)
     * @param Z Atomic number
     * @param A Atomic mass (amu)
     * @param E_particle Particle energy (MeV)
     * @param M_particle Particle mass (MeV/c²)
     * @param rng Random number generator (for LossType::SAMPLED)
     * @param type Type of energy loss calculation
     * @return Energy loss (MeV)
     * 
     * Example:
     * @code
     * RandomGenerator rng;
     * double loss = EnergyLoss::Calculate(
     *     10.0,      // 10 cm path
     *     0.0708,    // LH2 density
     *     1.0,       // Hydrogen Z
     *     1.008,     // Hydrogen A
     *     8000.0,    // 8 GeV electron
     *     0.511,     // Electron mass
     *     rng,
     *     EnergyLoss::LossType::SAMPLED
     * );
     * @endcode
     */
    static double Calculate(double length, double density, 
                           double Z, double A,
                           double E_particle, double M_particle,
                           RandomGenerator& rng,
                           LossType type = LossType::SAMPLED);
    
  /**
     * @brief Calculate most probable energy loss (Bethe-Bloch)
     * 
     * @param thickness Material thickness (g/cm²)
     * @param Z Atomic number
     * @param A Atomic mass (amu)
     * @param beta Particle velocity (v/c)
     * @param gamma Lorentz factor
     * @param density Material density (g/cm³)
     * @return Most probable energy loss (MeV)
     */
    static double MostProbableEnergyLoss(double thickness, 
                                        double Z, double A,
                                        double beta, double gamma,
                                        double density);
    
  /**
     * @brief Calculate Bethe-Bloch dE/dx
     * 
     * @param beta Particle velocity (v/c)
     * @param gamma Lorentz factor
     * @param Z Atomic number
     * @param A Atomic mass
     * @param I Mean excitation energy (MeV)
     * @param density Material density (g/cm³)
     * @return dE/dx in MeV/(g/cm²)
     */
    static double BetheBloch(double beta, double gamma,
                            double Z, double A, double I,
                            double density);
    
  /**
     * @brief Get mean excitation energy for element
     * 
     * @param Z Atomic number
     * @return Mean excitation energy I (MeV)
     * 
     * Uses formula:
     * - I = 21.8 eV for Hydrogen (Z=1)
     * - I = 16*Z^0.9 eV for other elements
     */
    static double GetIonizationPotential(double Z);
    
  /**
     * @brief Calculate density effect correction
     * 
     * @param beta Particle velocity (v/c)
     * @param gamma Lorentz factor
     * @param Z Atomic number
     * @param A Atomic mass
     * @param density Material density (g/cm³)
     * @return Density correction delta
     */
    static double DensityCorrection(double beta, double gamma,
                                   double Z, double A, double density);
    
  /**
     * @brief Calculate plasma frequency parameter
     * 
     * @param density Material density (g/cm³)
     * @param Z Atomic number
     * @param A Atomic mass
     * @return Plasma frequency hbar*omega_p (MeV)
     */
    static double PlasmaFrequency(double density, double Z, double A);
    
  /**
     * @brief Calculate xi parameter for Landau distribution
     * 
     * @param thickness Material thickness (g/cm²)
     * @param Z Atomic number
     * @param A Atomic mass
     * @param beta Particle velocity (v/c)
     * @return Xi parameter (MeV)
     */
    static double CalcXi(double thickness, double Z, double A, double beta);
    
  /**
     * @brief Sample from Landau distribution
     * 
     * @param rng Random number generator
     * @return Lambda parameter for Landau distribution
     * 
     * The Landau parameter lambda is sampled from:
     * lambda = -2*ln(x) where x is from Gaussian tail
     */
    static double SampleLandau(RandomGenerator& rng);

private:
  /**
     * @brief Sample from Gaussian tail (x > 0)
     * 
     * @param rng Random number generator
     * @param nsigma Number of standard deviations for cutoff
     * @return Random value from Gaussian tail
     */
    static double GaussianTail(RandomGenerator& rng, double nsigma = 10.0);
};

  /**
 * @class Material
 * @brief Properties of common materials
 * 
 * Provides standard material properties for energy loss calculations.
 */
class Material {
public:
    std::string name;
    double Z;           ///< Atomic number
    double A;           ///< Atomic mass (amu)
    double density;     ///< Density (g/cm³)
    double X0;          ///< Radiation length (g/cm²)
    
    // Standard materials
    static Material LiquidHydrogen();
    static Material LiquidDeuterium();
    static Material Helium3();
    static Material Helium4();
    static Material Carbon();
    static Material Aluminum();
    static Material Iron();
    static Material Copper();
    static Material Gold();
    
  /**
     * @brief Get material by name
     * @param name Material name (case-insensitive)
     * @return Material properties
     * @throws std::runtime_error if material not found
     */
    static Material GetByName(const std::string& name);
};

} // namespace simc

#endif // SIMC_ENERGY_LOSS_H
