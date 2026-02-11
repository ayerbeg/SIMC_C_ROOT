// include/simc/Bremsstrahlung.h
// Exact QED bremsstrahlung calculations

#ifndef SIMC_BREMSSTRAHLUNG_H
#define SIMC_BREMSSTRAHLUNG_H

#include "simc/core/SimcConstants.h"
#include "simc/physics/MathUtils.h"
#include <array>

namespace simc {

/**
 * @struct BremResult
 * @brief Result of bremsstrahlung calculation
 */
struct BremResult {
    double bsoft = 0.0;        ///< Soft photon correction δ_soft
    double bhard = 0.0;        ///< Hard photon correction δ_hard
    double dbsoft = 0.0;       ///< Derivative dδ_soft/dE_γ
    double rad_factor = 1.0;   ///< Radiative cross section factor
};

/**
 * @class Bremsstrahlung
 * @brief Exact QED bremsstrahlung calculations
 * 
 * Implements exact soft photon bremsstrahlung calculations following
 * Tsai and Maximon-Tjon prescriptions. Handles both on-shell and 
 * off-shell kinematics.
 * 
 * References:
 * - Tsai, Phys. Rev. 122, 1898 (1961)
 * - Maximon & Tjon, Phys. Rev. C 62, 054320 (2000)
 */
class Bremsstrahlung {
public:
    /**
     * @brief Calculate bremsstrahlung with on-shell protons
     * @param ein Incident electron energy (MeV)
     * @param eout Scattered electron energy (MeV)
     * @param egamma Photon energy (MeV)
     * @param radiate_proton Include proton radiation?
     * @param exponentiate Use exp(-δ)?
     * @param include_hard Include hard corrections?
     * @param calculate_spence Calculate full Spence function?
     * @return BremResult structure
     */
    static BremResult Calculate(double ein, double eout, double egamma,
                                bool radiate_proton,
                                bool exponentiate = true,
                                bool include_hard = true,
                                bool calculate_spence = false);
    
    /**
     * @brief Calculate bremsstrahlung with off-shell protons
     * @param egamma Photon energy (MeV)
     * @param k_i Initial electron 3-momentum [px, py, pz] (MeV/c)
     * @param k_f Final electron 3-momentum [px, py, pz] (MeV/c)
     * @param p_i Initial proton 4-momentum [E, px, py, pz] (MeV)
     * @param p_f Final proton 4-momentum [E, px, py, pz] (MeV)
     * @param radiate_proton Include proton radiation?
     * @param exponentiate Use exp(-δ)?
     * @param include_hard Include hard corrections?
     * @param calculate_spence Calculate full Spence function?
     * @return BremResult structure
     */
    static BremResult CalculateOffShell(double egamma,
                                       const std::array<double, 3>& k_i,
                                       const std::array<double, 3>& k_f,
                                       const std::array<double, 4>& p_i,
                                       const std::array<double, 4>& p_f,
                                       bool radiate_proton,
                                       bool exponentiate = true,
                                       bool include_hard = true,
                                       bool calculate_spence = false);

private:
    /**
     * @brief Calculate interference integral
     * @param calculate_spence Include Spence function terms?
     * @param alpha Kinematic parameter α = 2m₁² - 2·(E₁E₂ - p₁·p₂)
     * @param ar1 Root 1 from quadratic equation
     * @param ar2 Root 2 from quadratic equation
     * @param e1 Energy 1 (GeV)
     * @param e2 Energy 2 (GeV)
     * @param de Photon energy (GeV)
     * @return Integral value
     */
    static double Inter(bool calculate_spence, double alpha,
                       double ar1, double ar2,
                       double e1, double e2, double de);
    
    /**
     * @brief Derivative of interference integral with respect to de
     * @param alpha Kinematic parameter
     * @param ar1 Root 1
     * @param ar2 Root 2
     * @param de Photon energy (GeV)
     * @return Derivative dI/d(de)
     */
    static double InterPrime(double alpha, double ar1, double ar2, double de);
    
    /**
     * @brief Ultra-relativistic limit calculation
     * @param ein Incident energy (GeV)
     * @param eout Scattered energy (GeV)
     * @param de Photon energy (GeV)
     * @param theta Scattering angle (rad)
     * @return Soft photon correction in ultra-relativistic limit
     * 
     * This is used as a check/approximation for the full calculation.
     */
    static double Srad(double ein, double eout, double de, double theta);
    
    /**
     * @brief Vacuum polarization correction
     * @param Q2 Four-momentum transfer squared (GeV²)
     * @return Vacuum polarization correction factor
     */
    static double VacuumPolarization(double Q2);
    
    /**
     * @brief Calculate roots for interference integrals
     * @param alpha Kinematic parameter
     * @param adot Dot product term
     * @param m1_sq Mass₁² (GeV²)
     * @param m2_sq Mass₂² (GeV²)
     * @param ar1 Root 1 (OUTPUT)
     * @param ar2 Root 2 (OUTPUT)
     * @return True if roots are real
     */
    static bool CalculateRoots(double alpha, double adot,
                              double m1_sq, double m2_sq,
                              double& ar1, double& ar2);
};

} // namespace simc

#endif // SIMC_BREMSSTRAHLUNG_H
