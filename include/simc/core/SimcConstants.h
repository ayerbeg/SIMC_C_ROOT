// include/simc/SimcConstants.h
// Physical constants for SIMC Monte Carlo
// All energies in MeV, distances in cm, angles in radians

#ifndef SIMC_CONSTANTS_H
#define SIMC_CONSTANTS_H

namespace simc {
namespace constants {

// ============================================================================
// Particle Masses (MeV/c^2)
// ============================================================================
constexpr double Me     = 0.51099906;      ///< Electron mass
constexpr double Mp     = 938.27231;       ///< Proton mass
constexpr double Mn     = 939.56563;       ///< Neutron mass
constexpr double Mpi    = 139.57018;       ///< Charged pion mass
constexpr double Mpi0   = 134.9766;        ///< Neutral pion mass
constexpr double Mmu    = 105.6583755;     ///< Muon mass
constexpr double Mk     = 493.677;         ///< Kaon mass
constexpr double Mrho   = 769.3;           ///< Rho meson mass
constexpr double Md     = 1875.613;        ///< Deuteron mass
constexpr double Mlambda = 1115.68;        ///< Lambda hyperon mass
constexpr double Msigma0 = 1192.64;        ///< Sigma0 hyperon mass
constexpr double Msigma_minus = 1197.45;   ///< Sigma- hyperon mass
constexpr double MDelta = 1232.0;          ///< Delta resonance mass

// ============================================================================
// Squared Masses (MeV^2) - Pre-calculated for efficiency
// ============================================================================
constexpr double Me2    = Me * Me;
constexpr double Mp2    = Mp * Mp;
constexpr double Mn2    = Mn * Mn;
constexpr double Mpi2   = Mpi * Mpi;
constexpr double Mpi02  = Mpi0 * Mpi0;
constexpr double Mmu2   = Mmu * Mmu;
constexpr double Mk2    = Mk * Mk;
constexpr double Mrho2  = Mrho * Mrho;
constexpr double Md2    = Md * Md;

// ============================================================================
// Physical Constants
// ============================================================================
constexpr double amu    = 931.49432;       ///< Atomic mass unit (MeV/c^2)
constexpr double hbarc  = 197.327053;      ///< hbar*c (MeV*fm)
constexpr double alpha  = 1.0 / 137.0359895; ///< Fine structure constant
constexpr double alpi   = alpha / 3.141592653589793; ///< alpha/pi
constexpr double euler  = 0.577215665;     ///< Euler's constant

// ============================================================================
// Mathematical Constants
// ============================================================================
constexpr double pi     = 3.141592653589793;
constexpr double twopi  = 2.0 * pi;
constexpr double degrad = 180.0 / pi;      ///< Degrees per radian (rad to deg)

// ============================================================================
// Unit Conversions
// ============================================================================
constexpr double GEV_TO_MEV = 1000.0;      ///< GeV to MeV conversion
constexpr double MEV_TO_GEV = 0.001;       ///< MeV to GeV conversion
constexpr double RAD_TO_MR  = 1000.0;      ///< Radians to milliradians
constexpr double MR_TO_RAD  = 0.001;       ///< Milliradians to radians
constexpr double RAD_TO_DEG = degrad;      ///< Radians to degrees
constexpr double DEG_TO_RAD = pi / 180.0;  ///< Degrees to radians

// ============================================================================
// Resonance Widths (MeV)
// ============================================================================
constexpr double Delta_width = 117.0;      ///< Delta resonance width

} // namespace constants
} // namespace simc

#endif // SIMC_CONSTANTS_H
