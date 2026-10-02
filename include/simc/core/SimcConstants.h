// include/simc/core/SimcConstants.h
// Physical constants of SIMC (constants.inc). Units: simc/core/Units.h (MeV, cm, rad).

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

// Unit conversions: see simc/core/Units.h (the single unit convention of SIMC++).

// ============================================================================
// Resonance Widths (MeV)
// ============================================================================
constexpr double Delta_width = 117.0;      ///< Delta resonance width

} // namespace constants
} // namespace simc

#endif // SIMC_CONSTANTS_H
