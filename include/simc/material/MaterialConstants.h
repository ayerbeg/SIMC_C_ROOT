/**
 * @file
 * @ingroup material
 * @brief Window and wall materials and geometry constants of the target (target.inc PARAMETERs).
 *
 * Window / wall materials and geometry constants. Fortran: target.inc (PARAMETERs, l.1-36).
 * All literals are doubles (-fdefault-real-8); the X0_cm values are the same double divisions
 * the Fortran compiler folds at compile time.
 */
#ifndef SIMC_MATERIAL_MATERIALCONSTANTS_H
#define SIMC_MATERIAL_MATERIALCONSTANTS_H

namespace simc {
namespace material {

/// @addtogroup material
/// @{

// Aluminium
constexpr double rho_Al = 2.70;  ///< density of aluminium [g/cm^3].
constexpr double Z_Al = 13.;  ///< effective Z of aluminium.
constexpr double A_Al = 26.98;  ///< effective A of aluminium.
constexpr double X0_Al = 24.01;  ///< radiation length of aluminium [g/cm^2].
constexpr double X0_cm_Al = X0_Al / rho_Al;  ///< radiation length of aluminium [cm].
// Mylar
constexpr double rho_mylar = 1.39;  ///< density of mylar [g/cm^3].
constexpr double Z_mylar = 4.545;  ///< effective Z of mylar.
constexpr double A_mylar = 8.735;  ///< effective A of mylar.
constexpr double X0_mylar = 39.95;  ///< radiation length of mylar [g/cm^2].
constexpr double X0_cm_mylar = X0_mylar / rho_mylar;  ///< radiation length of mylar [cm].
// Kevlar
constexpr double rho_kevlar = 0.74;  ///< density of kevlar [g/cm^3].
constexpr double Z_kevlar = 2.67;  ///< effective Z of kevlar.
constexpr double A_kevlar = 4.67;  ///< effective A of kevlar.
constexpr double X0_kevlar = 55.2;  ///< radiation length of kevlar [g/cm^2].
constexpr double X0_cm_kevlar = X0_kevlar / rho_kevlar;  ///< radiation length of kevlar [cm].
// Air (80% N2, 20% O2)
constexpr double rho_air = 0.00121;  ///< density of air [g/cm^3].
constexpr double Z_air = 7.2;  ///< effective Z of air.
constexpr double A_air = 14.4;  ///< effective A of air.
constexpr double X0_air = 36.66;  ///< radiation length of air [g/cm^2].
constexpr double X0_cm_air = X0_air / rho_air;  ///< radiation length of air [cm].

constexpr double inch_cm = 2.54;  ///< centimetres per inch.
/// target.inc: "stupid name to avoid conflict with definition in constants.inc" — note it is
/// 3.14159265358979, NOT constants.inc pi = 3.141592653589793.
constexpr double target_pi = 3.14159265358979;

/// enerloss_new.f: parameter(me=0.51099906)
constexpr double me_eloss = 0.51099906;

/// @}

}  // namespace material
}  // namespace simc

#endif  // SIMC_MATERIAL_MATERIALCONSTANTS_H
