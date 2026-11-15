/**
 * @file
 * @ingroup core
 * @brief THE unit convention of SIMC++ and the named conversion factors.
 *
 * Everything inside the framework (simc::RunSetup, the events, the stage interfaces, the kernels) uses the internal
 * units of the original SIMC:
 *
 * | quantity | unit | quantity | unit | quantity | unit |
 * |---|---|---|---|---|---|
 * | energy, mass | MeV | momentum | MeV/c | mass^2 | MeV^2 |
 * | length, position | cm | angle | rad | slopes | dx/dz, dy/dz |
 * | target thickness | g/cm^2 | density | g/cm^3 | delta | dp/p in % |
 *
 * Cross sections are in the units the physics routines return (documented per model).
 *
 * Exceptions are the INPUT units of the `.inp` / JSON parameters (simc::SimcParameters: degrees, mg/cm^2, mrad, ...),
 * converted once by the set-up, and the OUTPUT units of the ntuple (GeV, mrad, ...), converted by the output code.
 * Both conversions use the named factors of simc::units. Kernels keep their own literals where the original writes
 * them (e.g. the bremsstrahlung code works in GeV with "/1000."), so that the ports stay line-by-line comparable.
 *
 * Every factor is exact and is used with the SAME operation as the original (x / 1000. stays a division by 1000), so
 * replacing a literal by its name does not change a single bit. Task P7.4.
 */
#ifndef SIMC_CORE_UNITS_H
#define SIMC_CORE_UNITS_H

#include "simc/core/SimcConstants.h"

namespace simc {
/**
 * @ingroup core
 * @brief Named, exact unit-conversion factors between the input / output units and the internal units (see the file
 *        description of core/Units.h).
 */
namespace units {

constexpr double kMeVPerGeV = 1000.0;     ///< x[GeV] = x[MeV] / kMeVPerGeV
constexpr double kMradPerRad = 1000.0;    ///< x[mrad] = x[rad] * kMradPerRad
constexpr double kMgPerG = 1000.0;        ///< t[g/cm^2] = t[mg/cm^2] / kMgPerG
constexpr double kMmPerCm = 10.0;         ///< x[mm] = x[cm] * kMmPerCm
constexpr double kPercent = 100.0;        ///< delta[%] = dp/p * kPercent
constexpr double kPermillePerPercent = 10.0;   ///< delta[10^-3] = delta[%] * kPermillePerPercent
constexpr double kDegPerRad = constants::degrad;   ///< angle[rad] = angle[deg] / kDegPerRad (Fortran degrad)

}  // namespace units
}  // namespace simc

#endif  // SIMC_CORE_UNITS_H
