// include/simc/core/Units.h
//
// THE unit convention of SIMC++ (task P7.4). Everything inside the framework — RunSetup, events, stage interfaces,
// kernels — uses the Fortran SIMC internal units:
//
//   energy, mass        MeV            momentum   MeV/c          mass^2   MeV^2
//   length, position    cm             angle      rad            slopes   dx/dz, dy/dz (dimensionless)
//   target thickness    g/cm^2         density    g/cm^3         delta    dp/p in %
//   cross sections      as the Fortran physics routines return them (documented per model)
//
// Exceptions are the INPUT units of the .inp / JSON parameters (SimcParameters, e.g. degrees, mg/cm^2, mrad), converted
// once by the set-up (dbase.f), and the OUTPUT units of the ntuple (GeV, mrad, ...; results_write.f), converted by the
// output code. Both conversions use the named factors below. Fortran kernels keep their own literals where the Fortran
// writes them (e.g. brem.f works in GeV with "/1000."), so that the ports stay line-by-line comparable.
//
// Every factor is exact and is used with the SAME operation as the Fortran (x / 1000. stays a division by 1000), so
// replacing a literal by its name does not change a single bit.
#ifndef SIMC_CORE_UNITS_H
#define SIMC_CORE_UNITS_H

#include "simc/core/SimcConstants.h"

namespace simc {
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
