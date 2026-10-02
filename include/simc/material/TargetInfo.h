// include/simc/material/TargetInfo.h
//
// Fortran: target.inc type target_info / COMMON /target_info/ targ (l.38-60), in the units the
// Fortran uses after dbase.f: thick g/cm^2 (dbase.f:451), length cm (452), angle rad (453),
// X0 g/cm^2, X0_cm cm, Eloss MeV, teff radiation lengths, musc_max rad, Coulomb MeV.
#ifndef SIMC_MATERIAL_TARGETINFO_H
#define SIMC_MATERIAL_TARGETINFO_H

#include "simc/core/RunSetup.h"

namespace simc {
namespace material {

/// Fortran: target.inc type target_function / target_info. Since P7.4 the target model works directly on the core
/// RunSetup types (they had a field-by-field duplicate here): index [i-1] for the Fortran (i), 1 beam, 2 e-, 3 hadron.
using TargetFunction = simc::TargetFunction;
using TargetInfo = simc::TargetInfo;

/// Fortran: target.f:316 type limits (min, max)
struct Limits {
    double min{0}, max{0};
};

}  // namespace material
}  // namespace simc

#endif  // SIMC_MATERIAL_TARGETINFO_H
