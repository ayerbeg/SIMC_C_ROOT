/**
 * @file
 * @ingroup material
 * @brief Target description types of the target model (target.inc type target_info, Limits).
 *
 * Fortran: target.inc type target_info / COMMON /target_info/ targ (l.38-60), in the units the
 * Fortran uses after dbase.f: thick g/cm^2 (dbase.f:451), length cm (452), angle rad (453),
 * X0 g/cm^2, X0_cm cm, Eloss MeV, teff radiation lengths, musc_max rad, Coulomb MeV.
 */
#ifndef SIMC_MATERIAL_TARGETINFO_H
#define SIMC_MATERIAL_TARGETINFO_H

#include "simc/core/RunSetup.h"

namespace simc {
namespace material {

/// @addtogroup material
/// @{

/// Fortran: target.inc type target_function / target_info. Since P7.4 the target model works directly on the core
/// RunSetup types (they had a field-by-field duplicate here): index [i-1] for the Fortran (i), 1 beam, 2 e-, 3 hadron.
using TargetFunction = simc::TargetFunction;
/// The target description (simc::TargetInfo, COMMON /target_info/ targ).
using TargetInfo = simc::TargetInfo;

/// Fortran: target.f:316 type limits (min, max)
struct Limits {
    double min{0} /**< @copydoc max */, max{0};  ///< lower and upper limit (units of the quantity).
};

/// @}

}  // namespace material
}  // namespace simc

#endif  // SIMC_MATERIAL_TARGETINFO_H
