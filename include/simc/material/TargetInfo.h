// include/simc/material/TargetInfo.h
//
// Fortran: target.inc type target_info / COMMON /target_info/ targ (l.38-60), in the units the
// Fortran uses after dbase.f: thick g/cm^2 (dbase.f:451), length cm (452), angle rad (453),
// X0 g/cm^2, X0_cm cm, Eloss MeV, teff radiation lengths, musc_max rad, Coulomb MeV.
#ifndef SIMC_MATERIAL_TARGETINFO_H
#define SIMC_MATERIAL_TARGETINFO_H

#include <array>

namespace simc {
namespace material {

/// Fortran: target.inc type target_function
struct TargetFunction {
    double ave{0}, min{0}, max{0};
};

/// Fortran: target.f:316 type limits (min, max)
struct Limits {
    double min{0}, max{0};
};

/// Fortran: target.inc type target_info. Index [i-1] for the Fortran (i), i = 1 beam, 2 e-, 3 hadron.
struct TargetInfo {
    double A{0}, Z{0}, N{0}, mass_amu{0}, M{0}, mrec_amu{0}, Mrec{0}, rho{0};
    double thick{0};    ///< g/cm^2
    double angle{0};    ///< rad
    double abundancy{0};
    double length{0};   ///< cm
    double zoffset{0}, X0{0}, X0_cm{0}, L1{0}, L2{0};
    double fr1{0}, fr2{0}, xoffset{0}, yoffset{0};
    int fr_pattern{0}, can{0};
    std::array<TargetFunction, 3> Eloss{}, teff{};
    TargetFunction Coulomb;
    double Coulomb_constant{0};
    std::array<double, 3> musc_max{};
    double musc_nsig_max{0};
    double Mtar_struck{0}, Mrec_struck{0};
};

}  // namespace material
}  // namespace simc

#endif  // SIMC_MATERIAL_TARGETINFO_H
