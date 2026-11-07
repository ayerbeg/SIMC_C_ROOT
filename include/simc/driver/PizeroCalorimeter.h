// include/simc/driver/PizeroCalorimeter.h
//
// The pi0 -> 2 gamma detection in a calorimeter hadron arm (task P14.3), as a driver component of its own (it uses only
// the ISpectrometer interface, so any calorimeter plug-in can be used): the photons of the
// pi0 decay (ntup%gamma1/2, from pizero_decay in complete_ev) are rotated to the calorimeter frame and sent one by one
// through the arm's ISpectrometer (normally "calo"); the arm accepts the event if both photons (pizero_ngamma = 2) or at
// least one (pizero_ngamma = 1) hit.
// Fortran: simc.f:montecarlo l.1489-1556 (rotation l.1495-1515, the two mc_calo calls l.1517-1544, ok_P_arm l.1546-1556)
// and l.1597-1609 (FP / RECON copies).
//
// Fortran behaviour (docs/FORTRAN_DIFFERENCES.md, D-062):
//   * FORTRAN-BUG(D-062 G1), corrected unless fortran_compat: the photons are passed m2 = Mh2 = Mpi0^2 and decay_flag =
//     doing_decay (simc.f:1518-1522/1535-1539), so every pi0 run with doing_decay = 1 STOPs in project.f:76-90 ('cannot
//     decay particle with mass= 134.98'). Corrected: a photon is massless and never decays (m2 = 0, decay_flag = false).
//   * FORTRAN-BUG(D-062 G3), corrected unless fortran_compat: pathlen is set to 0 once before the hadron arm
//     (simc.f:1462); the second photon starts from the first photon's path length, so FP%p%path is the sum of both
//     (~600 cm for drift_to_cal 300). Corrected: each photon starts at 0, and FP%p%path is the path of the photon whose
//     position is reported in FP%p.
//   * FORTRAN-QUIRK(D-062 G4): FP%p (ntuple *fp of the hadron arm) = the LAST photon that hit (mc_calo writes x_fp only
//     on success); with pizero_ngamma = 1 and photon 2 missed, photon 1. xcal/ycal_gamma = -1e10 for a missed photon.
//   * the photon calls use the same vertex position (x, y, z at z = 0) and the pi0 delta; dx/dz = ex/ez, dy/dz = ey/ez.
#ifndef SIMC_DRIVER_PIZEROCALORIMETER_H
#define SIMC_DRIVER_PIZEROCALORIMETER_H

#include <array>
#include <string>

#include "simc/stages/ISpectrometer.h"

namespace simc {
namespace driver {

/// One photon in the calorimeter frame (simc.f:1495-1515): ex, ey, ez and the slopes ex/ez, ey/ez.
struct PhotonSlopes {
    double ex{0}, ey{0}, ez{0};
    double dxdz{0}, dydz{0};
};
/// simc.f:1495-1515: rotate the lab photon (E, px, py, pz; z along the beam) by the arm angle theta to the frame whose z
/// axis points to the calorimeter centre. arm 8 (SHMS side): ey = py cos - pz sin, ez = py sin + pz cos; arm 7 (HMS side):
/// ey = py cos + pz sin, ez = -py sin + pz cos.
PhotonSlopes rotateToCalorimeter(const std::array<double, 4>& gamma, double theta, int arm);

struct PizeroCaloResult {
    bool ok{false};                    ///< ok_P_arm
    bool ok_gamma1{false}, ok_gamma2{false};
    double xcal_gamma1{-1.0e10}, ycal_gamma1{-1.0e10}, xcal_gamma2{-1.0e10}, ycal_gamma2{-1.0e10};
    SpecResult fp;                     ///< the photon reported in FP%p (last hit) — x_fp..dy_fp
    double pathlen{0};                 ///< FP%p%path
    double dxdz_last{0}, dydz_last{0}; ///< dx_P_arm/dy_P_arm after the calls (photon 2 slopes)
    std::string stop_reason;           ///< "gamma1:<reason>" / "gamma2:<reason>" of the photons that missed
};

class PizeroCalorimeter {
public:
    /// `arm` = the hadron arm code (7 or 8), `ngamma` = pizero_ngamma, `fortran_compat` = keep D-062 G1/G3.
    PizeroCalorimeter(int arm, int ngamma, bool fortran_compat) : arm_(arm), ngamma_(ngamma), compat_(fortran_compat) {}

    /// `base` = the hadron-arm input of simc.f (p_spec, th_spec, delta, x, y, z at z = 0, m2 = Mh2, flags, raster);
    /// gamma1/gamma2 = ntup%gamma1/2 (E, px, py, pz, MeV). Throws FortranStop for pizero_ngamma not 1 or 2 (simc.f:1555).
    PizeroCaloResult detect(ISpectrometer& calo, const SpecInput& base, const std::array<double, 4>& gamma1,
                            const std::array<double, 4>& gamma2, IRandom& rng) const;

private:
    int arm_, ngamma_;
    bool compat_;
};

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_PIZEROCALORIMETER_H
