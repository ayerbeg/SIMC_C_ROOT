// src/driver/PizeroCalorimeter.cpp — see include/simc/driver/PizeroCalorimeter.h.
#include "simc/driver/PizeroCalorimeter.h"

#include <cmath>

#include "simc/core/Errors.h"

namespace simc {
namespace driver {

// Fortran: simc.f:1495-1515
PhotonSlopes rotateToCalorimeter(const std::array<double, 4>& g, double theta, int arm) {
    PhotonSlopes s;
    if (arm == 8) {
        s.ex = g[1];
        s.ey = g[2] * std::cos(theta) - g[3] * std::sin(theta);
        s.ez = g[2] * std::sin(theta) + g[3] * std::cos(theta);
    } else {
        s.ex = g[1];
        s.ey = g[2] * std::cos(theta) + g[3] * std::sin(theta);
        s.ez = -g[2] * std::sin(theta) + g[3] * std::cos(theta);
    }
    s.dxdz = s.ex / s.ez;   // simc.f:1518-1519 / 1532-1533
    s.dydz = s.ey / s.ez;
    return s;
}

// Fortran: simc.f:1489-1556
PizeroCaloResult PizeroCalorimeter::detect(ISpectrometer& calo, const SpecInput& base,
                                           const std::array<double, 4>& gamma1, const std::array<double, 4>& gamma2,
                                           IRandom& rng) const {
    PizeroCaloResult res;
    const PhotonSlopes s1 = rotateToCalorimeter(gamma1, base.th_spec, arm_);
    const PhotonSlopes s2 = rotateToCalorimeter(gamma2, base.th_spec, arm_);

    SpecInput in = base;
    if (!compat_) {
        // FORTRAN-BUG(D-062 G1): original Fortran simc.f:1518-1522/1535-1539: `call mc_calo(..., m2, mc_smear,
        // mc_smear, doing_decay, ...)` with m2 = Mh2 = Mpi0**2: the "photon" can decay and project.f STOPs; corrected
        // because a photon is massless and stable (m2 = 0, decay_flag = .false.).
        in.m2 = 0.0;
        in.decay_flag = false;
    }
    // first photon
    in.dxdz = s1.dxdz;
    in.dydz = s1.dydz;
    in.pathlen_start = base.pathlen_start;
    const SpecResult r1 = calo.transport(in, rng);
    res.ok_gamma1 = r1.ok;
    if (r1.ok) {
        res.xcal_gamma1 = r1.x_fp;
        res.ycal_gamma1 = r1.y_fp;
    } else {
        res.stop_reason = "gamma1:" + r1.stop_reason;
    }
    // second photon
    in.dxdz = s2.dxdz;
    in.dydz = s2.dydz;
    // FORTRAN-BUG(D-062 G3): original Fortran simc.f:1462 sets pathlen = 0 once; the second mc_calo call starts from the
    // first photon's path length (FP%p%path ~ twice the drift); corrected: every photon starts at 0.
    in.pathlen_start = compat_ ? r1.pathlen : 0.0;
    const SpecResult r2 = calo.transport(in, rng);
    res.ok_gamma2 = r2.ok;
    if (r2.ok) {
        res.xcal_gamma2 = r2.x_fp;
        res.ycal_gamma2 = r2.y_fp;
    } else {
        res.stop_reason += (res.stop_reason.empty() ? "" : " ") + std::string("gamma2:") + r2.stop_reason;
    }
    res.dxdz_last = s2.dxdz;
    res.dydz_last = s2.dydz;

    // simc.f:1546-1556
    if (ngamma_ == 2) {
        res.ok = r1.ok && r2.ok;   // require both photons
    } else if (ngamma_ == 1) {
        res.ok = r1.ok || r2.ok;   // require only one photon
    } else {
        throw FortranStop("pizero_ngamma not set correctly (should be 1 or 2), stopping");
    }
    // FP%p (simc.f:1599-1603): the montecarlo locals xfp..dyfp hold the last photon that hit (D-062 G4)
    const SpecResult& last = r2.ok ? r2 : r1;
    res.fp = last;
    res.pathlen = compat_ ? r2.pathlen : last.pathlen;   // simc.f:1604 FP%p%path = pathlen (both calls, G3)
    return res;
}

}  // namespace driver
}  // namespace simc
