// src/spectrometers/calo/McCalo.cpp — see include/simc/spectrometers/calo/kernel/McCalo.h (kernel).
#include "simc/spectrometers/calo/kernel/McCalo.h"

#include <cmath>
#include <stdexcept>

#include "simc/stages/IRandom.h"

namespace simc {
namespace calo {

namespace {
// mc_calo.f:31-34
constexpr double h_entr = 30.75;   // calorimeter 1/2 width: (30 blocks * 2.05 cm)/2.0
constexpr double v_entr = 36.9;    // calorimeter 1/2 height: (36 blocks * 2.05 cm)/2.0
const char* const kNames[kNumCounters + 1] = {"", "trials", "slit_hor", "slit_vert", "successes"};
}  // namespace

const char* counterName(int idx) {
    if (idx < 1 || idx > kNumCounters) throw std::out_of_range("calo::counterName: bad index");
    return kNames[idx];
}

// Fortran: calo/mc_calo_recon.f (l.1-60)
void mcCaloRecon(const optics::Track& t, double fry, double delta_y, double delta_z, double drift_to_cal,
                 double& delta_p, double& delta_t, double& delta_phi, double& y_tgt) {
    delta_p = t.dpps;
    delta_phi = (t.xs - fry) / (drift_to_cal - delta_z);
    delta_t = (t.ys - delta_y) / (drift_to_cal - delta_z);
    y_tgt = -100.;
}

// Fortran: calo/mc_calo.f:mc_calo (l.1-191)
McCaloDetail McCalo::mcCalo(McCaloArgs& a, IRandom& rng) {
    McCaloDetail det;
    optics::Track& t = track_;   // COMMON /track/
    // ctheta = cos(th_spec), stheta = sin(th_spec) (l.111-112): unused
    a.ok_spec = false;
    optics::DecayState d;
    d.decay_flag = a.decay_flag;
    d.dflag = false;   // particle has not decayed yet
    d.m2 = a.m2;
    d.ctau = a.ctau;
    d.decdist = a.decdist;
    d.Mh2_final = a.Mh2_final;
    d.rng = &rng;
    d.fortran_compat = compat_;
    cnt_[kTrials] = cnt_[kTrials] + 1;

    t.xs = a.x;
    t.ys = a.y;
    t.zs = a.z;
    t.dxdzs = a.dxdz;
    t.dydzs = a.dydz;
    t.dpps = a.dpp;
    d.ph = a.p_spec * (1. + t.dpps / 100.);   // p

    const double zdrift = a.drift_to_cal;
    optics::project(t, zdrift, a.pathlen, &d);
    if (std::fabs(t.ys) > h_entr) {
        ++cnt_[kSlitHor];
        det.stop_code = kSlitHor;
    } else if (std::fabs(t.xs) > v_entr) {
        ++cnt_[kSlitVert];
        det.stop_code = kSlitVert;
    } else {
        a.x_fp = t.xs;
        a.y_fp = t.ys;
        a.dx_fp = t.dxdzs;
        a.dy_fp = t.dydzs;
        // (mc_calo.f:167: the call of mc_calo_recon is commented out; dpp, dxdz, dydz, y are returned unchanged)
        a.ok_spec = true;
        cnt_[kSuccesses] = cnt_[kSuccesses] + 1;
    }
    // label 500
    a.m2 = d.m2;
    a.decdist = d.decdist;
    a.Mh2_final = d.Mh2_final;
    return det;
}

}  // namespace calo
}  // namespace simc
