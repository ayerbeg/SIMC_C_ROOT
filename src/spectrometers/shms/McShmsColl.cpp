// src/spectrometers/shms/McShmsColl.cpp
// Fortran: shms/mc_shms_coll.f:mc_shms_coll (l.1-149) — 20-step collimator stepping for pions/muons
// (using_SHMScoll), the SHMS version of mc_hms_coll: octagon 8.50/12.50 -> 8.65/12.85, 6.35 cm thick,
// coll_radl = 0.42084 (HMS value scaled by 6.35/6.30). Hadronic absorption uses hms/pion_coll_absorb.f.
#include <cmath>

#include "simc/spectrometers/hms/HmsSpectrometer.h"   // hms::pionCollAbsorb
#include "simc/spectrometers/shms/ShmsSpectrometer.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace shms {

bool ShmsSpectrometer::mcShmsColl(optics::DecayState& d, double p_spec, double& pathlen, IRandom& rng) {
    constexpr double h_entr = 8.50;
    constexpr double v_entr = 12.50;
    constexpr double h_exit = 8.65;
    constexpr double v_exit = 12.85;
    constexpr double x_off = 0.000;
    constexpr double y_off = +0.000;
    constexpr int nstep = 20;
    constexpr double coll_dens = 17.0;
    constexpr double zcoll = 69.45;
    constexpr double acoll = 171.56797;
    constexpr double coll_radl = 0.42084;

    optics::Track& t = track_;
    double& m2 = d.m2;
    double& p = d.ph;

    double eloss = 0.;
    const double step_size = 6.35 / nstep;
    double h_step = h_entr;
    double v_step = v_entr;
    double epart = std::sqrt(p * p + m2);

    for (int n = 1; n <= nstep; ++n) {
        bool step_flag = false;
        if (std::fabs(t.ys - y_off) > h_step) {
            step_flag = true;
            cnt_[kSlitHor] = cnt_[kSlitHor] + 1;
            thick_temp_ = step_size;
        }
        if (std::fabs(t.xs - x_off) > v_step) {
            step_flag = true;
            cnt_[kSlitVert] = cnt_[kSlitVert] + 1;
            thick_temp_ = step_size;
        }
        if (std::fabs(t.xs - x_off) > (-v_step / h_step * std::fabs(t.ys - y_off) + 3 * v_step / 2)) {
            step_flag = true;
            cnt_[kSlitOct] = cnt_[kSlitOct] + 1;
            thick_temp_ = step_size;
        }
        const double thick = thick_temp_;

        if (step_flag) {
            if (m2 > 12000. && m2 < 20000.) {   // pion: hadronic interaction
                const double trans = hms::pionCollAbsorb(p, thick);
                const double rantemp = rng.uniform();
                if (rantemp > trans) return false;
            }
            const double coll_radw = thick / coll_radl;
            mat_.musc(m2, p, coll_radw, t.dydzs, t.dxdzs, rng);
            epart = std::sqrt(p * p + m2);
            eloss = mat_.energyLoss(thick, coll_dens, zcoll, acoll, epart, std::sqrt(m2), 1, rng);
            epart = epart - eloss;
            if (epart < std::sqrt(m2)) return false;
            p = std::sqrt(epart * epart - m2);
            t.dpps = 100. * (p / p_spec - 1.);
        }

        optics::project(t, step_size, pathlen, &d);

        h_step = h_step + (h_exit - h_entr) / nstep;
        v_step = v_step + (v_exit - v_entr) / nstep;
    }
    return true;
}

}  // namespace shms
}  // namespace simc
