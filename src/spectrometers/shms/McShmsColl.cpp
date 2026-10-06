// src/spectrometers/shms/McShmsColl.cpp
// Fortran: shms/mc_shms_coll.f:mc_shms_coll (l.1-149) — 20-step collimator stepping for pions/muons
// (using_SHMScoll), the SHMS version of mc_hms_coll: octagon 8.50/12.50 -> 8.65/12.85, 6.35 cm thick,
// coll_radl = 0.42084 (HMS value scaled by 6.35/6.30). Hadronic absorption uses hms/pion_coll_absorb.f.
#include <cmath>

#include "simc/spectrometers/hms/kernel/McHms.h"   // hms::pionCollAbsorb
#include "simc/spectrometers/shms/kernel/McShms.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace shms {

bool McShms::mcShmsColl(optics::DecayState& d, double p_spec, double& pathlen, IRandom& rng) {
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

    // FORTRAN-BUG(D-021, D-037 F1): original Fortran shms/mc_shms_coll.f:79/86/93: `shmsSTOP_slit_hor =
    // shmsSTOP_slit_hor + 1` (and _vert, _oct) for EVERY 1/20 step outside the opening, also for rays that pass;
    // corrected as in mc_hms_coll (McHmsColl.cpp): each counter at most once per ray (the unused slit_*_flag flags).
    // fortran_compat counts per step.
    const bool perStep = fortranCompat();
    bool slit_h_flag = false, slit_v_flag = false, slit_o_flag = false;
    auto count = [&](int code, bool& flag) {
        if (perStep || !flag) cnt_[code] = cnt_[code] + 1;
        flag = true;
    };

    for (int n = 1; n <= nstep; ++n) {
        bool step_flag = false;
        if (std::fabs(t.ys - y_off) > h_step) {
            step_flag = true;
            count(kSlitHor, slit_h_flag);
            thick_temp_ = step_size;
        }
        if (std::fabs(t.xs - x_off) > v_step) {
            step_flag = true;
            count(kSlitVert, slit_v_flag);
            thick_temp_ = step_size;
        }
        if (std::fabs(t.xs - x_off) > (-v_step / h_step * std::fabs(t.ys - y_off) + 3 * v_step / 2)) {
            step_flag = true;
            count(kSlitOct, slit_o_flag);
            thick_temp_ = step_size;
        }
        const double thick = thick_temp_;

        if (step_flag) {
            if (m2 > 12000. && m2 < 20000.) {   // pion: hadronic interaction
                const double trans = hms::pionCollAbsorb(p, thick);
                const double rantemp = rng.uniform();
                if (collTrace_) hms::traceAbsorb(*collTrace_, p, thick, trans, rantemp);
                if (rantemp > trans) return false;
            }
            // FORTRAN-DIFF-CANDIDATE(D-037 F3): shms/mc_shms_coll.f:117-120: musc (2 gauss1) and the Landau-sampled enerloss_new
            // (typeflag 1, 1 gauss1) are called whatever ms_flag / mc_smear say, so the collimator scatters and
            // straggles even with mc_smear = 0. Question: is the collimator scattering intended physics (it is the
            // purpose of the stepping) or should it follow ms_flag? Kept in both modes.
            const double coll_radw = thick / coll_radl;
            mat_.musc(m2, p, coll_radw, t.dydzs, t.dxdzs, rng);
            epart = std::sqrt(p * p + m2);
            eloss = mat_.energyLoss(thick, coll_dens, zcoll, acoll, epart, std::sqrt(m2), 1, rng);
            if (collTrace_) hms::traceEloss(*collTrace_, thick, epart, std::sqrt(m2), eloss);
            epart = epart - eloss;
            if (epart < std::sqrt(m2)) return false;
            // FORTRAN-BUG(D-037 F4): original Fortran shms/mc_shms_coll.f:123-126 `if(epart.lt.sqrt(m2)) goto 500` is
            // unreachable (enerloss_new clamps the loss); corrected: a clamped loss stops the particle (sSTOP_coll),
            // see McHmsColl.cpp. fortran_compat keeps the Fortran.
            if (!fortranCompat() && mat_.lastEnergyLossClamped()) {
                if (collTrace_) collTrace_->ranged_out = true;
                return false;
            }
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
