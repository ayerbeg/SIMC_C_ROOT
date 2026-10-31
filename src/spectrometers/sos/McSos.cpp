// src/spectrometers/sos/McSos.cpp — see include/simc/spectrometers/sos/kernel/McSos.h (kernel).
// mc_sos (sos/mc_sos.f). The detector hut is in McSosHut.cpp, the ISpectrometer plug-in in SosSpectrometer.cpp.
#include "simc/spectrometers/sos/kernel/McSos.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "simc/stages/IRandom.h"

#ifndef SIMC_SOS_DEFAULT_MATRIX_DIR
#define SIMC_SOS_DEFAULT_MATRIX_DIR "data/matrices/simc_fortran"
#endif

namespace simc {
namespace sos {

namespace {

constexpr int spectr = 2;   // mc_sos.f:36-37, sos is spec #2

// sos/apertures_sos.inc (l.7-19)
constexpr double r2_quad = 163.84;      // 12.7 cm radius^2
constexpr double w_bm01 = 8.0;          // BM01 half-gap in cm
constexpr double w_bm02 = 8.0;          // BM02 half-gap in cm
constexpr double t_bm01_in = 41.32;     // BM01-IN YP top limit
constexpr double b_bm01_in = -30.52;    // BM01-IN YP bot limit
constexpr double t_bm01_out = 53.05;    // BM01-OUT YP top limit
constexpr double b_bm01_out = -65.07;   // BM01-OUT YP bot limit
constexpr double t_bm02_in = 51.73;     // BM02-IN YP top limit
constexpr double b_bm02_in = -66.38;    // BM02-IN YP bot limit
constexpr double t_bm02_out = 51.52;    // BM02-OUT YP top limit
constexpr double b_bm02_out = -55.85;   // BM02-OUT YP bot limit
constexpr double w_exit = 8.57;         // Horiz 1/2-gap exit flng
constexpr double t_exit = 51.52;        // exit flange top limit
constexpr double b_exit = -55.85;       // exit flange bot limit

// mc_sos.f:58-67: 'large' collimator, z-positions of the fixed slit
constexpr double h_entr = 7.201;
constexpr double v_entr = 4.696;
constexpr double h_exit = 7.567;
constexpr double v_exit = 4.935;
constexpr double z_entr = 126.3e0;          // nominally 1.263 m
constexpr double z_exit = z_entr + 6.3e0;   // 6.3 cm thick

constexpr double pi = 3.141592653589793;   // constants.inc (used by the mock sieve only)

const char* const kNames[kNumCounters + 1] = {"",        "trials",   "slit_hor", "slit_vert", "slit_oct",
                                              "quad_in", "quad_mid", "quad_out", "bm01_in",   "bm01_out",
                                              "bm02_in", "bm02_out", "exit",     "hut",       "dc1",
                                              "dc2",     "scin",     "cal",      "successes"};

}  // namespace

const char* counterName(int idx) {
    if (idx < 1 || idx > kNumCounters) throw std::out_of_range("sos::counterName: bad index");
    return kNames[idx];
}

std::string McSos::defaultMatrixDir() { return SIMC_SOS_DEFAULT_MATRIX_DIR; }

McSos::McSos() : eng_(std::make_shared<optics::CosyEngine>()), matrixDir_(defaultMatrixDir()), log_(&std::cout) {}

void McSos::setLog(std::ostream* os) {
    log_ = os;
    mat_.setLog(os);
    eng_->setLog(os);
}

// Fortran: mc_sos.f:155-160 (first_time_here: transp_init, n_classes must be 10)
void McSos::ensureForward() {
    if (eng_->hasForward(spectr)) return;
    const int n_classes = eng_->initForward(spectr, matrixDir_);
    if (n_classes != 10) throw simc::FortranStop("MC_SOS, wrong number of transport classes");
}

// Fortran: mc_sos_recon.f:70-102 (firsttime block)
void McSos::ensureRecon() {
    if (!eng_->hasRecon(spectr)) eng_->initRecon(spectr, matrixDir_);
}

// Fortran: mc_sos.f:mc_sos (l.1-394)
McSosDetail McSos::mcSos(McSosArgs& a, IRandom& rng) {
    McSosDetail det;
    // Initialize ok_spec to .false., reset decay flag (mc_sos.f:115-117)
    a.ok_spec = false;
    optics::DecayState d;
    d.decay_flag = a.decay_flag;
    d.dflag = false;   // particle has not decayed yet
    d.m2 = a.m2;
    d.ctau = a.ctau;
    d.decdist = a.decdist;
    d.Mh2_final = a.Mh2_final;
    d.rng = &rng;
    d.fortran_compat = fortranCompat();   // D-037 F7/F9 (project.f)
    d.trace = transportTrace_;
    if (transportTrace_) transportTrace_->reset();
    cnt_[kTrials] = cnt_[kTrials] + 1;
    det.stop_code = body(a, d, rng, det);
    // label 500: copy back the COMMON / argument variables the Fortran modifies in place
    a.m2 = d.m2;
    a.decdist = d.decdist;
    a.Mh2_final = d.Mh2_final;
    return det;
}

int McSos::body(McSosArgs& a, optics::DecayState& d, IRandom& rng, McSosDetail& det) {
    optics::Track& t = track_;   // COMMON /track/
    double xt, yt;
    xt = a.th_spec;   // avoid 'unused variable' error for th_spec (mc_sos.f:118)

    // Force particles to go through the sieve slit holes, for mock sieve option (mc_sos.f:122-137).
    if (useSieve_) {
        xt = a.x + z_entr * a.dxdz;   // project to collimator
        yt = a.y + z_entr * a.dydz;
        xt = 1.524 * std::round(xt / 1.524);   // shift to nearest hole (nint: half away from zero)
        if (std::fabs(yt) < 1.524) {
            yt = 1.016 * std::round(yt / 1.016);   // middle 3 holes
        } else {
            yt = yt / std::fabs(yt) * (2.54 * std::round((std::fabs(yt) - 1.016) / 2.54) + 1.016);   // outer holes
        }
        const double rt = 0.254 * std::sqrt(rng.uniform());   // distance from center of hole (r=2.54mm)
        const double tht = 2 * pi * rng.uniform();              // angle of offset
        xt = xt + rt * std::cos(tht);
        yt = yt + rt * std::sin(tht);
        a.dxdz = (xt - a.x) / z_entr;   // force to correct angle
        a.dydz = (yt - a.y) / z_entr;
    }

    // Save spectrometer coordinates.
    t.xs = a.x;
    t.ys = a.y;
    t.zs = a.z;
    t.dxdzs = a.dxdz;
    t.dydzs = a.dydz;

    // particle momentum
    t.dpps = a.dpp;
    d.ph = a.p_spec * (1. + t.dpps / 100.);   // p

    // Read in transport coefficients.
    ensureForward();

    optics::CosyEngine& eng = *eng_;
    double& pathlen = a.pathlen;
    double zdrift;

    // Check front of fixed slit (mc_sos.f:167-180).
    zdrift = z_entr;
    optics::project(t, zdrift, pathlen, &d);
    if (std::fabs(t.ys) > h_entr) return stop(kSlitHor);
    if (std::fabs(t.xs) > v_entr) return stop(kSlitVert);
    if (std::fabs(t.xs) > (-v_entr / h_entr * std::fabs(t.ys) + 3 * v_entr / 2)) return stop(kSlitOct);

    // Check back of fixed slit, at about 1.325 meter (mc_sos.f:184-197)
    zdrift = z_exit - z_entr;
    optics::project(t, zdrift, pathlen, &d);
    if (std::fabs(t.ys) > h_exit) return stop(kSlitHor);
    if (std::fabs(t.xs) > v_exit) return stop(kSlitVert);
    if (std::fabs(t.xs) > (-v_exit / h_exit * std::fabs(t.ys) + 3 * v_exit / 2)) return stop(kSlitOct);

    // Go to quad IN mag bound. (150 cm from target) (mc_sos.f:202-208)
    if (!eng.adrift(spectr, 1) && log_) *log_ << " Transformation #1 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 1) - z_exit;
    optics::project(t, zdrift, pathlen, &d);   // project and decay
    if ((t.xs * t.xs + t.ys * t.ys) > r2_quad) return stop(kQuadIn);

    // Go to quad midplane.
    eng.transp(spectr, 2, t, 35.0e0, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r2_quad) return stop(kQuadMid);

    // Go to quad OUT mag boundary.
    eng.transp(spectr, 3, t, 35.0e0, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r2_quad) return stop(kQuadOut);

    // Go to BM01 IN magnetic boundary. Find intersection with rotated aperture plane (mc_sos.f:228-238).
    eng.transp(spectr, 4, t, 80.0e0, pathlen, &d);
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, -45.0e0, xt, yt, fortranCompat());
    if ((yt > w_bm01) || (-yt > (w_bm01 - 0.05 * 2.54)) || (-xt > t_bm01_in) || (-xt < b_bm01_in))
        return stop(kBm01In);

    // Go to BM01 OUT magnetic boundary (mc_sos.f:243-251).
    eng.transp(spectr, 5, t, 169.52e0, pathlen, &d);
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, 45.0e0, xt, yt, fortranCompat());
    if ((std::fabs(yt) > w_bm01) || (-xt > t_bm01_out) || (-xt < b_bm01_out)) return stop(kBm01Out);

    // Go to BM02 IN magnetic boundary (mc_sos.f:256-264).
    eng.transp(spectr, 6, t, 80.80e0, pathlen, &d);
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, 49.0e0, xt, yt, fortranCompat());
    if ((std::fabs(yt) > w_bm02) || (-xt > t_bm02_in) || (-xt < b_bm02_in)) return stop(kBm02In);

    // Go to BM02 OUT magnetic boundary (mc_sos.f:269-277).
    eng.transp(spectr, 7, t, 77.06e0, pathlen, &d);
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, 57.0e0, xt, yt, fortranCompat());
    if ((std::fabs(yt) > w_bm02) || (-xt > t_bm02_out) || (-xt < b_bm02_out)) return stop(kBm02Out);

    // Current (may 99) SOS model: MEs to the old exit flange (mc_sos.f:307-316). The aperture plane is rotated by 45
    // degrees (a 57-degree call is commented out, l.310).
    eng.transp(spectr, 8, t, 43.82e0, pathlen, &d);   // old flange
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, 45.0e0, xt, yt, fortranCompat());
    if ((std::fabs(yt) > w_exit) || (-xt > t_exit) || (-xt < b_exit)) return stop(kExit);

    // Go to aperture at -10.84 cm (from -55.18 to -10.84 is 44.34 cm), mc_sos.f:326-333.
    // Aperture is +/- 18.6944cm wide at top, +/- 10.998 cm wide at bottom, and +/- 37.694 cm tall.
    eng.transp(spectr, 9, t, 44.34e0, pathlen, &d);   // new aperture
    const double tmpwidth = 10.998 + 0.10209 * (t.xs + 37.694);
    if ((std::fabs(t.xs) > 37.694) || std::fabs(t.ys) > tmpwidth) return stop(kExit);

    // Go to window at -3.22 cm (from -10.84 to -3.22 is 7.62 cm), mc_sos.f:339-347.
    // FORTRAN-DIFF-CANDIDATE(D-059 P14.1): a straight line WITHOUT project (the project call is commented out, l.337):
    // no pathlength, no decay check over these 7.62 cm, zs unchanged; the hut then starts at -3.206267 (mc_sos.f:349-360
    // "WE REALLY NEED TO FIX THIS UP"). Question: should pathlen include the 7.62 cm (FP path / time of flight)?
    t.xs = t.xs + 7.62 * t.dxdzs;
    t.ys = t.ys + 7.62 * t.dydzs;
    // DC1 is +/-32cm (offset by ~8cm(at the moment), so this is almost non-cut.
    if ((std::fabs(t.xs) > 38.1) || std::fabs(t.ys) > 12.7) return stop(kExit);

    // If we get this far, the particle is in the hut.
    cnt_[kHut] = cnt_[kHut] + 1;
    det.reached_hut = true;

    // and track through the detector hut (mc_sos.f:368-370)
    const int hutCode = mcSosHut(d, a, -3.206267e0, rng, det);
    if (hutCode != kSuccess) return hutCode;   // (.not.ok) goto 500

    // replace xs,ys,... with 'tracked' quantities.
    t.xs = a.x_fp;
    t.ys = a.y_fp;
    t.dxdzs = a.dx_fp;
    t.dydzs = a.dy_fp;

    // Reconstruct target quantities.
    ensureRecon();
    double dpp_recon, dth_recon, dph_recon, y_recon;
    eng.recon(spectr, t, a.fry, dpp_recon, dth_recon, dph_recon, y_recon);

    // Fill output to return to main code
    a.dpp = dpp_recon;
    a.dxdz = dph_recon;
    a.dydz = dth_recon;
    a.y = y_recon;

    a.ok_spec = true;
    cnt_[kSuccesses] = cnt_[kSuccesses] + 1;
    return kSuccess;
}

}  // namespace sos
}  // namespace simc
