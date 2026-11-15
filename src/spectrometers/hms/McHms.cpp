// src/spectrometers/hms/McHms.cpp — see include/simc/spectrometers/hms/kernel/McHms.h (kernel).
// mc_hms (mc_hms.f) and hit_dipole. The collimator is in McHmsColl.cpp, the detector hut in McHmsHut.cpp, the
// ISpectrometer plug-in in HmsSpectrometer.cpp.
#include "simc/spectrometers/hms/kernel/McHms.h"

#include "simc/core/DataPath.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "simc/stages/IRandom.h"


namespace simc {
namespace hms {

namespace {

constexpr int spectr = 1;   // mc_hms.f:41, HMS is spec #1

// apertures_hms.inc
constexpr double r_Q1 = 20.50;
constexpr double r_Q2 = 30.22;
constexpr double r_Q3 = 30.22;
constexpr double x_d1 = 34.29, y_d1 = 12.07;
constexpr double x_d2 = 27.94, y_d2 = 18.42;
constexpr double x_d3 = 13.97, y_d3 = 18.95;
constexpr double x_d4 = 1.956, y_d4 = 20.32;
constexpr double x_d5 = 27.94, y_d5 = 12.065, r_d5 = 6.35;
constexpr double a_d6 = -0.114, b_d6 = 20.54;

// mc_hms.f:36, 70-94 (HMS-100 collimator)
constexpr double h_entr = 4.575;
constexpr double v_entr = 11.646;
constexpr double h_exit = 4.759;
constexpr double v_exit = 12.114;
constexpr double x_off = +0.000;
constexpr double y_off = +0.028;
constexpr double z_off = +40.17;
constexpr double z_entr = 126.2e0 + z_off;   // nominally 1.262 m
constexpr double z_exit = z_entr + 6.3e0;    // 6.3 cm thick
constexpr double z_dip1 = 64.77e0;           // end of 26.65 inch pipe
constexpr double z_dip2 = z_dip1 + 297.18e0;   // ~117 inch pipe
constexpr double z_dip3 = z_dip2 + 115.57e0;   // 45.5 inch pipe
// real*8 x_offset_pipes/2.8/, y_offset_pipes/0.0/ (SAVEd variables, never changed)
constexpr double x_offset_pipes = 2.8;
constexpr double y_offset_pipes = 0.0;
constexpr double pi = 3.141592653589793;   // constants.inc

const char* const kNames[kNumCounters + 1] = {"",       "trials", "slit_hor", "slit_vert", "slit_oct", "Q1_in",
                                              "Q1_mid", "Q1_out", "Q2_in",    "Q2_mid",    "Q2_out",   "Q3_in",
                                              "Q3_mid", "Q3_out", "D1_in",    "D1_out",    "hut",      "dc1",
                                              "dc2",    "scin",   "cal",      "successes", "coll"};

}  // namespace

const char* counterName(int idx) {
    if (idx < 1 || idx > kNumCounters) throw std::out_of_range("hms::counterName: bad index");
    return kNames[idx];
}

// Fortran: mc_hms.f:hit_dipole (l.445-492). Returns .true. if the point is OUTSIDE the dipole aperture.
bool hitDipole(double x, double y) {
    const double x_local = std::fabs(x);   // symmetry: first quadrant only
    const double y_local = std::fabs(y);
    const bool check1 = (x_local <= x_d1) && (y_local <= y_d1);
    const bool check2 = (x_local <= x_d2) && (y_local <= y_d2);
    const bool check3 = (x_local <= x_d3) && (y_local <= y_d3);
    const bool check4 = (x_local <= x_d4) && (y_local <= y_d4);
    // the rounded corner ((x-x_d5)**2 + (y-y_d5)**2 .le. r_d5**2)
    const double dx5 = x_local - x_d5, dy5 = y_local - y_d5;
    const bool check5 = (dx5 * dx5 + dy5 * dy5) <= r_d5 * r_d5;
    // the slanted piece
    const bool check6 = (x_local >= x_d4) && (x_local <= x_d3) && ((y_local - a_d6 * x_local - b_d6) <= 0.0);
    const bool miss_dipole = check1 || check2 || check3 || check4 || check5 || check6;
    return !miss_dipole;
}

std::string McHms::defaultMatrixDir() { return simc::DataDir("matrices/simc_fortran"); }

McHms::McHms()
    : eng_(std::make_shared<optics::CosyEngine>()), matrixDir_(defaultMatrixDir()), log_(&std::cout) {}

void McHms::setLog(std::ostream* os) {
    log_ = os;
    mat_.setLog(os);
    eng_->setLog(os);
}

// Fortran: mc_hms.f:185-190 (first_time_here: transp_init, n_classes must be 12)
void McHms::ensureForward() {
    if (eng_->hasForward(spectr)) return;
    const int n_classes = eng_->initForward(spectr, matrixDir_);
    if (n_classes != 12) throw simc::FortranStop("MC_HMS, wrong number of transport classes");
}

// Fortran: mc_hms_recon.f firsttime block
void McHms::ensureRecon() {
    if (!eng_->hasRecon(spectr)) eng_->initRecon(spectr, matrixDir_);
}

// Fortran: mc_hms.f:mc_hms (l.1-441)
McHmsDetail McHms::mcHms(McHmsArgs& a, IRandom& rng) {
    McHmsDetail det;
    // Initialize ok_spec to .false., reset decay flag (mc_hms.f:150-152)
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
    if (collTrace_) collTrace_->reset();
    cnt_[kTrials] = cnt_[kTrials] + 1;
    det.stop_code = body(a, d, rng, det);
    // label 500: copy back the COMMON / argument variables the Fortran modifies in place
    a.m2 = d.m2;
    a.decdist = d.decdist;
    a.Mh2_final = d.Mh2_final;
    return det;
}

int McHms::body(McHmsArgs& a, optics::DecayState& d, IRandom& rng, McHmsDetail& det) {
    optics::Track& t = track_;   // COMMON /track/
    double xt, yt;
    xt = a.th_spec;   // avoid 'unused variable' error for th_spec (mc_hms.f:153)

    // Force particles to go through the sieve slit holes, for mock sieve option (mc_hms.f:157-168).
    if (useSieve_) {
        xt = a.x + z_entr * a.dxdz;   // project to collimator
        yt = a.y + z_entr * a.dydz;
        xt = 2.540 * std::round(xt / 2.54);   // shift to nearest hole (nint: half away from zero)
        yt = 1.524 * std::round(yt / 1.524);
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

    // Check front of fixed slit.
    zdrift = z_entr;
    optics::project(t, zdrift, pathlen, &d);   // project and decay

    // SECTION ADDED FOR HMS COLLIMATOR (mc_hms.f:206-254)
    if (a.using_HMScoll && (d.m2 > 100.0 * 100.0) && (d.m2 < 200.0 * 200.0)) {   // pions/muons
        const bool Hcollsuccess = mcHmsColl(d, a.p_spec, pathlen, rng);
        if (!Hcollsuccess) {
            saveCollSuccess_ = 1.5;   // COMMON /save_coll/, never read in SIMC (D-037 F10)
            return stop(kColl);
        }
    } else {   // just apertures at front and back
        if (std::fabs(t.ys - y_off) > h_entr) return stop(kSlitHor);
        if (std::fabs(t.xs - x_off) > v_entr) return stop(kSlitVert);
        if (std::fabs(t.xs - x_off) > (-v_entr / h_entr * std::fabs(t.ys - y_off) + 3 * v_entr / 2))
            return stop(kSlitOct);

        // Check back of fixed slit.
        zdrift = z_exit - z_entr;
        optics::project(t, zdrift, pathlen, &d);
        if (std::fabs(t.ys - y_off) > h_exit) return stop(kSlitHor);
        if (std::fabs(t.xs - x_off) > v_exit) return stop(kSlitVert);
        if (std::fabs(t.xs - x_off) > (-v_exit / h_exit * std::fabs(t.ys - y_off) + 3 * v_exit / 2))
            return stop(kSlitOct);
    }

    // Go to Q1 IN mag bound. Drift rather than using COSY matrices
    if (!eng.adrift(spectr, 1) && log_) *log_ << " Transformation #1 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 1) - z_exit;
    optics::project(t, zdrift, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q1 * r_Q1) return stop(kQ1In);

    // Check aperture at 2/3 of Q1.
    eng.transp(spectr, 2, t, 125.233e0, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q1 * r_Q1) return stop(kQ1Mid);

    // Go to Q1 OUT mag boundary.
    eng.transp(spectr, 3, t, 62.617e0, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q1 * r_Q1) return stop(kQ1Out);

    // Go to Q2 IN mag bound.
    if (!eng.adrift(spectr, 4) && log_) *log_ << " Transformation #4 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 4);
    optics::project(t, zdrift, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q2 * r_Q2) return stop(kQ2In);

    // Check aperture at 2/3 of Q2.
    eng.transp(spectr, 5, t, 143.90e0, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q2 * r_Q2) return stop(kQ2Mid);

    // Go to Q2 OUT mag boundary.
    eng.transp(spectr, 6, t, 71.95e0, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q2 * r_Q2) return stop(kQ2Out);

    // Go to Q3 IN mag bound.
    if (!eng.adrift(spectr, 7) && log_) *log_ << " Transformation #7 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 7);
    optics::project(t, zdrift, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q3 * r_Q3) return stop(kQ3In);

    // Check aperture at 2/3 of Q3.
    eng.transp(spectr, 8, t, 143.8e0, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q3 * r_Q3) return stop(kQ3Mid);

    // Go to Q3 OUT mag boundary.
    eng.transp(spectr, 9, t, 71.9e0, pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q3 * r_Q3) return stop(kQ3Out);

    // Go to D1 IN magnetic boundary, find intersection with rotated aperture plane.
    if (!eng.adrift(spectr, 10) && log_) *log_ << " Transformation #10 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 10);
    optics::project(t, zdrift, pathlen, &d);
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, -6.0e0, xt, yt, fortranCompat());
    if (hitDipole(xt, yt)) return stop(kD1In);

    // Go to D1 OUT magnetic boundary.
    eng.transp(spectr, 11, t, 526.053e0, pathlen, &d);
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, 6.0e0, xt, yt, fortranCompat());
    if (hitDipole(xt, yt)) return stop(kD1Out);

    // Vacuum pipes following the dipole. First the odd piece interfacing with the dipole itself (mc_hms.f:363)
    {
        const double ddx = xt - x_offset_pipes, ddy = yt - y_offset_pipes;
        if (((ddx * ddx + ddy * ddy) > 30.48 * 30.48) || (std::fabs(yt - y_offset_pipes) > 20.5232))
            return stop(kD1Out);
    }

    // Check the exit of the 26.65 inch pipe
    zdrift = z_dip1;
    optics::project(t, zdrift, pathlen, &d);
    {
        const double ddx = t.xs - x_offset_pipes, ddy = t.ys - y_offset_pipes;
        if ((ddx * ddx + ddy * ddy) > 1145.518) return stop(kD1Out);
    }

    // check exit of long (117 inch) pipe
    zdrift = z_dip2 - z_dip1;
    optics::project(t, zdrift, pathlen, &d);
    {
        const double ddx = t.xs - x_offset_pipes, ddy = t.ys - y_offset_pipes;
        if ((ddx * ddx + ddy * ddy) > 1512.2299) return stop(kD1Out);
    }

    // lastly check the exit of the last piece of pipe. 45.5 inches long, 30.62 dia.
    zdrift = z_dip3 - z_dip2;
    optics::project(t, zdrift, pathlen, &d);
    {
        const double ddx = t.xs - x_offset_pipes, ddy = t.ys - y_offset_pipes;
        if ((ddx * ddx + ddy * ddy) > 2162.9383) return stop(kD1Out);
    }

    // If we get this far, the particle is in the hut.
    cnt_[kHut] = cnt_[kHut] + 1;
    det.reached_hut = true;

    // and track through the detector hut
    if (!eng.adrift(spectr, 12) && log_) *log_ << " Transformation #12 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 12) - z_dip3;   // distance left to go
    const int hutCode = mcHmsHut(d, a, -zdrift, rng, det);
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

}  // namespace hms
}  // namespace simc

