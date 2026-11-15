// src/spectrometers/hrs/McHrs.cpp — see include/simc/spectrometers/hrs/kernel/McHrs.h (kernel).
// mc_hrsl (hrsl/mc_hrsl.f) and mc_hrsr (hrsr/mc_hrsr.f), one implementation with a side parameter. The detector hut is in
// McHrsHut.cpp, the ISpectrometer plug-ins in HrsSpectrometer.cpp.
#include "simc/spectrometers/hrs/kernel/McHrs.h"

#include "simc/core/DataPath.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include "simc/stages/IRandom.h"


namespace simc {
namespace hrs {

namespace {

using quad = Quad;

// apertures_hrs*.inc (identical for both sides)
constexpr double r_Q1 = 15.0;    // 15.0 cm radius
constexpr double r_Q2 = 30.22;   // 30.22 cm radius
constexpr double r_Q3 = 30.22;

// mc_hrs*.f:56-63: large collimator (h_exit differs per side, HrsSideParams)
constexpr double h_entr = 3.145;
constexpr double v_entr = 6.090;
constexpr double v_exit = 6.485;
constexpr double y_off = 0.0;
constexpr double x_off = 0.0;
constexpr double z_off = 0.0;

constexpr double pi = 3.141592653589793;   // constants.inc (mock sieve only)

// z_entr = <lit>d0 + z_off, z_exit = z_entr + 8.0d0: REAL(16) parameter expressions rounded to the real*8 PARAMETERs.
double zEntr(unsigned long long tenths) {
    const quad q = static_cast<quad>(tenths) / static_cast<quad>(10) + static_cast<quad>(z_off);
    return static_cast<double>(q);
}
double zExit(double z_entr) { return static_cast<double>(static_cast<quad>(z_entr) + static_cast<quad>(8)); }

// a decimal literal m * 10^-e as the binary128 constant (exact numerator and denominator, one correctly rounded division)
quad quadLiteral(unsigned long long mant, int e) {
    quad den = 1;
    for (int i = 0; i < e; ++i) den = den * static_cast<quad>(10);
    return static_cast<quad>(mant) / den;
}

const char* const kNames[kNumCounters + 1] = {"",       "trials", "slit_hor", "slit_vert", "Q1_in",  "Q1_mid",
                                              "Q1_out", "Q2_in",  "Q2_mid",   "Q2_out",    "Q3_in",  "Q3_mid",
                                              "Q3_out", "D1_in",  "D1_out",   "hut",       "dc1",    "dc2",
                                              "s1",     "s2",     "cal",      "successes"};

}  // namespace

double low8(quad q) {
    static_assert(sizeof(quad) == 16, "binary128 expected");
    double d;
    std::memcpy(&d, &q, sizeof d);   // x86-64 little endian: the low 64 mantissa bits
    return d;
}

const HrsSideParams& sideParams(Side s) {
    // mc_hrsr.f:58, 67, 524 / mc_hrsl.f:58, 67, 525
    static const HrsSideParams R{3, 3.340, zEntr(1100), zExit(zEntr(1100)), 0.48, "hrsr"};
    static const HrsSideParams L{4, 3.335, zEntr(1109), zExit(zEntr(1109)), 0.78, "hrsl"};
    return s == Side::Right ? R : L;
}

const std::array<TranspLength, 7>& transpLengths() {
    // mc_hrs*.f:248, 259, 307, 318, 388, 455, 466
    static const std::array<TranspLength, 7> t = {{
        {2, 62.75333333, low8(quadLiteral(6275333333ULL, 8))},
        {3, 31.37666667, low8(quadLiteral(3137666667ULL, 8))},
        {5, 121.77333333, low8(quadLiteral(12177333333ULL, 8))},
        {6, 60.88666667, low8(quadLiteral(6088666667ULL, 8))},
        {8, 659.73445725, low8(quadLiteral(65973445725ULL, 8))},
        {10, 121.7866667, low8(quadLiteral(1217866667ULL, 7))},
        {11, 60.89333333, low8(quadLiteral(6089333333ULL, 8))},
    }};
    return t;
}

const char* counterName(int idx) {
    if (idx < 1 || idx > kNumCounters) throw std::out_of_range("hrs::counterName: bad index");
    return kNames[idx];
}

std::string McHrs::defaultMatrixDir() { return simc::DataDir("matrices/simc_fortran"); }

McHrs::McHrs(Side side)
    : side_(side), eng_(std::make_shared<optics::CosyEngine>()), matrixDir_(defaultMatrixDir()), log_(&std::cout) {}

void McHrs::setLog(std::ostream* os) {
    log_ = os;
    mat_.setLog(os);
    eng_->setLog(os);
}

double McHrs::transpLength(int cls) const {
    for (const auto& t : transpLengths())
        if (t.cls == cls) return fortranCompat() ? t.compat : t.intended;
    throw std::out_of_range("McHrs::transpLength: class without a transp call");
}

// Fortran: mc_hrs*.f:148-153 (first_time_here: transp_init, 12 classes). The single repository copy of the (identical)
// hrsl/ and hrsr/ files is read for both sides.
void McHrs::ensureForward() {
    const int spectr = params().spectr;
    if (eng_->hasForward(spectr)) return;
    const int n_classes = eng_->loadForward(spectr, matrixDir_ + "/hrs/hrs_forward_cosy.dat");
    if (n_classes != 12)
        throw simc::FortranStop(side_ == Side::Left ? "MC_HRSL, wrong number of transport classes"
                                                    : "MC_HRSR, wrong number of transport classes");
}

// Fortran: mc_hrs*_recon.f firsttime block
void McHrs::ensureRecon() {
    const int spectr = params().spectr;
    if (!eng_->hasRecon(spectr)) eng_->loadRecon(spectr, matrixDir_ + "/hrs/hrs_recon_cosy.dat");
}

// Fortran: mc_hrs*.f (l.1-551)
McHrsDetail McHrs::mcHrs(McHrsArgs& a, IRandom& rng) {
    McHrsDetail det;
    a.ok_spec = false;
    optics::DecayState d;
    d.decay_flag = a.decay_flag;
    d.dflag = false;   // particle has not decayed yet
    d.m2 = a.m2;
    d.ctau = a.ctau;
    d.decdist = a.decdist;
    d.Mh2_final = a.Mh2_final;
    d.rng = &rng;
    d.fortran_compat = fortranCompat();
    d.trace = transportTrace_;
    if (transportTrace_) transportTrace_->reset();
    cnt_[kTrials] = cnt_[kTrials] + 1;
    det.stop_code = body(a, d, rng, det);
    a.m2 = d.m2;
    a.decdist = d.decdist;
    a.Mh2_final = d.Mh2_final;
    return det;
}

int McHrs::body(McHrsArgs& a, optics::DecayState& d, IRandom& rng, McHrsDetail& det) {
    const HrsSideParams& P = params();
    const int spectr = P.spectr;
    const double h_exit = P.h_exit, z_entr = P.z_entr, z_exit = P.z_exit;
    optics::Track& t = track_;
    double xt, yt;
    xt = a.th_spec;   // avoid 'unused variable' error for th_spec (mc_hrs*.f:117)

    // mock sieve slit (mc_hrs*.f:122-131). FORTRAN-QUIRK: rt and tht are drawn but never applied to xt, yt.
    if (useSieve_) {
        xt = a.x + z_exit * a.dxdz;
        yt = a.y + z_exit * a.dydz;
        xt = 2.50 * std::round(xt / 2.50);
        yt = 1.25 * std::round(yt / 1.25);
        const double rt = 0.1 * std::sqrt(rng.uniform());
        const double tht = 2 * pi * rng.uniform();
        (void)rt;
        (void)tht;
        a.dxdz = (xt - a.x) / z_exit;
        a.dydz = (yt - a.y) / z_exit;
    }

    t.xs = a.x;
    t.ys = a.y;
    t.zs = a.z;
    t.dxdzs = a.dxdz;
    t.dydzs = a.dydz;
    t.dpps = a.dpp;
    d.ph = a.p_spec * (1. + t.dpps / 100.);

    ensureForward();
    optics::CosyEngine& eng = *eng_;
    double& pathlen = a.pathlen;
    double zdrift, ztmp;
    auto project = [&](double z) { optics::project(t, z, pathlen, &d); };
    auto rad = [&]() { return std::sqrt(t.xs * t.xs + t.ys * t.ys); };
    auto stopHere = [&](int code, double where) { return stop(code, where, t.xs, t.ys); };

    // Circular apertures before slitbox (only important for no collimator), l.159-180
    zdrift = 65.686;
    ztmp = zdrift;
    project(zdrift);
    if (rad() > 7.3787) return stopHere(kSlitHor, 20.);
    zdrift = 80.436 - ztmp;
    ztmp = 80.436;
    project(zdrift);
    if (rad() > 7.4092) return stopHere(kSlitHor, 21.);

    // Check front of fixed slit (l.184-200).
    zdrift = z_entr - ztmp;
    project(zdrift);
    if (std::fabs(t.ys - y_off) > h_entr) return stopHere(kSlitHor, 1.);
    if (std::fabs(t.xs - x_off) > v_entr) return stopHere(kSlitVert, 2.);

    // Check back of fixed slit (l.203-218).
    zdrift = z_exit - z_entr;
    project(zdrift);
    if (std::fabs(t.ys - y_off) > h_exit) return stopHere(kSlitHor, 3.);
    if (std::fabs(t.xs - x_off) > v_exit) return stopHere(kSlitVert, 4.);

    // Aperture before Q1 (l.222-231)
    ztmp = 135.064;
    zdrift = ztmp - z_exit;
    project(zdrift);
    if (rad() > 12.5222) return stopHere(kQ1In, 22.);

    // Go to Q1 IN mag bound (l.235-244)
    if (!eng.adrift(spectr, 1) && log_) *log_ << " Transformation #1 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 1) - ztmp;
    project(zdrift);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q1 * r_Q1) return stopHere(kQ1In, 5.);

    // Check aperture at 2/3 of Q1 / Q1 OUT (l.248-266). FORTRAN-BUG(D-061 F1): REAL(16) length literals.
    eng.transp(spectr, 2, t, transpLength(2), pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q1 * r_Q1) return stopHere(kQ1Mid, 6.);
    eng.transp(spectr, 3, t, transpLength(3), pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q1 * r_Q1) return stopHere(kQ1Out, 7.);

    // Apertures after Q1, before Q2 (l.270-290)
    zdrift = 300.464 - 253.16;   // Q1 exit is z=253.16
    ztmp = zdrift;               // distance from Q1 exit
    project(zdrift);
    if (rad() > 14.9225) return stopHere(kQ1Out, 23.);
    zdrift = 314.464 - 300.464;
    ztmp = ztmp + zdrift;
    project(zdrift);
    if (rad() > 20.9550) return stopHere(kQ2In, 24.);

    // Go to Q2 IN mag bound (l.294-303)
    if (!eng.adrift(spectr, 4) && log_) *log_ << " Transformation #4 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 4) - ztmp;
    project(zdrift);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q2 * r_Q2) return stopHere(kQ2In, 8.);

    // Q2 2/3 and OUT (l.307-325)
    eng.transp(spectr, 5, t, transpLength(5), pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q2 * r_Q2) return stopHere(kQ2Mid, 9.);
    eng.transp(spectr, 6, t, transpLength(6), pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q2 * r_Q2) return stopHere(kQ2Out, 10.);

    // Apertures after Q2, before D1 (l.329-360)
    zdrift = 609.664 - 553.020;   // Q2 exit is z=553.02
    ztmp = zdrift;
    project(zdrift);
    if (rad() > 30.0073) return stopHere(kQ2Out, 25.);
    zdrift = 641.800 - 609.664;
    ztmp = ztmp + zdrift;
    project(zdrift);
    if (rad() > 30.0073) return stopHere(kQ2Out, 26.);
    zdrift = 819.489 - 641.800;
    ztmp = ztmp + zdrift;
    project(zdrift);
    if (std::fabs(t.xs) > 50.0 || std::fabs(t.ys) > 15.0) return stopHere(kD1In, 27.);

    // Go to D1 IN magnetic boundary, intersection with the rotated aperture plane (l.364-383)
    if (!eng.adrift(spectr, 7) && log_) *log_ << " Transformation #7 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 7) - ztmp;
    project(zdrift);
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, -30.0, xt, yt, fortranCompat());
    if (std::fabs(xt - 2.500) > 52.5) return stop(kD1In, 11., xt, yt);   // -50 < x < +55
    if ((std::fabs(yt) + 0.01861 * xt) > 12.5) return stop(kD1In, 12., xt, yt);   // tan(1.066) ~ 0.01861

    // Go to D1 OUT magnetic boundary (l.388-408)
    eng.transp(spectr, 8, t, transpLength(8), pathlen, &d);
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, 30.0, xt, yt, fortranCompat());
    if (std::fabs(xt - 2.500) > 52.5) return stop(kD1Out, 13., xt, yt);
    if ((std::fabs(yt) + 0.01861 * xt) > 12.5) return stop(kD1Out, 14., xt, yt);

    // Apertures after D1, before Q3 (l.413-445)
    zdrift = 1745.33546 - 1655.83446;   // D1 exit is z=1655.83446
    ztmp = zdrift;
    project(zdrift);
    if (rad() > 30.3276) return stopHere(kD1Out, 28.);
    if (std::fabs(t.xs) > 50.0 || std::fabs(t.ys) > 15.0) return stopHere(kD1Out, 29.);
    zdrift = 1759.00946 - 1745.33546;
    ztmp = ztmp + zdrift;
    project(zdrift);
    if (rad() > 30.3276) return stopHere(kQ3In, 30.);

    // Go to Q3 IN mag bound (l.449-458)
    if (!eng.adrift(spectr, 9) && log_) *log_ << " Transformation #9 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 9) - ztmp;
    project(zdrift);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q3 * r_Q3) return stopHere(kQ3In, 15.);

    // Q3 2/3 and OUT (l.462-480)
    eng.transp(spectr, 10, t, transpLength(10), pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q3 * r_Q3) return stopHere(kQ3Mid, 16.);
    eng.transp(spectr, 11, t, transpLength(11), pathlen, &d);
    if ((t.xs * t.xs + t.ys * t.ys) > r_Q3 * r_Q3) return stopHere(kQ3Out, 17.);

    // Apertures after Q3 (l.484-493)
    zdrift = 2080.38746 - 1997.76446;   // Q3 exit is z=1997.76446
    ztmp = zdrift;
    project(zdrift);
    if (std::fabs(t.xs) > 35.56 || std::fabs(t.ys) > 17.145) return stopHere(kQ3Out, 31.);

    // Vacuum window is 15.522cm before FP (which is at VDC1), l.497-506
    zdrift = 2327.47246 - 2080.38746;
    ztmp = ztmp + zdrift;
    project(zdrift);
    if (std::fabs(t.xs) > 99.76635 || std::fabs(t.ys) > 17.145) return stopHere(kQ3Out, 32.);

    // If we get this far, the particle is in the hut.
    cnt_[kHut] = cnt_[kHut] + 1;
    det.reached_hut = true;

    if (!eng.adrift(spectr, 12) && log_) *log_ << " Transformation #12 is NOT a drift\n";
    zdrift = eng.driftdist(spectr, 12) - ztmp;   // distance left to go
    const int hutCode = mcHrsHut(d, a, -zdrift, rng, det);
    if (hutCode != kSuccess) return hutCode;

    // replace xs,ys,... with 'tracked' quantities.
    t.xs = a.x_fp;
    t.ys = a.y_fp;
    t.dxdzs = a.dx_fp;
    t.dydzs = a.dy_fp;

    // Apply offset to y_fp (detectors offset w.r.t optical axis), not to ys (recon), l.524-525.
    a.y_fp = a.y_fp - P.vdc_offset;

    ensureRecon();
    double dpp_recon, dth_recon, dph_recon, y_recon;
    eng.recon(spectr, t, a.fry, dpp_recon, dth_recon, dph_recon, y_recon);
    // (mc_hrs*.f:529-532: 'if (.not.ok)' after the recon tests the hut flag, always true here)

    a.dpp = dpp_recon;
    a.dxdz = dph_recon;
    a.dydz = dth_recon;
    a.y = y_recon;

    a.ok_spec = true;
    cnt_[kSuccesses] = cnt_[kSuccesses] + 1;
    return kSuccess;
}

}  // namespace hrs
}  // namespace simc
