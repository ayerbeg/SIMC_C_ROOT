// src/spectrometers/shms/ShmsSpectrometer.cpp — see include/simc/spectrometers/shms/ShmsSpectrometer.h.
// mc_shms (mc_shms.f), the ISpectrometer glue and the registration. Collimator: McShmsColl.cpp, hut: McShmsHut.cpp.
#include "simc/spectrometers/shms/ShmsSpectrometer.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "simc/core/Registry.h"
#include "simc/core/SimcParameters.h"
#include "simc/stages/IRandom.h"
#include "ShmsHutConstants.h"

#ifndef SIMC_SHMS_DEFAULT_MATRIX_DIR
#define SIMC_SHMS_DEFAULT_MATRIX_DIR "data/matrices/simc_fortran"
#endif

namespace simc {
namespace shms {

namespace {

// apertures_shms.inc
constexpr double r_HBx = 11.2;
constexpr double r_HBfym = -4.13, r_HBfyp = 11.75;
constexpr double r_HBmenym = -5.45, r_HBmenyp = 11.74;
constexpr double r_HBmexym = -10.25, r_HBmexyp = 11.71;
constexpr double r_HBbym = -11.71, r_HBbyp = 11.70;
constexpr double r_Q1 = 20.00;
constexpr double r_Q2 = 30.00;
constexpr double r_Q3 = 30.00;
constexpr double r_D1 = 30.00;

// mc_shms.f:89-94 collimator, 176-210 drift lengths (2017 MEs)
constexpr double h_entr = 8.5;
constexpr double v_entr = 12.5;
constexpr double h_exit = 8.65;
constexpr double v_exit = 12.85;
constexpr double x_off = +0.00;
constexpr double y_off = +0.00;
constexpr double zd_hbin = 118.39;
constexpr double zd_hbmen = 17.61;
constexpr double zd_hbmex = 80.0;
constexpr double zd_hbout = 17.61;
constexpr double z_entr = 25.189;   // 82.789 cm from center of HB
constexpr double z_thick = 6.35;
constexpr double zd_q1in = 58.39;
constexpr double zd_q1men = 28.35;
constexpr double zd_q1mid = 93.65;
constexpr double zd_q1mex = 93.65;
constexpr double zd_q1out = 28.35;
constexpr double zd_q2in = 25.55;
constexpr double zd_q2men = 39.1;
constexpr double zd_q2mid = 79.35;
constexpr double zd_q2mex = 79.35;
constexpr double zd_q2out = 39.1;
constexpr double zd_q3in = 28.10;
constexpr double zd_q3men = 39.1;
constexpr double zd_q3mid = 79.35;
constexpr double zd_q3mex = 79.35;
constexpr double zd_q3out = 39.1;
constexpr double zd_q3d1trans = 18.00;
constexpr double zd_d1flare = 30.10;
constexpr double zd_d1men = 39.47;   // FORTRAN-QUIRK (D-015 F5): COSY class 22 drift is 39.4715
constexpr double zd_d1mid = 36.406263;   // zd_d1mid1..7
constexpr double zd_d1mex = 36.406263;
constexpr double zd_d1out = 60.68;   // COSY class 31: 60.68058
constexpr double zd_fp = 307.95;     // COSY class 32: 307.9514
constexpr double pi = 3.141592653589793;   // constants.inc

const char* const kNamesArr[kNumSlots + 1] = {
    "",        "trials",  "successes", "targ_hor", "targ_vert", "targ_oct", "slit_hor", "slit_vert", "slit_oct",
    "hb_in",   "hb_men",  "hb_mex",    "hb_out",   "q1_in",     "q1_men",   "q1_mid",   "q1_mex",    "q1_out",
    "q2_in",   "q2_men",  "q2_mid",    "q2_mex",   "q2_out",    "q3_in",    "q3_men",   "q3_mid",    "q3_mex",
    "q3_out",  "D1_in",   "D1_flr",    "D1_men",   "D1_mid1",   "D1_mid2",  "D1_mid3",  "D1_mid4",   "D1_mid5",
    "D1_mid6", "D1_mid7", "D1_mex",    "D1_out",   "BP_in",     "BP_out",   "hut",      "dc1",       "dc2",
    "cer1",    "cer2",    "s1",        "s2",       "s3",        "cal",      "stop_id",  "coll"};

// indices into ShmsTrace
enum Tr : int {
    tHbIn, tHbMen, tHbMex, tHbOut, tQ1In, tQ1Men, tQ1Mid, tQ1Mex, tQ1Out, tQ2In, tQ2Men, tQ2Mid, tQ2Mex, tQ2Out,
    tQ3In, tQ3Men, tQ3Mid, tQ3Mex, tQ3Out, tDIn, tDFlr, tDMen, tDM1, tDMex = tDM1 + 7, tDOut
};

}  // namespace

const char* const ShmsTrace::kNames[ShmsTrace::kN] = {
    "hb_in", "hb_men", "hb_mex", "hb_out", "q1_in", "q1_men", "q1_mid", "q1_mex", "q1_out", "q2_in", "q2_men",
    "q2_mid", "q2_mex", "q2_out", "q3_in", "q3_men", "q3_mid", "q3_mex", "q3_out", "d_in", "d_flr", "d_men",
    "d_m1", "d_m2", "d_m3", "d_m4", "d_m5", "d_m6", "d_m7", "d_mex", "d_out"};

const char* counterName(int idx) {
    if (idx < 1 || idx > kNumSlots) throw std::out_of_range("shms::counterName: bad index");
    return kNamesArr[idx];
}

std::string ShmsSpectrometer::defaultMatrixDir() { return SIMC_SHMS_DEFAULT_MATRIX_DIR; }

ShmsSpectrometer::ShmsSpectrometer()
    : eng_(std::make_shared<optics::CosyEngine>()), matrixDir_(defaultMatrixDir()), log_(&std::cout) {}

void ShmsSpectrometer::setLog(std::ostream* os) {
    log_ = os;
    mat_.setLog(os);
    eng_->setLog(os);
}

void ShmsSpectrometer::configure(const Json& config, const SimcParameters& params) {
    usingColl_ = params.using_SHMScoll;
    ctau_ = params.ctau;
    if (!config.is_object()) return;
    if (config.contains("matrix_dir")) matrixDir_ = config.at("matrix_dir").get<std::string>();
    if (config.contains("spectr")) spectr_ = config.at("spectr").get<int>();
    if (config.contains("using_SHMScoll")) usingColl_ = config.at("using_SHMScoll").get<bool>();
    if (config.contains("ctau")) ctau_ = config.at("ctau").get<double>();
    if (config.contains("use_sieve")) useSieve_ = config.at("use_sieve").get<bool>();
    if (config.contains("use_2sieve")) use2Sieve_ = config.at("use_2sieve").get<bool>();
    if (config.contains("use_coll")) useColl_ = config.at("use_coll").get<bool>();
    if (config.contains("skip_hb")) skipHb_ = config.at("skip_hb").get<bool>();
    if (config.contains("ms_flag")) msOverride_ = config.at("ms_flag").get<bool>();
    if (config.contains("wcs_flag")) wcsOverride_ = config.at("wcs_flag").get<bool>();
    if (config.contains("decay")) decayOverride_ = config.at("decay").get<bool>();
    if (config.contains("verbose") && !config.at("verbose").get<bool>()) setLog(nullptr);
}

void ShmsSpectrometer::init(const StageContext&) {
    ensureForward(spectr_);
    ensureRecon(spectr_);
}

// Fortran: mc_shms.f:335-347 (first_time_here)
void ShmsSpectrometer::ensureForward(int spectr) {
    if (spectr != 5) {
        if (log_) *log_ << " LSA tune no longer in use!\n";
        throw optics::FortranStop("mc_shms: spectr = " + std::to_string(spectr) + " (LSA tune no longer in use!)");
    }
    if (eng_->hasForward(spectr)) return;
    const int n_classes = eng_->initForward(spectr, matrixDir_);
    if (n_classes != 32) throw optics::FortranStop("Bender-SHMS, wrong number of transport classes");
}

// Fortran: mc_shms_recon.f firsttime block
void ShmsSpectrometer::ensureRecon(int spectr) {
    if (!eng_->hasRecon(spectr)) eng_->initRecon(spectr, matrixDir_);
}

const std::map<std::string, long>& ShmsSpectrometer::stopCounters() const {
    cntMap_.clear();
    for (int i = 1; i <= kNumSlots; ++i)
        if (i != kStopIdSlot) cntMap_[kNamesArr[i]] = cnt_[i];
    cntMap_["q2_mid_dup"] = cnt_[kQ2Mid];   // simc.f:519-520 prints q2_mid in the Q3 'mid' slot
    return cntMap_;
}

// Fortran: mc_shms.f:mc_shms (l.1-1108)
McShmsDetail ShmsSpectrometer::mcShms(McShmsArgs& a, IRandom& rng) {
    McShmsDetail det;
    cerFlag_ = true;    // mc_shms.f:257-258
    vacFlag_ = false;
    cnt_[kStopIdSlot] = 0;   // stop_id = 0
    a.ok_spec = false;
    optics::DecayState d;
    d.decay_flag = a.decay_flag;
    d.dflag = false;
    d.m2 = a.m2;
    d.ctau = a.ctau;
    d.decdist = a.decdist;
    d.Mh2_final = a.Mh2_final;
    d.rng = &rng;
    cnt_[kTrials] = cnt_[kTrials] + 1;
    det.stop_code = body(a, d, rng, det);
    det.stop_id = static_cast<int>(cnt_[kStopIdSlot]);
    a.m2 = d.m2;
    a.decdist = d.decdist;
    a.Mh2_final = d.Mh2_final;
    return det;
}

int ShmsSpectrometer::body(McShmsArgs& a, optics::DecayState& d, IRandom& rng, McShmsDetail& det) {
    optics::Track& t = track_;
    ShmsTrace& tr = trace_;
    double xt, yt, zdrift;
    xt = a.th_spec;   // avoid 'unused variable' error for th_spec (mc_shms.f:265)
    // spec(1:58) = -100000. (spec_ntuple): the spec array is never output — not ported.

    // Sieve in front of HB (mc_shms.f:276-288)
    if (use2Sieve_) {
        const double z_slit = zd_hbin;
        xt = a.x + z_slit * a.dxdz;
        yt = a.y + z_slit * a.dydz;
        xt = 2.5 * std::round(xt / 2.5);   // nint: half away from zero
        yt = 1.64 * std::round(yt / 1.64);
        const double rt = 0.3 * std::sqrt(rng.uniform());
        const double tht = 2 * pi * rng.uniform();
        xt = xt + rt * std::cos(tht);
        yt = yt + rt * std::sin(tht);
        a.dxdz = (xt - a.x) / z_slit;
        a.dydz = (yt - a.y) / z_slit;
    }
    // Mock sieve slit at the collimator (mc_shms.f:293-315)
    if (useSieve_) {
        const double z_slit = zd_hbin + zd_hbmen + zd_hbmex + zd_hbout + z_entr;
        xt = a.x + z_slit * a.dxdz;
        yt = a.y + z_slit * a.dydz;
        xt = 2.5 * std::round(xt / 2.5);
        yt = 1.64 * std::round(yt / 1.64);
        const double rt = 0.3 * std::sqrt(rng.uniform());
        const double tht = 2 * pi * rng.uniform();
        xt = xt + rt * std::cos(tht);
        yt = yt + rt * std::sin(tht);
        a.dxdz = (xt - a.x) / z_slit;
        a.dydz = (yt - a.y) / z_slit;
    }

    // particle momentum
    t.dpps = a.dpp;
    d.ph = a.p_spec * (1. + t.dpps / 100.);

    t.xs = a.x;
    t.ys = a.y;
    t.zs = a.z;
    t.dxdzs = a.dxdz;
    t.dydzs = a.dydz;

    // Read in transport coefficients (first_time_here; spectr must be 5).
    const int spectr = a.spectr;
    ensureForward(spectr);
    optics::CosyEngine& eng = *eng_;
    double& pathlen = a.pathlen;
    auto project = [&](double z) { optics::project(t, z, pathlen, &d); };
    auto transp = [&](int cls, double zd) { eng.transp(spectr, cls, t, zd, pathlen, &d); };

    // HB (mc_shms.f:414-498)
    if (skipHb_) {
        zdrift = zd_hbin + zd_hbmen + zd_hbmex + zd_hbout;
        project(zdrift);
    } else {
        // Go to HB Mech entrance.
        zdrift = zd_hbin;
        project(zdrift);
        xt = t.xs;
        yt = t.ys;
        optics::rotateVaxis(t, 1.5, xt, yt);
        yt = yt + 1.51;
        tr.x[tHbIn] = xt;
        tr.y[tHbIn] = yt;
        if ((xt * xt > r_HBx * r_HBx) || (yt > r_HBfyp) || (yt < r_HBfym)) return stop(kHBIn, 1);

        // Go to HB Mag entrance
        zdrift = zd_hbmen;
        project(zdrift);
        xt = t.xs;
        yt = t.ys;
        optics::rotateVaxis(t, 1.5, xt, yt);
        yt = yt + 0.98;
        tr.x[tHbMen] = xt;
        tr.y[tHbMen] = yt;
        if ((xt * xt > r_HBx * r_HBx) || (yt > r_HBmenyp) || (yt < r_HBmenym)) return stop(kHBMen, 2);

        // Go to HB Mag exit
        transp(3, zd_hbmex);
        xt = t.xs;
        yt = t.ys;
        optics::rotateVaxis(t, -1.5, xt, yt);
        yt = yt + 0.98;
        tr.x[tHbMex] = xt;
        tr.y[tHbMex] = yt;
        if ((xt * xt > r_HBx * r_HBx) || (yt > r_HBmexyp) || (yt < r_HBmexym)) return stop(kHBMex, 3);

        // Go to HB exit
        zdrift = zd_hbout;
        project(zdrift);
        xt = t.xs;
        yt = t.ys;
        optics::rotateVaxis(t, -1.5, xt, yt);
        yt = yt + 1.51;
        tr.x[tHbOut] = xt;
        tr.y[tHbOut] = yt;
        if ((xt * xt > r_HBx * r_HBx) || (yt > r_HBbyp) || (yt < r_HBbym)) return stop(kHBOut, 4);
    }

    // GH collimator routine including pion/muon scattering (mc_shms.f:503-568)
    if (a.using_SHMScoll && (d.m2 > 100.0 * 100.0) && (d.m2 < 200.0 * 200.0)) {
        zdrift = z_entr;
        project(zdrift);
        const bool Scollsuccess = mcShmsColl(d, a.p_spec, pathlen, rng);
        if (!Scollsuccess) {
            // save_Scollsuccess = 1.5 (COMMON /save_coll/, never read)
            return stop(kColl);   // stop_id stays 0
        }
    } else {
        if (useColl_) {
            zdrift = z_entr;
            project(zdrift);
            xt = t.xs;
            yt = t.ys;
            // FORTRAN-QUIRK: the collimator stops leave stop_id = 0
            if (std::fabs(yt - y_off) > h_entr) return stop(kSlitHor);
            if (std::fabs(xt - x_off) > v_entr) return stop(kSlitVert);
            if (std::fabs(xt - x_off) > (-v_entr / h_entr * std::fabs(yt - y_off) + 3 * v_entr / 2))
                return stop(kSlitOct);
            zdrift = z_thick;
            project(zdrift);
            xt = t.xs;
            yt = t.ys;
            if (std::fabs(yt - y_off) > (h_exit)) return stop(kSlitHor);
            if (std::fabs(xt - x_off) > (v_exit)) return stop(kSlitVert);
            if (std::fabs(xt - x_off) > ((-v_exit) / (h_exit)*std::fabs(yt - y_off) + 3 * (v_exit) / 2))
                return stop(kSlitOct);
        }
    }

    auto circ = [&](double r) { return (t.xs * t.xs + t.ys * t.ys) > r * r; };
    auto rec = [&](int i) {
        tr.x[i] = t.xs;
        tr.y[i] = t.ys;
    };

    // Q1 (mc_shms.f:573-638)
    if (useColl_) {
        zdrift = zd_q1in - z_entr - z_thick;
    } else {
        zdrift = zd_q1in;
    }
    project(zdrift);
    rec(tQ1In);
    if (circ(r_Q1)) return stop(kQ1In, 5);
    zdrift = zd_q1men;
    project(zdrift);
    rec(tQ1Men);
    if (circ(r_Q1)) return stop(kQ1Men, 6);
    transp(7, zd_q1mid);
    rec(tQ1Mid);
    if (circ(r_Q1)) return stop(kQ1Mid, 7);
    transp(8, zd_q1mex);
    rec(tQ1Mex);
    if (circ(r_Q1)) return stop(kQ1Mex, 8);
    zdrift = zd_q1out;
    project(zdrift);
    rec(tQ1Out);
    if (circ(r_Q1)) return stop(kQ1Out, 9);

    // Q2 (mc_shms.f:641-701)
    zdrift = zd_q2in;
    project(zdrift);
    rec(tQ2In);
    if (circ(r_Q2)) return stop(kQ2In, 10);
    zdrift = zd_q2men;
    project(zdrift);
    rec(tQ2Men);
    if (circ(r_Q2)) return stop(kQ2Men, 11);
    transp(12, zd_q2mid);
    rec(tQ2Mid);
    if (circ(r_Q2)) return stop(kQ2Mid, 12);
    transp(13, zd_q2mex);
    rec(tQ2Mex);
    if (circ(r_Q2)) return stop(kQ2Mex, 13);
    zdrift = zd_q2out;
    project(zdrift);
    rec(tQ2Out);
    if (circ(r_Q2)) return stop(kQ2Out, 14);

    // Q3 (mc_shms.f:704-764)
    zdrift = zd_q3in;
    project(zdrift);
    rec(tQ3In);
    if (circ(r_Q3)) return stop(kQ3In, 15);
    zdrift = zd_q3men;
    project(zdrift);
    rec(tQ3Men);
    if (circ(r_Q3)) return stop(kQ3Men, 16);
    transp(17, zd_q3mid);
    rec(tQ3Mid);
    if (circ(r_Q3)) return stop(kQ3Mid, 17);
    transp(18, zd_q3mex);
    rec(tQ3Mex);
    if (circ(r_Q3)) return stop(kQ3Mex, 18);
    zdrift = zd_q3out;
    project(zdrift);
    rec(tQ3Out);
    if (circ(r_Q3)) return stop(kQ3Out, 19);

    // D1 mechanical entrance (mc_shms.f:815-825)
    zdrift = zd_q3d1trans;
    project(zdrift);
    rec(tDIn);
    if (circ(r_D1)) return stop(kD1In, 20);

    // Rotated/offset apertures: xt,yt = rotate_haxis(ang) of (xs,ys), then xt += off.
    auto dip = [&](double ang, double off, int i) {
        xt = t.xs;
        yt = t.ys;
        optics::rotateHaxis(t, ang, xt, yt);
        xt = xt + off;
        tr.x[i] = xt;
        tr.y[i] = yt;
        return (xt * xt + yt * yt) > r_D1 * r_D1;
    };

    // D1 entrance flare (mc_shms.f:832-846): xt = xt - 3.5
    zdrift = zd_d1flare;
    project(zdrift);
    if (dip(9.200, -3.5, tDFlr)) return stop(kD1Flr, 21);
    // D1 magnetic entrance
    zdrift = zd_d1men;
    project(zdrift);
    if (dip(9.200, 2.82, tDMen)) return stop(kD1Men, 22);
    // dipole interior apertures 1..7 (mc_shms.f:874-993)
    transp(23, zd_d1mid);
    if (dip(6.9, 8.05, tDM1)) return stop(kD1Mid1, 23);
    transp(24, zd_d1mid);
    if (dip(4.6, 11.75, tDM1 + 1)) return stop(kD1Mid1 + 1, 24);
    transp(25, zd_d1mid);
    if (dip(2.3, 13.96, tDM1 + 2)) return stop(kD1Mid1 + 2, 25);
    transp(26, zd_d1mid);
    if (dip(0.0, 14.70, tDM1 + 3)) return stop(kD1Mid1 + 3, 26);
    transp(27, zd_d1mid);
    if (dip(-2.3, 13.96, tDM1 + 4)) return stop(kD1Mid1 + 4, 27);
    transp(28, zd_d1mid);
    if (dip(-4.6, 11.75, tDM1 + 5)) return stop(kD1Mid1 + 5, 28);
    transp(29, zd_d1mid);
    if (dip(-6.9, 8.05, tDM1 + 6)) return stop(kD1Mid1 + 6, 29);
    // D1 magnetic exit
    transp(30, zd_d1mex);
    if (dip(-9.2, 2.82, tDMex)) return stop(kD1Mex, 30);
    // D1 OUT mechanical boundary: xt = xt - 6.88
    zdrift = zd_d1out;
    project(zdrift);
    if (dip(-9.20, -6.88, tDOut)) return stop(kD1Out, 31);

    // Transport to the hut (not all the way to the focal plane), mc_shms.f:1037-1047
    if (cerFlag_) {
        zdrift = zd_fp + hut::hcer_1_zentrance;   // 307.95 cm to fp - 291.7 cm to Cer entrance
    } else {
        if (vacFlag_) {
            zdrift = zd_fp + hut::hcer_1_zexit;
        } else {
            zdrift = zd_fp + hut::hcer_1_zentrance;
        }
    }
    project(zdrift);

    // If we get this far, the particle is in the hut.
    cnt_[kHut] = cnt_[kHut] + 1;
    det.reached_hut = true;
    tr.xs_fp = t.xs;
    tr.ys_fp = t.ys;
    tr.dxdzs_fp = t.dxdzs;
    tr.dydzs_fp = t.dydzs;

    // and track through the detector hut (zinit = 0.0, unused by mc_shms_hut)
    const int hutCode = mcShmsHut(d, a, 0.0, rng, det);
    if (hutCode != kSuccess) return hutCode;

    t.xs = a.x_fp;
    t.ys = a.y_fp;
    t.dxdzs = a.dx_fp;
    t.dydzs = a.dy_fp;

    // Reconstruct target quantities.
    ensureRecon(spectr);
    double dpp_recon, dth_recon, dph_recon, y_recon;
    eng.recon(spectr, t, a.fry, dpp_recon, dth_recon, dph_recon, y_recon);

    a.dpp = dpp_recon;
    a.dxdz = dph_recon;
    a.dydz = dth_recon;
    a.y = y_recon;
    a.ok_spec = true;
    cnt_[kSuccesses] = cnt_[kSuccesses] + 1;
    return kSuccess;
}

// ISpectrometer glue: simc.f:montecarlo l.1482-1489 (hadron arm) / l.1747-1755 (electron arm).
SpecResult ShmsSpectrometer::transport(const SpecInput& in, IRandom& rng) {
    McShmsArgs a;
    a.p_spec = in.p_spec;
    a.th_spec = in.th_spec;
    a.dpp = in.delta;
    a.x = in.x;
    a.y = in.y;
    a.z = in.z;
    a.dxdz = in.dxdz;
    a.dydz = in.dydz;
    a.m2 = in.m2;
    a.ms_flag = msOverride_ ? *msOverride_ : in.ms_flag;
    a.wcs_flag = wcsOverride_ ? *wcsOverride_ : in.wcs_flag;
    a.decay_flag = decayOverride_ ? *decayOverride_ : in.decay_flag;
    a.fry = in.fry;
    a.pathlen = 0.0;
    a.spectr = spectr_;
    a.using_SHMScoll = usingColl_;
    a.ctau = ctau_;
    a.decdist = 0.0;
    a.Mh2_final = in.m2;
    const McShmsDetail det = mcShms(a, rng);

    SpecResult r;
    r.ok = a.ok_spec;
    r.stop_reason = det.stop_code == kSuccess ? std::string() : std::string(counterName(det.stop_code));
    if (det.fit_done) {
        r.x_fp = a.x_fp;
        r.dx_fp = a.dx_fp;
        r.y_fp = a.y_fp;
        r.dy_fp = a.dy_fp;
    }
    if (r.ok) {
        r.delta_recon = a.dpp;
        r.dxdz_recon = a.dxdz;
        r.dydz_recon = a.dydz;
        r.y_recon = a.y;
    }
    r.pathlen = a.pathlen;
    r.resmult = det.reached_hut ? a.resmult : 0.0;
    r.m2_final = a.m2;
    r.resmult_written = det.reached_hut;
    r.decdist = a.decdist;
    r.Mh2_final = a.Mh2_final;
    return r;
}

}  // namespace shms
}  // namespace simc

SIMC_REGISTER(simc::ISpectrometer, "shms", simc::shms::ShmsSpectrometer)
