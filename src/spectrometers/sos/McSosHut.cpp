// src/spectrometers/sos/McSosHut.cpp
// Fortran: sos/mc_sos_hut.f:mc_sos_hut (l.1-548) — the SOS detector hut: exit window, two drift chambers (6 planes
// each, wire-chamber smearing, CERNLIB lfit of the focal-plane track), hodoscopes (3/4 trigger), Cherenkov,
// calorimeter; multiple scattering in every layer if ms_flag.
// All PARAMETERs are copied verbatim; with -fdefault-real-8 they are real*8, and the parameter-only
// expressions are folded in double, like the constexpr ones here.
#include <cmath>

#include "simc/optics/Lfit.h"
#include "simc/spectrometers/sos/kernel/McSos.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace sos {

namespace {

// Math constants (mc_sos_hut.f:26-30; this routine has its own, shorter pi)
constexpr double pi = 3.141592654;
constexpr double d_r = pi / 180.;

// spectrometer exit window, 15 mil Kevlar, 5 mil Mylar (use 20 mil,X0=53.3) (l.43-48)
constexpr double sfoil_exit_radlen = 53.3;
constexpr double sfoil_exit_thick = 0.020 * 2.54;
constexpr double sfoil_exit_zpos = -3.22;
constexpr double sfoil_exit_ang = 0. * d_r;   // 0 deg. AFTER installation
// spectrometer air gaps
constexpr double sair_radlen = 30420.;
// chamber gas
constexpr double sdc_radlen = 16700.0;
constexpr double sdc_thick = 0.61775;
// effective wire plane: 60 micron field wires + 30 micron sense wires, Tungsten
constexpr double sdc_wire_radlen = 0.35;
constexpr double sdc_wire_thick = 0.0000354;
// chamber cathode foils, 1/2 mil of Mylar
constexpr double sdc_cath_radlen = 28.7;
constexpr double sdc_cath_thick = 0.0005 * 2.54;
// hodoscopes
constexpr double sscin_radlen = 42.4;
// Cherenkov entrance foil, 0.5 mm Al
constexpr double scer_entr_radlen = 8.90;
constexpr double scer_entr_thick = 0.050;
// Cherenkov, 1 atm Freon 12
constexpr double scer_radlen = 4810.0;
// Cherenkov mirror, 75 mu plus 2 cm Rotacell
constexpr double scer_mir_radlen = 400.0;
constexpr double scer_mir_thick = 2.0;
// Cherenkov exit foil, 0.5 mm Al
constexpr double scer_exit_radlen = 8.90;
constexpr double scer_exit_thick = 0.050;
// (aerogel and shower-counter radiation lengths, l.92-109, are declared but not used)

// Wire chamber resolutions (sigma): real*8 sdc_sigma(1:12) / 12*0.030 /
constexpr double sdc_sigma[12] = {0.030, 0.030, 0.030, 0.030, 0.030, 0.030,
                                  0.030, 0.030, 0.030, 0.030, 0.030, 0.030};

// Wire plane positions (l.116-143)
constexpr int sdc_nr_cham = 2;
constexpr int sdc_nr_plan = 6;
constexpr double sdc_1_zpos = 6.25;
constexpr double sdc_2_zpos = 55.77;
constexpr double sdc_del_plane = sdc_thick + sdc_wire_thick + sdc_cath_thick;
constexpr double sdc_1_left = 24.0;
constexpr double sdc_1_right = -24.0;
constexpr double sdc_1y_offset = -1.822;
constexpr double sdc_1_top = -32.0;
constexpr double sdc_1_bot = 32.0;
constexpr double sdc_1x_offset = 8.649;
constexpr double sdc_2_left = 24.0;
constexpr double sdc_2_right = -24.0;
constexpr double sdc_2y_offset = -1.976;
constexpr double sdc_2_top = -32.0;
constexpr double sdc_2_bot = 32.0;
constexpr double sdc_2x_offset = -1.532;

// Scintillator positions and thicknesses (l.146-173)
constexpr double sscin_1y_zpos = 73.61;
constexpr double sscin_1x_zpos = 97.11;
constexpr double sscin_2y_zpos = 249.51;
constexpr double sscin_2x_zpos = 290.81;
constexpr double sscin_1x_thick = 1.040;   // 1.040 for overlap
constexpr double sscin_1y_thick = 1.098;
constexpr double sscin_2x_thick = 1.040;
constexpr double sscin_2y_thick = 1.098;
constexpr double sscin_1_left = 18.25;
constexpr double sscin_1_right = -18.25;
constexpr double sscin_1x_offset = 2.8;
constexpr double sscin_1_top = -31.75;
constexpr double sscin_1_bot = 31.75;
constexpr double sscin_1y_offset = 2.25;
constexpr double sscin_2_left = 18.25;
constexpr double sscin_2_right = -18.25;
constexpr double sscin_2x_offset = 4.9;
constexpr double sscin_2_top = -56.25;
constexpr double sscin_2_bot = 56.25;
constexpr double sscin_2y_offset = 2.9;

// Cherenkov position
constexpr double scer_zentrance = 130.000;
constexpr double scer_zmirror = 155.000;
constexpr double scer_zexit = 160.000;

// Calorimeter position (only the back plane is used: the drift to the back of the calorimeter; the aperture
// limits scal_left/right/top/bottom, l.188-191, belong to the commented-out cuts)
constexpr double scal_4ta_zpos = 346.01;

constexpr int scintrig = 3;           // set trigger to 3/4
constexpr double nsig_max = 99.0e0;   // max #/sigma for gaussian ran #s

bool inHodo1(const optics::Track& t) {   // mc_sos_hut.f:422-425, 433-436
    return t.ys < (sscin_1_left + sscin_1y_offset) && t.ys > (sscin_1_right + sscin_1y_offset) &&
           t.xs < (sscin_1_bot + sscin_1x_offset) && t.xs > (sscin_1_top + sscin_1x_offset);
}
bool inHodo2(const optics::Track& t) {   // mc_sos_hut.f:472-475, 483-486
    return t.ys < (sscin_2_left + sscin_2y_offset) && t.ys > (sscin_2_right + sscin_2y_offset) &&
           t.xs < (sscin_2_bot + sscin_2x_offset) && t.xs > (sscin_2_top + sscin_2x_offset);
}

}  // namespace

int McSos::mcSosHut(optics::DecayState& d, McSosArgs& a, double zinit, IRandom& rng, McSosDetail& det) {
    optics::Track& t = track_;
    double& pathlen = a.pathlen;
    // m2 = d.m2, p = d.ph (updated by project() on decay)
    auto musc = [&](double radw) {
        if (a.ms_flag) mat_.musc(d.m2, d.ph, radw, t.dydzs, t.dxdzs, rng);
    };
    auto muscExt = [&](double radw, double drift) {
        if (a.ms_flag) mat_.muscExt(d.m2, d.ph, radw, drift, t.dydzs, t.dxdzs, t.ys, t.xs, rng);
    };
    auto project = [&](double drift) { optics::project(t, drift, pathlen, &d); };

    // Resolution multiplier: 'more conservative values' (mc_sos_hut.f:251-257). Always drawn.
    const double tmpran = rng.uniform();
    const double tmplim = 0.15;   // fraction of events with >6 hits
    double resmult;
    if (tmpran < tmplim) {
        resmult = 2.0;   // resolution is 3.0x worse (both DCs) [sic]
    } else {
        resmult = 1.0;
    }
    a.resmult = resmult;

    // Initialize the xdc and ydc arrays to zero (real*4 for lfit)
    float xdc[12], ydc[12], zdc[12];
    for (int i = 0; i < 12; ++i) {
        xdc[i] = 0.f;
        ydc[i] = 0.f;
    }
    int scincount = 0;
    double radw, drift;

    // Go to spectrometer exit foil. (Drift forwards from zinit.) (mc_sos_hut.f:282-287)
    // Approximation is used in calculating xt = x at foil crossing.
    // FORTRAN-QUIRK(D-061 F2): zinit = -3.206267 lies after the foil at -3.22, so drift < 0 and is clamped to
    // 0.001 cm (l.284).
    const double xt = t.xs + t.dxdzs * sfoil_exit_zpos;
    drift = (sfoil_exit_zpos + xt * std::tan(sfoil_exit_ang)) - zinit;
    if (drift <= 0.001) drift = 0.001;   // avoid drift<=0
    project(drift);                      // drift and decay
    radw = sfoil_exit_thick / sfoil_exit_radlen / std::cos(sfoil_exit_ang);
    musc(radw);

    // Go to first drift chamber set (air MS before the DCs instead of 1/2 way through), l.297-301
    drift = (sdc_1_zpos - 0.5 * sdc_nr_plan * sdc_del_plane) - (sfoil_exit_zpos + xt * std::tan(sfoil_exit_ang));
    radw = drift / sair_radlen;
    project(drift);
    muscExt(radw, drift);

    // One chamber: mc_sos_hut.f:303-334 (jchamber=1) and 352-384 (jchamber=2), identical code.
    auto chamber = [&](int jchamber) {
        const int npl_off = (jchamber - 1) * sdc_nr_plan;
        for (int iplane = 1; iplane <= sdc_nr_plan; ++iplane) {
            radw = sdc_cath_thick / sdc_cath_radlen;
            musc(radw);
            drift = 0.5 * sdc_thick;
            radw = drift / sdc_radlen;
            musc(radw);
            drift = drift + sdc_cath_thick;
            project(drift);
            radw = sdc_wire_thick / sdc_wire_radlen;
            musc(radw);
            double tmpran1, tmpran2;
            if (a.wcs_flag) {
                tmpran1 = rng.gauss1(nsig_max);   // gaussian, truncated at 99 sigma
                tmpran2 = rng.gauss1(nsig_max);
            } else {
                tmpran1 = 0.;
                tmpran2 = 0.;
            }
            const int k = npl_off + iplane - 1;
            xdc[k] = static_cast<float>(t.xs + sdc_sigma[k] * tmpran1 * resmult);
            ydc[k] = static_cast<float>(t.ys + sdc_sigma[k] * tmpran2 * resmult);
            if (iplane == 1 || iplane == 3 || iplane == 5) {
                xdc[k] = 0.f;   // assume 3 planes are x, 3 are y
            } else {
                ydc[k] = 0.f;
            }
            drift = 0.5 * sdc_thick;
            radw = drift / sdc_radlen;
            musc(radw);
            drift = drift + sdc_wire_thick;
            project(drift);
        }
    };

    chamber(1);
    if (t.xs > (sdc_1_bot - sdc_1x_offset) || t.xs < (sdc_1_top - sdc_1x_offset) ||
        t.ys > (sdc_1_left - sdc_1y_offset) || t.ys < (sdc_1_right - sdc_1y_offset)) {
        return stop(kDc1);
    }
    radw = sdc_cath_thick / sdc_cath_radlen;
    musc(radw);

    // at last cathode foil of first drift chamber set, drift to next (l.347-350)
    drift = sdc_2_zpos - sdc_1_zpos - sdc_nr_plan * sdc_del_plane;
    radw = drift / sair_radlen;
    project(drift);
    muscExt(radw, drift);

    chamber(2);
    if (t.xs > (sdc_2_bot - sdc_2x_offset) || t.xs < (sdc_2_top - sdc_2x_offset) ||
        t.ys > (sdc_2_left - sdc_2y_offset) || t.ys < (sdc_2_right - sdc_2y_offset)) {
        return stop(kDc2);
    }
    radw = sdc_cath_thick / sdc_cath_radlen;
    musc(radw);

    // fit track to give focal plane values, use LFIT from GENLIB (mc_sos_hut.f:399-414)
    for (int jchamber = 1; jchamber <= sdc_nr_cham; ++jchamber) {
        const int npl_off = (jchamber - 1) * sdc_nr_plan;
        for (int iplane = 1; iplane <= sdc_nr_plan; ++iplane) {
            if (jchamber == 1)
                zdc[npl_off + iplane - 1] =
                    static_cast<float>(sdc_1_zpos + (iplane - 0.5 - 0.5 * sdc_nr_plan) * sdc_del_plane);
            if (jchamber == 2)
                zdc[npl_off + iplane - 1] =
                    static_cast<float>(sdc_2_zpos + (iplane - 0.5 - 0.5 * sdc_nr_plan) * sdc_del_plane);
        }
    }
    float badf, xfp4, yfp4, dxfp4, dyfp4;
    cern::lfit(zdc, xdc, 12, 0, dxfp4, xfp4, badf);
    cern::lfit(zdc, ydc, 12, 0, dyfp4, yfp4, badf);
    a.x_fp = static_cast<double>(xfp4);
    a.y_fp = static_cast<double>(yfp4);
    a.dx_fp = static_cast<double>(dxfp4);
    a.dy_fp = static_cast<double>(dyfp4);
    det.fit_done = true;
    det.xfp4 = xfp4;
    det.dxfp4 = dxfp4;
    det.yfp4 = yfp4;
    det.dyfp4 = dyfp4;

    // at last cathode foil of second drift chamber set, drift to hodoscopes (l.418-427)
    drift = sscin_1y_zpos - sdc_2_zpos - 0.5 * sdc_nr_plan * sdc_del_plane;
    radw = drift / sair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (inHodo1(t)) scincount = scincount + 1;
    radw = sscin_1y_thick / sscin_radlen;
    musc(radw);

    drift = sscin_1x_zpos - sscin_1y_zpos;
    radw = drift / sair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (inHodo1(t)) scincount = scincount + 1;
    radw = sscin_1x_thick / sscin_radlen;
    musc(radw);

    // finished first hodoscope, drift to cherenkov (l.442-464)
    drift = scer_zentrance - sscin_1x_zpos;
    radw = drift / sair_radlen;
    project(drift);
    muscExt(radw, drift);

    radw = scer_entr_thick / scer_entr_radlen;
    musc(radw);

    drift = scer_zmirror - scer_zentrance;
    radw = drift / scer_radlen;
    project(drift);
    muscExt(radw, drift);

    radw = scer_mir_thick / scer_mir_radlen;
    musc(radw);

    drift = scer_zexit - scer_zmirror;
    radw = drift / scer_radlen;
    project(drift);
    muscExt(radw, drift);

    radw = scer_exit_thick / scer_exit_radlen;
    musc(radw);

    // drift to second hodoscope (l.468-488)
    drift = sscin_2y_zpos - scer_zexit;
    radw = drift / sair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (inHodo2(t)) scincount = scincount + 1;
    radw = sscin_2y_thick / sscin_radlen;
    musc(radw);

    drift = sscin_2x_zpos - sscin_2y_zpos;
    radw = drift / sair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (inHodo2(t)) scincount = scincount + 1;
    radw = sscin_2x_thick / sscin_radlen;
    musc(radw);

    // Test on scintillator trigger (DJG 8/18/98)
    if (scincount < scintrig) return stop(kScin);

    // Drift to back of calorimeter (so that pathlen is calculated to the back of the spectrometer), l.500-503.
    drift = scal_4ta_zpos - sscin_2x_zpos;
    radw = drift / sair_radlen;
    project(drift);
    muscExt(radw, drift);

    // FORTRAN-QUIRK(D-061 F7): the calorimeter aperture cut (l.519-525) and the fiducial cut on the fitted track
    // (l.535-539) are commented out, so sSTOP_cal is never incremented; xcal/ycal (l.533-534) are computed but unused.
    return kSuccess;   // ok_hut = .true.
}

}  // namespace sos
}  // namespace simc
