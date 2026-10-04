// src/spectrometers/hms/McHmsHut.cpp
// Fortran: hms/mc_hms_hut.f:mc_hms_hut (l.1-608) — the HMS detector hut: exit window, two drift chambers
// (6 planes each, wire-chamber smearing, CERNLIB lfit of the focal-plane track), aerogel, hodoscopes (3/4
// trigger), Cherenkov, calorimeter; multiple scattering in every layer if ms_flag.
// All PARAMETERs are copied verbatim; with -fdefault-real-8 they are real*8, and the parameter-only
// expressions are folded in double, like the constexpr ones here.
#include <cmath>

#include "simc/spectrometers/hms/kernel/McHms.h"
#include "simc/optics/Lfit.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace hms {

namespace {

// spectrometer exit window: 0.011 inch preformed aluminum (mc_hms_hut.f:57-59)
constexpr double hfoil_exit_radlen = 8.90;
constexpr double hfoil_exit_thick = 0.011 * 2.54;
// spectrometer air gaps
constexpr double hair_radlen = 30420.;
// chamber entrance foil, 1 mil of Mylar
constexpr double hdc_entr_radlen = 28.7;
constexpr double hdc_entr_thick = 0.001 * 2.54;
// chamber gas, 50/50 ethane/argon
constexpr double hdc_radlen = 16700.0;
constexpr double hdc_thick = 1.8;
// effective wire plane
constexpr double hdc_wire_radlen = 0.35;
constexpr double hdc_wire_thick = 0.0000049;
// effective cathode plane
constexpr double hdc_cath_radlen = 7.2;
constexpr double hdc_cath_thick = 0.000177;
// chamber exit foil
constexpr double hdc_exit_radlen = 28.7;
constexpr double hdc_exit_thick = 0.001 * 2.54;
// Aerogel entrance foil / aerogel / air / exit foil
constexpr double haer_entr_radlen = 8.90;
constexpr double haer_entr_thick = 0.15;
constexpr double haer_radlen = 150.0;
constexpr double haer_thick = 9.0;
constexpr double haer_air_radlen = 30420.0;
constexpr double haer_air_thick = 16.0;
constexpr double haer_exit_radlen = 8.90;
constexpr double haer_exit_thick = 0.1;
// hodoscopes
constexpr double hscin_radlen = 42.4;
// Cherenkov entrance foil, gas (0.5 atm Freon), mirror, exit foil
constexpr double hcer_entr_radlen = 8.90;
constexpr double hcer_entr_thick = 0.040 * 2.54;
constexpr double hcer_radlen = 9620.0;
constexpr double hcer_mir_radlen = 400.0;
constexpr double hcer_mir_thick = 2.0;
constexpr double hcer_exit_radlen = 8.90;
constexpr double hcer_exit_thick = 0.040 * 2.54;
// shower counter (only used by the commented-out cut)
// Wire chamber resolutions (sigma): real*8 hdc_sigma(1:12) / 12*0.030 /
constexpr double hdc_sigma[12] = {0.030, 0.030, 0.030, 0.030, 0.030, 0.030,
                                  0.030, 0.030, 0.030, 0.030, 0.030, 0.030};
// Wire plane positions
constexpr int hdc_nr_cham = 2;
constexpr int hdc_nr_plan = 6;
constexpr double hdc_1_zpos = -52.1084;
constexpr double hdc_2_zpos = 29.2608;
constexpr double hdc_del_plane = hdc_thick + hdc_wire_thick + hdc_cath_thick;
constexpr double hdc_1_left = 26.0;
constexpr double hdc_1_right = -26.0;
constexpr double hdc_1y_offset = 1.443;
constexpr double hdc_1_top = -56.5;
constexpr double hdc_1_bot = 56.5;
constexpr double hdc_1x_offset = 1.670;
constexpr double hdc_2_left = 26.0;
constexpr double hdc_2_right = -26.0;
constexpr double hdc_2y_offset = 2.753;
constexpr double hdc_2_top = -56.5;
constexpr double hdc_2_bot = 56.5;
constexpr double hdc_2x_offset = 2.758;
// Aerogel positions
constexpr double haer_zentrance = 35.699;
constexpr double haer_zexit = 60.949;
// Scintillator positions and thicknesses
constexpr double hscin_1x_zpos = 77.830;
constexpr double hscin_1y_zpos = 97.520;
constexpr double hscin_2x_zpos = 298.820;
constexpr double hscin_2y_zpos = 318.510;
constexpr double hscin_1x_thick = 1.067;
constexpr double hscin_1y_thick = 1.067;
constexpr double hscin_2x_thick = 1.067;
constexpr double hscin_2y_thick = 1.067;
constexpr double hscin_1x_left = 37.75;
constexpr double hscin_1x_right = -37.75;
constexpr double hscin_1x_offset = -1.3;
constexpr double hscin_1y_top = -60.25;
constexpr double hscin_1y_bot = 60.25;
constexpr double hscin_1y_offset = -1.3;
// hscin_2* limits/offsets exist in the Fortran (2x_offset -0.60, 2y_offset -2.40) but are never used:
// FORTRAN-DIFF-CANDIDATE(D-015 F7): mc_hms_hut.f:537-551 test the second hodoscope with the hscin_1* values.
// Question: copy-paste, or the same active area on purpose? (offsets differ: -0.60/-2.40 vs hodoscope 1)
// Cherenkov position
constexpr double hcer_zentrance = 110.000;
constexpr double hcer_zmirror = 230.000;
constexpr double hcer_zexit = 265.000;
// Calorimeter position
constexpr double hcal_4ta_zpos = 371.69;

constexpr int scintrig = 3;          // set trigger to 3/4
constexpr double nsig_max = 99.0e0;  // max #/sigma for gaussian ran #s

bool inHodo1(const optics::Track& t) {
    return t.ys < (hscin_1x_left + hscin_1y_offset) && t.ys > (hscin_1x_right + hscin_1y_offset) &&
           t.xs < (hscin_1y_bot + hscin_1x_offset) && t.xs > (hscin_1y_top + hscin_1x_offset);
}

}  // namespace

int McHms::mcHmsHut(optics::DecayState& d, McHmsArgs& a, double zinit, IRandom& rng, McHmsDetail& det) {
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

    // Resolution multiplier: 'more conservative values' (mc_hms_hut.f:292-298). Always drawn.
    const double tmpran = rng.uniform();
    const double tmplim = 0.15;   // fraction of events with >6 hits
    double resmult;
    if (tmpran < tmplim) {
        resmult = 2.0;
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

    // Go to spectrometer exit foil, 25cm before DC1 (Drift forwards from zinit).
    drift = (hdc_1_zpos - 25.000) - zinit;
    project(drift);
    radw = hfoil_exit_thick / hfoil_exit_radlen;
    musc(radw);

    // Go to first drift chamber set (air mult. scattering before drifting... applied after project)
    drift = (hdc_1_zpos - 0.5 * hdc_nr_plan * hdc_del_plane) - (hdc_1_zpos - 25.000);
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);

    // One chamber: mc_hms_hut.f:337-370 (jchamber=1) and 391-425 (jchamber=2), identical code.
    auto chamber = [&](int jchamber) {
        radw = hdc_entr_thick / hdc_entr_radlen;
        musc(radw);
        const int npl_off = (jchamber - 1) * hdc_nr_plan;
        for (int iplane = 1; iplane <= hdc_nr_plan; ++iplane) {
            radw = hdc_cath_thick / hdc_cath_radlen;
            musc(radw);
            drift = 0.5 * hdc_thick;
            radw = drift / hdc_radlen;
            musc(radw);
            drift = drift + hdc_cath_thick;
            project(drift);
            radw = hdc_wire_thick / hdc_wire_radlen;
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
            xdc[k] = static_cast<float>(t.xs + hdc_sigma[k] * tmpran1 * resmult);
            ydc[k] = static_cast<float>(t.ys + hdc_sigma[k] * tmpran2 * resmult);
            if (iplane == 2 || iplane == 5) {
                xdc[k] = 0.f;   // y plane, no x information
            } else {
                ydc[k] = 0.f;   // x-like plane, no y info
            }
            drift = 0.5 * hdc_thick;
            radw = drift / hdc_radlen;
            musc(radw);
            drift = drift + hdc_wire_thick;
            project(drift);
        }
        radw = hdc_exit_thick / hdc_exit_radlen;
        musc(radw);
    };

    chamber(1);
    if (t.xs > (hdc_1_bot - hdc_1x_offset) || t.xs < (hdc_1_top - hdc_1x_offset) ||
        t.ys > (hdc_1_left - hdc_1y_offset) || t.ys < (hdc_1_right - hdc_1y_offset)) {
        return stop(kDc1);
    }
    radw = hdc_cath_thick / hdc_cath_radlen;
    musc(radw);

    // at last cathode foil of first drift chamber set, drift to next
    drift = (hdc_2_zpos - 0.5 * hdc_nr_plan * hdc_del_plane) - (hdc_1_zpos + 0.5 * hdc_nr_plan * hdc_del_plane);
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);

    chamber(2);
    if (t.xs > (hdc_2_bot - hdc_2x_offset) || t.xs < (hdc_2_top - hdc_2x_offset) ||
        t.ys > (hdc_2_left - hdc_2y_offset) || t.ys < (hdc_2_right - hdc_2y_offset)) {
        return stop(kDc2);
    }
    radw = hdc_cath_thick / hdc_cath_radlen;
    musc(radw);

    // fit track to give new focal plane values, use LFIT from GENLIB (mc_hms_hut.f:439-455)
    for (int jchamber = 1; jchamber <= hdc_nr_cham; ++jchamber) {
        const int npl_off = (jchamber - 1) * hdc_nr_plan;
        for (int iplane = 1; iplane <= hdc_nr_plan; ++iplane) {
            if (jchamber == 1)
                zdc[npl_off + iplane - 1] =
                    static_cast<float>(hdc_1_zpos + (iplane - 0.5 - 0.5 * hdc_nr_plan) * hdc_del_plane);
            if (jchamber == 2)
                zdc[npl_off + iplane - 1] =
                    static_cast<float>(hdc_2_zpos + (iplane - 0.5 - 0.5 * hdc_nr_plan) * hdc_del_plane);
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

    // at last cathode foil of second drift chamber set, drift to aerogel
    drift = haer_zentrance - hdc_2_zpos - 0.5 * hdc_nr_plan * hdc_del_plane;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);

    radw = haer_entr_thick / haer_entr_radlen;
    musc(radw);

    drift = haer_thick;
    radw = drift / haer_radlen;
    project(drift);
    muscExt(radw, drift);

    drift = haer_air_thick;
    radw = drift / haer_air_radlen;
    project(drift);
    muscExt(radw, drift);

    radw = haer_exit_thick / haer_exit_radlen;
    musc(radw);

    // at aerogel exit, drift to hodoscopes
    drift = hscin_1x_zpos - haer_zexit;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (inHodo1(t)) scincount = scincount + 1;
    radw = hscin_1x_thick / hscin_radlen;
    musc(radw);

    drift = hscin_1y_zpos - hscin_1x_zpos;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (inHodo1(t)) scincount = scincount + 1;
    radw = hscin_1y_thick / hscin_radlen;
    musc(radw);

    // finished first hodoscope, drift to cherenkov
    drift = hcer_zentrance - hscin_1y_zpos;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);

    radw = hcer_entr_thick / hcer_entr_radlen;
    musc(radw);

    drift = hcer_zmirror - hcer_zentrance;
    radw = drift / hcer_radlen;
    project(drift);
    muscExt(radw, drift);

    radw = hcer_mir_thick / hcer_mir_radlen;
    musc(radw);

    // FORTRAN-DIFF-CANDIDATE(D-015 F7): mc_hms_hut.f:523-529 apply the exit-foil scattering before the drift from the
    // mirror to the exit AND again after it. Question: double counting, or a second foil?
    drift = hcer_zexit - hcer_zmirror;
    radw = hcer_exit_thick / hcer_exit_radlen;
    musc(radw);
    project(drift);

    radw = hcer_exit_thick / hcer_exit_radlen;
    musc(radw);

    // drift to second hodoscope
    drift = hscin_2x_zpos - hcer_zexit;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (inHodo1(t)) scincount = scincount + 1;   // FORTRAN-DIFF-CANDIDATE(D-015 F7): hodoscope-1 limits (l.537-540)
    radw = hscin_2x_thick / hscin_radlen;
    musc(radw);

    drift = hscin_2y_zpos - hscin_2x_zpos;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (inHodo1(t)) scincount = scincount + 1;   // FORTRAN-DIFF-CANDIDATE(D-015 F7): hodoscope-1 limits (l.548-551)
    radw = hscin_2y_thick / hscin_radlen;
    musc(radw);

    // Test on scintillator trigger (DJG 8/18/98)
    if (scincount < scintrig) return stop(kScin);

    // Drift to back of calorimeter (so that pathlen is calculated to the back of the spectrometer).
    drift = hcal_4ta_zpos - hscin_2y_zpos;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);

    // FORTRAN-QUIRK (mc_hms_hut.f:582-588): the calorimeter aperture cut is commented out, so hSTOP_cal is
    // never incremented. xcal/ycal (l.596-597) are computed but unused (the fiducial cut is commented out too).
    return kSuccess;   // ok_hut = .true.
}

}  // namespace hms
}  // namespace simc
