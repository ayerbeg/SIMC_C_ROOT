// src/spectrometers/hrs/McHrsHut.cpp
// Fortran: hrsl/mc_hrsl_hut.f (l.1-498) = hrsr/mc_hrsr_hut.f (names only differ) — the HRS detector hut: exit foil, two
// VDCs (6 planes each, wire-chamber smearing, CERNLIB lfit), S1, Cherenkov, S2, drift to the back of the calorimeter;
// multiple scattering in every layer if ms_flag. PARAMETERs copied verbatim (real*8 under -fdefault-real-8).
#include <cmath>

#include "simc/optics/Lfit.h"
#include "simc/spectrometers/hrs/kernel/McHrs.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace hrs {

namespace {

// spectrometer exit window, 0.1 mm of titanium (X0=3.56cm)
constexpr double hfoil_exit_radlen = 3.56;
constexpr double hfoil_exit_thick = 0.01;
// spectrometer air gaps
constexpr double hair_radlen = 30420.;
// chamber entrance foil, .18 mm of Kapton
constexpr double hdc_entr_radlen = 34.4;
constexpr double hdc_entr_thick = 0.00018 * 2.54;
// chamber gas, 35/65 ethane/argon
constexpr double hdc_radlen = 16700.0;
constexpr double hdc_thick = 1.5;
// effective wire plane
constexpr double hdc_wire_radlen = 0.35;
constexpr double hdc_wire_thick = 0.0000049;
// effective cathode plane, Be/Cu
constexpr double hdc_cath_radlen = 7.2;
constexpr double hdc_cath_thick = 0.000177;
// chamber exit foil
constexpr double hdc_exit_radlen = 34.4;
constexpr double hdc_exit_thick = 0.00018 * 2.54;
// hodoscopes
constexpr double hscin_radlen = 42.4;
// Cherenkov entrance foil, 0.5 atm CO2, mirror, exit foil
constexpr double hcer_entr_radlen = 8.90;
constexpr double hcer_entr_thick = 0.040 * 2.54;
constexpr double hcer_radlen = 36620.0;
constexpr double hcer_mir_radlen = 400.0;
constexpr double hcer_mir_thick = 2.0;
constexpr double hcer_exit_radlen = 8.90;
constexpr double hcer_exit_thick = 0.040 * 2.54;

constexpr double hdc_sigma[12] = {0.0225, 0.0225, 0.0225, 0.0225, 0.0225, 0.0225,
                                  0.0225, 0.0225, 0.0225, 0.0225, 0.0225, 0.0225};
constexpr int hdc_nr_cham = 2;
constexpr int hdc_nr_plan = 6;
// Drift chamber 1 is the focal plane, so shift all zpos values by 25cm
constexpr double hdc_1_zpos = -25.0 + 25.0;
constexpr double hdc_2_zpos = 25.0 + 25.0;
constexpr double hdc_del_plane = hdc_thick + hdc_wire_thick + hdc_cath_thick;
constexpr double hdc_1_left = 14.4;
constexpr double hdc_1_right = -14.4;
constexpr double hdc_1y_offset = 0.000;
constexpr double hdc_1_top = -105.6;
constexpr double hdc_1_bot = 105.6;
constexpr double hdc_1x_offset = 0.000;
constexpr double hdc_2_left = 14.4;
constexpr double hdc_2_right = -14.4;
constexpr double hdc_2y_offset = 0.000;
constexpr double hdc_2_top = -105.6;
constexpr double hdc_2_bot = 105.6;
constexpr double hdc_2x_offset = 0.000;
// Scintillators
constexpr double hscin_1x_zpos = 95.0 + 25.0;
constexpr double hscin_2x_zpos = 288.3 + 25.0;
constexpr double hscin_1x_thick = 0.5 * 1.067;
constexpr double hscin_2x_thick = 0.5 * 1.067;
constexpr double hscin_1x_left = 18.0;
constexpr double hscin_1x_right = -18.0;
constexpr double hscin_2x_left = 30.0;
constexpr double hscin_2x_right = -30.0;
// Cherenkov
constexpr double hcer_zentrance = 137.0 + 25.0;
constexpr double hcer_zmirror = 197.0 + 25.0;
constexpr double hcer_zexit = 237.0 + 25.0;
// Calorimeter back plane (the aperture/fiducial cuts are commented out, l.463-470, 480-487)
constexpr double hcal_4ta_zpos = 407.3 + 25.0;

constexpr double nsig_max = 99.0;   // parameter(nsig_max=99.0d0): REAL(16) literal, real*8 PARAMETER -> 99 exactly

}  // namespace

int McHrs::mcHrsHut(optics::DecayState& d, McHrsArgs& a, double zinit, IRandom& rng, McHrsDetail& det) {
    optics::Track& t = track_;
    double& pathlen = a.pathlen;
    auto musc = [&](double radw) {
        if (a.ms_flag) mat_.musc(d.m2, d.ph, radw, t.dydzs, t.dxdzs, rng);
    };
    auto muscExt = [&](double radw, double drift) {
        if (a.ms_flag) mat_.muscExt(d.m2, d.ph, radw, drift, t.dydzs, t.dxdzs, t.ys, t.xs, rng);
    };
    auto project = [&](double drift) { optics::project(t, drift, pathlen, &d); };
    auto stopHere = [&](int code, double where) { return stop(code, where, t.xs, t.ys); };
    // rotate_haxis(45.0d0,xt,yt) (l.291-293, 354-356). FORTRAN-BUG(P14.2 H2): 45.0d0 is REAL(16); the real*8 dummy rotang
    // reads its low 8 bytes, +0.0, so the Fortran applies NO rotation (rotate_haxis with 0 degrees returns xt, yt
    // unchanged). Corrected: 45 degrees as written. fortran_compat: 0 degrees.
    const double vdcRot = fortranCompat() ? low8(static_cast<Quad>(45)) : 45.0;

    // Initialize the xdc and ydc arrays to zero (real*4 for lfit)
    float xdc[12], ydc[12], zdc[12];
    for (int i = 0; i < 12; ++i) {
        xdc[i] = 0.f;
        ydc[i] = 0.f;
    }
    const double resmult = 1.0;   // l.243: the HRS has no resolution multiplier
    a.resmult = resmult;
    double radw, drift, xt, yt;

    // Scatter in spectrometer exit foil, which is located at zinit (l.252-253).
    radw = hfoil_exit_thick / hfoil_exit_radlen;
    musc(radw);

    // Go to first drift chamber set (l.259-262)
    drift = (hdc_1_zpos - 0.5 * hdc_nr_plan * hdc_del_plane) - zinit;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);

    // One chamber: l.264-302 (jchamber=1) and 316-353 (jchamber=2), identical code.
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
                tmpran1 = rng.gauss1(nsig_max);
                tmpran2 = rng.gauss1(nsig_max);
            } else {
                tmpran1 = 0.;
                tmpran2 = 0.;
            }
            const int k = npl_off + iplane - 1;
            xdc[k] = static_cast<float>(t.xs + hdc_sigma[k] * tmpran1 * resmult);
            ydc[k] = static_cast<float>(t.ys + hdc_sigma[k] * tmpran2 * resmult);
            if (iplane == 1 || iplane == 3 || iplane == 5) {
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
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, vdcRot, xt, yt, fortranCompat());
    if (xt > (hdc_1_bot - hdc_1x_offset) || xt < (hdc_1_top - hdc_1x_offset) || yt > (hdc_1_left - hdc_1y_offset) ||
        yt < (hdc_1_right - hdc_1y_offset))
        return stopHere(kDc1, 17.);
    radw = hdc_cath_thick / hdc_cath_radlen;
    musc(radw);

    // at last cathode foil of first drift chamber set, drift to next (l.310-314)
    drift = hdc_2_zpos - hdc_1_zpos - hdc_nr_plan * hdc_del_plane;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);

    chamber(2);
    xt = t.xs;
    yt = t.ys;
    optics::rotateHaxis(t, vdcRot, xt, yt, fortranCompat());
    if (xt > (hdc_2_bot - hdc_2x_offset) || xt < (hdc_2_top - hdc_2x_offset) || yt > (hdc_2_left - hdc_2y_offset) ||
        yt < (hdc_2_right - hdc_2y_offset))
        return stopHere(kDc2, 19.);
    radw = hdc_cath_thick / hdc_cath_radlen;
    musc(radw);

    // fit track to give new focal plane values, use LFIT from GENLIB (l.372-389)
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

    // drift to S1 (l.393-405)
    drift = hscin_1x_zpos - hdc_2_zpos - 0.5 * hdc_nr_plan * hdc_del_plane;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (t.ys > hscin_1x_left || t.ys < hscin_1x_right) return stopHere(kS1, 20.);
    radw = hscin_1x_thick / hscin_radlen;
    musc(radw);

    // finished first hodoscope, drift to cherenkov (l.409-430)
    drift = hcer_zentrance - hscin_1x_zpos;
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
    drift = hcer_zexit - hcer_zmirror;
    radw = drift / hcer_radlen;
    project(drift);
    muscExt(radw, drift);
    radw = hcer_exit_thick / hcer_exit_radlen;
    musc(radw);

    // drift to second hodoscope (l.434-448)
    drift = hscin_2x_zpos - hcer_zexit;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (t.ys > hscin_2x_left || t.ys < hscin_2x_right) return stopHere(kS2, 21.);
    radw = hscin_2x_thick / hscin_radlen;
    musc(radw);

    // Drift to back of calorimeter (l.452-455)
    drift = hcal_4ta_zpos - hscin_2x_zpos;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);

    // FORTRAN-QUIRK(P14.2): the calorimeter aperture and fiducial cuts are commented out (l.463-470, 480-487): the cal
    // counter is never incremented; xcal/ycal (l.477-478) are unused.
    return kSuccess;   // ok_hut = .true.
}

}  // namespace hrs
}  // namespace simc
