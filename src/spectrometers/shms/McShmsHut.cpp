// src/spectrometers/shms/McShmsHut.cpp
// Fortran: shms/mc_shms_hut.f:mc_shms_hut (l.1-458) — the SHMS detector hut: heavy-gas Cherenkov (or vacuum pipe /
// helium bag, unreachable: mc_shms sets cer_flag = .true. on every call), two drift chambers (lfit), S1 (x,y),
// noble-gas Cherenkov, S2 (x,y), calorimeter; use_det_cut = .true. (every detector cut stops the ray).
// The final `close (unit=chan)` (l.455) has no counterpart. The real*4 `spec` array is not ported (never output).
#include <cmath>

#include "ShmsHutConstants.h"
#include "simc/optics/Lfit.h"
#include "simc/spectrometers/shms/ShmsSpectrometer.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace shms {

using namespace hut;

namespace {
constexpr double nsig_max = 99.0e0;
constexpr bool use_det_cut = true;   // mc_shms_hut.f:66
}  // namespace

int ShmsSpectrometer::mcShmsHut(optics::DecayState& d, McShmsArgs& a, double /*zinit (unused)*/, IRandom& rng,
                                McShmsDetail& det) {
    optics::Track& t = track_;
    double& pathlen = a.pathlen;
    auto musc = [&](double radw) {
        if (a.ms_flag) mat_.musc(d.m2, d.ph, radw, t.dydzs, t.dxdzs, rng);
    };
    auto muscExt = [&](double radw, double drift) {
        if (a.ms_flag) mat_.muscExt(d.m2, d.ph, radw, drift, t.dydzs, t.dxdzs, t.ys, t.xs, rng);
    };
    auto project = [&](double drift) { optics::project(t, drift, pathlen, &d); };
    // stop with stop_id; with use_det_cut = .false. the Fortran would only count and continue
    auto hutStop = [&](int code, int id) {
        cnt_[code] = cnt_[code] + 1;
        cnt_[kStopIdSlot] = id;
        return code;
    };

    // Initialize (mc_shms_hut.f:74-88); hdc_del_plane is a local variable here (not a PARAMETER)
    const double hdc_del_plane = hdc_thick + hdc_wire_thick + hdc_cath_thick;
    const double resmult = 1.0;
    a.resmult = resmult;
    float xdc[12], ydc[12], zdc[12];
    for (int i = 0; i < 12; ++i) {
        xdc[i] = 0.f;
        ydc[i] = 0.f;
    }
    double radw, drift;

    if (cerFlag_) {
        // (the drift to the Cherenkov entrance is done in mc_shms)
        radw = hfoil_exit_thick / hfoil_exit_radlen;
        musc(radw);
        radw = hcer_entr_thick / hcer_entr_radlen;
        musc(radw);
        drift = hcer_1_zmirror - hcer_1_zentrance - hcer_mirglass_thick / 2;
        radw = drift / hcer_1_radlen;
        project(drift);
        muscExt(radw, drift);

        radw = hcer_mirglass_thick / hcer_mirglass_radlen;
        musc(radw);

        drift = hcer_1_zexit - hcer_1_zmirror - hcer_mirglass_thick / 2;
        radw = drift / hcer_1_radlen;
        project(drift);
        muscExt(radw, drift);
        radw = hcer_exit_thick / hcer_exit_radlen;
        musc(radw);
    } else {
        if (vacFlag_) {   // vacuum pipe
            radw = hfoil_exit_thick / hfoil_exit_radlen;
            musc(radw);
        } else {   // helium bag
            radw = hfoil_exit_thick / hfoil_exit_radlen;
            musc(radw);
            radw = helbag_al_thick / helbag_al_radlen;
            musc(radw);
            radw = helbag_mylar_thick / helbag_mylar_radlen;
            musc(radw);
            drift = hcer_1_zexit - hcer_1_zentrance;
            radw = drift / helbag_hel_radlen;
            project(drift);
            muscExt(radw, drift);
            radw = helbag_al_thick / helbag_al_radlen;
            musc(radw);
            radw = helbag_mylar_thick / helbag_mylar_radlen;
            musc(radw);
        }
    }

    // Go to first drift chamber set
    drift = (hdc_1_zpos - 0.5 * hdc_nr_plan * hdc_del_plane) - hcer_1_zexit;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);

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
        if (use_det_cut) return hutStop(kDc1, 18);
    }
    radw = hdc_cath_thick / hdc_cath_radlen;
    musc(radw);

    // at last cathode foil of first drift chamber set, drift to next
    drift = hdc_2_zpos - hdc_1_zpos - hdc_nr_plan * hdc_del_plane;
    project(drift);
    radw = drift / hair_radlen;
    muscExt(radw, drift);

    chamber(2);
    if (t.xs > (hdc_2_bot - hdc_2x_offset) || t.xs < (hdc_2_top - hdc_2x_offset) ||
        t.ys > (hdc_2_left - hdc_2y_offset) || t.ys < (hdc_2_right - hdc_2y_offset)) {
        if (use_det_cut) return hutStop(kDc2, 19);
    }
    radw = hdc_cath_thick / hdc_cath_radlen;
    musc(radw);

    // drift to the 1st hodoscope
    drift = hscin_1x_zpos - hdc_2_zpos - 0.5 * hdc_nr_plan * hdc_del_plane;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (t.ys > (hscin_1x_left + hscin_1y_offset) || t.ys < (hscin_1x_right + hscin_1y_offset)) {
        if (use_det_cut) return hutStop(kS1, 20);
    }
    radw = hscin_1x_thick / hscin_radlen;
    musc(radw);
    drift = hscin_1y_zpos - hscin_1x_zpos;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (t.xs > (hscin_1y_bot + hscin_1x_offset) || t.xs < (hscin_1y_top + hscin_1x_offset)) {
        if (use_det_cut) return hutStop(kS1, 21);   // FORTRAN-QUIRK: both S1 planes count in s1
    }
    radw = hscin_1y_thick / hscin_radlen;
    musc(radw);

    // drift to the second cherenkov
    drift = hcer_2_zentrance - hscin_1y_zpos;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);

    radw = hcer_2_entr_thick / hcer_2_entr_radlen;
    musc(radw);

    drift = hcer_2_zmirror - hcer_2_zentrance;
    radw = drift / hcer_2_radlen;
    project(drift);
    muscExt(radw, drift);

    radw = hcer_mir_thick / hcer_mir_radlen;
    musc(radw);

    drift = hcer_2_zexit - hcer_2_zmirror;
    radw = drift / hcer_2_radlen;
    project(drift);
    muscExt(radw, drift);

    radw = hcer_2_exit_thick / hcer_2_exit_radlen;
    musc(radw);

    // drift to 2nd hodoscope
    drift = hscin_2x_zpos - hcer_2_zexit;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (t.ys > (hscin_2x_left + hscin_2y_offset) || t.ys < (hscin_2x_right + hscin_2y_offset)) {
        if (use_det_cut) return hutStop(kS3, 22);   // FORTRAN-QUIRK: the 2x plane counts in s3
    }
    radw = hscin_2x_thick / hscin_radlen;
    musc(radw);
    drift = hscin_2y_zpos - hscin_2x_zpos;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (t.xs > (hscin_2y_bot + hscin_2x_offset) || t.xs < (hscin_2y_top + hscin_2x_offset)) {
        if (use_det_cut) return hutStop(kS2, 23);   // FORTRAN-QUIRK: the 2y plane counts in s2
    }
    radw = hscin_2y_thick / hscin_radlen;
    musc(radw);

    // drift to the back of the calorimeter; true position cut
    drift = hcal_4ta_zpos - hscin_2y_zpos;
    radw = drift / hair_radlen;
    project(drift);
    muscExt(radw, drift);
    if (t.ys > hcal_left || t.ys < hcal_right || t.xs > hcal_bottom || t.xs < hcal_top) {
        if (use_det_cut) return hutStop(kCal, 24);
    }

    // fit track to give new focal plane values, use LFIT from GENLIB
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
    float x0, y0, dx, dy, badf;
    cern::lfit(zdc, xdc, 12, 0, dx, x0, badf);
    cern::lfit(zdc, ydc, 12, 0, dy, y0, badf);
    a.x_fp = static_cast<double>(x0);
    a.y_fp = static_cast<double>(y0);
    a.dx_fp = static_cast<double>(dx);
    a.dy_fp = static_cast<double>(dy);
    det.fit_done = true;
    det.xfp4 = x0;
    det.dxfp4 = dx;
    det.yfp4 = y0;
    det.dyfp4 = dy;

    // FORTRAN-QUIRK: `cal` also counts the fiducial cut on the FITTED track 5 cm inside the calorimeter
    const double xcal = a.x_fp + a.dx_fp * hcal_4ta_zpos;
    const double ycal = a.y_fp + a.dy_fp * hcal_4ta_zpos;
    if (ycal > (hcal_left - 5.0) || ycal < (hcal_right + 5.0) || xcal > (hcal_bottom - 5.0) ||
        xcal < (hcal_top + 5.0)) {
        if (use_det_cut) return hutStop(kCal, 25);
    }
    return kSuccess;   // ok_hut = .true.
}

}  // namespace shms
}  // namespace simc
