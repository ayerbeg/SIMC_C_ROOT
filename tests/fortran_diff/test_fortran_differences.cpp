// tests/fortran_diff/test_fortran_differences.cpp — task P7.5 (D-031, D-032).
//
// One sub-test per item of docs/FORTRAN_DIFFERENCES.md that changes behaviour. For each corrected clear Fortran bug
// (class a) the test runs the C++ code twice: with fortran_compat = true it must give the Fortran result (taken from
// the Fortran source/reference where possible), with the default (corrected) the corrected result.
//
//   usage: test_fortran_differences <item> <reference dir> <simc_gfortran dir>
// Exit 0 pass, 1 fail, 77 skip (reference data missing).
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "simc/core/Json.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/SpectrometerId.h"
#include "simc/core/random/Mt64Random.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/driver/RunReport.h"
#include "simc/driver/SimcRun.h"
#include "simc/input/InpReader.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/optics/CosyEngine.h"
#include "simc/output/SummaryJson.h"
#include "simc/physics/proton/ProtonCrossSections.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/setup/SetupBuilder.h"
#include "simc/spectrometers/hms/kernel/McHms.h"

using namespace simc;

namespace {

int g_checks = 0, g_fail = 0;
void expect(bool ok, const std::string& what) {
    ++g_checks;
    if (!ok) ++g_fail;
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what.c_str());
}

std::string g_ref, g_fortran;

// ------------------------------------------------------------------------------------------------ D-032
// One number per spectrometer; the retired code 6 is rejected at input (.inp and JSON) and by the set-up.
int testSpectrometerId() {
    const auto& t = SpectrometerTable();
    expect(t.size() == 7, "table has the 7 spectrometers");
    for (const auto& r : t) {
        expect(SpectrometerFromCode(r.code) == r.id && Code(r.id) == r.code, std::string("code round trip ") + r.key);
        expect(SpectrometerFromKey(r.key) == r.id, std::string("key round trip ") + r.key);
        expect(ArmCodeError(r.code, "electron").empty(), std::string("code accepted ") + r.key);
    }
    // dbase.f:247-263: phi = 3pi/2 for HMS, HRS-R, calo HMS side
    expect(IsBeamRight(SpectrometerId::HMS) && IsBeamRight(SpectrometerId::HRSR) &&
               IsBeamRight(SpectrometerId::CALO_HMS_SIDE) && !IsBeamRight(SpectrometerId::SHMS) &&
               !IsBeamRight(SpectrometerId::SOS) && !IsBeamRight(SpectrometerId::HRSL) &&
               !IsBeamRight(SpectrometerId::CALO_SHMS_SIDE),
           "beam side of every spectrometer as dbase.f:247-263");
    expect(!SpectrometerFromCode(6) && !SpectrometerFromCode(0) && !SpectrometerFromCode(9), "6, 0, 9 are not ids");
    expect(ArmCodeError(6, "electron") == "electron_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS",
           "message for electron_arm = 6");
    expect(ArmCodeError(6, "hadron") == "hadron_arm = 6 was the retired SHMS LSA tune; use 5 = SHMS",
           "message for hadron_arm = 6");

    // JSON input
    for (const char* which : {"electron_arm", "hadron_arm"}) {
        SimcParameters p;
        Json j = {{"electron_arm", 5}, {"hadron_arm", 1}};
        j[which] = 6;
        std::string what;
        try {
            FromJson(j, p);
        } catch (const std::invalid_argument& e) {
            what = e.what();
        }
        expect(what.find(std::string(which) + " = 6 was the retired SHMS LSA tune; use 5 = SHMS") != std::string::npos,
               std::string("FromJson rejects ") + which + " = 6: '" + what + "'");
        SimcParameters ok;
        bool threw = false;
        try {
            FromJson(Json{{"electron_arm", 5}, {"hadron_arm", 1}}, ok);
        } catch (const std::exception&) {
            threw = true;
        }
        expect(!threw && ok.electron_arm == 5 && ok.hadron_arm == 1, "FromJson accepts 5/1");
    }

    // .inp input
    const std::string dir = "fortran_diff_d032";
    std::system(("mkdir -p " + dir).c_str());
    for (const char* which : {"electron_arm", "hadron_arm"}) {
        const std::string f = dir + "/" + which + "6.inp";
        {
            std::ofstream o(f);
            o << "begin parm x\n  electron_arm = " << (std::string(which) == "electron_arm" ? 6 : 5)
              << "\n  hadron_arm = " << (std::string(which) == "hadron_arm" ? 6 : 1) << "\nend parm x\n";
        }
        SimcParameters p;
        InpReader rd(p);
        const bool ok = rd.loadRun(f);
        expect(!ok && rd.errors().size() == 1 &&
                   rd.errors()[0] == std::string(which) + " = 6 was the retired SHMS LSA tune; use 5 = SHMS",
               std::string("InpReader rejects ") + which + " = 6");
    }
    // set-up from parameters built in code (both modes)
    for (bool compat : {false, true}) {
        SimcParameters p;
        p.electron_arm = 5;
        p.hadron_arm = 6;
        SetupOptions so;
        so.fortran_compat = compat;
        const RunSetup s = BuildRunSetup(p, so);
        expect(!s.ok() && s.stops.front().find("hadron_arm = 6 was the retired SHMS LSA tune") != std::string::npos,
               std::string("BuildRunSetup stops for hadron_arm = 6, fortran_compat=") + (compat ? "true" : "false"));
    }
    return 0;
}


bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

// ------------------------------------------------------------------------------------------------ D-006 pass_cuts
// simc.f:237 tests the hadron delta against SPedge%e%delta%max.
int testPassCuts() {
    RunSetup s;
    s.SPedge.e.delta = {-10.0, 22.0};   // SHMS-like electron arm
    s.SPedge.p.delta = {-10.0, 10.0};   // HMS-like hadron arm
    Event r;                            // every other recon value is 0: inside the +-1e10 defaults
    r.p.delta = 15.0;                   // above the hadron limit, below the electron limit
    s.fortran_compat = true;
    expect(driver::PassCuts(s, r), "compat: p delta 15% passes (tested against SPedge%e%delta%max = 22%)");
    s.fortran_compat = false;
    expect(!driver::PassCuts(s, r), "corrected: p delta 15% fails (SPedge%p%delta%max = 10%)");
    r.p.delta = 5.0;
    for (bool c : {true, false}) {
        s.fortran_compat = c;
        expect(driver::PassCuts(s, r), std::string("p delta 5% passes in both modes, compat=") + (c ? "1" : "0"));
    }
    r.p.delta = 25.0;
    for (bool c : {true, false}) {
        s.fortran_compat = c;
        expect(!driver::PassCuts(s, r), std::string("p delta 25% fails in both modes, compat=") + (c ? "1" : "0"));
    }
    return 0;
}

// ------------------------------------------------------------------------------------------------ D-017 target.f:501
// extreme_trip_thru_target, hadron, pudding can / cryo17: zz from the%min (Fortran) vs thp%min (corrected).
int testExtremeTrip() {
    const std::string dir = g_ref + "/functions/extreme_trip_thru_target";
    simcl1::Table in, out;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !out.load(dir + "/output.txt", err)) {
        std::printf("SKIP: %s\n", err.c_str());
        return 77;
    }
    // Run extreme_trip_thru_target for one table row; `thpShift` moves the hadron angle range (synthetic variant).
    auto run = [&](const simcl1::Row& r, bool compat, double thpMin, double thpMax, double& zzHad) {
        material::HallCCryoTarget T;
        std::ostringstream log;
        T.setLog(&log);
        auto& targ = T.targ();
        targ.A = r["targ_A"];
        targ.Z = r["targ_Z"];
        targ.N = targ.A - targ.Z;
        targ.rho = r["targ_rho"];
        targ.thick = r["targ_thick_mg"] / 1000.;
        targ.length = targ.thick / targ.rho;
        targ.angle = r["targ_angle_deg"] / constants::degrad;
        targ.can = static_cast<int>(std::lround(r["targ_can"]));
        T.setArms(static_cast<int>(std::lround(r["electron_arm"])), static_cast<int>(std::lround(r["hadron_arm"])));
        const double Mh = r["Mh"];
        T.targetInit(std::lround(r["using_Eloss"]) > 0, std::lround(r["using_Coulomb"]) > 0,
                     std::fabs(r["the_cent_deg"]) / constants::degrad, std::fabs(r["thp_cent_deg"]) / constants::degrad,
                     r["Pp_cent"], Mh * Mh, r["Ebeam"], r["Pe_cent"]);
        const material::Limits the{r["the_min"], r["the_max"]}, thp{thpMin, thpMax};
        const material::Limits pe{r["pe_min"], r["pe_max"]}, pp{r["pp_min"], r["pp_max"]};
        const material::Limits z{-0.5 * targ.length, 0.5 * targ.length};
        T.extremeTripThruTarget(r["Ebeam"], the, thp, pe, pp, z, Mh, compat);
        const auto result = T.targ();
        // independent recomputation of the hadron maximum (target.f:497-505) from z = max(-L/2/tan(angle), -L/2), with
        // the angle the%min (Fortran, compat) or thp%min (corrected)
        double zz = -(targ.length / 2.) / std::tan(compat ? the.min : thp.min);
        zz = std::max(zz, (-targ.length / 2.));
        zzHad = zz;
        ScriptedRandom none;
        const double m = Mh;
        const double emin = std::sqrt(pp.min * pp.min + m * m), emax = std::sqrt(pp.max * pp.max + m * m);
        double E1, E2, t1, t2;
        T.trip(3, zz, emin, thp.min, E1, t1, m, 3, none);
        T.trip(3, zz, emax, thp.min, E2, t2, m, 3, none);
        expect(sameBits(result.Eloss[2].max, std::max(E1, E2)) && sameBits(result.teff[2].max, std::max(t1, t2)),
               std::string(compat ? "  compat" : "  corrected") + ": hadron Eloss/teff max = trip from z = " +
                   std::to_string(zz) + " cm (angle " + (compat ? "the_min" : "thp_min") + ")");
        return result;
    };
    int ncan23 = 0, nEffect = 0;
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        const simcl1::Row r(in, i), o(out, i);
        const int can = static_cast<int>(std::lround(r["targ_can"]));
        if (can != 2 && can != 3) continue;
        ++ncan23;
        std::printf("row %zu (can %d): the_min %.4f rad, thp_min %.4f rad\n", i + 1, can, r["the_min"], r["thp_min"]);
        double zc, zr;
        const auto tc = run(r, true, r["thp_min"], r["thp_max"], zc);
        const auto tr = run(r, false, r["thp_min"], r["thp_max"], zr);
        expect(sameBits(tc.Eloss[2].max, o["Eloss_max_3"]) && sameBits(tc.teff[2].max, o["teff_max_3"]),
               "  compat Eloss_max_3/teff_max_3 = Fortran table bitwise");
        // both angles < 45 deg: -L/2/tan(angle) < -L/2, so z is clamped to -L/2 with either angle -> no effect here
        expect(zc == zr && sameBits(tr.Eloss[2].max, tc.Eloss[2].max),
               "  both angles < 45 deg -> z = -L/2 either way, no effect on this row");
        // synthetic variant: hadron at 52-60 deg (thp_min > 45 deg) -> the bug matters
        const double a0 = 52. / constants::degrad, a1 = 60. / constants::degrad;
        const auto sc = run(r, true, a0, a1, zc);
        const auto sr = run(r, false, a0, a1, zr);
        expect(zc != zr, "  hadron at 52-60 deg: z of the maximum-loss trip " + std::to_string(zc) + " (compat) vs " +
                             std::to_string(zr) + " cm (corrected); Eloss_max_3 " + std::to_string(sc.Eloss[2].max) +
                             " vs " + std::to_string(sr.Eloss[2].max) + " MeV");
        if (!sameBits(sc.Eloss[2].max, sr.Eloss[2].max)) ++nEffect;
        bool restSame = true;
        for (int k = 0; k < 3; ++k) {
            restSame = restSame && sameBits(sc.Eloss[k].min, sr.Eloss[k].min) && sameBits(sc.teff[k].min, sr.teff[k].min);
            if (k < 2)
                restSame = restSame && sameBits(sc.Eloss[k].max, sr.Eloss[k].max) &&
                           sameBits(sc.teff[k].max, sr.teff[k].max) && sameBits(sc.musc_max[k], sr.musc_max[k]);
        }
        expect(restSame, "  all other extreme values identical in both modes");
    }
    expect(ncan23 > 0, "the table has pudding-can/cryo17 rows (" + std::to_string(ncan23) + ")");
    expect(nEffect > 0, "the corrected hadron maximum loss differs for a hadron angle > 45 deg (" +
                            std::to_string(nEffect) + " of " + std::to_string(ncan23) + " rows)");
    return 0;
}

// ------------------------------------------------------------------------------------------------ D-014 deForest
// physics_proton.f:84-87: the cos_phi warning is tested after the clamp (never prints).
int testDeForestWarning() {
    Event ev;
    ev.Ein = 4000.;
    ev.e.E = 3000.;
    ev.e.theta = 0.3;
    ev.nu = 1000.;
    ev.q = 1500.;
    ev.Q2 = 1.2e6;
    ev.Pm = 100.;
    ev.p.E = 1500.;
    ev.uq = {0.6, 0.0, 0.8};
    ev.up = {-0.5, 0.1, 0.9};   // not a unit vector: the cos_phi formula then exceeds 1 (1.0358)
    double out[2];
    std::string logs[2];
    for (int compat = 0; compat < 2; ++compat) {
        proton::DeForestState st;
        st.deForest_flag = 0;
        st.Mh2 = constants::Mp * constants::Mp;
        st.nevent = 42;
        st.fortran_compat = compat == 1;
        std::ostringstream log;
        st.log = &log;
        out[compat] = proton::deForest(ev, st);
        logs[compat] = log.str();
    }
    expect(logs[1].find("deForest give cos_phi") == std::string::npos, "compat: no cos_phi warning (as the Fortran)");
    expect(logs[0].find("WARNING: deForest give cos_phi =") != std::string::npos &&
               logs[0].find("1.035") != std::string::npos && logs[0].find("at event") != std::string::npos,
           "corrected: the warning prints the unclamped value: '" + logs[0].substr(0, logs[0].find('\n')) + "'");
    expect(sameBits(out[0], out[1]) && std::isfinite(out[0]), "the cross section is identical in both modes");
    return 0;
}

// ------------------------------------------------------------------------------------------------ D-028 F2
// event.f:1641 spectrometer_angles: sqrt of a rounding-negative argument -> NaN.
int testSpectrometerAngles() {
    const std::string dir = g_ref + "/functions/spectrometer_angles";
    simcl1::Table in, out;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !out.load(dir + "/output.txt", err)) {
        std::printf("SKIP: %s\n", err.c_str());
        return 77;
    }
    std::size_t nanF = 0, compatNanOk = 0, corrFinite = 0, otherSame = 0, other = 0;
    double maxCorrDy = 0;
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        const simcl1::Row r(in, i), o(out, i);
        double dxa, dya, dxb, dyb;
        SpectrometerAngles(r["theta0"], r["phi0"], dxa, dya, r["theta"], r["phi"], true);
        SpectrometerAngles(r["theta0"], r["phi0"], dxb, dyb, r["theta"], r["phi"], false);
        if (std::isnan(o["dy"])) {
            ++nanF;
            if (std::isnan(dya) && sameBits(dxa, o["dx"])) ++compatNanOk;
            if (std::isfinite(dyb) && sameBits(dxb, dxa)) {
                ++corrFinite;
                maxCorrDy = std::max(maxCorrDy, std::fabs(dyb));
            }
        } else {
            ++other;
            if (sameBits(dxa, dxb) && sameBits(dya, dyb) && sameBits(dya, o["dy"])) ++otherSame;
        }
    }
    expect(nanF > 0, "the Fortran table has NaN dy rows: " + std::to_string(nanF));
    expect(compatNanOk == nanF, "compat: NaN in every Fortran NaN row (" + std::to_string(compatNanOk) + ")");
    expect(corrFinite == nanF && maxCorrDy == 0.0,
           "corrected: dy = 0 (finite) in every Fortran NaN row (" + std::to_string(corrFinite) + "), max |dy| = " +
               std::to_string(maxCorrDy));
    expect(otherSame == other, "all other rows bitwise identical to the Fortran in both modes (" +
                                   std::to_string(otherSame) + "/" + std::to_string(other) + ")");
    return 0;
}

// ------------------------------------------------------------------------------------------------ D-015 F9 raddeg
int testRotateRaddeg() {
    optics::Track t;
    t.dxdzs = 0.03;
    t.dydzs = -0.02;
    const double ang = -6.0;   // HMS dipole rotation (mc_hms.f)
    for (int compat = 0; compat < 2; ++compat) {
        const double rd = compat ? 0.017453292 : constants::pi / 180.;
        const double th = ang * rd;
        const double tan_th = std::tan(th), sin_th = std::sin(th), cos_th = std::cos(th);
        const double alpha_p = (t.dxdzs + tan_th) / (1. - t.dxdzs * tan_th);
        const double beta_p = t.dydzs / (cos_th - t.dxdzs * sin_th);
        const double x0 = 30.0, y0 = 5.0;
        double x = x0, y = y0;
        optics::rotateHaxis(t, ang, x, y, compat == 1);
        expect(sameBits(x, x0 * (cos_th + alpha_p * sin_th)) && sameBits(y, y0 + x0 * beta_p * sin_th),
               std::string("rotate_haxis with raddeg = ") + (compat ? "0.017453292 (compat)" : "pi/180 (corrected)"));
    }
    double xa = 30, ya = 5, xb = 30, yb = 5;
    optics::rotateVaxis(t, 1.5, xa, ya, true);
    optics::rotateVaxis(t, 1.5, xb, yb, false);
    expect(!sameBits(ya, yb) && std::fabs(ya - yb) < 1e-6, "rotate_vaxis: the modes differ at the 1e-8 level (" +
                                                               std::to_string(ya - yb) + " cm)");
    return 0;
}

// ------------------------------------------------------------------------------------------------ D-021 slit counters
// mc_hms_coll.f increments hSTOP_slit_* once per 1/20 step; corrected: once per ray.
int testHmsCollSlitCounters() {
    const double mpi2 = constants::Mpi * constants::Mpi;
    long counts[2][3] = {{0}};
    long ok[2] = {0, 0};
    long maxPerRay[2] = {0, 0};
    long totalRays = 2000;
    std::vector<int> stops[2];
    for (int compat = 0; compat < 2; ++compat) {
        hms::McHms h;
        h.setFortranCompat(compat == 1);
        std::ostringstream log;
        h.setLog(&log);
        Mt64Random rng(2025);
        Mt64Random geo(7);
        for (long i = 0; i < totalRays; ++i) {
            hms::McHmsArgs a;
            a.p_spec = 2000.;
            a.dpp = 0.;
            a.dxdz = (geo.uniform() - 0.5) * 0.16;
            a.dydz = (geo.uniform() - 0.5) * 0.08;
            a.m2 = mpi2;
            a.using_HMScoll = true;
            a.ms_flag = true;
            const hms::HmsCounters before = h.counters();
            const hms::McHmsDetail d = h.mcHms(a, rng);
            stops[compat].push_back(d.stop_code);
            if (a.ok_spec) ++ok[compat];
            for (int k = 0; k < 3; ++k) {
                const long inc = h.counters()[hms::kSlitHor + k] - before[hms::kSlitHor + k];
                maxPerRay[compat] = std::max(maxPerRay[compat], inc);
            }
        }
        for (int k = 0; k < 3; ++k) counts[compat][k] = h.counters()[hms::kSlitHor + k];
    }
    std::printf("  slit hor/vert/oct: compat %ld %ld %ld, corrected %ld %ld %ld; successes %ld / %ld\n", counts[1][0],
                counts[1][1], counts[1][2], counts[0][0], counts[0][1], counts[0][2], ok[1], ok[0]);
    expect(stops[0] == stops[1] && ok[0] == ok[1], "the transport (stop reasons, successes, RNG use) is identical");
    expect(maxPerRay[1] > 1, "compat: some ray is counted more than once per counter (max " +
                                 std::to_string(maxPerRay[1]) + ", as mc_hms_coll.f)");
    expect(maxPerRay[0] == 1, "corrected: at most one count per ray and counter");
    expect(counts[0][0] < counts[1][0] || counts[0][1] < counts[1][1] || counts[0][2] < counts[1][2],
           "corrected totals are smaller");
    return 0;
}

// ------------------------------------------------------------------------------------------------ D-006 transp.f:158
// K -> pi pi decay in the first half of a transp element: m_final = Mk (Fortran) vs Mpi (corrected).
int testKaonDecayMass() {
    double m2[2], ph[2];
    for (int compat = 0; compat < 2; ++compat) {
        optics::CosyEngine eng;
        eng.setFortranCompat(compat == 1);
        std::ostringstream log;
        eng.setLog(&log);
        eng.initForward(1, hms::McHms::defaultMatrixDir());
        ScriptedRandom rng;
        // z_decay = -dlen*ln(1-u) tiny -> first half; rph, rth1; kaon branch u = 0.9 >= 0.7 -> two pions
        rng.load({}, {1.e-9, 0.25, 0.3, 0.9});
        optics::DecayState d;
        d.decay_flag = true;
        d.m2 = constants::Mk * constants::Mk;
        d.ph = 1500.;
        d.ctau = 371.2;
        d.Mh2_final = d.m2;
        d.rng = &rng;
        optics::Track t;
        double path = 0;
        eng.transp(1, 1, t, 10.0, path, &d);
        m2[compat] = d.m2;
        ph[compat] = d.ph;
        expect(d.dflag && rng.uniformsLeft() == 0, std::string("decay in the first half, 2-pion mode, compat=") +
                                                       (compat ? "1" : "0"));
    }
    expect(sameBits(m2[1], constants::Mk * constants::Mk), "compat: the decay product gets the kaon mass (transp.f:158)");
    expect(sameBits(m2[0], constants::Mpi * constants::Mpi), "corrected: the decay product is a pion (Mpi^2)");
    expect(ph[0] != ph[1], "the decay momenta differ (" + std::to_string(ph[1]) + " vs " + std::to_string(ph[0]) + " MeV/c)");
    return 0;
}

// ------------------------------------------------------------------------------------------------ D-011/D-018 report
// .geni SHMS labels (simc.f:519-520, 526) and .hist formats (simc.f:815-818, 897-898, 907, 927, 727-732).
int testReportFormats() {
    const std::string inp = g_ref + "/runs/eep_h_elastic_rad_A/eep_h_elastic_rad_A.inp";
    if (!std::ifstream(inp)) {
        std::printf("SKIP: %s missing\n", inp.c_str());
        return 77;
    }
    SimcParameters p;
    InpReader rd(p);
    rd.loadRun(inp);
    p.spec.e.offset.xptar = 1.5;   // mrad (input units)
    p.spec.p.offset.yptar = -2.0;
    p.use_first_cer = 1;
    driver::RunReport r;
    StopCounters sh;
    sh.spectrometer = "shms";
    sh.role = "electron";
    const std::vector<std::string> keys = StopCounterKeys("shms");
    for (const std::string& k : keys) sh.counts[k] = 0;
    sh.counts["q3_mid"] = 0;
    sh.counts["trials"] = 100;
    sh.counts["q2_mid"] = 11;
    sh.counts["q3_mid"] = 33;
    sh.counts["dc1"] = 4;
    sh.counts["dc2"] = 5;
    sh.counts["successes"] = 50;
    r.stops.push_back(sh);
    std::string hist[2], geni[2];
    for (int compat = 0; compat < 2; ++compat) {
        SetupOptions so;
        so.fortran_compat = compat == 1;
        so.data_dir = g_fortran;
        const RunSetup s = BuildRunSetup(p, so);
        geni[compat] = driver::FormatGeni(s, r);
        hist[compat] = driver::FormatHist(s, r);
    }
    auto line = [](const std::string& t, const std::string& key) {
        const auto i = t.find(key);
        if (i == std::string::npos) return std::string("<missing>");
        const auto b = t.rfind('\n', i) + 1;
        return t.substr(b, t.find('\n', i) - b);
    };
    // .geni
    const std::string q3c = line(geni[1], "Q3 phys"), q3 = line(geni[0], "Q3 phys");
    expect(q3c.find(" 11 ") != std::string::npos && q3c.find(" 33") == std::string::npos,
           "compat .geni Q3 mid column = shmsSTOP_q2_mid: '" + q3c + "'");
    expect(q3.find(" 33 ") != std::string::npos && q3.find(" 11 ") == std::string::npos,
           "corrected .geni Q3 mid column = shmsSTOP_q3_mid: '" + q3 + "'");
    expect(geni[1].find("DC1, DC2, Scin, Cal") != std::string::npos && line(geni[0], "DC1, DC2").find("Scin") == std::string::npos,
           "SHMS DC line: compat label 'DC1, DC2, Scin, Cal', corrected 'DC1, DC2': '" + line(geni[0], "DC1, DC2") + "'");
    // .hist
    const std::string offc = line(hist[1], "xptar offset"), off = line(hist[0], "xptar offset");
    expect(offc.find("0.0015") != std::string::npos, "compat xptar offset printed in rad under 'mr': '" + offc + "'");
    expect(off.find("1.5000") != std::string::npos, "corrected xptar offset printed in mr: '" + off + "'");
    expect(line(hist[1], "deForest_flag").find("deForest_flag= F") != std::string::npos &&
               line(hist[0], "deForest_flag").find("deForest_flag= 0") != std::string::npos,
           "deForest_flag: compat L format 'F', corrected integer '0'");
    expect(line(hist[1], "use_first_cer").find("use_first_cer= T") != std::string::npos &&
               line(hist[0], "use_first_cer").find("use_first_cer= 1") != std::string::npos,
           "use_first_cer: compat 'T', corrected '1'");
    expect(line(hist[1], "use_expon").find("use_expon=  F") != std::string::npos &&
               line(hist[0], "use_expon").find("use_expon=  0") != std::string::npos,
           "use_expon: compat 'F', corrected '0'");
    // simc.f:727-732: the (e,e'p)pi title goes to unit 6 in the Fortran, to the .hist when corrected
    for (int compat = 0; compat < 2; ++compat) {
        SetupOptions so;
        so.fortran_compat = compat == 1;
        so.data_dir = g_fortran;
        RunSetup s = BuildRunSetup(p, so);
        s.flags.doing_eep = s.flags.doing_hyd_elast = false;
        s.flags.doing_delta = s.flags.doing_hyddelta = true;
        std::ostringstream cap;
        std::streambuf* old = std::cout.rdbuf(cap.rdbuf());
        const std::string h = driver::FormatHist(s, r);
        std::cout.rdbuf(old);
        const bool inHist = h.find("H(e,e'p)pi") != std::string::npos;
        const bool onStdout = cap.str().find("H(e,e'p)pi") != std::string::npos;
        expect(compat ? (!inHist && onStdout) : (inHist && !onStdout),
               std::string("(e,e'p)pi title: ") + (compat ? "compat on unit 6 only" : "corrected in the .hist only"));
    }
    // everything else identical (the lines that differ are only the ones above)
    std::istringstream a(hist[1]), b(hist[0]);
    std::string la, lb;
    int ndiff = 0;
    while (std::getline(a, la) && std::getline(b, lb))
        if (la != lb) ++ndiff;
    expect(ndiff == 5, ".hist: exactly 5 lines differ between the modes (xptar/yptar offsets, deForest_flag, "
                       "use_first_cer, use_expon): " + std::to_string(ndiff));
    return 0;
}

const std::map<std::string, std::function<int()>>& items() {
    static const std::map<std::string, std::function<int()>> m = {
        {"D032_spectrometer_id", testSpectrometerId},
        {"D006_pass_cuts", testPassCuts},
        {"D017_extreme_trip_hadron_angle", testExtremeTrip},
        {"D014_deforest_cos_phi_warning", testDeForestWarning},
        {"D028_spectrometer_angles_nan", testSpectrometerAngles},
        {"D015_rotate_raddeg", testRotateRaddeg},
        {"D021_hms_coll_slit_counters", testHmsCollSlitCounters},
        {"D006_kaon_decay_mass", testKaonDecayMass},
        {"D011_D018_report_formats", testReportFormats},
    };
    return m;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <item> [reference dir] [simc_gfortran dir]\n", argv[0]);
        return 2;
    }
    g_ref = argc > 2 ? argv[2] : "";
    g_fortran = argc > 3 ? argv[3] : "";
    const auto it = items().find(argv[1]);
    if (it == items().end()) {
        std::fprintf(stderr, "unknown item %s\n", argv[1]);
        return 2;
    }
    const int rc = it->second();
    if (rc == 77) return 77;
    std::printf("%s: %d checks, %d failures\n", argv[1], g_checks, g_fail);
    return g_fail ? 1 : 0;
}
