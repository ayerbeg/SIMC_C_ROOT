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
#include "simc/core/Registry.h"
#include "simc/core/SpectrometerId.h"
#include "simc/core/random/Mt64Random.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/driver/FortranWriter.h"
#include "simc/driver/RunReport.h"
#include "simc/driver/SimcRun.h"
#include "simc/input/InpReader.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/optics/CosyEngine.h"
#include "simc/output/NtupleLayout.h"
#include "simc/output/SummaryJson.h"
#include "simc/physics/delta/DeltaModel.h"
#include "simc/physics/nuclear/SpectralFunctionTables.h"
#include "simc/radiation/SimcRadiator.h"
#include "simc/reactions/NuclearEep.h"
#include "simc/physics/proton/ProtonCrossSections.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/setup/SetupBuilder.h"
#include "simc/spectrometers/hms/kernel/McHms.h"
#include "simc/spectrometers/shms/kernel/McShms.h"
#include "simc/spectrometers/hrs/kernel/McHrs.h"
#include "simc/stages/IOutputSink.h"

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


// ------------------------------------------------------------------------------------------------ D-033 (P8.1/P8.2)
std::string spectralDir() { return std::string(SIMC_SOURCE_DIR) + "/data/spectral"; }

// F3, sf_lookup.f:93: SFd = SF/4/3.1415926535/(Pm*Pm)/5.0/20.0 (Fortran) vs the per-node differential spectral function
// with the widths of the file (corrected).
int testBenharSfDiff() {
    const double pi = constants::pi;
    for (const char* file : {"benharsf_56.dat", "benharsf_4.dat"}) {
        nuclear::SfTable c, k;   // c = compat, k = corrected
        std::string err;
        if (!c.load(spectralDir() + "/" + file, true, nullptr, &err) || !k.load(spectralDir() + "/" + file, true, nullptr, &err)) {
            std::printf("skip: %s\n", err.c_str());
            return 77;
        }
        c.setFortranCompat(true);
        k.setFortranCompat(false);
        const std::string f = file;
        // a node with strength
        int iPm = 5, iEm = 3;
        for (int j = 1; j <= c.numEm(); ++j)
            if (c.sfval(j, iPm) > c.sfval(iEm, iPm)) iEm = j;
        const double Pm = c.Pmval(iPm), Em = c.Emval(iEm);
        const double fortran = c.sfval(iEm, iPm) / 4 / 3.1415926535 / (Pm * Pm) / 5.0 / 20.0;
        expect(c.lookupDiff(Em, Pm) == fortran, f + " compat: SFd = SF/4/3.1415926535/Pm^2/5/20 at a node");
        const double node = c.sfval(iEm, iPm) / (4. * pi * Pm * Pm * c.dPm(iPm) * c.dEm(iEm));
        expect(std::abs(k.lookupDiff(Em, Pm) / node - 1) < 1e-14,
               f + " corrected: SFd = sfval/(4 pi Pmval^2 dPm dEm) at a node, widths from the file (dPm " +
                   std::to_string(c.dPm(iPm)) + ", dEm " + std::to_string(c.dEm(iEm)) + ")");
        // Pm = 0: +Inf in the Fortran, finite (the first node's density) when corrected
        const double c0 = c.lookupDiff(Em, 0.0), k0 = k.lookupDiff(Em, 0.0);
        expect(std::isinf(c0) || std::isnan(c0), f + " compat: Pm = 0 gives Inf/NaN (" + std::to_string(c0) + ")");
        expect(std::isfinite(k0) && k0 == k.lookupDiff(Em, c.Pmval(1)),
               f + " corrected: Pm = 0 gives the first node's value (" + std::to_string(k0) + ")");
        // no 1/Pm^2 divergence below the first node
        expect(k.lookupDiff(Em, 0.5 * c.Pmval(1)) == k0 && c.lookupDiff(Em, 0.5 * c.Pmval(1)) > 3.9 * c.lookupDiff(Em, c.Pmval(1)),
               f + " below the first node: compat grows as 1/Pm^2, corrected flat");
        // normalisation of the corrected density over the nodes: sum S * 4 pi P^2 dP dE = sum of the (non-negative)
        // bin probabilities = 1 (benharsf_4.dat has one negative node, -1.454e-6, set to 0 by the SF < 1e-20 test)
        double sum = 0, pos = 0;
        for (int i = 1; i <= k.numPm(); ++i)
            for (int j = 1; j <= k.numEm(); ++j) {
                sum += k.lookupDiff(k.Emval(j), k.Pmval(i)) * 4. * pi * k.Pmval(i) * k.Pmval(i) * k.dPm(i) * k.dEm(j);
                pos += k.sfval(j, i) < 1e-20 ? 0. : k.sfval(j, i);
            }
        char buf[160];
        std::snprintf(buf, sizeof buf, "%.15f (sum of the non-negative bin probabilities %.15f)", sum, pos);
        expect(std::abs(sum / pos - 1) < 1e-12, f + " corrected density integrates over the table to " + buf);
    }
    return 0;
}

// F6, sf_lookup.f:160-184: Em == Emval(numEm) uses the SAVEd Em1/Em2/sf1/sf2 of the previous call (NaN on the first).
int testSfLookupStale() {
    nuclear::SfTable c, k;
    std::string err;
    if (!c.load(spectralDir() + "/benharsf_12.dat", true, nullptr, &err) ||
        !k.load(spectralDir() + "/benharsf_12.dat", true, nullptr, &err))
        return 77;
    c.setFortranCompat(true);
    k.setFortranCompat(false);
    const int n = c.numEm(), iPm = 4;
    const double Pm = c.Pmval(iPm), Elast = c.Emval(n);
    expect(std::isnan(c.lookup(Elast, Pm)), "compat: first call at Em = Emval(numEm) is NaN (unset SAVEd locals)");
    expect(k.lookup(Elast, Pm) == (c.sfval(n, iPm) < 1e-20 ? 0. : c.sfval(n, iPm)),
           "corrected: Em = Emval(numEm) gives the table value of the last node");
    // after a call in the first interval, the compat value is the extrapolation of THAT interval to Elast
    const double E1 = 0.5 * (c.Emval(1) + c.Emval(2));
    c.lookup(E1, Pm);
    k.lookup(E1, Pm);
    const double s1 = c.sfval(1, iPm), s2 = c.sfval(2, iPm);
    double stale = s1 + (Elast - c.Emval(1)) * (s2 - s1) / (c.Emval(2) - c.Emval(1));
    if (stale < 1e-20) stale = 0;
    expect(c.lookup(Elast, Pm) == stale, "compat: a later call at Emval(numEm) reuses the previous interval");
    expect(k.lookup(Elast, Pm) == (c.sfval(n, iPm) < 1e-20 ? 0. : c.sfval(n, iPm)),
           "corrected: no state is carried between calls");
    // the reference table (gen_src 6 rows) is reproduced in compat mode by l1_sf_lookup
    return 0;
}

/// A D(e,e'p)/A(e,e'p) "process" (set-up, target, radiator, reaction) for one P8.F1 configuration input.
struct EepProc {
    SimcParameters p;
    RunSetup s;
    material::HallCCryoTarget target;
    radiation::SimcRadiator rad;
    std::unique_ptr<IPhysicsModel> model;
    std::unique_ptr<reactions::NuclearEep> reaction;
};
std::unique_ptr<EepProc> eepProc(const std::string& inp, bool compat) {
    auto P = std::make_unique<EepProc>();
    InpReader rd(P->p);
    if (!rd.loadRun(inp)) return nullptr;
    SetupOptions so;
    so.fortran_compat = compat;
    so.data_dir = spectralDir();
    P->s = BuildRunSetup(P->p, so);
    P->target.setLog(nullptr);
    P->target.configure(Json{{"Mh", P->s.Mh}}, P->p);
    P->target.setLog(nullptr);
    P->target.init(P->s);
    P->rad.setLog(nullptr);
    P->rad.init(P->s);
    StageContext sc(P->s);
    sc.target = &P->target;
    sc.radiator = &P->rad;
    P->model = Registry<IPhysicsModel>::create("deforest");
    P->model->configure(Json{{"Mh", P->s.Mh}}, P->p);
    P->model->init(sc);
    sc.physics = P->model.get();
    P->reaction = std::make_unique<reactions::NuclearEep>();
    P->reaction->configure(Json{{"verbose", false}}, P->p);
    P->reaction->setLog(nullptr);
    P->reaction->init(sc);
    return P;
}

/// complete_ev of one complete_ev_eepA input row in one mode.
bool completeEvRow(EepProc& P, const simcl1::Row& r, Event& v, MainEvent& m) {
    v = Event{};
    m = MainEvent{};
    v.Ein = r["Ein"];
    v.e.xptar = r["e_xptar"];
    v.e.yptar = r["e_yptar"];
    v.e.theta = r["e_theta"];
    v.e.phi = r["e_phi"];
    v.e.E = r["e_E"];
    v.e.P = r["e_P"];
    v.e.delta = r["e_delta"];
    v.p.xptar = r["p_xptar"];
    v.p.yptar = r["p_yptar"];
    v.p.theta = r["p_theta"];
    v.p.phi = r["p_phi"];
    v.p.E = r["p_E"];
    v.p.P = r["p_P"];
    v.p.delta = r["p_delta"];
    v.Pm = r["in_v_Pm"];
    if (std::lround(r["calc_angles"]) == 1) {
        PhysicsAngles(P.s.spec.e.theta, P.s.spec.e.phi, v.e.xptar, v.e.yptar, v.e.theta, v.e.phi);
        PhysicsAngles(P.s.spec.p.theta, P.s.spec.p.phi, v.p.xptar, v.p.yptar, v.p.theta, v.p.phi);
    }
    m.target.z = r["t_z"];
    m.target.teff[0] = r["t_teff_1"];
    std::vector<double> g;
    for (int i = 1; i <= 12; ++i) g.push_back(r["g" + std::to_string(i)]);
    ScriptedRandom rng;
    rng.load(g, {});
    return P.reaction->completeEv(m, v, rng);
}

const char* kEepInp[] = {"", "p8f1_d_rad_cap", "p8f1_fe56_rad_cap", "l1_d_norad", "l1_fe56_norad", "l1_c12_rad",
                         "l1_au197_rad", "l1_he4_rad", "l1_fe56_bh_rad"};

// F5 (event.f:926, PmPer from the previous try's Pm) and F7 (event.f:944, heavy Mrec without a guard), on the rows of
// the complete_ev_eepA table: compat must equal the Fortran output, corrected must differ only as documented.
int testCompleteEvEep(bool f5) {
    const std::string dir = g_ref + "/functions/complete_ev_eepA";
    simcl1::Table in, out;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !out.load(dir + "/output.txt", err)) return 77;
    std::map<int, std::unique_ptr<EepProc>> pc, pk;
    int nCompat = 0, nCorr = 0, nRows = 0;
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        simcl1::Row r(in, i), o(out, i);
        const int cfg = static_cast<int>(std::lround(r["cfg"]));
        if (cfg < 1 || cfg > 8) continue;
        const bool fortranNaNEm = std::lround(o["success"]) == 1 && std::isnan(o["v_Em"]);
        const bool stalePm = std::lround(o["success"]) == 1 && r["in_v_Pm"] != o["v_Pm"];
        if (f5 ? !stalePm : !fortranNaNEm) continue;
        if (f5 && nRows >= 50) continue;   // enough
        ++nRows;
        const std::string inp = g_ref + "/functions/eep_state/inp/" + kEepInp[cfg] + ".inp";
        auto& C = pc[cfg];
        auto& K = pk[cfg];
        if (!C) C = eepProc(inp, true);
        if (!K) K = eepProc(inp, false);
        if (!C || !K) return 77;
        Event vc, vk;
        MainEvent mc, mk;
        const bool okc = completeEvRow(*C, r, vc, mc);
        const bool okk = completeEvRow(*K, r, vk, mk);
        if (f5) {
            const double fresh = std::sqrt(std::max(0., vk.Pm * vk.Pm - vk.PmPar * vk.PmPar - vk.PmOop * vk.PmOop));
            const double stale = std::sqrt(std::max(0., r["in_v_Pm"] * r["in_v_Pm"] - vc.PmPar * vc.PmPar - vc.PmOop * vc.PmOop));
            const double ref = o["v_PmPer"];
            if (okc && std::memcmp(&vc.PmPer, &ref, 8) == 0 && vc.PmPer == stale) ++nCompat;
            if (okk && vk.PmPer == fresh && vk.Pm == vc.Pm && vk.PmPar == vc.PmPar) ++nCorr;
        } else {
            if (okc && std::isnan(vc.Em) && std::isnan(vc.Mrec)) ++nCompat;
            if (!okk) ++nCorr;
        }
    }
    if (f5) {
        expect(nRows > 0 && nCompat == nRows, "compat: PmPer from the stale in_v_Pm, equal to the Fortran (" +
                                                  std::to_string(nCompat) + "/" + std::to_string(nRows) + " rows)");
        expect(nRows > 0 && nCorr == nRows, "corrected: PmPer from this event's Pm (" + std::to_string(nCorr) + "/" +
                                                std::to_string(nRows) + " rows)");
    } else {
        expect(nRows > 0 && nCompat == nRows, "compat: Emiss^2 < Pmiss^2 gives NaN Em/Mrec with success (" +
                                                  std::to_string(nCompat) + "/" + std::to_string(nRows) + " rows)");
        expect(nRows > 0 && nCorr == nRows, "corrected: the try fails (" + std::to_string(nCorr) + "/" +
                                                std::to_string(nRows) + " rows)");
    }
    return 0;
}
int testPmPer() { return testCompleteEvEep(true); }
int testMrecGuard() { return testCompleteEvEep(false); }

// ------------------------------------------------------------------------------------------------ D-037 (P6.4)
// F1 for the SHMS: mc_shms_coll.f:79/86/93 count slit_* per 1/20 step (Fortran) vs once per ray (corrected), as D-021.
int testShmsCollSlitCounters() {
    const double mpi2 = constants::Mpi * constants::Mpi;
    long counts[2][3] = {{0}};
    long ok[2] = {0, 0}, maxPerRay[2] = {0, 0};
    std::vector<int> stops[2];
    for (int compat = 0; compat < 2; ++compat) {
        shms::McShms sh;
        sh.setFortranCompat(compat == 1);
        std::ostringstream log;
        sh.setLog(&log);
        Mt64Random rng(2026);
        Mt64Random geo(11);
        for (long i = 0; i < 2000; ++i) {
            shms::McShmsArgs a;
            a.p_spec = 2500.;
            a.dxdz = (geo.uniform() - 0.5) * 0.10;
            a.dydz = (geo.uniform() - 0.5) * 0.08;
            a.m2 = mpi2;
            a.using_SHMScoll = true;
            a.ms_flag = true;
            const shms::ShmsCounters before = sh.counters();
            const shms::McShmsDetail d = sh.mcShms(a, rng);
            stops[compat].push_back(d.stop_code);
            if (a.ok_spec) ++ok[compat];
            for (int k = 0; k < 3; ++k)
                maxPerRay[compat] = std::max(maxPerRay[compat], sh.counters()[shms::kSlitHor + k] - before[shms::kSlitHor + k]);
        }
        for (int k = 0; k < 3; ++k) counts[compat][k] = sh.counters()[shms::kSlitHor + k];
    }
    std::printf("  SHMS slit hor/vert/oct: compat %ld %ld %ld, corrected %ld %ld %ld; successes %ld / %ld\n",
                counts[1][0], counts[1][1], counts[1][2], counts[0][0], counts[0][1], counts[0][2], ok[1], ok[0]);
    expect(stops[0] == stops[1] && ok[0] == ok[1], "the transport (stop reasons, successes, RNG use) is identical");
    expect(maxPerRay[1] > 1, "compat: some ray is counted more than once per counter (max " +
                                 std::to_string(maxPerRay[1]) + ", as mc_shms_coll.f)");
    expect(maxPerRay[0] == 1, "corrected: at most one count per ray and counter");
    return 0;
}

// F2: hSTOP_coll / sSTOP_coll are written nowhere (Fortran) vs a .geni line "Collimator stops (pi/mu)" (corrected).
int testCollStopsGeni() {
    const std::string inp = g_ref + "/runs/eep_h_elastic_norad_A/eep_h_elastic_norad_A.inp";
    if (!std::ifstream(inp)) {
        std::printf("SKIP: %s missing\n", inp.c_str());
        return 77;
    }
    SimcParameters p;
    InpReader rd(p);
    rd.loadRun(inp);   // SHMS electron arm, HMS hadron arm
    driver::RunReport r;
    for (const char* spec : {"hms", "shms"}) {
        StopCounters sc;
        sc.spectrometer = spec;
        sc.role = std::string(spec) == "hms" ? "hadron" : "electron";
        const std::string sp(spec);
        const std::vector<std::string>& keys = StopCounterKeys(sp);
        for (const std::string& k : keys) sc.counts[k] = 0;
        sc.counts["trials"] = 3932;
        sc.counts["coll"] = std::string(spec) == "hms" ? 1253 : 258;
        r.stops.push_back(sc);
    }
    for (int on = 0; on < 2; ++on) {
        p.using_HMScoll = p.using_SHMScoll = (on == 1);
        std::string geni[2];
        for (int compat = 0; compat < 2; ++compat) {
            SetupOptions so;
            so.fortran_compat = compat == 1;
            so.data_dir = g_fortran;
            geni[compat] = driver::FormatGeni(BuildRunSetup(p, so), r);
        }
        expect(geni[1].find("Collimator") == std::string::npos,
               std::string("compat: no collimator line in the .geni (simc.f:462-532), using_*coll=") + (on ? "1" : "0"));
        const auto i = geni[0].find("Collimator stops (pi/mu)");
        const auto j = i == std::string::npos ? i : geni[0].find("Collimator stops (pi/mu)", i + 1);
        const std::string l1 = i == std::string::npos ? "" : geni[0].substr(i, geni[0].find('\n', i) - i);
        const std::string l2 = j == std::string::npos ? "" : geni[0].substr(j, geni[0].find('\n', j) - j);
        expect(l1.find("1253") != std::string::npos && l2.find("258") != std::string::npos,
               "corrected: HMS and SHMS collimator stops written: '" + l1 + "' / '" + l2 + "'");
        // the corrected .geni parses back (summary/hist2json reader) with the coll counters
        std::vector<std::string> warn;
        const Json parsed = ParseFortranGeni(geni[0], &warn);
        const long hc = parsed.contains("hms") ? parsed["hms"].value("coll", -1L) : -1L;
        const long sc = parsed.contains("shms") ? parsed["shms"].value("coll", -1L) : -1L;
        long collWarn = 0;
        for (const std::string& w : warn) {
            std::printf("  parser warning: %s\n", w.c_str());
            if (w.find("Collimator") != std::string::npos) ++collWarn;
        }
        expect(hc == 1253 && sc == 258 && collWarn == 0, "the corrected .geni parses back (coll 1253 / 258)");
    }
    return 0;
}

// Runs one row of an mc_hms P6.F2 table through McHms with the row's scripted random numbers.
struct HmsRowRun {
    hms::McHmsArgs a;
    hms::McHmsDetail det;
    optics::CollTrace ct;
    long dcoll{0};
};
HmsRowRun runHmsRow(const simcl1::Row& r, hms::McHms& h) {
    HmsRowRun o;
    h.setCollTrace(&o.ct);
    const long coll0 = h.counters()[hms::kColl];
    ScriptedRandom rng;
    std::vector<double> g, u;
    for (int k = 1; r.has("g" + std::to_string(k)); ++k) g.push_back(r["g" + std::to_string(k)]);
    for (int k = 1; r.has("u" + std::to_string(k)); ++k) u.push_back(r["u" + std::to_string(k)]);
    rng.load(g, u);
    hms::McHmsArgs& a = o.a;
    a.p_spec = r["p_spec"];
    a.th_spec = r["th_spec"];
    a.dpp = r["dpp"];
    a.x = r["x"];
    a.y = r["y"];
    a.z = r["z"];
    a.dxdz = r["dxdz"];
    a.dydz = r["dydz"];
    a.m2 = r["m2"];
    a.fry = r["fry"];
    a.ms_flag = r["ms_flag"] != 0;
    a.wcs_flag = r["wcs_flag"] != 0;
    a.decay_flag = r["decay_flag"] != 0;
    a.using_HMScoll = r["using_coll"] != 0;
    a.ctau = r["ctau"];
    a.decdist = r["decdist_in"];
    a.Mh2_final = r["Mh2_final_in"];
    a.resmult = r["resmult_in"];
    a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = -999999.;
    o.det = h.mcHms(a, rng);
    o.dcoll = h.counters()[hms::kColl] - coll0;
    h.setCollTrace(nullptr);
    return o;
}

// F4: a particle that ranges out in the collimator continues with p ~ 5e-3 MeV/c (Fortran, enerloss_new clamps the
// loss) vs stops in the collimator (corrected). Rows of reference/functions/mc_hms_decay (e.g. row 720: a backward
// muon from a kaon decay).
int testCollRangedOut() {
    const std::string dir = g_ref + "/functions/mc_hms_decay";
    simcl1::Table in, out;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !out.load(dir + "/output.txt", err)) {
        std::printf("SKIP: %s\n", err.c_str());
        return 77;
    }
    long nRanged = 0, compatOk = 0, corrOk = 0, otherSame = 0, other = 0;
    hms::McHms hc, hd;   // one Fortran-mode and one corrected kernel, as one process each
    hc.setFortranCompat(true);
    hd.setFortranCompat(false);
    std::ostringstream log;
    hc.setLog(&log);
    hd.setLog(&log);
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        const simcl1::Row r(in, i), ref(out, i);
        const HmsRowRun c = runHmsRow(r, hc), d = runHmsRow(r, hd);
        if (!d.ct.ranged_out) {
            if (ref["decayed"] != 0) continue;   // decays also change with F7 (D037_F7_project_decay_drift)
            ++other;
            if (d.det.stop_code == c.det.stop_code && d.dcoll == c.dcoll) ++otherSame;
            continue;
        }
        ++nRanged;
        const bool cm = c.det.stop_code == std::lround(ref["stop_code"]) && sameBits(c.a.pathlen, ref["pathlen"]) &&
                        c.dcoll == 0;
        const bool co = d.det.stop_code == hms::kColl && d.dcoll == 1 && !d.a.ok_spec;
        compatOk += cm;
        corrOk += co;
        std::printf("  row %zu: Fortran/compat stop %s, pathlen %.6g cm (ref %.6g); corrected stop %s after %ld "
                    "collimator energy-loss steps\n",
                    i + 1, c.det.stop_code ? hms::counterName(c.det.stop_code) : "success", c.a.pathlen,
                    ref["pathlen"], d.det.stop_code ? hms::counterName(d.det.stop_code) : "success",
                    d.ct.n_eloss_calls);
    }
    expect(nRanged >= 1, "rows of mc_hms_decay that range out in the collimator: " + std::to_string(nRanged));
    expect(compatOk == nRanged, "compat: those rows continue as in the Fortran table (stop reason, pathlen bitwise)");
    expect(corrOk == nRanged, "corrected: they stop in the collimator (hSTOP_coll + 1)");
    expect(otherSame == other, "no other row changes its stop reason or coll count (" + std::to_string(otherSame) +
                                   "/" + std::to_string(other) + ")");
    return 0;
}

// F7: project.f:113 finishes the drift with the post-decay slopes (Fortran) vs ends the track at z_drift (corrected).
// All decaying rows of reference/functions/project_decay.
int testProjectDecayDrift() {
    const std::string dir = g_ref + "/functions/project_decay";
    simcl1::Table in, out;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !out.load(dir + "/output.txt", err)) {
        std::printf("SKIP: %s\n", err.c_str());
        return 77;
    }
    long n = 0, compatBit = 0, corrExact = 0;
    double maxDx[2] = {0, 0}, maxDpath[2] = {0, 0};   // [synthetic, capture]
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        const simcl1::Row r(in, i), ref(out, i);
        if (ref["decayed"] == 0 || std::lround(r["routine"]) != 1) continue;
        double xo[2], yo[2], po[2], zdec = 0, xdec = 0, s_after = 0, norm_old = 0;
        for (int compat = 0; compat < 2; ++compat) {
            ScriptedRandom rng;
            std::vector<double> u;
            for (int k = 1; k <= 4; ++k) u.push_back(r["u" + std::to_string(k)]);
            rng.load({}, u);
            optics::Track t;
            t.xs = r["xs"];
            t.dxdzs = r["dxdzs"];
            t.ys = r["ys"];
            t.dydzs = r["dydzs"];
            t.dpps = r["dpps"];
            optics::DecayState d;
            d.decay_flag = true;
            d.m2 = r["m2"];
            d.ph = r["ph"];
            d.ctau = r["ctau"];
            d.decdist = r["decdist_in"];
            d.Mh2_final = r["Mh2_final_in"];
            d.rng = &rng;
            d.fortran_compat = compat == 1;
            optics::TransportTrace tt;
            d.trace = &tt;
            double path = r["pathlen_in"];
            optics::project(t, r["zd"], path, &d);
            xo[compat] = t.xs;
            yo[compat] = t.ys;
            po[compat] = path;
            zdec = tt.dec.z_decay;
            xdec = tt.dec.x;
            s_after = tt.dec.dxdz_after;
            norm_old = std::sqrt(1 + r["dxdzs"] * r["dxdzs"] + r["dydzs"] * r["dydzs"]);
        }
        ++n;
        if (sameBits(xo[1], ref["xs_out"]) && sameBits(yo[1], ref["ys_out"]) && sameBits(po[1], ref["pathlen_out"]))
            ++compatBit;
        // corrected: the remaining z is z_drift - z_decay/sqrt(1+s_parent**2), so the track ends at z_drift
        if (sameBits(xo[0], xdec + s_after * (r["zd"] - zdec / norm_old))) ++corrExact;
        const int k = r["gen_src"] == 1 ? 1 : 0;
        maxDx[k] = std::max(maxDx[k], std::max(std::fabs(xo[0] - xo[1]), std::fabs(yo[0] - yo[1])));
        maxDpath[k] = std::max(maxDpath[k], std::fabs(po[0] - po[1]));
    }
    std::printf("  %ld decaying project rows; corrected - Fortran: max |dx|,|dy| %.4g cm (synthetic) %.4g cm (capture "
                "steps), max |dpathlen| %.4g / %.4g cm\n",
                n, maxDx[0], maxDx[1], maxDpath[0], maxDpath[1]);
    expect(n > 0 && compatBit == n, "compat reproduces project_decay bitwise (" + std::to_string(compatBit) + "/" +
                                        std::to_string(n) + ")");
    expect(corrExact == n, "corrected: the drift after the decay ends at z_drift (" + std::to_string(corrExact) + "/" +
                               std::to_string(n) + ")");
    expect(maxDx[0] > 0, "the two modes differ");
    return 0;
}

// F9: a decay of a particle that is neither a pion nor a kaon executes a bare STOP (Fortran) vs is never checked for
// (corrected) and the configuration doing_decay + such a hadron is rejected at start-up.
int testDecayNonPiK() {
    // kernel: a proton with decay_flag and ctau = 0 decays at once in the Fortran (z_decay = -0)
    for (int compat = 0; compat < 2; ++compat) {
        ScriptedRandom rng;
        rng.load({}, {0.3, 0.2, 0.1});
        optics::Track t;
        t.dxdzs = 0.01;
        optics::DecayState d;
        d.decay_flag = true;
        d.m2 = constants::Mp * constants::Mp;
        d.ph = 2000.;
        d.ctau = 0.;
        d.rng = &rng;
        d.fortran_compat = compat == 1;
        double path = 0;
        bool stopped = false;
        std::string what;
        try {
            optics::project(t, 100., path, &d);
        } catch (const FortranStop& e) {
            stopped = true;
            what = e.what();
        }
        if (compat)
            expect(stopped && what.find("cannot decay particle with mass=  938.27") != std::string::npos,
                   "compat: FortranStop '" + what + "' (project.f:87-91)");
        else
            expect(!stopped && rng.nUniform() == 0 && !d.dflag && sameBits(t.xs, 0.01 * 100.),
                   "corrected: no decay check (no random number), plain drift");
    }
    // set-up: doing_decay = 1 with the proton hadron of H(e,e'p)
    const std::string inp = g_ref + "/runs/eep_h_elastic_norad_A/eep_h_elastic_norad_A.inp";
    if (!std::ifstream(inp)) {
        std::printf("SKIP: %s missing\n", inp.c_str());
        return 77;
    }
    SimcParameters p;
    InpReader rd(p);
    rd.loadRun(inp);
    p.doing_decay = true;
    for (int compat = 0; compat < 2; ++compat) {
        SetupOptions so;
        so.fortran_compat = compat == 1;
        so.data_dir = g_fortran;
        bool rejected = false;
        std::string what;
        try {
            BuildRunSetup(p, so);
        } catch (const std::invalid_argument& e) {
            rejected = true;
            what = e.what();
        }
        if (compat)
            expect(!rejected, "compat: the set-up accepts doing_decay with a proton (dbase.f)");
        else
            expect(rejected && what.find("neither a pion nor a kaon") != std::string::npos,
                   "corrected: rejected at start-up: '" + what + "'");
    }
    // a pion hadron is accepted in both modes
    for (int compat = 0; compat < 2; ++compat) {
        SimcParameters q = p;
        q.doing_pion = true;
        SetupOptions so;
        so.fortran_compat = compat == 1;
        so.data_dir = g_fortran;
        bool rejected = false;
        try {
            BuildRunSetup(q, so);
        } catch (const std::invalid_argument&) {
            rejected = true;
        } catch (const std::exception&) {
        }
        expect(!rejected, std::string("doing_decay with a pion hadron is accepted, compat=") + (compat ? "1" : "0"));
    }
    return 0;
}

// ------------------------------------------------------------------------------------------------ D-038 (P9.2, pion)
std::string pionInp(const std::string& name) { return g_ref + "/functions/pion_state/inp/" + name + ".inp"; }
bool loadPionInp(const std::string& name, SimcParameters& p) {
    const std::string inp = pionInp(name);
    if (!std::ifstream(inp)) {
        std::printf("SKIP: %s missing\n", inp.c_str());
        return false;
    }
    InpReader rd(p);
    return rd.loadRun(inp);
}
std::string specDir() { return std::string(SIMC_SOURCE_DIR) + "/data/spectral"; }

// D-038 F1: dbase.f:823-826 stops charged H(e,e'pi-)Delta++ (which_pion = 3 without doing_pizero) too.
int testDeltaPPStop() {
    SimcParameters p, q;
    if (!loadPionInp("l1s_h_delta3", p) || !loadPionInp("l1s_h_pi0delta3", q)) return 77;
    for (int compat = 0; compat < 2; ++compat) {
        SetupOptions so;
        so.fortran_compat = compat == 1;
        so.data_dir = specDir();
        const RunSetup s = BuildRunSetup(p, so);
        const RunSetup s0 = BuildRunSetup(q, so);
        const bool stopped = !s.ok() && s.stops.front() == "STOP Can't do Delta+ production from neutron for pi0";
        const bool stopped0 = !s0.ok() && s0.stops.front() == "STOP Can't do Delta+ production from neutron for pi0";
        if (compat)
            expect(stopped, "compat: charged H(e,e'pi-)Delta++ (which_pion = 3) STOPs as the Fortran (dbase.f:823-826)");
        else
            expect(!stopped && s.ok() && s.targ.Mtar_struck == constants::Mp && s.targ.Mrec_struck == constants::MDelta,
                   "corrected: charged H(e,e'pi-)Delta++ is set up (Mtar_struck = Mp, Mrec_struck = MDelta)");
        expect(stopped0, std::string("pi0 Delta0 from the neutron on H STOPs in both modes, compat=") +
                             (compat ? "1" : "0"));
    }
    // corrected: the configuration runs (50 events) and the weight carries the 0.55 factor of event.f:1479-1484
    driver::RunOptions opt;
    opt.quiet_stages = true;
    opt.override_ngen = true;
    opt.ngen = 50;
    opt.run_name = "d038_f1";
    opt.data_dir = specDir();
    std::ostringstream log;
    opt.log = &log;
    driver::SimcRun run(p, Json{{"sinks", Json::array()}}, opt);
    bool ran = false;
    std::string what;
    try {
        run.run();
        ran = run.totals().nevent == 50 && run.totals().wtcontribute > 0;
    } catch (const std::exception& e) {
        what = e.what();
    }
    expect(ran, "corrected: H(e,e'pi-)Delta++ runs 50 events" + (what.empty() ? "" : " (" + what + ")"));
    return 0;
}

// D-038 F2: pi0 'from the neutron' (which_pion = 1, doing_pizero) accepted on hydrogen by dbase.f:141-142.
int testPi0NeutronOnH() {
    SimcParameters p;
    if (!loadPionInp("l1_h_k2_pi0n", p)) return 77;
    for (int compat = 0; compat < 2; ++compat) {
        SetupOptions so;
        so.fortran_compat = compat == 1;
        so.data_dir = specDir();
        bool rejected = false;
        std::string what;
        RunSetup s;
        try {
            s = BuildRunSetup(p, so);
        } catch (const std::invalid_argument& e) {
            rejected = true;
            what = e.what();
        }
        if (compat)
            expect(!rejected && s.ok() && s.targ.Mtar_struck == constants::Mn,
                   "compat: accepted, Mtar_struck = Mn on a proton target (dbase.f:141-142, 333)");
        else
            expect(rejected && what.find("hydrogen has no neutron") != std::string::npos,
                   "corrected: rejected at start-up: '" + what + "'");
    }
    return 0;
}

// D-038 F5: CENTRAL.sigcc of pion runs computed with main0%epsilon = t = phi_pq = 0 (simc.f:1236-1290).
int testCentralSigcc() {
    SimcParameters p;
    if (!loadPionInp("p9f1_pip_h_dec1_cap", p)) return 77;
    double sig[2] = {0, 0};
    MainEvent m0[2];
    Event v0[2];
    for (int compat = 0; compat < 2; ++compat) {
        driver::RunOptions opt;
        opt.quiet_stages = true;
        opt.fortran_compat = compat == 1;
        opt.data_dir = specDir();
        std::ostringstream log;
        opt.log = &log;
        driver::SimcRun run(p, Json{{"sinks", Json::array()}}, opt);
        run.initialise();
        sig[compat] = run.central().central.sigcc;
        m0[compat] = run.central().main0;
        v0[compat] = run.central().vertex0;
    }
    // the Fortran value (reference/functions/pion_state cfg 1 c_sigcc; .hist CENTRAL.sigcc 0.839106E-06)
    expect(sameBits(sig[1], 8.3910606773484538e-07) || std::fabs(sig[1] / 0.839106e-6 - 1) < 1e-6,
           "compat: CENTRAL.sigcc = " + fortran::E(sig[1], 13, 6) + " (Fortran 0.839106E-06)");
    expect(m0[1].epsilon == 0 && m0[1].t == 0 && m0[1].phi_pq == 0, "compat: main0 epsilon, t, phi_pq stay 0");
    const Event& v = v0[0];
    const double tt = v.Q2 - constants::Mpi * constants::Mpi + 2 * v.nu * v.p.E -
                      2 * v.p.P * v.q * std::cos(m0[0].theta_pq);
    expect(sameBits(m0[0].epsilon, 1. / (1. + 2. * (1 + v.nu * v.nu / v.Q2) * (std::tan(v.e.theta / 2.) * std::tan(v.e.theta / 2.)))) ||
               std::fabs(m0[0].epsilon - v.epsilon) < 1e-12,
           "corrected: main0 epsilon = the central epsilon (" + std::to_string(m0[0].epsilon) + ")");
    expect(std::fabs(m0[0].t - tt) <= 1e-9 * std::fabs(tt) && m0[0].t > 0,   // SIMC t = -t_Mandelstam
           "corrected: main0 t = the central t (" + std::to_string(m0[0].t / 1e6) + " GeV^2)");
    char b[160];
    std::snprintf(b, sizeof b, "corrected CENTRAL.sigcc %.6E differs from the Fortran %.6E (ratio %.4f)", sig[0], sig[1],
                  sig[0] / sig[1]);
    expect(sig[0] != sig[1] && sig[0] > 0, b);
    return 0;
}

// D-043 F13 (= D-038 F5 for kaons): calculate_central evaluates peeK with main0 epsilon = t = theta_pq = phi_pq = 0
// (simc.f:1236-1290). Kaon capture input p10f1_kl_h_dec1_cap (K+ Lambda H, K1).
int testCentralSigccKaon() {
    SimcParameters p;
    const std::string inp = g_ref + "/functions/kaon_state/inp/p10f1_kl_h_dec1_cap.inp";
    if (!std::ifstream(inp)) {
        std::printf("SKIP: %s missing\n", inp.c_str());
        return 77;
    }
    InpReader rd(p);
    if (!rd.loadRun(inp)) return 1;
    double sig[2] = {0, 0};
    MainEvent m0[2];
    for (int compat = 0; compat < 2; ++compat) {
        driver::RunOptions opt;
        opt.quiet_stages = true;
        opt.fortran_compat = compat == 1;
        opt.data_dir = specDir();
        std::ostringstream log;
        opt.log = &log;
        driver::SimcRun run(p, Json{{"sinks", Json::array()}}, opt);
        run.initialise();
        sig[compat] = run.central().central.sigcc;
        m0[compat] = run.central().main0;
    }
    // the Fortran value (reference/functions/kaon_state cfg 1 c_sigcc; .hist CENTRAL.sigcc 0.143536E-06)
    expect(std::fabs(sig[1] / 0.143536e-6 - 1) < 5e-6, "compat: CENTRAL.sigcc = " + fortran::E(sig[1], 13, 6) +
                                                            " (Fortran 0.143536E-06)");
    expect(m0[1].epsilon == 0 && m0[1].t == 0 && m0[1].theta_pq == 0 && m0[1].phi_pq == 0,
           "compat: main0 epsilon, t, theta_pq, phi_pq stay 0");
    expect(m0[0].epsilon > 0.8 && m0[0].epsilon < 1 && m0[0].t > 0, "corrected: main0 epsilon, t = the central values");
    char b[160];
    std::snprintf(b, sizeof b, "corrected CENTRAL.sigcc %.6E differs from the Fortran %.6E (ratio %.4f)", sig[0], sig[1],
                  sig[0] / sig[1]);
    expect(sig[0] != sig[1] && sig[0] > 0, b);
    return 0;
}

// D-038 F15: the .hist final-state line is missing for which_pion = 11 (simc.f:750-755).
int testHistWhichPion11() {
    SimcParameters p;
    if (!loadPionInp("l1_c12coh_k2_pim", p)) return 77;
    driver::RunReport r;
    for (int compat = 0; compat < 2; ++compat) {
        SetupOptions so;
        so.fortran_compat = compat == 1;
        so.data_dir = specDir();
        const RunSetup s = BuildRunSetup(p, so);
        const std::string h = driver::FormatHist(s, r);
        const bool has = h.find("Final State is A + pi") != std::string::npos;
        if (compat)
            expect(!has, "compat: no final-state line for which_pion = 11 (simc.f:750-755)");
        else
            expect(has, "corrected: '****----  Final State is A + pi ----****' for which_pion = 11");
    }
    return 0;
}


// ------------------------------------------------------------------------------------------------ D-048 (P11.1, Delta)
std::string deltaInp(const std::string& name) {
    const std::string a = g_ref + "/functions/delta_state/inp/" + name + ".inp";
    if (std::ifstream(a)) return a;
    return g_ref + "/runs/" + name + "/" + name + ".inp";
}
bool loadDeltaInp(const std::string& name, SimcParameters& p) {
    const std::string inp = deltaInp(name);
    if (!std::ifstream(inp)) {
        std::printf("SKIP: %s missing\n", inp.c_str());
        return false;
    }
    InpReader rd(p);
    return rd.loadRun(inp);
}

/// Output sink that checks the complete_ev energy-momentum balance of every contributing event at the vertex:
/// (nu + efer - E_p)^2 - |q + |pfer| pfer_hat - p|^2 = Mrec_struck^2 (event.f:633-679, both roots), and collects
/// plausibility figures (weights, root, Fermi momentum, the radiating tail).
struct DeltaCheckSink : IOutputSink {
    double mrec2{0};
    long n{0}, bad{0}, nonfinite{0}, tail1{0}, pferNonZero{0};
    double maxrel{0}, sumw{0}, wmin{1e300}, wmax{0}, pmin{1e300}, pmax{0};
    std::string name() const override { return "delta_check"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext&) override {}
    void beginRun(const RunSetup& s) override { mrec2 = s.targ.Mrec_struck * s.targ.Mrec_struck; }
    void fill(const MainEvent& m, const Event& v, const Event&, const Event&, bool success) override {
        if (!success) return;
        ++n;
        const FermiMomentum& f = m.pfermi;
        const double ap = std::fabs(f.pfer);
        const double rx = v.q * v.uq.x + ap * f.pferx - v.p.P * v.up.x;
        const double ry = v.q * v.uq.y + ap * f.pfery - v.p.P * v.up.y;
        const double rz = v.q * v.uq.z + ap * f.pferz - v.p.P * v.up.z;
        const double er = v.nu + f.efer - v.p.E;
        const double m2 = er * er - (rx * rx + ry * ry + rz * rz);
        const double rel = std::fabs(m2 - mrec2) / (er * er);   // relative to the terms that cancel
        maxrel = std::max(maxrel, rel);
        if (!(rel < 1e-9)) ++bad;
        if (!std::isfinite(m.weight) || !(m.weight > 0)) ++nonfinite;
        sumw += m.weight;
        wmin = std::min(wmin, m.weight);
        wmax = std::max(wmax, m.weight);
        pmin = std::min(pmin, v.p.P);
        pmax = std::max(pmax, v.p.P);
        if (m.ntup.radarm == 1) ++tail1;
        if (f.pfer != 0) ++pferNonZero;
    }
    void endRun(const RunSummary&) override {}
};

struct DeltaRunResult {
    bool ran{false};
    std::string what;
    long ntried{0}, nevent{0};
    double wtcontr{0};
    DeltaCheckSink* sink{nullptr};
    std::vector<std::string> messages;
};
DeltaRunResult runDelta(const SimcParameters& p, bool compat, long ngen, std::unique_ptr<driver::SimcRun>& keep) {
    DeltaRunResult r;
    driver::RunOptions opt;
    opt.quiet_stages = true;
    opt.fortran_compat = compat;
    opt.override_ngen = true;
    opt.ngen = ngen;
    opt.override_seed = true;
    opt.seed = 4242;
    opt.run_name = "d048";
    opt.data_dir = specDir();
    static std::ostringstream log;
    opt.log = &log;
    keep = std::make_unique<driver::SimcRun>(p, Json{{"sinks", Json::array()}}, opt);
    auto sk = std::make_unique<DeltaCheckSink>();
    r.sink = sk.get();
    keep->addSink(std::move(sk));
    try {
        keep->run();
        r.ran = true;
        r.ntried = keep->totals().ntried;
        r.nevent = keep->totals().nevent;
        r.wtcontr = keep->totals().wtcontribute;
        r.messages = keep->setup().messages;
    } catch (const std::exception& e) {
        r.what = e.what();
    }
    return r;
}

// D-048 F1 (D-018): dbase.f:779-788 tests the pion flags doing_deutpi/doing_hepi, so every D/A (e,e'p)pi run STOPs;
// corrected: D/A Delta runs (the D/A branches of event.f/init.f, never run in the Fortran: self-consistency only).
int testDeltaDAStop() {
    const char* names[] = {"l1s_d_d1", "l1s_he3_d1", "l1s_c12_d1"};
    for (const char* nm : names) {
        SimcParameters p;
        if (!loadDeltaInp(nm, p)) return 77;
        for (int compat = 0; compat < 2; ++compat) {
            SetupOptions so;
            so.fortran_compat = compat == 1;
            so.data_dir = specDir();
            const RunSetup s = BuildRunSetup(p, so);
            const std::string msg = "STOP I don't have ANY idea what (e,e'p)pi we're doing!!!";
            if (compat) {
                expect(!s.ok() && s.stops.front() == msg, std::string(nm) + " compat: STOPs as the Fortran: " +
                                                              (s.ok() ? "<no STOP>" : s.stops.front()));
            } else {
                const bool d = std::lround(p.targ.A) == 2;
                expect(s.ok() && (d ? s.flags.doing_deutdelta : s.flags.doing_hedelta) && !s.pdist.pval.empty() &&
                           s.VERTEXedge.Pm.max == s.pdist.pval.back() && (d || s.VERTEXedge.Em.max > 0),
                       std::string(nm) + " corrected: set up (" + (d ? "doing_deutdelta" : "doing_hedelta") +
                           ", momentum distribution " + s.pdist.file + ", nump " + std::to_string(s.pdist.pval.size()) +
                           (d ? "" : ", spectral function " + s.sf.file + ", Em max " + std::to_string(s.VERTEXedge.Em.max)) +
                           ")");
            }
        }
        // corrected: 300 events; energy-momentum balance of every event at the vertex, finite positive weights
        std::unique_ptr<driver::SimcRun> keep;
        const DeltaRunResult r = runDelta(p, false, 300, keep);
        expect(r.ran && r.nevent == 300, std::string(nm) + " corrected: runs 300 events" + (r.what.empty() ? "" : " (" + r.what + ")"));
        if (!r.ran) continue;
        const DeltaCheckSink& k = *r.sink;
        char b[400];
        std::snprintf(b, sizeof b,
                      "%s corrected: %ld events, Ntried %ld, (nu+efer-Ep)^2-|q+pfer-p|^2 = Mpi^2 for all (max rel %.2e, %ld bad), "
                      "weights finite > 0 (%ld bad; %.3g..%.3g), Fermi momentum non-zero in %ld, p %.0f..%.0f MeV/c, "
                      "tail-1 events %ld",
                      nm, k.n, r.ntried, k.maxrel, k.bad, k.nonfinite, k.wmin, k.wmax, k.pferNonZero, k.pmin, k.pmax, k.tail1);
        expect(k.n == 300 && k.bad == 0 && k.nonfinite == 0 && k.pferNonZero == k.n, b);
    }
    return 0;
}

// D-048 F2: radc.f:288-293 sets no tail-1 Egamma limits for doing_delta: every ntail = 1 try fails.
int testDeltaTail1() {
    SimcParameters p;
    if (!loadDeltaInp("p11f1_dh_d1_rad_cap", p)) return 77;
    DeltaRunResult r[2];
    std::unique_ptr<driver::SimcRun> keep[2];
    for (int compat = 0; compat < 2; ++compat) r[compat] = runDelta(p, compat == 1, 400, keep[compat]);
    if (!r[0].ran || !r[1].ran) {
        expect(false, "runs: " + r[0].what + " / " + r[1].what);
        return 0;
    }
    char b[300];
    std::snprintf(b, sizeof b, "compat: 0 of %ld events radiate in tail 1 (radarm = 1), Ntried %ld", r[1].sink->n, r[1].ntried);
    expect(r[1].sink->tail1 == 0, b);
    std::snprintf(b, sizeof b, "corrected: %ld of %ld events radiate in tail 1, Ntried %ld (compat %ld, ratio %.3f)",
                  r[0].sink->tail1, r[0].sink->n, r[0].ntried, r[1].ntried,
                  static_cast<double>(r[0].ntried) / static_cast<double>(r[1].ntried));
    expect(r[0].sink->tail1 > 0 && r[0].ntried < r[1].ntried, b);
    expect(r[0].sink->bad == 0 && r[1].sink->bad == 0, "both modes: vertex energy-momentum balance holds for every event");
    return 0;
}

// D-048 F3: physics_delta.f:99-100 passes binary128 quotients to the REAL*8 dummies of sig_blok.
int testDeltaSigblok() {
    simcl1::Table in, out;
    std::string err;
    const std::string dir = g_ref + "/functions/peedelta";
    if (!in.load(dir + "/input.txt", err) || !out.load(dir + "/output.txt", err)) {
        std::printf("SKIP: %s\n", err.c_str());
        return 77;
    }
    long nrow = 0, compatBad = 0, corrBad = 0, valBad = 0, nonfiniteCompat = 0, nonfiniteCorr = 0, diffModes = 0;
    for (std::size_t i = 0; i < out.rows.size(); ++i) {
        simcl1::Row r(out, i);
        Event v;
        v.Ein = r["v_Ein"];
        v.nu = r["v_nu"];
        v.q = r["v_q"];
        v.Q2 = r["v_Q2"];
        v.e.E = r["v_e_E"];
        v.p.P = r["v_p_P"];
        v.p.E = r["v_p_E"];
        v.up = {r["v_up_x"], r["v_up_y"], r["v_up_z"]};
        v.uq = {r["v_uq_x"], r["v_uq_y"], r["v_uq_z"]};
        MainEvent m;
        m.epsilon = r["in_m_epsilon"];
        m.t = r["in_m_t"];
        m.phi_pq = r["in_m_phi_pq"];
        m.pfermi = {r["pfer"], r["pferx"], r["pfery"], r["pferz"], r["efer"]};
        delta::DeltaConfig c;
        c.which_pion = static_cast<int>(std::lround(r["eff_which_pion"]));
        c.Mtar_struck = r["eff_Mtar_struck"];
        c.sigblok_binary128_quirk = true;
        MainEvent m1 = m, m2 = m;
        const double v1 = delta::peedelta(v, m1, c);
        c.sigblok_binary128_quirk = false;
        const double v2 = delta::peedelta(v, m2, c);
        ++nrow;
        auto same = [](double a, double b) { return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, 8) == 0; };
        if (!same(m1.ntup.sigcm1, r["ntup_sigcm1"])) ++compatBad;
        if (!same(m2.ntup.sigcm1, r["hx_sigblok_intended"])) ++corrBad;
        if (!same(v1, r["peedelta"]) || !same(v2, v1)) ++valBad;
        if (!std::isfinite(m1.ntup.sigcm1)) ++nonfiniteCompat;
        if (!std::isfinite(m2.ntup.sigcm1)) ++nonfiniteCorr;
        if (!same(m1.ntup.sigcm1, m2.ntup.sigcm1)) ++diffModes;
    }
    char b[300];
    std::snprintf(b, sizeof b, "compat: ntup%%sigcm = the Fortran (low 8 bytes of the binary128 arguments) on %ld of %ld rows "
                  "(%ld non-finite)", nrow - compatBad, nrow, nonfiniteCompat);
    expect(compatBad == 0 && nrow > 3000, b);
    std::snprintf(b, sizeof b, "corrected: ntup%%sigcm = sig_blok(t, Q2, s in GeV^2) = hx_sigblok_intended on %ld of %ld rows "
                  "(%ld non-finite; %ld rows differ from compat)", nrow - corrBad, nrow, nonfiniteCorr, diffModes);
    expect(corrBad == 0 && nonfiniteCorr == 0, b);
    std::snprintf(b, sizeof b, "both modes: the peedelta value (the weight) is the Fortran one on all rows (%ld bad)", valBad);
    expect(valBad == 0, b);
    return 0;
}

// D-048 F6: dbase.f:184-186 quote typo in the A >= 2 Delta warning.
int testDeltaWarningText() {
    SimcParameters p;
    if (!loadDeltaInp("l1s_d_d1", p)) return 77;
    std::string fortranLine;   // the verbatim line of the Fortran run (reference/runs/delta_d_d1_rad_A/*.out)
    {
        std::ifstream f(g_ref + "/runs/delta_d_d1_rad_A/delta_d_d1_rad_A.out");
        for (std::string l; std::getline(f, l);)
            if (l.find("WARNING: Delta cross section model") != std::string::npos) fortranLine = l;
    }
    for (int compat = 0; compat < 2; ++compat) {
        SetupOptions so;
        so.fortran_compat = compat == 1;
        so.data_dir = specDir();
        const RunSetup s = BuildRunSetup(p, so);
        std::string got;
        for (const auto& m : s.messages)
            if (m.find("WARNING: Delta cross section model") != std::string::npos) got = m;
        if (compat)
            expect(!fortranLine.empty() && " " + got == fortranLine,
                   "compat: the Fortran text (list-directed, leading blank) '" + got.substr(0, 40) + "...'");
        else
            expect(got == "WARNING: Delta cross section model only set up for proton target!", "corrected: '" + got + "'");
    }
    return 0;
}

// D-048 F9 (= D-038 F5): the central peedelta sees main0 epsilon = t = phi_pq = 0 (simc.f:1236-1290).
int testCentralSigccDelta() {
    SimcParameters p;
    if (!loadDeltaInp("p11f1_dh_d1_rad_cap", p)) return 77;
    double sig[2] = {0, 0};
    MainEvent m0[2];
    for (int compat = 0; compat < 2; ++compat) {
        driver::RunOptions opt;
        opt.quiet_stages = true;
        opt.fortran_compat = compat == 1;
        opt.data_dir = specDir();
        std::ostringstream log;
        opt.log = &log;
        driver::SimcRun run(p, Json{{"sinks", Json::array()}}, opt);
        run.initialise();
        sig[compat] = run.central().central.sigcc;
        m0[compat] = run.central().main0;
    }
    // the Fortran value (.hist of delta_h_d1_rad_A / p11f1_dh_d1_rad_cap: CENTRAL.sigcc = 0.227051E+01)
    expect(std::fabs(sig[1] / 2.27051 - 1) < 5e-6, "compat: CENTRAL.sigcc = " + fortran::E(sig[1], 13, 6) +
                                                       " (Fortran 0.227051E+01)");
    expect(m0[1].epsilon == 0 && m0[1].t == 0 && m0[1].phi_pq == 0, "compat: main0 epsilon, t, phi_pq stay 0");
    expect(m0[0].epsilon > 0.8 && m0[0].epsilon < 1, "corrected: main0 epsilon = the central epsilon");
    char b[160];
    std::snprintf(b, sizeof b, "corrected CENTRAL.sigcc %.6E vs Fortran %.6E (ratio %.4f)", sig[0], sig[1], sig[0] / sig[1]);
    expect(sig[0] != sig[1] && sig[0] > 0, b);
    return 0;
}


// ------------------------------------------------------------------------------------------------ D-051 (P12.2, SIDIS)
bool loadSidisInp(const std::string& name, SimcParameters& p) {
    std::string inp = g_ref + "/functions/sidis/sidis_state/inp/" + name + ".inp";
    if (!std::ifstream(inp)) inp = g_ref + "/runs/" + name + "/" + name + ".inp";
    if (!std::ifstream(inp)) {
        std::printf("SKIP: %s missing\n", inp.c_str());
        return false;
    }
    InpReader rd(p);
    return rd.loadRun(inp);
}

/// Sink: the SIDIS ntuple of every contributing event (FillSemiNtuple of the run's mode) and the vertex PmPer test.
struct SemiSink : IOutputSink {
    const RunSetup* s{nullptr};
    long n{0}, pfermiNan{0}, pfermiZero{0}, pmperCurrent{0};
    std::string name() const override { return "semi_check"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext&) override {}
    void beginRun(const RunSetup& rs) override { s = &rs; }
    void fill(const MainEvent& m, const Event& v, const Event& o, const Event& r, bool success) override {
        if (!success) return;
        ++n;
        std::vector<double> ntu(kSemiNtupleSize, 0.0);
        FillSemiNtuple(*s, m, v, o, r, ntu);
        if (std::isnan(ntu[53])) ++pfermiNan;
        if (ntu[53] == 0.0) ++pfermiZero;
        const double pp = std::sqrt(std::max(0.e0, v.Pm * v.Pm - v.PmPar * v.PmPar - v.PmOop * v.PmOop));
        if (sameBits(pp, v.PmPer)) ++pmperCurrent;
    }
    void endRun(const RunSummary&) override {}
};

struct SemiRun {
    bool ran{false};
    std::string what;
    SemiSink* sink{nullptr};
};
SemiRun runSemi(const SimcParameters& p, bool compat, long ngen, std::unique_ptr<driver::SimcRun>& keep) {
    SemiRun r;
    driver::RunOptions opt;
    opt.quiet_stages = true;
    opt.fortran_compat = compat;
    opt.override_ngen = true;
    opt.ngen = ngen;
    opt.override_seed = true;
    opt.seed = 4242;
    opt.run_name = "d051";
    opt.data_dir = specDir();
    static std::ostringstream log;
    opt.log = &log;
    keep = std::make_unique<driver::SimcRun>(p, Json{{"sinks", Json::array()}}, opt);
    auto sk = std::make_unique<SemiSink>();
    r.sink = sk.get();
    keep->addSink(std::move(sk));
    try {
        keep->run();
        r.ran = true;
    } catch (const std::exception& e) {
        r.what = e.what();
    }
    return r;
}

// D-051 F3 (results_write.f:208-209, pfermi 0/0 for pfer = 0) and F9 (event.f:922, PmPer from the previous Pm), H pi+.
int testSemiPfermiPmPer() {
    SimcParameters p;
    if (!loadSidisInp("l1_h_k1_pip_norad", p)) return 77;
    SemiRun r[2];
    std::unique_ptr<driver::SimcRun> keep[2];
    for (int compat = 0; compat < 2; ++compat) r[compat] = runSemi(p, compat == 1, 300, keep[compat]);
    if (!r[0].ran || !r[1].ran) {
        expect(false, "runs: " + r[0].what + " / " + r[1].what);
        return 0;
    }
    const SemiSink &c = *r[1].sink, &k = *r[0].sink;
    char b[300];
    std::snprintf(b, sizeof b, "F3 compat: pfermi = 0/0 = NaN in %ld of %ld H events (as the Fortran .bin)", c.pfermiNan, c.n);
    expect(c.n == 300 && c.pfermiNan == c.n, b);
    std::snprintf(b, sizeof b, "F3 corrected: pfermi = 0 in %ld of %ld H events", k.pfermiZero, k.n);
    expect(k.n == 300 && k.pfermiZero == k.n, b);
    std::snprintf(b, sizeof b, "F9 compat: vertex PmPer from this event's Pm in %ld of %ld events (stale Pm otherwise)",
                  c.pmperCurrent, c.n);
    expect(c.pmperCurrent < c.n / 10, b);
    std::snprintf(b, sizeof b, "F9 corrected: vertex PmPer from this event's Pm in %ld of %ld events", k.pmperCurrent, k.n);
    expect(k.pmperCurrent == k.n, b);
    return 0;
}

// D-051 F12 (simc.f:679-700: no pi0 title) and F14 (simc.f:701-717: A >= 3 kaon SIDIS STOPs at the report).
int testSemiHistTitles() {
    {
        SimcParameters p;
        if (!loadSidisInp("l1_h_k1_pip_norad", p)) return 77;
        p.doing_pizero = true;
        p.pizero_ngamma = 2;
        driver::RunReport rr;
        for (int compat = 0; compat < 2; ++compat) {
            SetupOptions so;
            so.fortran_compat = compat == 1;
            so.data_dir = specDir();
            const RunSetup s = BuildRunSetup(p, so);
            const std::string h = driver::FormatHist(s, rr);
            if (compat)
                expect(h.find("H(e,e'pi+)X") != std::string::npos && h.find("pi0)X") == std::string::npos,
                       "F12 compat: a pi0 SIDIS run is titled H(e,e'pi+)X (doing_hplus = 1)");
            else
                expect(h.find("H(e,e'pi0)X") != std::string::npos, "F12 corrected: title H(e,e'pi0)X");
        }
    }
    {
        SimcParameters p;
        if (!loadSidisInp("l1_c12_k1_kp_rad", p)) return 77;
        driver::RunReport rr;
        for (int compat = 0; compat < 2; ++compat) {
            SetupOptions so;
            so.fortran_compat = compat == 1;
            so.data_dir = specDir();
            const RunSetup s = BuildRunSetup(p, so);
            std::string h, err;
            try {
                h = driver::FormatHist(s, rr);
            } catch (const FortranStop& e) {
                err = e.what();
            }
            if (compat)
                expect(err.find("I don't have ANY idea what A(e,e'k)X we're doing!!!") != std::string::npos,
                       "F14 compat: C-12 K+ SIDIS STOPs when the .hist is written: " + err);
            else
                expect(err.empty() && h.find("A(e,e'k+)X") != std::string::npos,
                       "F14 corrected: the .hist is written with A(e,e'k+)X" + (err.empty() ? "" : " (" + err + ")"));
        }
        // corrected: the C-12 K+ run itself (100 events, finite positive weights)
        std::unique_ptr<driver::SimcRun> keep;
        const SemiRun r = runSemi(p, false, 100, keep);
        expect(r.ran && r.sink->n == 100, "F14 corrected: C-12 K+ SIDIS runs 100 events" + (r.what.empty() ? "" : " (" + r.what + ")"));
    }
    return 0;
}


// ------------------------------------------------------------------------------------------------ D-061 (P14.2)
// F1: mc_hrs*.f pass REAL(16) transport lengths to transp's real*8 zd: FP path = -3.134165696542769E+232 and the HRS kaon
// survival weight +Inf (reference/runs/p14f1_kaon_hrs_probe: Weight = +Inf for every event). Corrected: the real lengths.
struct WeightSink : IOutputSink {
    std::vector<double> w, path;
    std::string name() const override { return "weights"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext&) override {}
    void beginRun(const RunSetup&) override {}
    void fill(const MainEvent& m, const Event&, const Event&, const Event&, bool success) override {
        if (success) {
            w.push_back(m.weight);
            path.push_back(m.FP.p.path);
        }
    }
    void endRun(const RunSummary&) override {}
};

int testHrsTranspLengths() {
    // kernel: one near-central ray through HRS-R and HRS-L
    for (hrs::Side side : {hrs::Side::Right, hrs::Side::Left}) {
        double path[2] = {0, 0};
        bool ok[2] = {false, false};
        for (int compat = 0; compat < 2; ++compat) {
            hrs::McHrs k(side);
            k.setLog(nullptr);
            k.setFortranCompat(compat == 1);
            hrs::McHrsArgs a;
            a.p_spec = 2500.;
            a.dpp = 0.3;
            a.x = 0.02;
            a.y = -0.01;
            a.dxdz = 0.002;
            a.dydz = -0.003;
            a.m2 = 938.27231 * 938.27231;
            a.Mh2_final = a.m2;
            Mt64Random rng(1);
            k.mcHrs(a, rng);
            path[compat] = a.pathlen;
            ok[compat] = a.ok_spec;
        }
        const std::string nm = hrs::sideParams(side).name;
        expect(ok[1] && sameBits(path[1], -3.134165696542769E+232),
               nm + " compat: FP path = -3.134165696542769E+232 (the REAL(16) 121.77333333d0 read as real*8)");
        expect(ok[0] && path[0] > 2000. && path[0] < 3500., nm + " corrected: FP path " + std::to_string(path[0]) + " cm");
    }
    // full run: HRS kaon (K+ Lambda, HRS-L e-, HRS-R K+, doing_decay = 0)
    const std::string inp = g_ref + "/runs/p14f1_kaon_hrs_probe/p14f1_kaon_hrs_probe.inp";
    if (!std::ifstream(inp)) {
        std::printf("SKIP: %s missing\n", inp.c_str());
        return 77;
    }
    SimcParameters p;
    InpReader rd(p);
    if (!rd.loadRun(inp)) return 1;
    p.ngen = 30;
    for (int compat = 0; compat < 2; ++compat) {
        driver::RunOptions opt;
        opt.quiet_stages = true;
        opt.fortran_compat = compat == 1;
        std::ostringstream log;
        opt.log = &log;
        driver::SimcRun run(p, Json{{"sinks", Json::array()}}, opt);
        auto* sink = static_cast<WeightSink*>(&run.addSink(std::make_unique<WeightSink>()));
        run.run();
        long ninf = 0, nfin = 0;
        for (double x : sink->w) (std::isinf(x) && x > 0 ? ninf : nfin) += std::isfinite(x) || std::isinf(x) ? 1 : 0;
        const std::string what = std::to_string(sink->w.size()) + " events: " + std::to_string(ninf) + " +Inf, " +
                                 std::to_string(nfin) + " finite";
        if (compat)
            expect(!sink->w.empty() && ninf == static_cast<long>(sink->w.size()),
                   "compat: HRS kaon Weight = +Inf for every event, as the Fortran probe run (" + what + ")");
        else
            expect(!sink->w.empty() && nfin == static_cast<long>(sink->w.size()),
                   "corrected: HRS kaon Weight finite for every event (" + what + ")");
    }
    return 0;
}

// F1e: rotate_haxis(45.0d0) receives +0.0 (REAL(16) literal): the VDC test is not rotated. A ray starting at y = 12.8 cm (13.2 at VDC1, 14.2 at VDC2) with
// dy/dz = 0.02 at x = 100 cm passes the unrotated VDC tests (|y| <= 14.4) and fails the 45 degree one (yt = y + 2.0).
int testHrsVdcRotation() {
    int code[2] = {-1, -1};
    for (int compat = 0; compat < 2; ++compat) {
        hrs::McHrs k(hrs::Side::Left);
        k.setLog(nullptr);
        k.setFortranCompat(compat == 1);
        optics::Track t;
        t.xs = 100.;
        t.ys = 12.8;
        t.dxdzs = 0.;
        t.dydzs = 0.02;
        k.setTrack(t);
        hrs::McHrsArgs a;
        optics::DecayState d;
        d.m2 = 0.51099906 * 0.51099906;
        d.ph = 2700.;
        hrs::McHrsDetail det;
        ScriptedRandom rng;
        code[compat] = k.mcHrsHutDirect(d, a, -15.522000000000105, rng, det);
    }
    expect(code[1] == hrs::kSuccess, "compat: no rotation, the ray passes the VDCs (code " + std::to_string(code[1]) + ")");
    expect(code[0] == hrs::kDc1, "corrected: 45 degree rotation, the ray stops at VDC1 (code " + std::to_string(code[0]) + ")");
    return 0;
}


// ------------------------------------------------------------------------------------------------ D-062 (P14.3)
struct PathSink : IOutputSink {
    std::vector<double> path, xcal1, xcal2;
    std::string name() const override { return "paths"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext&) override {}
    void beginRun(const RunSetup&) override {}
    void fill(const MainEvent& m, const Event&, const Event&, const Event&, bool success) override {
        if (success) {
            path.push_back(m.FP.p.path);
            xcal1.push_back(m.ntup.xcal_gamma1);
            xcal2.push_back(m.ntup.xcal_gamma2);
        }
    }
    void endRun(const RunSummary&) override {}
};

/// Runs `inp` (after `edit`) in one mode; returns "" or the FortranStop message; `sink` gets the events.
std::string runInp(SimcParameters p, bool compat, PathSink* sink, long* nevents = nullptr) {
    driver::RunOptions opt;
    opt.quiet_stages = true;
    opt.fortran_compat = compat;
    opt.data_dir = specDir();
    std::ostringstream log;
    opt.log = &log;
    try {
        driver::SimcRun run(p, Json{{"sinks", Json::array()}}, opt);
        PathSink* own = sink;
        if (own) {
            auto s = std::make_unique<PathSink>();
            PathSink* raw = s.get();
            run.addSink(std::move(s));
            run.run();
            *own = *raw;
        } else {
            run.run();
        }
        if (nevents) *nevents = run.totals().nevent;
    } catch (const FortranStop& e) {
        return e.what();
    }
    return "";
}

bool loadRunInp(const std::string& run, SimcParameters& p) {
    const std::string inp = g_ref + "/runs/" + run + "/" + run + ".inp";
    if (!std::ifstream(inp)) {
        std::printf("SKIP: %s missing\n", inp.c_str());
        return false;
    }
    InpReader rd(p);
    return rd.loadRun(inp);
}

// G1: pi0 photons are passed m2 = Mpi0^2 and doing_decay (simc.f:1518-1522/1535-1539): any doing_decay = 1 pi0 run STOPs
// (reference/runs/p14f2_excl_h_dec1_probe: 'cannot decay particle with mass=  134.98' in the first event); corrected:
// massless stable photons, the run completes.
int testPi0PhotonDecay() {
    SimcParameters p;
    if (!loadRunInp("p14f2_excl_h_dec1_probe", p)) return 77;
    p.ngen = 20;
    const std::string c = runInp(p, true, nullptr);
    expect(c.find("cannot decay particle with mass=  134.98") != std::string::npos, "compat: STOP '" + c + "'");
    long n = 0;
    const std::string k = runInp(p, false, nullptr, &n);
    expect(k.empty() && n == 20, "corrected: doing_decay = 1 pi0 run completes (" + std::to_string(n) + " events" +
                                     (k.empty() ? "" : ", STOP " + k) + ")");
    return 0;
}

// G2: the calorimeter electron arm passes decay_flag = doing_decay (simc.f:1761-1765). pi+ H (pion_h_k1_rad_dec1_A) with
// the electron in the calorimeter (arm 7, drift_to_cal 300 cm) and ctau = 0.001 cm: compat STOPs at the first electron
// "decay" (mass 0.51); corrected: the electron never decays, the run completes.
int testCaloElectronDecay() {
    SimcParameters p;
    if (!loadRunInp("pion_h_k1_rad_dec1_A", p)) return 77;
    p.electron_arm = 7;
    p.drift_to_cal = 300.;
    p.ctau = 0.001;
    p.ngen = -3000;
    const std::string c = runInp(p, true, nullptr);
    expect(c.find("cannot decay particle with mass=    0.51") != std::string::npos, "compat: STOP '" + c + "'");
    const std::string k = runInp(p, false, nullptr);
    expect(k.empty(), "corrected: the calorimeter electron does not decay, the run completes" + (k.empty() ? "" : " (" + k + ")"));
    return 0;
}

// G3: pathlen is not reset for the second photon (simc.f:1462): FP%p%path = both photons (~2 x 300 cm); corrected: the
// path of the photon reported in FP%p (~300 cm).
int testPi0PathLength() {
    SimcParameters p;
    if (!loadRunInp("p14f2_excl_h_ng2_cap", p)) return 77;
    p.ngen = 30;
    PathSink a, b;
    runInp(p, true, &a);
    runInp(p, false, &b);
    auto range = [](const std::vector<double>& v) {
        double lo = 1e99, hi = -1e99;
        for (double x : v) lo = std::min(lo, x), hi = std::max(hi, x);
        return std::make_pair(lo, hi);
    };
    const auto ra = range(a.path), rb = range(b.path);
    char m[200];
    std::snprintf(m, sizeof m, "compat: FP%%p%%path in [%.2f, %.2f] cm (two photons, drift 300 cm)", ra.first, ra.second);
    expect(!a.path.empty() && ra.first > 600. && ra.second < 620., m);
    std::snprintf(m, sizeof m, "corrected: FP%%p%%path in [%.2f, %.2f] cm (one photon)", rb.first, rb.second);
    expect(!b.path.empty() && rb.first > 300. && rb.second < 310., m);
    return 0;
}

const std::map<std::string, std::function<int()>>& items() {
    static const std::map<std::string, std::function<int()>> m = {
        {"D032_spectrometer_id", testSpectrometerId},
        {"D061_F1_hrs_transp_lengths", testHrsTranspLengths},
        {"D061_F1e_hrs_vdc_rotation", testHrsVdcRotation},
        {"D062_G1_pi0_photon_decay", testPi0PhotonDecay},
        {"D062_G2_calo_electron_decay", testCaloElectronDecay},
        {"D062_G3_pi0_path_length", testPi0PathLength},
        {"D006_pass_cuts", testPassCuts},
        {"D017_extreme_trip_hadron_angle", testExtremeTrip},
        {"D014_deforest_cos_phi_warning", testDeForestWarning},
        {"D028_spectrometer_angles_nan", testSpectrometerAngles},
        {"D015_rotate_raddeg", testRotateRaddeg},
        {"D021_hms_coll_slit_counters", testHmsCollSlitCounters},
        {"D006_kaon_decay_mass", testKaonDecayMass},
        {"D011_D018_report_formats", testReportFormats},
        {"D033_F3_benhar_sf_diff", testBenharSfDiff},
        {"D033_F5_pmper_stale_pm", testPmPer},
        {"D033_F6_sf_lookup_stale", testSfLookupStale},
        {"D033_F7_heavy_mrec_guard", testMrecGuard},
        {"D037_F1_shms_coll_slit_counters", testShmsCollSlitCounters},
        {"D037_F2_coll_stops_geni", testCollStopsGeni},
        {"D037_F4_coll_ranged_out", testCollRangedOut},
        {"D037_F7_project_decay_drift", testProjectDecayDrift},
        {"D037_F9_decay_non_pi_k", testDecayNonPiK},
        {"D038_F1_delta_pp_stop", testDeltaPPStop},
        {"D038_F2_pi0_neutron_on_h", testPi0NeutronOnH},
        {"D038_F5_central_sigcc", testCentralSigcc},
        {"D043_F13_central_sigcc_kaon", testCentralSigccKaon},
        {"D038_F15_hist_which_pion_11", testHistWhichPion11},
        {"D048_F1_delta_DA_stop", testDeltaDAStop},
        {"D048_F2_delta_tail1_egamma", testDeltaTail1},
        {"D048_F3_delta_sigblok_args", testDeltaSigblok},
        {"D048_F6_delta_warning_text", testDeltaWarningText},
        {"D048_F9_central_sigcc_delta", testCentralSigccDelta},
        {"D051_F3_F9_semi_pfermi_pmper", testSemiPfermiPmPer},
        {"D051_F12_F14_semi_hist_titles", testSemiHistTitles},
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
