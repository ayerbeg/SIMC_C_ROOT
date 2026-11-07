// tests/l1/test_l1_calo.cpp
//
// L1 tests of the calorimeter arm (task P14.3) against the P14.F2 tables reference/functions/calo/ (README.txt there):
//   mc       <dir> [tol] : mc_calo (calo/mc_calo.f), scripted grnd/gauss1, every output column, counter deltas d_1..d_4
//   stops    <dir>       : mc_calo_stops — the rows that end in the project.f STOP (simc::FortranStop, same message)
//   recon    <dir> [tol] : mc_calo_recon (calo/mc_calo_recon.f, not called by SIMC)
//   rotation <dir>       : calo_frame_rotation — simc.f:1495-1515 photon rotation (driver::rotateToCalorimeter), bitwise
//   selftest             : plug-in, pi0 two-photon component in both modes (D-062 G1/G3), counters
// Fortran tables -> fortran_compat. Exit: 0 pass, 1 fail, 77 table missing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "simc/core/Errors.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/driver/MonteCarlo.h"
#include "simc/driver/PizeroCalorimeter.h"
#include "simc/spectrometers/calo/CaloSpectrometer.h"
#include "simc/spectrometers/calo/kernel/McCalo.h"

using namespace simc;
using simc::calo::McCalo;
using simc::calo::McCaloArgs;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

constexpr double kUnset = -999999.;

std::vector<ColumnSpec> present(const std::string& dir, const std::vector<ColumnSpec>& known) {
    simcl1::Table ref;
    std::string err;
    std::vector<ColumnSpec> c;
    if (!ref.load(dir + "/output.txt", err)) return c;
    for (const auto& k : known)
        if (ref.col(k.name) >= 0) c.push_back(k);
    return c;
}

McCaloArgs argsOf(const Row& r) {
    McCaloArgs a;
    a.spect_side = static_cast<int>(std::lround(r["spect_side"]));
    a.p_spec = r["p_spec"];
    a.th_spec = r["th_spec"];
    a.dpp = r["dpp"];
    a.x = r["x"];
    a.y = r["y"];
    a.z = r["z"];
    a.dxdz = r["dxdz"];
    a.dydz = r["dydz"];
    a.m2 = r["m2"];
    a.frx = r["frx"];
    a.fry = r["fry"];
    a.ms_flag = std::lround(r["ms_flag"]) != 0;
    a.wcs_flag = std::lround(r["wcs_flag"]) != 0;
    a.decay_flag = std::lround(r["decay_flag"]) != 0;
    a.ctau = r["ctau"];
    a.decdist = r["decdist_in"];
    a.Mh2_final = r["Mh2_final_in"];
    a.resmult = r["resmult_in"];
    a.pathlen = r["pathlen_in"];
    a.using_tgt_field = std::lround(r["using_tgt_field"]) != 0;
    a.prot_zbeam = r["prot_zbeam"];
    a.drift_to_cal = r["drift_to_cal"];
    a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
    return a;
}

void loadQueues(const Row& r, ScriptedRandom& rng) {
    std::vector<double> g, u;
    for (int k = 1; r.has("g" + std::to_string(k)); ++k) g.push_back(r["g" + std::to_string(k)]);
    for (int k = 1; r.has("u" + std::to_string(k)); ++k) u.push_back(r["u" + std::to_string(k)]);
    rng.load(g, u);
}

int runMc(const std::string& dir, double tol) {
    McCalo k;
    k.setFortranCompat(true);
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::vector<ColumnSpec> known;
    for (const char* n : {"ok", "stop_code", "x_fp", "dx_fp", "y_fp", "dy_fp", "resmult"}) known.push_back({n, 0.0, 0.0});
    for (const char* n : {"dpp_out", "dxdz_out", "dydz_out", "y_out", "pathlen", "m2_out", "xs_end", "ys_end", "dxdzs_end",
                          "dydzs_end", "dpps_end", "decdist", "Mh2_final"})
        known.push_back({n, tol, 0.0});
    for (const char* n : {"n_gauss", "n_unif", "rng_err", "row", "d_1", "d_2", "d_3", "d_4"}) known.push_back({n, 0.0, 0.0});
    return simcl1::run("mc_calo (" + dir + ")", dir, present(dir, known), [&](const Row& r, Outputs& o) {
        loadQueues(r, rng);
        McCaloArgs a = argsOf(r);
        const calo::CaloCounters before = k.counters();
        const calo::McCaloDetail det = k.mcCalo(a, rng);
        for (int i = 1; i <= 4; ++i) o["d_" + std::to_string(i)] = static_cast<double>(k.counters()[i] - before[i]);
        o["ok"] = a.ok_spec ? 1. : 0.;
        o["stop_code"] = det.stop_code;
        o["x_fp"] = a.x_fp;
        o["dx_fp"] = a.dx_fp;
        o["y_fp"] = a.y_fp;
        o["dy_fp"] = a.dy_fp;
        o["resmult"] = a.resmult;
        o["dpp_out"] = a.dpp;
        o["dxdz_out"] = a.dxdz;
        o["dydz_out"] = a.dydz;
        o["y_out"] = a.y;
        o["pathlen"] = a.pathlen;
        o["m2_out"] = a.m2;
        const optics::Track& t = k.track();
        o["xs_end"] = t.xs;
        o["ys_end"] = t.ys;
        o["dxdzs_end"] = t.dxdzs;
        o["dydzs_end"] = t.dydzs;
        o["dpps_end"] = t.dpps;
        o["decdist"] = a.decdist;
        o["Mh2_final"] = a.Mh2_final;
        o["n_gauss"] = static_cast<double>(rng.nGauss());
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["rng_err"] = rng.error() ? 1. : 0.;
        o["row"] = static_cast<double>(r.number());
    });
}

int runStops(const std::string& dir) {
    simcl1::Table in;
    std::string err;
    if (!in.load(dir + "/input.txt", err)) {
        std::printf("mc_calo_stops: SKIP (%s)\n", err.c_str());
        return 77;
    }
    std::ifstream f(dir + "/output.txt");
    std::string line;
    std::getline(f, line);   // header
    std::vector<std::string> msg;
    while (std::getline(f, line)) {
        std::istringstream is(line);
        int row, ex;
        is >> row >> ex;
        std::string rest;
        std::getline(is, rest);
        rest = rest.substr(rest.find_first_not_of(' '));
        const auto bar = rest.find(" | ");
        msg.push_back(bar == std::string::npos ? rest : rest.substr(0, bar));
    }
    int fail = 0;
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        Row r(in, i);
        McCalo k;
        k.setFortranCompat(true);
        ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
        loadQueues(r, rng);
        McCaloArgs a = argsOf(r);
        std::string got = "(no STOP)";
        try {
            k.mcCalo(a, rng);
        } catch (const FortranStop& e) {
            got = e.what();
        }
        const bool ok = i < msg.size() && got.find(msg[i]) != std::string::npos;
        std::printf("  row %zu: Fortran '%s' C++ '%s' %s\n", i + 1, i < msg.size() ? msg[i].c_str() : "?", got.c_str(),
                    ok ? "ok" : "FAIL");
        if (!ok) ++fail;
    }
    std::printf("mc_calo_stops: %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

int runRecon(const std::string& dir, double tol) {
    std::vector<ColumnSpec> known;
    for (const char* n : {"delta_p", "delta_t", "delta_phi", "y_tgt"}) known.push_back({n, tol, 0.0});
    return simcl1::run("mc_calo_recon (" + dir + ")", dir, present(dir, known), [&](const Row& r, Outputs& o) {
        optics::Track t;
        t.xs = r["xs"];
        t.ys = r["ys"];
        t.dxdzs = r["dxdzs"];
        t.dydzs = r["dydzs"];
        t.dpps = r["dpps"];
        double dp, dt, dphi, y;
        calo::mcCaloRecon(t, r["fry"], r["delta_y"], r["delta_z"], r["drift_to_cal"], dp, dt, dphi, y);
        o["delta_p"] = dp;
        o["delta_t"] = dt;
        o["delta_phi"] = dphi;
        o["y_tgt"] = y;
    });
}

bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

int runRotation(const std::string& dir) {
    simcl1::Table t;
    std::string err;
    if (!t.load(dir + "/output.txt", err)) {
        std::printf("calo_frame_rotation: SKIP (%s)\n", err.c_str());
        return 77;
    }
    long nbad = 0;
    for (std::size_t i = 0; i < t.rows.size(); ++i) {
        Row r(t, i);
        const std::array<double, 4> g = {r["g_E"], r["g_px"], r["g_py"], r["g_pz"]};
        const driver::PhotonSlopes s =
            driver::rotateToCalorimeter(g, r["spec_p_theta"], static_cast<int>(std::lround(r["hadron_arm"])));
        if (!sameBits(s.dxdz, r["dxdz"]) || !sameBits(s.dydz, r["dydz"])) {
            if (nbad < 5) std::printf("  row %zu differs: %.17g %.17g vs %.17g %.17g\n", i + 1, s.dxdz, s.dydz, r["dxdz"], r["dydz"]);
            ++nbad;
        }
    }
    std::printf("calo_frame_rotation: %zu rows, %ld differ: %s\n", t.rows.size(), nbad, nbad ? "FAIL" : "PASS");
    return nbad ? 1 : 0;
}

// ---------------------------------------------------------------------------------------------- selftest
int nfail = 0;
void check(bool ok, const std::string& what) {
    std::printf("  %-100s %s\n", what.c_str(), ok ? "ok" : "FAIL");
    if (!ok) ++nfail;
}

int runSelftest() {
    constexpr double Mpi0 = 134.9766;
    SimcParameters p;
    p.drift_to_cal = 300.;
    p.hadron_arm = 8;
    RunSetup setup;
    // two photons of a pi0 along the calorimeter axis (theta 0.2 rad, arm 8), one missing
    const double th = 0.2;
    auto lab = [&](double ey_cal, double ez_cal) {   // inverse of the arm-8 rotation
        return std::array<double, 4>{std::hypot(ey_cal, ez_cal), 0.0, ey_cal * std::cos(th) + ez_cal * std::sin(th),
                                     -ey_cal * std::sin(th) + ez_cal * std::cos(th)};
    };
    const std::array<double, 4> g1 = lab(10., 2000.), g2 = lab(400., 1500.);   // g2 at dy/dz 0.27: |y| = 80 cm, misses
    for (int compat = 0; compat < 2; ++compat) {
        for (int ng : {1, 2}) {
            auto sp = Registry<ISpectrometer>::create("calo");
            sp->configure(Json::object(), p);
            RunSetup s;
            s.fortran_compat = compat == 1;
            sp->init(StageContext(s));
            SpecInput in;
            in.p_spec = 4000.;
            in.th_spec = th;
            in.m2 = Mpi0 * Mpi0;
            in.decay_flag = false;
            ScriptedRandom rng;
            rng.load({}, {});
            const driver::PizeroCalorimeter pc(8, ng, compat == 1);
            const driver::PizeroCaloResult r = pc.detect(*sp, in, g1, g2, rng);
            const std::string m = std::string(compat ? " (compat" : " (corrected") + ", ngamma " + std::to_string(ng) + ")";
            check(r.ok_gamma1 && !r.ok_gamma2 && r.ok == (ng == 1), "photon 1 hits, photon 2 misses, ok_P_arm" + m);
            check(r.xcal_gamma2 == -1.0e10 && r.ycal_gamma2 == -1.0e10 && r.xcal_gamma1 == r.fp.x_fp,
                  "xcal/ycal_gamma2 = -1e10, FP = photon 1" + m);
            const double p1 = 300. * std::sqrt(1. + std::pow(10. / 2000., 2));
            const double p2 = 300. * std::sqrt(1. + std::pow(400. / 1500., 2));
            if (compat)
                check(std::fabs(r.pathlen - (p1 + p2)) < 1e-9, "D-062 G3 compat: FP path = sum of both photons" + m);
            else
                check(std::fabs(r.pathlen - p1) < 1e-9, "D-062 G3 corrected: FP path = the reported photon's" + m);
            const auto& c = sp->stopCounters();
            check(c.at("trials") == 2 && c.at("successes") == 1 && c.at("slit_hor") == 1,
                  "Calo Trials counts photons (2 trials, 1 success, 1 slit_hor)" + m);
        }
    }
    // D-062 G1: compat photons with m2 = Mpi0^2 and decay on STOP (project.f), corrected photons never decay
    {
        SimcParameters q = p;
        q.ctau = 780.45;
        for (int compat = 0; compat < 2; ++compat) {
            auto sp = Registry<ISpectrometer>::create("calo");
            sp->configure(Json::object(), q);
            RunSetup s;
            s.fortran_compat = compat == 1;
            sp->init(StageContext(s));
            SpecInput in;
            in.p_spec = 4000.;
            in.th_spec = th;
            in.m2 = Mpi0 * Mpi0;
            in.decay_flag = true;
            ScriptedRandom rng;
            rng.load({}, std::vector<double>(20, 1e-9));   // decay right away
            bool stopped = false;
            try {
                driver::PizeroCalorimeter(8, 2, compat == 1).detect(*sp, in, g1, g1, rng);
            } catch (const FortranStop& e) {
                stopped = std::string(e.what()).find("cannot decay particle with mass=  134.98") != std::string::npos;
            }
            if (compat)
                check(stopped, "D-062 G1 compat: a pi0 'photon' with doing_decay STOPs (mass 134.98)");
            else
                check(!stopped && rng.nUniform() == 0, "D-062 G1 corrected: photons massless, no decay, no random numbers");
        }
    }
    check(driver::MonteCarlo::DefaultArm(7).type == "calo" && driver::MonteCarlo::DefaultArm(8).type == "calo",
          "MonteCarlo::DefaultArm(7/8) = \"calo\"");
    std::printf("selftest: %s (%d failures)\n", nfail ? "FAIL" : "PASS", nfail);
    return nfail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s mc|stops|recon|rotation <dir> [tol] | selftest\n", argv[0]);
        return 2;
    }
    const std::string mode = argv[1];
    try {
        if (mode == "selftest") return runSelftest();
        const double tol = argc > 3 ? std::atof(argv[3]) : 0.0;
        if (argc > 2) {
            if (mode == "mc") return runMc(argv[2], tol);
            if (mode == "stops") return runStops(argv[2]);
            if (mode == "recon") return runRecon(argv[2], tol);
            if (mode == "rotation") return runRotation(argv[2]);
        }
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
