// tests/l1/test_l1_shms.cpp
//
// L1 tests of the ISpectrometer "shms" (task P6.3) against the Fortran harness tables (task P6.F1):
//   mc   <table_dir> [rel_tol]  : reference/functions/mc_shms (flags off) and mc_shms_smear (ms=wcs=1, scripted
//                                 g1..g210): every drv_mc.f column incl. stop code, stop_id and RNG consumption,
//                                 plus the COMMON /spectrometer_shms/ counters summed over the table.
//   p12  <table_dir>            : reference/functions/mc_p12_check: the SHMS (electron-arm) predictions ss* of
//                                 chk_p12.f for the solved target points (both modes).
//   dist <ref_dir> <config> [nray] [seed] : mc_distributions shms_nosmear | shms_smear, 1e6 Mt64 rays.
//   selftest                    : registry, transport(), counter keys, LSA stop, collimator/decay smoke tests.
// Exit: 0 pass, 1 fail, 77 table missing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "SpecDistCompare.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/Mt64Random.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/spectrometers/shms/ShmsSpectrometer.h"

using namespace simc;
using simc::shms::McShmsArgs;
using simc::shms::McShmsDetail;
using simc::shms::ShmsSpectrometer;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

constexpr double kPi = 3.141592653589793;   // constants.inc
constexpr double kDegrad = 180. / kPi;
constexpr double kMe = 0.51099906;
constexpr double kMe2 = kMe * kMe;           // Me2 = Me**2 (PARAMETER)
constexpr double kUnset = -999999.;

std::string codeName(int c) { return c == 0 ? "success" : shms::counterName(c); }

// ---------------------------------------------------------------------------------------------- mc / mc_smear
int runMc(const std::string& dir, double tol) {
    ShmsSpectrometer sp;
    std::ostringstream log;
    sp.setLog(&log);
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::vector<std::pair<int, int>> codes;

    std::vector<ColumnSpec> c;
    for (const char* n : {"ok", "stop_code", "stop_id", "reached_hut"}) c.push_back({n, 0.0, 0.0});
    for (const char* n : {"x_fp", "dx_fp", "y_fp", "dy_fp"}) c.push_back({n, 0.0, 0.0});   // real*4: exact
    for (const char* n : {"dpp_out", "dxdz_out", "dydz_out", "y_out", "pathlen", "m2_out"}) c.push_back({n, tol, 0.0});
    c.push_back({"resmult", 0.0, 0.0});
    for (const char* n : {"xs_end", "ys_end", "dxdzs_end", "dydzs_end", "dpps_end"}) c.push_back({n, tol, 0.0});
    for (const char* n : {"n_gauss", "n_unif", "rng_err", "nsig_last", "row"}) c.push_back({n, 0.0, 0.0});

    int rc = simcl1::run("mc_shms (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        std::vector<double> g, u;
        for (int k = 1; r.has("g" + std::to_string(k)); ++k) g.push_back(r["g" + std::to_string(k)]);
        for (int k = 1; r.has("u" + std::to_string(k)); ++k) u.push_back(r["u" + std::to_string(k)]);
        rng.load(g, u);
        McShmsArgs a;
        a.spectr = static_cast<int>(std::lround(r["spectr"]));
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
        a.ms_flag = std::lround(r["ms_flag"]) != 0;
        a.wcs_flag = std::lround(r["wcs_flag"]) != 0;
        a.decay_flag = std::lround(r["decay_flag"]) != 0;
        a.using_SHMScoll = std::lround(r["using_coll"]) != 0;
        a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
        a.resmult = kUnset;
        a.pathlen = 0.;
        const McShmsDetail det = sp.mcShms(a, rng);
        codes.push_back({det.stop_code, det.stop_id});
        o["ok"] = a.ok_spec ? 1. : 0.;
        o["stop_code"] = det.stop_code;
        o["stop_id"] = det.stop_id;
        o["reached_hut"] = det.reached_hut ? 1. : 0.;
        o["x_fp"] = a.x_fp;
        o["dx_fp"] = a.dx_fp;
        o["y_fp"] = a.y_fp;
        o["dy_fp"] = a.dy_fp;
        o["dpp_out"] = a.dpp;
        o["dxdz_out"] = a.dxdz;
        o["dydz_out"] = a.dydz;
        o["y_out"] = a.y;
        o["pathlen"] = a.pathlen;
        o["m2_out"] = a.m2;
        o["resmult"] = a.resmult;
        const optics::Track& t = sp.track();
        o["xs_end"] = t.xs;
        o["ys_end"] = t.ys;
        o["dxdzs_end"] = t.dxdzs;
        o["dydzs_end"] = t.dydzs;
        o["dpps_end"] = t.dpps;
        o["n_gauss"] = static_cast<double>(rng.nGauss());
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["rng_err"] = rng.error() ? 1. : 0.;
        o["nsig_last"] = rng.lastNsig();
        o["row"] = static_cast<double>(r.number());
    });
    if (rc == 77) return rc;

    simcl1::Table ref;
    std::string err;
    if (!ref.load(dir + "/output.txt", err)) return 1;
    const int ic = ref.col("stop_code"), isid = ref.col("stop_id"), ih = ref.col("reached_hut");
    std::map<std::pair<int, int>, long> nref, ncpp, nagree;
    std::map<int, long> byCode;
    long hutRef = 0;
    for (std::size_t i = 0; i < ref.rows.size(); ++i) {
        const std::pair<int, int> k{static_cast<int>(std::lround(ref.rows[i][ic])),
                                    static_cast<int>(std::lround(ref.rows[i][isid]))};
        ++nref[k];
        ++byCode[k.first];
        hutRef += std::lround(ref.rows[i][ih]);
        if (i < codes.size()) {
            ++ncpp[codes[i]];
            if (codes[i] == k) ++nagree[k];
        }
    }
    std::printf("  stop reason agreement (row by row):\n  %4s %7s %-10s %8s %8s %8s\n", "code", "stop_id", "name", "ref",
                "cpp", "agree");
    bool okc = true;
    for (const auto& kv : nref) {
        std::printf("  %4d %7d %-10s %8ld %8ld %8ld\n", kv.first.first, kv.first.second, codeName(kv.first.first).c_str(),
                    kv.second, ncpp[kv.first], nagree[kv.first]);
        if (nagree[kv.first] != kv.second || ncpp[kv.first] != kv.second) okc = false;
    }
    const shms::ShmsCounters& cn = sp.counters();
    std::printf("  counters: trials %ld (rows %zu), hut %ld (ref %ld), successes %ld (ref %ld)\n", cn[shms::kTrials],
                ref.rows.size(), cn[shms::kHut], hutRef, cn[shms::kSuccesses], byCode[0]);
    if (cn[shms::kTrials] != static_cast<long>(ref.rows.size()) || cn[shms::kHut] != hutRef ||
        cn[shms::kSuccesses] != byCode[0])
        okc = false;
    for (int k = shms::kTargHor; k <= shms::kColl; ++k) {
        if (k == shms::kHut || k == shms::kStopIdSlot) continue;
        if (cn[k] != byCode[k]) {
            std::printf("  counter %s = %ld, reference rows with that stop code = %ld\n", shms::counterName(k), cn[k],
                        byCode[k]);
            okc = false;
        }
    }
    const auto& m = sp.stopCounters();
    if (m.size() != 52 || m.at("q2_mid_dup") != cn[shms::kQ2Mid] || m.count("stop_id")) okc = false;
    if (!log.str().empty()) {
        std::printf("  unexpected messages:\n%s", log.str().c_str());
        okc = false;
    }
    std::printf("  stop reasons / counters: %s\n", okc ? "PASS" : "FAIL");
    return (rc == 0 && okc) ? 0 : 1;
}

// ---------------------------------------------------------------------------------------------- p12
// chk_p12.f:arms, electron arm = SHMS (simc.f:1647-1755).
int runP12(const std::string& dir, double tol) {
    simcl1::Table in, out;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !out.load(dir + "/output.txt", err)) {
        std::printf("p12: SKIP (%s)\n", err.c_str());
        return 77;
    }
    ShmsSpectrometer sp;
    Mt64Random rng(12345);
    const char* cols[8] = {"ssdelta", "ssyptar", "ssxptar", "ssytar", "ssxfp", "ssxpfp", "ssyfp", "ssypfp"};
    std::vector<simcl1::ColumnResult> res(8);
    for (int k = 0; k < 8; ++k) res[k].name = cols[k];
    int fail = 0;
    for (std::size_t i = 0; i < out.rows.size(); ++i) {
        Row o(out, i);
        const long ev = std::lround(o["event"]);
        std::size_t j = 0;
        while (j < in.rows.size() && std::lround(Row(in, j)["event"]) != ev) ++j;
        if (j == in.rows.size()) return 1;
        Row a(in, j);
        const double tx = o["tx"], ty = o["ty"], tz = o["tz"];
        const double pE = 8360.989;
        const double thE = std::fabs(12.5) / kDegrad;
        const double cE = std::cos(thE), sE = std::sin(thE);
        const double phiE = kPi / 2.;
        const double pe_ = pE * (1. + a["ssdeltai"] / 100.);
        McShmsArgs m;
        m.p_spec = pE;
        m.th_spec = thE;
        m.dpp = 100 * (pe_ - 0. - 0. - pE) / pE;
        m.dydz = a["ssyptari"];
        m.dxdz = a["ssxptari"];
        m.x = -ty;
        m.y = -tx * cE - tz * sE * std::sin(phiE);
        m.z = tz * cE + tx * sE * std::sin(phiE);
        m.x = m.x - m.z * m.dxdz;
        m.y = m.y - m.z * m.dydz;
        m.z = 0.0;
        m.fry = m.x;
        m.m2 = kMe2;
        m.pathlen = 0.;
        m.spectr = 5;
        sp.mcShms(m, rng);
        const double v[8] = {m.dpp, m.dydz, m.dxdz, m.y, m.x_fp, m.dx_fp, m.y_fp, m.dy_fp};
        const bool oks = std::lround(o["ok_shms"]) == 1;
        if (m.ok_spec != oks) {
            std::printf("  event %ld mode %ld: ok_shms C++ %d, Fortran %d\n", ev, std::lround(o["mode"]), m.ok_spec,
                        oks);
            ++fail;
        }
        for (int k = 0; k < 8; ++k) {
            const double f = o[cols[k]];
            const double d = simcl1::relDiff(v[k], f);
            auto& rr = res[k];
            ++rr.rows;
            if (v[k] == f) ++rr.bitwise_equal;
            if (d > (k >= 4 ? 0.0 : tol)) ++rr.failed;
            if (d > rr.max_rel) {
                rr.max_rel = d;
                rr.max_abs = std::fabs(v[k] - f);
                rr.worst_row = i + 1;
                rr.worst_ref = f;
                rr.worst_cpp = v[k];
            }
        }
    }
    std::vector<ColumnSpec> spec;
    for (int k = 0; k < 8; ++k) spec.push_back({cols[k], k >= 4 ? 0.0 : tol, 0.0});
    const bool ok = simcl1::report("mc_p12_check SHMS arm (" + dir + ")", res, spec);
    return (ok && fail == 0) ? 0 : 1;
}

// ---------------------------------------------------------------------------------------------- dist
int runDist(const std::string& refdir, const std::string& cfg, long nray, std::uint64_t seed) {
    const bool smear = cfg == "shms_smear";
    if (cfg != "shms_nosmear" && !smear) return 2;
    // drv_mcdist.f source for spectr 5
    const double p_spec = 8360.989, th = std::fabs(12.5) / kDegrad, m2 = kMe2;
    const double dlo = -12., dhi = 25., ax = 0.06, ay = 0.04, xw = 0.1, yw = 1.0;
    ShmsSpectrometer sp;
    Mt64Random rng(seed);
    std::map<std::string, std::vector<double>> v;
    std::map<simcl1::dist::CodeKey, long> ncode;
    for (long i = 0; i < nray; ++i) {
        const double d0 = dlo + (dhi - dlo) * rng.uniform();
        const double a0 = ax * (2. * rng.uniform() - 1.);
        const double b0 = ay * (2. * rng.uniform() - 1.);
        const double x0 = xw * (2. * rng.uniform() - 1.);
        const double y0 = yw * (2. * rng.uniform() - 1.);
        McShmsArgs a;
        a.p_spec = p_spec;
        a.th_spec = th;
        a.dpp = d0;
        a.dxdz = a0;
        a.dydz = b0;
        a.x = x0;
        a.y = y0;
        a.z = 0.;
        a.fry = x0;
        a.m2 = m2;
        a.ms_flag = a.wcs_flag = smear;
        a.pathlen = 0.;
        a.spectr = 5;
        const McShmsDetail det = sp.mcShms(a, rng);
        ++ncode[{det.stop_code, det.stop_id}];
        if (!a.ok_spec) continue;
        v["x_fp"].push_back(a.x_fp);
        v["dx_fp"].push_back(a.dx_fp);
        v["y_fp"].push_back(a.y_fp);
        v["dy_fp"].push_back(a.dy_fp);
        v["dpp_rec"].push_back(a.dpp);
        v["dxdz_rec"].push_back(a.dxdz);
        v["dydz_rec"].push_back(a.dydz);
        v["y_rec"].push_back(a.y);
        v["pathlen"].push_back(a.pathlen);
        v["resmult"].push_back(a.resmult);
        v["res_dpp"].push_back(a.dpp - d0);
        v["res_dxdz"].push_back(a.dxdz - a0);
        v["res_dydz"].push_back(a.dydz - b0);
        v["res_y"].push_back(a.y - y0);
    }
    std::printf("  Mt64 seed %llu\n", static_cast<unsigned long long>(seed));
    const int f = simcl1::dist::compare(refdir, cfg, nray, ncode, v,
                                        {"x_fp", "dx_fp", "y_fp", "dy_fp", "dpp_rec", "dxdz_rec", "dydz_rec", "y_rec",
                                         "pathlen", "resmult", "res_dpp", "res_dxdz", "res_dydz", "res_y"},
                                        codeName);
    return f == 77 ? 77 : (f ? 1 : 0);
}

// ---------------------------------------------------------------------------------------------- selftest
int runSelftest() {
    int fail = 0;
    auto expect = [&](bool ok, const std::string& what) {
        std::printf("  %-78s %s\n", what.c_str(), ok ? "ok" : "FAIL");
        if (!ok) ++fail;
    };
    {
        auto s = Registry<ISpectrometer>::create("shms");
        expect(s && s->name() == "shms", "Registry<ISpectrometer>::create(\"shms\")");
        SimcParameters p;
        s->configure(Json{{"ms_flag", false}}, p);
        Mt64Random rng(7);
        SpecInput in;   // mc_shms table row 1 (x = fry = -0.5, electron)
        in.p_spec = 8360.989;
        in.x = -0.5;
        in.fry = -0.5;
        in.m2 = kMe2;
        in.ms_flag = true;   // overridden by the JSON
        SpecResult r = s->transport(in, rng);
        expect(r.ok && r.stop_reason.empty() && r.x_fp == 7.3325866460800171E-001 &&
                   r.dx_fp == 2.6422363589517772E-004 && r.pathlen == 2.1517928545169530E+003 &&
                   r.delta_recon == 6.2609139619365581E-002 && r.resmult == 1.0 &&
                   s->stopCounters().at("successes") == 1 && s->stopCounters().at("q2_mid_dup") == 0,
               "transport(): mc_shms row 1 reproduced exactly (fp, pathlen, recon delta), keys");
    }
    {
        ShmsSpectrometer sp;
        std::ostringstream log;
        sp.setLog(&log);
        Mt64Random rng(1);
        McShmsArgs a;
        a.p_spec = 8000;
        a.m2 = kMe2;
        a.spectr = 6;
        bool thr = false;
        try {
            sp.mcShms(a, rng);
        } catch (const optics::FortranStop&) {
            thr = true;
        }
        expect(thr, "spectr = 6 (LSA) -> FortranStop");
    }
    {   // collimator (using_SHMScoll, pions): smoke test, no Fortran table yet
        ShmsSpectrometer sp;
        std::ostringstream log;
        sp.setLog(&log);
        Mt64Random rng(11);
        const double mpi2 = 139.57018 * 139.57018;
        long ok = 0;
        for (int i = 0; i < 20000; ++i) {
            McShmsArgs a;
            a.p_spec = 3000.;
            a.dpp = -10 + 30 * rng.uniform();
            a.dxdz = 0.06 * (2 * rng.uniform() - 1);
            a.dydz = 0.04 * (2 * rng.uniform() - 1);
            a.m2 = mpi2;
            a.using_SHMScoll = true;
            sp.mcShms(a, rng);
            ok += a.ok_spec;
        }
        std::printf("    coll run: successes %ld, coll stops %ld, slit_hor %ld (per-step counts)\n", ok,
                    sp.counters()[shms::kColl], sp.counters()[shms::kSlitHor]);
        expect(ok > 0 && sp.counters()[shms::kColl] > 0, "collimator stepping runs (pions)");
    }
    {
        ShmsSpectrometer sp;
        Mt64Random rng(13);
        const double mpi2 = 139.57018 * 139.57018;
        long decayed = 0;
        for (int i = 0; i < 2000; ++i) {
            McShmsArgs a;
            a.p_spec = 1000.;
            a.m2 = mpi2;
            a.decay_flag = true;
            a.ctau = 780.45;
            a.Mh2_final = mpi2;
            sp.mcShms(a, rng);
            if (a.m2 != mpi2) ++decayed;
        }
        std::printf("    decay run: %ld of 2000 central pions decayed\n", decayed);
        expect(decayed > 100 && decayed < 1900, "decay branch reachable, m2 updated");
    }
    std::printf("selftest: %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s mc|p12|dist|selftest ...\n", argv[0]);
        return 2;
    }
    const std::string mode = argv[1];
    try {
        if (mode == "selftest") return runSelftest();
        if (argc < 3) return 2;
        const double tol = argc > 3 && mode != "dist" ? std::atof(argv[3]) : 1e-12;
        if (mode == "mc") return runMc(argv[2], tol);
        if (mode == "p12") return runP12(argv[2], tol);
        if (mode == "dist") {
            if (argc < 4) return 2;
            const long n = argc > 4 ? std::atol(argv[4]) : 1000000L;
            const std::uint64_t seed = argc > 5 ? std::strtoull(argv[5], nullptr, 10) : 20260932ULL;
            return runDist(argv[2], argv[3], n, seed);
        }
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
