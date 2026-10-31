// tests/l1/test_l1_sos.cpp
//
// Tests of the ISpectrometer "sos" (task P14.1: sos/mc_sos.f, mc_sos_hut.f, mc_sos_recon.f).
//   mc <table_dir> [rel_tol] : a P14.F1 ray table in the drv_mc.f format of reference/functions/mc_hms (README there):
//                              input.txt = spectr(=2) p_spec th_spec dpp x y z dxdz dydz m2 fry ms_flag wcs_flag
//                              decay_flag [using_coll] [ctau] u1.. [g1..]; output.txt = the inputs + the output columns.
//                              Every output column of the reference that this adapter knows is compared (real*8:
//                              |rel| <= rel_tol; real*4 focal plane, codes, RNG counts: exact), plus the stop reasons
//                              row by row and the COMMON /spectrometer2/ counters summed over the table.
//                              stop_code = index of the sSTOP_* counter in struct_sos.inc order (2 slit_hor ... 17 cal).
//   selftest                 : no Fortran table needed: maps, RNG consumption per stop (no smearing / scripted zero
//                              gaussians), counter bookkeeping on a flat ray scan, recon sanity, the mc_sos_recon zero
//                              guard, mock sieve, decay, the plug-in adapter (raster fry), the default arm choice.
// Exit: 0 pass, 1 fail, 77 table missing.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/Mt64Random.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/driver/MonteCarlo.h"
#include "simc/optics/CosyEngine.h"
#include "simc/spectrometers/sos/SosSpectrometer.h"
#include "simc/spectrometers/sos/kernel/McSos.h"

using namespace simc;
using simc::sos::McSos;
using simc::sos::McSosArgs;
using simc::sos::McSosDetail;
using simc::sos::SosSpectrometer;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

constexpr double kUnset = -999999.;   // drv_mc.f preset of x_fp..dy_fp, resmult
constexpr double kMe2 = 0.51099906 * 0.51099906;
constexpr double kMp2 = 938.27231 * 938.27231;
constexpr double kMpi = 139.57018;
constexpr double kMmu = 105.6583755;   // constants.inc:25 (project.f:77 decay to a muon)

// ---------------------------------------------------------------------------------------------- mc (P14.F1 tables)
int runMc(const std::string& dir, double tol) {
    simcl1::Table ref;
    std::string err;
    if (!ref.load(dir + "/output.txt", err)) {
        std::printf("mc_sos (%s): SKIP (%s)\n", dir.c_str(), err.c_str());
        return 77;
    }
    SosSpectrometer plugin;
    McSos& sos = plugin.kernel();
    sos.setFortranCompat(true);   // D-031: Fortran tables -> exact Fortran behaviour
    std::ostringstream log;
    sos.setLog(&log);
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::vector<int> codes;

    // the output columns this adapter can produce, compared if present in the reference
    std::vector<ColumnSpec> known;
    for (const char* n : {"ok", "stop_code", "reached_hut"}) known.push_back({n, 0.0, 0.0});
    for (const char* n : {"x_fp", "dx_fp", "y_fp", "dy_fp"}) known.push_back({n, 0.0, 0.0});   // real*4: exact
    for (const char* n : {"dpp_out", "dxdz_out", "dydz_out", "y_out", "pathlen", "m2_out"}) known.push_back({n, tol, 0.0});
    known.push_back({"resmult", 0.0, 0.0});
    for (const char* n : {"xs_end", "ys_end", "dxdzs_end", "dydzs_end", "dpps_end"}) known.push_back({n, tol, 0.0});
    for (const char* n : {"decdist", "Mh2_final"}) known.push_back({n, tol, 0.0});
    for (const char* n : {"n_gauss", "n_unif", "rng_err", "nsig_last", "row"}) known.push_back({n, 0.0, 0.0});
    // stop_where/x_stop/y_stop: -1 for the SOS (no such COMMON); d_k = change of word k of COMMON /spectrometer2/
    for (const char* n : {"stop_where", "x_stop", "y_stop"}) known.push_back({n, 0.0, 0.0});
    for (int k = 1; k <= 21; ++k) known.push_back({"d_" + std::to_string(k), 0.0, 0.0});
    std::vector<ColumnSpec> c;
    for (const auto& k : known)
        if (ref.col(k.name) >= 0) c.push_back(k);
    std::printf("  comparing %zu output columns\n", c.size());

    int rc = simcl1::run("mc_sos (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        std::vector<double> g, u;
        for (int k = 1; r.has("g" + std::to_string(k)); ++k) g.push_back(r["g" + std::to_string(k)]);
        for (int k = 1; r.has("u" + std::to_string(k)); ++k) u.push_back(r["u" + std::to_string(k)]);
        rng.load(g, u);
        if (std::lround(r["spectr"]) != 2) throw std::runtime_error("not an SOS row");
        McSosArgs a;
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
        if (r.has("ctau")) a.ctau = r["ctau"];
        a.decdist = r.has("decdist_in") ? r["decdist_in"] : 0.;
        a.Mh2_final = r.has("Mh2_final_in") ? r["Mh2_final_in"] : a.m2;
        a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
        a.resmult = r.has("resmult_in") ? r["resmult_in"] : kUnset;   // drv_mc_sh.f: resmult = resmult_in
        a.pathlen = 0.;
        const sos::SosCounters before = sos.counters();
        const McSosDetail det = sos.mcSos(a, rng);
        for (int k = 1; k <= 21; ++k)
            o["d_" + std::to_string(k)] =
                k <= sos::kNumCounters ? static_cast<double>(sos.counters()[k] - before[k]) : 0.;
        o["stop_where"] = o["x_stop"] = o["y_stop"] = -1.;
        codes.push_back(det.stop_code);
        o["ok"] = a.ok_spec ? 1. : 0.;
        o["stop_code"] = det.stop_code;
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
        const optics::Track& t = sos.track();
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
        o["nsig_last"] = rng.lastNsig();
        o["row"] = static_cast<double>(r.number());
    });
    if (rc == 77) return rc;

    // Per-reason agreement and the COMMON /spectrometer2/ counters summed over the table.
    bool okc = true;
    const int ic = ref.col("stop_code"), ih = ref.col("reached_hut");
    if (ic >= 0) {
        std::map<int, long> nref, ncpp, nagree;
        long hutRef = 0;
        for (std::size_t i = 0; i < ref.rows.size(); ++i) {
            const int cr = static_cast<int>(std::lround(ref.rows[i][ic]));
            ++nref[cr];
            if (ih >= 0) hutRef += std::lround(ref.rows[i][ih]);
            if (i < codes.size()) {
                ++ncpp[codes[i]];
                if (codes[i] == cr) ++nagree[cr];
            }
        }
        std::printf("  stop reason agreement (row by row):\n  %4s %-10s %8s %8s %8s\n", "code", "name", "ref", "cpp",
                    "agree");
        for (const auto& kv : nref) {
            const char* nm = kv.first == 0 ? "success" : sos::counterName(kv.first);
            std::printf("  %4d %-10s %8ld %8ld %8ld\n", kv.first, nm, kv.second, ncpp[kv.first], nagree[kv.first]);
            if (nagree[kv.first] != kv.second || ncpp[kv.first] != kv.second) okc = false;
        }
        const sos::SosCounters& cn = sos.counters();
        std::printf("  counters: trials %ld (rows %zu), hut %ld (ref %ld), successes %ld (ref %ld)\n", cn[sos::kTrials],
                    ref.rows.size(), cn[sos::kHut], hutRef, cn[sos::kSuccesses], nref[0]);
        if (cn[sos::kTrials] != static_cast<long>(ref.rows.size()) || (ih >= 0 && cn[sos::kHut] != hutRef) ||
            cn[sos::kSuccesses] != nref[0])
            okc = false;
        for (int k = sos::kSlitHor; k <= sos::kCal; ++k) {
            if (k == sos::kHut) continue;
            if (cn[k] != nref[k]) {
                std::printf("  counter %s = %ld, reference rows with that stop code = %ld\n", sos::counterName(k), cn[k],
                            nref[k]);
                okc = false;
            }
        }
    }
    if (!log.str().empty()) {
        std::printf("  unexpected messages:\n%s", log.str().c_str());
        okc = false;
    }
    std::printf("  stop reasons / counters: %s\n", okc ? "PASS" : "FAIL");
    return (rc == 0 && okc) ? 0 : 1;
}

// ---------------------------------------------------------------------------------------------- transp / chain / recon / hut
std::vector<ColumnSpec> present(const std::string& dir, const std::vector<ColumnSpec>& known) {
    simcl1::Table ref;
    std::string err;
    std::vector<ColumnSpec> c;
    if (!ref.load(dir + "/output.txt", err)) return c;
    for (const auto& k : known)
        if (ref.col(k.name) >= 0) c.push_back(k);
    return c;
}

// reference/functions/sos/transp_sos: one COSY class alone (shared/transp.f, no decay)
int runTransp(const std::string& dir, double tol) {
    optics::CosyEngine eng;
    eng.setFortranCompat(true);
    std::ostringstream log;
    eng.setLog(&log);
    eng.initForward(2, McSos::defaultMatrixDir());
    std::vector<ColumnSpec> known;
    for (const char* n : {"xs_out", "dxdzs_out", "ys_out", "dydzs_out", "dpps_out", "pathlen", "delta_z", "driftdist"})
        known.push_back({n, tol, 0.0});
    known.push_back({"adrift", 0.0, 0.0});
    const int rc = simcl1::run("transp_sos (" + dir + ")", dir, present(dir, known), [&](const Row& r, Outputs& o) {
        const int sp = static_cast<int>(std::lround(r["spectr"])), cls = static_cast<int>(std::lround(r["class"]));
        optics::Track t;
        t.xs = r["xs"];
        t.dxdzs = r["dxdzs"];
        t.ys = r["ys"];
        t.dydzs = r["dydzs"];
        t.dpps = r["dpps"];
        double pathlen = 0.;
        const double dz = eng.transp(sp, cls, t, r["zd"], pathlen, nullptr);
        o["xs_out"] = t.xs;
        o["dxdzs_out"] = t.dxdzs;
        o["ys_out"] = t.ys;
        o["dydzs_out"] = t.dydzs;
        o["dpps_out"] = t.dpps;
        o["pathlen"] = pathlen;
        (void)dz;
        o["delta_z"] = pathlen - r["zd"];   // the table's definition (README: delta_z = pathlen - zd)
        o["adrift"] = eng.adrift(sp, cls) ? 1. : 0.;
        o["driftdist"] = eng.driftdist(sp, cls);
    });
    if (rc == 0 && !log.str().empty()) std::printf("  messages:\n%s", log.str().c_str());
    return rc;
}

// reference/functions/sos/transp_sos_chain: the aperture-free project/transp/rotate_haxis sequence of mc_sos.f
int runChain(const std::string& dir, double tol) {
    optics::CosyEngine eng;
    eng.setFortranCompat(true);
    std::ostringstream log;
    eng.setLog(&log);
    eng.initForward(2, McSos::defaultMatrixDir());
    const char* steps[] = {"coll_ent", "coll_ext", "quad_in", "quad_mid", "quad_out", "bm01_in", "bm01_out",
                           "bm02_in", "bm02_out", "exit_old", "exit_new", "exit_win", "hut_in", "fp"};
    std::vector<ColumnSpec> known;
    for (const char* st : steps)
        for (const char* q : {"_x", "_y", "_dxdz", "_dydz", "_path", "_xt", "_yt"})
            known.push_back({std::string(st) + q, tol, 0.0});
    for (const char* n : {"cosy_x", "cosy_dxdz", "cosy_y", "cosy_dydz", "cosy_dz"}) known.push_back({n, tol, 0.0});
    constexpr double z_entr = 126.3e0, z_exit = z_entr + 6.3e0;   // mc_sos.f:66-67
    return simcl1::run("transp_sos_chain (" + dir + ")", dir, present(dir, known), [&](const Row& r, Outputs& o) {
        optics::Track t;
        t.xs = r["x"];
        t.ys = r["y"];
        t.dxdzs = r["dxdz"];
        t.dydzs = r["dydz"];
        t.dpps = r["dpp"];
        double path = 0.;
        auto rec = [&](const std::string& st, double xt, double yt) {
            o[st + "_x"] = t.xs;
            o[st + "_y"] = t.ys;
            o[st + "_dxdz"] = t.dxdzs;
            o[st + "_dydz"] = t.dydzs;
            o[st + "_path"] = path;
            o[st + "_xt"] = xt;
            o[st + "_yt"] = yt;
        };
        auto rot = [&](const std::string& st, double ang) {
            double xt = t.xs, yt = t.ys;
            optics::rotateHaxis(t, ang, xt, yt, true);
            rec(st, xt, yt);
        };
        optics::project(t, z_entr, path);
        rec("coll_ent", t.xs, t.ys);
        optics::project(t, z_exit - z_entr, path);
        rec("coll_ext", t.xs, t.ys);
        optics::project(t, eng.driftdist(2, 1) - z_exit, path);
        rec("quad_in", t.xs, t.ys);
        eng.transp(2, 2, t, 35.0e0, path);
        rec("quad_mid", t.xs, t.ys);
        eng.transp(2, 3, t, 35.0e0, path);
        rec("quad_out", t.xs, t.ys);
        eng.transp(2, 4, t, 80.0e0, path);
        rot("bm01_in", -45.0e0);
        eng.transp(2, 5, t, 169.52e0, path);
        rot("bm01_out", 45.0e0);
        eng.transp(2, 6, t, 80.80e0, path);
        rot("bm02_in", 49.0e0);
        eng.transp(2, 7, t, 77.06e0, path);
        rot("bm02_out", 57.0e0);
        eng.transp(2, 8, t, 43.82e0, path);
        rot("exit_old", 45.0e0);
        eng.transp(2, 9, t, 44.34e0, path);
        rec("exit_new", t.xs, t.ys);
        t.xs = t.xs + 7.62 * t.dxdzs;
        t.ys = t.ys + 7.62 * t.dydzs;
        rec("exit_win", t.xs, t.ys);
        rec("hut_in", -3.206267e0, 0.0);   // the harness records xt = zinit, yt = 0 (no rotated plane)
        optics::project(t, 3.22, path);
        rec("fp", t.xs, t.ys);
        optics::Track c;
        c.xs = r["x"];
        c.ys = r["y"];
        c.dxdzs = r["dxdz"];
        c.dydzs = r["dydz"];
        c.dpps = r["dpp"];
        double p2 = 0., dz = 0.;
        for (int k = 1; k <= eng.nClasses(2); ++k) dz = dz + eng.transp(2, k, c, 0.0, p2);
        o["cosy_x"] = c.xs;
        o["cosy_dxdz"] = c.dxdzs;
        o["cosy_y"] = c.ys;
        o["cosy_dydz"] = c.dydzs;
        o["cosy_dz"] = dz;
    });
}

// reference/functions/sos/recon_sos: mc_sos_recon.f
int runRecon(const std::string& dir, double tol) {
    optics::CosyEngine eng;
    eng.initRecon(2, McSos::defaultMatrixDir());
    std::vector<ColumnSpec> known;
    for (const char* n : {"dpp_recon", "dph_recon", "dth_recon", "y_recon"}) known.push_back({n, tol, 0.0});
    return simcl1::run("recon_sos (" + dir + ")", dir, present(dir, known), [&](const Row& r, Outputs& o) {
        optics::Track t;
        t.xs = r["xfp"];
        t.dxdzs = r["xpfp"];
        t.ys = r["yfp"];
        t.dydzs = r["ypfp"];
        double dpp, dth, dph, y;
        eng.recon(2, t, r["fry"], dpp, dth, dph, y);
        o["dpp_recon"] = dpp;
        o["dph_recon"] = dph;
        o["dth_recon"] = dth;
        o["y_recon"] = y;
    });
}

// reference/functions/sos/hut_sos: mc_sos_hut.f alone (scripted RNG)
int runHut(const std::string& dir, double tol) {
    McSos sos;
    sos.setFortranCompat(true);
    std::ostringstream log;
    sos.setLog(&log);
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::vector<ColumnSpec> known;
    for (const char* n : {"ok_hut", "stop_code", "x_fp", "dx_fp", "y_fp", "dy_fp", "resmult", "dflag_out"})
        known.push_back({n, 0.0, 0.0});
    for (const char* n : {"pathlen", "m2_out", "p_out", "xs_end", "ys_end", "dxdzs_end", "dydzs_end", "dpps_end", "decdist",
                          "Mh2_final"})
        known.push_back({n, tol, 0.0});
    for (const char* n : {"stop_where", "x_stop", "y_stop", "n_gauss", "n_unif", "rng_err", "nsig_last", "row"})
        known.push_back({n, 0.0, 0.0});
    for (int k = 1; k <= 21; ++k) known.push_back({"d_" + std::to_string(k), 0.0, 0.0});
    const int rc = simcl1::run("hut_sos (" + dir + ")", dir, present(dir, known), [&](const Row& r, Outputs& o) {
        std::vector<double> g, u;
        for (int k = 1; r.has("g" + std::to_string(k)); ++k) g.push_back(r["g" + std::to_string(k)]);
        for (int k = 1; r.has("u" + std::to_string(k)); ++k) u.push_back(r["u" + std::to_string(k)]);
        rng.load(g, u);
        optics::Track t;
        t.xs = r["xs"];
        t.ys = r["ys"];
        t.dxdzs = r["dxdzs"];
        t.dydzs = r["dydzs"];
        t.dpps = r["dpps"];
        sos.setTrack(t);
        McSosArgs a;
        a.ms_flag = std::lround(r["ms_flag"]) != 0;
        a.wcs_flag = std::lround(r["wcs_flag"]) != 0;
        a.decay_flag = std::lround(r["decay_flag"]) != 0;
        a.resmult = r["resmult_in"];
        a.pathlen = r["pathlen_in"];
        a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
        optics::DecayState d;
        d.decay_flag = a.decay_flag;
        d.dflag = std::lround(r["dflag_in"]) != 0;
        d.m2 = r["m2"];
        d.ph = r["p"];
        d.ctau = r["ctau"];
        d.decdist = r["decdist_in"];
        d.Mh2_final = r["Mh2_final_in"];
        d.rng = &rng;
        d.fortran_compat = true;
        McSosDetail det;
        const sos::SosCounters before = sos.counters();
        const int code = sos.mcSosHutDirect(d, a, r["zinit"], rng, det);
        o["ok_hut"] = code == 0 ? 1. : 0.;
        o["stop_code"] = code == 0 ? -1. : static_cast<double>(code);
        o["x_fp"] = a.x_fp;
        o["dx_fp"] = a.dx_fp;
        o["y_fp"] = a.y_fp;
        o["dy_fp"] = a.dy_fp;
        o["resmult"] = a.resmult;
        o["pathlen"] = a.pathlen;
        o["m2_out"] = d.m2;
        o["p_out"] = d.ph;
        o["dflag_out"] = d.dflag ? 1. : 0.;
        const optics::Track& e = sos.track();
        o["xs_end"] = e.xs;
        o["ys_end"] = e.ys;
        o["dxdzs_end"] = e.dxdzs;
        o["dydzs_end"] = e.dydzs;
        o["dpps_end"] = e.dpps;
        o["stop_where"] = o["x_stop"] = o["y_stop"] = -1.;
        o["decdist"] = d.decdist;
        o["Mh2_final"] = d.Mh2_final;
        o["n_gauss"] = static_cast<double>(rng.nGauss());
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["rng_err"] = rng.error() ? 1. : 0.;
        o["nsig_last"] = rng.lastNsig();
        for (int k = 1; k <= 21; ++k)
            o["d_" + std::to_string(k)] =
                k <= sos::kNumCounters ? static_cast<double>(sos.counters()[k] - before[k]) : 0.;
        o["row"] = static_cast<double>(r.number());
    });
    if (rc != 77 && !log.str().empty()) {
        std::printf("  unexpected messages:\n%s", log.str().c_str());
        return 1;
    }
    return rc;
}

// ---------------------------------------------------------------------------------------------- selftest
int nfail = 0;
void check(bool ok, const std::string& what) {
    std::printf("  %-96s %s\n", what.c_str(), ok ? "ok" : "FAIL");
    if (!ok) ++nfail;
}

McSosArgs ray(double p_spec, double dpp, double x, double y, double dxdz, double dydz, double m2, bool smear = false) {
    McSosArgs a;
    a.p_spec = p_spec;
    a.th_spec = 0.5;
    a.dpp = dpp;
    a.x = x;
    a.y = y;
    a.z = 0.;
    a.dxdz = dxdz;
    a.dydz = dydz;
    a.m2 = m2;
    a.fry = x;
    a.ms_flag = smear;
    a.wcs_flag = smear;
    a.decay_flag = false;
    a.Mh2_final = m2;
    a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
    a.resmult = kUnset;
    return a;
}

int runSelftest() {
    // 1. maps
    {
        McSos s;
        std::ostringstream log;
        s.setLog(&log);
        bool loaded = true;
        try {
            s.loadMaps();
        } catch (const std::exception& e) {
            loaded = false;
            std::printf("  %s\n", e.what());
        }
        check(loaded, "SOS COSY maps load (" + McSos::defaultMatrixDir() + "/sos)");
        if (!loaded) return 1;
        optics::CosyEngine& e = s.engine();
        check(e.nClasses(2) == 10, "sos/forward_cosy.dat: 10 classes (mc_sos.f:158)");
        check(e.adrift(2, 1) && e.driftdist(2, 1) == 150., "class 1 is a 150 cm drift (mc_sos.f:202-203)");
        check(e.reconTerms(2).size() == 461, "sos/recon_cosy.dat: 461 terms (got " +
                                                 std::to_string(e.reconTerms(2).size()) + ")");
        check(log.str().empty(), "no transp_init length warnings");
    }

    // 2. RNG consumption: no smearing -> 1 grnd for a ray reaching the hut, 0 gauss1; nothing before the hut
    McSos s;
    std::ostringstream log;
    s.setLog(&log);
    s.setFortranCompat(true);
    ScriptedRandom rng;
    {
        McSosArgs a = ray(1000., 1.0, 0.1, -0.05, 0.002, -0.003, kMe2);
        rng.load({}, {0.5});
        const McSosDetail d = s.mcSos(a, rng);
        check(a.ok_spec && d.stop_code == 0 && d.reached_hut && d.fit_done, "near-central electron passes");
        check(rng.nUniform() == 1 && rng.nGauss() == 0 && !rng.error(), "no smearing: 1 grnd (resmult), 0 gauss1");
        check(a.resmult == 1.0, "u = 0.5 -> resmult 1");
        check(std::fabs(a.dpp - 1.0) < 0.02 && std::fabs(a.dxdz - 0.002) < 5e-4 && std::fabs(a.dydz + 0.003) < 5e-4,
              "recon close to the generated ray (delta " + std::to_string(a.dpp) + ", xptar " + std::to_string(a.dxdz) +
                  ", yptar " + std::to_string(a.dydz) + ")");
        std::printf("    pathlen = %.6f cm, x_fp %.5f dx_fp %.6f y_fp %.5f dy_fp %.6f\n", a.pathlen, a.x_fp, a.dx_fp,
                    a.y_fp, a.dy_fp);
        check(a.pathlen > 1000. && a.pathlen < 1300., "pathlen to the back of the calorimeter");
        // the same ray with ms = wcs = on and all gaussians 0: same result bitwise, exactly 180 gauss1(99)
        McSosArgs b = ray(1000., 1.0, 0.1, -0.05, 0.002, -0.003, kMe2, true);
        rng.load(std::vector<double>(400, 0.0), {0.1});
        const McSosDetail db = s.mcSos(b, rng);
        check(b.ok_spec && db.stop_code == 0, "smeared (zero gaussians) passes");
        check(rng.nGauss() == 180 && rng.nUniform() == 1 && rng.lastNsig() == 99.0,
              "ms=wcs=on success: 180 gauss1(99) + 1 grnd (got " + std::to_string(rng.nGauss()) + ")");
        check(b.resmult == 2.0, "u = 0.1 -> resmult 2");
        check(b.dpp == a.dpp && b.dxdz == a.dxdz && b.dydz == a.dydz && b.y == a.y && b.x_fp == a.x_fp &&
                  b.dy_fp == a.dy_fp && b.pathlen == a.pathlen,
              "zero gaussians reproduce the unsmeared ray bitwise");
    }
    // a ray stopped at the slit: no random numbers
    {
        McSosArgs a = ray(1000., 0., 0., 0., 0.0, 0.2, kMp2);
        rng.load({}, {});
        const McSosDetail d = s.mcSos(a, rng);
        check(!a.ok_spec && d.stop_code == sos::kSlitHor && rng.nUniform() == 0 && rng.nGauss() == 0,
              "dydz = 0.2 stops at slit_hor, no random numbers");
        check(a.x_fp == kUnset && a.resmult == kUnset && a.dpp == 0. && a.dydz == 0.2,
              "outputs not written for a stopped ray");
        McSosArgs b = ray(1000., 0., 0., 0., 0.1, 0.0, kMp2);
        const McSosDetail db = s.mcSos(b, rng);
        check(db.stop_code == sos::kSlitVert, "dxdz = 0.1 stops at slit_vert");
        McSosArgs c = ray(1000., 0., 0., 0., 0.033, 0.05, kMp2);
        const McSosDetail dc = s.mcSos(c, rng);
        check(dc.stop_code == sos::kSlitOct, "corner ray stops at slit_oct (got " + std::to_string(dc.stop_code) + ")");
    }

    // 3. flat scan: counter bookkeeping, every stop code consistent with the counters, recon sanity
    {
        McSos k;
        k.setLog(&log);
        k.setFortranCompat(true);
        Mt64Random mt(20261030);
        std::map<int, long> tally;
        double s2d = 0, s2x = 0, s2y = 0;
        long nres = 0;
        const long N = 200000;
        for (long i = 0; i < N; ++i) {
            const bool point = (i % 2) == 0;   // half of the rays from a point target (recon sanity)
            const double dpp = -30. + 60. * mt.uniform();
            const double x = point ? 0. : -1. + 2. * mt.uniform();
            const double y = point ? 0. : -2. + 4. * mt.uniform();
            const double dxdz = -0.07 + 0.14 * mt.uniform();
            const double dydz = -0.09 + 0.18 * mt.uniform();
            McSosArgs a = ray(1200., dpp, x, y, dxdz, dydz, (i % 3) ? kMp2 : kMe2);
            a.fry = point ? 0. : x;
            const McSosDetail d = k.mcSos(a, mt);
            ++tally[d.stop_code];
            if (d.stop_code == 0 && point) {
                s2d += (a.dpp - dpp) * (a.dpp - dpp);
                s2x += (a.dxdz - dxdz) * (a.dxdz - dxdz);
                s2y += (a.dydz - dydz) * (a.dydz - dydz);
                ++nres;
            }
        }
        const sos::SosCounters& c = k.counters();
        std::printf("    %-10s %8s %8s\n", "counter", "value", "rays");
        long sumStops = 0;
        bool same = true;
        for (int i = 1; i <= sos::kNumCounters; ++i) {
            const long rays = i == sos::kSuccesses ? tally[0] : (i == sos::kTrials || i == sos::kHut) ? -1 : tally[i];
            std::printf("    %-10s %8ld %8ld\n", sos::counterName(i), c[i], rays);
            if (rays >= 0 && rays != c[i]) same = false;
            if (i != sos::kTrials && i != sos::kHut && i != sos::kSuccesses) sumStops += c[i];
        }
        check(same, "every stop code matches the counter it incremented");
        check(c[sos::kTrials] == N && c[sos::kTrials] == sumStops + c[sos::kSuccesses], "trials = stops + successes");
        check(c[sos::kHut] == c[sos::kDc1] + c[sos::kDc2] + c[sos::kScin] + c[sos::kCal] + c[sos::kSuccesses],
              "hut = dc1 + dc2 + scin + cal + successes");
        check(c[sos::kCal] == 0, "sSTOP_cal never incremented (mc_sos_hut.f:519-525 commented out)");
        int nhit = 0;
        for (int i : {sos::kSlitHor, sos::kSlitVert, sos::kSlitOct, sos::kQuadIn, sos::kQuadMid, sos::kQuadOut,
                      sos::kBm01In, sos::kBm01Out, sos::kBm02In, sos::kBm02Out, sos::kExit, sos::kDc1, sos::kDc2,
                      sos::kScin})
            if (c[i] > 0) ++nhit;
        std::printf("    %d of 14 stop reasons reached by the flat scan\n", nhit);
        check(c[sos::kSuccesses] > 1000 && c[sos::kSlitHor] > 0 && c[sos::kSlitVert] > 0 && c[sos::kSlitOct] > 0,
              "successes and the three slit stops are populated");
        const double rd = std::sqrt(s2d / nres), rx = std::sqrt(s2x / nres), ry = std::sqrt(s2y / nres);
        std::printf("    point-target recon rms (no smearing, %ld rays): delta %.4g %%, xptar %.4g, yptar %.4g\n", nres,
                    rd, rx, ry);
        check(nres > 1000 && rd < 0.1 && rx < 2e-3 && ry < 2e-3, "recon rms small (forward and recon maps agree)");
        check(log.str().empty(), "no messages");
    }

    // 4. mc_sos_recon zero guard: EVERY hut value |v| <= 1e-30 -> 1e-30 (mc_sos_recon.f:120-122)
    {
        optics::CosyEngine& e = s.engine();
        optics::Track t0, t1;
        t1.xs = 1e-28;      // xs/100 = 1e-30
        t1.dxdzs = 1e-30;
        t1.ys = 1e-28;
        t1.dydzs = 1e-30;
        double a0, b0, c0, d0, a1, b1, c1, d1;
        e.recon(2, t0, 0.0, a0, b0, c0, d0);
        e.recon(2, t1, 1e-28, a1, b1, c1, d1);
        check(a0 == a1 && b0 == b1 && c0 == c1 && d0 == d1, "recon(2) of a zero track == recon of the 1e-30 track");
        check(std::isfinite(a0) && std::isfinite(d0), "recon(2) of a zero track is finite");
    }

    // 5. mock sieve (mc_sos.f:122-137): 2 grnd per call before the transport
    {
        McSos k;
        k.setLog(&log);
        k.setUseSieve(true);
        McSosArgs a = ray(1000., 0.5, 0.0, 0.0, 0.004, 0.003, kMe2);
        rng.load({}, {0.25, 0.5, 0.5});
        const McSosDetail d = k.mcSos(a, rng);
        check(d.stop_code == 0 && rng.nUniform() == 3, "use_sieve: 2 grnd + 1 (hut) for a success");
    }

    // 6. decay in flight (pions, short ctau): decays happen, m2 changes to the muon mass, decdist accumulated
    {
        McSos k;
        k.setLog(&log);
        Mt64Random mt(7);
        long ndec = 0, nok = 0;
        for (int i = 0; i < 20000; ++i) {
            McSosArgs a = ray(800., -10. + 20. * mt.uniform(), 0., 0., -0.04 + 0.08 * mt.uniform(),
                              -0.06 + 0.12 * mt.uniform(), kMpi * kMpi);
            a.decay_flag = true;
            a.ctau = 780.45;
            k.mcSos(a, mt);
            if (a.ok_spec) ++nok;
            if (a.m2 != kMpi * kMpi) {
                ++ndec;
                if (std::fabs(std::sqrt(a.m2) - kMmu) > 0.01 || a.Mh2_final != a.m2) {
                    check(false, "decay product is a muon with Mh2_final = m2");
                    break;
                }
            }
        }
        std::printf("    pions: %ld decays, %ld successes of 20000\n", ndec, nok);
        check(ndec > 100 && nok > 100, "pion decays in flight occur in the SOS");
    }

    // 7. the plug-in: registry, counters, raster fry (simc.f:1472) as the recon argument
    {
        auto sp = Registry<ISpectrometer>::create("sos");
        SimcParameters p;
        sp->configure(Json::object(), p);
        RunSetup setup;
        sp->init(StageContext(setup));
        SpecInput in;
        in.p_spec = 1000.;
        in.th_spec = 0.5;
        in.delta = 1.0;
        in.x = 0.1;
        in.y = -0.05;
        in.dxdz = 0.002;
        in.dydz = -0.003;
        in.m2 = kMe2;
        in.fry = 0.1;          // xtar_init: NOT used by the SOS
        in.raster_fry = 0.2;   // simc.f fry
        Mt64Random mt(3);
        const SpecResult r = sp->transport(in, mt);
        McSos k;
        k.setLog(nullptr);
        McSosArgs a = ray(1000., 1.0, 0.1, -0.05, 0.002, -0.003, kMe2);
        a.fry = 0.2;
        Mt64Random mt2(3);
        k.mcSos(a, mt2);
        check(sp->name() == "sos" && r.ok && r.resolution_factor.has_value(), "plug-in \"sos\": transport ok, resfac set");
        check(r.delta_recon == a.dpp && r.dxdz_recon == a.dxdz && r.dydz_recon == a.dydz && r.y_recon == a.y &&
                  r.x_fp == a.x_fp && r.pathlen == a.pathlen,
              "plug-in passes SpecInput::raster_fry as the recon fry (kernel with fry = 0.2, bitwise)");
        const auto& m = sp->stopCounters();
        check(m.size() == 18 && m.at("trials") == 1 && m.at("successes") == 1 && m.count("bm02_out") == 1,
              "stopCounters: 18 sSTOP_* keys");
    }

    // 8. the default arm of code 2 is the SOS plug-in
    {
        bool ok = false;
        try {
            ok = driver::MonteCarlo::DefaultArm(2).type == "sos";
        } catch (const std::exception& e) {
            std::printf("  %s\n", e.what());
        }
        check(ok, "MonteCarlo::DefaultArm(2) = \"sos\"");
        bool np = false;
        try {
            driver::MonteCarlo::DefaultArm(3);
        } catch (const NotPorted&) {
            np = true;
        }
        check(np, "MonteCarlo::DefaultArm(3) (HRS-R) still NotPorted (P14.2)");
    }

    std::printf("selftest: %s (%d failures)\n", nfail ? "FAIL" : "PASS", nfail);
    return nfail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s mc|transp|chain|recon|hut <table_dir> [rel_tol] | selftest\n", argv[0]);
        return 2;
    }
    const std::string mode = argv[1];
    try {
        if (mode == "selftest") return runSelftest();
        const double tol = argc > 3 ? std::atof(argv[3]) : 1e-12;
        if (mode == "mc" && argc > 2) return runMc(argv[2], tol);
        if (mode == "transp" && argc > 2) return runTransp(argv[2], tol);
        if (mode == "chain" && argc > 2) return runChain(argv[2], tol);
        if (mode == "recon" && argc > 2) return runRecon(argv[2], tol);
        if (mode == "hut" && argc > 2) return runHut(argv[2], tol);
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
