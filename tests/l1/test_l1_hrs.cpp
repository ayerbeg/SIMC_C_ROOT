// tests/l1/test_l1_hrs.cpp
//
// Tests of the ISpectrometers "hrsr" / "hrsl" (task P14.2: hrsr/mc_hrsr.f, hrsl/mc_hrsl.f, *_hut.f, *_recon.f) against the
// P14.F1 tables reference/functions/hrs{r,l}/ (drv_*_sh.f formats, as the SOS tables of test_l1_sos.cpp). The side
// follows from the spectr column (3 = HRS-R, 4 = HRS-L). The Fortran tables are compared with fortran_compat (the
// REAL(16) quirks D-061 F1 transp lengths, H2 rotate_haxis(45.0d0)).
//   mc|transp|chain|recon|hut <table_dir> [rel_tol]
//   selftest : the REAL(16) garbage lengths (README values), side constants, RNG use, counters, corrected mode, plug-ins.
// Exit: 0 pass, 1 fail, 77 table missing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
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
#include "simc/spectrometers/hrs/HrsSpectrometer.h"
#include "simc/spectrometers/hrs/kernel/McHrs.h"

using namespace simc;
using simc::hrs::McHrs;
using simc::hrs::McHrsArgs;
using simc::hrs::McHrsDetail;
using simc::hrs::Side;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

constexpr double kUnset = -999999.;
constexpr double kMe2 = 0.51099906 * 0.51099906;
constexpr double kMp2 = 938.27231 * 938.27231;

Side sideOf(double spectr) {
    const long s = std::lround(spectr);
    if (s == 3) return Side::Right;
    if (s == 4) return Side::Left;
    throw std::runtime_error("not an HRS row (spectr " + std::to_string(s) + ")");
}

std::vector<ColumnSpec> present(const std::string& dir, const std::vector<ColumnSpec>& known) {
    simcl1::Table ref;
    std::string err;
    std::vector<ColumnSpec> c;
    if (!ref.load(dir + "/output.txt", err)) return c;
    for (const auto& k : known)
        if (ref.col(k.name) >= 0) c.push_back(k);
    return c;
}

/// One kernel per side, Fortran mode, messages captured.
struct Kernels {
    std::ostringstream log;
    std::map<Side, std::unique_ptr<McHrs>> k;
    McHrs& get(Side s) {
        auto& p = k[s];
        if (!p) {
            p = std::make_unique<McHrs>(s);
            p->setFortranCompat(true);
            p->setLog(&log);
        }
        return *p;
    }
};

void loadQueues(const Row& r, ScriptedRandom& rng) {
    std::vector<double> g, u;
    for (int k = 1; r.has("g" + std::to_string(k)); ++k) g.push_back(r["g" + std::to_string(k)]);
    for (int k = 1; r.has("u" + std::to_string(k)); ++k) u.push_back(r["u" + std::to_string(k)]);
    rng.load(g, u);
}

void counterDeltas(const hrs::HrsCounters& before, const hrs::HrsCounters& after, Outputs& o) {
    for (int k = 1; k <= 21; ++k) o["d_" + std::to_string(k)] = static_cast<double>(after[k] - before[k]);
}

// ---------------------------------------------------------------------------------------------- mc
int runMc(const std::string& dir, double tol) {
    simcl1::Table ref;
    std::string err;
    if (!ref.load(dir + "/output.txt", err)) {
        std::printf("mc_hrs (%s): SKIP (%s)\n", dir.c_str(), err.c_str());
        return 77;
    }
    Kernels K;
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::vector<int> codes;
    Side side = Side::Right;
    std::vector<ColumnSpec> known;
    for (const char* n : {"ok", "stop_code", "reached_hut", "x_fp", "dx_fp", "y_fp", "dy_fp", "resmult"})
        known.push_back({n, 0.0, 0.0});
    for (const char* n : {"dpp_out", "dxdz_out", "dydz_out", "y_out", "pathlen", "m2_out", "xs_end", "ys_end",
                          "dxdzs_end", "dydzs_end", "dpps_end", "decdist", "Mh2_final", "x_stop", "y_stop"})
        known.push_back({n, tol, 0.0});
    for (const char* n : {"stop_where", "n_gauss", "n_unif", "rng_err", "nsig_last", "row"}) known.push_back({n, 0.0, 0.0});
    for (int k = 1; k <= 21; ++k) known.push_back({"d_" + std::to_string(k), 0.0, 0.0});
    const std::vector<ColumnSpec> c = present(dir, known);
    std::printf("  comparing %zu output columns\n", c.size());

    int rc = simcl1::run("mc_hrs (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        loadQueues(r, rng);
        side = sideOf(r["spectr"]);
        McHrs& k = K.get(side);
        McHrsArgs a;
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
        a.resmult = r.has("resmult_in") ? r["resmult_in"] : kUnset;
        a.pathlen = 0.;
        k.setDump({kUnset, kUnset, kUnset});   // drv_mc_sh.f presets the COMMON /g_dump_all_*/
        const hrs::HrsCounters before = k.counters();
        const McHrsDetail det = k.mcHrs(a, rng);
        counterDeltas(before, k.counters(), o);
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
        const optics::Track& t = k.track();
        o["xs_end"] = t.xs;
        o["ys_end"] = t.ys;
        o["dxdzs_end"] = t.dxdzs;
        o["dydzs_end"] = t.dydzs;
        o["dpps_end"] = t.dpps;
        o["stop_where"] = k.dump().stop_where;
        o["x_stop"] = k.dump().x_stop;
        o["y_stop"] = k.dump().y_stop;
        o["decdist"] = a.decdist;
        o["Mh2_final"] = a.Mh2_final;
        o["n_gauss"] = static_cast<double>(rng.nGauss());
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["rng_err"] = rng.error() ? 1. : 0.;
        o["nsig_last"] = rng.lastNsig();
        o["row"] = static_cast<double>(r.number());
    });
    if (rc == 77) return rc;

    bool okc = true;
    const int ic = ref.col("stop_code");
    std::map<int, long> nref, ncpp, nagree;
    for (std::size_t i = 0; i < ref.rows.size() && ic >= 0; ++i) {
        const int cr = static_cast<int>(std::lround(ref.rows[i][ic]));
        ++nref[cr];
        if (i < codes.size()) {
            ++ncpp[codes[i]];
            if (codes[i] == cr) ++nagree[cr];
        }
    }
    std::printf("  stop reason agreement (row by row):\n  %4s %-10s %8s %8s %8s\n", "code", "name", "ref", "cpp", "agree");
    for (const auto& kv : nref) {
        std::printf("  %4d %-10s %8ld %8ld %8ld\n", kv.first, kv.first == 0 ? "success" : hrs::counterName(kv.first),
                    kv.second, ncpp[kv.first], nagree[kv.first]);
        if (nagree[kv.first] != kv.second || ncpp[kv.first] != kv.second) okc = false;
    }
    const hrs::HrsCounters& cn = K.get(side).counters();
    std::printf("  counters: trials %ld (rows %zu), successes %ld (ref %ld)\n", cn[hrs::kTrials], ref.rows.size(),
                cn[hrs::kSuccesses], nref[0]);
    if (cn[hrs::kTrials] != static_cast<long>(ref.rows.size()) || cn[hrs::kSuccesses] != nref[0]) okc = false;
    // transp prints 'PROBLEM WITH TRANSFORMATION' for every REAL(16) length (D-061 F1): expected, not checked here
    std::printf("  stop reasons / counters: %s\n", okc ? "PASS" : "FAIL");
    return (rc == 0 && okc) ? 0 : 1;
}

// ---------------------------------------------------------------------------------------------- transp
int runTransp(const std::string& dir, double tol) {
    optics::CosyEngine eng;
    eng.setFortranCompat(true);
    std::ostringstream log;
    eng.setLog(&log);
    for (int sp : {3, 4}) eng.loadForward(sp, McHrs::defaultMatrixDir() + "/hrs/hrs_forward_cosy.dat");
    std::vector<ColumnSpec> known;
    for (const char* n : {"xs_out", "dxdzs_out", "ys_out", "dydzs_out", "dpps_out", "pathlen", "delta_z", "driftdist"})
        known.push_back({n, tol, 0.0});
    known.push_back({"adrift", 0.0, 0.0});
    return simcl1::run("transp_hrs (" + dir + ")", dir, present(dir, known), [&](const Row& r, Outputs& o) {
        const int sp = static_cast<int>(std::lround(r["spectr"])), cls = static_cast<int>(std::lround(r["class"]));
        optics::Track t;
        t.xs = r["xs"];
        t.dxdzs = r["dxdzs"];
        t.ys = r["ys"];
        t.dydzs = r["dydzs"];
        t.dpps = r["dpps"];
        double pathlen = 0.;
        eng.transp(sp, cls, t, r["zd"], pathlen, nullptr);
        o["xs_out"] = t.xs;
        o["dxdzs_out"] = t.dxdzs;
        o["ys_out"] = t.ys;
        o["dydzs_out"] = t.dydzs;
        o["dpps_out"] = t.dpps;
        o["pathlen"] = pathlen;
        o["delta_z"] = pathlen - r["zd"];   // the table's definition
        o["adrift"] = eng.adrift(sp, cls) ? 1. : 0.;
        o["driftdist"] = eng.driftdist(sp, cls);
    });
}

// ---------------------------------------------------------------------------------------------- chain
// The aperture-free project/transp/rotate_haxis sequence of mc_hrs*.f (fortran_compat: REAL(16) transp lengths).
int runChain(const std::string& dir, double tol) {
    Kernels K;
    const char* steps[] = {"circ1",    "circ2",    "coll_ent", "coll_ext", "pre_q1",   "q1_in",    "q1_mid",
                           "q1_out",   "post_q1a", "post_q1b", "q2_in",    "q2_mid",   "q2_out",   "post_q2a",
                           "post_q2b", "post_q2c", "d1_in",    "d1_out",   "post_d1a", "post_d1b", "q3_in",
                           "q3_mid",   "q3_out",   "post_q3",  "vac_win",  "hut_in",   "fp"};
    std::vector<ColumnSpec> known;
    for (const char* st : steps)
        for (const char* q : {"_x", "_y", "_dxdz", "_dydz", "_path", "_xt", "_yt"})
            known.push_back({std::string(st) + q, tol, 0.0});
    for (const char* n : {"cosy_x", "cosy_dxdz", "cosy_y", "cosy_dydz", "cosy_dz"}) known.push_back({n, tol, 0.0});
    return simcl1::run("transp_hrs_chain (" + dir + ")", dir, present(dir, known), [&](const Row& r, Outputs& o) {
        McHrs& k = K.get(sideOf(r["spectr"]));
        k.loadMaps();
        optics::CosyEngine& eng = k.engine();
        const hrs::HrsSideParams& P = k.params();
        const int sp = P.spectr;
        optics::Track t;
        t.xs = r["x"];
        t.ys = r["y"];
        t.dxdzs = r["dxdz"];
        t.dydzs = r["dydz"];
        t.dpps = r["dpp"];
        double path = 0., ztmp;
        auto rec = [&](const std::string& st, double xt, double yt) {
            o[st + "_x"] = t.xs;
            o[st + "_y"] = t.ys;
            o[st + "_dxdz"] = t.dxdzs;
            o[st + "_dydz"] = t.dydzs;
            o[st + "_path"] = path;
            o[st + "_xt"] = xt;
            o[st + "_yt"] = yt;
        };
        auto pr = [&](const std::string& st, double z) {
            optics::project(t, z, path);
            rec(st, t.xs, t.ys);
        };
        auto tr = [&](const std::string& st, int cls) {
            eng.transp(sp, cls, t, k.transpLength(cls), path);
            rec(st, t.xs, t.ys);
        };
        auto rot = [&](const std::string& st, double ang) {
            double xt = t.xs, yt = t.ys;
            optics::rotateHaxis(t, ang, xt, yt, true);
            rec(st, xt, yt);
        };
        pr("circ1", 65.686);
        ztmp = 65.686;
        pr("circ2", 80.436 - ztmp);
        ztmp = 80.436;
        pr("coll_ent", P.z_entr - ztmp);
        pr("coll_ext", P.z_exit - P.z_entr);
        ztmp = 135.064;
        pr("pre_q1", ztmp - P.z_exit);
        pr("q1_in", eng.driftdist(sp, 1) - ztmp);
        tr("q1_mid", 2);
        tr("q1_out", 3);
        double zd = 300.464 - 253.16;
        ztmp = zd;
        pr("post_q1a", zd);
        zd = 314.464 - 300.464;
        ztmp = ztmp + zd;
        pr("post_q1b", zd);
        pr("q2_in", eng.driftdist(sp, 4) - ztmp);
        tr("q2_mid", 5);
        tr("q2_out", 6);
        zd = 609.664 - 553.020;
        ztmp = zd;
        pr("post_q2a", zd);
        zd = 641.800 - 609.664;
        ztmp = ztmp + zd;
        pr("post_q2b", zd);
        zd = 819.489 - 641.800;
        ztmp = ztmp + zd;
        pr("post_q2c", zd);
        optics::project(t, eng.driftdist(sp, 7) - ztmp, path);
        rot("d1_in", -30.0);
        eng.transp(sp, 8, t, k.transpLength(8), path);
        rot("d1_out", 30.0);
        zd = 1745.33546 - 1655.83446;
        ztmp = zd;
        pr("post_d1a", zd);
        zd = 1759.00946 - 1745.33546;
        ztmp = ztmp + zd;
        pr("post_d1b", zd);
        pr("q3_in", eng.driftdist(sp, 9) - ztmp);
        tr("q3_mid", 10);
        tr("q3_out", 11);
        zd = 2080.38746 - 1997.76446;
        ztmp = zd;
        pr("post_q3", zd);
        zd = 2327.47246 - 2080.38746;
        ztmp = ztmp + zd;
        pr("vac_win", zd);
        const double zinit = -(eng.driftdist(sp, 12) - ztmp);
        rec("hut_in", zinit, 0.0);
        pr("fp", -zinit);
        optics::Track c;
        c.xs = r["x"];
        c.ys = r["y"];
        c.dxdzs = r["dxdz"];
        c.dydzs = r["dydz"];
        c.dpps = r["dpp"];
        double p2 = 0., dz = 0.;
        for (int cl = 1; cl <= eng.nClasses(sp); ++cl) dz = dz + eng.transp(sp, cl, c, 0.0, p2);
        o["cosy_x"] = c.xs;
        o["cosy_dxdz"] = c.dxdzs;
        o["cosy_y"] = c.ys;
        o["cosy_dydz"] = c.dydzs;
        o["cosy_dz"] = dz;
    });
}

// ---------------------------------------------------------------------------------------------- recon
int runRecon(const std::string& dir, double tol) {
    optics::CosyEngine eng;
    for (int sp : {3, 4}) eng.loadRecon(sp, McHrs::defaultMatrixDir() + "/hrs/hrs_recon_cosy.dat");
    std::vector<ColumnSpec> known;
    for (const char* n : {"dpp_recon", "dph_recon", "dth_recon", "y_recon"}) known.push_back({n, tol, 0.0});
    return simcl1::run("recon_hrs (" + dir + ")", dir, present(dir, known), [&](const Row& r, Outputs& o) {
        optics::Track t;
        t.xs = r["xfp"];
        t.dxdzs = r["xpfp"];
        t.ys = r["yfp"];
        t.dydzs = r["ypfp"];
        double dpp, dth, dph, y;
        eng.recon(static_cast<int>(std::lround(r["spectr"])), t, r["fry"], dpp, dth, dph, y);
        o["dpp_recon"] = dpp;
        o["dph_recon"] = dph;
        o["dth_recon"] = dth;
        o["y_recon"] = y;
    });
}

// ---------------------------------------------------------------------------------------------- hut
int runHut(const std::string& dir, double tol) {
    Kernels K;
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::vector<ColumnSpec> known;
    for (const char* n : {"ok_hut", "stop_code", "x_fp", "dx_fp", "y_fp", "dy_fp", "resmult", "dflag_out"})
        known.push_back({n, 0.0, 0.0});
    for (const char* n : {"pathlen", "m2_out", "p_out", "xs_end", "ys_end", "dxdzs_end", "dydzs_end", "dpps_end", "decdist",
                          "Mh2_final", "x_stop", "y_stop"})
        known.push_back({n, tol, 0.0});
    for (const char* n : {"stop_where", "n_gauss", "n_unif", "rng_err", "nsig_last", "row"}) known.push_back({n, 0.0, 0.0});
    for (int k = 1; k <= 21; ++k) known.push_back({"d_" + std::to_string(k), 0.0, 0.0});
    return simcl1::run("hut_hrs (" + dir + ")", dir, present(dir, known), [&](const Row& r, Outputs& o) {
        loadQueues(r, rng);
        McHrs& k = K.get(sideOf(r["spectr"]));
        optics::Track t;
        t.xs = r["xs"];
        t.ys = r["ys"];
        t.dxdzs = r["dxdzs"];
        t.dydzs = r["dydzs"];
        t.dpps = r["dpps"];
        k.setTrack(t);
        k.setDump({kUnset, kUnset, kUnset});
        McHrsArgs a;
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
        McHrsDetail det;
        const hrs::HrsCounters before = k.counters();
        const int code = k.mcHrsHutDirect(d, a, r["zinit"], rng, det);
        counterDeltas(before, k.counters(), o);
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
        const optics::Track& e = k.track();
        o["xs_end"] = e.xs;
        o["ys_end"] = e.ys;
        o["dxdzs_end"] = e.dxdzs;
        o["dydzs_end"] = e.dydzs;
        o["dpps_end"] = e.dpps;
        o["stop_where"] = k.dump().stop_where;
        o["x_stop"] = k.dump().x_stop;
        o["y_stop"] = k.dump().y_stop;
        o["decdist"] = d.decdist;
        o["Mh2_final"] = d.Mh2_final;
        o["n_gauss"] = static_cast<double>(rng.nGauss());
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["rng_err"] = rng.error() ? 1. : 0.;
        o["nsig_last"] = rng.lastNsig();
        o["row"] = static_cast<double>(r.number());
    });
}

// ---------------------------------------------------------------------------------------------- selftest
int nfail = 0;
void check(bool ok, const std::string& what) {
    std::printf("  %-100s %s\n", what.c_str(), ok ? "ok" : "FAIL");
    if (!ok) ++nfail;
}

McHrsArgs ray(double p_spec, double dpp, double x, double y, double dxdz, double dydz, double m2, bool smear = false) {
    McHrsArgs a;
    a.p_spec = p_spec;
    a.th_spec = 0.3;
    a.dpp = dpp;
    a.x = x;
    a.y = y;
    a.dxdz = dxdz;
    a.dydz = dydz;
    a.m2 = m2;
    a.fry = 0.;
    a.ms_flag = smear;
    a.wcs_flag = smear;
    a.Mh2_final = m2;
    a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
    a.resmult = kUnset;
    return a;
}

int runSelftest() {
    // 1. the REAL(16) garbage lengths (reference/functions/hrsl/transp_hrsl_chain/README.txt)
    {
        const double readme[7] = {-1.3355790489890208E-041, -1.2493079844832989E-127, -3.1341656965427691E+232,
                                  3.0220523071482735E+189,  -1.4594567076032913E-071, -1.0091251298040999E-065,
                                  -1.2355194031115920E+107};
        bool same = true;
        for (int i = 0; i < 7; ++i) {
            const double c = hrs::transpLengths()[i].compat;
            if (std::fabs(c - readme[i]) > 1e-15 * std::fabs(readme[i])) same = false;
            std::printf("    class %2d: intended %.10g, REAL(16) read as real*8 %.16e\n", hrs::transpLengths()[i].cls,
                        hrs::transpLengths()[i].intended, c);
        }
        check(same, "REAL(16) transp lengths read as real*8 = the P14.F1 values (D-061 F1)");
        check(hrs::low8(static_cast<hrs::Quad>(45)) == 0.0, "rotate_haxis(45.0d0): the real*8 dummy reads +0.0 (D-061 F1e)");
        const auto& R = hrs::sideParams(Side::Right);
        const auto& L = hrs::sideParams(Side::Left);
        check(R.spectr == 3 && L.spectr == 4 && R.z_entr == 110.0 && L.z_entr == 110.9 && R.z_exit == 118.0 &&
                  L.z_exit == 110.9 + 8.0 && R.h_exit == 3.340 && L.h_exit == 3.335 && R.vdc_offset == 0.48 &&
                  L.vdc_offset == 0.78,
              "side constants (mc_hrsr.f / mc_hrsl.f l.58, 67-68, 524-525)");
    }
    for (Side side : {Side::Right, Side::Left}) {
        const std::string nm = hrs::sideParams(side).name;
        // 2. maps, a near-central ray, no RNG without smearing; smearing with zero gaussians: 176 gauss1, same result
        McHrs k(side);
        std::ostringstream log;
        k.setLog(&log);
        k.setFortranCompat(false);
        k.loadMaps();
        check(k.engine().nClasses(k.params().spectr) == 12, nm + ": 12 forward classes");
        ScriptedRandom rng;
        McHrsArgs a = ray(1500., 0.5, 0.05, -0.02, 0.003, -0.004, kMe2);
        rng.load({}, {});
        const McHrsDetail d = k.mcHrs(a, rng);
        check(a.ok_spec && d.stop_code == 0 && rng.nUniform() == 0 && rng.nGauss() == 0 && a.resmult == 1.0,
              nm + ": near-central electron passes, no random numbers, resmult 1");
        std::printf("    %s corrected: pathlen %.6f cm, delta %.6f xptar %.6f yptar %.6f\n", nm.c_str(), a.pathlen, a.dpp,
                    a.dxdz, a.dydz);
        check(std::fabs(a.dpp - 0.5) < 0.02 && std::fabs(a.dxdz - 0.003) < 5e-4 && std::fabs(a.dydz + 0.004) < 5e-4,
              nm + ": recon close to the generated ray");
        check(a.pathlen > 2000. && a.pathlen < 3500., nm + ": corrected pathlen is physical");
        check(log.str().empty(), nm + ": corrected mode prints no transp length warnings");
        McHrsArgs b = ray(1500., 0.5, 0.05, -0.02, 0.003, -0.004, kMe2, true);
        rng.load(std::vector<double>(300, 0.0), {});
        k.mcHrs(b, rng);
        check(b.ok_spec && rng.nGauss() == 176 && b.dpp == a.dpp && b.y_fp == a.y_fp && b.pathlen == a.pathlen,
              nm + ": ms=wcs=on success: 176 gauss1, zero gaussians reproduce the ray bitwise");
        // y_fp carries the VDC offset, the recon input (ys) not
        check(b.y_fp == k.track().ys - k.params().vdc_offset, nm + ": y_fp = fit (recon input ys) - VDC offset");
        // 3. fortran_compat: garbage pathlen past Q2, same acceptance decision for this ray
        McHrs kc(side);
        std::ostringstream logc;
        kc.setLog(&logc);
        kc.setFortranCompat(true);
        McHrsArgs c = ray(1500., 0.5, 0.05, -0.02, 0.003, -0.004, kMe2);
        rng.load({}, {});
        kc.mcHrs(c, rng);
        char buf[64];
        std::snprintf(buf, sizeof buf, "%.6e", c.pathlen);
        check(c.ok_spec && c.pathlen < -1e232, nm + ": fortran_compat pathlen = REAL(16) garbage (" + buf + ")");
        check(logc.str().find("PROBLEM WITH TRANSFORMATION") != std::string::npos,
              nm + ": fortran_compat prints the transp.f length warning");
        // 4. flat scan: counter bookkeeping
        McHrs ks(side);
        ks.setLog(nullptr);
        Mt64Random mt(side == Side::Right ? 33 : 44);
        std::map<int, long> tally;
        const long N = 100000;
        for (long i = 0; i < N; ++i) {
            McHrsArgs s = ray(1500., -8. + 16. * mt.uniform(), -0.3 + 0.6 * mt.uniform(), -0.3 + 0.6 * mt.uniform(),
                              -0.07 + 0.14 * mt.uniform(), -0.04 + 0.08 * mt.uniform(), (i % 2) ? kMp2 : kMe2);
            ++tally[ks.mcHrs(s, mt).stop_code];
        }
        const hrs::HrsCounters& cn = ks.counters();
        long sum = 0;
        bool same = true;
        for (int i = 2; i <= hrs::kCal; ++i) {
            if (i == hrs::kHut) continue;
            sum += cn[i];
            if (cn[i] != tally[i]) same = false;
        }
        check(same && cn[hrs::kTrials] == N && N == sum + cn[hrs::kSuccesses] && cn[hrs::kSuccesses] == tally[0] &&
                  cn[hrs::kHut] == cn[hrs::kDc1] + cn[hrs::kDc2] + cn[hrs::kS1] + cn[hrs::kS2] + cn[hrs::kSuccesses] &&
                  cn[hrs::kCal] == 0,
              nm + ": counters consistent (trials " + std::to_string(N) + ", successes " +
                  std::to_string(cn[hrs::kSuccesses]) + "), cal never incremented");
    }
    // 5. plug-ins and default arms
    {
        auto r = Registry<ISpectrometer>::create("hrsr");
        auto l = Registry<ISpectrometer>::create("hrsl");
        check(r->name() == "hrsr" && l->name() == "hrsl", "plug-ins \"hrsr\", \"hrsl\" registered");
        SimcParameters p;
        r->configure(Json::object(), p);
        RunSetup setup;
        r->init(StageContext(setup));
        SpecInput in;
        in.p_spec = 1500.;
        in.delta = 0.5;
        in.x = 0.05;
        in.y = -0.02;
        in.dxdz = 0.003;
        in.dydz = -0.004;
        in.m2 = kMe2;
        in.fry = 0.3;
        in.raster_fry = 0.1;
        Mt64Random mt(1);
        const SpecResult res = r->transport(in, mt);
        McHrs k(Side::Right);
        k.setLog(nullptr);
        McHrsArgs a = ray(1500., 0.5, 0.05, -0.02, 0.003, -0.004, kMe2);
        a.fry = 0.1;
        Mt64Random mt2(1);
        k.mcHrs(a, mt2);
        check(res.ok && res.y_recon == a.y && res.delta_recon == a.dpp && res.resolution_factor == 1.0,
              "plug-in uses the raster fry (bitwise vs kernel)");
        check(r->stopCounters().size() == 21, "stopCounters: 21 keys");
        check(driver::MonteCarlo::DefaultArm(3).type == "hrsr" && driver::MonteCarlo::DefaultArm(4).type == "hrsl",
              "MonteCarlo::DefaultArm(3/4)");
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
        if (argc > 2) {
            if (mode == "mc") return runMc(argv[2], tol);
            if (mode == "transp") return runTransp(argv[2], tol);
            if (mode == "chain") return runChain(argv[2], tol);
            if (mode == "recon") return runRecon(argv[2], tol);
            if (mode == "hut") return runHut(argv[2], tol);
        }
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
