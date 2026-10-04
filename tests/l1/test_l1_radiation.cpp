// tests/l1/test_l1_radiation.cpp
//
// L1 tests of the radiative-correction chain (task P5.1) against the P5.F1 Fortran tables in
// reference/functions/<table>/ (see reference/P5.F1_report.md):
//
//   deterministic (tolerance given on the command line, 1e-12; NaN/Inf rows must match):
//     gamma spen spence schwinger lambda_dave brem bremos inter inter_prime srad extrad_friedrich extrad_phi
//     radc_init radc_init_ev basicrad_init_ev peaked_rad_weight peaked_rad_weight_r1x3
//   fixed RNG (ScriptedRandom, consumption counts compared):
//     basicrad_fixedrng generate_rad_fixedrng
//   generate_rad_capture <capture_generate_rad.txt> <radc_state dir> <fortran dir>
//     every generate_rad call (21258) of the P1.2 rad_A stream, recorded in runs/p5f1_radA_capture
//   radc_state <dir> <tol> <runs dir>: BuildRunSetup of the 22 configurations (rad flags, limits, edges) and
//     radc_init_ev at the recorded central vertex (the radc COMMON after calculate_central, central%rad)
//   rad_distributions <dir> <cfgid> <radc_state dir>: real RNG (Mt64Random, 1e6 samples): quantiles, KS, chi2
//   selftest: registry / interface
//
// generate_rad calls complete_ev after the incoming-electron photon (radc.f:324). Since P4.1 the callback is the
// real port, reactions::HydrogenElastic::completeEv (event.f:complete_ev), wired to the same target and radiator.
// Its output is also checked against the recorded pre-state columns of the table (pre_e_E, pre_p_E, pre_teff_2/3,
// pre_hardcorfac, pre_frac_*, ...).
//
// Exit: 0 pass, 1 fail, 77 table missing.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "TMath.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/random/Mt64Random.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/input/InpReader.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/radiation/kernel/Brem.h"
#include "simc/radiation/kernel/Radc.h"
#include "simc/radiation/SimcRadiator.h"
#include "simc/reactions/HydrogenElastic.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/setup/SetupBuilder.h"

using namespace simc;
using namespace simc::radiation;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

bool b(const Row& r, const char* c) { return r[c] != 0; }
int iv(const Row& r, const char* c) { return static_cast<int>(std::lround(r[c])); }

std::vector<ColumnSpec> specs(const std::vector<std::string>& names, double tol) {
    std::vector<ColumnSpec> v;
    for (const auto& n : names) v.push_back({n, tol, 0.0});
    return v;
}

int countLines(const std::string& text, const std::string& needle) {
    std::istringstream is(text);
    std::string l;
    int n = 0;
    while (std::getline(is, l))
        if (l.find(needle) != std::string::npos) ++n;
    return n;
}
std::string readFile(const std::string& path) {
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

BremFlags bremFlags(const Row& r) {
    BremFlags f;
    f.exponentiate = b(r, "exponentiate");
    f.include_hard = b(r, "include_hard");
    f.calculate_spence = b(r, "calculate_spence");
    f.produce_output = r.has("produce_output") && b(r, "produce_output");
    return f;
}

// ---------------------------------------------------------------------------------------------------------------
// deterministic tables
// ---------------------------------------------------------------------------------------------------------------
int runGamma(const std::string& dir, double tol) {
    return simcl1::run("gamma (radc.f:92-116 dead code + intrinsic = std::tgamma)", dir,
                       specs({"gamma_radc", "gamma_intrinsic"}, tol), [&](const Row& r, Outputs& o) {
                           o["gamma_radc"] = gammaRadc(r["x"]);
                           o["gamma_intrinsic"] = std::tgamma(r["x"]);
                       });
}

int runSpen(const std::string& dir, double tol) {
    return simcl1::run("spen (radc.f:743-764)", dir, specs({"spen"}, tol),
                       [&](const Row& r, Outputs& o) { o["spen"] = spen(r["x"]); });
}

int runSpence(const std::string& dir, double tol) {
    return simcl1::run("spence (brem.f:287-320)", dir, specs({"spence"}, tol),
                       [&](const Row& r, Outputs& o) { o["spence"] = spence(r["ax"]); });
}

int runSchwinger(const std::string& dir, double tol) {
    return simcl1::run("schwinger (radc.f:711-739)", dir, specs({"schwinger", "dsoft", "dhard"}, tol),
                       [&](const Row& r, Outputs& o) {
                           Event v;
                           v.Ein = r["Ein"];
                           v.e.E = r["e_E"];
                           v.e.theta = r["e_theta"];
                           v.Q2 = r["Q2"];
                           v.nu = r["nu"];
                           double ds = 0, dh = 0;
                           o["schwinger"] = schwinger(r["Ecutoff"], v, b(r, "include_hard"), r["targ_A"], r["etta"],
                                                      iv(r, "use_expon"), ds, dh);
                           o["dsoft"] = ds;
                           o["dhard"] = dh;
                       });
}

int runLambdaDave(const std::string& dir, double tol) {
    bool warned = false;
    std::ostringstream log;
    int rc = simcl1::run("lambda_dave (radc.f:768-821)", dir, specs({"lambda_dave"}, tol),
                         [&](const Row& r, Outputs& o) {
                             o["lambda_dave"] = lambdaDave(iv(r, "itail"), iv(r, "plus_flag"), b(r, "doing_proton"),
                                                           r["e1"], r["e2"], r["e3"], r["p3"], r["th"], warned, &log);
                         });
    const std::string w = "WARNING: lambda_dave(3)<0";
    const int nf = countLines(readFile(dir + "/driver.log"), w), nc = countLines(log.str(), w);
    std::printf("  warning box '%s': Fortran %d, C++ %d -> %s\n", w.c_str(), nf, nc, nf == nc ? "OK" : "MISMATCH");
    return rc != 0 ? rc : (nf == nc ? 0 : 1);
}

int runBrem(const std::string& dir, double tol) {
    std::ostringstream log;
    return simcl1::run("brem (brem.f:6-214)", dir, specs({"brem", "bsoft", "bhard", "dbsoft"}, tol),
                       [&](const Row& r, Outputs& o) {
                           double bs = 0, bh = 0, dbs = 0;
                           o["brem"] = brem(bremFlags(r), r["ein"], r["eout"], r["egamma"], b(r, "radiate_proton"),
                                            bs, bh, dbs, &log);
                           o["bsoft"] = bs;
                           o["bhard"] = bh;
                           o["dbsoft"] = dbs;
                       });
}

int runBremos(const std::string& dir, double tol) {
    std::ostringstream log;
    return simcl1::run("bremos (brem.f:344-574)", dir, specs({"bremos", "bsoft", "bhard", "dbsoft"}, tol),
                       [&](const Row& r, Outputs& o) {
                           double bs = 0, bh = 0, dbs = 0;
                           o["bremos"] = bremos(bremFlags(r), r["egamma"], r["k_ix"], r["k_iy"], r["k_iz"], r["k_fx"],
                                                r["k_fy"], r["k_fz"], r["p_ix"], r["p_iy"], r["p_iz"], r["p_fx"],
                                                r["p_fy"], r["p_fz"], r["p_fe"], b(r, "radiate_proton"), bs, bh, dbs,
                                                &log);
                           o["bsoft"] = bs;
                           o["bhard"] = bh;
                           o["dbsoft"] = dbs;
                       });
}

int runInter(const std::string& dir, double tol) {
    return simcl1::run("inter (brem.f:218-247)", dir, specs({"inter"}, tol), [&](const Row& r, Outputs& o) {
        o["inter"] = inter(b(r, "calculate_spence"), r["alpha"], r["ar1"], r["ar2"], r["e1"], r["e2"], r["de"]);
    });
}

int runInterPrime(const std::string& dir, double tol) {
    return simcl1::run("inter_prime (brem.f:581-596)", dir, specs({"inter_prime"}, tol),
                       [&](const Row& r, Outputs& o) {
                           o["inter_prime"] = interPrime(r["alpha"], r["ar1"], r["ar2"], r["de"]);
                       });
}

int runSrad(const std::string& dir, double tol) {
    return simcl1::run("srad (brem.f:251-283)", dir, specs({"dsoft", "dsoftz", "dsoftzz"}, tol),
                       [&](const Row& r, Outputs& o) {
                           std::ostringstream log;
                           const SradValues v = srad(r["ak"], r["akp"], r["eang"], r["q2"], r["ame"], r["am"], r["ap"],
                                                     r["de"], true, &log);
                           o["dsoft"] = v.dsoft;
                           o["dsoftz"] = v.dsoftz;
                           o["dsoftzz"] = v.dsoftzz;
                       });
}

int runExtradFriedrich(const std::string& dir, double tol) {
    return simcl1::run("extrad_friedrich (radc.f:650-664)", dir, specs({"dbrem", "dbrem_prime"}, tol),
                       [&](const Row& r, Outputs& o) {
                           double d = 0, dp = 0;
                           extradFriedrich(r["etatzai"], r["Ei"], r["Ecutoff"], r["trad"], d, dp);
                           o["dbrem"] = d;
                           o["dbrem_prime"] = dp;
                       });
}

int runExtradPhi(const std::string& dir, double tol) {
    return simcl1::run("extrad_phi (radc.f:668-707)", dir, specs({"extrad_phi"}, tol), [&](const Row& r, Outputs& o) {
        RadcRun run;
        RadcState st;
        run.extrad_flag = iv(r, "extrad_flag");
        run.etatzai = r["etatzai"];
        st.bt[1] = r["bt_1"];
        st.bt[2] = r["bt_2"];
        st.g[1] = r["g_1"];
        st.g[2] = r["g_2"];
        o["extrad_phi"] = extradPhi(run, st, iv(r, "itail"), r["E1"], r["E2"], r["Egamma"]);
    });
}

int runRadcInit(const std::string& dir, double tol) {
    return simcl1::run("radc_init (init.f:576-651)", dir,
                       specs({"extrad_flag", "etatzai", "produce_output", "exponentiate", "include_hard",
                              "calculate_spence"},
                             tol),
                       [&](const Row& r, Outputs& o) {
                           const RadcInitResult x = radcInit(r["targ_Z"], r["targ_L1"], r["targ_L2"], iv(r, "rad_flag"),
                                                             iv(r, "extrad_flag_in"), iv(r, "use_expon"),
                                                             b(r, "debug_1"));
                           o["extrad_flag"] = x.extrad_flag;
                           o["etatzai"] = x.etatzai;
                           o["produce_output"] = x.brem.produce_output;
                           o["exponentiate"] = x.brem.exponentiate;
                           o["include_hard"] = x.brem.include_hard;
                           o["calculate_spence"] = x.brem.calculate_spence;
                       });
}

const std::vector<std::string> kRadcCommonCols = {
    "etta",    "bt_1",    "bt_2",    "lambda_1", "lambda_2", "lambda_3", "rad_proton_this_ev", "hardcorfac",
    "g_0",     "g_1",     "g_2",     "g_3",      "g_4",      "c_int_0",  "c_int_1",            "c_int_2",
    "c_int_3", "g_int",   "c_ext_0", "c_ext_1",  "c_ext_2",  "c_ext_3",  "g_ext",              "c_0",
    "c_1",     "c_2",     "c_3",     "c_4",      "frac_1",   "frac_2",   "frac_3"};

void putRadcCommon(const RadcState& st, Outputs& o, const std::string& prefix = "") {
    o[prefix + "etta"] = st.etta;
    o[prefix + "bt_1"] = st.bt[1];
    o[prefix + "bt_2"] = st.bt[2];
    for (int i = 1; i <= 3; ++i) o[prefix + "lambda_" + std::to_string(i)] = st.lambda[i];
    o[prefix + "rad_proton_this_ev"] = st.rad_proton_this_ev;
    o[prefix + "hardcorfac"] = st.hardcorfac;
    for (int i = 0; i <= 4; ++i) o[prefix + "g_" + std::to_string(i)] = st.g[i];
    for (int i = 0; i <= 4; ++i) o[prefix + "c_" + std::to_string(i)] = st.c[i];
    for (int i = 0; i <= 3; ++i) o[prefix + "c_int_" + std::to_string(i)] = st.c_int[i];
    for (int i = 0; i <= 3; ++i) o[prefix + "c_ext_" + std::to_string(i)] = st.c_ext[i];
    o[prefix + "g_int"] = st.g_int;
    o[prefix + "g_ext"] = st.g_ext;
    for (int i = 1; i <= 3; ++i) o[prefix + "frac_" + std::to_string(i)] = st.frac[i];
}

int runRadcInitEv(const std::string& dir, double tol) {
    // One state for the whole table = one Fortran process: the SAVEd dsoft_prime of the intcor_mode=0 rows (first
    // in the table) is 0, as in a SIMC run.
    RadcState st;
    std::ostringstream log;
    return simcl1::run("radc_init_ev (init.f:655-728)", dir, specs(kRadcCommonCols, tol),
                       [&](const Row& r, Outputs& o) {
                           RadcRun run;
                           run.etatzai = r["etatzai"];
                           run.intcor_mode = iv(r, "intcor_mode");
                           run.use_offshell_rad = b(r, "use_offshell_rad");
                           run.use_expon = iv(r, "use_expon");
                           run.doing_tail[3] = b(r, "doing_tail_3");
                           run.brem = bremFlags(r);
                           run.targ_A = r["targ_A"];
                           MainEvent m;
                           m.target.teff[0] = r["teff_1"];
                           m.target.teff[1] = r["teff_2"];
                           Event v;
                           v.Ein = r["Ein"];
                           v.e.E = r["e_E"];
                           v.e.P = r["e_P"];
                           v.e.theta = r["e_theta"];
                           v.ue = {r["ue_x"], r["ue_y"], r["ue_z"]};
                           v.p.E = r["p_E"];
                           v.p.P = r["p_P"];
                           v.up = {r["up_x"], r["up_y"], r["up_z"]};
                           v.Q2 = r["Q2"];
                           v.nu = r["nu"];
                           radcInitEv(run, st, m, v, &log);
                           putRadcCommon(st, o);
                       });
}

int runBasicradInitEv(const std::string& dir, double tol) {
    std::vector<std::string> cols = {"g_0",     "g_1",     "g_2",     "g_3",     "g_4",   "c_int_0", "c_int_1",
                                     "c_int_2", "c_int_3", "g_int",   "c_ext_0", "c_ext_1", "c_ext_2", "c_ext_3",
                                     "g_ext",   "c_0",     "c_1",     "c_2",     "c_3",   "c_4"};
    return simcl1::run("basicrad_init_ev (init.f:732-813)", dir, specs(cols, tol), [&](const Row& r, Outputs& o) {
        RadcState st;
        for (int i = 1; i <= 3; ++i) st.lambda[i] = r["lambda_" + std::to_string(i)];
        st.bt[1] = r["bt_1"];
        st.bt[2] = r["bt_2"];
        st.g[4] = r["g_4"];
        basicradInitEv(st, r["e1"], r["e2"], r["e3"]);
        Outputs all;
        putRadcCommon(st, all);
        for (const auto& c : cols) o[c] = all[c];
    });
}

int runPeakedRadWeight(const std::string& dir, double tol, bool fresh_per_row) {
    RadcState shared;   // one process for the main table (SAVEd dsoft_intmin/max carried over)
    std::ostringstream log;
    const std::string title = fresh_per_row ? "peaked_rad_weight_r1x3 (fresh state per row)"
                                            : "peaked_rad_weight (radc.f:523-646)";
    return simcl1::run(title, dir, specs({"peaked_rad_weight"}, tol), [&](const Row& r, Outputs& o) {
        RadcState fresh;
        RadcState& st = fresh_per_row ? fresh : shared;
        RadcRun run;
        run.rad_flag = iv(r, "rad_flag");
        run.extrad_flag = iv(r, "extrad_flag");
        run.use_offshell_rad = b(r, "use_offshell_rad");
        run.brem = bremFlags(r);
        run.etatzai = r["etatzai"];
        st.bt[1] = r["bt_1"];
        st.bt[2] = r["bt_2"];
        st.g[1] = r["g_1"];
        st.g[2] = r["g_2"];
        st.g[4] = r["g_4"];
        st.g_int = r["g_int"];
        st.c_int[0] = r["c_int_0"];
        st.g_ext = r["g_ext"];
        st.c_ext[0] = r["c_ext_0"];
        st.c[4] = r["c_4"];
        st.rad_proton_this_ev = b(r, "rad_proton_this_ev");
        Event v;
        v.Ein = r["Ein"];
        v.e.E = r["e_E"];
        v.ue = {r["ue_x"], r["ue_y"], r["ue_z"]};
        v.p.P = r["p_P"];
        v.p.E = r["p_E"];
        v.up = {r["up_x"], r["up_y"], r["up_z"]};
        o["peaked_rad_weight"] = peakedRadWeight(run, st, v, r["Egamma"], r["emin"], r["emax"],
                                                 r["basicrad_val_reciprocal"], r["basicrad_weight"], &log);
    });
}

// ---------------------------------------------------------------------------------------------------------------
// fixed RNG
// ---------------------------------------------------------------------------------------------------------------
int runBasicradFixed(const std::string& dir, double tol) {
    return simcl1::run("basicrad_fixedrng (radc.f:3-88, scripted grnd)", dir,
                       specs({"itail_out", "Egamma", "weight", "val_reciprocal", "n_unif", "n_gauss", "rng_err"}, tol),
                       [&](const Row& r, Outputs& o) {
                           RadcState st;
                           for (int i = 0; i <= 4; ++i) {
                               st.g[i] = r["g_" + std::to_string(i)];
                               st.c[i] = r["c_" + std::to_string(i)];
                           }
                           ScriptedRandom rng;
                           rng.load({}, {r["u1"], r["u2"]});
                           int itail = iv(r, "itail");
                           double eg = 0, w = 0, vr = 0;
                           basicrad(st, itail, r["Egamma_lo"], r["Egamma_hi"], eg, w, vr, rng);
                           o["itail_out"] = itail;
                           o["Egamma"] = eg;
                           o["weight"] = w;
                           o["val_reciprocal"] = vr;
                           o["n_unif"] = static_cast<double>(rng.nUniform());
                           o["n_gauss"] = static_cast<double>(rng.nGauss());
                           o["rng_err"] = rng.error();
                       });
}

// ---- configuration of one Fortran process (rs_init) ------------------------------------------------------------
struct CfgProcess {
    SimcParameters params;
    RunSetup setup;
    material::HallCCryoTarget target;
    SimcRadiator rad;
    std::unique_ptr<reactions::HydrogenElastic> reaction;   // event.f:complete_ev (P4.1)
    double ntup_radphot{0}, ntup_radarm{0};   // COMMON /ntuple_variables/ (not part of main in the Fortran)
    std::string error;
};

std::string cfgInp(const std::string& radc_state_dir, int cfg) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "/inp/cfg%02d.inp", cfg);
    return radc_state_dir + buf;
}

/// Fortran rs_init: min_max_init + dbase_read + radc_init (BuildRunSetup), plus the target and radiator.
std::unique_ptr<CfgProcess> makeProcess(const std::string& inp, const std::string& fortran_dir) {
    auto P = std::make_unique<CfgProcess>();
    InpReader rd(P->params);
    if (!rd.loadRun(inp)) {
        P->error = "cannot read " + inp;
        return P;
    }
    SetupOptions opt;
    opt.fortran_compat = true;   // D-031: Fortran references -> exact Fortran behaviour
    opt.data_dir = fortran_dir;
    P->setup = BuildRunSetup(P->params, opt);
    if (!P->setup.ok()) {
        P->error = "setup STOP: " + P->setup.stops.front();
        return P;
    }
    P->target.setLog(nullptr);
    Json tcfg = Json::object();
    tcfg["Mh"] = P->setup.Mh;
    P->target.configure(tcfg, P->params);
    P->target.setLog(nullptr);
    P->target.init(P->setup);
    P->rad.setLog(nullptr);
    P->rad.init(P->setup);
    P->reaction = std::make_unique<reactions::HydrogenElastic>();
    P->reaction->configure(Json::object(), P->params);
    P->reaction->setLog(nullptr);
    StageContext sc(P->setup);
    sc.target = &P->target;
    sc.radiator = &P->rad;
    P->reaction->init(sc);
    return P;
}

const std::vector<std::string> kGenRadOut = {
    "success",      "n_unif",       "n_gauss",      "rng_err",      "ntail",        "radphot",      "radarm",
    "gen_weight",   "Egamma_used_1", "Egamma_used_2", "Egamma_used_3", "Egamma_min_1", "Egamma_max_1", "Egamma_min_2",
    "Egamma_max_2", "Egamma_min_3", "Egamma_max_3", "lim1",         "lim2",         "lim3",         "hardcorfac",
    "g_4",          "frac_1",       "frac_2",       "vertex_Ein",   "vertex_e_E",   "vertex_p_E",   "vertex_p_P",
    "vertex_Q2",    "teff_2",       "teff_3",       "Eloss_2",      "Eloss_3",      "orig_Ein",     "orig_e_E",
    "orig_e_P",     "orig_e_delta", "orig_p_E",     "orig_p_P",     "orig_p_delta"};
const std::vector<std::string> kGenRadPre = {"ok_pre",     "n_gauss_pre", "n_unif_pre", "rng_err_pre", "pre_e_E",
                                             "pre_p_E",    "pre_teff_2",  "pre_teff_3", "pre_hardcorfac",
                                             "pre_frac_1", "pre_frac_2",  "pre_frac_3", "pre_g_4",     "pre_c_4",
                                             "pre_rad_proton_this_ev"};

/// One drv_generate_rad row (see fortran_harness/radiation/drv_generate_rad.f) or one capture row.
void generateRadRow(CfgProcess& P, const Row& r, Outputs& o, bool capture, bool with_pre) {
    MainEvent m;
    Event v, orig;
    m.target.x = r["t_x"];
    m.target.y = r["t_y"];
    m.target.z = r["t_z"];
    m.target.rasterx = r["t_rasterx"];
    m.target.rastery = r["t_rastery"];
    m.target.Eloss[0] = r["t_Eloss_1"];
    m.target.teff[0] = r["t_teff_1"];
    m.target.Coulomb = r["t_Coulomb"];
    m.gen_weight = capture ? r["cap_pre_gen_weight"] : 1.0;
    v.Ein = r["Ein"];
    v.e.xptar = r["e_xptar"];
    v.e.yptar = r["e_yptar"];
    v.e.theta = r["e_theta"];
    v.e.phi = r["e_phi"];

    ScriptedRandom rng;
    std::vector<double> pg;
    for (int i = 1; i <= iv(r, "npg"); ++i) pg.push_back(r["pg" + std::to_string(i)]);
    rng.load(pg, {});
    const CompleteEvFn complete_ev = P.reaction->completeEvFn();
    const bool okpre = complete_ev(m, v, rng);
    m.Trec = v.Trec;
    m.sigcc = 1.0;
    RadcState& st = P.rad.state();
    if (with_pre) {
        o["ok_pre"] = okpre;
        o["n_gauss_pre"] = static_cast<double>(rng.nGauss());
        o["n_unif_pre"] = static_cast<double>(rng.nUniform());
        o["rng_err_pre"] = rng.error();
        o["pre_e_E"] = v.e.E;
        o["pre_p_E"] = v.p.E;
        o["pre_teff_2"] = m.target.teff[1];
        o["pre_teff_3"] = m.target.teff[2];
        o["pre_hardcorfac"] = st.hardcorfac;
        o["pre_frac_1"] = st.frac[1];
        o["pre_frac_2"] = st.frac[2];
        o["pre_frac_3"] = st.frac[3];
        o["pre_g_4"] = st.g[4];
        o["pre_c_4"] = st.c[4];
        o["pre_rad_proton_this_ev"] = st.rad_proton_this_ev;
    }
    if (!capture) {   // the driver zeroes them; in the capture run they are the stale COMMON values (not compared)
        for (int i = 1; i <= 3; ++i) st.Egamma_min[i] = st.Egamma_max[i] = 0.;
    }
    std::vector<double> u, g;
    for (int i = 1; i <= 6; ++i) u.push_back(r["u" + std::to_string(i)]);
    for (int i = 1; i <= 12; ++i) g.push_back(r["g" + std::to_string(i)]);
    rng.load(g, u);
    bool success = false;
    m.ntup.radphot = P.ntup_radphot;   // COMMON /ntuple_variables/: stale if generate_rad is not called
    m.ntup.radarm = P.ntup_radarm;
    if (okpre) success = P.rad.radiate(m, v, orig, rng, complete_ev);
    P.ntup_radphot = m.ntup.radphot;
    P.ntup_radarm = m.ntup.radarm;

    o["success"] = success;
    o["n_unif"] = static_cast<double>(rng.nUniform());
    o["n_gauss"] = static_cast<double>(rng.nGauss());
    o["rng_err"] = rng.error();
    o["ntail"] = st.ntail;
    o["radphot"] = m.ntup.radphot;
    o["radarm"] = m.ntup.radarm;
    o["gen_weight"] = m.gen_weight;
    for (int i = 1; i <= 3; ++i) {
        const std::string k = std::to_string(i);
        o["Egamma_used_" + k] = st.Egamma_used[i];
        o["Egamma_min_" + k] = st.Egamma_min[i];
        o["Egamma_max_" + k] = st.Egamma_max[i];
        o["lim" + k] = st.lim[i];
    }
    o["hardcorfac"] = st.hardcorfac;
    o["g_4"] = st.g[4];
    o["frac_1"] = st.frac[1];
    o["frac_2"] = st.frac[2];
    o["vertex_Ein"] = v.Ein;
    o["vertex_e_E"] = v.e.E;
    o["vertex_p_E"] = v.p.E;
    o["vertex_p_P"] = v.p.P;
    o["vertex_Q2"] = v.Q2;
    o["teff_2"] = m.target.teff[1];
    o["teff_3"] = m.target.teff[2];
    o["Eloss_2"] = m.target.Eloss[1];
    o["Eloss_3"] = m.target.Eloss[2];
    o["orig_Ein"] = orig.Ein;
    o["orig_e_E"] = orig.e.E;
    o["orig_e_P"] = orig.e.P;
    o["orig_e_delta"] = orig.e.delta;
    o["orig_p_E"] = orig.p.E;
    o["orig_p_P"] = orig.p.P;
    o["orig_p_delta"] = orig.p.delta;
}

int runGenerateRadFixed(const std::string& dir, double tol, const std::string& radc_state_dir,
                        const std::string& fortran_dir) {
    std::map<int, std::unique_ptr<CfgProcess>> procs;
    std::vector<std::string> cols = kGenRadPre;
    cols.insert(cols.end(), kGenRadOut.begin(), kGenRadOut.end());
    std::map<int, std::array<long, 3>> perSrc;   // gen_src -> rows, success
    int rc = simcl1::run("generate_rad_fixedrng (radc.f:120-519, scripted grnd/gauss1, full init per cfg)", dir,
                         specs(cols, tol), [&](const Row& r, Outputs& o) {
                             const int cfg = iv(r, "cfg");
                             auto& P = procs[cfg];
                             if (!P) {
                                 P = makeProcess(cfgInp(radc_state_dir, cfg), fortran_dir);
                                 if (!P->error.empty()) std::printf("  cfg %d: %s\n", cfg, P->error.c_str());
                             }
                             generateRadRow(*P, r, o, false, true);
                             auto& c = perSrc[iv(r, "gen_src")];
                             ++c[0];
                             if (o["success"] != 0) ++c[1];
                         });
    for (const auto& kv : perSrc)
        std::printf("  gen_src %d: %ld rows, %ld successes\n", kv.first, kv.second[0], kv.second[1]);
    return rc;
}

/// Every generate_rad call of the capture run (reference/runs/p5f1_radA_capture/capture_generate_rad.txt).
/// Compared as check_vs_capture does: all columns; Egamma_min/max_i only where lim_i != 0 (they are stale COMMON
/// values otherwise), orig_* only on success.
int runGenerateRadCapture(const std::string& path, double tol, const std::string& radc_state_dir,
                          const std::string& fortran_dir) {
    simcl1::Table t;
    std::string err;
    if (!t.load(path, err)) {
        std::printf("generate_rad_capture: SKIP (%s)\n", err.c_str());
        return 77;
    }
    auto P = makeProcess(cfgInp(radc_state_dir, 0), fortran_dir);
    if (!P->error.empty()) {
        std::printf("generate_rad_capture: %s\n", P->error.c_str());
        return 1;
    }
    std::map<std::string, std::size_t> bad;
    std::size_t rows_bad = 0, nsucc = 0;
    std::map<std::string, double> maxrel;
    for (std::size_t i = 0; i < t.rows.size(); ++i) {
        Row r(t, i);
        Outputs o;
        generateRadRow(*P, r, o, true, false);
        const bool succ = r["success"] != 0;
        nsucc += succ;
        bool any = false;
        for (const auto& c : kGenRadOut) {
            if (c.rfind("orig_", 0) == 0 && !succ) continue;
            if (c.rfind("Egamma_m", 0) == 0) {
                const std::string k = c.substr(c.size() - 1);
                if (r["lim" + k] == 0) continue;
            }
            if (c == "n_unif" || c == "n_gauss") {   // the capture counts are cap_n_unif/cap_n_gauss
                if (o[c] != r["cap_" + c]) {
                    ++bad[c];
                    any = true;
                }
                continue;
            }
            if (c == "rng_err") continue;
            const double f = r[c], x = o[c];
            const double rel = simcl1::relDiff(x, f);
            if (!(std::isnan(f) && std::isnan(x)) && rel > tol) {
                ++bad[c];
                any = true;
            }
            if (!std::isnan(rel)) maxrel[c] = std::max(maxrel[c], rel);
        }
        if (any) {
            ++rows_bad;
            if (rows_bad <= 5) std::printf("  mismatch at capture row %zu (ntried %.0f)\n", i + 1, r["cap_ntried"]);
        }
    }
    double worst = 0;
    std::string wc;
    for (const auto& kv : maxrel)
        if (kv.second > worst) {
            worst = kv.second;
            wc = kv.first;
        }
    std::printf("L1 generate_rad_capture (%s)\n  rows %zu (successes %zu), rows with a mismatch: %zu, max|rel| %.3e (%s)\n",
                path.c_str(), t.rows.size(), nsucc, rows_bad, worst, wc.empty() ? "-" : wc.c_str());
    for (const auto& kv : bad) std::printf("  column %s: %zu mismatches\n", kv.first.c_str(), kv.second);
    std::printf("  => %s\n", rows_bad ? "FAIL" : "PASS");
    return rows_bad ? 1 : 0;
}

// ---------------------------------------------------------------------------------------------------------------
// radc_state: BuildRunSetup rad quantities + radc_init_ev at the recorded central vertex
// ---------------------------------------------------------------------------------------------------------------
int runRadcState(const std::string& dir, double tol, const std::string& fortran_dir) {
    std::vector<std::string> cols = {"extrad_flag",   "doing_tail_1",  "doing_tail_2",   "doing_tail_3",
                                     "hardwired_rad", "produce_output", "exponentiate",  "include_hard",
                                     "calculate_spence", "etatzai",    "targ_L1",        "targ_L2",
                                     "Egamma_tot_max", "Egamma1_max",  "Egamma2_max",    "Egamma3_max",
                                     "Ebeam_vertex_ave", "teff_ave_1", "teff_ave_2",     "teff_ave_3",
                                     "edge_e_E_min",  "edge_e_E_max",  "edge_p_E_min",   "edge_p_E_max",
                                     "edge_Em_min",   "edge_Em_max",   "VERTEXedge_Trec_min", "VERTEXedge_Trec_max",
                                     "gen_e_E_min",   "gen_e_E_max",   "gen_p_E_min",    "gen_p_E_max",
                                     "spec_e_P",      "spec_p_P",      "targ_A"};
    for (const auto& c : kRadcCommonCols) cols.push_back(c);
    for (const auto& c : kRadcCommonCols)
        if (c != "etta" && c != "rad_proton_this_ev" && c != "g_4" && c != "c_4") cols.push_back("cr_" + c);
    cols.push_back("cr_etatzai");
    // input.txt carries the RAW extrad_flag, output.txt the value after the radc_init default under the same name:
    // rename the input column so that the runner's echo check does not mistake it for a stale table.
    simcl1::Table in, ref;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !ref.load(dir + "/output.txt", err)) {
        std::printf("L1 radc_state: SKIP (%s)\n", err.c_str());
        return 77;
    }
    for (auto& n : in.names)
        if (n == "extrad_flag") n = "extrad_flag_in";
    const auto cs = specs(cols, tol);
    auto res = simcl1::compare(in, ref, cs, [&](const Row& r, Outputs& o) {
                           const Row rr(ref, r.index());   // the central vertex is a Fortran OUTPUT (calculate_central)
                           const int cfg = iv(r, "cfg");
                           auto P = makeProcess(cfgInp(dir, cfg), fortran_dir);
                           if (!P->error.empty()) std::printf("  cfg %d: %s\n", cfg, P->error.c_str());
                           const RunSetup& s = P->setup;
                           const RadcRun& run = P->rad.run();
                           o["extrad_flag"] = run.extrad_flag;
                           for (int i = 1; i <= 3; ++i) o["doing_tail_" + std::to_string(i)] = run.doing_tail[i];
                           o["hardwired_rad"] = run.hardwired_rad;
                           o["produce_output"] = run.brem.produce_output;
                           o["exponentiate"] = run.brem.exponentiate;
                           o["include_hard"] = run.brem.include_hard;
                           o["calculate_spence"] = run.brem.calculate_spence;
                           o["etatzai"] = run.etatzai;
                           o["targ_L1"] = s.targ.L1;
                           o["targ_L2"] = s.targ.L2;
                           o["Egamma_tot_max"] = run.Egamma_tot_max;
                           o["Egamma1_max"] = run.Egamma1_max;
                           o["Egamma2_max"] = run.Egamma2_max;
                           o["Egamma3_max"] = run.Egamma3_max;
                           o["Ebeam_vertex_ave"] = s.Ebeam_vertex_ave;
                           for (int i = 1; i <= 3; ++i) o["teff_ave_" + std::to_string(i)] = s.targ.teff[i - 1].ave;
                           o["edge_e_E_min"] = run.edge_e_E.min;
                           o["edge_e_E_max"] = run.edge_e_E.max;
                           o["edge_p_E_min"] = run.edge_p_E.min;
                           o["edge_p_E_max"] = run.edge_p_E.max;
                           o["edge_Em_min"] = run.edge_Em.min;
                           o["edge_Em_max"] = run.edge_Em.max;
                           o["VERTEXedge_Trec_min"] = run.VERTEXedge_Trec.min;
                           o["VERTEXedge_Trec_max"] = run.VERTEXedge_Trec.max;
                           o["gen_e_E_min"] = s.gen.e.E.min;
                           o["gen_e_E_max"] = s.gen.e.E.max;
                           o["gen_p_E_min"] = s.gen.p.E.min;
                           o["gen_p_E_max"] = s.gen.p.E.max;
                           o["spec_e_P"] = run.spec_e_P;
                           o["spec_p_P"] = run.spec_p_P;
                           o["targ_A"] = run.targ_A;
                           // simc.f:1257-1261: main0%target%teff(2:3) = targ%teff%ave, radc_init_ev(main0, vertex0).
                           // teff(1) of main0 is the teff of the beam from the central calculation = targ%teff(1)%ave.
                           MainEvent m;
                           for (int i = 0; i < 3; ++i) m.target.teff[i] = s.targ.teff[i].ave;
                           Event v;
                           v.Ein = rr["c_Ein"];
                           v.e.E = rr["c_e_E"];
                           v.e.P = rr["c_e_P"];
                           v.e.theta = rr["c_e_theta"];
                           v.e.phi = rr["c_e_phi"];
                           v.p.E = rr["c_p_E"];
                           v.p.P = rr["c_p_P"];
                           v.Q2 = rr["c_Q2"];
                           v.nu = rr["c_nu"];
                           v.ue = {rr["c_ue_x"], rr["c_ue_y"], rr["c_ue_z"]};
                           v.up = {rr["c_up_x"], rr["c_up_y"], rr["c_up_z"]};
                           P->rad.initEvent(m, v);
                           const RadcState& st = P->rad.state();
                           putRadcCommon(st, o);
                           // central%rad (simc.f:1263-1277), filled only if using_rad
                           Outputs cr;
                           RadcState zero;
                           putRadcCommon(s.sw.using_rad ? st : zero, cr, "cr_");
                           for (const auto& kv : cr) o[kv.first] = kv.second;
                           o["cr_etatzai"] = s.sw.using_rad ? run.etatzai : 0.0;
                       }, err);
    if (!err.empty()) {
        std::printf("L1 radc_state: ERROR %s\n", err.c_str());
        return 1;
    }
    return simcl1::report("radc_state (BuildRunSetup rad part + radc_init_ev at the central vertex)", res, cs) ? 0 : 1;
}

// ---------------------------------------------------------------------------------------------------------------
// real-RNG distributions
// ---------------------------------------------------------------------------------------------------------------
struct Hist {
    std::vector<double> lo, hi, cnt;
    double under{0}, over{0};
};

bool loadHist(const std::string& path, Hist& h) {
    simcl1::Table t;
    std::string err;
    if (!t.load(path, err)) return false;
    for (const auto& r : t.rows) {
        h.lo.push_back(r[0]);
        h.hi.push_back(r[1]);
        h.cnt.push_back(r[2]);
    }
    std::ifstream f(path);
    std::string l;
    while (std::getline(f, l)) {
        auto p = l.find("underflow");
        if (!l.empty() && l[0] == '#' && p != std::string::npos) {
            std::sscanf(l.c_str() + p, "underflow %lf overflow %lf", &h.under, &h.over);
            break;
        }
    }
    return true;
}

double quantileInv(const std::vector<double>& s, double p) {   // numpy 'inverted_cdf'
    std::size_t k = static_cast<std::size_t>(std::ceil(p * static_cast<double>(s.size())));
    if (k < 1) k = 1;
    return s[k - 1];
}
double countBelow(const std::vector<double>& s, double x) {
    return static_cast<double>(std::lower_bound(s.begin(), s.end(), x) - s.begin());
}
double countAtOrBelow(const std::vector<double>& s, double x) {
    return static_cast<double>(std::upper_bound(s.begin(), s.end(), x) - s.begin());
}

/// Binned KS + chi2 of `s` (already transformed to the histogram variable) against a Fortran histogram.
int compareHist(const std::string& label, const Hist& h, const std::vector<double>& s, double nf, double& pKS,
                double& pChi2) {
    const double n = static_cast<double>(s.size());
    const std::size_t nb = h.cnt.size();
    double cumF = h.under, D = 0, chi2 = 0, prevC = 0;
    int ndf = -1;
    const double scale = nf / n;
    auto addChi = [&](double a, double bb) {
        if (a + bb > 0) {
            chi2 += (a - bb) * (a - bb) / (a + bb);
            ++ndf;
        }
    };
    for (std::size_t i = 0; i <= nb; ++i) {
        const double cC = i < nb ? countBelow(s, h.lo[i]) : countAtOrBelow(s, h.hi[nb - 1]);
        if (i == 0) {
            addChi(h.under, cC * scale);
        } else {
            addChi(h.cnt[i - 1], (cC - prevC) * scale);
            cumF += h.cnt[i - 1];
        }
        prevC = cC;
        D = std::max(D, std::fabs(cumF / nf - cC / n));
    }
    addChi(h.over, (n - prevC) * scale);
    const double ne = n * nf / (n + nf);
    pKS = TMath::KolmogorovProb(std::sqrt(ne) * D);
    pChi2 = ndf > 0 ? TMath::Prob(chi2, ndf) : 1.0;
    std::printf("    %-22s KS D=%.3e p=%.4f | chi2=%.1f/%d p=%.4f\n", label.c_str(), D, pKS, chi2, ndf, pChi2);
    return (pKS < 1e-3 || pChi2 < 1e-3) ? 1 : 0;
}

int compareDistribution(const std::string& dir, int cfgid, const std::string& var, std::vector<double> s,
                        const std::map<std::string, double>& fsum) {
    std::sort(s.begin(), s.end());
    const double n = static_cast<double>(s.size());
    const double nf = fsum.at("nsel");
    int fail = 0;
    std::printf("  cfg %d %s: n_cpp=%.0f n_f=%.0f\n", cfgid, var.c_str(), n, nf);
    std::printf("    %-5s %22s %22s %10s %8s\n", "q", "Fortran", "C++", "rel diff", "z");
    const double Q[] = {0.001, 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99, 0.999};
    const char* QN[] = {"q0p1", "q1", "q5", "q25", "q50", "q75", "q95", "q99", "q99p9"};
    double zmax = 0;
    for (int i = 0; i < 9; ++i) {
        const double qf = fsum.at(QN[i]);
        const double qc = quantileInv(s, Q[i]);
        // ECDF of the C++ sample at the Fortran quantile. Atoms (e.g. gen_weight: one value per tail at fixed
        // kinematics) make the ECDF jump at qf: z = 0 if p lies inside the jump [F(qf-), F(qf)].
        const double Flo = countBelow(s, qf) / n, Fhi = countAtOrBelow(s, qf) / n;
        const double sig = std::sqrt(Q[i] * (1 - Q[i]) * (1. / n + 1. / nf));
        const double z = Q[i] < Flo ? (Flo - Q[i]) / sig : (Q[i] > Fhi ? (Fhi - Q[i]) / sig : 0.0);
        zmax = std::max(zmax, std::fabs(z));
        if (std::fabs(z) > 4.) ++fail;
        std::printf("    %-5s %22.14e %22.14e %10.2e %8.3f\n", QN[i], qf, qc, qf != 0 ? (qc - qf) / std::fabs(qf) : 0.,
                    z);
    }
    double pKS, pChi2;
    char hp[128];
    std::snprintf(hp, sizeof hp, "/hist_%02d_%s.txt", cfgid, var.c_str());
    Hist h;
    if (!loadHist(dir + hp, h)) {
        std::printf("    missing %s\n", hp);
        return 1;
    }
    fail += compareHist("linear histogram", h, s, nf, pKS, pChi2);
    std::snprintf(hp, sizeof hp, "/hist_%02d_%s_log.txt", cfgid, var.c_str());
    Hist hl;
    if (loadHist(dir + hp, hl)) {
        std::vector<double> ls;
        for (double x : s)
            if (x > 0) ls.push_back(std::log10(x));
        double nfl = hl.under + hl.over;
        for (double c : hl.cnt) nfl += c;
        fail += compareHist("log10 histogram", hl, ls, nfl, pKS, pChi2);
    }
    double m = 0, m2 = 0;
    for (double x : s) m += x;
    m /= n;
    for (double x : s) m2 += (x - m) * (x - m);
    std::printf("    info: mean F %.6e C %.6e | rms F %.6e C %.6e | min F %.6e C %.6e | max F %.6e C %.6e | max|z_q| %.2f\n",
                fsum.at("mean"), m, fsum.at("rms"), std::sqrt(m2 / n), fsum.at("min"), s.front(), fsum.at("max"),
                s.back(), zmax);
    std::printf("    => %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

int runRadDistributions(const std::string& dir, int cfgid, const std::string& radc_state_dir,
                        const std::string& fortran_dir) {
    simcl1::Table cfg;
    std::string err;
    if (!cfg.load(dir + "/config.txt", err) || !std::ifstream(dir + "/summary.txt")) {
        std::printf("rad_distributions: SKIP (%s)\n", err.c_str());
        return 77;
    }
    int irow = -1;
    for (std::size_t i = 0; i < cfg.rows.size(); ++i)
        if (std::lround(cfg.rows[i][0]) == cfgid) irow = static_cast<int>(i);
    if (irow < 0) {
        std::printf("rad_distributions: cfgid %d not found\n", cfgid);
        return 1;
    }
    Row c(cfg, static_cast<std::size_t>(irow));
    const int kind = iv(c, "kind"), nsamp = iv(c, "nsamp");
    const long seed = std::lround(c["seed"]);
    // summary rows of this cfgid: var -> columns
    std::map<std::string, std::map<std::string, double>> fs;
    std::ifstream fsf(dir + "/summary.txt");
    std::string line;
    std::vector<std::string> names;
    while (std::getline(fsf, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line);
        std::vector<std::string> w;
        std::string x;
        while (is >> x) w.push_back(x);
        if (names.empty()) {
            names = w;
            continue;
        }
        if (std::stol(w[0]) != cfgid) continue;
        std::map<std::string, double> m;
        for (std::size_t k = 0; k < names.size(); ++k)
            if (names[k] != "var") m[names[k]] = std::strtod(w[k].c_str(), nullptr);
        fs[w[3]] = m;
    }

    auto P = makeProcess(cfgInp(radc_state_dir, 0), fortran_dir);
    if (!P->error.empty()) {
        std::printf("rad_distributions: %s\n", P->error.c_str());
        return 1;
    }
    Mt64Random rng(static_cast<std::uint64_t>(seed));
    std::printf("L1 rad_distributions cfgid %d kind %d: Mt64Random seed %ld, %d samples\n", cfgid, kind, seed, nsamp);
    int fail = 0;
    const RunSetup& s = P->setup;
    const CompleteEvFn complete_ev = P->reaction->completeEvFn();
    // Central pre-state (drv_raddist.f): trip_thru_target(1,0,Ebeam,0,...,Me,4), vertex at the central angles,
    // complete_ev (gauss1 from the real stream), main%Trec.
    MainEvent m0;
    Event v0;
    double el = 0, tf = 0;
    P->target.trip(1, 0.0, s.Ebeam, 0.0, el, tf, constants::Me, 4, rng);
    m0.target.Eloss[0] = el;
    m0.target.teff[0] = tf;
    m0.gen_weight = 1.0;
    v0.Ein = s.Ebeam_vertex_ave;
    v0.e.theta = s.spec.e.theta;
    v0.e.phi = s.spec.e.phi;
    if (!complete_ev(m0, v0, rng)) {
        std::printf("  complete_ev failed for the central pre-state\n");
        return 1;
    }
    m0.Trec = v0.Trec;
    m0.sigcc = 1.0;
    const RadcState st0 = P->rad.state();
    std::printf("  pre: Ein teff1 Eloss1 e_E p_E  %.16e %.16e %.16e %.16e %.16e\n", v0.Ein, tf, el, v0.e.E, v0.p.E);
    std::printf("  pre: frac1..3 g4 hardcorfac   %.16e %.16e %.16e %.16e %.16e\n", st0.frac[1], st0.frac[2], st0.frac[3],
                st0.g[4], st0.hardcorfac);
    {   // the pre-state printed by drv_raddist (driver.log) must be reproduced exactly
        const std::string log = readFile(dir + "/driver.log");
        const auto p1 = log.find("pre: Ein teff1 Eloss1 e_E p_E");
        const auto p2 = log.find("pre: frac1..3 g4 hardcorfac");
        if (p1 != std::string::npos && p2 != std::string::npos) {
            double a[5], bb[5];
            std::sscanf(log.c_str() + p1 + 29, "%lf %lf %lf %lf %lf", &a[0], &a[1], &a[2], &a[3], &a[4]);
            std::sscanf(log.c_str() + p2 + 28, "%lf %lf %lf %lf %lf", &bb[0], &bb[1], &bb[2], &bb[3], &bb[4]);
            const double mine[10] = {v0.Ein, tf, el, v0.e.E, v0.p.E, st0.frac[1], st0.frac[2], st0.frac[3], st0.g[4],
                                     st0.hardcorfac};
            int nb = 0;
            for (int i = 0; i < 10; ++i)
                if (mine[i] != (i < 5 ? a[i] : bb[i - 5])) ++nb;
            std::printf("  pre-state vs driver.log: %d of 10 values differ -> %s\n", nb, nb ? "FAIL" : "OK");
            fail += nb ? 1 : 0;
        }
    }

    if (kind == 1) {
        const int itail = iv(c, "itail");
        std::vector<double> eg(static_cast<std::size_t>(nsamp));
        for (int i = 0; i < nsamp; ++i) {
            int it = itail;
            double e = 0, w = 0, vr = 0;
            basicrad(st0, it, c["Egamma_lo"], c["Egamma_hi"], e, w, vr, rng);
            eg[static_cast<std::size_t>(i)] = e;
        }
        fail += compareDistribution(dir, cfgid, "Egamma", eg, fs.at("Egamma"));
    } else {
        std::vector<double> e1, e2, e3, rp, gw;
        long nt[4] = {0, 0, 0, 0}, nsucc = 0;
        for (int i = 0; i < nsamp; ++i) {
            MainEvent m = m0;
            Event v = v0, orig;
            P->rad.state() = st0;
            const bool ok = P->rad.radiate(m, v, orig, rng, complete_ev);
            const RadcState& st = P->rad.state();
            nsucc += ok;
            ++nt[st.ntail];
            if (ok && st.ntail == 1) e1.push_back(st.Egamma_used[1]);
            if (ok && st.ntail == 2) e2.push_back(st.Egamma_used[2]);
            if (ok && st.ntail == 3) e3.push_back(st.Egamma_used[3]);
            rp.push_back(m.ntup.radphot);
            gw.push_back(m.gen_weight);
        }
        // fractions (fractions.txt: Fortran counts) -> binomial z
        const std::string fr = readFile(dir + "/fractions.txt");
        std::printf("  success fraction C++ %.10f (%ld)\n", nsucc / double(nsamp), nsucc);
        for (int k = 1; k <= 3; ++k) {
            const std::string key = "ntail=" + std::to_string(k) + ": fraction of calls ";
            const auto p = fr.find(key);
            double ff = 0;
            if (p != std::string::npos) std::sscanf(fr.c_str() + p + key.size(), "%lf", &ff);
            const double fc = nt[k] / double(nsamp);
            const double z = (fc - ff) / std::sqrt(ff * (1 - ff) * 2. / nsamp);
            std::printf("  ntail=%d fraction: Fortran %.6f C++ %.6f z=%.2f\n", k, ff, fc, z);
            if (std::fabs(z) > 4) ++fail;
        }
        fail += compareDistribution(dir, cfgid, "Egamma1_ntail1", e1, fs.at("Egamma1_ntail1"));
        fail += compareDistribution(dir, cfgid, "Egamma2_ntail2", e2, fs.at("Egamma2_ntail2"));
        fail += compareDistribution(dir, cfgid, "Egamma3_ntail3", e3, fs.at("Egamma3_ntail3"));
        fail += compareDistribution(dir, cfgid, "radphot", rp, fs.at("radphot"));
        fail += compareDistribution(dir, cfgid, "gen_weight", gw, fs.at("gen_weight"));
    }
    std::printf("  => %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

// ---------------------------------------------------------------------------------------------------------------
int runSelftest() {
    int fail = 0;
    if (!Registry<IRadiator>::contains("simc_rad")) {
        std::printf("IRadiator 'simc_rad' not registered\n");
        ++fail;
    }
    auto r = Registry<IRadiator>::create("simc_rad");
    if (r->name() != "simc_rad") ++fail;
    // FortranStop paths
    RadcRun run;
    RadcState st;
    run.extrad_flag = 3;
    try {
        extradPhi(run, st, 0, 1., 1., 1.);
        ++fail;
        std::printf("extrad_phi(0) with extrad_flag 3 did not STOP\n");
    } catch (const FortranStop&) {
    }
    try {
        radcInit(1, 5.31, 6.144, 0, -1, 0, false);
        ++fail;
        std::printf("radc_init with extrad_flag < 0 did not STOP\n");
    } catch (const FortranStop&) {
    }
    std::printf("radiation selftest: %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mode> <dir> [tol] [...]\n", argv[0]);
        return 1;
    }
    const std::string what = argv[1];
    if (what == "selftest") return runSelftest();
    if (argc < 3) return 1;
    const std::string dir = argv[2];
    const double tol = argc > 3 ? std::atof(argv[3]) : 1e-12;
    try {
        if (what == "gamma") return runGamma(dir, tol);
        if (what == "spen") return runSpen(dir, tol);
        if (what == "spence") return runSpence(dir, tol);
        if (what == "schwinger") return runSchwinger(dir, tol);
        if (what == "lambda_dave") return runLambdaDave(dir, tol);
        if (what == "brem") return runBrem(dir, tol);
        if (what == "bremos") return runBremos(dir, tol);
        if (what == "inter") return runInter(dir, tol);
        if (what == "inter_prime") return runInterPrime(dir, tol);
        if (what == "srad") return runSrad(dir, tol);
        if (what == "extrad_friedrich") return runExtradFriedrich(dir, tol);
        if (what == "extrad_phi") return runExtradPhi(dir, tol);
        if (what == "radc_init") return runRadcInit(dir, tol);
        if (what == "radc_init_ev") return runRadcInitEv(dir, tol);
        if (what == "basicrad_init_ev") return runBasicradInitEv(dir, tol);
        if (what == "peaked_rad_weight") return runPeakedRadWeight(dir, tol, false);
        if (what == "peaked_rad_weight_r1x3") return runPeakedRadWeight(dir, tol, true);
        if (what == "basicrad_fixedrng") return runBasicradFixed(dir, tol);
        if (what == "generate_rad_fixedrng" && argc > 5) return runGenerateRadFixed(dir, tol, argv[4], argv[5]);
        if (what == "generate_rad_capture" && argc > 5) return runGenerateRadCapture(dir, tol, argv[4], argv[5]);
        if (what == "radc_state" && argc > 4) return runRadcState(dir, tol, argv[4]);
        if (what == "rad_distributions" && argc > 5)
            return runRadDistributions(dir, std::atoi(argv[3]), argv[4], argv[5]);
    } catch (const std::exception& e) {
        std::printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "unknown mode %s (or missing arguments)\n", what.c_str());
    return 1;
}
