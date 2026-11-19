// tests/l1/test_l1_reactions.cpp
//
// L1 tests of the event generation (task P4.1) against the P4.F1 Fortran tables and captures
// (reference/functions/{generate_fixedrng, complete_ev, complete_recon_ev, complete_main, transform_to_cm},
// reference/functions/event_capture/<cfg>/, reference/P4.F1_report.md).
//
//   table <name> <dir> <tol> <runs dir> <fortran dir>
//       generate_fixedrng / complete_ev / complete_recon_ev / complete_main: every output column of the table,
//       one "process" (setup + target + radiator + beam + reaction) per cfg, rows in file order;
//       complete_main also checks the product identity weight = ((SF_weight*jacobian)*gen_weight)*sigcc.
//   transform_to_cm <dir> <tol>
//   capture <capture dir> <inp> <tol> <fortran dir>
//       replay of the real calls of one capture run: generate (ScriptedRandom fed with the captured phase-1 and
//       phase-6 grnd/gauss1 values of every try, persistent main/vertex/orig/radiator state as in SIMC),
//       complete_ev (first call of every try), complete_recon_ev and complete_main (every call).
//   selftest
// Exit: 0 pass, 1 fail, 77 reference missing.
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/input/InpReader.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/radiation/SimcRadiator.h"
#include "simc/radiation/kernel/RadcState.h"
#include "simc/reactions/EventChain.h"
#include "simc/reactions/HydrogenElastic.h"
#include "simc/reactions/StandardBeam.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/setup/SetupBuilder.h"

using namespace simc;
using namespace simc::reactions;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

int iv(const Row& r, const char* c) { return static_cast<int>(std::lround(r[c])); }

// ---- output helpers (column names of the P4.F1 drivers) ----------------------------------------------------------
void putArm(const std::string& p, const ArmFull& a, Outputs& o) {
    o[p + "delta"] = a.delta;
    o[p + "xptar"] = a.xptar;
    o[p + "yptar"] = a.yptar;
    o[p + "z"] = a.z;
    o[p + "theta"] = a.theta;
    o[p + "phi"] = a.phi;
    o[p + "E"] = a.E;
    o[p + "P"] = a.P;
}
void putEvent(const std::string& p, const Event& v, Outputs& o) {
    o[p + "Ein"] = v.Ein;
    o[p + "Em"] = v.Em;
    o[p + "Pm"] = v.Pm;
    o[p + "Emiss"] = v.Emiss;
    o[p + "Pmiss"] = v.Pmiss;
    o[p + "Pmx"] = v.Pmx;
    o[p + "Pmy"] = v.Pmy;
    o[p + "Pmz"] = v.Pmz;
    o[p + "PmPar"] = v.PmPar;
    o[p + "PmPer"] = v.PmPer;
    o[p + "PmOop"] = v.PmOop;
    o[p + "nu"] = v.nu;
    o[p + "q"] = v.q;
    o[p + "Q2"] = v.Q2;
    o[p + "Trec"] = v.Trec;
    o[p + "W"] = v.W;
    o[p + "Mrec"] = v.Mrec;
    o[p + "epsilon"] = v.epsilon;
    o[p + "theta_pq"] = v.theta_pq;
    o[p + "theta_tarq"] = v.theta_tarq;
    o[p + "phi_pq"] = v.phi_pq;
    o[p + "phi_targ"] = v.phi_targ;
    o[p + "beta"] = v.beta;
    o[p + "phi_s"] = v.phi_s;
    o[p + "phi_c"] = v.phi_c;
    o[p + "zhad"] = v.zhad;
    o[p + "pt2"] = v.pt2;
    o[p + "xbj"] = v.xbj;
    putArm(p + "e_", v.e, o);
    putArm(p + "p_", v.p, o);
    o[p + "ue_x"] = v.ue.x;
    o[p + "ue_y"] = v.ue.y;
    o[p + "ue_z"] = v.ue.z;
    o[p + "up_x"] = v.up.x;
    o[p + "up_y"] = v.up.y;
    o[p + "up_z"] = v.up.z;
    o[p + "uq_x"] = v.uq.x;
    o[p + "uq_y"] = v.uq.y;
    o[p + "uq_z"] = v.uq.z;
}
void putMain(const MainEvent& m, Outputs& o) {
    o["m_weight"] = m.weight;
    o["m_SF_weight"] = m.SF_weight;
    o["m_gen_weight"] = m.gen_weight;
    o["m_jacobian"] = m.jacobian;
    o["m_Ein_shift"] = m.Ein_shift;
    o["m_Ee_shift"] = m.Ee_shift;
    o["m_sigcc"] = m.sigcc;
    o["m_sigcc_recon"] = m.sigcc_recon;
    o["m_sigcent"] = m.sigcent;
    o["m_epsilon"] = m.epsilon;
    o["m_theta_pq"] = m.theta_pq;
    o["m_theta_tarq"] = m.theta_tarq;
    o["m_phi_pq"] = m.phi_pq;
    o["m_phi_targ"] = m.phi_targ;
    o["m_beta"] = m.beta;
    o["m_w"] = m.w;
    o["m_t"] = m.t;
    o["m_tmin"] = m.tmin;
    o["m_q2"] = m.q2;
    o["m_pcm"] = m.pcm;
    o["m_thetacm"] = m.thetacm;
    o["m_phicm"] = m.phicm;
    o["m_wcm"] = m.wcm;
    o["m_davejac"] = m.davejac;
    o["m_johnjac"] = m.johnjac;
    o["m_Trec"] = m.Trec;
}
void putTarget(const EventTarget& t, Outputs& o) {
    o["t_x"] = t.x;
    o["t_y"] = t.y;
    o["t_z"] = t.z;
    o["t_rasterx"] = t.rasterx;
    o["t_rastery"] = t.rastery;
    for (int i = 0; i < 3; ++i) {
        o["t_teff_" + std::to_string(i + 1)] = t.teff[i];
        o["t_Eloss_" + std::to_string(i + 1)] = t.Eloss[i];
    }
    o["t_Coulomb"] = t.Coulomb;
}
void putRad(const radiation::RadcState& st, Outputs& o, const std::string& p = "r_") {
    o[p + "ntail"] = st.ntail;
    for (int i = 1; i <= 3; ++i) o[p + "Egamma_used_" + std::to_string(i)] = st.Egamma_used[i];
    o[p + "hardcorfac"] = st.hardcorfac;
    for (int i = 1; i <= 3; ++i) o[p + "frac_" + std::to_string(i)] = st.frac[i];
    o[p + "g_4"] = st.g[4];
    o[p + "c_4"] = st.c[4];
    o[p + "rad_proton_ev"] = st.rad_proton_this_ev;
}

// ---- one Fortran process (es_init) --------------------------------------------------------------------------------
/// IRadiator decorator: records whether generate_rad was called (the table's rad_called / fail_stage).
struct RecordingRadiator : IRadiator {
    radiation::SimcRadiator inner;
    bool called{false}, result{false};
    std::string name() const override { return "recording"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext& c) override { inner.init(c); }
    void initEvent(MainEvent& m, const Event& v) override { inner.initEvent(m, v); }
    bool radiate(MainEvent& m, Event& v, Event& o, IRandom& r, const CompleteEvFn& c) override {
        called = true;
        result = inner.radiate(m, v, o, r, c);
        return result;
    }
};

struct Process {
    SimcParameters params;
    RunSetup setup;
    material::HallCCryoTarget target;
    RecordingRadiator rad;
    StandardBeam beam;
    std::unique_ptr<HydrogenElastic> reaction;
    std::unique_ptr<IPhysicsModel> model;
    MainContext ctx;
    std::string error;
};

std::unique_ptr<Process> makeProcess(const std::string& inp, const std::string& fortran_dir) {
    auto P = std::make_unique<Process>();
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
    Json tcfg = Json::object();
    tcfg["Mh"] = P->setup.Mh;
    P->target.setLog(nullptr);
    P->target.configure(tcfg, P->params);
    P->target.setLog(nullptr);
    P->target.init(P->setup);
    P->rad.inner.setLog(nullptr);
    P->rad.init(P->setup);
    StageContext sc(P->setup);
    sc.target = &P->target;
    sc.radiator = &P->rad;
    P->beam.configure(Json::object(), P->params);
    P->beam.init(sc);
    P->model = Registry<IPhysicsModel>::create("sigep");
    P->model->configure(Json::object(), P->params);
    P->model->init(P->setup);
    sc.physics = P->model.get();
    P->reaction = std::make_unique<HydrogenElastic>(P->setup.flags.doing_phsp ? "phase_space" : "hydrogen_elastic");
    P->reaction->configure(Json::object(), P->params);
    P->reaction->setLog(nullptr);
    P->reaction->init(sc);
    P->ctx = MainContext::FromSetup(P->setup);
    return P;
}

std::string cfgInp(const std::string& runs, int cfg) {
    static const char* stems[] = {"", "eep_h_elastic_norad_A", "eep_h_elastic_rad_A", "acc_phsp_hms_shms_A"};
    if (cfg < 1 || cfg > 3) return "";
    return runs + "/" + stems[cfg] + "/" + stems[cfg] + ".inp";
}

/// One generate call as the SIMC main loop does it (simc.f:174-183 resets; event.f:generate).
bool runGenerate(Process& P, MainEvent& m, Event& v, Event& o, ScriptedRandom& rng, Outputs& out) {
    m.ntup.radphot = 0.0;   // simc.f:175-183
    m.ntup.radarm = 0.0;
    m.ntup.resfac = 0.0;
    m.ntup.sigcm = 0.0;
    m.ntup.krel = 0.0;
    m.ntup.mm = 0.0;
    m.ntup.mmA = 0.0;
    m.ntup.t = 0.0;
    m.ntup.decdist = 0.0;
    P.rad.called = false;
    const bool ok = generateEvent(P.beam, *P.reaction, m, v, o, rng);
    out["success"] = ok;
    out["fail_stage"] = ok ? 0 : (P.rad.called ? 2 : 1);
    out["rad_called"] = P.rad.called;
    out["n_unif"] = static_cast<double>(rng.nUniform());
    out["n_gauss"] = static_cast<double>(rng.nGauss());
    out["rng_err"] = rng.error();
    out["radphot"] = m.ntup.radphot;
    out["radarm"] = m.ntup.radarm;
    putMain(m, out);
    putTarget(m.target, out);
    putEvent("v_", v, out);
    putEvent("o_", o, out);
    putRad(P.rad.inner.state(), out);
    return ok;
}

/// Run a table: every column of output.txt that is not an input column is compared.
int runAuto(const std::string& title, const std::string& dir, double tol, const simcl1::Adapter& adapter) {
    simcl1::Table in, ref;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !ref.load(dir + "/output.txt", err)) {
        std::printf("L1 %s: SKIP (%s)\n", title.c_str(), err.c_str());
        return 77;
    }
    std::set<std::string> inNames(in.names.begin(), in.names.end());
    std::vector<ColumnSpec> cols;
    for (const auto& n : ref.names)
        if (!inNames.count(n)) cols.push_back({n, tol, 0.0});
    auto res = simcl1::compare(in, ref, cols, adapter, err);
    if (!err.empty()) {
        std::printf("L1 %s: ERROR %s\n", title.c_str(), err.c_str());
        return 1;
    }
    std::size_t rows = res.empty() ? 0 : res[0].rows, nb = 0, nc = 0;
    for (const auto& r : res) {
        nb += r.bitwise_equal;
        nc += r.rows;
    }
    const bool ok = simcl1::report(title, res, cols);
    std::printf("  summary: %zu rows x %zu columns, %zu of %zu values bitwise equal\n", rows, cols.size(), nb, nc);
    return ok ? 0 : 1;
}

// ---- tables -------------------------------------------------------------------------------------------------------
int runTable(const std::string& name, const std::string& dir, double tol, const std::string& runs,
             const std::string& fortran_dir) {
    std::map<int, std::unique_ptr<Process>> procs;
    auto proc = [&](int cfg) -> Process& {
        auto& P = procs[cfg];
        if (!P) {
            P = makeProcess(cfgInp(runs, cfg), fortran_dir);
            if (!P->error.empty()) throw std::runtime_error("cfg " + std::to_string(cfg) + ": " + P->error);
        }
        return *P;
    };
    std::size_t identity_bad = 0, identity_rows = 0;
    int rc = 1;
    if (name == "generate_fixedrng") {
        rc = runAuto("generate_fixedrng (event.f:126-428, scripted grnd/gauss1)", dir, tol,
                     [&](const Row& r, Outputs& o) {
                         Process& P = proc(iv(r, "cfg"));
                         MainEvent m;
                         Event v, orig;
                         for (int i = 1; i <= 3; ++i) P.rad.inner.state().Egamma_min[i] = P.rad.inner.state().Egamma_max[i] = 0.;
                         std::vector<double> u, g;
                         for (int i = 1; i <= 10; ++i) u.push_back(r["u" + std::to_string(i)]);
                         for (int i = 1; i <= 24; ++i) g.push_back(r["g" + std::to_string(i)]);
                         ScriptedRandom rng;
                         rng.load(g, u);
                         runGenerate(P, m, v, orig, rng, o);
                     });
    } else if (name == "complete_ev") {
        rc = runAuto("complete_ev (event.f:432-1052)", dir, tol, [&](const Row& r, Outputs& o) {
            Process& P = proc(iv(r, "cfg"));
            const RunSetup& s = P.setup;
            MainEvent m;
            Event v;
            v.Ein = r["Ein"];
            v.e.xptar = r["e_xptar"];
            v.e.yptar = r["e_yptar"];
            v.e.theta = r["e_theta"];
            v.e.phi = r["e_phi"];
            if (iv(r, "calc_angles") == 1)
                PhysicsAngles(s.spec.e.theta, s.spec.e.phi, v.e.xptar, v.e.yptar, v.e.theta, v.e.phi);
            o["in_e_theta"] = v.e.theta;
            o["in_e_phi"] = v.e.phi;
            m.target.z = r["t_z"];
            m.target.teff[0] = r["t_teff_1"];
            std::vector<double> g;
            for (int i = 1; i <= 12; ++i) g.push_back(r["g" + std::to_string(i)]);
            ScriptedRandom rng;
            rng.load(g, {});
            const bool ok = P.reaction->completeEv(m, v, rng);
            o["success"] = ok;
            o["n_gauss"] = static_cast<double>(rng.nGauss());
            o["n_unif"] = static_cast<double>(rng.nUniform());
            o["rng_err"] = rng.error();
            putEvent("v_", v, o);
            putMain(m, o);
            o["t_teff_2"] = m.target.teff[1];
            o["t_teff_3"] = m.target.teff[2];
            o["t_Eloss_2"] = m.target.Eloss[1];
            o["t_Eloss_3"] = m.target.Eloss[2];
            o["ntup_krel"] = m.ntup.krel;
            putRad(P.rad.inner.state(), o);
        });
    } else if (name == "complete_recon_ev") {
        rc = runAuto("complete_recon_ev (event.f:1056-1359)", dir, tol, [&](const Row& r, Outputs& o) {
            Process& P = proc(iv(r, "cfg"));
            Event rec;
            for (const char* a : {"e", "p"}) {
                ArmFull& x = (a[0] == 'e') ? rec.e : rec.p;
                const std::string p = std::string("in_") + a + "_";
                x.delta = r[p + "delta"];
                x.yptar = r[p + "yptar"];
                x.xptar = r[p + "xptar"];
                x.z = r[p + "z"];
                x.theta = r[p + "theta"];
                x.phi = r[p + "phi"];
                x.E = r[p + "E"];
                x.P = r[p + "P"];
            }
            o["success"] = P.reaction->completeRecon(rec);
            putEvent("r_", rec, o);
            o["ntup_mm"] = P.reaction->reconNtup().mm;
            o["ntup_mmA"] = P.reaction->reconNtup().mmA;
            o["ntup_t"] = P.reaction->reconNtup().t;
        });
    } else if (name == "complete_main") {
        rc = runAuto("complete_main (event.f:1363-1569)", dir, tol, [&](const Row& r, Outputs& o) {
            Process& P = proc(iv(r, "cfg"));
            MainEvent m;
            Event v, rec, v0;
            v.Ein = r["v_Ein"];
            v.e.E = r["v_e_E"];
            v.e.theta = r["v_e_theta"];
            v.Q2 = r["v_Q2"];
            rec.Ein = r["r_Ein"];
            rec.e.E = r["r_e_E"];
            rec.e.theta = r["r_e_theta"];
            rec.Q2 = r["r_Q2"];
            m.jacobian = r["in_jacobian"];
            m.gen_weight = r["in_gen_weight"];
            MainContext ctx = P.ctx;
            if (iv(r, "ov_coulomb") == 1) {   // harness override (not a SIMC state)
                ctx.using_Coulomb = true;
                ctx.Coulomb_ave = r["ov_coulomb_ave"];
            }
            P.reaction->mainContext() = ctx;   // the harness override of the COMMON state for this row
            const bool ok = P.reaction->completeMain(iv(r, "force_sigcc") == 1, m, v, v0, rec);
            P.reaction->mainContext() = P.ctx;
            o["success"] = ok;
            o["SF_weight"] = m.SF_weight;
            o["sigcc"] = m.sigcc;
            o["sigcc_recon"] = m.sigcc_recon;
            o["jacobian"] = m.jacobian;
            o["gen_weight"] = m.gen_weight;
            o["weight"] = m.weight;
            const double c = 1.0 + ctx.Coulomb_ave / ctx.Ebeam;
            o["hx_coul_fac"] = ctx.using_Coulomb ? c * c : 1.0;
            o["hx_prod"] = m.SF_weight * m.jacobian * m.gen_weight * m.sigcc;
            ++identity_rows;
            if (std::memcmp(&o["hx_prod"], &m.weight, sizeof(double)) != 0) ++identity_bad;
            putMain(m, o);
        });
        std::printf("  product identity weight == ((SF_weight*jacobian)*gen_weight)*sigcc: %zu of %zu rows violate\n",
                    identity_bad, identity_rows);
        if (identity_bad) rc = 1;
    } else {
        std::printf("unknown table %s\n", name.c_str());
    }
    return rc;
}

int runTransformToCm(const std::string& dir, double tol) {
    return runAuto("transform_to_cm (jacobians.f:1-282)", dir, tol, [&](const Row& r, Outputs& o) {
        Event v;
        MainEvent m;
        v.Ein = r["Ein"];
        v.nu = r["nu"];
        v.q = r["q"];
        v.p.P = r["p_P"];
        v.p.E = r["p_E"];
        v.uq = {r["uq_x"], r["uq_y"], r["uq_z"]};
        v.up = {r["up_x"], r["up_y"], r["up_z"]};
        m.phi_pq = r["phi_pq"];
        const CmTransform c = transformToCm(v, m, r["pfer"], r["pferx"], r["pfery"], r["pferz"], r["efer"]);
        o["gstar"] = c.gstar; o["bstar"] = c.bstar; o["bstarx"] = c.bstarx; o["bstary"] = c.bstary; o["bstarz"] = c.bstarz;
        o["nustar"] = c.nustar; o["qstar"] = c.qstar; o["qstarx"] = c.qstarx; o["qstary"] = c.qstary; o["qstarz"] = c.qstarz;
        o["ehadcm"] = c.ehadcm; o["phadcm"] = c.phadcm; o["phadcmx"] = c.phadcmx; o["phadcmy"] = c.phadcmy;
        o["phadcmz"] = c.phadcmz;
        o["ebeamcm"] = c.ebeamcm; o["pbeamcm"] = c.pbeamcm; o["pbeamcmx"] = c.pbeamcmx; o["pbeamcmy"] = c.pbeamcmy;
        o["pbeamcmz"] = c.pbeamcmz;
        o["etarcm"] = c.etarcm; o["ptarcm"] = c.ptarcm; o["ptarcmx"] = c.ptarcmx; o["ptarcmy"] = c.ptarcmy;
        o["ptarcmz"] = c.ptarcmz;
        o["thetacm"] = c.thetacm; o["phicm"] = c.phicm; o["phiqn"] = c.phiqn; o["jacobian"] = c.jacobian;
        o["jac_old"] = c.jac_old;
    });
}

// ---- captures -----------------------------------------------------------------------------------------------------
/// Load a whitespace table, plain or .xz (through `xz -dc`).
bool loadAny(const std::string& base, simcl1::Table& t, std::string& err) {
    std::string path = base;
    bool xz = false;
    if (!std::ifstream(path)) {
        path = base + ".xz";
        xz = true;
        if (!std::ifstream(path)) {
            err = "missing " + base + "[.xz]";
            return false;
        }
    }
    t.names.clear();
    t.rows.clear();
    t.path = path;
    FILE* f = xz ? popen(("xz -dc '" + path + "'").c_str(), "r") : std::fopen(path.c_str(), "r");
    if (!f) {
        err = "cannot open " + path;
        return false;
    }
    std::string line;
    char buf[65536];
    while (std::fgets(buf, sizeof buf, f)) {
        line = buf;
        const auto p = line.find_first_not_of(" \t\r\n");
        if (p == std::string::npos || line[p] == '#') continue;
        std::istringstream ss(line);
        std::string w;
        if (t.names.empty()) {
            while (ss >> w) t.names.push_back(w);
            continue;
        }
        std::vector<double> r;
        r.reserve(t.names.size());
        while (ss >> w) r.push_back(std::strtod(w.c_str(), nullptr));
        if (r.size() != t.names.size()) {
            err = path + ": bad row";
            xz ? pclose(f) : std::fclose(f);
            return false;
        }
        t.rows.push_back(std::move(r));
    }
    xz ? pclose(f) : std::fclose(f);
    return true;
}

struct CaptureCompare {
    CaptureCompare(std::string w, double t) : what(std::move(w)), tol(t) {}
    std::string what;
    double tol;
    std::size_t rows{0}, rowsBad{0}, values{0}, bitwise{0};
    double maxrel{0};
    std::string worstCol;
    std::map<std::string, std::size_t> badCols;
    void row(const Row& r, const simcl1::Table& t, const Outputs& o, const std::set<std::string>& exclude) {
        ++rows;
        bool any = false;
        for (const auto& n : t.names) {
            if (exclude.count(n)) continue;
            auto it = o.find(n);
            if (it == o.end()) continue;
            const double f = r[n], c = it->second;
            ++values;
            if (std::memcmp(&f, &c, sizeof f) == 0 || (std::isnan(f) && std::isnan(c)) || simcl1::u1UpstreamChange(n, f, c)) {
                ++bitwise;
                continue;
            }
            const double rel = simcl1::relDiff(c, f);
            if (rel > maxrel) {
                maxrel = rel;
                worstCol = n;
            }
            if (!(rel <= tol)) {
                ++badCols[n];
                any = true;
            }
        }
        if (any) {
            ++rowsBad;
            if (rowsBad <= 3) std::printf("    %s: mismatch at row %zu (ntried %.0f)\n", what.c_str(), r.number(), r["ntried"]);
        }
    }
    bool report() const {
        std::printf("  %-18s rows %6zu, values %9zu, bitwise %9zu, rows with a mismatch %zu, max|rel| %.3e %s\n",
                    what.c_str(), rows, values, bitwise, rowsBad, maxrel, worstCol.c_str());
        for (const auto& kv : badCols) std::printf("      column %s: %zu\n", kv.first.c_str(), kv.second);
        return rowsBad == 0 && rows > 0;
    }
};

int runCapture(const std::string& capdir, const std::string& inp, double tol, const std::string& fortran_dir) {
    simcl1::Table rngT, genT, cevT, recT, mainT;
    std::string err;
    if (!loadAny(capdir + "/cap_rng.txt", rngT, err) || !loadAny(capdir + "/cap_generate.txt", genT, err) ||
        !loadAny(capdir + "/cap_complete_ev.txt", cevT, err) || !loadAny(capdir + "/cap_recon.txt", recT, err) ||
        !loadAny(capdir + "/cap_main.txt", mainT, err)) {
        std::printf("capture %s: SKIP (%s)\n", capdir.c_str(), err.c_str());
        return 77;
    }
    std::printf("L1 capture replay %s (setup %s)\n", capdir.c_str(), inp.c_str());
    // per try: the grnd and gauss1 values of phases 1 (generate) and 6 (generate_rad), in order
    std::map<long, std::pair<std::vector<double>, std::vector<double>>> stream;   // ntried -> (u, g)
    {
        const int cN = rngT.col("ntried"), cP = rngT.col("phase"), cK = rngT.col("kind"), cV = rngT.col("value");
        for (const auto& r : rngT.rows) {
            const int ph = static_cast<int>(r[cP]);
            if (ph != 1 && ph != 6) continue;
            auto& q = stream[std::lround(r[cN])];
            (static_cast<int>(r[cK]) == 1 ? q.first : q.second).push_back(r[cV]);
        }
    }
    int fail = 0;

    // generate: one persistent process, as the SIMC main loop
    {
        auto P = makeProcess(inp, fortran_dir);
        if (!P->error.empty()) throw std::runtime_error(P->error);
        MainEvent m;
        Event v, orig;
        CaptureCompare cc{"generate", tol};
        std::size_t leftover = 0;
        const std::set<std::string> excl = {"ntried", "nevent", "n_unif_rad", "n_gauss_rad", "n_int"};
        for (std::size_t i = 0; i < genT.rows.size(); ++i) {
            Row r(genT, i);
            const auto& q = stream[std::lround(r["ntried"])];
            ScriptedRandom rng;
            rng.load(q.second, q.first);
            Outputs o;
            runGenerate(*P, m, v, orig, rng, o);
            if (rng.uniformsLeft() || rng.gaussiansLeft() || rng.error()) ++leftover;
            // capture column names
            o["Ein_shift"] = m.Ein_shift;
            o["Ee_shift"] = m.Ee_shift;
            o["gen_weight"] = m.gen_weight;
            o["jacobian"] = m.jacobian;
            o["ntail"] = P->rad.inner.state().ntail;
            o["hardcorfac"] = P->rad.inner.state().hardcorfac;
            for (int k = 1; k <= 3; ++k) o["Egamma_used_" + std::to_string(k)] = P->rad.inner.state().Egamma_used[k];
            cc.row(r, genT, o, excl);
        }
        fail += cc.report() ? 0 : 1;
        std::printf("      tries whose captured random numbers were not consumed exactly: %zu\n", leftover);
        if (leftover) ++fail;
    }
    // complete_ev (first call of every try)
    {
        auto P = makeProcess(inp, fortran_dir);
        CaptureCompare cc{"complete_ev", tol};
        for (std::size_t i = 0; i < cevT.rows.size(); ++i) {
            Row r(cevT, i);
            MainEvent m;
            Event v;
            v.Ein = r["Ein"];
            v.e.xptar = r["e_xptar"];
            v.e.yptar = r["e_yptar"];
            v.e.theta = r["e_theta"];
            v.e.phi = r["e_phi"];
            m.target.z = r["t_z"];
            m.target.teff[0] = r["t_teff_1"];
            std::vector<double> g;
            for (int k = 1; k <= iv(r, "ng"); ++k) g.push_back(r["g" + std::to_string(k)]);
            ScriptedRandom rng;
            rng.load(g, {});
            Outputs o;
            o["success"] = P->reaction->completeEv(m, v, rng);
            putEvent("v_", v, o);
            o["jacobian"] = m.jacobian;
            o["t_teff_2"] = m.target.teff[1];
            o["t_teff_3"] = m.target.teff[2];
            o["t_Eloss_2"] = m.target.Eloss[1];
            o["t_Eloss_3"] = m.target.Eloss[2];
            putRad(P->rad.inner.state(), o);
            cc.row(r, cevT, o, {"ntried"});
        }
        fail += cc.report() ? 0 : 1;
    }
    // complete_recon_ev and complete_main (deterministic)
    {
        auto P = makeProcess(inp, fortran_dir);
        CaptureCompare cr{"complete_recon_ev", tol}, cm{"complete_main", tol};
        for (std::size_t i = 0; i < recT.rows.size(); ++i) {
            Row r(recT, i);
            Event rec;
            for (const char* a : {"e", "p"}) {
                ArmFull& x = (a[0] == 'e') ? rec.e : rec.p;
                const std::string p = std::string("in_") + a + "_";
                x.delta = r[p + "delta"]; x.yptar = r[p + "yptar"]; x.xptar = r[p + "xptar"]; x.z = r[p + "z"];
                x.theta = r[p + "theta"]; x.phi = r[p + "phi"]; x.E = r[p + "E"]; x.P = r[p + "P"];
            }
            Outputs o;
            o["success"] = P->reaction->completeRecon(rec);
            putEvent("r_", rec, o);
            cr.row(r, recT, o, {"ntried"});
        }
        std::size_t identity_bad = 0;
        for (std::size_t i = 0; i < mainT.rows.size(); ++i) {
            Row r(mainT, i);
            MainEvent m;
            Event v, rec, v0;
            v.Ein = r["v_Ein"]; v.e.E = r["v_e_E"]; v.e.theta = r["v_e_theta"]; v.Q2 = r["v_Q2"];
            rec.Ein = r["r_Ein"]; rec.e.E = r["r_e_E"]; rec.e.theta = r["r_e_theta"]; rec.Q2 = r["r_Q2"];
            m.jacobian = r["in_jacobian"];
            m.gen_weight = r["in_gen_weight"];
            m.SF_weight = r["in_SF_weight"];   // stale values before the call (overwritten)
            m.sigcc = r["in_sigcc"];
            m.weight = r["in_weight"];
            Outputs o;
            o["success"] = P->reaction->completeMain(iv(r, "force_sigcc") == 1, m, v, v0, rec);
            o["SF_weight"] = m.SF_weight; o["sigcc"] = m.sigcc; o["sigcc_recon"] = m.sigcc_recon;
            o["jacobian"] = m.jacobian; o["gen_weight"] = m.gen_weight; o["weight"] = m.weight;
            const double prod = m.SF_weight * m.jacobian * m.gen_weight * m.sigcc;
            if (std::memcmp(&prod, &m.weight, sizeof prod) != 0) ++identity_bad;
            cm.row(r, mainT, o, {"ntried"});
        }
        fail += cr.report() ? 0 : 1;
        fail += cm.report() ? 0 : 1;
        std::printf("      complete_main product identity violated in %zu of %zu rows\n", identity_bad, mainT.rows.size());
        if (identity_bad) ++fail;
    }
    std::printf("  => %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

int runSelftest() {
    int fail = 0;
    for (const char* n : {"hydrogen_elastic", "phase_space"})
        if (!Registry<IReaction>::contains(n)) {
            std::printf("IReaction '%s' not registered\n", n);
            ++fail;
        }
    if (!Registry<IBeam>::contains("standard")) ++fail;
    std::printf("reactions selftest: %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    const std::string what = argv[1];
    try {
        if (what == "selftest") return runSelftest();
        if (what == "table" && argc > 6) return runTable(argv[2], argv[3], std::atof(argv[4]), argv[5], argv[6]);
        if (what == "transform_to_cm" && argc > 3) return runTransformToCm(argv[2], std::atof(argv[3]));
        if (what == "capture" && argc > 5) return runCapture(argv[2], argv[3], std::atof(argv[4]), argv[5]);
    } catch (const std::exception& e) {
        std::printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "usage: %s table|transform_to_cm|capture|selftest ...\n", argv[0]);
    return 1;
}
