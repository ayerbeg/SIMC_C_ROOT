// tests/l1/test_l1_pion_event.cpp
//
// L1 tests of the exclusive pion reaction (task P9.2, reactions::PionProduction) against the P9.F1 Fortran tables
// (reference/functions/, reference/P9.F1_report.md). fortran_compat mode (D-031), bitwise (tolerance argument, 0 =
// bitwise). One "process" (set-up + stages) per configuration, rows in table order, as the Fortran drivers.
//
//   pion_state <dir> <inp dir> <data/spectral> <data/pion> <tol>   start-up state + calculate_central per cfg 1-15,
//                                                                  and the 8 STOP configurations (stops.txt)
//   generate <dir> <inp dir> <data/spectral> <data/pion> <tol>     event.f:generate (pion branches), scripted RNG
//   complete_ev <dir> ...                                          event.f:complete_ev (pion branches), scripted RNG
//   complete_recon_ev <dir> ...                                    event.f:complete_recon_ev (+ ntup%mm/mmA/t)
//   complete_main <dir> ...                                        event.f:complete_main (pion branch)
//   capture <cap dir> <inp> <data/spectral> <data/pion> <tol>      every real generate / complete_ev /
//                                                                  complete_recon_ev / complete_main call of a capture
// Kaon (task P10.2): the same modes with two extra arguments `kaon <data/kaon>` run reactions::KaonProduction and the
// "kaon" model against the P10.F1 tables (kaon_state, generate_kaon_fixedrng, complete_ev_kaon, complete_recon_ev_kaon,
// complete_main_kaon, event_capture_kaon); kaon_state is mode pion_state with the kaon columns.
// Delta (task P11.1): `delta` as argument 7 (no data directory) runs reactions::DeltaProduction and the "delta" model
// against the P11.F1 tables (delta_state, generate_delta_fixedrng, complete_ev_delta, complete_recon_ev_delta,
// complete_main_delta, event_capture_delta) and adds the mode `peedelta <dir> ...` (physics_delta.f:peedelta with the
// driver's re-evaluation columns hx_*).
// SIDIS (task P12.2): `semi <data/sidis>` as arguments 7-8 runs reactions::SidisProduction and the "sidis" model
// against the P12.F1 tables (sidis/generate_semi_fixedrng, complete_ev_semi, complete_recon_ev_semi, complete_main_semi,
// event_capture_semi).
// Exit: 0 pass, 1 fail, 77 reference missing.
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
#include "simc/core/Errors.h"
#include "simc/core/Registry.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/driver/SimcRun.h"
#include "simc/driver/Central.h"
#include "simc/input/InpReader.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/radiation/SimcRadiator.h"
#include "simc/radiation/kernel/RadcState.h"
#include "simc/physics/delta/DeltaModel.h"
#include "simc/physics/kaon/KaonModel.h"
#include "simc/physics/pion/PionCrossSections.h"
#include "simc/reactions/DeltaProduction.h"
#include "simc/reactions/EventChain.h"
#include "simc/reactions/KaonProduction.h"
#include "simc/reactions/PionProduction.h"
#include "simc/reactions/SidisProduction.h"
#include "simc/reactions/RhoProduction.h"
#include "simc/core/SimcConstants.h"
#include "simc/physics/sidis/SidisModel.h"
#include "simc/reactions/StandardBeam.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/setup/SetupBuilder.h"

using namespace simc;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

int iv(const Row& r, const char* c) { return static_cast<int>(std::lround(r[c])); }
bool same(double a, double b) { return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof a) == 0; }

std::string gSpectral, gPion;   // data directories
std::string gKaon;              // data/kaon: non-empty = the kaon tables (P10.2)
bool isKaon() { return !gKaon.empty(); }
bool gDelta = false;            // the Delta tables (P11.1)
bool isDelta() { return gDelta; }
std::string gSidis;             // data/sidis: non-empty = the SIDIS tables (P12.2)
bool isSemi() { return !gSidis.empty(); }
bool gRho = false;              // the rho tables (P13.2)
bool isRho() { return gRho; }

/// cfg -> SIMC input file (reference/functions/pion_state/configs.txt)
std::string cfgInp(const std::string& inpdir, int cfg) {
    if (isRho()) {   // reference/functions/rho/rho_state/configs.txt
        static const char* d[] = {"",          "p13f1_rh_k1_rad_cap", "p13f1_rd_k1_rad_cap", "l1_h_k1_norad",
                                  "l1_h_k1_pim", "l1_h_k2_rad",        "l1_he3_k1_rad",       "l1_h_k1_dec1"};
        if (cfg < 1 || cfg > 7) throw std::runtime_error("bad cfg " + std::to_string(cfg));
        return inpdir + "/" + d[cfg] + ".inp";
    }
    if (isSemi()) {   // reference/functions/sidis/sidis_state/configs.txt
        static const char* d[] = {"",
                                  "p12f1_sh_k1_pip_rad_cap",
                                  "p12f1_sd_k1_pip_rad_cap",
                                  "l1_h_k1_pip_norad",
                                  "l1_h_k1_pim_rad",
                                  "l1_d_k1_pip_fermi",
                                  "l1_d_k1_pim_rad",
                                  "l1_h_k1_kp_rad",
                                  "l1_h_k1_km_rad",
                                  "l1_d_k1_kp_rad",
                                  "l1_h_k2_pip_rad",
                                  "l1_c12_k1_pip_rad"};
        if (cfg < 1 || cfg > 11) throw std::runtime_error("bad cfg " + std::to_string(cfg));
        return inpdir + "/" + d[cfg] + ".inp";
    }
    if (isDelta()) {   // reference/functions/delta_state/configs.txt
        static const char* d[] = {"", "p11f1_dh_d1_rad_cap", "p11f1_dh_d2_rad_cap", "l1_h_d1_norad",
                                  "l1_h_d2_norad", "l1_h_d1_wp1", "l1_h_d1_wp2"};
        if (cfg < 1 || cfg > 6) throw std::runtime_error("bad cfg " + std::to_string(cfg));
        return inpdir + "/" + d[cfg] + ".inp";
    }
    if (isKaon()) {   // reference/functions/kaon_state/configs.txt
        static const char* k[] = {"",
                                  "p10f1_kl_h_dec1_cap",
                                  "p10f1_ksm_d_dec1_cap",
                                  "l1_h_k1_norad",
                                  "l1_h_k2_dec0",
                                  "l1_h_k1_sig0_dec0",
                                  "l1_d_k1_lam_dec1",
                                  "l1_d_k2_sig0_dec0",
                                  "l1_he3_k1_dec1",
                                  "l1_he3coh_k2c_dec1",
                                  "l1_c12_k1_dec0",
                                  "l1_c12coh_k2c_dec1",
                                  "l1_he3coh11_k2c_dec1",
                                  "l1_he3coh12_k2c_dec1",
                                  "l1_he4coh_k2c_dec1",
                                  "l1_d_k2_sigm_dec0"};
        if (cfg < 1 || cfg > 15) throw std::runtime_error("bad cfg " + std::to_string(cfg));
        return inpdir + "/" + k[cfg] + ".inp";
    }
    static const char* f[] = {"",
                              "p9f1_pip_h_dec1_cap",
                              "p9f1_pim_d_dec1_cap",
                              "l1_h_k1_norad",
                              "l1_h_k2_rad",
                              "l1_he3coh_k2",
                              "l1_d_k2_pip",
                              "l1_he3qf_k2_pip",
                              "l1_h_k2_pi0",
                              "l1_d_k2_pi0n",
                              "l1_c12coh_k2_pim",
                              "l1_h_k2_delta2",
                              "l1_d_k2_delta3",
                              "l1_d_k2_pi0delta3",
                              "l1_h_k2_pi0delta2",
                              "l1_h_k2_pi0n"};
    if (cfg < 1 || cfg > 15) throw std::runtime_error("bad cfg " + std::to_string(cfg));
    return inpdir + "/" + f[cfg] + ".inp";
}

// ------------------------------------------------------------------------------------------------- event fields
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
void getArm(const std::string& p, const Row& r, ArmFull& a) {
    a.delta = r[p + "delta"];
    a.xptar = r[p + "xptar"];
    a.yptar = r[p + "yptar"];
    a.z = r[p + "z"];
    a.theta = r[p + "theta"];
    a.phi = r[p + "phi"];
    a.E = r[p + "E"];
    a.P = r[p + "P"];
}
#define EVFIELDS(X)                                                                                                    \
    X(Ein) X(Em) X(Pm) X(Emiss) X(Pmiss) X(Pmx) X(Pmy) X(Pmz) X(PmPar) X(PmPer) X(PmOop) X(nu) X(q) X(Q2) X(Trec)     \
        X(W) X(Mrec) X(epsilon) X(theta_pq) X(theta_tarq) X(phi_pq) X(phi_targ) X(beta) X(phi_s) X(phi_c) X(zhad)    \
            X(pt2) X(xbj)
void putEvent(const std::string& p, const Event& v, Outputs& o) {
#define PUTF(f) o[p + #f] = v.f;
    EVFIELDS(PUTF)
#undef PUTF
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
void getEvent(const std::string& p, const Row& r, Event& v) {
#define GETF(f) v.f = r[p + #f];
    EVFIELDS(GETF)
#undef GETF
    getArm(p + "e_", r, v.e);
    getArm(p + "p_", r, v.p);
    v.ue = {r[p + "ue_x"], r[p + "ue_y"], r[p + "ue_z"]};
    v.up = {r[p + "up_x"], r[p + "up_y"], r[p + "up_z"]};
    v.uq = {r[p + "uq_x"], r[p + "uq_y"], r[p + "uq_z"]};
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
        o["t_teff_" + std::to_string(i + 1)] = t.teff[static_cast<std::size_t>(i)];
        o["t_Eloss_" + std::to_string(i + 1)] = t.Eloss[static_cast<std::size_t>(i)];
    }
    o["t_Coulomb"] = t.Coulomb;
}
void putRad(const radiation::RadcState& st, Outputs& o) {
    o["r_ntail"] = st.ntail;
    for (int i = 1; i <= 3; ++i) o["r_Egamma_used_" + std::to_string(i)] = st.Egamma_used[static_cast<std::size_t>(i)];
    o["r_hardcorfac"] = st.hardcorfac;
    for (int i = 1; i <= 3; ++i) o["r_frac_" + std::to_string(i)] = st.frac[static_cast<std::size_t>(i)];
    o["r_g_4"] = st.g[4];
    o["r_c_4"] = st.c[4];
    o["r_rad_proton_ev"] = st.rad_proton_this_ev;
}
void putPfermi(const FermiMomentum& f, Outputs& o) {
    o["pfer"] = f.pfer;
    o["pferx"] = f.pferx;
    o["pfery"] = f.pfery;
    o["pferz"] = f.pferz;
    o["efer"] = f.efer;
}
void getPfermi(const Row& r, FermiMomentum& f) {
    f.pfer = r["pfer"];
    f.pferx = r["pferx"];
    f.pfery = r["pfery"];
    f.pferz = r["pferz"];
    f.efer = r["efer"];
}
void putGammas(const NtupVars& n, Outputs& o) {
    for (int i = 0; i < 4; ++i) {
        o["ntup_g1_" + std::to_string(i + 1)] = n.gamma1[static_cast<std::size_t>(i)];
        o["ntup_g2_" + std::to_string(i + 1)] = n.gamma2[static_cast<std::size_t>(i)];
    }
}

/// IRadiator decorator: records whether generate_rad was called.
struct RecordingRadiator : IRadiator {
    radiation::SimcRadiator inner;
    bool called{false};
    std::string name() const override { return "recording"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(const StageContext& c) override { inner.init(c); }
    void initEvent(MainEvent& m, const Event& v) override { inner.initEvent(m, v); }
    bool radiate(MainEvent& m, Event& v, Event& o, IRandom& r, const CompleteEvFn& c) override {
        called = true;
        return inner.radiate(m, v, o, r, c);
    }
};

/// One Fortran process (es_init of estate9.f): set-up, target, radiator, beam, the "pion" model, the reaction.
struct Process {
    SimcParameters params;
    RunSetup setup;
    material::HallCCryoTarget target;
    RecordingRadiator rad;
    reactions::StandardBeam beam;
    std::unique_ptr<IPhysicsModel> model;
    std::unique_ptr<reactions::PionProduction> reaction;
    reactions::MainContext ctx;
    NtupVars centralNtup;   ///< COMMON /ntuple_variables/ as calculate_central leaves it (semi: ntup%xfermi)
};

std::unique_ptr<Process> makeProcess(const std::string& inp, bool compat = true) {
    auto P = std::make_unique<Process>();
    InpReader rd(P->params);
    if (!rd.loadRun(inp)) throw std::runtime_error("cannot read " + inp);
    SetupOptions opt;
    opt.fortran_compat = compat;
    opt.data_dir = gSpectral;
    P->setup = BuildRunSetup(P->params, opt);
    if (!P->setup.ok()) throw std::runtime_error(inp + ": setup STOP " + P->setup.stops.front());
    P->setup.data_dir = gSpectral;
    P->target.setLog(nullptr);
    P->target.configure(Json{{"Mh", P->setup.Mh}}, P->params);
    P->target.setLog(nullptr);
    P->target.init(P->setup);
    P->rad.inner.setLog(nullptr);
    P->rad.init(P->setup);
    StageContext sc(P->setup);
    sc.target = &P->target;
    sc.radiator = &P->rad;
    P->beam.configure(Json::object(), P->params);
    P->beam.init(sc);
    P->model = Registry<IPhysicsModel>::create(
        isRho() ? "rho" : (isSemi() ? "sidis" : (isDelta() ? "delta" : (isKaon() ? "kaon" : "pion"))));
    P->model->configure(Json{{"data_dir", isSemi() ? gSidis : (isKaon() ? gKaon : gPion)}, {"verbose", false}}, P->params);
    P->model->init(sc);
    sc.physics = P->model.get();
    if (isRho())
        P->reaction = std::make_unique<reactions::RhoProduction>();
    else if (isSemi())
        P->reaction = std::make_unique<reactions::SidisProduction>();
    else if (isDelta())
        P->reaction = std::make_unique<reactions::DeltaProduction>();
    else if (isKaon())
        P->reaction = std::make_unique<reactions::KaonProduction>();
    else
        P->reaction = std::make_unique<reactions::PionProduction>();
    P->reaction->configure(Json{{"verbose", false}, {"data_dir", gPion}}, P->params);
    P->reaction->init(sc);
    P->ctx = P->reaction->mainContext();
    // delta (P11.1): es_init ends with calculate_central, whose radc_init_ev leaves the COMMON /radccom/ values that a
    // row whose complete_ev fails before radc_init_ev reports (cfg 6, which_pion 2: every row fails)
    // semi (P12.2): es_init's calculate_central also runs the central peepiX (SAVEd zaero, survivalprob, PDF warnings)
    // rho (P13.2): es_init's calculate_central also leaves the radiative COMMON of its radc_init_ev
    if (isDelta() || isSemi() || isRho())
        P->centralNtup = driver::calculateCentral(P->setup, *P->reaction, P->rad, *P->model).main0.ntup;
    return P;
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
    // rho (P13.F1): the drivers group the output rows by cfg; output.txt echoes every input column -> it is the input
    if (isRho() && in.rows.size() == ref.rows.size()) in = ref;
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

/// One generate call as the SIMC main loop does it (simc.f:174-183 resets; event.f:generate).
bool runGenerate(Process& P, long ntried, MainEvent& m, Event& v, Event& o, ScriptedRandom& rng, Outputs& out) {
    m.ntup.radphot = 0.0;
    m.ntup.radarm = 0.0;
    m.ntup.resfac = 0.0;
    m.ntup.sigcm = 0.0;
    m.ntup.krel = 0.0;
    m.ntup.mm = 0.0;
    m.ntup.mmA = 0.0;
    m.ntup.t = 0.0;
    m.ntup.decdist = 0.0;
    P.rad.called = false;
    P.reaction->beginTry(ntried);
    const bool ok = reactions::generateEvent(P.beam, *P.reaction, m, v, o, rng);
    out["success"] = ok;
    out["fail_stage"] = ok ? 0 : (P.rad.called ? 2 : 1);
    if (isRho()) {
        auto& R = dynamic_cast<reactions::RhoProduction&>(*P.reaction);
        // fail_stage as drv_generate13: 1 = generate_rad not called (complete_ev failed, or rho_decay failed without
        // radiation), 2 = generate_rad called (and it or rho_decay failed)
        out["rad_success"] = R.lastRadSuccess();
        out["Mh"] = m.gnrl.Mh;
        out["Mh2"] = m.gnrl.Mh2;
        out["Mh2_final"] = m.gnrl.Mh2_final;
        out["ntup_rhomass"] = m.ntup.rhomass;
        out["ntup_rhotheta"] = m.ntup.rhotheta;
        out["rho_decay_called"] = R.lastDecayCalled();
        out["rho_decay_success"] = R.lastDecaySuccess();
        out["rho_decay_n_unif"] = R.lastDecayNUnif();
    }
    out["rad_called"] = P.rad.called;
    out["n_unif"] = static_cast<double>(rng.nUniform());
    out["n_gauss"] = static_cast<double>(rng.nGauss());
    out["rng_err"] = rng.error();
    out["radphot"] = m.ntup.radphot;
    out["radarm"] = m.ntup.radarm;
    out["ntup_krel"] = m.ntup.krel;
    out["ntup_xfermi"] = m.ntup.xfermi;
    out["targ_Mrec_struck"] = P.reaction->mrecStruck();
    putPfermi(m.pfermi, out);
    putGammas(m.ntup, out);
    putMain(m, out);
    putTarget(m.target, out);
    putEvent("v_", v, out);
    putEvent("o_", o, out);
    putRad(P.rad.inner.state(), out);
    return ok;
}

/// complete_ev inputs as generate leaves them (drv_complete_ev9.f / the capture).
void loadCompleteEvInputs(const Process& P, const Row& r, MainEvent& m, Event& v, bool calcAngles) {
    const RunSetup& s = P.setup;
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
    v.Em = r["v_Em_in"];
    if (isSemi()) {   // the hadron generated by generate (event.f:289-298) and the stale vertex%Pm (event.f:922)
        v.p.E = r["p_E"];
        v.p.P = r["p_P"];
        v.p.delta = r["p_delta"];
        v.Pm = r["in_v_Pm"];
    }
    if (calcAngles) {
        PhysicsAngles(s.spec.e.theta, s.spec.e.phi, v.e.xptar, v.e.yptar, v.e.theta, v.e.phi);
        PhysicsAngles(s.spec.p.theta, s.spec.p.phi, v.p.xptar, v.p.yptar, v.p.theta, v.p.phi);
    }
    getPfermi(r, m.pfermi);
    m.target.z = r["t_z"];
    m.target.teff[0] = r["t_teff_1"];
}

void putRhoGnrl(const MainEvent& m, Outputs& o) {
    o["Mh"] = m.gnrl.Mh;
    o["Mh2"] = m.gnrl.Mh2;
    o["Mh2_final"] = m.gnrl.Mh2_final;
    o["ntup_rhomass"] = m.ntup.rhomass;
}

void putCompleteEv(Process& P, bool ok, const MainEvent& m, const Event& v, const ScriptedRandom& rng, Outputs& o) {
    o["success"] = ok;
    o["n_gauss"] = static_cast<double>(rng.nGauss());
    o["n_unif"] = static_cast<double>(rng.nUniform());
    o["rng_err"] = rng.error();
    putEvent("v_", v, o);
    putMain(m, o);
    o["jacobian"] = m.jacobian;
    o["t_teff_2"] = m.target.teff[1];
    o["t_teff_3"] = m.target.teff[2];
    o["t_Eloss_2"] = m.target.Eloss[1];
    o["t_Eloss_3"] = m.target.Eloss[2];
    o["ntup_krel"] = m.ntup.krel;
    o["targ_Mrec_struck"] = P.reaction->mrecStruck();
    putGammas(m.ntup, o);
    putRad(P.rad.inner.state(), o);
}

/// complete_main inputs (drv_complete_main9.f / the capture) and the call; outputs as the tables.
bool doCompleteMain(Process& P, const Row& r, Outputs& o) {
    MainEvent m;
    Event v, rec, v0;
    getEvent("v_", r, v);
    getEvent("r_", r, rec);
    m.jacobian = r["in_jacobian"];
    m.gen_weight = r["in_gen_weight"];
    m.SF_weight = r["in_SF_weight"];
    m.sigcc = r["in_sigcc"];
    m.sigcc_recon = r["in_sigcc_recon"];
    m.weight = r["in_weight"];
    if (r.has("in_m_epsilon")) m.epsilon = r["in_m_epsilon"];
    if (r.has("in_m_t")) m.t = r["in_m_t"];
    if (r.has("in_m_phi_pq")) m.phi_pq = r["in_m_phi_pq"];
    if (isKaon() || isSemi()) {   // the fields peeK's / peepiX's survival-probability block reads
        m.theta_pq = r["in_m_theta_pq"];
        m.FP.p.path = r["in_m_FP_p_path"];
        m.FP.p.dx = r["in_m_FP_p_dx"];
        m.FP.p.dy = r["in_m_FP_p_dy"];
        m.SP.p.delta = r["in_m_SP_p_delta"];
        m.ntup.decdist = r["decdist_in"];
    }
    getPfermi(r, m.pfermi);
    if (r.has("in_ntup_sigcm1")) m.ntup.sigcm1 = r["in_ntup_sigcm1"];
    if (r.has("in_ntup_sigcm2")) m.ntup.sigcm2 = r["in_ntup_sigcm2"];
    if (r.has("in_ntup_sigcm")) m.ntup.sigcm = r["in_ntup_sigcm"];
    reactions::MainContext ctx = P.ctx;
    if (r.has("ov_coulomb") && iv(r, "ov_coulomb") == 1) {
        ctx.using_Coulomb = iv(r, "ov_using_coulomb") != 0;
        ctx.Coulomb_ave = r["ov_coulomb_ave"];
    }
    P.reaction->mainContext() = ctx;
    // kaon: the row overrides of doing_decay / hadron_arm / ctau (ov_state, drv_complete_main10)
    kaon::KaonModel* km = isKaon() ? dynamic_cast<kaon::KaonModel*>(P.model.get()) : nullptr;
    kaon::PeekConfig savedPeek;
    SimcParameters savedParams;
    if (km) {
        savedPeek = km->config();
        savedParams = P.setup.params;
        if (r.has("ov_state") && iv(r, "ov_state") == 1) {
            km->config().doing_decay = iv(r, "ov_doing_decay") != 0;
            km->config().hadron_arm = iv(r, "ov_hadron_arm");
            km->config().ctau = r["ov_ctau"];
            P.setup.params.doing_decay = km->config().doing_decay;
        }
    }
    // delta: the row overrides of which_pion / targ%Mtar_struck (ov_state, drv_complete_main11)
    delta::DeltaModel* dm = isDelta() ? dynamic_cast<delta::DeltaModel*>(P.model.get()) : nullptr;
    delta::DeltaConfig savedDelta;
    if (dm) {
        savedDelta = dm->config();
        if (r.has("ov_state") && iv(r, "ov_state") == 1) {
            dm->config().which_pion = iv(r, "ov_which_pion");
            dm->config().Mtar_struck = r["ov_Mtar_struck"];
        }
        m.ntup.decdist = 0.0;
    }
    // semi (P12.2): main fields at entry and the row overrides of doing_hplus/doing_pizero/do_fermi/sigc_* (ov_state)
    sidis::SidisModel* sm = isSemi() ? dynamic_cast<sidis::SidisModel*>(P.model.get()) : nullptr;
    sidis::PeepiXConfig savedSidis;
    if (sm) {
        m.sigcent = r["in_sigcent"];
        m.ntup.dilu = r["in_ntup_dilu"];
        m.davejac = r["in_m_davejac"];
        m.johnjac = r["in_m_johnjac"];
        m.ntup.xfermi = r["in_ntup_xfermi"];
        savedSidis = sm->config();
        if (r.has("ov_state") && iv(r, "ov_state") == 1) {
            sm->config().doing_hplus = iv(r, "ov_hplus") != 0;
            sm->config().doing_pizero = iv(r, "ov_pizero") != 0;
            sm->config().do_fermi = iv(r, "ov_do_fermi") != 0;
            sm->config().sigc_flag = iv(r, "ov_sigc_flag");
            sm->config().sigc_nbin = iv(r, "ov_sigc_nbin");
            sm->config().sigc_kin_min = r["ov_sigc_kin_min"];
            sm->config().sigc_kin_max = r["ov_sigc_kin_max"];
            sm->config().sigc_kin_ind = r["ov_sigc_kin_ind"];
        }
    }
    const bool ok = P.reaction->completeMain(iv(r, "force_sigcc") == 1, m, v, v0, rec);
    if (sm) {
        o["decdist_out"] = m.ntup.decdist;
        o["sigcent"] = m.sigcent;
        o["ntup_xfermi"] = m.ntup.xfermi;
        o["ntup_dilu"] = m.ntup.dilu;
        o["eff_doing_hplus"] = sm->config().doing_hplus ? 1. : 0.;
        o["eff_doing_pizero"] = sm->config().doing_pizero ? 1. : 0.;
        o["eff_do_fermi"] = sm->config().do_fermi ? 1. : 0.;
        sm->config() = savedSidis;
    }
    P.reaction->mainContext() = P.ctx;
    if (dm) {
        o["decdist_out"] = m.ntup.decdist;
        o["eff_which_pion"] = dm->config().which_pion;
        o["eff_Mtar_struck"] = dm->config().Mtar_struck;
        dm->config() = savedDelta;
    }
    if (km) {
        o["decdist_out"] = m.ntup.decdist;
        o["eff_doing_decay"] = km->config().doing_decay ? 1. : 0.;
        o["eff_hadron_arm"] = km->config().hadron_arm;
        o["eff_ctau"] = km->config().ctau;
        km->config() = savedPeek;
        P.setup.params = savedParams;
    }
    o["success"] = ok;
    o["SF_weight"] = m.SF_weight;
    o["sigcc"] = m.sigcc;
    o["sigcc_recon"] = m.sigcc_recon;
    o["jacobian"] = m.jacobian;
    o["gen_weight"] = m.gen_weight;
    o["weight"] = m.weight;
    o["m_thetacm"] = m.thetacm;
    o["m_phicm"] = m.phicm;
    o["m_pcm"] = m.pcm;
    o["m_davejac"] = m.davejac;
    o["m_johnjac"] = m.johnjac;
    o["m_wcm"] = m.wcm;
    o["ntup_sigcm1"] = m.ntup.sigcm1;
    o["ntup_sigcm2"] = m.ntup.sigcm2;
    o["ntup_sigcm"] = m.ntup.sigcm;
    o["efer_out"] = m.pfermi.efer;   // rho: peerho writes the COMMON efer (D-057 F6)
    o["m_w"] = m.w;
    o["m_t"] = m.t;
    o["eff_using_Coulomb"] = ctx.using_Coulomb ? 1. : 0.;
    o["eff_coulomb_ave"] = ctx.Coulomb_ave;
    const double c = 1.0 + ctx.Coulomb_ave / ctx.Ebeam;
    o["hx_coul_fac"] = ctx.using_Coulomb ? c * c : 1.0;
    const double hp = m.SF_weight * m.jacobian * m.gen_weight * m.sigcc;
    o["hx_prod"] = hp;
    o[(isKaon() || isDelta() || isSemi() || isRho()) ? "hx_ratio" : "hx_tgtweight"] = hp != 0.0 ? m.weight / hp : 0.0;
    return ok;
}

// rho (P13.2): the stale state a try of the SIMC main loop starts from (drv_generate13 / the capture): the previous rho's
// lab theta/phi in vertex%p%xptar/yptar, COMMON /gnrl/ Mh, Mh2, ntup%rhomass, rhotheta and (table only) orig.
void setRhoStale(const Row& r, MainEvent& m, Event& v, Event& orig) {
    v.p.xptar = r["in_v_p_xptar"];
    v.p.yptar = r["in_v_p_yptar"];
    m.gnrl = HadronMass{r["in_Mh"], r["in_Mh2"], r.has("in_Mh2_final") ? r["in_Mh2_final"] : r["in_Mh2"]};
    m.gnrl_per_event = true;
    m.ntup.rhomass = r["in_rhomass"];
    m.ntup.rhotheta = r["in_rhotheta"];
    if (r.has("in_o_Ein")) getEvent("in_o_", r, orig);
}

/// rho complete_ev inputs as generate leaves them (drv_complete_ev13.f / the capture).
void loadCompleteEvInputsRho(const Process& P, const Row& r, MainEvent& m, Event& v, bool calcAngles) {
    const RunSetup& s = P.setup;
    v.Ein = r["Ein"];
    v.e.xptar = r["e_xptar"];
    v.e.yptar = r["e_yptar"];
    v.e.theta = r["e_theta"];
    v.e.phi = r["e_phi"];
    v.e.E = r["e_E"];
    v.e.P = r["e_P"];
    v.e.delta = r["e_delta"];
    v.p.xptar = r["in_p_xptar"];
    v.p.yptar = r["in_p_yptar"];
    v.p.theta = r["p_theta"];
    v.p.phi = r["p_phi"];
    v.p.E = r["in_p_E"];
    v.p.P = r["in_p_P"];
    v.p.delta = r["in_p_delta"];
    v.up = {r["in_up_x"], r["in_up_y"], r["in_up_z"]};
    v.Em = r["v_Em_in"];
    v.Pm = r["in_v_Pm"];
    v.Mrec = r["in_v_Mrec"];
    if (calcAngles) {
        PhysicsAngles(s.spec.e.theta, s.spec.e.phi, v.e.xptar, v.e.yptar, v.e.theta, v.e.phi);
        PhysicsAngles(s.spec.p.theta, s.spec.p.phi, v.p.xptar, v.p.yptar, v.p.theta, v.p.phi);
    }
    getPfermi(r, m.pfermi);
    m.target.z = r["t_z"];
    m.target.teff[0] = r["t_teff_1"];
    m.gnrl = HadronMass{r["in_Mh"], r["in_Mh2"], r["in_Mh2"]};
    m.gnrl_per_event = true;
    m.ntup.rhomass = r["in_rhomass"];
}

int runEventTable(const std::string& name, const std::string& dir, const std::string& inpdir, double tol) {
    std::map<int, std::unique_ptr<Process>> procs;   // one Fortran process per cfg
    auto proc = [&](int cfg) -> Process& {
        auto& P = procs[cfg];
        if (!P) P = makeProcess(cfgInp(inpdir, cfg));
        return *P;
    };
    if (name == "generate") {
        return runAuto(std::string(isDelta() ? "[delta tables] " : "") + "generate_pion_fixedrng (event.f:126-428, pion branches, scripted grnd/gauss1)", dir, tol,
                       [&](const Row& r, Outputs& o) {
                           Process& P = proc(iv(r, "cfg"));
                           MainEvent m;   // drv_generate9: main, vertex, orig = 0 per row
                           if (isSemi()) m.ntup.xfermi = P.centralNtup.xfermi;   // COMMON ntup%xfermi of the central peepiX
                           Event v, orig;
                           v.Pm = r["in_v_Pm"];
                           if (isRho()) setRhoStale(r, m, v, orig);
                           for (int i = 1; i <= 3; ++i)
                               P.rad.inner.state().Egamma_min[static_cast<std::size_t>(i)] =
                                   P.rad.inner.state().Egamma_max[static_cast<std::size_t>(i)] = 0.;
                           std::vector<double> u, g;
                           for (int i = 1; r.has("u" + std::to_string(i)); ++i) u.push_back(r["u" + std::to_string(i)]);
                           for (int i = 1; r.has("g" + std::to_string(i)); ++i) g.push_back(r["g" + std::to_string(i)]);
                           ScriptedRandom rng;
                           rng.load(g, u);
                           runGenerate(P, iv(r, "gen_ntried"), m, v, orig, rng, o);
                       });
    }
    if (name == "complete_ev") {
        return runAuto(std::string(isDelta() ? "[delta tables] " : "") + "complete_ev_pion (event.f:432-1052, pion branches)", dir, tol, [&](const Row& r, Outputs& o) {
            const int cfg = iv(r, "cfg");
            Process& P = proc(cfg);
            P.reaction->setMrecStruck(P.setup.targ.Mrec_struck);   // drv_complete_ev9: restored before every row
            MainEvent m;
            Event v;
            if (isRho())
                loadCompleteEvInputsRho(P, r, m, v, iv(r, "calc_angles") == 1);
            else
                loadCompleteEvInputs(P, r, m, v, iv(r, "calc_angles") == 1);
            o["in_e_theta"] = v.e.theta;
            o["in_e_phi"] = v.e.phi;
            o["in_p_theta"] = v.p.theta;
            o["in_p_phi"] = v.p.phi;
            std::vector<double> g, u;
            for (int i = 1; i <= 12; ++i) g.push_back(r["g" + std::to_string(i)]);
            for (int i = 1; i <= 4; ++i) u.push_back(r["u" + std::to_string(i)]);
            ScriptedRandom rng;
            rng.load(g, u);
            P.reaction->beginTry(iv(r, "gen_ntried"));
            const bool ok = P.reaction->completeEv(m, v, rng);
            putCompleteEv(P, ok, m, v, rng, o);
            if (isRho()) putRhoGnrl(m, o);
        });
    }
    if (name == "complete_recon_ev") {
        return runAuto(std::string(isDelta() ? "[delta tables] " : "") + "complete_recon_ev_pion (event.f:1056-1359)", dir, tol, [&](const Row& r, Outputs& o) {
            Process& P = proc(iv(r, "cfg"));
            Event rec;
            getArm("in_e_", r, rec.e);
            getArm("in_p_", r, rec.p);
            if (isRho()) dynamic_cast<reactions::RhoProduction&>(*P.reaction).setReconMh2(r["in_Mh2"]);
            o["success"] = P.reaction->completeRecon(rec);
            putEvent("r_", rec, o);
            o["ntup_mm"] = P.reaction->reconNtup().mm;
            o["ntup_mmA"] = P.reaction->reconNtup().mmA;
            o["ntup_t"] = P.reaction->reconNtup().t;
        });
    }
    if (name == "complete_main") {
        std::size_t identity_bad = 0;
        const int rc = runAuto(std::string(isDelta() ? "[delta tables] " : "") + "complete_main_pion (event.f:1363-1569, pion branch)", dir, tol,
                               [&](const Row& r, Outputs& o) {
                                   Process& P = proc(iv(r, "cfg"));
                                   doCompleteMain(P, r, o);
                               });
        (void)identity_bad;
        return rc;
    }
    if (name == "peedelta" && isDelta()) {
        return runAuto("peedelta (physics_delta.f:1-135 + driver re-evaluation)", dir, tol, [&](const Row& r, Outputs& o) {
            Process& P = proc(iv(r, "cfg"));
            MainEvent m;
            Event v;
            getEvent("v_", r, v);
            m.epsilon = r["in_m_epsilon"];
            m.t = r["in_m_t"];
            m.phi_pq = r["in_m_phi_pq"];
            getPfermi(r, m.pfermi);
            m.ntup.sigcm1 = r["in_ntup_sigcm1"];
            m.ntup.sigcm2 = r["in_ntup_sigcm2"];
            m.ntup.sigcm = r["in_ntup_sigcm"];
            delta::DeltaConfig cfg = dynamic_cast<delta::DeltaModel&>(*P.model).config();
            if (r.has("ov_state") && iv(r, "ov_state") == 1) {
                cfg.which_pion = iv(r, "ov_which_pion");
                cfg.Mtar_struck = r["ov_Mtar_struck"];
            }
            delta::PeedeltaTrace tr;
            const double val = delta::peedelta(v, m, cfg, nullptr, &tr);
            o["peedelta"] = val;
            o["m_thetacm"] = m.thetacm;
            o["m_phicm"] = m.phicm;
            o["m_pcm"] = m.pcm;
            o["m_davejac"] = m.davejac;
            o["m_johnjac"] = m.johnjac;
            o["m_wcm"] = m.wcm;
            o["ntup_sigcm1"] = m.ntup.sigcm1;
            o["ntup_sigcm2"] = m.ntup.sigcm2;
            o["ntup_sigcm"] = m.ntup.sigcm;
            o["eff_which_pion"] = cfg.which_pion;
            o["eff_Mtar_struck"] = cfg.Mtar_struck;
            o["hx_tfcos"] = tr.tfcos;
            o["hx_s"] = tr.s;
            o["hx_k_eq"] = tr.k_eq;
            o["hx_fac"] = tr.fac;
            o["hx_gtpr"] = tr.gtpr;
            o["hx_t_seen"] = tr.t_arg;
            o["hx_q2_seen"] = tr.q2_arg;
            o["hx_s_seen"] = tr.s_arg;
            o["hx_t_gev"] = m.t / 1.e6;
            o["hx_q2_gev"] = v.Q2 / 1.e6;
            o["hx_s_gev"] = tr.s / 1.e6;
            o["hx_sigblok_seen"] = pion::sigBlok(m.thetacm, m.phicm, tr.t_arg, tr.q2_arg, tr.s_arg, m.epsilon,
                                                 cfg.Mtar_struck / 1000., cfg.which_pion);
            // the corrected path (D-048 F3): peedelta with the doubles must give the intended value
            delta::DeltaConfig cc = cfg;
            cc.sigblok_binary128_quirk = false;
            MainEvent m2 = m;
            delta::peedelta(v, m2, cc);
            o["hx_sigblok_intended"] = m2.ntup.sigcm1;
            o["hx_ratio"] = val / (m.davejac * (tr.gtpr * tr.fac));
        });
    }
    std::printf("unknown table %s\n", name.c_str());
    return 1;
}

// ------------------------------------------------------------------------------------------------- pion_state
struct Cmp {
    Cmp(std::string w, double t) : what(std::move(w)), tol(t) {}
    std::string what;
    double tol;
    std::map<std::string, std::size_t> bad;
    std::size_t n{0}, eq{0};
    double maxrel{0};
    std::string worst;
    void check(const std::string& col, double cpp, double ref, std::size_t row) {
        ++n;
        if (same(cpp, ref)) {
            ++eq;
            return;
        }
        const double rel = simcl1::relDiff(cpp, ref);
        if (rel > maxrel || std::isnan(rel)) {
            maxrel = rel;
            worst = col;
        }
        if (!(rel <= tol)) {
            if (bad[col]++ == 0)
                std::printf("    %s %s row %zu: ref %.17g cpp %.17g\n", what.c_str(), col.c_str(), row, ref, cpp);
        }
    }
    bool report() const {
        std::size_t nb = 0;
        for (const auto& kv : bad) nb += kv.second;
        std::printf("  %-40s %9zu values, %9zu bitwise, %zu beyond tol, max|rel| %.3e %s\n", what.c_str(), n, eq, nb,
                    maxrel, worst.c_str());
        return nb == 0 && n > 0;
    }
};

int runPionState(const std::string& dir, const std::string& inpdir, double tol) {
    simcl1::Table out;
    std::string err;
    if (!out.load(dir + "/output.txt", err)) {
        std::printf("pion_state: SKIP (%s)\n", err.c_str());
        return 77;
    }
    Cmp c{isDelta() ? "delta_state (start-up + calculate_central)"
                    : (isKaon() ? "kaon_state (start-up + calculate_central)" : "pion_state (start-up + calculate_central)"),
          tol};
    for (std::size_t i = 0; i < out.rows.size(); ++i) {
        Row r(out, i);
        const int cfg = iv(r, "cfg");
        SimcParameters p;
        InpReader rd(p);
        if (!rd.loadRun(cfgInp(inpdir, cfg))) throw std::runtime_error("cannot read cfg " + std::to_string(cfg));
        driver::RunOptions opt;
        opt.fortran_compat = true;
        opt.data_dir = gSpectral;
        opt.quiet_stages = true;
        opt.run_name = "pion_state_cfg" + std::to_string(cfg);
        std::ostringstream log;
        opt.log = &log;
        Json pipe{{"sinks", Json::array()},
                  {"physics", {{"type", "pion"}, {"data_dir", gPion}}},
                  {"reaction", {{"type", "pion_production"}, {"data_dir", gPion}}}};
        if (isKaon())
            pipe = Json{{"sinks", Json::array()},
                        {"physics", {{"type", "kaon"}, {"data_dir", gKaon}}},
                        {"reaction", {{"type", "kaon_production"}, {"data_dir", gPion}}}};
        if (isDelta())
            pipe = Json{{"sinks", Json::array()},
                        {"physics", {{"type", "delta"}}},
                        {"reaction", {{"type", "delta_production"}, {"data_dir", gPion}}}};
        if (isRho())
            pipe = Json{{"sinks", Json::array()}, {"physics", {{"type", "rho"}}}, {"reaction", {{"type", "rho_production"}}}};
        // pi0 configurations use the HMS/SHMS as hadron arm in the L1 inputs (the calorimeter is P14)
        driver::SimcRun run(p, pipe, opt);
        run.initialise();
        const RunSetup& s = run.setup();
        const TargetInfo& t = s.targ;
        const ReactionFlags& f = s.flags;
        Outputs o;
        o["doing_pion"] = f.doing_pion;
        o["doing_rho"] = f.doing_rho;
        o["doing_hydrho"] = f.doing_hydrho;
        o["doing_deutrho"] = f.doing_deutrho;
        o["doing_herho"] = f.doing_herho;
        o["doing_hplus"] = f.doing_hplus;
        o["doing_semi"] = f.doing_semi;
        o["do_fermi"] = f.do_fermi;
        o["doing_delta"] = f.doing_delta;
        o["doing_hyddelta"] = f.doing_hyddelta;
        o["doing_deutdelta"] = f.doing_deutdelta;
        o["doing_hedelta"] = f.doing_hedelta;
        o["doing_hyd_elast"] = f.doing_hyd_elast;
        o["doing_deuterium"] = f.doing_deuterium;
        o["doing_heavy"] = f.doing_heavy;
        o["doing_phsp"] = f.doing_phsp;
        o["doing_pizero"] = f.doing_pizero;
        o["doing_kaon"] = f.doing_kaon;
        o["which_kaon"] = s.params.which_kaon;
        o["doing_hydkaon"] = f.doing_hydkaon;
        o["doing_deutkaon"] = f.doing_deutkaon;
        o["doing_hekaon"] = f.doing_hekaon;
        o["doing_semika"] = f.doing_semika;
        o["which_pion"] = s.params.which_pion;
        o["pizero_ngamma"] = s.params.pizero_ngamma;
        o["doing_hydpi"] = f.doing_hydpi;
        o["doing_deutpi"] = f.doing_deutpi;
        o["doing_hepi"] = f.doing_hepi;
        o["doing_eep"] = f.doing_eep;
        o["doing_decay"] = s.params.doing_decay;
        o["ctau"] = s.params.ctau;
        o["electron_arm"] = s.params.electron_arm;
        o["hadron_arm"] = s.params.hadron_arm;
        o["using_rad"] = s.sw.using_rad;
        o["using_Eloss"] = s.sw.using_Eloss;
        o["using_Coulomb"] = s.sw.using_Coulomb;
        o["mc_smear"] = s.sw.mc_smear;
        o["sign_hadron"] = s.sign_hadron;
        o["targ_A"] = t.A;
        o["targ_Z"] = t.Z;
        o["targ_N"] = t.N;
        o["targ_mass_amu"] = t.mass_amu;
        o["targ_M"] = t.M;
        o["targ_mrec_amu"] = t.mrec_amu;
        o["targ_Mrec"] = t.Mrec;
        o["targ_Mtar_struck"] = t.Mtar_struck;
        o["targ_Mrec_struck"] = t.Mrec_struck;
        o["targ_rho"] = t.rho;
        o["targ_thick"] = t.thick;
        o["targ_length"] = t.length;
        o["targ_X0"] = t.X0;
        o["targ_X0_cm"] = t.X0_cm;
        o["targ_Coulomb_ave"] = t.Coulomb.ave;
        o["targ_Coulomb_constant"] = t.Coulomb_constant;
        for (int k = 0; k < 3; ++k) {
            const std::string ks = std::to_string(k + 1);
            const std::size_t kk = static_cast<std::size_t>(k);
            o["targ_Eloss_ave_" + ks] = t.Eloss[kk].ave;
            o["targ_teff_ave_" + ks] = t.teff[kk].ave;
            o["targ_musc_max_" + ks] = t.musc_max[kk];
        }
        o["Mh"] = s.Mh;
        o["Mh2"] = s.Mh2;
        o["Mh2_final"] = s.Mh2_final;
        o["Ebeam"] = s.Ebeam;
        o["dEbeam"] = s.dEbeam;
        o["Ebeam_vertex_ave"] = s.Ebeam_vertex_ave;
        o["nump"] = static_cast<double>(s.pdist.pval.size());
        o["pval_nump"] = s.pdist.pval.empty() ? 0.0 : s.pdist.pval.back();
        const EdgeTrue& V = s.VERTEXedge;
        auto cut = [&](const std::string& n, const Cut& x) {
            o[n + "_min"] = x.min;
            o[n + "_max"] = x.max;
        };
        cut("VE_Em", V.Em);
        cut("VE_Pm", V.Pm);
        cut("VE_Mrec", V.Mrec);
        cut("VE_Trec", V.Trec);
        cut("VE_Trec_struck", V.Trec_struck);
        cut("edge_e_E", s.edge.e.E);
        cut("edge_p_E", s.edge.p.E);
        cut("cuts_Em", s.cuts.Em);
        cut("gen_sumEgen", s.gen.sumEgen);
        cut("gen_e_E", s.gen.e.E);
        cut("gen_e_delta", s.gen.e.delta);
        cut("gen_e_yptar", s.gen.e.yptar);
        cut("gen_e_xptar", s.gen.e.xptar);
        cut("gen_p_E", s.gen.p.E);
        cut("gen_p_delta", s.gen.p.delta);
        cut("gen_p_yptar", s.gen.p.yptar);
        cut("gen_p_xptar", s.gen.p.xptar);
        o["Egamma_tot_max"] = s.rad.Egamma_tot_max;
        o["Egamma1_max"] = s.rad.Egamma1_max;
        o["Egamma2_max"] = s.rad.Egamma2_max;
        o["Egamma3_max"] = s.rad.Egamma3_max;
        const driver::CentralEvent& C = run.central().central;
        o["c_sigcc"] = C.sigcc;
        o["c_nu"] = C.nu;
        o["c_q"] = C.q;
        o["c_Q2"] = C.Q2;
        o["c_Em"] = C.Em;
        o["c_Pm"] = C.Pm;
        o["c_W"] = C.W;
        o["c_MM"] = C.MM;
        o["c_rad_hardcorfac"] = C.rad.hardcorfac;
        const MainEvent& m0 = run.central().main0;
        o["c_ntup_sigcm"] = m0.ntup.sigcm;
        o["c_ntup_sigcm1"] = m0.ntup.sigcm1;
        o["c_ntup_sigcm2"] = m0.ntup.sigcm2;
        o["c_efer"] = m0.pfermi.efer;
        putEvent("v0_", run.central().vertex0, o);
        o["c_decdist"] = m0.ntup.decdist;   // COMMON decdist after the central peeK (kaon)
        o["c_ntup_krel"] = m0.ntup.krel;
        o["c_ntup_rhomass"] = m0.ntup.rhomass;
        o["c_ntup_rhotheta"] = m0.ntup.rhotheta;
        o["c_ntup_t"] = m0.ntup.t;
        o["c_ntup_mm"] = m0.ntup.mm;
        o["c_ntup_mmA"] = m0.ntup.mmA;
        o["c_pfer"] = m0.pfermi.pfer;
        // rho (P13.2): no 3He spectral function / Em table is loaded (dbase.f:594, D-057 F4)
        o["nume"] = 0.;
        o["sf_numEm"] = 0.;
        o["sf_numPm"] = 0.;
        o["spec_e_P"] = s.spec.e.P;
        o["spec_e_theta"] = s.spec.e.theta;
        o["spec_e_phi"] = s.spec.e.phi;
        o["spec_p_P"] = s.spec.p.P;
        o["spec_p_theta"] = s.spec.p.theta;
        o["spec_p_phi"] = s.spec.p.phi;
        o["SPedge_p_delta_min"] = s.SPedge.p.delta.min;
        o["SPedge_p_delta_max"] = s.SPedge.p.delta.max;
        o["SPedge_e_delta_min"] = s.SPedge.e.delta.min;
        o["SPedge_e_delta_max"] = s.SPedge.e.delta.max;
        o["dE_edge_test"] = s.params.dE_edge_test;
        if (isKaon()) {   // spot values of /sigkaon/ (the kaon model's tables)
            const kaon::SaghaiTables& T = dynamic_cast<const kaon::KaonModel&>(run.model()).tables();
            o["weightc_1_1"] = T.wc(1, 1);
            o["weightc_20_50"] = T.wc(20, 50);
            o["weightd_1_1_1"] = T.wd(1, 1, 1);
            o["weightd_8_40_30"] = T.wd(8, 40, 30);
            o["zrff1_1_1_1"] = T.zrff[0][kaon::SaghaiTables::idxP(1, 1, 1)];
            o["ziff6_10_11_19"] = T.ziff[5][kaon::SaghaiTables::idxP(10, 11, 19)];
            o["zsrff1_1_1_1"] = T.zsrff[0][kaon::SaghaiTables::idxS(1, 1, 1)];
            o["zsiff6_20_10_19"] = T.zsiff[5][kaon::SaghaiTables::idxS(20, 10, 19)];
        }
        for (const auto& n : out.names) {
            if (n == "cfg") continue;
            auto it = o.find(n);
            if (it == o.end()) {
                std::printf("    column %s not produced\n", n.c_str());
                c.bad[n]++;
                continue;
            }
            c.check(n, it->second, r[n], r.number());
        }
    }
    bool ok = c.report();

    // momentum distributions (pdist.txt: cfg i pval mprob)
    simcl1::Table pd;
    if (pd.load(dir + "/pdist.txt", err)) {
        Cmp cp{"pion_state pdist.txt (pval, mprob)", tol};
        std::map<int, RunSetup> st;
        // configurations with a momentum distribution according to output.txt (nump > 0). pdist.txt also contains
        // 2 x 1050 rows labelled cfg 4 and 5 (H and He-3 coherent, nump = 0 in output.txt, no file read by dbase.f):
        // a labelling inconsistency of the reference, reported and not compared.
        std::map<int, double> nump;
        for (std::size_t i = 0; i < out.rows.size(); ++i) {
            Row r(out, i);
            nump[iv(r, "cfg")] = r["nump"];
        }
        std::map<int, long> skipped;
        for (std::size_t i = 0; i < pd.rows.size(); ++i) {
            Row r(pd, i);
            const int cfg = iv(r, "cfg");
            if (!(nump[cfg] > 0)) {
                ++skipped[cfg];
                continue;
            }
            auto it = st.find(cfg);
            if (it == st.end()) {
                SimcParameters p;
                InpReader rd(p);
                rd.loadRun(cfgInp(inpdir, cfg));
                SetupOptions so;
                so.fortran_compat = true;
                so.data_dir = gSpectral;
                it = st.emplace(cfg, BuildRunSetup(p, so)).first;
            }
            const std::size_t k = static_cast<std::size_t>(iv(r, "i") - 1);
            const MomentumDistribution& d = it->second.pdist;
            cp.check("pval", k < d.pval.size() ? d.pval[k] : -1., r["pval"], r.number());
            cp.check("mprob", k < d.mprob.size() ? d.mprob[k] : -1., r["mprob"], r.number());
        }
        ok = cp.report() && ok;
        for (const auto& kv : skipped)
            std::printf("    pdist.txt: %ld rows labelled cfg %d not compared (output.txt nump = 0 for that cfg)\n",
                        kv.second, kv.first);
    }

    if (isRho()) {
        // stops.txt (4He and C-12 (e,e'rho): STOP at dbase.f:876): refused at input in both modes (D-057 F16)
        int nref = 0;
        for (const double A : {4.0, 12.0}) {
            for (const bool compat : {true, false}) {
                SimcParameters p;
                InpReader rd(p);
                rd.loadRun(cfgInp(inpdir, 1));
                p.targ.A = A;
                SetupOptions so;
                so.fortran_compat = compat;
                so.data_dir = gSpectral;
                try {
                    BuildRunSetup(p, so);
                } catch (const std::invalid_argument& e) {
                    if (std::string(e.what()).find("only set up for H, D and 3He") != std::string::npos) ++nref;
                }
            }
        }
        std::printf("  STOP configurations A = 4, 12 (dbase.f:876): refused at input in both modes: %d of 4\n", nref);
        ok = ok && nref == 4;
        std::printf("pion_state: %s\n", ok ? "PASS" : "FAIL");
        return ok ? 0 : 1;
    }
    // STOP configurations (stops.txt): the set-up must stop with the Fortran message in fortran_compat mode
    std::ifstream sf(dir + "/stops.txt");
    int nstop = 0, badstop = 0;
    for (std::string line; std::getline(sf, line);) {
        if (line.empty() || line[0] == '#' || line.rfind("cfg", 0) == 0) continue;
        std::vector<std::string> fld;
        std::stringstream ss(line);
        for (std::string x; std::getline(ss, x, '|');) {
            const auto a = x.find_first_not_of(' '), b = x.find_last_not_of(' ');
            fld.push_back(a == std::string::npos ? "" : x.substr(a, b - a + 1));
        }
        if (fld.size() < 5) continue;
        ++nstop;
        const std::string inp = inpdir + "/" + fld[1].substr(fld[1].find('/') + 1) + ".inp";
        SimcParameters p;
        InpReader rd(p);
        std::string got;
        if (!rd.loadRun(inp)) {
            got = "<cannot read " + inp + ">";
        } else {
            SetupOptions so;
            so.fortran_compat = true;
            so.data_dir = gSpectral;
            const RunSetup s = BuildRunSetup(p, so);
            got = s.stops.empty() ? "<no STOP>" : s.stops.front();
        }
        const bool eq = got == fld[4];
        std::printf("    stop cfg %s: %-60s %s\n", fld[0].c_str(), got.c_str(), eq ? "ok" : ("REF: " + fld[4]).c_str());
        if (!eq) ++badstop;
    }
    std::printf("  STOP configurations: %d, %d differ\n", nstop, badstop);
    ok = ok && nstop == (isDelta() ? 3 : (isKaon() ? 4 : 8)) && badstop == 0;
    std::printf("pion_state: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

// ------------------------------------------------------------------------------------------------- capture
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
    std::string acc;
    char buf[1 << 16];
    while (std::fgets(buf, sizeof buf, f)) {
        acc += buf;
        if (acc.back() != '\n') continue;
        std::string line;
        line.swap(acc);
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
    void row(const Row& r, const simcl1::Table& t, const Outputs& o, const std::set<std::string>& exclude,
             const std::set<std::string>* skip = nullptr) {
        ++rows;
        bool any = false;
        for (const auto& n : t.names) {
            if (exclude.count(n)) continue;
            if (skip && skip->count(n)) continue;
            auto it = o.find(n);
            if (it == o.end()) continue;
            const double f = r[n], c = it->second;
            ++values;
            if (same(f, c)) {
                ++bitwise;
                continue;
            }
            const double rel = simcl1::relDiff(c, f);
            if (rel > maxrel) {
                maxrel = rel;
                worstCol = n;
            }
            if (!(rel <= tol)) {
                if (badCols[n]++ == 0)
                    std::printf("    %s: %s first mismatch at row %zu (ntried %.0f): ref %.17g cpp %.17g\n", what.c_str(),
                                n.c_str(), r.number(), r["ntried"], f, c);
                any = true;
            }
        }
        if (any) ++rowsBad;
    }
    bool report() const {
        std::printf("  %-18s rows %6zu, values %9zu, bitwise %9zu, rows with a mismatch %zu, max|rel| %.3e %s\n",
                    what.c_str(), rows, values, bitwise, rowsBad, maxrel, worstCol.c_str());
        for (const auto& kv : badCols) std::printf("      column %s: %zu\n", kv.first.c_str(), kv.second);
        return rowsBad == 0 && rows > 0;
    }
};

int runCapture(const std::string& capdir, const std::string& inp, double tol) {
    simcl1::Table rngT, genT, cevT, recT, mainT;
    std::string err;
    if (!loadAny(capdir + "/cap_rng.txt", rngT, err) || !loadAny(capdir + "/cap_generate.txt", genT, err) ||
        !loadAny(capdir + "/cap_complete_ev.txt", cevT, err) || !loadAny(capdir + "/cap_recon.txt", recT, err) ||
        !loadAny(capdir + "/cap_main.txt", mainT, err)) {
        std::printf("capture %s: SKIP (%s)\n", capdir.c_str(), err.c_str());
        return 77;
    }
    std::printf("L1 capture replay %s (setup %s)\n", capdir.c_str(), inp.c_str());
    std::map<long, std::pair<std::vector<double>, std::vector<double>>> stream;   // ntried -> (u, g), phases 1 and 6
    {
        const int cN = rngT.col("ntried"), cP = rngT.col("phase"), cK = rngT.col("kind"), cV = rngT.col("value");
        for (const auto& r : rngT.rows) {
            const int ph = static_cast<int>(r[static_cast<std::size_t>(cP)]);
            // semi (P12.F1 captures): phase 9 = the gauss1 of the tail-1 complete_ev re-call inside generate_rad
            // rho (P13.F1): 12 generate_rho, 13 rho_decay, 9/15 the tail-1 complete_ev re-call and its generate_rho
            if (ph != 1 && ph != 6 && !(isSemi() && ph == 9) &&
                !(isRho() && (ph == 9 || ph == 12 || ph == 13 || ph == 15)))
                continue;
            auto& q = stream[std::lround(r[static_cast<std::size_t>(cN)])];
            (static_cast<int>(r[static_cast<std::size_t>(cK)]) == 1 ? q.first : q.second)
                .push_back(r[static_cast<std::size_t>(cV)]);
        }
    }
    int fail = 0;
    // generate: one persistent process and persistent main/vertex/orig, as the SIMC main loop
    {
        auto P = makeProcess(inp);
        MainEvent m;
        Event v, orig;
        CaptureCompare cc{"generate", tol};
        std::size_t leftover = 0;
        const std::set<std::string> excl = {"ntried", "nevent", "n_unif_rad", "n_gauss_rad", "n_int", "in_v_Pm",
                                            "in_v_p_xptar", "in_v_p_yptar", "in_Mh", "in_Mh2", "in_Mh2_final",
                                            "in_rhomass", "in_rhotheta"};
        // rho (two-pass captures: gaps between the full-logged tries): the stale state comes from the record; orig is not
        // recorded at entry, so the o_* fields of failed tries are not compared, nor (fail_stage 1) the fields that
        // complete_ev / generate_rho leave unassigned (checks.txt of the capture)
        std::set<std::string> rhoSkipFail, rhoSkipCev;
        for (const auto& n : genT.names)
            if (n.rfind("o_", 0) == 0) rhoSkipFail.insert(n);
        rhoSkipCev = rhoSkipFail;
        for (const char* n : {"v_Pmx", "v_Pmy", "v_Pmz", "v_Pmiss", "v_Emiss", "v_PmPar", "v_PmOop", "v_PmPer", "v_Trec",
                              "v_p_E", "v_p_P", "v_p_delta", "t_Eloss_2", "t_Eloss_3", "t_teff_2", "t_teff_3", "m_Trec",
                              "m_w", "m_t", "m_tmin", "m_q2", "m_epsilon", "m_theta_pq", "m_phi_pq", "ntup_krel", "v_up_x",
                              "v_up_y", "v_up_z", "v_p_xptar", "v_p_yptar", "ntup_rhomass", "Mh", "Mh2", "Mh2_final",
                              "hardcorfac", "ntup_rhotheta"})
            rhoSkipCev.insert(n);
        for (std::size_t i = 0; i < genT.rows.size(); ++i) {
            Row r(genT, i);
            const long nt = std::lround(r["ntried"]);
            if (isRho()) setRhoStale(r, m, v, orig);
            const auto& q = stream[nt];
            ScriptedRandom rng;
            rng.load(q.second, q.first);
            Outputs o;
            runGenerate(*P, nt, m, v, orig, rng, o);
            if (rng.uniformsLeft() || rng.gaussiansLeft() || rng.error()) ++leftover;
            o["Ein_shift"] = m.Ein_shift;
            o["Ee_shift"] = m.Ee_shift;
            o["gen_weight"] = m.gen_weight;
            o["jacobian"] = m.jacobian;
            o["ntail"] = P->rad.inner.state().ntail;
            o["hardcorfac"] = P->rad.inner.state().hardcorfac;
            for (int k = 1; k <= 3; ++k)
                o["Egamma_used_" + std::to_string(k)] = P->rad.inner.state().Egamma_used[static_cast<std::size_t>(k)];
            if (isRho() && r["success"] == 0)
                cc.row(r, genT, o, excl, r["fail_stage"] == 1 ? &rhoSkipCev : &rhoSkipFail);
            else
                cc.row(r, genT, o, excl);
        }
        fail += cc.report() ? 0 : 1;
        std::printf("      tries whose captured random numbers were not consumed exactly: %zu\n", leftover);
        if (leftover) ++fail;
    }
    // complete_ev (first call of every try); on failures the fields complete_ev does not assign before returning keep
    // the previous try's values in SIMC (checks.txt of the capture): not compared.
    {
        auto P = makeProcess(inp);
        CaptureCompare cc{"complete_ev", tol};
        const std::set<std::string> staleOnFail = {
            "v_Pmx", "v_Pmy", "v_Pmz", "v_Pmiss", "v_Emiss", "v_PmPar", "v_PmOop", "v_PmPer", "v_Trec", "v_p_E",
            "v_p_P", "v_p_delta", "t_Eloss_2", "t_Eloss_3", "t_teff_2", "t_teff_3", "m_Trec", "m_w", "m_t", "m_tmin",
            "m_q2", "m_epsilon", "m_theta_pq", "m_phi_pq", "ntup_krel", "v_up_x", "v_up_y", "v_up_z", "r_hardcorfac",
            "r_frac_1", "r_frac_2", "r_frac_3", "r_g_4", "r_c_4", "r_rad_proton_ev",
            "v_zhad", "v_pt2"};   // semi: assigned after the two-pion threshold test (event.f:990-992)
        for (std::size_t i = 0; i < cevT.rows.size(); ++i) {
            Row r(cevT, i);
            MainEvent m;
            Event v;
            if (isRho())
                loadCompleteEvInputsRho(*P, r, m, v, false);
            else
                loadCompleteEvInputs(*P, r, m, v, false);
            std::vector<double> g, u;
            for (int k = 1; k <= iv(r, "ng"); ++k) g.push_back(r["g" + std::to_string(k)]);
            for (int k = 1; k <= 4; ++k) u.push_back(r["u" + std::to_string(k)]);
            ScriptedRandom rng;
            rng.load(g, u);
            P->reaction->beginTry(std::lround(r["ntried"]));
            Outputs o;
            const bool ok = P->reaction->completeEv(m, v, rng);
            putCompleteEv(*P, ok, m, v, rng, o);
            if (isRho()) putRhoGnrl(m, o);
            std::set<std::string> skipFail = staleOnFail;
            if (isRho() && !ok)   // checks.txt: v_*, r_* and the stale-prone fields of failed complete_ev calls
                for (const auto& n : cevT.names)
                    if (n.rfind("v_", 0) == 0 || n.rfind("r_", 0) == 0 || n == "ntup_rhomass" || n == "Mh" ||
                        n == "Mh2" || n == "Mh2_final" || n == "jacobian")
                        skipFail.insert(n);
            cc.row(r, cevT, o, {"ntried", "in_v_Pm", "n_unif_ph1", "gen_weight"}, ok ? nullptr : &skipFail);
        }
        fail += cc.report() ? 0 : 1;
    }
    // complete_recon_ev and complete_main (deterministic)
    {
        auto P = makeProcess(inp);
        CaptureCompare cr{"complete_recon_ev", tol}, cm{"complete_main", tol};
        for (std::size_t i = 0; i < recT.rows.size(); ++i) {
            Row r(recT, i);
            Event rec;
            getArm("in_e_", r, rec.e);
            getArm("in_p_", r, rec.p);
            if (isRho())   // COMMON Mh2: Mrho2 for calculate_central (ntried 0), the decay pion's Mpi2 for every event
                dynamic_cast<reactions::RhoProduction&>(*P->reaction)
                    .setReconMh2(std::lround(r["ntried"]) == 0 ? P->setup.Mh2 : constants::Mpi2);
            Outputs o;
            o["success"] = P->reaction->completeRecon(rec);
            putEvent("r_", rec, o);
            cr.row(r, recT, o, {"ntried"});
        }
        for (std::size_t i = 0; i < mainT.rows.size(); ++i) {
            Row r(mainT, i);
            Outputs o;
            doCompleteMain(*P, r, o);
            cm.row(r, mainT, o, {"ntried"});
        }
        fail += cr.report() ? 0 : 1;
        fail += cm.report() ? 0 : 1;
    }
    std::printf("  => %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 7) {
        std::fprintf(stderr, "usage: %s pion_state|generate|complete_ev|complete_recon_ev|complete_main|capture <dir> "
                             "<inp dir|inp> <data/spectral> <data/pion> <tol>\n", argv[0]);
        return 2;
    }
    const std::string what = argv[1];
    gSpectral = argv[4];
    gPion = argv[5];
    if (argc > 8 && std::string(argv[7]) == "kaon") gKaon = argv[8];
    if (argc > 7 && std::string(argv[7]) == "delta") gDelta = true;
    if (argc > 8 && std::string(argv[7]) == "semi") gSidis = argv[8];
    if (argc > 7 && std::string(argv[7]) == "rho") gRho = true;
    const double tol = std::atof(argv[6]);
    try {
        if (what == "pion_state") return runPionState(argv[2], argv[3], tol);
        if (what == "generate" || what == "complete_ev" || what == "complete_recon_ev" || what == "complete_main" ||
            what == "peedelta")
            return runEventTable(what, argv[2], argv[3], tol);
        if (what == "capture") return runCapture(argv[2], argv[3], tol);
    } catch (const std::exception& e) {
        std::printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "unknown mode %s\n", what.c_str());
    return 2;
}
