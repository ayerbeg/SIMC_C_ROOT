// tests/l1/test_l1_eep.cpp
//
// L1 tests of D(e,e'p) / A(e,e'p) (tasks P8.1, P8.2) against the P8.F1 Fortran tables (reference/functions/, see
// reference/P8.F1_report.md). Everything in fortran_compat mode (D-031) and bitwise (tolerance argument, 0 = bitwise).
//
//   files <data/spectral> <fortran dir>              the data/spectral copies are byte-identical to the Fortran files
//   theory_init <dir> <inp dir> <data dir> <tol>     init.f:theory_init through BuildRunSetup (shells + theory.txt)
//   sf_lookup_init <dir> <inp dir> <data dir> <tol>  sf_lookup_init per file (per_file.txt) and after the Benhar set-ups
//   sf_lookup <dir> <data dir> <tol>                 sf_lookup + sf_lookup_diff, rows in order (SAVEd state)
//   generate_em <dir> <data dir> <tol>               generate_em with scripted grnd
//   eep_state <dir> <inp dir> <data dir> <tol>       the whole start-up state per configuration (dbase/init/central)
//   generate <dir> <inp dir> <data dir> <tol>        event.f:generate (D, heavy), scripted grnd/gauss1
//   complete_ev <dir> <inp dir> <data dir> <tol>     event.f:complete_ev (D, heavy), scripted gauss1
//   complete_recon_ev <dir> <inp dir> <data dir> <tol>
//   complete_main <dir> <inp dir> <data dir> <tol>   event.f:complete_main incl. deForest/Coulomb overrides
//   capture <cap dir> <inp> <data dir> <tol>         replay of the real calls of a capture run (generate with the
//                                                    captured stream and persistent state, complete_ev, recon, main)
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
#include "simc/core/Registry.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/driver/SimcRun.h"
#include "simc/input/InpReader.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/radiation/SimcRadiator.h"
#include "simc/radiation/kernel/RadcState.h"
#include "simc/reactions/EventChain.h"
#include "simc/reactions/NuclearEep.h"
#include "simc/reactions/StandardBeam.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/physics/nuclear/SpectralFunctionTables.h"
#include "simc/physics/nuclear/SpectralFunctions.h"
#include "simc/setup/SetupBuilder.h"

using namespace simc;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

int iv(const Row& r, const char* c) { return static_cast<int>(std::lround(r[c])); }

bool same(double a, double b) { return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof a) == 0; }

/// Bitwise / tolerance bookkeeping of free-form comparisons.
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

/// cfg -> SIMC input file (reference/functions/eep_state/configs.txt)
std::string cfgInp(const std::string& inpdir, int cfg) {
    static const char* f[] = {"",
                              "p8f1_d_rad_cap",
                              "p8f1_fe56_rad_cap",
                              "l1_d_norad",
                              "l1_fe56_norad",
                              "l1_c12_rad",
                              "l1_au197_rad",
                              "l1_he4_rad",
                              "l1_fe56_bh_rad",
                              "l1_he3_rad",
                              "l1_ca40_rad",
                              "l1_he3_bh_rad",
                              "l1_he4_bh_rad",
                              "l1_c12_bh_rad",
                              "l1_au197_bh_rad"};
    if (cfg < 1 || cfg > 14) throw std::runtime_error("bad cfg " + std::to_string(cfg));
    return inpdir + "/" + f[cfg] + ".inp";
}

RunSetup setupFor(const std::string& inpdir, int cfg, const std::string& data_dir, SimcParameters* pout = nullptr) {
    SimcParameters p;
    InpReader rd(p);
    const std::string inp = cfgInp(inpdir, cfg);
    if (!rd.loadRun(inp)) throw std::runtime_error("cannot read " + inp);
    SetupOptions opt;
    opt.fortran_compat = true;
    opt.data_dir = data_dir;
    RunSetup s = BuildRunSetup(p, opt);
    if (!s.ok()) throw std::runtime_error("cfg " + std::to_string(cfg) + ": setup STOP " + s.stops.front());
    if (pout) *pout = p;
    return s;
}

// ------------------------------------------------------------------------------------------------- files
int runFiles(const std::string& dataDir, const std::string& fortranDir) {
    std::ifstream sums(dataDir + "/MD5SUMS");
    if (!sums) return 77;
    int bad = 0, n = 0;
    std::string md5, name;
    while (sums >> md5 >> name) {
        ++n;
        auto slurp = [](const std::string& p, std::string& out) {
            std::ifstream f(p, std::ios::binary);
            if (!f) return false;
            std::ostringstream s;
            s << f.rdbuf();
            out = s.str();
            return true;
        };
        std::string a, b;
        if (!slurp(dataDir + "/" + name, a)) {
            std::printf("  missing copy %s\n", name.c_str());
            ++bad;
            continue;
        }
        if (!slurp(fortranDir + "/" + name, b)) {
            std::printf("  %s: Fortran file missing (skip)\n", name.c_str());
            continue;
        }
        const bool eq = a == b;
        std::printf("  %-20s %s\n", name.c_str(), eq ? "identical to the Fortran file" : "DIFFERS from the Fortran file");
        if (!eq) ++bad;
    }
    std::printf("spectral data files: %d listed, %d bad -> %s\n", n, bad, bad || !n ? "FAIL" : "PASS");
    return (bad || !n) ? 1 : 0;
}

// ------------------------------------------------------------------------------------------------- theory_init
int runTheoryInit(const std::string& dir, const std::string& inpdir, const std::string& dataDir, double tol) {
    simcl1::Table out, th;
    std::string err;
    if (!out.load(dir + "/output.txt", err) || !th.load(dir + "/theory.txt", err)) {
        std::printf("theory_init: SKIP (%s)\n", err.c_str());
        return 77;
    }
    std::map<int, RunSetup> setups;
    auto S = [&](int cfg) -> const RunSetup& {
        auto it = setups.find(cfg);
        if (it == setups.end()) it = setups.emplace(cfg, setupFor(inpdir, cfg, dataDir)).first;
        return it->second;
    };
    Cmp c{"theory_init shells (output.txt)", tol};
    for (std::size_t i = 0; i < out.rows.size(); ++i) {
        Row r(out, i);
        const RunSetup& s = S(iv(r, "cfg"));
        const TheoryData& t = s.theory;
        const std::size_t m = static_cast<std::size_t>(iv(r, "m") - 1);
        c.check("nrhoPm", t.nrhoPm, r["nrhoPm"], r.number());
        c.check("E_Fermi", t.E_Fermi, r["E_Fermi"], r.number());
        if (m >= t.nprot.size()) {
            c.check("m", -1, r["m"], r.number());
            continue;
        }
        c.check("nprot_theory", t.nprot[m], r["nprot_theory"], r.number());
        c.check("Em_theory", t.Em[m], r["Em_theory"], r.number());
        c.check("Emsig_theory", t.Emsig[m], r["Emsig_theory"], r.number());
        c.check("bs_norm_theory", t.bs_norm[m], r["bs_norm_theory"], r.number());
        c.check("Em_int_theory", t.Em_int[m], r["Em_int_theory"], r.number());
        c.check("Pm_n", t.Pm[m].n, r["Pm_n"], r.number());
        c.check("Pm_bin", t.Pm[m].bin, r["Pm_bin"], r.number());
        c.check("Pm_min", t.Pm[m].min, r["Pm_min"], r.number());
        c.check("Pm_max", t.Pm[m].max, r["Pm_max"], r.number());
        c.check("doing_deuterium", s.flags.doing_deuterium ? 1. : 0., r["doing_deuterium"], r.number());
    }
    Cmp ct{"theory_init theory(m,n) (theory.txt)", tol};
    for (std::size_t i = 0; i < th.rows.size(); ++i) {
        Row r(th, i);
        const TheoryData& t = S(iv(r, "cfg")).theory;
        const std::size_t m = static_cast<std::size_t>(iv(r, "m") - 1);
        const std::size_t n = static_cast<std::size_t>(iv(r, "n"));
        ct.check("Pm_center", t.Pm[m].min + (static_cast<double>(n) - 0.5) * t.Pm[m].bin, r["Pm_center"], r.number());
        ct.check("theory", t.theory[m][n], r["theory"], r.number());
    }
    // the theory file actually chosen (D-033 F1: A not in {2,12,56,197} -> c12.theory, with the explicit warning)
    int fb = 0;
    for (const auto& kv : setups) {
        const RunSetup& s = kv.second;
        std::printf("    cfg %2d: A = %5.1f -> %-14s%s\n", kv.first, s.targ.A, s.sf.file.c_str(),
                    s.sf.fallback ? "  (Fortran fallback, warning printed)" : "");
        const int A = static_cast<int>(std::lround(s.targ.A));
        if (s.sf.fallback != (A != 2 && A != 12 && A != 56 && A != 197)) ++fb;
    }
    const bool ok = c.report() && ct.report() && fb == 0;
    std::printf("theory_init: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

// ------------------------------------------------------------------------------------------------- sf_lookup_init
const char* kSfFile[] = {"", "benharsf_3mod.dat", "benharsf_4.dat", "benharsf_12.dat", "benharsf_56.dat",
                         "benharsf_197.dat"};

int runSfLookupInit(const std::string& dir, const std::string& inpdir, const std::string& dataDir, double tol) {
    simcl1::Table pf, out;
    std::string err;
    if (!pf.load(dir + "/per_file.txt", err) || !out.load(dir + "/output.txt", err)) {
        std::printf("sf_lookup_init: SKIP (%s)\n", err.c_str());
        return 77;
    }
    // direct sf_lookup_init per (file, proton flag)
    std::map<std::pair<int, int>, std::unique_ptr<nuclear::SfTable>> tabs;
    Cmp c{"sf_lookup_init per file (per_file.txt)", tol};
    std::vector<std::string> msgs;
    for (std::size_t i = 0; i < pf.rows.size(); ++i) {
        Row r(pf, i);
        const auto key = std::make_pair(iv(r, "sf_id"), iv(r, "proton"));
        auto& T = tabs[key];
        if (!T) {
            T = std::make_unique<nuclear::SfTable>();
            T->setFortranCompat(true);
            // driver.log holds the messages of the proton-column processes only
            if (!T->load(dataDir + "/" + kSfFile[key.first], key.second == 1, key.second == 1 ? &msgs : nullptr, &err)) {
                std::printf("  %s\n", err.c_str());
                return 1;
            }
        }
        const int iPm = iv(r, "iPm"), iEm = iv(r, "iEm");
        c.check("numPm", T->numPm(), r["numPm"], r.number());
        c.check("numEm", T->numEm(), r["numEm"], r.number());
        c.check("sftotnorm", T->sftotnorm(), r["sftotnorm"], r.number());
        c.check("Pmval", T->Pmval(iPm), r["Pmval"], r.number());
        c.check("dPm", T->dPm(iPm), r["dPm"], r.number());
        c.check("Emval", T->Emval(iEm), r["Emval"], r.number());
        c.check("dEm", T->dEm(iEm), r["dEm"], r.number());
        c.check("sfnorm_iPm", T->sfnorm(iPm), r["sfnorm_iPm"], r.number());
        c.check("sfval", T->sfval(iEm, iPm), r["sfval"], r.number());
    }
    // the normalisation message, as printed by the Fortran ('Uncorrected S.F. has normalization=')
    {
        std::ifstream lf(dir + "/driver.log");
        std::set<std::string> fl;
        for (std::string l; std::getline(lf, l);)
            if (l.find("Uncorrected S.F.") != std::string::npos) fl.insert(l);
        int miss = 0;
        for (const auto& m : msgs)
            if (m.find("Uncorrected") != std::string::npos && !fl.count(m)) {
                std::printf("    message not in driver.log: |%s|\n", m.c_str());
                ++miss;
            }
        std::printf("  normalisation messages: %zu C++, all in the Fortran driver.log: %s\n", msgs.size(),
                    miss ? "NO" : "yes");
        if (miss) c.bad["message"] += static_cast<std::size_t>(miss);
    }
    // /sftable/ after the full set-up of the Benhar configurations (the "benhar" plug-in reads the file)
    Cmp cs{"sf_lookup_init after set-up (output.txt)", tol};
    std::map<int, std::unique_ptr<nuclear::BenharSpectralFunction>> sfs;
    std::map<int, RunSetup> setups;
    for (std::size_t i = 0; i < out.rows.size(); ++i) {
        Row r(out, i);
        const int cfg = iv(r, "cfg");
        auto& B = sfs[cfg];
        if (!B) {
            SimcParameters p;
            setups[cfg] = setupFor(inpdir, cfg, dataDir, &p);
            B = std::make_unique<nuclear::BenharSpectralFunction>();
            B->setLog(nullptr);
            B->configure(Json::object(), p);
            B->init(setups[cfg]);
        }
        nuclear::SfTable& T = B->table();
        const int iPm = iv(r, "iPm"), iEm = iv(r, "iEm");
        cs.check("numPm", T.numPm(), r["numPm"], r.number());
        cs.check("numEm", T.numEm(), r["numEm"], r.number());
        cs.check("sftotnorm", T.sftotnorm(), r["sftotnorm"], r.number());
        cs.check("Pmval", T.Pmval(iPm), r["Pmval"], r.number());
        cs.check("dPm", T.dPm(iPm), r["dPm"], r.number());
        cs.check("Emval", T.Emval(iEm), r["Emval"], r.number());
        cs.check("dEm", T.dEm(iEm), r["dEm"], r.number());
        cs.check("sfval", T.sfval(iEm, iPm), r["sfval"], r.number());
    }
    for (const auto& kv : sfs)
        std::printf("    cfg %2d: A = %5.1f -> %s (transparency %.17g)\n", kv.first, setups[kv.first].targ.A,
                    kv.second->file().c_str(), kv.second->transparency());
    const bool ok = c.report() && cs.report();
    std::printf("sf_lookup_init: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

// ------------------------------------------------------------------------------------------------- sf_lookup
int runSfLookup(const std::string& dir, const std::string& dataDir, double tol) {
    std::map<std::pair<int, int>, std::unique_ptr<nuclear::SfTable>> tabs;   // one "process" per (file, proton)
    std::string err;
    return simcl1::run("sf_lookup + sf_lookup_diff (sf_lookup.f:85-190, rows in order)", dir,
                       {{"SF", tol}, {"SFd", tol}}, [&](const Row& r, Outputs& o) {
                           const auto key = std::make_pair(iv(r, "sf_id"), iv(r, "proton"));
                           auto& T = tabs[key];
                           if (!T) {
                               T = std::make_unique<nuclear::SfTable>();
                               T->setFortranCompat(true);
                               if (!T->load(dataDir + "/" + kSfFile[key.first], key.second == 1, nullptr, &err))
                                   throw std::runtime_error(err);
                           }
                           o["SF"] = T->lookup(r["Em"], r["Pm"]);
                           o["SFd"] = T->lookupDiff(r["Em"], r["Pm"]);
                       });
}

// ------------------------------------------------------------------------------------------------- generate_em
int runGenerateEm(const std::string& dir, const std::string& dataDir, double tol) {
    std::map<std::pair<int, int>, std::unique_ptr<nuclear::SfTable>> tabs;
    std::string err;
    return simcl1::run("generate_em (sf_lookup.f:195-245, scripted grnd)", dir,
                       {{"Em", tol}, {"n_unif", 0}, {"rng_err", 0}}, [&](const Row& r, Outputs& o) {
                           const auto key = std::make_pair(iv(r, "sf_id"), iv(r, "proton"));
                           auto& T = tabs[key];
                           if (!T) {
                               T = std::make_unique<nuclear::SfTable>();
                               T->setFortranCompat(true);
                               if (!T->load(dataDir + "/" + kSfFile[key.first], key.second == 1, nullptr, &err))
                                   throw std::runtime_error(err);
                           }
                           ScriptedRandom rng;
                           rng.load({}, {r["u1"], r["u2"]});
                           o["Em"] = T->generateEm(r["Pm"], rng);
                           o["n_unif"] = static_cast<double>(rng.nUniform());
                           o["rng_err"] = rng.error() ? 1. : 0.;
                       });
}


// =====================================================================================================================
// Event-level tables (P8.2)
// =====================================================================================================================
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

/// IRadiator decorator: records whether generate_rad was called (the table's rad_called / fail_stage).
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

/// One Fortran process (es_init): set-up, target, radiator, beam, the deForest models, the reaction.
struct Process {
    SimcParameters params;
    RunSetup setup;
    material::HallCCryoTarget target;
    RecordingRadiator rad;
    reactions::StandardBeam beam;
    std::unique_ptr<IPhysicsModel> model;                  // "deforest" with the input's deForest_flag
    std::map<int, std::unique_ptr<IPhysicsModel>> fixed;   // deForest_flag -1, 0, 1 (harness overrides)
    std::unique_ptr<reactions::NuclearEep> reaction;
    reactions::MainContext ctx;
};

std::unique_ptr<Process> makeProcess(const std::string& inp, const std::string& dataDir, bool compat = true) {
    auto P = std::make_unique<Process>();
    InpReader rd(P->params);
    if (!rd.loadRun(inp)) throw std::runtime_error("cannot read " + inp);
    SetupOptions opt;
    opt.fortran_compat = compat;
    opt.data_dir = dataDir;
    P->setup = BuildRunSetup(P->params, opt);
    if (!P->setup.ok()) throw std::runtime_error(inp + ": setup STOP " + P->setup.stops.front());
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
    const Json mc = Json{{"Mh", P->setup.Mh}, {"verbose", false}};
    P->model = Registry<IPhysicsModel>::create("deforest");
    P->model->configure(mc, P->params);
    P->model->init(sc);
    const char* names[3] = {"deforest_cc1_onshell", "deforest_cc1", "deforest_cc2"};
    for (int f = -1; f <= 1; ++f) {
        auto& m = P->fixed[f];
        m = Registry<IPhysicsModel>::create(names[f + 1]);
        m->configure(mc, P->params);
        m->init(sc);
    }
    sc.physics = P->model.get();
    P->reaction = std::make_unique<reactions::NuclearEep>(P->setup.flags.doing_deuterium ? "deuterium_eep" : "nucleus_eep");
    P->reaction->configure(Json{{"verbose", false}}, P->params);
    P->reaction->setLog(nullptr);
    P->reaction->init(sc);
    P->ctx = P->reaction->mainContext();
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

/// One generate call as the SIMC main loop does it (simc.f:174-183 resets, ntried; event.f:generate).
bool runGenerate(Process& P, long ntried, MainEvent& m, Event& v, Event& o, ScriptedRandom& rng, Outputs& out) {
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
    P.reaction->beginTry(ntried);
    const bool ok = reactions::generateEvent(P.beam, *P.reaction, m, v, o, rng);
    out["success"] = ok;
    out["fail_stage"] = ok ? 0 : (P.rad.called ? 2 : 1);
    out["rad_called"] = P.rad.called;
    out["n_unif"] = static_cast<double>(rng.nUniform());
    out["n_gauss"] = static_cast<double>(rng.nGauss());
    out["rng_err"] = rng.error();
    out["radphot"] = m.ntup.radphot;
    out["radarm"] = m.ntup.radarm;
    out["ntup_krel"] = m.ntup.krel;
    putMain(m, out);
    putTarget(m.target, out);
    putEvent("v_", v, out);
    putEvent("o_", o, out);
    putRad(P.rad.inner.state(), out);
    return ok;
}

int runEventTable(const std::string& name, const std::string& dir, const std::string& inpdir,
                  const std::string& dataDir, double tol) {
    std::map<int, std::unique_ptr<Process>> procs;   // one Fortran process per cfg
    auto proc = [&](int cfg) -> Process& {
        auto& P = procs[cfg];
        if (!P) P = makeProcess(cfgInp(inpdir, cfg), dataDir);
        return *P;
    };
    if (name == "generate") {
        return runAuto("generate_eepA_fixedrng (event.f:126-428, D/heavy, scripted grnd/gauss1)", dir, tol,
                       [&](const Row& r, Outputs& o) {
                           Process& P = proc(iv(r, "cfg"));
                           MainEvent m;
                           Event v, orig;
                           v.Pm = r["in_v_Pm"];   // stale vertex%Pm of the previous try (D-033 F5)
                           for (int i = 1; i <= 3; ++i)
                               P.rad.inner.state().Egamma_min[static_cast<std::size_t>(i)] =
                                   P.rad.inner.state().Egamma_max[static_cast<std::size_t>(i)] = 0.;
                           std::vector<double> u, g;
                           for (int i = 1; i <= 16; ++i) u.push_back(r["u" + std::to_string(i)]);
                           for (int i = 1; i <= 32; ++i) g.push_back(r["g" + std::to_string(i)]);
                           ScriptedRandom rng;
                           rng.load(g, u);
                           runGenerate(P, iv(r, "gen_ntried"), m, v, orig, rng, o);
                       });
    }
    if (name == "complete_ev") {
        return runAuto("complete_ev_eepA (event.f:432-1052, D/heavy)", dir, tol, [&](const Row& r, Outputs& o) {
            Process& P = proc(iv(r, "cfg"));
            const RunSetup& s = P.setup;
            MainEvent m;
            Event v;
            v.Ein = r["Ein"];
            v.p.xptar = r["p_xptar"];
            v.p.yptar = r["p_yptar"];
            v.p.theta = r["p_theta"];
            v.p.phi = r["p_phi"];
            v.p.E = r["p_E"];
            v.p.P = r["p_P"];
            v.p.delta = r["p_delta"];
            v.e.xptar = r["e_xptar"];
            v.e.yptar = r["e_yptar"];
            v.e.theta = r["e_theta"];
            v.e.phi = r["e_phi"];
            v.e.E = r["e_E"];
            v.e.P = r["e_P"];
            v.e.delta = r["e_delta"];
            v.Pm = r["in_v_Pm"];
            if (iv(r, "calc_angles") == 1) {
                PhysicsAngles(s.spec.e.theta, s.spec.e.phi, v.e.xptar, v.e.yptar, v.e.theta, v.e.phi);
                PhysicsAngles(s.spec.p.theta, s.spec.p.phi, v.p.xptar, v.p.yptar, v.p.theta, v.p.phi);
            }
            o["in_e_theta"] = v.e.theta;
            o["in_e_phi"] = v.e.phi;
            o["in_p_theta"] = v.p.theta;
            o["in_p_phi"] = v.p.phi;
            m.target.z = r["t_z"];
            m.target.teff[0] = r["t_teff_1"];
            std::vector<double> g;
            for (int i = 1; i <= 12; ++i) g.push_back(r["g" + std::to_string(i)]);
            ScriptedRandom rng;
            rng.load(g, {});
            P.reaction->beginTry(iv(r, "gen_ntried"));
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
    }
    if (name == "complete_recon_ev") {
        return runAuto("complete_recon_ev_eepA (event.f:1056-1359)", dir, tol, [&](const Row& r, Outputs& o) {
            Process& P = proc(iv(r, "cfg"));
            Event rec;
            getArm("in_e_", r, rec.e);
            getArm("in_p_", r, rec.p);
            o["success"] = P.reaction->completeRecon(rec);
            putEvent("r_", rec, o);
            o["ntup_mm"] = P.reaction->reconNtup().mm;
            o["ntup_mmA"] = P.reaction->reconNtup().mmA;
            o["ntup_t"] = P.reaction->reconNtup().t;
        });
    }
    if (name == "complete_main") {
        std::size_t identity_bad = 0;
        const int rc = runAuto("complete_main_eepA (event.f:1363-1569)", dir, tol, [&](const Row& r, Outputs& o) {
            Process& P = proc(iv(r, "cfg"));
            MainEvent m;
            Event v, rec, v0;
            getEvent("v_", r, v);
            getEvent("r_", r, rec);
            m.jacobian = r["in_jacobian"];
            m.SF_weight = r["in_SF_weight"];
            m.sigcc = r["in_sigcc"];
            m.sigcc_recon = r["in_sigcc_recon"];
            m.weight = r["in_weight"];
            m.gen_weight = r["in_gen_weight"];
            // harness overrides of the COMMON state for this row
            reactions::MainContext ctx = P.ctx;
            IPhysicsModel* model = P.model.get();
            int flag = P.params.deForest_flag;
            if (iv(r, "ov_deforest") == 1) {
                flag = iv(r, "ov_deforest_flag");
                model = P.fixed.at(flag).get();
            }
            if (iv(r, "ov_coulomb") == 1) {
                ctx.using_Coulomb = iv(r, "ov_using_coulomb") != 0;
                ctx.Coulomb_ave = r["ov_coulomb_ave"];
            }
            model->beginEvent(iv(r, "gen_ntried"));
            P.reaction->mainContext() = ctx;
            P.reaction->setPhysics(model);
            const bool ok = P.reaction->completeMain(iv(r, "force_sigcc") == 1, m, v, v0, rec);
            P.reaction->mainContext() = P.ctx;
            P.reaction->setPhysics(P.model.get());
            o["success"] = ok;
            o["SF_weight"] = m.SF_weight;
            o["sigcc"] = m.sigcc;
            o["sigcc_recon"] = m.sigcc_recon;
            o["jacobian"] = m.jacobian;
            o["gen_weight"] = m.gen_weight;
            o["weight"] = m.weight;
            o["eff_deForest_flag"] = flag;
            o["eff_using_Coulomb"] = ctx.using_Coulomb ? 1. : 0.;
            o["eff_coulomb_ave"] = ctx.Coulomb_ave;
            const double c = 1.0 + ctx.Coulomb_ave / ctx.Ebeam;
            o["hx_coul_fac"] = ctx.using_Coulomb ? c * c : 1.0;
            o["hx_prod"] = m.SF_weight * m.jacobian * m.gen_weight * m.sigcc;
            o["hx_deforest_v"] = model->sigcc(v, m);
            o["hx_deforest_r"] = model->sigccRecon(rec, m);
            if (ok && !same(o["hx_prod"], m.weight)) ++identity_bad;
        });
        std::printf("  product identity weight == ((SF_weight*jacobian)*gen_weight)*sigcc violated in %zu rows\n",
                    identity_bad);
        return identity_bad ? 1 : rc;
    }
    std::printf("unknown table %s\n", name.c_str());
    return 1;
}

// ------------------------------------------------------------------------------------------------- eep_state
int runEepState(const std::string& dir, const std::string& inpdir, const std::string& dataDir, double tol) {
    simcl1::Table out;
    std::string err;
    if (!out.load(dir + "/output.txt", err)) {
        std::printf("eep_state: SKIP (%s)\n", err.c_str());
        return 77;
    }
    Cmp c{"eep_state (start-up + calculate_central)", tol};
    for (std::size_t i = 0; i < out.rows.size(); ++i) {
        Row r(out, i);
        const int cfg = iv(r, "cfg");
        SimcParameters p;
        InpReader rd(p);
        if (!rd.loadRun(cfgInp(inpdir, cfg))) throw std::runtime_error("cannot read cfg " + std::to_string(cfg));
        driver::RunOptions opt;
        opt.fortran_compat = true;
        opt.data_dir = dataDir;
        opt.quiet_stages = true;
        opt.run_name = "eep_state_cfg" + std::to_string(cfg);
        std::ostringstream log;
        opt.log = &log;
        driver::SimcRun run(p, Json{{"sinks", Json::array()}}, opt);
        run.initialise();
        const RunSetup& s = run.setup();
        const TargetInfo& t = s.targ;
        const ReactionFlags& f = s.flags;
        Outputs o;
        o["doing_eep"] = f.doing_eep;
        o["doing_hyd_elast"] = f.doing_hyd_elast;
        o["doing_deuterium"] = f.doing_deuterium;
        o["doing_heavy"] = f.doing_heavy;
        o["use_benhar_sf"] = s.params.use_benhar_sf;
        o["transparency"] = s.params.transparency;
        o["using_rad"] = s.sw.using_rad;
        o["using_Eloss"] = s.sw.using_Eloss;
        o["using_Coulomb"] = s.sw.using_Coulomb;
        o["deForest_flag"] = s.params.deForest_flag;
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
        o["targ_angle"] = t.angle;
        o["targ_X0"] = t.X0;
        o["targ_X0_cm"] = t.X0_cm;
        o["targ_L1"] = t.L1;
        o["targ_L2"] = t.L2;
        o["targ_Coulomb_ave"] = t.Coulomb.ave;
        o["targ_Coulomb_min"] = t.Coulomb.min;
        o["targ_Coulomb_max"] = t.Coulomb.max;
        o["targ_Coulomb_constant"] = t.Coulomb_constant;
        for (int k = 0; k < 3; ++k) {
            const std::string ks = std::to_string(k + 1);
            const std::size_t kk = static_cast<std::size_t>(k);
            o["targ_Eloss_ave_" + ks] = t.Eloss[kk].ave;
            o["targ_Eloss_min_" + ks] = t.Eloss[kk].min;
            o["targ_Eloss_max_" + ks] = t.Eloss[kk].max;
            o["targ_teff_ave_" + ks] = t.teff[kk].ave;
            o["targ_musc_max_" + ks] = t.musc_max[kk];
        }
        o["Mh"] = s.Mh;
        o["Mh2"] = s.Mh2;
        o["Ebeam"] = s.Ebeam;
        o["Ebeam_vertex_ave"] = s.Ebeam_vertex_ave;
        o["E_Fermi"] = s.theory.E_Fermi;
        o["nrhoPm"] = s.theory.nrhoPm;
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
        cut("edge_Em", s.edge.Em);
        cut("edge_e_E", s.edge.e.E);
        cut("edge_p_E", s.edge.p.E);
        o["slop_total_Em_used"] = s.slop.total.Em.used;
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
        o["c_rad_hardcorfac"] = C.rad.hardcorfac;
        putEvent("v0_", run.central().vertex0, o);
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
    const bool ok = c.report();
    std::printf("eep_state: %s\n", ok ? "PASS" : "FAIL");
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
             const std::set<std::string>* only = nullptr) {
        ++rows;
        bool any = false;
        for (const auto& n : t.names) {
            if (exclude.count(n)) continue;
            if (only && !only->count(n)) continue;
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

int runCapture(const std::string& capdir, const std::string& inp, const std::string& dataDir, double tol) {
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
            if (ph != 1 && ph != 6) continue;
            auto& q = stream[std::lround(r[static_cast<std::size_t>(cN)])];
            (static_cast<int>(r[static_cast<std::size_t>(cK)]) == 1 ? q.first : q.second)
                .push_back(r[static_cast<std::size_t>(cV)]);
        }
    }
    int fail = 0;
    // generate: one persistent process and persistent main/vertex/orig, as the SIMC main loop
    {
        auto P = makeProcess(inp, dataDir);
        MainEvent m;
        Event v, orig;
        CaptureCompare cc{"generate", tol};
        std::size_t leftover = 0;
        const std::set<std::string> excl = {"ntried", "nevent", "n_unif_rad", "n_gauss_rad", "n_int"};
        for (std::size_t i = 0; i < genT.rows.size(); ++i) {
            Row r(genT, i);
            const long nt = std::lround(r["ntried"]);
            const auto& q = stream[nt];
            ScriptedRandom rng;
            rng.load(q.second, q.first);
            Outputs o;
            o["in_v_Pm"] = v.Pm;   // the stale vertex%Pm at entry (D-033 F5)
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
            cc.row(r, genT, o, excl);
        }
        fail += cc.report() ? 0 : 1;
        std::printf("      tries whose captured random numbers were not consumed exactly: %zu\n", leftover);
        if (leftover) ++fail;
    }
    // complete_ev (first call of every try); on failures only the fields complete_ev assigned before returning
    {
        auto P = makeProcess(inp, dataDir);
        CaptureCompare cc{"complete_ev", tol};
        const std::set<std::string> onFail = {"success", "v_nu", "v_Q2", "v_q", "v_xbj", "v_uq_x", "v_uq_y", "v_uq_z",
                                              "v_ue_x", "v_ue_y", "v_ue_z", "v_up_x", "v_up_y", "v_up_z", "v_Em",
                                              "v_Mrec", "jacobian"};
        for (std::size_t i = 0; i < cevT.rows.size(); ++i) {
            Row r(cevT, i);
            MainEvent m;
            Event v;
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
            m.target.z = r["t_z"];
            m.target.teff[0] = r["t_teff_1"];
            std::vector<double> g;
            for (int k = 1; k <= iv(r, "ng"); ++k) g.push_back(r["g" + std::to_string(k)]);
            ScriptedRandom rng;
            rng.load(g, {});
            P->reaction->beginTry(std::lround(r["ntried"]));
            Outputs o;
            const bool ok = P->reaction->completeEv(m, v, rng);
            o["success"] = ok;
            putEvent("v_", v, o);
            o["jacobian"] = m.jacobian;
            o["ntup_krel"] = m.ntup.krel;
            o["t_teff_2"] = m.target.teff[1];
            o["t_teff_3"] = m.target.teff[2];
            o["t_Eloss_2"] = m.target.Eloss[1];
            o["t_Eloss_3"] = m.target.Eloss[2];
            putRad(P->rad.inner.state(), o);
            cc.row(r, cevT, o, {"ntried"}, ok ? nullptr : &onFail);
        }
        fail += cc.report() ? 0 : 1;
    }
    // complete_recon_ev and complete_main (deterministic; one process: the sf_lookup state would be shared)
    {
        auto P = makeProcess(inp, dataDir);
        CaptureCompare cr{"complete_recon_ev", tol}, cm{"complete_main", tol};
        for (std::size_t i = 0; i < recT.rows.size(); ++i) {
            Row r(recT, i);
            Event rec;
            getArm("in_e_", r, rec.e);
            getArm("in_p_", r, rec.p);
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
            getEvent("v_", r, v);
            getEvent("r_", r, rec);
            m.jacobian = r["in_jacobian"];
            m.gen_weight = r["in_gen_weight"];
            m.SF_weight = r["in_SF_weight"];
            m.sigcc = r["in_sigcc"];
            m.sigcc_recon = r["in_sigcc_recon"];
            m.weight = r["in_weight"];
            Outputs o;
            const bool ok = P->reaction->completeMain(iv(r, "force_sigcc") == 1, m, v, v0, rec);
            o["success"] = ok;
            o["SF_weight"] = m.SF_weight;
            o["sigcc"] = m.sigcc;
            o["sigcc_recon"] = m.sigcc_recon;
            o["jacobian"] = m.jacobian;
            o["gen_weight"] = m.gen_weight;
            o["weight"] = m.weight;
            const double prod = m.SF_weight * m.jacobian * m.gen_weight * m.sigcc;
            if (ok && !same(prod, m.weight)) ++identity_bad;
            cm.row(r, mainT, o, {"ntried"});
        }
        fail += cr.report() ? 0 : 1;
        fail += cm.report() ? 0 : 1;
        std::printf("      complete_main product identity violated in %zu of %zu rows\n", identity_bad,
                    mainT.rows.size());
        if (identity_bad) ++fail;
    }
    std::printf("  => %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    const std::string what = argv[1];
    try {
        if (what == "files" && argc > 3) return runFiles(argv[2], argv[3]);
        if (what == "theory_init" && argc > 5) return runTheoryInit(argv[2], argv[3], argv[4], std::atof(argv[5]));
        if (what == "sf_lookup_init" && argc > 5) return runSfLookupInit(argv[2], argv[3], argv[4], std::atof(argv[5]));
        if (what == "sf_lookup" && argc > 4) return runSfLookup(argv[2], argv[3], std::atof(argv[4]));
        if (what == "generate_em" && argc > 4) return runGenerateEm(argv[2], argv[3], std::atof(argv[4]));
        if (what == "eep_state" && argc > 5) return runEepState(argv[2], argv[3], argv[4], std::atof(argv[5]));
        if ((what == "generate" || what == "complete_ev" || what == "complete_recon_ev" || what == "complete_main") &&
            argc > 5)
            return runEventTable(what, argv[2], argv[3], argv[4], std::atof(argv[5]));
        if (what == "capture" && argc > 5) return runCapture(argv[2], argv[3], argv[4], std::atof(argv[5]));
    } catch (const std::exception& e) {
        std::printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "usage: %s files|theory_init|sf_lookup_init|sf_lookup|generate_em ...\n", argv[0]);
    return 1;
}
