// tests/l1/test_l1_kaon.cpp
//
// L1 tests of the kaon physics (task P10.1) against the P10.F1 Fortran tables (reference/functions/, see
// reference/P10.F1_report.md). The tolerance argument is |rel| (0 = bitwise); the runner reports the bitwise-equal rows
// per column; NaN/Inf rows must match (NaN == NaN).
//
//   files <data/kaon> <fortran dir>                        md5 of the copies (MD5SUMS) + byte equality with the Fortran
//   saghai_tables <dir> <data/kaon>                        COMMON /sigkaon/ after dbase.f:624-679 (fortran_compat
//                                                          reader), all 4 arrays, bitwise
//   sig_factorized <dir> <tol>                             physics_kaon.f:173-237, rows in order (SAVEd fact_w)
//   eekeek|eekeeks <dir> <data/kaon> <tol> <fint 0|1>      physics_kaon.f:241-477; fint 0 = as built (compat, the
//                                                          eekeek/eekeeks tables), 1 = the fint value (D-043 F1
//                                                          corrected, the *_fintabi tables; F3 grid kept)
//   phspwght <dir> <data/kaon> <tol> <fint 0|1>            physics_kaon.f:479-544 (+ the HMS/SHMS STOP of stops.txt)
//   peek <dir> <kaon_state dir> <data/kaon> <tol> <fint>   physics_kaon.f:1-168 with the per-configuration state and
//                                                          the SAVEd zaero (rows in order per cfg)
//   plugin <peek dir> <kaon_state dir> <data/kaon>         the IPhysicsModel "kaon" from the Registry on the peek rows
//   fortran_diff_D043 <peek dir> <kaon_state dir> <data/kaon> <eekeek dir> <saghai dir>
//                                                          D-043 F1/F2/F3/F5/F6 in both modes
// Exit: 0 pass, 1 fail, 77 reference missing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "Md5.h"
#include "simc/core/Errors.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/physics/kaon/KaonCrossSections.h"
#include "simc/physics/kaon/KaonModel.h"
#include "simc/stages/IPhysicsModel.h"

using namespace simc;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

int iv(const Row& r, const char* c) { return static_cast<int>(std::lround(r[c])); }
bool same(double a, double b) { return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof a) == 0; }

std::vector<ColumnSpec> specs(const std::vector<std::string>& names, double tol) {
    std::vector<ColumnSpec> v;
    for (const auto& n : names) v.push_back(ColumnSpec{n, tol, 0.0});
    return v;
}

bool loadTables(const std::string& dataDir, bool corrected, kaon::SaghaiTables& t) {
    std::string err;
    if (!t.load(dataDir, corrected, nullptr, &err)) {
        std::printf("kaon tables: %s\n", err.c_str());
        return false;
    }
    return true;
}

// ------------------------------------------------------------------------------------------------- files
int runFiles(const std::string& dataDir, const std::string& fortranDir) {
    std::ifstream sums(dataDir + "/MD5SUMS");
    if (!sums) return 77;
    int bad = 0, n = 0;
    std::string md5, name;
    while (sums >> md5 >> name) {
        ++n;
        std::string a, b;
        if (!simcl1::slurp(dataDir + "/" + name, a)) {
            std::printf("  missing copy %s\n", name.c_str());
            ++bad;
            continue;
        }
        const std::string got = simcl1::md5hex(a);
        const bool okmd5 = got == md5;
        std::string cmp = "Fortran file missing (not compared)";
        bool okF = true;
        if (simcl1::slurp(fortranDir + "/" + name, b)) {
            okF = a == b;
            cmp = okF ? "identical to the Fortran file" : "DIFFERS from the Fortran file";
        }
        std::printf("  %-18s md5 %s %s, %s\n", name.c_str(), got.c_str(), okmd5 ? "ok" : "MISMATCH", cmp.c_str());
        if (!okmd5 || !okF) ++bad;
    }
    std::printf("kaon data files: %d listed, %d bad -> %s\n", n, bad, bad || n != 4 ? "FAIL" : "PASS");
    return (bad || n != 4) ? 1 : 0;
}

// ------------------------------------------------------------------------------------------------- saghai_tables
int runSaghaiTables(const std::string& dir, const std::string& dataDir) {
    kaon::SaghaiTables t;
    if (!loadTables(dataDir, false, t)) return 1;
    int fail = 0;
    auto cmpFile = [&](const std::string& file, const std::function<double(const Row&, const std::string&)>& get,
                       const std::vector<std::string>& cols) {
        simcl1::Table ref;
        std::string err;
        if (!ref.load(dir + "/" + file, err)) {
            std::printf("  %s: %s\n", file.c_str(), err.c_str());
            ++fail;
            return;
        }
        std::size_t n = 0, eq = 0;
        for (std::size_t i = 0; i < ref.rows.size(); ++i) {
            Row r(ref, i);
            for (const auto& c : cols) {
                ++n;
                const double c1 = get(r, c), f = r[c];
                if (same(c1, f))
                    ++eq;
                else if (n - eq <= 3)
                    std::printf("    %s row %zu %s: ref %.9g cpp %.9g\n", file.c_str(), i + 1, c.c_str(), f, c1);
            }
        }
        std::printf("  %-20s %6zu rows, %7zu values, %7zu bitwise -> %s\n", file.c_str(), ref.rows.size(), n, eq,
                    n == eq && n > 0 ? "ok" : "MISMATCH");
        if (n != eq || n == 0) ++fail;
    };
    cmpFile("weightc.txt", [&](const Row& r, const std::string&) { return t.wc(iv(r, "j"), iv(r, "i")); }, {"weightc"});
    cmpFile("weightd.txt", [&](const Row& r, const std::string&) { return t.wd(iv(r, "k"), iv(r, "j"), iv(r, "i")); },
            {"weightd"});
    std::vector<std::string> pcols, scols;
    for (int k = 1; k <= 6; ++k) {
        pcols.push_back("zrff" + std::to_string(k));
        pcols.push_back("ziff" + std::to_string(k));
        scols.push_back("zsrff" + std::to_string(k));
        scols.push_back("zsiff" + std::to_string(k));
    }
    cmpFile("saghai_proton.txt",
            [&](const Row& r, const std::string& c) {
                const std::size_t i = kaon::SaghaiTables::idxP(iv(r, "iread"), iv(r, "iq2"), iv(r, "iang"));
                const int k = c[c.size() - 1] - '1';
                return static_cast<double>(c[1] == 'r' ? t.zrff[static_cast<std::size_t>(k)][i]
                                                       : t.ziff[static_cast<std::size_t>(k)][i]);
            },
            pcols);
    cmpFile("saghai_sigma0.txt",
            [&](const Row& r, const std::string& c) {
                const std::size_t i = kaon::SaghaiTables::idxS(iv(r, "iread"), iv(r, "iq2"), iv(r, "iang"));
                const int k = c[c.size() - 1] - '1';
                return static_cast<double>(c[2] == 'r' ? t.zsrff[static_cast<std::size_t>(k)][i]
                                                       : t.zsiff[static_cast<std::size_t>(k)][i]);
            },
            scols);
    std::printf("saghai_tables (dbase.f:624-679, fortran_compat reader): %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

// ------------------------------------------------------------------------------------------------- sig_factorized
int runSigFactorized(const std::string& dir, double tol) {
    kaon::SigFactorizedState st;   // one process, rows in order (SAVEd fact_w)
    std::ostringstream log;
    const int rc = simcl1::run("sig_factorized (physics_kaon.f:173-237)", dir, specs({"sig"}, tol),
                               [&](const Row& r, Outputs& o) {
                                   o["sig"] = kaon::sigFactorized(r["q2"], r["w"], r["t"], r["pk"], r["mrec"], st, &log);
                               });
    std::size_t warn = 0;
    for (std::size_t p = log.str().find("does not agree"); p != std::string::npos; p = log.str().find("does not agree", p + 1))
        ++warn;
    std::printf("  'Kaon C.M. momentum passed to sig_factorized does not agree' warnings: %zu\n", warn);
    return rc;
}

// ------------------------------------------------------------------------------------------------- Saghai
int runSaghai(const std::string& which, const std::string& dir, const std::string& dataDir, double tol, bool fintValue) {
    kaon::SaghaiTables t;
    if (!loadTables(dataDir, false, t)) return 1;   // the Fortran tables (the _fintabi references keep F2/F3)
    kaon::SaghaiOptions opt;
    opt.fint_value = fintValue;
    const bool s = which == "eekeeks";
    std::vector<std::string> cols = {"sigma_eep", "ss_out", "q22_out", "angl_out", "hx_px1", "hx_px2", "hx_px3"};
    for (int k = 1; k <= 6; ++k) cols.push_back("hx_zr" + std::to_string(k));
    for (int k = 1; k <= 6; ++k) cols.push_back("hx_zi" + std::to_string(k));
    return simcl1::run(which + (s ? " (physics_kaon.f:357-477)" : " (physics_kaon.f:241-355)") +
                           (fintValue ? ", fint value (D-043 F1 corrected)" : ", as built (fortran_compat)"),
                       dir, specs(cols, tol), [&](const Row& r, Outputs& o) {
                           kaon::SaghaiTrace tr;
                           const double v = s ? kaon::eekeeks(r["ss"], r["q22"], r["angl"], r["theta"], r["phi"], r["epsi"],
                                                              r["mrec_struck"], t, opt, &tr)
                                              : kaon::eekeek(r["ss"], r["q22"], r["angl"], r["theta"], r["phi"], r["epsi"],
                                                             r["mrec_struck"], t, opt, &tr);
                           o["sigma_eep"] = v;
                           o["ss_out"] = r["ss"];
                           o["q22_out"] = r["q22"];
                           o["angl_out"] = r["angl"];
                           for (int k = 0; k < 3; ++k) o["hx_px" + std::to_string(k + 1)] = tr.px[k];
                           for (int k = 0; k < 6; ++k) {
                               o["hx_zr" + std::to_string(k + 1)] = tr.zr[k];
                               o["hx_zi" + std::to_string(k + 1)] = tr.zi[k];
                           }
                       });
}

int runPhspwght(const std::string& dir, const std::string& dataDir, double tol, bool fintValue) {
    kaon::SaghaiTables t;
    if (!loadTables(dataDir, false, t)) return 1;
    kaon::SaghaiOptions opt;
    opt.fint_value = fintValue;
    // hx_fint = the fint value the routine sees: 0 as built, the REAL*4 value with the ABI wrapper (_fintabi)
    int rc = simcl1::run(std::string("phspwght (physics_kaon.f:479-544, dead code)") +
                             (fintValue ? ", fint value (D-043 F1)" : ", as built"),
                         dir, specs({"weight", "hx_fint"}, tol), [&](const Row& r, Outputs& o) {
                             const kaon::PhspwghtResult p = kaon::phspwght(r["delta"], r["theta"], r["phi"],
                                                                           iv(r, "electron_arm"), iv(r, "hadron_arm"), t, opt);
                             o["weight"] = p.weight;
                             o["hx_fint"] = p.fint_seen;   // the fint value the routine sees (0 as built)
                         });
    // HMS/SHMS: three lines and a bare STOP (stops.txt)
    std::ostringstream log;
    bool stopped = false;
    try {
        kaon::phspwght(0., 0., 0., 1, 5, t, opt, &log);
    } catch (const FortranStop&) {
        stopped = true;
    }
    std::ifstream sf(dir + "/stops.txt");
    std::string ref, line;
    while (std::getline(sf, line))
        if (!line.empty() && line[0] == ' ' && line.find("NOTE") == std::string::npos) ref += line + "\n";
    const bool okStop = stopped && (ref.empty() || ref == log.str());
    std::printf("  phspwght(electron_arm 1, hadron_arm 5): %s, message %s\n", stopped ? "STOP" : "no STOP",
                ref.empty() ? "(no reference)" : (ref == log.str() ? "identical to stops.txt" : "DIFFERS"));
    if (!okStop) {
        std::printf("%s---\n%s", log.str().c_str(), ref.c_str());
        rc = 1;
    }
    return rc;
}

// ------------------------------------------------------------------------------------------------- peeK
/// Per-configuration state from kaon_state/output.txt (the Fortran start-up of each configuration).
struct CfgState {
    double Mtar_struck{0}, Mrec_struck{0}, ctau{0}, spec_p_P{0}, Mh2{0};
    int hadron_arm{0};
    bool doing_decay{false};
};
bool loadStates(const std::string& stateDir, std::map<int, CfgState>& st) {
    simcl1::Table t;
    std::string err;
    if (!t.load(stateDir + "/output.txt", err)) {
        std::printf("kaon_state: %s\n", err.c_str());
        return false;
    }
    for (std::size_t i = 0; i < t.rows.size(); ++i) {
        Row r(t, i);
        CfgState c;
        c.Mtar_struck = r["targ_Mtar_struck"];
        c.Mrec_struck = r["targ_Mrec_struck"];
        c.ctau = r["ctau"];
        c.spec_p_P = r["v0_p_P"];   // vertex0%p%P = spec%p%P*(1+0/100)
        c.Mh2 = r["Mh2"];
        c.hadron_arm = iv(r, "hadron_arm");
        c.doing_decay = iv(r, "doing_decay") != 0;
        st[iv(r, "cfg")] = c;
    }
    return !st.empty();
}

kaon::PeekConfig peekConfig(const CfgState& c, bool corrected, bool fintValue) {
    kaon::PeekConfig p;
    p.Mtar_struck = c.Mtar_struck;
    p.Mrec_struck = c.Mrec_struck;
    p.doing_decay = c.doing_decay;
    p.hadron_arm = c.hadron_arm;
    p.ctau = c.ctau;
    p.spec_p_P = c.spec_p_P;
    p.Mh2 = c.Mh2;
    p.survival_p_kaon = corrected;
    p.survival_path_fixed = corrected;
    p.saghai.fint_value = fintValue || corrected;
    p.saghai.proton_grid_fixed = corrected;
    return p;
}

/// the row overrides of drv_peek / drv_complete_main10 (ov_state)
void applyOverrides(const Row& r, kaon::PeekConfig& p) {
    if (r.has("ov_state") && iv(r, "ov_state") == 1) {
        p.doing_decay = iv(r, "ov_doing_decay") != 0;
        p.hadron_arm = iv(r, "ov_hadron_arm");
        p.ctau = r["ov_ctau"];
    }
}

void readEvent(const Row& r, const std::string& pre, Event& v) {
    auto g = [&](const char* n) { return r[pre + n]; };
    v.Ein = g("Ein"); v.Em = g("Em"); v.Pm = g("Pm"); v.Emiss = g("Emiss"); v.Pmiss = g("Pmiss");
    v.Pmx = g("Pmx"); v.Pmy = g("Pmy"); v.Pmz = g("Pmz"); v.PmPar = g("PmPar"); v.PmPer = g("PmPer");
    v.PmOop = g("PmOop"); v.nu = g("nu"); v.q = g("q"); v.Q2 = g("Q2"); v.Trec = g("Trec"); v.W = g("W");
    v.Mrec = g("Mrec"); v.epsilon = g("epsilon"); v.theta_pq = g("theta_pq"); v.theta_tarq = g("theta_tarq");
    v.phi_pq = g("phi_pq"); v.phi_targ = g("phi_targ"); v.beta = g("beta"); v.phi_s = g("phi_s"); v.phi_c = g("phi_c");
    v.zhad = g("zhad"); v.pt2 = g("pt2"); v.xbj = g("xbj");
    auto arm = [&](const std::string& a, ArmFull& x) {
        x.delta = r[pre + a + "_delta"]; x.xptar = r[pre + a + "_xptar"]; x.yptar = r[pre + a + "_yptar"];
        x.z = r[pre + a + "_z"]; x.theta = r[pre + a + "_theta"]; x.phi = r[pre + a + "_phi"];
        x.E = r[pre + a + "_E"]; x.P = r[pre + a + "_P"];
    };
    arm("e", v.e);
    arm("p", v.p);
    v.ue = {g("ue_x"), g("ue_y"), g("ue_z")};
    v.up = {g("up_x"), g("up_y"), g("up_z")};
    v.uq = {g("uq_x"), g("uq_y"), g("uq_z")};
}

void readMain(const Row& r, MainEvent& m) {
    m.epsilon = r["in_m_epsilon"];
    m.t = r["in_m_t"];
    m.theta_pq = r["in_m_theta_pq"];
    m.phi_pq = r["in_m_phi_pq"];
    m.FP.p.path = r["in_m_FP_p_path"];
    m.FP.p.dx = r["in_m_FP_p_dx"];
    m.FP.p.dy = r["in_m_FP_p_dy"];
    m.SP.p.delta = r["in_m_SP_p_delta"];
    m.ntup.decdist = r["decdist_in"];
    m.ntup.sigcm1 = r["in_ntup_sigcm1"];
    m.ntup.sigcm2 = r["in_ntup_sigcm2"];
    m.ntup.sigcm = r["in_ntup_sigcm"];
    m.pfermi = FermiMomentum{r["pfer"], r["pferx"], r["pfery"], r["pferz"], r["efer"]};
}

void putPeek(double sig, double sp, const MainEvent& m, const kaon::PeekConfig& p, Outputs& o) {
    o["peek"] = sig;
    o["survivalprob"] = sp;
    o["decdist_out"] = m.ntup.decdist;
    o["m_thetacm"] = m.thetacm;
    o["m_phicm"] = m.phicm;
    o["m_pcm"] = m.pcm;
    o["m_davejac"] = m.davejac;
    o["m_johnjac"] = m.johnjac;
    o["m_wcm"] = m.wcm;
    o["ntup_sigcm1"] = m.ntup.sigcm1;
    o["ntup_sigcm2"] = m.ntup.sigcm2;
    o["ntup_sigcm"] = m.ntup.sigcm;
    o["eff_doing_decay"] = p.doing_decay ? 1. : 0.;
    o["eff_hadron_arm"] = p.hadron_arm;
    o["eff_ctau"] = p.ctau;
}

const std::vector<std::string> kPeekCols = {"peek",      "survivalprob", "decdist_out", "m_thetacm",   "m_phicm",
                                            "m_pcm",     "m_davejac",    "m_johnjac",   "m_wcm",       "ntup_sigcm1",
                                            "ntup_sigcm2", "ntup_sigcm", "eff_doing_decay", "eff_hadron_arm", "eff_ctau"};

int runPeek(const std::string& dir, const std::string& stateDir, const std::string& dataDir, double tol, bool fintValue) {
    std::map<int, CfgState> st;
    if (!loadStates(stateDir, st)) return 77;
    kaon::SaghaiTables tab;
    if (!loadTables(dataDir, false, tab)) return 1;
    std::map<int, kaon::PeekState> pst;   // one Fortran process per cfg (SAVEd zaero, phi, fact_w)
    std::map<int, std::size_t> perCfg;
    std::ostringstream log;
    const int rc = simcl1::run(std::string("peeK (physics_kaon.f:1-168), fortran_compat") +
                                   (fintValue ? " + fint value (D-043 F1)" : ""),
                               dir, specs(kPeekCols, tol), [&](const Row& r, Outputs& o) {
                                   const int cfg = iv(r, "cfg");
                                   ++perCfg[cfg];
                                   kaon::PeekConfig pc = peekConfig(st.at(cfg), false, fintValue);
                                   applyOverrides(r, pc);
                                   Event v;
                                   MainEvent m;
                                   readEvent(r, "v_", v);
                                   readMain(r, m);
                                   double sp = r["in_survivalprob"];
                                   const double sig = kaon::peeK(v, m, sp, pc, tab, pst[cfg], &log);
                                   putPeek(sig, sp, m, pc, o);
                               });
    std::printf("  rows per cfg:");
    for (const auto& kv : perCfg) std::printf(" %d:%zu", kv.first, kv.second);
    std::printf("\n");
    return rc;
}

// ------------------------------------------------------------------------------------------------- plug-in
std::unique_ptr<IPhysicsModel> makeModel(const CfgState& c, const std::string& dataDir, bool compat) {
    SimcParameters p;
    p.doing_kaon = true;
    p.doing_decay = c.doing_decay;
    p.hadron_arm = c.hadron_arm;
    p.ctau = c.ctau;
    RunSetup s = RunSetup::FromParameters(p);
    s.targ.Mtar_struck = c.Mtar_struck;
    s.targ.Mrec_struck = c.Mrec_struck;
    s.spec.p.P = c.spec_p_P;
    s.Mh2 = c.Mh2;
    s.fortran_compat = compat;
    auto m = Registry<IPhysicsModel>::create("kaon");
    m->configure(Json{{"data_dir", dataDir}, {"verbose", false}}, p);
    m->init(s);
    return m;
}

int runPlugin(const std::string& dir, const std::string& stateDir, const std::string& dataDir) {
    std::map<int, CfgState> st;
    if (!loadStates(stateDir, st)) return 77;
    simcl1::Table in, ref;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !ref.load(dir + "/output.txt", err)) return 77;
    int bad = 0;
    try {
        SimcParameters p;
        RunSetup s = RunSetup::FromParameters(p);
        auto m = Registry<IPhysicsModel>::create("kaon");
        m->configure(Json{{"data_dir", dataDir}}, p);
        m->init(s);
        std::printf("  plugin: init without targ%%Mtar_struck did not throw\n");
        ++bad;
    } catch (const std::invalid_argument& e) {
        std::printf("  plugin: non-kaon set-up refused: %s\n", e.what());
    }
    std::map<int, std::unique_ptr<IPhysicsModel>> models;
    std::size_t rows = 0, eq = 0;
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        Row r(in, i), f(ref, i);
        const int cfg = iv(r, "cfg");
        auto it = models.find(cfg);
        if (it == models.end()) it = models.emplace(cfg, makeModel(st.at(cfg), dataDir, true)).first;
        IPhysicsModel& model = *it->second;
        auto* km = dynamic_cast<kaon::KaonModel*>(&model);
        if (!km) throw std::runtime_error("the \"kaon\" model is not a kaon::KaonModel");
        const kaon::PeekConfig saved = km->config();
        applyOverrides(r, km->config());
        km->setSurvivalProbability(r["in_survivalprob"]);
        Event v;
        MainEvent m;
        readEvent(r, "v_", v);
        readMain(r, m);
        const MainEvent before = m;
        const double c1 = model.sigcc(v, m);   // const: main must not change
        const bool unchanged = std::memcmp(&before.thetacm, &m.thetacm, sizeof(double)) == 0 &&
                               std::memcmp(&before.ntup.sigcm, &m.ntup.sigcm, sizeof(double)) == 0 &&
                               std::memcmp(&before.ntup.decdist, &m.ntup.decdist, sizeof(double)) == 0;
        const double c2 = model.sigccMain(v, m);
        const double rr = model.sigccRecon(v, m);
        Outputs o;
        putPeek(c2, km->survivalProbability(), m, km->config(), o);
        km->config() = saved;
        bool ok = same(c1, c2) && unchanged && rr == 1.0;
        for (const auto& c : kPeekCols) ok = ok && same(o[c], f[c]);
        ++rows;
        if (ok) {
            ++eq;
        } else if (bad++ < 5) {
            std::printf("  plugin row %zu (cfg %d): sigcc %.17g sigccMain %.17g ref %.17g unchanged %d recon %g\n", i + 1, cfg,
                        c1, c2, f["peek"], unchanged ? 1 : 0, rr);
        }
    }
    std::printf("plugin \"kaon\" (Registry): %zu rows, %zu bitwise (all outputs) -> %s\n", rows, eq, bad ? "FAIL" : "PASS");
    return bad ? 1 : 0;
}

// ------------------------------------------------------------------------------------------------- D-043
int runFortranDiff(const std::string& peekDir, const std::string& stateDir, const std::string& dataDir,
                   const std::string& eekeekDir, const std::string& saghaiDir) {
    std::map<int, CfgState> st;
    if (!loadStates(stateDir, st)) return 77;
    simcl1::Table in, ref, ein, eref;
    std::string err;
    if (!in.load(peekDir + "/input.txt", err) || !ref.load(peekDir + "/output.txt", err)) return 77;
    kaon::SaghaiTables tc, td;
    std::vector<std::string> msgs;
    if (!loadTables(dataDir, false, tc)) return 1;
    if (!td.load(dataDir, true, &msgs, &err)) {
        std::printf("%s\n", err.c_str());
        return 1;
    }
    int fail = 0;

    // ---- F2: the corrected reader puts the true ff1..ff3 into zs*ff1..3 (the compat reader has them in zs*ff4..6 of the
    // same record, check_file_mapping.txt) and reads ff4..ff6 (lost in compat).
    {
        std::size_t n = 0, eq = 0, diff46 = 0;
        for (int k = 0; k < 3; ++k)
            for (std::size_t i = 0; i < tc.zsrff[0].size(); ++i) {
                n += 2;
                eq += (td.zsrff[static_cast<std::size_t>(k)][i] == tc.zsrff[static_cast<std::size_t>(k + 3)][i] ? 1 : 0) +
                      (td.zsiff[static_cast<std::size_t>(k)][i] == tc.zsiff[static_cast<std::size_t>(k + 3)][i] ? 1 : 0);
                diff46 += (td.zsrff[static_cast<std::size_t>(k + 3)][i] != tc.zsrff[static_cast<std::size_t>(k + 3)][i]) +
                          (td.zsiff[static_cast<std::size_t>(k + 3)][i] != tc.zsiff[static_cast<std::size_t>(k + 3)][i]);
            }
        // the proton tables, weightc and weightd are read identically in both modes
        bool same_other = tc.weightc == td.weightc && tc.weightd == td.weightd;
        for (int k = 0; k < 6; ++k)
            same_other = same_other && tc.zrff[static_cast<std::size_t>(k)] == td.zrff[static_cast<std::size_t>(k)] &&
                         tc.ziff[static_cast<std::size_t>(k)] == td.ziff[static_cast<std::size_t>(k)];
        const bool ok = n == eq && diff46 > 0 && same_other && msgs.size() == 1 &&
                        msgs[0].find("18 records") != std::string::npos;
        std::printf("D-043 F2 (dbase.f:662-677 saghai_sigma0.dat one line off): corrected zs*ff1..3 == compat zs*ff4..6 "
                    "for %zu of %zu values; %zu ff4..6 values differ from compat; other tables identical: %s; header "
                    "message: %s -> %s\n",
                    eq, n, diff46, same_other ? "yes" : "NO", msgs.empty() ? "(none)" : msgs[0].c_str(), ok ? "PASS" : "FAIL");
        if (!ok) ++fail;
    }

    // ---- F1 + F3 on the eekeek rows: compat = 0 for every finite row; F1 only = the fintabi value; F3 moves only the
    // rows with s > 4.7 GeV^2 (between/beyond the corrected nodes 5.1 and 5.4)
    if (ein.load(eekeekDir + "/input.txt", err) && eref.load(eekeekDir + "/output.txt", err)) {
        std::size_t rows = 0, zeroCompat = 0, f3changed = 0, f3changedLowS = 0, f1nonzero = 0;
        kaon::SaghaiOptions oc, o1, o13;
        o1.fint_value = true;
        o13.fint_value = true;
        o13.proton_grid_fixed = true;
        for (std::size_t i = 0; i < ein.rows.size(); ++i) {
            Row r(ein, i);
            const double vc = kaon::eekeek(r["ss"], r["q22"], r["angl"], r["theta"], r["phi"], r["epsi"], r["mrec_struck"], tc, oc);
            const double v1 = kaon::eekeek(r["ss"], r["q22"], r["angl"], r["theta"], r["phi"], r["epsi"], r["mrec_struck"], tc, o1);
            const double v13 = kaon::eekeek(r["ss"], r["q22"], r["angl"], r["theta"], r["phi"], r["epsi"], r["mrec_struck"], tc, o13);
            ++rows;
            if (!std::isfinite(vc) || vc == 0.0) ++zeroCompat;
            if (v1 != 0.0 && std::isfinite(v1)) ++f1nonzero;
            if (!same(v1, v13)) {
                ++f3changed;
                if (r["ss"] <= static_cast<float>(4.7)) ++f3changedLowS;
            }
        }
        const bool ok = zeroCompat == rows && f1nonzero > 0 && f3changed > 0 && f3changedLowS == 0;
        std::printf("D-043 F1/F3 (eekeek): %zu rows: compat 0 (or NaN edge) in %zu; with the fint value non-zero in %zu; "
                    "the corrected s nodes 5.1/5.4 change %zu rows (%zu of them with s <= 4.7 GeV^2) -> %s\n",
                    rows, zeroCompat, f1nonzero, f3changed, f3changedLowS, ok ? "PASS" : "FAIL");
        if (!ok) ++fail;
    } else {
        std::printf("D-043 F1/F3: eekeek table missing (%s)\n", err.c_str());
        ++fail;
    }

    // ---- F5/F6 and the weight-relevant outputs on the peek rows: compat reproduces the table; corrected differs only
    // in survivalprob/decdist (decay OFF rows) and sigcm1
    {
        std::map<int, kaon::PeekState> pc, pd;
        std::size_t rows = 0, bad = 0, decOff = 0, spChanged = 0, shmsStale = 0, shmsStaleFixed = 0;
        double maxRelSp = 0;
        for (std::size_t i = 0; i < in.rows.size(); ++i) {
            Row r(in, i), f(ref, i);
            const int cfg = iv(r, "cfg");
            kaon::PeekConfig cc = peekConfig(st.at(cfg), false, false), cd = peekConfig(st.at(cfg), true, true);
            applyOverrides(r, cc);
            applyOverrides(r, cd);
            Event v;
            MainEvent mc, md;
            readEvent(r, "v_", v);
            readMain(r, mc);
            readMain(r, md);
            double spc = r["in_survivalprob"], spd = spc;
            const double zc0 = pc[cfg].zaero;
            const double sc = kaon::peeK(v, mc, spc, cc, tc, pc[cfg]);
            const double sd = kaon::peeK(v, md, spd, cd, td, pd[cfg]);
            ++rows;
            bool ok = same(sc, f["peek"]) && same(spc, f["survivalprob"]) && same(mc.ntup.decdist, f["decdist_out"]) &&
                      same(sc, sd) && same(mc.ntup.sigcm2, md.ntup.sigcm2) && same(mc.wcm, md.wcm);
            if (!cc.doing_decay) {
                ++decOff;
                // corrected: beta*gamma of p_kaon = spec%p%P*(1+delta/100); zaero = 0 for the SHMS; sqrt
                const double pk = cd.spec_p_P * (1. + md.SP.p.delta / 100.);
                const double bk = pk / std::sqrt(pk * pk + cd.Mh2);
                const double gk = 1. / std::sqrt(1. - bk * bk);
                double za = 0.;
                if (cd.hadron_arm == 2) za = -82.8;
                if (cd.hadron_arm == 3 || cd.hadron_arm == 4) za = -183.;
                const double path = md.FP.p.path + za * std::sqrt(1 + md.FP.p.dx * md.FP.p.dx + md.FP.p.dy * md.FP.p.dy);
                const double expect = 1. / std::exp(path / (cd.ctau * bk * gk));
                ok = ok && std::fabs(spd - expect) <= 1e-15 * std::fabs(expect) && same(md.ntup.decdist, spd);
                if (!same(spc, spd)) {
                    ++spChanged;
                    maxRelSp = std::max(maxRelSp, std::fabs(spd / spc - 1.));
                }
                if ((cc.hadron_arm == 5 || cc.hadron_arm == 6) && zc0 != 0.) {
                    ++shmsStale;   // compat: the stale zaero of an earlier arm 2/3/4 row is used
                    if (pd[cfg].zaero == 0.) ++shmsStaleFixed;
                }
            } else {
                ok = ok && same(spd, spc) && same(md.ntup.decdist, mc.ntup.decdist);
            }
            if (!ok && bad++ < 5)
                std::printf("  peek row %zu (cfg %d): compat %.17g ref %.17g corrected %.17g sp %.17g/%.17g\n", i + 1, cfg, sc,
                            f["peek"], sd, spc, spd);
        }
        const bool ok = bad == 0 && spChanged > 0 && shmsStale > 0 && shmsStale == shmsStaleFixed;
        std::printf("D-043 F5/F6 (physics_kaon.f:148-165): %zu rows, %zu with doing_decay = 0: corrected survival "
                    "probability differs in %zu (max |rel| %.3e), %zu SHMS rows with a stale zaero in compat (0 in "
                    "corrected: %zu); peeK, sigcm2, wcm identical in both modes; %zu bad -> %s\n",
                    rows, decOff, spChanged, maxRelSp, shmsStale, shmsStaleFixed, bad, ok ? "PASS" : "FAIL");
        if (!ok) ++fail;
    }
    std::printf("fortran_diff D-043 (kaon model): %s\n", fail ? "FAIL" : "PASS");
    (void)saghaiDir;
    return fail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mode> ... (see the file header)\n", argv[0]);
        return 2;
    }
    const std::string what = argv[1];
    try {
        if (what == "files" && argc > 3) return runFiles(argv[2], argv[3]);
        if (what == "saghai_tables" && argc > 3) return runSaghaiTables(argv[2], argv[3]);
        if (what == "sig_factorized" && argc > 3) return runSigFactorized(argv[2], std::atof(argv[3]));
        if ((what == "eekeek" || what == "eekeeks") && argc > 5)
            return runSaghai(what, argv[2], argv[3], std::atof(argv[4]), std::atoi(argv[5]) != 0);
        if (what == "phspwght" && argc > 5) return runPhspwght(argv[2], argv[3], std::atof(argv[4]), std::atoi(argv[5]) != 0);
        if (what == "peek" && argc > 6) return runPeek(argv[2], argv[3], argv[4], std::atof(argv[5]), std::atoi(argv[6]) != 0);
        if (what == "plugin" && argc > 4) return runPlugin(argv[2], argv[3], argv[4]);
        if (what == "fortran_diff_D043" && argc > 6) return runFortranDiff(argv[2], argv[3], argv[4], argv[5], argv[6]);
    } catch (const std::exception& e) {
        std::printf("ERROR: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "bad arguments for mode %s\n", what.c_str());
    return 2;
}
