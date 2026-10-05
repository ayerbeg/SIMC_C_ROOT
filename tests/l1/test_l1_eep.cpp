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
#include "simc/core/random/ScriptedRandom.h"
#include "simc/input/InpReader.h"
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
    } catch (const std::exception& e) {
        std::printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "usage: %s files|theory_init|sf_lookup_init|sf_lookup|generate_em ...\n", argv[0]);
    return 1;
}
