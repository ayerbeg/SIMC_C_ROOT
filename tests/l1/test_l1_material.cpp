// tests/l1/test_l1_material.cpp
//
// L1 tests of the material / target routines (task P3.1) against the Fortran harness tables in
// reference/functions/<table>/ (conventions: reference/functions/README.md, reference/P3.F1_report.md).
//
// Usage: test_l1_material <mode> <table_dir> [rel_tol | cfgid]
//   gauss1               gauss1_fixedrng (real algorithm on scripted uniforms)
//   enerloss_new         enerloss_new (typeflags 4,2,3 + static-x rows)
//   enerloss_fixedrng    enerloss_new_fixedrng (typeflag 1, scripted gauss1)
//   musc | musc_ext | target_musc   *_fixedrng tables (scripted gauss1)
//   trip                 trip_thru_target_fixedrng
//   target_init          target_init
//   extreme              extreme_trip_thru_target
//   distributions <dir> <cfgid>     real-RNG distributions: Mt64Random, 1e6 samples, quantiles + KS + chi2
//   selftest             plug-in registry, exceptions, ScriptedRandom counting, quirks
// Deterministic and fixed-RNG tables: every output column |rel| <= rel_tol (default 1e-12); the RNG
// consumption columns (n_gauss, n_unif, nsig_arg, rng_err) must be identical. The WARNING lines printed by
// the routines must match the Fortran driver.log text.
// Exit: 0 pass, 1 fail, 77 table missing.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "TMath.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/Gauss1.h"
#include "simc/core/random/Mt64Random.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/material/HallCCryoTarget.h"
#include "simc/material/SimcMaterial.h"

using namespace simc;
using simc::material::HallCCryoTarget;
using simc::material::Limits;
using simc::material::SimcMaterial;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

double g_tol = 1e-12;

std::vector<ColumnSpec> cols(std::initializer_list<const char*> real, std::initializer_list<const char*> exact) {
    std::vector<ColumnSpec> v;
    for (auto n : real) v.push_back({n, g_tol, 0.0});
    for (auto n : exact) v.push_back({n, 0.0, 0.0});
    return v;
}

// ---- WARNING-line comparison with driver.log ---------------------------------------------------------
std::vector<std::string> lines(std::istream& in, bool dropDriver) {
    std::vector<std::string> v;
    std::string l;
    while (std::getline(in, l)) {
        if (dropDriver && (l.rfind(" drv_", 0) == 0 || l.find("RANLUX") != std::string::npos)) continue;
        v.push_back(l);
    }
    return v;
}

int compareLog(const std::string& dir, const std::string& cppLog) {
    std::ifstream f(dir + "/driver.log");
    if (!f) {
        std::printf("  driver.log missing: warning comparison skipped\n");
        return 0;
    }
    std::istringstream c(cppLog);
    auto fw = lines(f, true), cw = lines(c, false);
    std::size_t same = 0;
    for (std::size_t i = 0; i < fw.size() && i < cw.size(); ++i) same += (fw[i] == cw[i]);
    std::printf("  warning lines: Fortran %zu, C++ %zu, identical text %zu\n", fw.size(), cw.size(), same);
    if (fw.size() != cw.size() || same != fw.size()) {
        for (std::size_t i = 0; i < fw.size() || i < cw.size(); ++i) {
            std::printf("   F: [%s]\n   C: [%s]\n", i < fw.size() ? fw[i].c_str() : "-", i < cw.size() ? cw[i].c_str() : "-");
            if (i > 8) break;
        }
        std::printf("  => FAIL (warnings)\n");
        return 1;
    }
    return 0;
}

int withLog(int rc, const std::string& dir, const std::ostringstream& log) {
    if (rc != 0) return rc;
    return compareLog(dir, log.str());
}

void counters(const ScriptedRandom& rng, Outputs& o, bool unif = true) {
    o["n_gauss"] = static_cast<double>(rng.nGauss());
    if (unif) o["n_unif"] = static_cast<double>(rng.nUniform());
    o["nsig_arg"] = rng.lastNsig();
    o["rng_err"] = rng.error() ? 1.0 : 0.0;
}

// ---- target state as fortran_harness/target_eloss/tstate.f:set_target_state ---------------------------
void setTargetState(HallCCryoTarget& T, const Row& r, double* EbeamOut = nullptr, double* MhOut = nullptr) {
    auto& targ = T.targ();
    targ.A = r["targ_A"];
    targ.Z = r["targ_Z"];
    targ.N = targ.A - targ.Z;
    targ.rho = r["targ_rho"];
    targ.thick = r["targ_thick_mg"] / 1000.;
    targ.length = targ.thick / targ.rho;
    targ.angle = r["targ_angle_deg"] / constants::degrad;
    targ.can = static_cast<int>(std::lround(r["targ_can"]));
    T.setArms(static_cast<int>(std::lround(r["electron_arm"])), static_cast<int>(std::lround(r["hadron_arm"])));
    const double Ebeam = r["Ebeam"];
    const double the = std::fabs(r["the_cent_deg"]) / constants::degrad;
    const double thp = std::fabs(r["thp_cent_deg"]) / constants::degrad;
    const double Mh = r["Mh"];
    const double Mh2 = Mh * Mh;
    const bool ue = std::lround(r["using_Eloss"]) > 0;
    const bool uc = std::lround(r["using_Coulomb"]) > 0;
    T.targetInit(ue, uc, the, thp, r["Pp_cent"], Mh2, Ebeam, r["Pe_cent"]);
    if (EbeamOut) *EbeamOut = Ebeam;
    if (MhOut) *MhOut = Mh;
}

// ---- deterministic / fixed-RNG tables ------------------------------------------------------------------
int runGauss1(const std::string& dir) {
    ScriptedRandom rng(ScriptedRandom::GaussMode::Algorithm);
    return simcl1::run("gauss1_fixedrng (" + dir + ")", dir, cols({"gauss1"}, {"n_unif", "rng_err"}),
                       [&](const Row& r, Outputs& o) {
                           std::vector<double> u;
                           for (int i = 1; i <= 12; ++i) u.push_back(r["u" + std::to_string(i)]);
                           rng.load({}, u);
                           o["gauss1"] = rng.gauss1(r["nsigmax"]);
                           o["n_unif"] = static_cast<double>(rng.nUniform());
                           o["rng_err"] = rng.error() ? 1.0 : 0.0;
                       });
}

int runEnerloss(const std::string& dir, bool fixedrng) {
    std::ostringstream log;
    SimcMaterial mat;   // one instance = one Fortran process (static x, numerr)
    mat.setLog(&log);
    ScriptedRandom rng;
    int rc = simcl1::run(std::string(fixedrng ? "enerloss_new_fixedrng" : "enerloss_new") + " (" + dir + ")", dir,
                         cols({"Eloss"}, {"n_gauss", "n_unif", "nsig_arg", "rng_err"}), [&](const Row& r, Outputs& o) {
                             rng.load({r["g1"]});
                             o["Eloss"] = mat.energyLoss(r["len"], r["dens"], r["zeff"], r["aeff"], r["epart"], r["mpart"],
                                                         static_cast<int>(std::lround(r["typeflag"])), rng);
                             counters(rng, o);
                         });
    return withLog(rc, dir, log);
}

int runMusc(const std::string& dir) {
    std::ostringstream log;
    SimcMaterial mat;
    mat.setLog(&log);
    ScriptedRandom rng;
    int rc = simcl1::run("musc_fixedrng (" + dir + ")", dir, cols({"dth", "dph"}, {"n_gauss", "nsig_arg", "rng_err"}),
                         [&](const Row& r, Outputs& o) {
                             rng.load({r["g1"], r["g2"]});
                             double dth = r["dth_in"], dph = r["dph_in"];
                             mat.musc(r["m2"], r["p"], r["rad_len"], dth, dph, rng);
                             o["dth"] = dth;
                             o["dph"] = dph;
                             counters(rng, o, false);
                         });
    return withLog(rc, dir, log);
}

int runMuscExt(const std::string& dir) {
    std::ostringstream log;
    SimcMaterial mat;
    mat.setLog(&log);
    ScriptedRandom rng;
    int rc = simcl1::run("musc_ext_fixedrng (" + dir + ")", dir,
                         cols({"dph", "dth", "y", "x"}, {"n_gauss", "nsig_arg", "rng_err"}), [&](const Row& r, Outputs& o) {
                             rng.load({r["g1"], r["g2"], r["g3"], r["g4"]});
                             double dph = r["dph_in"], dth = r["dth_in"], y = r["y_in"], x = r["x_in"];
                             mat.muscExt(r["m2"], r["p"], r["rad_len"], r["x_len"], dph, dth, y, x, rng);
                             o["dph"] = dph;
                             o["dth"] = dth;
                             o["y"] = y;
                             o["x"] = x;
                             counters(rng, o, false);
                         });
    return withLog(rc, dir, log);
}

int runTargetMusc(const std::string& dir) {
    std::ostringstream log;
    HallCCryoTarget T;
    T.setLog(&log);
    ScriptedRandom rng;
    int rc = simcl1::run("target_musc_fixedrng (" + dir + ")", dir,
                         cols({"dangle_1", "dangle_2"}, {"n_gauss", "nsig_arg", "rng_err"}), [&](const Row& r, Outputs& o) {
                             rng.load({r["g1"], r["g2"]});
                             double d[2] = {0, 0};
                             T.targetMusc(r["p"], r["beta"], r["teff"], d, rng);
                             o["dangle_1"] = d[0];
                             o["dangle_2"] = d[1];
                             counters(rng, o, false);
                         });
    return withLog(rc, dir, log);
}

int runTrip(const std::string& dir) {
    std::ostringstream log;
    HallCCryoTarget T;   // one instance for all rows = one Fortran process (static th_can, enerloss state)
    T.setLog(&log);
    ScriptedRandom rng;
    int rc = simcl1::run("trip_thru_target_fixedrng (" + dir + ")", dir,
                         cols({"Eloss", "radlen", "X0_cm", "targ_length"}, {"n_gauss", "nsig_arg", "rng_err"}),
                         [&](const Row& r, Outputs& o) {
                             setTargetState(T, r);
                             rng.load({r["g1"], r["g2"], r["g3"], r["g4"], r["g5"]});
                             double el = 0, rl = 0;
                             T.trip(static_cast<int>(std::lround(r["narm"])), r["zpos"], r["energy"], r["theta"], el, rl,
                                    r["mass"], static_cast<int>(std::lround(r["typeflag"])), rng);
                             o["Eloss"] = el;
                             o["radlen"] = rl;
                             counters(rng, o, false);
                             o["X0_cm"] = T.targ().X0_cm;
                             o["targ_length"] = T.targ().length;
                         });
    return withLog(rc, dir, log);
}

int runTargetInit(const std::string& dir) {
    std::ostringstream log;
    HallCCryoTarget T;
    T.setLog(&log);
    int rc = simcl1::run("target_init (" + dir + ")", dir,
                         cols({"targ_thick", "targ_length", "L1", "L2", "X0", "X0_cm", "Eloss_ave_1", "Eloss_ave_2",
                               "Eloss_ave_3", "teff_ave_1", "teff_ave_2", "teff_ave_3", "Coulomb_ave", "Coulomb_const"},
                              {}),
                         [&](const Row& r, Outputs& o) {
                             setTargetState(T, r);
                             const auto& t = T.targ();
                             o["targ_thick"] = t.thick;
                             o["targ_length"] = t.length;
                             o["L1"] = t.L1;
                             o["L2"] = t.L2;
                             o["X0"] = t.X0;
                             o["X0_cm"] = t.X0_cm;
                             for (int i = 0; i < 3; ++i) {
                                 o["Eloss_ave_" + std::to_string(i + 1)] = t.Eloss[i].ave;
                                 o["teff_ave_" + std::to_string(i + 1)] = t.teff[i].ave;
                             }
                             o["Coulomb_ave"] = t.Coulomb.ave;
                             o["Coulomb_const"] = t.Coulomb_constant;
                         });
    return withLog(rc, dir, log);
}

int runExtreme(const std::string& dir) {
    std::ostringstream log;
    HallCCryoTarget T;
    T.setLog(&log);
    std::vector<ColumnSpec> c = cols({"targ_length", "X0_cm"}, {"musc_nsig_max"});
    for (int i = 1; i <= 3; ++i)
        for (const char* n : {"Eloss_min_", "Eloss_max_", "teff_min_", "teff_max_", "musc_max_"})
            c.push_back({n + std::to_string(i), g_tol, 0.0});
    int rc = simcl1::run("extreme_trip_thru_target (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        double Ebeam = 0, Mh = 0;
        setTargetState(T, r, &Ebeam, &Mh);
        const auto& t = T.targ();
        Limits the{r["the_min"], r["the_max"]}, thp{r["thp_min"], r["thp_max"]};
        Limits pe{r["pe_min"], r["pe_max"]}, pp{r["pp_min"], r["pp_max"]};
        Limits z{-0.5 * t.length, 0.5 * t.length};   // init.f:243-244
        T.extremeTripThruTarget(Ebeam, the, thp, pe, pp, z, Mh);
        o["targ_length"] = t.length;
        o["X0_cm"] = t.X0_cm;
        for (int i = 0; i < 3; ++i) {
            const std::string k = std::to_string(i + 1);
            o["Eloss_min_" + k] = t.Eloss[i].min;
            o["Eloss_max_" + k] = t.Eloss[i].max;
            o["teff_min_" + k] = t.teff[i].min;
            o["teff_max_" + k] = t.teff[i].max;
            o["musc_max_" + k] = t.musc_max[i];
        }
        o["musc_nsig_max"] = t.musc_nsig_max;
    });
    return withLog(rc, dir, log);
}

// ---- distributions (real RNG) ---------------------------------------------------------------------------
struct Hist {
    std::vector<double> lo, hi;
    std::vector<double> cnt;
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
    // under/overflow counts are in the comment line
    std::ifstream f(path);
    std::string l;
    while (std::getline(f, l)) {
        auto p = l.find("underflow");
        if (l.size() && l[0] == '#' && p != std::string::npos) {
            std::sscanf(l.c_str() + p, "underflow %lf overflow %lf", &h.under, &h.over);
            break;
        }
    }
    return true;
}

// numpy 'inverted_cdf': smallest sample x with ECDF(x) >= p
double quantileInv(const std::vector<double>& s, double p) {
    std::size_t n = s.size();
    std::size_t k = static_cast<std::size_t>(std::ceil(p * static_cast<double>(n)));
    if (k < 1) k = 1;
    return s[k - 1];
}

double countBelow(const std::vector<double>& s, double x) {   // #(s < x)
    return static_cast<double>(std::lower_bound(s.begin(), s.end(), x) - s.begin());
}
double countAtOrBelow(const std::vector<double>& s, double x) {   // #(s <= x)
    return static_cast<double>(std::upper_bound(s.begin(), s.end(), x) - s.begin());
}

/// Returns 0 pass, 1 fail. Prints one block per variable.
int compareDistribution(const std::string& dir, int cfgid, const std::string& var, std::vector<double>& s,
                        const std::map<std::string, double>& fsum) {
    std::sort(s.begin(), s.end());
    const double n = static_cast<double>(s.size());
    char hp[64];
    std::snprintf(hp, sizeof hp, "/hist_%02d_%s.txt", cfgid, var.c_str());
    Hist h;
    if (!loadHist(dir + hp, h)) {
        std::printf("  missing %s\n", hp);
        return 1;
    }
    const double nf = fsum.at("nsamp");
    // Quantiles: z = (F_cpp(q_fortran) - p) / sqrt(p(1-p)(1/n_cpp + 1/n_f))
    const double Q[] = {0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99};
    const char* QN[] = {"q01", "q05", "q25", "q50", "q75", "q95", "q99"};
    double zmax = 0;
    int fail = 0;
    std::printf("  cfg %d %s: n_cpp=%.0f n_f=%.0f\n", cfgid, var.c_str(), n, nf);
    std::printf("    %-4s %22s %22s %10s %8s\n", "q", "Fortran", "C++", "rel diff", "z");
    for (int i = 0; i < 7; ++i) {
        const double qf = fsum.at(QN[i]);
        const double qc = quantileInv(s, Q[i]);
        const double Fc = countAtOrBelow(s, qf) / n;   // ECDF of the C++ sample at the Fortran quantile
        const double sig = std::sqrt(Q[i] * (1 - Q[i]) * (1. / n + 1. / nf));
        const double z = (Fc - Q[i]) / sig;
        zmax = std::max(zmax, std::fabs(z));
        if (std::fabs(z) > 4.) ++fail;
        std::printf("    %-4s %22.14e %22.14e %10.2e %8.3f\n", QN[i], qf, qc, (qc - qf) / std::fabs(qf), z);
    }
    // KS on the 201 bin edges (numpy histogram: [lo_i, hi_i) except the last bin, which includes hi)
    const std::size_t nb = h.cnt.size();
    double cumF = h.under, D = 0;
    double chi2 = 0;
    int ndf = -1;
    auto addChi = [&](double a, double b) {   // a: Fortran count, b: C++ count scaled to the Fortran N
        if (a + b > 0) {
            chi2 += (a - b) * (a - b) / (a + b);
            ++ndf;
        }
    };
    const double scale = nf / n;
    double prevC = 0;
    for (std::size_t i = 0; i <= nb; ++i) {
        double cC;   // C++ count below edge i
        if (i < nb)
            cC = countBelow(s, h.lo[i]);
        else
            cC = countAtOrBelow(s, h.hi[nb - 1]);
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
    const double pKS = TMath::KolmogorovProb(std::sqrt(ne) * D);
    const double pChi2 = TMath::Prob(chi2, ndf);
    if (pKS < 1e-3) ++fail;
    if (pChi2 < 1e-3) ++fail;
    // informative (not tested, D-012): mean / rms / IQR
    double m = 0, m2 = 0;
    for (double x : s) m += x;
    m /= n;
    for (double x : s) m2 += (x - m) * (x - m);
    const double rms = std::sqrt(m2 / n);
    const double iqrC = quantileInv(s, 0.75) - quantileInv(s, 0.25), iqrF = fsum.at("q75") - fsum.at("q25");
    std::printf("    KS (binned, 201 edges): D=%.3e p=%.4f | chi2(%d bins)=%.1f/%d p=%.4f | max|z_q|=%.2f\n", D, pKS,
                ndf + 1, chi2, ndf, pChi2, zmax);
    std::printf("    info: mean F %.6e C %.6e | rms F %.6e C %.6e | IQR F %.6e C %.6e | max F %.6e C %.6e\n",
                fsum.at("mean"), m, fsum.at("rms"), rms, iqrF, iqrC, fsum.at("max"), s.back());
    std::printf("    => %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

int runDistributions(const std::string& dir, int cfgid) {
    simcl1::Table cfg;
    std::string err;
    if (!cfg.load(dir + "/config.txt", err) || cfg.rows.empty()) {
        std::printf("distributions: SKIP (%s)\n", err.c_str());
        return 77;
    }
    // summary.txt: parse manually (the 'var' column is text)
    std::map<std::string, std::map<std::string, double>> fsum;   // var -> column -> value
    {
        std::ifstream f(dir + "/summary.txt");
        std::string l;
        std::vector<std::string> names;
        while (std::getline(f, l)) {
            std::istringstream ss(l);
            std::vector<std::string> w;
            std::string t;
            while (ss >> t) w.push_back(t);
            if (w.empty()) continue;
            if (names.empty()) {
                names = w;
                continue;
            }
            if (std::lround(std::atof(w[0].c_str())) != cfgid) continue;
            std::map<std::string, double> m;
            std::string var;
            for (std::size_t i = 0; i < names.size() && i < w.size(); ++i) {
                if (names[i] == "var")
                    var = w[i];
                else
                    m[names[i]] = std::atof(w[i].c_str());
            }
            fsum[var] = m;
        }
    }
    const simcl1::Table& t = cfg;
    std::size_t irow = t.rows.size();
    for (std::size_t i = 0; i < t.rows.size(); ++i)
        if (std::lround(Row(t, i)["cfgid"]) == cfgid) irow = i;
    if (irow == t.rows.size() || fsum.empty()) {
        std::printf("distributions: cfgid %d not found\n", cfgid);
        return 1;
    }
    Row r(t, irow);
    const int kind = static_cast<int>(std::lround(r["kind"]));
    const std::size_t nsamp = static_cast<std::size_t>(std::lround(r["nsamp"]));
    // The C++ stream is independent of RANLUX (D-002); the seed is fixed (same number as the Fortran
    // configuration) so the test is reproducible.
    Mt64Random rng(static_cast<std::uint64_t>(std::lround(r["seed"])));
    std::printf("L1 distributions cfgid %d kind %d: Mt64Random seed %ld, %zu samples\n", cfgid, kind, std::lround(r["seed"]),
                nsamp);
    std::ostringstream log;
    int fail = 0;
    if (kind == 1) {
        SimcMaterial mat;
        mat.setLog(&log);
        std::vector<double> s(nsamp);
        for (auto& x : s) x = mat.energyLoss(r["len"], r["dens"], r["zeff"], r["aeff"], r["epart"], r["mpart"], 1, rng);
        fail += compareDistribution(dir, cfgid, "Eloss", s, fsum.at("Eloss"));
    } else if (kind == 2) {
        SimcMaterial mat;
        mat.setLog(&log);
        std::vector<double> a(nsamp), b(nsamp);
        for (std::size_t i = 0; i < nsamp; ++i) {
            double dth = 0., dph = 0.;
            mat.musc(r["m2"], r["p"], r["rad_len"], dth, dph, rng);
            a[i] = dth;
            b[i] = dph;
        }
        fail += compareDistribution(dir, cfgid, "dth", a, fsum.at("dth"));
        fail += compareDistribution(dir, cfgid, "dph", b, fsum.at("dph"));
    } else {
        HallCCryoTarget T;
        T.setLog(&log);
        setTargetState(T, r);
        std::vector<double> s(nsamp);
        double rl = 0;
        const int narm = static_cast<int>(std::lround(r["narm"]));
        for (auto& x : s) T.trip(narm, r["zpos"], r["energy"], r["theta"], x, rl, r["mass"], 1, rng);
        std::printf("  radlen = %.16e\n", rl);
        fail += compareDistribution(dir, cfgid, "Eloss", s, fsum.at("Eloss"));
    }
    return fail ? 1 : 0;
}

// ---- self test ----------------------------------------------------------------------------------------
int runSelftest() {
    int fail = 0;
    auto expect = [&](bool ok, const std::string& what) {
        std::printf("  %-70s %s\n", what.c_str(), ok ? "ok" : "FAIL");
        if (!ok) ++fail;
    };
    SimcParameters params;
    expect(Registry<IRandom>::contains("mt64"), "IRandom 'mt64' registered");
    expect(Registry<IMaterialModel>::contains("simc_material"), "IMaterialModel 'simc_material' registered");
    expect(Registry<ITarget>::contains("hallc_cryo"), "ITarget 'hallc_cryo' registered");

    // Mt64Random: range, reproducibility, truncation
    {
        Mt64Random a(42), b(42);
        bool inRange = true, same = true;
        for (int i = 0; i < 100000; ++i) {
            double u = a.uniform();
            inRange = inRange && u > 0. && u < 1.;
            same = same && (u == b.uniform());
        }
        expect(inRange, "Mt64Random uniform() in (0,1)");
        expect(same, "Mt64Random reproducible for equal seeds");
        bool trunc = true;
        for (int i = 0; i < 100000; ++i) trunc = trunc && std::fabs(a.gauss1(0.5)) <= 0.5;
        expect(trunc, "Mt64Random gauss1(0.5) truncated at 0.5");
        auto r = Registry<IRandom>::create("mt64");
        r->configure(Json{{"seed", 42}}, params);
        Mt64Random c(42);
        expect(r->uniform() == c.uniform(), "mt64 plug-in honours JSON seed");
    }
    // ScriptedRandom: counting and exhaustion
    {
        ScriptedRandom s;
        s.load({0.25}, {0.1});
        double g = s.gauss1(3.5), g2 = s.gauss1(99.);
        double u = s.uniform(), u2 = s.uniform();
        expect(g == 0.25 && std::isnan(g2) && u == 0.1 && std::isnan(u2) && s.error(), "ScriptedRandom NaN on exhaustion");
        expect(s.nGauss() == 2 && s.nUniform() == 2 && s.lastNsig() == 99., "ScriptedRandom counts and nsig");
        // gauss1 s == 1 edge (F6): u1 = 1, u2 = 0.5 -> v1 = 1, v2 = 0, s = 1 -> g = +0
        ScriptedRandom a(ScriptedRandom::GaussMode::Algorithm);
        a.load({}, {1.0, 0.5});
        double e = a.gauss1(10.);
        expect(e == 0.0 && a.nUniform() == 2, "gauss1 accepts s == 1 (returns 0, F6)");
    }
    // enerloss_new: invalid typeflag keeps the previous x (F3)
    {
        SimcMaterial m;
        m.setLog(nullptr);
        ScriptedRandom none;
        double e3 = m.energyLoss(1.0, 0.0723, 1., 1., 1000., constants::Me, 3, none);
        double e0 = m.energyLoss(1.0, 0.0723, 1., 1., 1000., constants::Me, 0, none);
        expect(e3 == e0, "enerloss_new typeflag 0 reuses x of the previous call (F3)");
        SimcMaterial fresh;
        fresh.setLog(nullptr);
        double big = fresh.energyLoss(1.0, 0.0723, 1., 1., 1e6, constants::Me, 7, none);   // x = 0 -> lambda = 1e5
        expect(big > 100., "enerloss_new x=0 -> lambda=1e5 (D-012)");
    }
    // musc_ext STOP conditions -> FortranStop; rad_len == 0 returns first
    {
        SimcMaterial m;
        m.setLog(nullptr);
        ScriptedRandom s;
        s.load({1, 1, 1, 1});
        double a = 0, b = 0, c = 0, d = 0;
        bool t1 = false, t2 = false;
        try {
            m.muscExt(0.26, 1000., 0.01, 0.0, a, b, c, d, s);
        } catch (const simc::FortranStop&) {
            t1 = true;
        }
        try {
            m.muscExt(0.26, 1000., -0.01, 1.0, a, b, c, d, s);
        } catch (const simc::FortranStop&) {
            t2 = true;
        }
        expect(t1 && t2, "musc_ext STOP (x_len<=0, rad_len<0) -> FortranStop");
        m.muscExt(0.26, 1000., 0.0, -1.0, a, b, c, d, s);
        expect(s.nGauss() == 0, "musc_ext rad_len==0 returns before the STOP test");
    }
    // Target: configure from parameters (dbase.f conversions) + init
    {
        SimcParameters p;
        p.targ.A = 1;
        p.targ.Z = 1;
        p.targ.rho = 0.0723;
        p.targ.thick = 100.;
        p.targ.can = 1;
        p.electron_arm = 5;
        p.hadron_arm = 1;
        p.Ebeam = 10600.;
        p.spec.e.theta = 12.5;
        p.spec.p.theta = 36.594;
        p.spec.e.P = 8360.989;
        p.spec.p.P = 3035.585;
        p.using_Eloss = true;
        p.using_Coulomb = true;   // forced off for Z=1
        auto T = Registry<ITarget>::create("hallc_cryo");
        T->configure(Json::object(), p);
        RunSetup setup;
        T->init(setup);
        auto* H = dynamic_cast<HallCCryoTarget*>(T.get());
        expect(H != nullptr, "hallc_cryo is a HallCCryoTarget");
        // target_init table cfg 1: Eloss_ave_1 = 2.3356216853792564E-01, teff_ave_2 = 1.2757750340418495E-02
        expect(H && H->targ().Eloss[0].ave == 2.3356216853792564E-01, "configure+init reproduce target_init cfg 1 Eloss_ave_1");
        expect(H && H->targ().teff[1].ave == 1.2757750340418495E-02, "configure+init reproduce target_init cfg 1 teff_ave_2");
        expect(H && H->targ().Coulomb.ave == 0.0, "using_Coulomb forced off for Z=1 (dbase.f:547)");
        Mt64Random r(1);
        bool zok = true;
        for (int i = 0; i < 10000; ++i) zok = zok && std::fabs(T->generateZ(r)) <= H->targ().length / 2.;
        expect(zok, "generateZ within +-length/2");
        // calorimeter arm -> throw
        H->setArms(7, 1);
        bool thr = false;
        double el, rl;
        ScriptedRandom none;
        try {
            H->trip(2, 0., 1000., 0.2, el, rl, constants::Me, 4, none);
        } catch (const std::invalid_argument&) {
            thr = true;
        }
        expect(thr, "trip_thru_target electron_arm 7 -> invalid_argument (out of scope)");
        // narm outside 1..3 behaves like the beam (computed GOTO fall-through)
        double e1, r1, e0, r0;
        H->trip(1, 0.1, 10600., 0., e1, r1, constants::Me, 4, none);
        H->trip(0, 0.1, 10600., 0.7, e0, r0, constants::Me, 4, none);
        expect(e1 == e0 && r1 == r0, "trip_thru_target narm=0 falls through to the beam branch");
        // bad can for a liquid -> FortranStop (dbase.f:461)
        p.targ.can = 4;
        bool thc = false;
        try {
            T->configure(Json::object(), p);
        } catch (const simc::FortranStop&) {
            thc = true;
        }
        expect(thc, "liquid target with can 4 -> FortranStop (dbase.f:461)");
    }
    // Pudding can: stale th_can (F4)
    {
        HallCCryoTarget T;
        std::ostringstream log;
        T.setLog(&log);
        auto& t = T.targ();
        t.A = 1;
        t.Z = 1;
        t.rho = 0.0723;
        t.thick = 0.2892;
        t.length = t.thick / t.rho;
        t.can = 2;
        T.setArms(1, 1);
        T.targetInit(true, false, 0.3, 0.6, 1000., constants::Mp * constants::Mp, 5000., 3000.);
        ScriptedRandom none;
        double el, rl;
        T.trip(2, 0.5, 3000., 0.4, el, rl, constants::Me, 4, none);
        const double th = T.thCan();
        T.trip(2, 3.9, 3000., 1.0, el, rl, constants::Me, 4, none);   // outside the can: z_can/R not <= 1
        expect(th != 0.0 && T.thCan() == th, "pudding can: th_can keeps its previous value (F4)");
        expect(log.str().find("z_can > can radius") != std::string::npos, "pudding can: warning printed");
    }
    std::printf("selftest: %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mode> <table_dir> [rel_tol|cfgid]\n", argv[0]);
        return 2;
    }
    const std::string mode = argv[1];
    if (mode == "selftest") return runSelftest();
    if (argc < 3) return 2;
    const std::string dir = argv[2];
    if (mode == "distributions") {
        if (argc < 4) return 2;
        return runDistributions(dir, std::atoi(argv[3]));
    }
    if (argc > 3) g_tol = std::atof(argv[3]);
    if (mode == "gauss1") return runGauss1(dir);
    if (mode == "enerloss_new") return runEnerloss(dir, false);
    if (mode == "enerloss_fixedrng") return runEnerloss(dir, true);
    if (mode == "musc") return runMusc(dir);
    if (mode == "musc_ext") return runMuscExt(dir);
    if (mode == "target_musc") return runTargetMusc(dir);
    if (mode == "trip") return runTrip(dir);
    if (mode == "target_init") return runTargetInit(dir);
    if (mode == "extreme") return runExtreme(dir);
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
