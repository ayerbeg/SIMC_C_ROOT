// tests/l1/test_l1_hms.cpp
//
// L1 tests of the ISpectrometer "hms" (task P6.2) against the Fortran harness tables (task P6.F1):
//   mc   <table_dir> [rel_tol]  : reference/functions/mc_hms (flags off) and mc_hms_smear (ms=wcs=1, scripted
//                                 g1..g210 / u1): every output column of drv_mc.f, incl. the stop code, the
//                                 random-number consumption and the aperture counters summed over the table.
//   p12  <table_dir>            : reference/functions/mc_p12_check: the HMS (hadron-arm) predictions hs* of chk_p12.f
//                                 for the solved target points of the 5 eep_h_elastic_norad_A events (both modes).
//   dist <ref_dir> <config> [nray] [seed] : reference/functions/mc_distributions, config hms_nosmear | hms_smear:
//                                 1e6 rays from the drv_mcdist.f flat source with Mt64Random (D-002: the sequence
//                                 differs, the distributions must agree): stop fractions (two-sample binomial pulls),
//                                 binned chi2 + binned KS against the 200-bin reference histograms, quantiles.
//   selftest                    : lfit, hit_dipole, pion_coll_absorb, registry, collimator/decay smoke tests.
// The COSY files come from HmsSpectrometer::defaultMatrixDir() (data/matrices/simc_fortran, md5-checked by
// l1_optics_cosy_files). Exit: 0 pass, 1 fail, 77 table missing.
#include <algorithm>
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
#include "simc/core/Registry.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/Mt64Random.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/spectrometers/hms/HmsSpectrometer.h"
#include "simc/spectrometers/hms/Lfit.h"

using namespace simc;
using simc::hms::HmsSpectrometer;
using simc::hms::McHmsArgs;
using simc::hms::McHmsDetail;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

constexpr double kPi = 3.141592653589793;   // constants.inc
constexpr double kDegrad = 180. / kPi;
constexpr double kMp = 938.27231;
constexpr double kUnset = -999999.;          // drv_mc.f / drv_mcdist.f preset

// ---------------------------------------------------------------------------------------------- mc / mc_smear
int runMc(const std::string& dir, double tol) {
    HmsSpectrometer hms;
    std::ostringstream log;
    hms.setLog(&log);
    ScriptedRandom rng(ScriptedRandom::GaussMode::Scripted);
    std::vector<int> codes;
    std::vector<int> huts;

    std::vector<ColumnSpec> c;
    for (const char* n : {"ok", "stop_code", "stop_id", "reached_hut"}) c.push_back({n, 0.0, 0.0});
    // real*4 focal plane: exact equality required
    for (const char* n : {"x_fp", "dx_fp", "y_fp", "dy_fp"}) c.push_back({n, 0.0, 0.0});
    for (const char* n : {"dpp_out", "dxdz_out", "dydz_out", "y_out", "pathlen", "m2_out"}) c.push_back({n, tol, 0.0});
    c.push_back({"resmult", 0.0, 0.0});
    for (const char* n : {"xs_end", "ys_end", "dxdzs_end", "dydzs_end", "dpps_end"}) c.push_back({n, tol, 0.0});
    for (const char* n : {"n_gauss", "n_unif", "rng_err", "nsig_last", "row"}) c.push_back({n, 0.0, 0.0});

    int rc = simcl1::run("mc_hms (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        std::vector<double> g, u;
        for (int k = 1; r.has("g" + std::to_string(k)); ++k) g.push_back(r["g" + std::to_string(k)]);
        for (int k = 1; r.has("u" + std::to_string(k)); ++k) u.push_back(r["u" + std::to_string(k)]);
        rng.load(g, u);
        if (std::lround(r["spectr"]) != 1) throw std::runtime_error("not an HMS row");
        McHmsArgs a;
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
        a.using_HMScoll = std::lround(r["using_coll"]) != 0;
        a.x_fp = a.dx_fp = a.y_fp = a.dy_fp = kUnset;
        a.resmult = kUnset;
        a.pathlen = 0.;
        const McHmsDetail det = hms.mcHms(a, rng);
        codes.push_back(det.stop_code);
        huts.push_back(det.reached_hut ? 1 : 0);
        o["ok"] = a.ok_spec ? 1. : 0.;
        o["stop_code"] = det.stop_code;
        o["stop_id"] = -1.;
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
        const optics::Track& t = hms.track();
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

    // Per-reason agreement and the COMMON /spectrometer1/ counters summed over the table.
    simcl1::Table ref;
    std::string err;
    if (!ref.load(dir + "/output.txt", err)) return 1;
    const int ic = ref.col("stop_code"), ih = ref.col("reached_hut");
    std::map<int, long> nref, ncpp, nagree;
    long hutRef = 0;
    for (std::size_t i = 0; i < ref.rows.size(); ++i) {
        const int cr = static_cast<int>(std::lround(ref.rows[i][ic]));
        ++nref[cr];
        hutRef += std::lround(ref.rows[i][ih]);
        if (i < codes.size()) {
            ++ncpp[codes[i]];
            if (codes[i] == cr) ++nagree[cr];
        }
    }
    std::printf("  stop reason agreement (row by row):\n  %4s %-10s %8s %8s %8s\n", "code", "name", "ref", "cpp", "agree");
    bool okc = true;
    for (const auto& kv : nref) {
        const char* nm = kv.first == 0 ? "success" : hms::counterName(kv.first);
        std::printf("  %4d %-10s %8ld %8ld %8ld\n", kv.first, nm, kv.second, ncpp[kv.first], nagree[kv.first]);
        if (nagree[kv.first] != kv.second || ncpp[kv.first] != kv.second) okc = false;
    }
    const hms::HmsCounters& cn = hms.counters();
    std::printf("  counters: trials %ld (rows %zu), hut %ld (ref %ld), successes %ld (ref %ld)\n", cn[hms::kTrials],
                ref.rows.size(), cn[hms::kHut], hutRef, cn[hms::kSuccesses], nref[0]);
    if (cn[hms::kTrials] != static_cast<long>(ref.rows.size()) || cn[hms::kHut] != hutRef ||
        cn[hms::kSuccesses] != nref[0])
        okc = false;
    for (int k = hms::kSlitHor; k <= hms::kColl; ++k) {
        if (k == hms::kHut || k == hms::kSuccesses) continue;
        if (cn[k] != nref[k]) {
            std::printf("  counter %s = %ld, reference rows with that stop code = %ld\n", hms::counterName(k), cn[k],
                        nref[k]);
            okc = false;
        }
    }
    const auto& m = hms.stopCounters();
    if (m.size() != 22 || m.at("trials") != cn[hms::kTrials]) okc = false;
    if (!log.str().empty()) {
        std::printf("  unexpected messages:\n%s", log.str().c_str());
        okc = false;
    }
    std::printf("  stop reasons / counters: %s\n", okc ? "PASS" : "FAIL");
    return (rc == 0 && okc) ? 0 : 1;
}

// ---------------------------------------------------------------------------------------------- p12
// chk_p12.f:arms (hadron arm = HMS part), inputs rebuilt exactly as the driver does.
int runP12(const std::string& dir, double tol) {
    simcl1::Table in, out;
    std::string err;
    if (!in.load(dir + "/input.txt", err) || !out.load(dir + "/output.txt", err)) {
        std::printf("p12: SKIP (%s)\n", err.c_str());
        return 77;
    }
    HmsSpectrometer hms;
    Mt64Random rng(12345);
    const char* cols[8] = {"hsdelta", "hsyptar", "hsxptar", "hsytar", "hsxfp", "hsxpfp", "hsyfp", "hsypfp"};
    std::vector<simcl1::ColumnResult> res(8);
    for (int k = 0; k < 8; ++k) res[k].name = cols[k];
    int fail = 0;
    for (std::size_t i = 0; i < out.rows.size(); ++i) {
        Row o(out, i);
        const long ev = std::lround(o["event"]);
        std::size_t j = 0;
        while (j < in.rows.size() && std::lround(Row(in, j)["event"]) != ev) ++j;
        if (j == in.rows.size()) {
            std::printf("  event %ld missing in input\n", ev);
            return 1;
        }
        Row a(in, j);
        const double tx = o["tx"], ty = o["ty"], tz = o["tz"];
        const double pP = 3035.585;
        const double thP = std::fabs(36.594) / kDegrad;
        const double cP = std::cos(thP), sP = std::sin(thP);
        const double phiP = 3 * kPi / 2.;
        McHmsArgs m;
        m.p_spec = pP;
        m.th_spec = thP;
        m.dpp = a["hsdeltai"];
        m.dydz = a["hsyptari"];
        m.dxdz = a["hsxptari"];
        m.x = -ty;
        m.y = -tx * cP - tz * sP * std::sin(phiP);
        m.z = tz * cP + tx * sP * std::sin(phiP);
        m.x = m.x - m.z * m.dxdz;
        m.y = m.y - m.z * m.dydz;
        m.z = 0.0;
        m.fry = m.x;
        m.m2 = kMp * kMp;
        m.pathlen = 0.;
        hms.mcHms(m, rng);
        const double v[8] = {m.dpp, m.dydz, m.dxdz, m.y, m.x_fp, m.dx_fp, m.y_fp, m.dy_fp};
        const bool okh = std::lround(o["ok_hms"]) == 1;
        if (m.ok_spec != okh) {
            std::printf("  event %ld mode %ld: ok_hms C++ %d, Fortran %d\n", ev, std::lround(o["mode"]), m.ok_spec, okh);
            ++fail;
        }
        for (int k = 0; k < 8; ++k) {
            const double f = o[cols[k]];
            const double d = simcl1::relDiff(v[k], f);
            auto& rr = res[k];
            ++rr.rows;
            if (v[k] == f) ++rr.bitwise_equal;
            const double tk = k >= 4 ? 0.0 : tol;   // focal plane: real*4, exact
            if (d > tk) ++rr.failed;
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
    const bool ok = simcl1::report("mc_p12_check HMS arm (" + dir + ")", res, spec);
    return (ok && fail == 0) ? 0 : 1;
}

// ---------------------------------------------------------------------------------------------- dist
// Regularized upper incomplete gamma Q(a,x) (series / continued fraction), for chi2 p-values.
double gammq(double a, double x) {
    if (x <= 0) return 1.0;
    const double gln = std::lgamma(a);
    if (x < a + 1) {
        double ap = a, sum = 1.0 / a, del = sum;
        for (int n = 0; n < 10000; ++n) {
            ap += 1;
            del *= x / ap;
            sum += del;
            if (std::fabs(del) < std::fabs(sum) * 1e-15) break;
        }
        return 1.0 - sum * std::exp(-x + a * std::log(x) - gln);
    }
    double b = x + 1 - a, cc = 1.0 / 1e-300, d = 1.0 / b, h = d;
    for (int i = 1; i < 10000; ++i) {
        const double an = -i * (i - a);
        b += 2;
        d = an * d + b;
        if (std::fabs(d) < 1e-300) d = 1e-300;
        cc = b + an / cc;
        if (std::fabs(cc) < 1e-300) cc = 1e-300;
        d = 1.0 / d;
        const double del = d * cc;
        h *= del;
        if (std::fabs(del - 1.0) < 1e-15) break;
    }
    return std::exp(-x + a * std::log(x) - gln) * h;
}

// Kolmogorov distribution Q_KS(lambda).
double qks(double lam) {
    if (lam < 0.2) return 1.0;
    double s = 0, sign = 1;
    for (int j = 1; j <= 200; ++j) {
        s += sign * 2 * std::exp(-2 * j * j * lam * lam);
        sign = -sign;
    }
    return std::min(1.0, std::max(0.0, s));
}

struct RefHist {
    std::vector<double> lo, hi;
    std::vector<double> n;
    double under{0}, over{0};
};

bool loadHist(const std::string& path, RefHist& h) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    std::getline(f, line);   // header
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            auto pu = line.find("underflow"), po = line.find("overflow");
            if (pu != std::string::npos) h.under = std::atof(line.c_str() + pu + 9);
            if (po != std::string::npos) h.over = std::atof(line.c_str() + po + 8);
            continue;
        }
        std::istringstream is(line);
        double a, b, c;
        is >> a >> b >> c;
        h.lo.push_back(a);
        h.hi.push_back(b);
        h.n.push_back(c);
    }
    return !h.n.empty();
}

int runDist(const std::string& refdir, const std::string& cfg, long nray, std::uint64_t seed) {
    const std::string base = refdir;
    {
        std::ifstream t(base + "/stopfrac.txt");
        if (!t) {
            std::printf("dist: SKIP (no %s/stopfrac.txt)\n", base.c_str());
            return 77;
        }
    }
    const bool smear = cfg == "hms_smear";
    if (cfg != "hms_nosmear" && !smear) return 2;
    // drv_mcdist.f source for spectr 1
    const double p_spec = 3035.585, th = std::fabs(36.594) / kDegrad, m2 = kMp * kMp;
    const double dlo = -12., dhi = 12., ax = 0.09, ay = 0.05, xw = 0.1, yw = 1.0;

    HmsSpectrometer hms;
    Mt64Random rng(seed);
    std::map<std::string, std::vector<double>> v;
    std::map<int, long> ncode;
    for (long i = 0; i < nray; ++i) {
        const double d0 = dlo + (dhi - dlo) * rng.uniform();
        const double a0 = ax * (2. * rng.uniform() - 1.);
        const double b0 = ay * (2. * rng.uniform() - 1.);
        const double x0 = xw * (2. * rng.uniform() - 1.);
        const double y0 = yw * (2. * rng.uniform() - 1.);
        McHmsArgs a;
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
        const McHmsDetail det = hms.mcHms(a, rng);
        ++ncode[det.stop_code];
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

    int fail = 0;
    // --- stop fractions
    std::map<int, long> nref;
    long Nref = 0;
    {
        std::ifstream f(base + "/stopfrac.txt");
        std::string line;
        std::getline(f, line);
        while (std::getline(f, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') {
                if (line.find("# " + cfg + ":") != std::string::npos) {
                    auto p = line.find("N=");
                    Nref = std::atol(line.c_str() + p + 2);
                }
                continue;
            }
            std::istringstream is(line);
            std::string c, nm;
            int code, sid;
            long cnt;
            double fr;
            is >> c >> code >> nm >> sid >> cnt >> fr;
            if (c == cfg) nref[code] += cnt;
        }
    }
    std::map<int, bool> all;
    for (auto& kv : nref) all[kv.first] = true;
    for (auto& kv : ncode) all[kv.first] = true;
    std::printf("dist %s: C++ %ld rays (Mt64 seed %llu), Fortran %ld rays (RANLUX)\n", cfg.c_str(), nray,
                static_cast<unsigned long long>(seed), Nref);
    std::printf("  %4s %-10s %10s %10s %12s %12s %8s\n", "code", "name", "n_ref", "n_cpp", "frac_ref", "frac_cpp", "pull");
    double maxpull = 0;
    for (auto& kv : all) {
        const int code = kv.first;
        const double n1 = static_cast<double>(nref[code]), n2 = static_cast<double>(ncode[code]);
        const double N1 = static_cast<double>(Nref), N2 = static_cast<double>(nray);
        const double p = (n1 + n2) / (N1 + N2);
        const double se = std::sqrt(p * (1 - p) * (1 / N1 + 1 / N2));
        const double pull = se > 0 ? (n2 / N2 - n1 / N1) / se : 0.0;
        maxpull = std::max(maxpull, std::fabs(pull));
        std::printf("  %4d %-10s %10.0f %10.0f %12.6e %12.6e %8.2f%s\n", code,
                    code == 0 ? "success" : hms::counterName(code), n1, n2, n1 / N1, n2 / N2, pull,
                    std::fabs(pull) > 5 ? "  <-- |pull|>5" : "");
        if (std::fabs(pull) > 5) ++fail;
    }
    std::printf("  max |pull| = %.2f\n", maxpull);

    // --- histograms: binned two-sample chi2 and KS, quantiles
    std::map<std::string, std::vector<double>> qref;
    {
        std::ifstream f(base + "/summary.txt");
        std::string line;
        std::getline(f, line);
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream is(line);
            std::string c, var;
            double n;
            is >> c >> var >> n;
            if (c != cfg) continue;
            std::vector<double> q(13);
            for (double& x : q) is >> x;   // mean rms min max q01 q05 q25 q50 q75 q95 q99 lo hi
            qref[var] = q;
        }
    }
    std::printf("  %-9s %9s %9s %8s %10s %10s  %s\n", "variable", "chi2/ndf", "p_chi2", "KS_D", "p_KS", "d_mean/se",
                "quantile diffs q01..q99 in units of rms");
    const char* vars[] = {"x_fp",    "dx_fp",   "y_fp",     "dy_fp",    "dpp_rec", "dxdz_rec", "dydz_rec",
                          "y_rec",   "pathlen", "resmult",  "res_dpp",  "res_dxdz", "res_dydz", "res_y"};
    for (const char* var : vars) {
        RefHist h;
        if (!loadHist(base + "/hist_" + cfg + "_" + var + ".txt", h)) {
            std::printf("  %s: reference histogram missing\n", var);
            ++fail;
            continue;
        }
        std::vector<double>& x = v[var];
        const std::size_t nb = h.n.size();
        const double lo = h.lo.front(), hi = h.hi.back();
        std::vector<double> c(nb, 0.0);
        double under = 0, over = 0;
        // dist_stats.py conventions: numpy.histogram bins over [lo,hi] (last bin closed), underflow = (x<lo),
        // overflow = (x>=hi) — a value equal to hi is counted in the last bin AND in the overflow (e.g. resmult=2).
        for (double xv : x) {
            if (xv >= hi) ++over;
            if (xv < lo) {
                ++under;
            } else if (xv > hi) {
            } else {
                std::size_t b = static_cast<std::size_t>((xv - lo) / (hi - lo) * static_cast<double>(nb));
                if (b >= nb) b = nb - 1;
                // guard against edge rounding: use the reference bin edges
                while (b > 0 && xv < h.lo[b]) --b;
                while (b + 1 < nb && xv >= h.hi[b]) ++b;
                ++c[b];
            }
        }
        // reference counts incl. under/overflow
        std::vector<double> r = h.n, k = c;
        r.push_back(h.under);
        r.push_back(h.over);
        k.push_back(under);
        k.push_back(over);
        double N1 = 0, N2 = 0;
        for (std::size_t i = 0; i < r.size(); ++i) {
            N1 += r[i];
            N2 += k[i];
        }
        double chi2 = 0;
        int ndf = -1;
        for (std::size_t i = 0; i < r.size(); ++i) {
            if (r[i] + k[i] <= 0) continue;
            const double t = std::sqrt(N2 / N1) * r[i] - std::sqrt(N1 / N2) * k[i];
            chi2 += t * t / (r[i] + k[i]);
            ++ndf;
        }
        const double pchi = ndf > 0 ? gammq(0.5 * ndf, 0.5 * chi2) : 1.0;
        // binned KS (order: underflow, bins, overflow)
        double c1 = h.under / N1, c2 = under / N2, D = std::fabs(c1 - c2);
        for (std::size_t i = 0; i < nb; ++i) {
            c1 += h.n[i] / N1;
            c2 += c[i] / N2;
            D = std::max(D, std::fabs(c1 - c2));
        }
        const double ne = N1 * N2 / (N1 + N2);
        const double pks = qks((std::sqrt(ne) + 0.12 + 0.11 / std::sqrt(ne)) * D);
        // mean and quantiles
        std::vector<double> s = x;
        std::sort(s.begin(), s.end());
        auto quant = [&](double q) {   // numpy linear interpolation
            const double pos = q * static_cast<double>(s.size() - 1);
            const std::size_t i0 = static_cast<std::size_t>(std::floor(pos));
            const std::size_t i1 = std::min(i0 + 1, s.size() - 1);
            return s[i0] + (s[i1] - s[i0]) * (pos - static_cast<double>(i0));
        };
        double mean = 0;
        for (double xv : x) mean += xv;
        mean /= static_cast<double>(x.size());
        const std::vector<double>& q = qref[var];
        std::string qd;
        double dmean = 0;
        if (q.size() == 13 && q[1] > 0) {
            dmean = (mean - q[0]) / (q[1] * std::sqrt(1.0 / N1 + 1.0 / N2));
            const double qs[7] = {0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99};
            char buf[32];
            for (int i = 0; i < 7; ++i) {
                std::snprintf(buf, sizeof buf, " %+.4f", (quant(qs[i]) - q[4 + i]) / q[1]);
                qd += buf;
            }
        }
        const bool bad = pchi < 1e-4 || pks < 1e-4;
        std::printf("  %-9s %9.3f %9.3g %8.5f %10.3g %10.2f %s%s\n", var, ndf > 0 ? chi2 / ndf : 0.0, pchi, D, pks, dmean,
                    qd.c_str(), bad ? "  <-- p < 1e-4" : "");
        if (bad) ++fail;
    }
    std::printf("dist %s: %s (%d failures; criteria |pull| <= 5, p_chi2 and p_KS >= 1e-4)\n", cfg.c_str(),
                fail ? "FAIL" : "PASS", fail);
    return fail ? 1 : 0;
}

// ---------------------------------------------------------------------------------------------- selftest
int runSelftest() {
    int fail = 0;
    auto expect = [&](bool ok, const std::string& what) {
        std::printf("  %-78s %s\n", what.c_str(), ok ? "ok" : "FAIL");
        if (!ok) ++fail;
    };
    {   // lfit: double accumulators, float outputs, y == 0 ignored
        const float x[4] = {1.f, 2.f, 3.f, 4.f};
        const float y[4] = {1.1f, 0.f, 3.2f, 3.9f};
        float a, b, e;
        cern::lfit(x, y, 4, 0, a, b, e);
        // hand computation over points 1,3,4 in double
        const double X[3] = {1., 3., 4.}, Y[3] = {1.1f, 3.2f, 3.9f};
        double sx = 0, sy = 0;
        for (int i = 0; i < 3; ++i) {
            sx += X[i];
            sy += Y[i];
        }
        const double xm = sx / 3., ym = sy / 3.;
        double sxy = 0, sxx = 0, syy = 0;
        for (int i = 0; i < 3; ++i) {
            sxy += (X[i] - xm) * (Y[i] - ym);
            sxx += (X[i] - xm) * (X[i] - xm);
            syy += (Y[i] - ym) * (Y[i] - ym);
        }
        const float ea = static_cast<float>(sxy / sxx);
        const float eb = static_cast<float>(ym - static_cast<double>(ea) * xm);
        const float ee = static_cast<float>((syy - sxy * static_cast<double>(ea)) / (3. - 2.));
        expect(a == ea && b == eb && e == ee, "lfit: y=0 point ignored, double sums, float results");
        float a2, b2, e2;
        cern::lfit(x, y, 4, 1, a2, b2, e2);
        expect(a2 != a, "lfit KEY=1 uses the y=0 point");
        const float z[3] = {0.f, 0.f, 5.f};
        cern::lfit(x, z, 3, 0, a, b, e);
        expect(a == 0.f && b == 0.f && e == 0.f, "lfit: one usable point -> A=B=E=0");
    }
    expect(!hms::hitDipole(0., 0.) && hms::hitDipole(40., 0.) && hms::hitDipole(0., 21.) &&
               !hms::hitDipole(30., 13.),
           "hit_dipole: centre inside, |x|>34.29 or |y|>20.54 outside, rounded corner inside");
    {
        const double t1 = hms::pionCollAbsorb(1000., 0.315), t0 = hms::pionCollAbsorb(1000., 0.);
        expect(t0 == 1.0 && t1 > 0.9 && t1 < 1.0, "pion_coll_absorb: transmission(0 cm) = 1, (0.315 cm) in (0.9,1)");
    }
    {
        auto s = Registry<ISpectrometer>::create("hms");
        expect(s && s->name() == "hms", "Registry<ISpectrometer>::create(\"hms\")");
        SimcParameters p;
        Json cfg = {{"using_HMScoll", true}, {"ms_flag", false}};
        s->configure(cfg, p);
        Mt64Random rng(7);
        SpecInput in;
        in.p_spec = 3035.585;
        in.m2 = kMp * kMp;
        in.ms_flag = true;   // overridden by JSON ms_flag=false
        in.wcs_flag = false;
        SpecResult r = s->transport(in, rng);
        expect(r.ok && r.stop_reason.empty() && r.x_fp == 0. && std::fabs(r.pathlen - 2677.493) < 1e-9 &&
                   s->stopCounters().at("successes") == 1,
               "transport(): central ray accepted, fp = 0, pathlen 2677.493 cm");
    }
    {   // collimator (using_HMScoll, pions): smoke test, no Fortran table yet
        HmsSpectrometer h;
        std::ostringstream log;
        h.setLog(&log);
        Mt64Random rng(11);
        const double mpi2 = 139.57018 * 139.57018;
        long ok = 0, coll = 0;
        for (int i = 0; i < 20000; ++i) {
            McHmsArgs a;
            a.p_spec = 2000.;
            a.dpp = -8 + 16 * rng.uniform();
            a.dxdz = 0.09 * (2 * rng.uniform() - 1);
            a.dydz = 0.05 * (2 * rng.uniform() - 1);
            a.m2 = mpi2;
            a.fry = 0;
            a.using_HMScoll = true;
            h.mcHms(a, rng);
            ok += a.ok_spec;
        }
        coll = h.counters()[hms::kColl];
        std::printf("    coll run: successes %ld, coll stops %ld, slit_hor %ld (per-step counts)\n", ok, coll,
                    h.counters()[hms::kSlitHor]);
        expect(ok > 0 && coll > 0 && h.counters()[hms::kSlitHor] > 0, "collimator stepping runs (pions)");
    }
    {   // decay: pions with a short ctau decay (optics::project decay branch) and change m2
        HmsSpectrometer h;
        Mt64Random rng(13);
        const double mpi2 = 139.57018 * 139.57018;
        long decayed = 0;
        for (int i = 0; i < 2000; ++i) {
            McHmsArgs a;
            a.p_spec = 1000.;
            a.m2 = mpi2;
            a.decay_flag = true;
            a.ctau = 780.45;
            a.Mh2_final = mpi2;
            h.mcHms(a, rng);
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
            const std::uint64_t seed = argc > 5 ? std::strtoull(argv[5], nullptr, 10) : 20260930ULL;
            return runDist(argv[2], argv[3], n, seed);
        }
    } catch (const std::exception& e) {
        std::printf("exception: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
