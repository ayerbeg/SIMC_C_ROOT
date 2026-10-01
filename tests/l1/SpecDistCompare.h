// tests/l1/SpecDistCompare.h — shared by test_l1_hms / test_l1_shms (tasks P6.2/P6.3).
//
// Compares a C++ run of N rays with the Fortran real-RNG run in reference/functions/mc_distributions
// (fortran_harness/spectrometers/dist_stats.py output):
//   * stop fractions per (stop code, stop_id): two-sample binomial pulls;
//   * per variable: binned two-sample chi2 and binned KS against the 200-bin reference histogram (dist_stats.py
//     conventions: numpy bins on [lo,hi] with the last bin closed, underflow x<lo, overflow x>=hi), the mean
//     difference in units of its standard error, and the quantiles q01..q99 in units of the rms.
// Pass: |pull| <= 5, p_chi2 >= 1e-4, p_KS >= 1e-4.
#ifndef SIMC_TESTS_L1_SPECDISTCOMPARE_H
#define SIMC_TESTS_L1_SPECDISTCOMPARE_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace simcl1 {
namespace dist {

/// Regularized upper incomplete gamma Q(a,x), for chi2 p-values.
inline double gammq(double a, double x) {
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

/// Kolmogorov distribution Q_KS(lambda).
inline double qks(double lam) {
    if (lam < 0.2) return 1.0;
    double s = 0, sign = 1;
    for (int j = 1; j <= 200; ++j) {
        s += sign * 2 * std::exp(-2 * j * j * lam * lam);
        sign = -sign;
    }
    return std::min(1.0, std::max(0.0, s));
}

struct RefHist {
    std::vector<double> lo, hi, n;
    double under{0}, over{0};
};

inline bool loadHist(const std::string& path, RefHist& h) {
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

using CodeKey = std::pair<int, int>;   ///< (stop code, stop_id)

/// Returns the number of failures (0 = pass); 77 if the reference is missing.
inline int compare(const std::string& base, const std::string& cfg, long nray, const std::map<CodeKey, long>& ncpp,
                   std::map<std::string, std::vector<double>>& v, const std::vector<std::string>& vars,
                   const std::function<std::string(int)>& codeName) {
    std::ifstream t(base + "/stopfrac.txt");
    if (!t) {
        std::printf("dist: SKIP (no %s/stopfrac.txt)\n", base.c_str());
        return 77;
    }
    int fail = 0;
    std::map<CodeKey, long> nref;
    long Nref = 0;
    {
        std::string line;
        std::getline(t, line);
        while (std::getline(t, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') {
                if (line.find("# " + cfg + ":") != std::string::npos) Nref = std::atol(line.c_str() + line.find("N=") + 2);
                continue;
            }
            std::istringstream is(line);
            std::string c, nm;
            int code, sid;
            long cnt;
            double fr;
            is >> c >> code >> nm >> sid >> cnt >> fr;
            if (c == cfg) nref[{code, sid}] += cnt;
        }
    }
    if (Nref == 0) {
        std::printf("dist: config %s not in %s/stopfrac.txt\n", cfg.c_str(), base.c_str());
        return 1;
    }
    std::map<CodeKey, bool> all;
    for (auto& kv : nref) all[kv.first] = true;
    for (auto& kv : ncpp) all[kv.first] = true;
    std::printf("dist %s: C++ %ld rays (Mt64), Fortran %ld rays (RANLUX)\n", cfg.c_str(), nray, Nref);
    std::printf("  %4s %7s %-10s %10s %10s %12s %12s %8s\n", "code", "stop_id", "name", "n_ref", "n_cpp", "frac_ref",
                "frac_cpp", "pull");
    double maxpull = 0;
    for (auto& kv : all) {
        const double n1 = nref.count(kv.first) ? static_cast<double>(nref.at(kv.first)) : 0.;
        const double n2 = ncpp.count(kv.first) ? static_cast<double>(ncpp.at(kv.first)) : 0.;
        const double N1 = static_cast<double>(Nref), N2 = static_cast<double>(nray);
        const double p = (n1 + n2) / (N1 + N2);
        const double se = std::sqrt(p * (1 - p) * (1 / N1 + 1 / N2));
        const double pull = se > 0 ? (n2 / N2 - n1 / N1) / se : 0.0;
        maxpull = std::max(maxpull, std::fabs(pull));
        std::printf("  %4d %7d %-10s %10.0f %10.0f %12.6e %12.6e %8.2f%s\n", kv.first.first, kv.first.second,
                    codeName(kv.first.first).c_str(), n1, n2, n1 / N1, n2 / N2, pull,
                    std::fabs(pull) > 5 ? "  <-- |pull|>5" : "");
        if (std::fabs(pull) > 5) ++fail;
    }
    std::printf("  max |pull| = %.2f\n", maxpull);

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
    for (const std::string& var : vars) {
        RefHist h;
        if (!loadHist(base + "/hist_" + cfg + "_" + var + ".txt", h)) {
            std::printf("  %s: reference histogram missing\n", var.c_str());
            ++fail;
            continue;
        }
        std::vector<double>& x = v[var];
        if (x.empty()) {
            std::printf("  %s: no C++ values\n", var.c_str());
            ++fail;
            continue;
        }
        const std::size_t nb = h.n.size();
        const double lo = h.lo.front(), hi = h.hi.back();
        std::vector<double> c(nb, 0.0);
        double under = 0, over = 0;
        // dist_stats.py: a value equal to hi is in the last bin AND in the overflow (e.g. resmult = 2)
        for (double xv : x) {
            if (xv >= hi) ++over;
            if (xv < lo) {
                ++under;
            } else if (xv <= hi) {
                std::size_t b = static_cast<std::size_t>((xv - lo) / (hi - lo) * static_cast<double>(nb));
                if (b >= nb) b = nb - 1;
                while (b > 0 && xv < h.lo[b]) --b;
                while (b + 1 < nb && xv >= h.hi[b]) ++b;
                ++c[b];
            }
        }
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
            const double tt = std::sqrt(N2 / N1) * r[i] - std::sqrt(N1 / N2) * k[i];
            chi2 += tt * tt / (r[i] + k[i]);
            ++ndf;
        }
        const double pchi = ndf > 0 ? gammq(0.5 * ndf, 0.5 * chi2) : 1.0;
        double c1 = h.under / N1, c2 = under / N2, D = std::fabs(c1 - c2);
        for (std::size_t i = 0; i < nb; ++i) {
            c1 += h.n[i] / N1;
            c2 += c[i] / N2;
            D = std::max(D, std::fabs(c1 - c2));
        }
        const double ne = N1 * N2 / (N1 + N2);
        const double pks = qks((std::sqrt(ne) + 0.12 + 0.11 / std::sqrt(ne)) * D);
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
        std::printf("  %-9s %9.3f %9.3g %8.5f %10.3g %10.2f %s%s\n", var.c_str(), ndf > 0 ? chi2 / ndf : 0.0, pchi, D,
                    pks, dmean, qd.c_str(), bad ? "  <-- p < 1e-4" : "");
        if (bad) ++fail;
    }
    std::printf("dist %s: %s (%d failures; criteria |pull| <= 5, p_chi2 and p_KS >= 1e-4)\n", cfg.c_str(),
                fail ? "FAIL" : "PASS", fail);
    return fail;
}

}  // namespace dist
}  // namespace simcl1

#endif  // SIMC_TESTS_L1_SPECDISTCOMPARE_H
