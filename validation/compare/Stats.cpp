// Stats.cpp - see Stats.h
#include "Stats.h"

#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace simcval {

double pullPValue(double z) {
    if (!std::isfinite(z)) return 0.0;
    return std::erfc(std::fabs(z) / std::sqrt(2.0));
}

double pull(double a, double sa, double b, double sb) {
    if (a == b) return 0.0;
    double s = std::sqrt(sa * sa + sb * sb);
    if (s <= 0) return (a > b) ? std::numeric_limits<double>::infinity()
                               : -std::numeric_limits<double>::infinity();
    return (a - b) / s;
}

bool binomialPull(int64_t kA, int64_t nA, int64_t kB, int64_t nB, double& fA, double& fB,
                  double& sA, double& sB, double& z) {
    fA = nA > 0 ? double(kA) / double(nA) : 0.0;
    fB = nB > 0 ? double(kB) / double(nB) : 0.0;
    sA = nA > 0 ? std::sqrt(std::max(0.0, fA * (1 - fA)) / nA) : 0.0;
    sB = nB > 0 ? std::sqrt(std::max(0.0, fB * (1 - fB)) / nB) : 0.0;
    if (kA == 0 && kB == 0) { z = 0; return false; }
    if (nA <= 0 || nB <= 0) { z = 0; return false; }
    double pp = double(kA + kB) / double(nA + nB);
    double s = std::sqrt(pp * (1 - pp) * (1.0 / nA + 1.0 / nB));
    if (fA == fB) z = 0;
    else z = (s > 0) ? (fA - fB) / s : std::numeric_limits<double>::infinity();
    return true;
}

bool contingencyChi2(const std::vector<int64_t>& a, const std::vector<int64_t>& b, double& chi2,
                     int& ndf, double& p) {
    double NA = 0, NB = 0;
    int K = 0;
    for (size_t k = 0; k < a.size(); ++k) {
        NA += a[k];
        NB += b[k];
        if (a[k] + b[k] > 0) ++K;
    }
    chi2 = 0;
    ndf = K - 1;
    p = 1;
    if (K < 2 || NA <= 0 || NB <= 0) return false;
    double N = NA + NB;
    for (size_t k = 0; k < a.size(); ++k) {
        double t = double(a[k] + b[k]);
        if (t <= 0) continue;
        double eA = t * NA / N, eB = t * NB / N;
        chi2 += (a[k] - eA) * (a[k] - eA) / eA + (b[k] - eB) * (b[k] - eB) / eB;
    }
    p = TMath::Prob(chi2, ndf);
    return true;
}

void perTriedMean(double S1, double S2, double N, double& Y, double& sY) {
    if (N <= 0) { Y = 0; sY = 0; return; }
    Y = S1 / N;
    sY = std::sqrt(std::max(0.0, S2 - S1 * S1 / N)) / N;
}

void weightedKS(const std::vector<double>& xa, const std::vector<double>& wa,
                const std::vector<double>& xb, const std::vector<double>& wb, double& D,
                double& p) {
    D = 0;
    p = 1;
    const size_t na = xa.size(), nb = xb.size();
    if (na == 0 || nb == 0) { p = (na == nb) ? 1.0 : 0.0; return; }
    auto W = [](const std::vector<double>& w, size_t i) { return w.empty() ? 1.0 : w[i]; };
    std::vector<size_t> ia(na), ib(nb);
    std::iota(ia.begin(), ia.end(), 0);
    std::iota(ib.begin(), ib.end(), 0);
    std::sort(ia.begin(), ia.end(), [&](size_t i, size_t j) { return xa[i] < xa[j]; });
    std::sort(ib.begin(), ib.end(), [&](size_t i, size_t j) { return xb[i] < xb[j]; });
    double SA = 0, SA2 = 0, SB = 0, SB2 = 0;
    for (size_t i = 0; i < na; ++i) { double w = W(wa, i); SA += w; SA2 += w * w; }
    for (size_t i = 0; i < nb; ++i) { double w = W(wb, i); SB += w; SB2 += w * w; }
    if (SA <= 0 || SB <= 0) { p = 0; return; }
    size_t i = 0, j = 0;
    double ca = 0, cb = 0;
    while (i < na || j < nb) {
        double x;
        if (j >= nb || (i < na && xa[ia[i]] <= xb[ib[j]])) x = xa[ia[i]];
        else x = xb[ib[j]];
        while (i < na && xa[ia[i]] == x) { ca += W(wa, ia[i]); ++i; }
        while (j < nb && xb[ib[j]] == x) { cb += W(wb, ib[j]); ++j; }
        D = std::max(D, std::fabs(ca / SA - cb / SB));
    }
    double neA = SA * SA / SA2, neB = SB * SB / SB2;
    double ne = neA * neB / (neA + neB);
    p = (D == 0) ? 1.0 : TMath::KolmogorovProb(D * std::sqrt(ne));
}

}  // namespace simcval
