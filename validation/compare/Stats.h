// Stats.h - statistical primitives used by simc_compare. All formulas are
// documented in validation/README.md ("Statistics").
#pragma once

#include <cstdint>
#include <vector>

namespace simcval {

// Two-sided p-value of a standard-normal pull.
double pullPValue(double pull);

// Pull of two independent measurements; returns 0 if a==b (also when both
// sigmas are 0), +-inf if the sigmas are 0 and a!=b.
double pull(double a, double sa, double b, double sb);

// Two-proportion z test with pooled proportion p = (kA+kB)/(nA+nB):
//   sigma^2 = p(1-p)(1/nA + 1/nB),  pull = (kA/nA - kB/nB)/sigma.
// Returns false if kA==kB==0 (test not applicable).
bool binomialPull(int64_t kA, int64_t nA, int64_t kB, int64_t nB, double& fA, double& fB,
                  double& sA, double& sB, double& z);

// 2xK homogeneity chi2 (contingency table) over categories with nA_k+nB_k>0.
// ndf = K'-1. Returns false if fewer than 2 populated categories.
bool contingencyChi2(const std::vector<int64_t>& a, const std::vector<int64_t>& b, double& chi2,
                     int& ndf, double& p);

// Mean per tried event of x_i (x_i=0 for events not in the tree):
//   Y = S1/N,  sigma^2 = (S2 - S1^2/N)/N^2  with S1=sum x, S2=sum x^2 over the tree.
void perTriedMean(double S1, double S2, double N, double& Y, double& sY);

// Unbinned two-sample weighted Kolmogorov-Smirnov test.
// D = sup |F_A(x) - F_B(x)| with weighted ECDFs; effective sizes
// n_eff = (sum w)^2 / sum w^2;  p = KolmogorovProb(D * sqrt(nA nB/(nA+nB))).
// Weights may be empty (unweighted). Ties are handled by stepping over equal x.
void weightedKS(const std::vector<double>& xa, const std::vector<double>& wa,
                const std::vector<double>& xb, const std::vector<double>& wb, double& D,
                double& p);

}  // namespace simcval
