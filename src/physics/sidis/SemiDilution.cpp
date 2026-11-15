// src/physics/sidis/SemiDilution.cpp — semi_dilution.f. DEAD CODE in SIMC: the only call, event.f:1524,
// is commented out. Ported for completeness (cheap) and covered by a unit test only; no plug-in calls it.
// REAL(8) arithmetic in the source order (D-024); sv**3 = sv*(sv*sv) (powi).
#include <cmath>
#include <ostream>

#include "simc/core/SimcConstants.h"
#include "simc/physics/sidis/SidisCrossSection.h"

namespace simc {
namespace sidis {

// Fortran: semi_dilution.f:semi_dilution (l.1-270)
double semiDilution(const Event& vertex, bool doing_semipi, bool doing_semika, bool doing_hplus, bool doing_deutsemi,
                    Cteq5Pdf& pdf, SemiDilutionState& st, std::ostream* log) {
    using constants::alpha;
    using constants::hbarc;
    using constants::Mp;
    using constants::pi;
    auto sq = [](double x) { return x * x; };
    constexpr double qu = 2. / 3., qd = -1. / 3.;
    constexpr double b = 3.76;          // GeV^-2
    constexpr double lambda = 0.227;    // 0.227 GeV for NLO
    constexpr double Q2zero = 2.0;      // GeV^2 for u,d,s,g

    if (st.first) {
        if (!pdf.loaded()) throw std::logic_error("semi_dilution: CTEQ5 table not loaded (SetCtq5)");
        st.first = false;
    }
    const double mp = Mp;
    const double s = (2. * vertex.Ein * mp + mp * mp) / 1.e6;   // convert to GeV2
    double xbj = vertex.Q2 / 2. / mp / vertex.nu;
    if (xbj > 1.0) {
        if (log) *log << " XBj is too large! " << xbj << "\n";
        xbj = 1.0;
    }
    const double y = vertex.nu / vertex.Ein;
    const double Q2gev = vertex.Q2 / 1.e6;
    double Qgev = std::sqrt(Q2gev);
    const double pt2gev = vertex.pt2 / 1.e6;

    const double u = pdf.pdf(1, xbj, Qgev, log);
    const double ubar = pdf.pdf(-1, xbj, Qgev, log);
    const double d = pdf.pdf(2, xbj, Qgev, log);
    const double dbar = pdf.pdf(-2, xbj, Qgev, log);

    // Binnewies et al. (PRD 52 p.4947): sv = ln( ln(Q2/Lambda^2)/ln(Q2_0/Lambda^2) )
    const double sv = std::log(std::log(Q2gev / (lambda * lambda)) / std::log(Q2zero / (lambda * lambda)));
    const double sv2 = sv * sv, sv3 = sv * sv2;
    if (doing_semipi) {
        st.N = 1.150 - 1.522 * sv + 1.378 * sv2 - 0.527 * sv3;
        st.a1 = -0.740 - 1.680 * sv + 1.546 * sv2 - 0.596 * sv3;
        st.a2 = 1.430 + 0.543 * sv - 0.023 * sv2;
    } else if (doing_semika) {
        if (st.kfirst) {
            if (log) *log << " Kaon FF in semi_dilution not iomplemented yet -using pion values\n";
            st.kfirst = false;
        }
        st.N = 1.150 - 1.522 * sv + 1.378 * sv2 - 0.527 * sv3;
        st.a1 = -0.740 - 1.680 * sv + 1.546 * sv2 - 0.596 * sv3;
        st.a2 = 1.430 + 0.543 * sv - 0.023 * sv2;
    }
    const double z = vertex.zhad;
    const double D_sum = st.N * std::pow(z, st.a1) * std::pow(1.0 - z, st.a2);
    // Ratio of D-/D+ from P. Geiger's thesis (HERMES)
    const double R_D = std::pow(1.0 - z, 0.083583) / std::pow(1.0 + z, 1.9838);
    const double D_fav = D_sum / (1.0 + R_D);
    const double D_unfav = D_sum / (1.0 + 1.0 / R_D);

    const double qu2 = sq(qu), qd2 = sq(qd);
    double sum_sq = qu2 * (u + ubar) + qd2 * (d + dbar);
    const double sum_sq_prot = sum_sq;
    const double sum_sq_neut = qu2 * (d + dbar) + qd2 * (u + ubar);
    if (doing_deutsemi) sum_sq = sum_sq + qu2 * (d + dbar) + qd2 * (u + ubar);

    double dsigdz;
    if (sum_sq > 0.) {
        if (doing_hplus)
            dsigdz = (qu2 * u * D_fav + qu2 * ubar * D_unfav + qd2 * d * D_unfav + qd2 * dbar * D_fav) / sum_sq;
        else
            dsigdz = (qu2 * u * D_unfav + qu2 * ubar * D_fav + qd2 * d * D_fav + qd2 * dbar * D_unfav) / sum_sq;
        if (doing_deutsemi) {
            if (doing_hplus)
                dsigdz = dsigdz + (qu2 * d * D_fav + qu2 * dbar * D_unfav + qd2 * u * D_unfav + qd2 * ubar * D_fav) / sum_sq;
            else
                dsigdz = dsigdz + (qu2 * d * D_unfav + qu2 * dbar * D_fav + qd2 * u * D_fav + qd2 * ubar * D_unfav) / sum_sq;
        }
    } else {
        dsigdz = 0.0;
    }
    double dsigdz_prot, dsigdz_neut;
    if (sum_sq_prot > 0.) {
        if (doing_hplus) {
            dsigdz_prot = (qu2 * u * D_fav + qu2 * ubar * D_unfav + qd2 * d * D_unfav + qd2 * dbar * D_fav) / sum_sq_prot;
            dsigdz_neut = (qu2 * d * D_fav + qu2 * dbar * D_unfav + qd2 * u * D_unfav + qd2 * ubar * D_fav) / sum_sq_neut;
        } else {
            dsigdz_prot = (qu2 * u * D_unfav + qu2 * ubar * D_fav + qd2 * d * D_fav + qd2 * dbar * D_unfav) / sum_sq_prot;
            dsigdz_neut = (qu2 * d * D_unfav + qu2 * dbar * D_fav + qd2 * u * D_fav + qd2 * ubar * D_unfav) / sum_sq_neut;
        }
    } else {
        dsigdz_prot = 0.0;
        dsigdz_neut = 0.0;
    }

    const double eE = vertex.e.E;
    const double nu = vertex.nu;
    const double hb2 = sq(hbarc / 1000.);
    const double sighad = dsigdz * b * std::exp(-b * pt2gev) / 2. / pi;
    const double sige = sq(alpha) * (1. + sq(1. - y)) * (eE / 1000.) / (s * xbj * sq(y) * (mp / 1000.) * (nu / 1000.)) * sum_sq;
    const double sigsemi = sige * sighad * hb2 * 10000.0;

    const double sighad_prot = dsigdz_prot * b * std::exp(-b * pt2gev) / 2. / pi;
    const double sige_prot =
        sq(alpha) * (1. + sq(1. - y)) * (eE / 1000.) / (s * xbj * sq(y) * (mp / 1000.) * (nu / 1000.)) * sum_sq_prot;
    const double sigsemi_prot = sige_prot * sighad_prot * hb2 * 10000.0;

    const double sighad_neut = dsigdz_neut * b * std::exp(-b * pt2gev) / 2. / pi;
    const double sige_neut =
        sq(alpha) * (1. + sq(1. - y)) * (eE / 1000.) / (s * xbj * sq(y) * (mp / 1000.) * (nu / 1000.)) * sum_sq_neut;
    const double sigsemi_neut = sige_neut * sighad_neut * hb2 * 10000.0;

    const double N15 = 8.0 * sigsemi_neut + 7.0 * sigsemi_prot;
    const double He4 = 2.0 * sigsemi_neut + 2.0 * sigsemi_prot;
    const double Al = 14.0 * sigsemi_neut + 13.0 * sigsemi_prot;
    const double Cu = 35. * sigsemi_neut + 29.0 * sigsemi_prot;
    const double Ni = 31.0 * sigsemi_neut + 28.0 * sigsemi_prot;

    double dilution;
    if (sigsemi_prot > 0.0 && sigsemi_neut > 0.0)
        dilution = sigsemi / (sigsemi + 0.3333 * N15 + 0.3165 * He4 + 0.0144 * Al + 0.0031 * Cu + 0.0013 * Ni);
    else
        dilution = 999;
    return dilution;
}

}  // namespace sidis
}  // namespace simc
