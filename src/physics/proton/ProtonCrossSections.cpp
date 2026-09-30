// src/physics/proton/ProtonCrossSections.cpp
//
// Fortran: simc_gfortran/physics_proton.f (sigep l.1-21, deForest l.25-133,
//          fofa_best_fit l.137-171, sigMott l.176-191); constants from constants.inc.
//
// Porting notes (D-010):
//  * SIMC is compiled with -fdefault-real-8: every literal (0.62, 1., 1.e-6 ...) is a double.
//  * Expressions are written in the Fortran order (left-associative); x**2 (integer power)
//    is x*x as gfortran emits it; x**3. (real power) is a libm pow() call -> std::pow.
//  * Division by 2., 4. (powers of two) is kept as written; results are identical to *0.5.
//  * No -ffast-math/FMA contraction may be used for this file (C++17 strict mode gives
//    -ffp-contract=off with GCC).
#include "simc/physics/proton/ProtonCrossSections.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "simc/core/SimcConstants.h"

namespace simc {
namespace proton {

using constants::alpha;
using constants::hbarc;
using constants::Mp2;

// Fortran: physics_proton.f:sigMott (l.176-191)
//   sig = (2.*alpha*hbarc*e0*cos(theta/2.)/Q2 )**2 ;  sigMott = sig*1.d4
// gfortran folds 2.*alpha*hbarc at compile time (objdump: ((K*e0)*cos(0.5*theta))/Q2);
// the left-associative C++ expression gives the same K and the same order.
double sigMott(double e0, double theta, double Q2) {
    const double t = 2. * alpha * hbarc * e0 * std::cos(theta / 2.) / Q2;
    const double sig = t * t;
    return sig * 1.e4;   // fm**2 --> microbarns
}

// Fortran: physics_proton.f:fofa_best_fit (l.137-171)
void fofa_best_fit(double qsquar, double& GE, double& GM) {
    const double mu_p = 2.793;

    // Q2 = -qsquar*(hbarc**2.)*1.e-6   (hbarc**2. is a folded PARAMETER expression = hbarc*hbarc)
    constexpr double hbarc2 = hbarc * hbarc;
    const double Q2 = -(qsquar * hbarc2 * 1.e-6);
    const double Q = std::sqrt(std::max(Q2, 0.));   // Fortran max(Q2,0.e0)

    const double Q3 = std::pow(Q, 3.);
    const double Q4 = std::pow(Q, 4.);
    const double Q5 = std::pow(Q, 5.);

    // Eqs 4, 5. NB: the Q2 terms use the SIGNED Q2 (timelike qsquar>0 gives Q=0 but Q2<0).
    double denom = 1. + 0.62 * Q + 0.68 * Q2 + 2.8 * Q3 + 0.83 * Q4;
    GE = 1. / denom;
    denom = 1. + 0.35 * Q + 2.44 * Q2 + 0.5 * Q3 + 1.04 * Q4 + 0.34 * Q5;
    GM = mu_p / denom;
    // (Eqs 6 are commented out in the Fortran and not ported.)
}

// Fortran: physics_proton.f:sigep (l.1-21). Units set by sigMott (microbarn/sr).
double sigep(const Event& vertex) {
    const double q4sq = vertex.Q2;
    double GE = 0, GM = 0;
    fofa_best_fit(-(q4sq / (hbarc * hbarc)), GE, GM);   // -q4sq/hbarc**2
    const double qmu4mp = q4sq / 4. / Mp2;
    const double W1p = GM * GM * qmu4mp;
    const double W2p = (GE * GE + GM * GM * qmu4mp) / (1.0 + qmu4mp);
    const double th = std::tan(vertex.e.theta / 2.);
    const double Wp = W2p + 2. * W1p * (th * th);
    // debug(5) print omitted (diagnostic only).
    return sigMott(vertex.e.E, vertex.e.theta, vertex.Q2) * vertex.e.E / vertex.Ein * Wp;
}

// Fortran: physics_proton.f:deForest (l.25-133)
// Flag 0: sigcc1; +1: sigcc2; -1: sigcc1 on-shell (Ebar -> E'-nu, qbar -> q).
// deForest's metric is (-1,1,1,1); the inner products below have regular signs and enter
// the structure functions with reversed signs (as in the Fortran comments).
double deForest(const Event& ev, const DeForestState& st) {
    const double Mh2 = st.Mh2;
    const int flag = st.deForest_flag;

    const double q4sq = -ev.Q2;
    const double q2 = ev.q * ev.q;
    double ebar, qbsq;
    if (flag >= 0) {
        ebar = std::sqrt(ev.Pm * ev.Pm + Mh2);
        const double d = ev.p.E - ebar;
        qbsq = d * d - q2;
    } else {
        ebar = ev.p.E - ev.nu;
        qbsq = q4sq;
    }

    const double uqup = ev.uq.x * ev.up.x + ev.uq.y * ev.up.y + ev.uq.z * ev.up.z;
    double sin_gamma = 1. - uqup * uqup;
    // Fortran l.74-78: clamp a (rounding-) negative value to 0 and warn.
    if (sin_gamma < 0) {
        char buf[128];
        std::snprintf(buf, sizeof buf, " WARNING: deForest came up with sin_gamma = %10.3f at event %10ld", sin_gamma,
                      st.nevent);
        (st.log ? *st.log : std::cout) << buf << '\n';
        sin_gamma = 0.0;
    }
    sin_gamma = std::sqrt(sin_gamma);
    double cos_phi = 0.0;
    if (sin_gamma != 0)
        cos_phi = (ev.uq.y * (ev.uq.y * ev.up.z - ev.uq.z * ev.up.y) - ev.uq.x * (ev.uq.z * ev.up.x - ev.uq.x * ev.up.z)) /
                  sin_gamma / std::sqrt(1. - ev.uq.z * ev.uq.z);
    if (std::abs(cos_phi) > 1.) {
        cos_phi = std::copysign(1.0, cos_phi);   // sign(1.0,cos_phi)
        // FORTRAN-QUIRK: physics_proton.f:84-85 tests (abs(cos_phi)-1.) > 1.e-10 AFTER cos_phi was
        // reset to +/-1, so the 'deForest give cos_phi' warning can never print. Reproduced: no warning.
    }

    double GE = 0, GM = 0;
    fofa_best_fit(q4sq / (hbarc * hbarc), GE, GM);
    const double qmu4mp = q4sq / 4. / Mp2;          // Mp2 (PARAMETER), not Mh2
    const double f1 = (GE - GM * qmu4mp) / (1.0 - qmu4mp);
    const double kf2 = (GM - GE) / (1.0 - qmu4mp);
    const double f1sq = f1 * f1;
    const double kf2_over_2m_allsq = kf2 * kf2 / 4. / Mh2;

    const double th = std::tan(ev.e.theta / 2.);
    const double tan2 = th * th;
    const double r = q4sq / q2;
    const double termC = r * r;
    const double termT = tan2 - q4sq / 2. / q2;
    const double termS = tan2 - r * (cos_phi * cos_phi);
    const double termI = (-q4sq / q2) * std::sqrt(tan2 - q4sq / q2) * cos_phi;

    double WC, WT, WS, WI;
    if (flag <= 0) {
        const double s = f1 + kf2;
        const double sumFF1 = s * s;
        const double sumFF2 = f1sq - qbsq * kf2 * kf2 / 4. / Mh2;
        const double eb_pE = ebar + ev.p.E;
        WC = (eb_pE * eb_pE) * sumFF2 - q2 * sumFF1;
        WT = -(2 * qbsq * sumFF1);
        WS = 4 * (ev.p.P * ev.p.P) * (sin_gamma * sin_gamma) * sumFF2;
        WI = -(4 * eb_pE * ev.p.P * sin_gamma * sumFF2);
    } else {
        const double pbarp = ebar * ev.p.E - ev.p.P * (ev.up.x * ev.Pmx + ev.up.y * ev.Pmy + ev.up.z * ev.Pmz);
        const double pbarq = ebar * ev.nu - ev.q * (ev.uq.x * ev.Pmx + ev.uq.y * ev.Pmy + ev.uq.z * ev.Pmz);
        const double pq = ev.p.E * ev.nu - ev.p.P * ev.q * (ev.up.x * ev.uq.x + ev.up.y * ev.uq.y + ev.up.z * ev.uq.z);
        const double qbarq = (ev.p.E - ebar) * ev.nu - q2;

        WC = (ebar * ev.p.E + (-pbarp + Mh2) / 2.) * f1sq - q2 * f1 * kf2 / 2. -
             ((-(pbarq * ev.p.E) - pq * ebar) * ev.nu + ebar * ev.p.E * q4sq + pbarq * pq - (-pbarp - Mh2) / 2. * q2) *
                 kf2_over_2m_allsq;
        WT = -((-pbarp + Mh2) * f1sq) - qbarq * f1 * kf2 + (2. * pbarq * pq + (-pbarp - Mh2) * q4sq) * kf2_over_2m_allsq;
        const double pPsg = ev.p.P * sin_gamma;
        WS = (pPsg * pPsg) * (f1sq - q4sq * kf2_over_2m_allsq);
        WI = pPsg * (-((ebar + ev.p.E) * f1sq) + ((-pbarq - pq) * ev.nu + (ebar + ev.p.E) * q4sq) * kf2_over_2m_allsq);
    }

    double allsum = termC * WC + termT * WT + termS * WS + termI * WI;
    if (flag <= 0) allsum = allsum / 4.0;
    // debug(5) print omitted (diagnostic only).
    return sigMott(ev.e.E, ev.e.theta, ev.Q2) * ev.p.P * allsum / ebar;   // microbarn*(MeV/sr)**2
}

}  // namespace proton
}  // namespace simc
