// src/physics/sidis/F1F2In21.cpp — the A = 1 part of F1F2IN21_v1.0.f (the Christy-Bosted proton F1/F2 fit, resonances +
// non-resonant background), used by the SIDIS cross section. See
// include/simc/physics/sidis/F1F2In21.h. REAL(8) everywhere; operation order of the source (D-024); x**2 and x**2.
// -> x*x, x**1. -> x, real exponents -> std::pow (the pow call sites of F1F2IN21_v1.0.o).
#include "simc/physics/sidis/F1F2In21.h"

#include <cmath>
#include <string>

#include "simc/core/Errors.h"

namespace simc {
namespace sidis {

namespace {

inline double max0(double v) { return v > 0.0 ? v : 0.0; }   // MAX(0.0, v) as maxsd (NaN -> 0)

enum class Nucleon { proton, neutron };

// Fortran: F1F2IN21_v1.0.f:RESMODP (l.373-670) and RESMODN (l.672-959). The two routines are the same code apart from
// the nucleon mass, the 4th-region height (l.183-187 / 470-474) and the non-resonant exponents t1/t2 (l.596-631 /
// 883-916); the differences are selected by `n`.
double resmod(Nucleon n, int sf, double w2, double q2, const double* xv) {
    auto xval = [&](int i) { return xv[i - 1]; };
    const bool P = n == Nucleon::proton;
    double mass[8], width[8], height[8], rescoef[7][5], sigr[8], wdif[3];
    double intwidth[8], kr[8], kcmr[8], ppicmr[8], ppi2cmr[8], petacmr[8], epicmr[8], epi2cmr[8];
    double eetacmr[8], br[8][4], ang[8], pgam[8], pwid[8][4], x0[8], xpr[3];

    const double mp = P ? 0.9382727 : 0.939565;
    const double mpi = 0.134977;
    const double meta = 0.547862;
    const double mp2 = mp * mp;
    const double W = std::sqrt(w2);
    wdif[1] = W - (mp + mpi);
    wdif[2] = W - (mp + meta);
    const double q20 = xval(50);

    // single pion branching ratios
    br[1][1] = 1.00;   // P33(1232)
    br[2][1] = 0.45;   // S11(1535)
    br[3][1] = 0.60;   // D13(1520)
    br[4][1] = 0.65;   // F15(1680)
    br[5][1] = 0.60;   // S11(1650)
    br[6][1] = 0.65;   // P11(1440) roper
    br[7][1] = 0.60;   // F37(1950)
    // eta branching ratios
    br[1][3] = 0.0;
    br[2][3] = 0.40;
    br[3][3] = 0.08;
    br[4][3] = 0.0;
    br[5][3] = 0.20;
    br[6][3] = 0.0;
    br[7][3] = 0.0;
    // 2-pion branching ratios
    for (int i = 1; i <= 7; ++i) br[i][2] = 1. - br[i][1] - br[i][3];
    // meson angular momentum
    ang[1] = 1.;
    ang[2] = 0.;
    ang[3] = 2.;
    ang[4] = 3.;
    ang[5] = 0.;
    ang[6] = 1.;
    ang[7] = 3.;
    for (int i = 1; i <= 7; ++i) x0[i] = 0.160;   // resonance damping parameter (0.16 in RESMODN: the same double)
    if (sf == 2) x0[1] = 0.07;                     // different Delta mass for sigL
    for (int i = 1; i <= 7; ++i) br[i][2] = 1. - br[i][1] - br[i][3];

    // dip = 1./(1.+q2/1.15)**2. is computed and never used
    const double mon = 1. / (1. + q2 / 1.5);   // 1./(1.+q2/1.5)**1.

    const double xb = q2 / (q2 + w2 - mp2);
    xpr[1] = 1.00 + (w2 - (mp + mpi) * (mp + mpi)) / (q2 + q20);
    xpr[1] = 1. / xpr[1];
    xpr[2] = 1. + (w2 - (mp + meta) * (mp + meta)) / (q2 + q20);
    xpr[2] = 1. / xpr[2];
    if (W <= (mp + mpi)) xpr[1] = 1.0;
    if (W <= (mp + meta)) xpr[2] = 1.0;

    // kinematics needed for the threshold relativistic Breit-Wigner
    const double k = (w2 - mp2) / 2. / mp;
    const double kcm = (w2 - mp2) / 2. / W;
    const double epicm = (w2 + mpi * mpi - mp2) / 2. / W;
    const double ppicm = std::sqrt(max0(epicm * epicm - mpi * mpi));
    const double epi2cm = (w2 + (2. * mpi) * (2. * mpi) - mp2) / 2. / W;
    const double ppi2cm = std::sqrt(max0(epi2cm * epi2cm - (2. * mpi) * (2. * mpi)));
    const double eetacm = (w2 + meta * meta - mp2) / 2. / W;
    const double petacm = std::sqrt(max0(eetacm * eetacm - meta * meta));

    int num = 0;
    for (int i = 1; i <= 6; ++i) {   // resonance masses
        num = num + 1;
        mass[i] = xval(i);
    }
    for (int i = 1; i <= 6; ++i) {   // resonance widths
        num = num + 1;
        intwidth[i] = xval(num);
        width[i] = intwidth[i];
    }
    mass[7] = xval(47);
    intwidth[7] = xval(48);
    width[7] = intwidth[7];

    for (int i = 1; i <= 7; ++i) {
        kr[i] = (mass[i] * mass[i] - mp2) / 2. / mp;
        kcmr[i] = (mass[i] * mass[i] - mp2) / 2. / mass[i];
        epicmr[i] = (mass[i] * mass[i] + mpi * mpi - mp2) / 2. / mass[i];
        ppicmr[i] = std::sqrt(max0(epicmr[i] * epicmr[i] - mpi * mpi));
        epi2cmr[i] = (mass[i] * mass[i] + (2. * mpi) * (2. * mpi) - mp2) / 2. / mass[i];
        ppi2cmr[i] = std::sqrt(max0(epi2cmr[i] * epi2cmr[i] - (2. * mpi) * (2. * mpi)));
        eetacmr[i] = (mass[i] * mass[i] + meta * meta - mp2) / 2. / mass[i];
        petacmr[i] = std::sqrt(max0(eetacmr[i] * eetacmr[i] - meta * meta));

        // partial widths: 1-pion decay mode
        pwid[i][1] = intwidth[i] * std::pow(ppicm / ppicmr[i], 2. * ang[i] + 1.) *
                     std::pow((ppicmr[i] * ppicmr[i] + x0[i] * x0[i]) / (ppicm * ppicm + x0[i] * x0[i]), ang[i]);
        // 2-pion decay mode
        pwid[i][2] = intwidth[i] * std::pow(ppi2cm / ppi2cmr[i], 2. * ang[i] + 4.) *
                     std::pow((ppi2cmr[i] * ppi2cmr[i] + x0[i] * x0[i]) / (ppi2cm * ppi2cm + x0[i] * x0[i]), ang[i] + 2);
        pwid[i][2] = W / mass[i] * pwid[i][2];
        // eta decay mode, only for the S11's
        pwid[i][3] = 0.;
        if (i == 2 || i == 5) {
            pwid[i][3] = intwidth[i] * std::pow(petacm / petacmr[i], 2. * ang[i] + 1.) *
                         std::pow((petacmr[i] * petacmr[i] + x0[i] * x0[i]) / (petacm * petacm + x0[i] * x0[i]), ang[i]);
        }
        const double rk = kcm / kcmr[i];
        pgam[i] = rk * rk * (kcmr[i] * kcmr[i] + x0[i] * x0[i]) / (kcm * kcm + x0[i] * x0[i]);
        pgam[i] = intwidth[i] * pgam[i];
        width[i] = br[i][1] * pwid[i][1] + br[i][2] * pwid[i][2] + br[i][3] * pwid[i][3];
    }

    // resonance Q^2 dependence
    for (int i = 1; i <= 6; ++i) {
        for (int j = 1; j <= 4; ++j) {
            num = num + 1;
            rescoef[i][j] = xval(num);
        }
        if (sf == 1) {
            height[i] = rescoef[i][1] * (1. + rescoef[i][2] * q2 / (1. + rescoef[i][3] * q2)) * std::pow(mon, rescoef[i][4]);
        } else {
            height[i] = (rescoef[i][1] + rescoef[i][2] * q2) * std::exp(-1. * rescoef[i][3] * q2);
        }
        height[i] = height[i] * height[i];
    }
    if (sf == 2) {   // 4th resonance region
        height[7] = P ? (xval(16) + xval(20) * q2) * std::exp(-1.0 * xval(24) * q2)
                      : (xval(44) + xval(45) * q2) * std::exp(-1.0 * xval(46) * q2);
    } else {
        height[7] = P ? xval(49) * std::pow(mon, xval(45)) : xval(49) * mon;
    }
    height[7] = height[7] * height[7];
    // the non-resonant coefficients nr_coef(3,4) = xval(37..48) are copied and never used (num ends at 48)

    // Breit-Wigners for all resonances
    double sig_res = 0.0;
    for (int i = 1; i <= 7; ++i) {
        const double dd = w2 - mass[i] * mass[i];
        sigr[i] = width[i] * pgam[i] / (dd * dd + (mass[i] * width[i]) * (mass[i] * width[i]));
        sigr[i] = height[i] * kr[i] / k * kcmr[i] / kcm * sigr[i] / intwidth[i];
        sig_res = sig_res + sigr[i];
    }
    sig_res = sig_res * W;
    if (sf == 2) sig_res = sig_res * q2;

    // non-resonant background
    double sig_nr = 0.;
    if (sf == 1 && xpr[1] < 1.0) {
        const double A0 = xval(37) / std::pow(1.0 + q2 / xval(42), xval(43));
        const double t1 = P ? xval(38) * std::log(1.06 + q2) + xval(39) / std::log(1.06 + q2)
                            : xval(38) * std::log(1.05 + q2) + xval(39) / (1.05 + q2);
        const double t2 = xval(40) * std::pow(1.0 + q2 / xval(41), xval(44));
        if (xpr[1] <= 1.0)   // 1-pi threshold
            sig_nr = 389.4 * A0 * std::pow(1. - xpr[1], t1) * std::pow(xpr[1], t2);
        if (xpr[2] <= 1.0)   // eta threshold
            sig_nr = sig_nr + xval(46) * 389.4 * A0 * std::pow(1. - xpr[2], t1) * std::pow(xpr[2], t2);
    } else if (sf == 2 && xpr[1] < 1.0) {
        const double a = 1.0 + q2 / xval(39);
        const double A0 = xval(37) / (a * a);   // **2.0
        const double t1 = xval(38) / (1.0 + q2 / (xval(40))) + xval(32) * std::log(q2 + xval(36));
        const double t2 = P ? xval(41) : xval(41) / std::pow(1.00 + q2 / xval(42), xval(43));
        if (xpr[1] <= 1.0) sig_nr = sig_nr + 389.4 * A0 * xb * std::pow(1. - xpr[1], t1) * std::pow(xpr[1], t2);
    }

    double sig = sig_res + sig_nr;
    // FORTRAN-QUIRK: (w-mp) < wdif(1) = w-(mp+mpi) is never true: the intended threshold cut is a no-op
    if ((W - mp) < wdif[1]) sig = 0.0;
    return sig;
}

// rescsp DATA xval (F1F2IN21_v1.0.f:205-225)
const double kXvalP[100] = {
    0.12291E+01, 0.15173E+01, 0.15044E+01, 0.17100E+01, 0.16801E+01, 0.14312E+01, 0.12616E+00, 0.23000E+00,
    0.92594E-01, 0.90606E-01, 0.75000E-01, 0.35067E+00, 0.75729E+01, 0.56091E+01, 0.94606E+01, 0.20156E+01,
    0.66190E+01, 0.41732E+00, 0.23980E-01, 0.53136E+01, 0.63752E+00, 0.11484E+02, 0.69949E-01, 0.26191E+01,
    0.53603E-01, 0.65000E+02, 0.15351E+00, 0.20624E+01, 0.23408E+01, 0.16100E+02, 0.62414E+02, 0.17201E+01,
    0.23261E+00, 0.65000E+02, 0.23292E+01, 0.14980E+01, 0.23000E+00, 0.63385E+00, 0.19093E-01, 0.61061E-01,
    0.29146E-02, 0.54388E+00, 0.77997E+00, 0.28783E+00, 0.10605E+01, 0.69793E+00, 0.20009E+01, 0.57000E+00,
    0.41632E+01, 0.38427E+00, 0.10000E+01, 0.99842E+00, 0.98719E+00, 0.10168E+01, 0.98945E+00, 0.99594E+00,
    0.98799E+00, 0.10271E+01, 0.10650E+01, 0.97920E+00, 0.10152E+01, 0.99622E+00, 0.81011E+01, 0.10070E-02,
    0.14857E+01, 0.33445E+01, 0.31641E-09, 0.69755E+02, 0.55228E+01, 0.14438E+00, 0.60474E+01, 0.65395E-07,
    0.14129E+01, 0.58609E+00, 0.36220E+01, 0.92699E+00, 0.14418E+01, 0.86403E-02, 0.10001E-03, 0.75106E+00,
    0.76077E+00, 0.42272E+00, 0.55511E-11, 0.52486E+00, 0.58153E+00, 0.15798E+01, 0.50105E+00, 0.89149E+02,
    0.72789E+00, 0.24813E-01, -.61906E+00, 0.10000E+01, 0.00000E+00, 0.00000E+00, 0.68158E+03, 0.12429E+01,
    0.00000E+00, 0.00000E+00, 0.00000E+00, 0.10000E-05};

// RESCSN DATA xvaln (F1F2IN21_v1.0.f:334-354)
const double kXvalN[100] = {
    0.12291E+01, 0.15173E+01, 0.15044E+01, 0.17100E+01, 0.16801E+01, 0.14312E+01, 0.12616E+00, 0.23000E+00,
    0.92594E-01, 0.90606E-01, 0.75000E-01, 0.35067E+00, 0.69500E+01, 0.86633E+01, 0.11557E+02, 0.22138E+01,
    0.44886E+01, 0.20500E+03, 0.84433E+03, 0.31167E+01, 0.96301E+00, 0.14956E+00, 0.20761E-07, 0.10440E+01,
    0.40143E-03, 0.90028E+02, 0.75248E-01, 0.20532E+00, 0.12444E-01, 0.34469E+03, 0.19948E+00, 0.26925E+01,
    0.48635E+01, 0.86000E+02, 0.67813E+04, 0.44281E+02, 0.29548E+00, 0.65421E+00, 0.23787E-09, 0.51967E-01,
    0.39926E-08, 0.29960E+00, 0.97516E+00, 0.46934E-01, 0.14246E+03, 0.55801E+00, 0.19349E+01, 0.27400E+00,
    0.38891E+00, 0.40000E-02, 0.10108E+01, 0.97020E+00, 0.98248E+00, 0.97768E+00, 0.10425E+01, 0.10198E+01,
    0.97822E+00, 0.98239E+00, 0.10103E+01, 0.10076E+01, 0.10044E+01, 0.99687E+00, 0.16696E+01, 0.10721E-06,
    0.54114E+00, 0.11923E+04, 0.55938E+02, 0.95000E+03, 0.39840E+02, 0.22026E+03, 0.30498E+01, 0.24459E+00,
    0.95574E+00, 0.35596E+00, 0.21228E-05, 0.96696E+01, 0.27563E+01, 0.93027E-01, 0.33559E+02, 0.31207E-01,
    0.29020E+02, 0.86417E+00, 0.36471E-08, 0.99167E+00, 0.68124E+00, 0.10000E-01, 0.90227E-01, 0.40115E+01,
    0.29915E+01, 0.45929E-01, -.16758E+01, 0.78493E+01, 0.78184E+01, 0.42074E+01, 0.41179E-05, 0.80597E+00,
    0.00000E+00, 0.00000E+00, 0.10045E+01, 0.62364E+00};

// The T and L parameter vectors of rescsp/RESCSN (l.227-234 / 357-363): xvalL = xval(51:100) except the masses and
// widths (1-12) and xval(47), xval(48), which are the T ones.
void splitXval(const double* xv, double* x1, double* xL) {
    for (int i = 1; i <= 50; ++i) {
        x1[i - 1] = xv[i - 1];
        xL[i - 1] = xv[50 + i - 1];
        if (i <= 12) xL[i - 1] = x1[i - 1];
        if (i == 47 || i == 48) xL[i - 1] = x1[i - 1];
    }
}

}  // namespace

// Fortran: F1F2IN21_v1.0.f:RESMODP (l.373-670)
double resmodp(int sf, double w2, double q2, const double* xval) { return resmod(Nucleon::proton, sf, w2, q2, xval); }
// Fortran: F1F2IN21_v1.0.f:RESMODN (l.672-959)
double resmodn(int sf, double w2, double q2, const double* xval) { return resmod(Nucleon::neutron, sf, w2, q2, xval); }

// Fortran: F1F2IN21_v1.0.f:rescsp (l.189-248). F1, FL and R (l.242-244) are computed and not returned.
void rescsp(double w2, double q2, double& sigT, double& sigL) {
    double x1[50], xL[50];
    splitXval(kXvalP, x1, xL);
    sigT = resmodp(1, w2, q2, x1);
    sigL = resmodp(2, w2, q2, xL);
}

// Fortran: F1F2IN21_v1.0.f:RESCSN (l.325-368)
void rescsn(double w2, double q2, double& sigtn, double& sigln) {
    double x1[50], xL[50];
    splitXval(kXvalN, x1, xL);
    sigtn = resmodn(1, w2, q2, x1);
    sigln = resmodn(2, w2, q2, xL);
}

// Fortran: F1F2IN21_v1.0.f:SF (l.2345-2371)
NucleonSF sf(double w2, double q2) {
    const double mp = 0.938272;
    const double mp2 = mp * mp;
    const double pi = 3.14159;
    const double pi2 = pi * pi;
    const double alpha = 1 / 137.03599;
    const double x = q2 / (q2 + w2 - mp2);
    double sigTp, sigLp, sigTn, sigLn;
    rescsp(w2, q2, sigTp, sigLp);
    rescsn(w2, q2, sigTn, sigLn);
    NucleonSF r;
    r.F1p = sigTp / 0.3894e3 / pi2 / alpha / 8.0 * std::abs(w2 - mp2);
    r.F1n = sigTn / 0.3894e3 / pi2 / alpha / 8.0 * std::abs(w2 - mp2);
    r.FLp = sigLp * 2.0 * x / 0.3894e3 / pi2 / alpha / 8.0 * std::abs(w2 - mp2);
    r.FLn = sigLn * 2.0 * x / 0.3894e3 / pi2 / alpha / 8.0 * std::abs(w2 - mp2);
    r.F2p = (2. * x * r.F1p + r.FLp) / (1. + 4. * mp2 * x * x / q2);
    r.F2n = (2. * x * r.F1n + r.FLn) / (1. + 4. * mp2 * x * x / q2);
    return r;
}

// Fortran: F1F2IN21_v1.0.f:F1F2IN21 (l.24-89)
void f1f2in21(double Z, double A, double qsq, double wsq, double& F1, double& F2) {
    const int IA = static_cast<int>(A);   // int(): truncation
    const int IZ = static_cast<int>(Z);
    if (IA < 2) {
        const NucleonSF s = sf(wsq, qsq);
        if (IZ < 1) {   // neutron
            F1 = s.F1n;
            F2 = s.F2n;
        } else {        // proton
            F1 = s.F1p;
            F2 = s.F2p;
        }
        return;
    }
    // l.60-86: deuteron (rescsd) and A > 2 (SFCROSS); for IA == 2 with IZ != 1 the Fortran leaves F1, F2 unchanged.
    throw NotPorted("F1F2IN21 with A = " + std::to_string(A) + ", Z = " + std::to_string(Z) +
                    ": only the free nucleon branch (int(A) < 2) is ported; SIMC (semi_physics.f) never calls the others");
}

}  // namespace sidis
}  // namespace simc
