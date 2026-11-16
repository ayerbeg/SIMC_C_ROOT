/**
 * @file
 * @ingroup physics
 * @brief Inelastic nucleon structure functions of M. E. Christy (F1F221), the free-nucleon part used by SIMC.
 *
 * Inelastic nucleon structure functions of M. E. Christy (F1F221 version 1.0 prerelease, 2021-2022), the part of
 * F1F2IN21_v1.0.f that SIMC reaches (task P12.1). semi_physics.f:peepiX (l.531-540) calls F1F2IN21 with A = 1 only
 * (Z = 1 proton, Z = 0 neutron), so the reachable code is:
 *   f1f2in21  F1F2IN21_v1.0.f:24-89     IA = int(A) < 2 branch (the deuteron and A > 2 branches throw NotPorted)
 *   sf        F1F2IN21_v1.0.f:2345-2371 reduced cross sections -> F1, FL, F2 of proton and neutron
 *   rescsp    F1F2IN21_v1.0.f:189-248   proton sigma_T, sigma_L from the fit parameters (DATA xval)
 *   rescsn    F1F2IN21_v1.0.f:325-368   neutron sigma_T, sigma_L (DATA xvaln)
 *   resmodp   F1F2IN21_v1.0.f:373-670   proton resonance + non-resonant fit, sf = 1 (T) or 2 (L)   [ub/sr/GeV]
 *   resmodn   F1F2IN21_v1.0.f:672-959   neutron
 * peepiX also calls SQESUB once (semi_physics.f:288-296), which only reads the deuteron smearing tables
 * f1f2tables/f1*-cdbonn.dat into SAVEd arrays used by the deuteron routines (never by sf): no numerical effect on SIMC,
 * not ported (documented in FORTRAN_DIFFERENCES.md).
 *
 * Numerics: F1F2IN21_v1.0.f is compiled with -fdefault-real-8; the routines above contain no D-exponent literal, so
 * everything is REAL(8) (no binary128 in resmodp/resmodn/rescsp/rescsn/sf in F1F2IN21_v1.0.o). x**2 and x**2. are x*x,
 * x**1. is x, every other real exponent is libm pow; operation order as in the source. Units: GeV, GeV^2.
 */
#ifndef SIMC_PHYSICS_SIDIS_F1F2IN21_H
#define SIMC_PHYSICS_SIDIS_F1F2IN21_H

namespace simc {
namespace sidis {

/// @addtogroup physics
/// @{

/// Fortran: RESMODP(sf, w2, q2, xval, sig); xval(1:50) = xval[0..49]
/// @param sf 1 = transverse, 2 = longitudinal.
/// @param w2 W^2 [GeV^2].
/// @param q2 Q2 [GeV^2].
/// @param xval the 50 fit parameters.
/// @return the cross section sig of the routine.
double resmodp(int sf, double w2, double q2, const double* xval);
/// Fortran: RESMODN(sf, w2, q2, xval, sig)
/// @param sf 1 = transverse, 2 = longitudinal.
/// @param w2 W^2 [GeV^2].
/// @param q2 Q2 [GeV^2].
/// @param xval the 50 fit parameters.
/// @return the cross section sig of the routine.
double resmodn(int sf, double w2, double q2, const double* xval);
/// Fortran: rescsp(W2, Q2, sigT, sigL)
/// @param w2 W^2 [GeV^2].
/// @param q2 Q2 [GeV^2].
/// @param[out] sigT transverse proton cross section.
/// @param[out] sigL longitudinal proton cross section.
void rescsp(double w2, double q2, double& sigT, double& sigL);
/// Fortran: RESCSN(w2, q2, sigtn, sigln)
/// @param w2 W^2 [GeV^2].
/// @param q2 Q2 [GeV^2].
/// @param[out] sigtn transverse neutron cross section.
/// @param[out] sigln longitudinal neutron cross section.
void rescsn(double w2, double q2, double& sigtn, double& sigln);

/// Proton and neutron structure functions F1, FL, F2.
struct NucleonSF {
    double F1p{0} /**< @copydoc F2n */, FLp{0} /**< @copydoc F2n */, F2p{0} /**< @copydoc F2n */, F1n{0} /**< @copydoc F2n */, FLn{0} /**< @copydoc F2n */, F2n{0};  ///< F1, FL, F2 of the proton and of the neutron.
};
/// Fortran: SF(w2, q2, F1p, FLp, F2p, F1n, FLn, F2n)
/// @param w2 W^2 [GeV^2].
/// @param q2 Q2 [GeV^2].
/// @return the structure functions.
NucleonSF sf(double w2, double q2);

/// Fortran: F1F2IN21(Z, A, QSQ, WSQ, F1, F2) for int(A) < 2 (free proton for int(Z) >= 1, free neutron otherwise).
/// int(A) >= 2 (deuteron rescsd, A > 2 SFCROSS) is not reachable from SIMC and throws simc::NotPorted.
/// @param Z target Z.
/// @param A target A.
/// @param qsq Q2 [GeV^2].
/// @param wsq W^2 [GeV^2].
/// @param[out] F1 structure function F1.
/// @param[out] F2 structure function F2.
/// @throws simc::NotPorted for int(A) >= 2.
void f1f2in21(double Z, double A, double qsq, double wsq, double& F1, double& F2);

/// @}

}  // namespace sidis
}  // namespace simc

#endif  // SIMC_PHYSICS_SIDIS_F1F2IN21_H
