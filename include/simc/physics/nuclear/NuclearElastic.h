/**
 * @file
 * @ingroup physics
 * @brief A(e,e')A nuclear elastic cross section (upstream SIMC nuclear_elastic.f, commit 4df7351).
 *
 * sig_nuc_elastic = (W2 + 2 W1 tan^2(theta/2)) * sigMott(E', theta, Q2) * M/(M + Ein (1 - cos theta)), with W1 = 0
 * and W2 = (Z Fgauss(Q2, A))^2 (spin-0 nucleus, Gaussian charge form factor; finding U-B2).
 * The IPhysicsModel "nuc_elastic" (event.f:complete_main at 20d46fc, l.1457-1459) calls it for vertex and recon.
 * @note Origin: upstream SIMC nuclear_elastic.f (sig_nuc_elastic l.1-27, NUC_FORM_FACTOR l.30-45, Fgauss l.48-66).
 */
#ifndef SIMC_PHYSICS_NUCLEAR_NUCLEARELASTIC_H
#define SIMC_PHYSICS_NUCLEAR_NUCLEARELASTIC_H

namespace simc {
struct Event;
namespace nuclear {

/// @addtogroup physics
/// @{

/// Fortran: nuclear_elastic.f:Fgauss (l.48-66). Gaussian nuclear charge form factor.
/// @param T Q2 [MeV^2].
/// @param A mass number.
/// @return the form factor (0 when (T/hbarc^2) 2.4^2/6 >= 80).
double Fgauss(double T, double A);

/// Fortran: nuclear_elastic.f:NUC_FORM_FACTOR (l.30-45): W1 = 0, FF = Fgauss(QSQ, A), W2 = (Z FF)^2.
/// @param QSQ Q2 [MeV^2].
/// @param W1 receives the transverse structure function (always 0).
/// @param W2 receives the longitudinal structure function (Z FF)^2.
/// @param FF receives the form factor Fgauss(QSQ, A).
/// @param Z target charge.
/// @param A mass number.
void NucFormFactor(double QSQ, double& W1, double& W2, double& FF, double Z, double A);

/// Fortran: nuclear_elastic.f:sig_nuc_elastic (l.1-27).
/// @param Ein beam energy at the vertex [MeV].
/// @param Eprime scattered electron energy vertex%e%E [MeV].
/// @param theta electron angle [rad].
/// @param Q2 four-momentum transfer squared [MeV^2].
/// @param Z target charge.
/// @param A mass number.
/// @param M target mass [MeV].
/// @return the cross section [microbarn/sr].
double SigNucElastic(double Ein, double Eprime, double theta, double Q2, double Z, double A, double M);

/// sig_nuc_elastic(ev) with ev%Ein, ev%e%E, ev%e%theta, ev%Q2 (as called from event.f:complete_main).
/// @param ev the vertex or reconstructed event.
/// @param Z target charge.
/// @param A mass number.
/// @param M target mass [MeV].
/// @return the cross section [microbarn/sr].
double SigNucElastic(const Event& ev, double Z, double A, double M);

/// @}

}  // namespace nuclear
}  // namespace simc

#endif  // SIMC_PHYSICS_NUCLEAR_NUCLEARELASTIC_H
