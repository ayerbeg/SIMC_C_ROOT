/**
 * @file
 * @ingroup radiation
 * @brief Exact soft-photon bremsstrahlung (brem.f), ported operation by operation (kernel header).
 *
 * KERNEL header (not installed): included only by src/radiation, the radiation plug-in and the L1 tests.
 *
 * Exact soft-photon bremsstrahlung (brem.f), ported operation by operation (task P5.1).
 * All internal kinematics are in GeV, as in the Fortran; dbsoft is returned per MeV.
 *
 *   brem         brem.f:6-214     on-shell ep, kinematics rebuilt from (ein, eout) assuming ELASTIC ep
 *   inter        brem.f:218-247   interference integral
 *   srad         brem.f:251-283   ultra-relativistic limit (diagnostic print only)
 *   spence       brem.f:287-320   the approximation actually used (0 for |x|<=1, -ln^2|x|/2 otherwise)
 *   bremos       brem.f:344-574   off-shell version (p_i%e and ami forced to mp)
 *   interPrime   brem.f:581-596   derivative of inter w.r.t. de
 * Not ported (never called in SIMC): ssum (brem.f:324-340), vac/del (brem.f:600-653).
 *
 * Numerics (D-010, P5.F1 Q9): every literal is a double; x**2 of a variable is x*x; x**0.5 is libm pow (NOT
 * sqrt); PARAMETER powers such as ame**4 are folded by gfortran with correct rounding (kAme4, kAm4).
 * produce_output prints to `log` (if non-null) in the Fortran order; the values are unchanged.
 */
#ifndef SIMC_RADIATION_BREM_H
#define SIMC_RADIATION_BREM_H

#include <iosfwd>

#include "simc/radiation/kernel/RadcState.h"

namespace simc {
namespace radiation {

/// @addtogroup radiation
/// @{

/// Fortran: brem.f:brem (l.6-214). Returns brem; bsoft, bhard, dbsoft (per MeV) are output arguments.
/// @param f the /bremcom/ flags.
/// @param ein incoming electron energy [MeV].
/// @param eout scattered electron energy [MeV].
/// @param egamma photon energy [MeV].
/// @param radiate_proton include the proton radiation.
/// @param[out] bsoft soft-photon correction.
/// @param[out] bhard hard-photon correction.
/// @param[out] dbsoft derivative of bsoft [1/MeV].
/// @param log stream for produce_output (not owned), or nullptr.
/// @return the correction factor brem.
double brem(const BremFlags& f, double ein, double eout, double egamma, bool radiate_proton, double& bsoft,
            double& bhard, double& dbsoft, std::ostream* log = nullptr);

/// Fortran: brem.f:inter (l.218-247).
/// @param calculate_spence include the Spence-function terms.
/// @param alpha,ar1,ar2 integral parameters (brem.f notation).
/// @param e1,e2 energies of the two particles [GeV].
/// @param de photon energy cut [GeV].
/// @return the interference integral.
double inter(bool calculate_spence, double alpha, double ar1, double ar2, double e1, double e2, double de);

/// Values srad computes (it has no output arguments; they are only printed).
struct SradValues {
    double dsoft{0} /**< @copydoc dsoftzz */, dsoftz{0} /**< @copydoc dsoftzz */, dsoftzz{0};  ///< the three soft corrections of the ultra-relativistic limit.
};
/// Fortran: brem.f:srad (l.251-283). Prints if produce_output; returns the three computed values.
/// @param ak,akp incoming and scattered electron momenta [GeV].
/// @param eang electron scattering angle [rad].
/// @param q2 four-momentum transfer squared [GeV^2].
/// @param ame electron mass [GeV].
/// @param am proton mass [GeV].
/// @param ap proton momentum [GeV].
/// @param de photon energy cut [GeV].
/// @param produce_output print the values.
/// @param log stream for the print (not owned), or nullptr.
/// @return the computed values.
SradValues srad(double ak, double akp, double eang, double q2, double ame, double am, double ap, double de,
                bool produce_output, std::ostream* log = nullptr);

/// Fortran: brem.f:spence (l.287-320).
/// @param ax the argument.
/// @return the approximated Spence function.
double spence(double ax);

/// Fortran: brem.f:bremos (l.344-574). Momenta in MeV/c, energies in MeV.
/// @param f the /bremcom/ flags.
/// @param egamma photon energy [MeV].
/// @param k_ix,k_iy,k_iz incoming electron 3-momentum [MeV/c].
/// @param k_fx,k_fy,k_fz scattered electron 3-momentum [MeV/c].
/// @param p_ix,p_iy,p_iz incoming proton 3-momentum (pm) [MeV/c].
/// @param p_fx,p_fy,p_fz scattered proton 3-momentum [MeV/c].
/// @param p_fe scattered proton energy [MeV].
/// @param radiate_proton include the proton radiation.
/// @param[out] bsoft soft-photon correction.
/// @param[out] bhard hard-photon correction.
/// @param[out] dbsoft derivative of bsoft [1/MeV].
/// @param log stream for produce_output (not owned), or nullptr.
/// @return the correction factor.
double bremos(const BremFlags& f, double egamma, double k_ix, double k_iy, double k_iz, double k_fx, double k_fy,
              double k_fz, double p_ix, double p_iy, double p_iz, double p_fx, double p_fy, double p_fz, double p_fe,
              bool radiate_proton, double& bsoft, double& bhard, double& dbsoft, std::ostream* log = nullptr);

/// Fortran: brem.f:inter_prime (l.581-596).
/// @param alpha,ar1,ar2 integral parameters (brem.f notation).
/// @param de photon energy cut [GeV].
/// @return the derivative of inter() with respect to de.
double interPrime(double alpha, double ar1, double ar2, double de);

/// @}

}  // namespace radiation
}  // namespace simc

#endif  // SIMC_RADIATION_BREM_H
