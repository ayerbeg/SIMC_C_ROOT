// include/simc/radiation/kernel/Brem.h
// KERNEL header (not installed): included only by src/radiation, the radiation plug-in and the L1 tests.
//
// Exact soft-photon bremsstrahlung (brem.f), ported operation by operation (task P5.1).
// All internal kinematics are in GeV, as in the Fortran; dbsoft is returned per MeV.
//
//   brem         brem.f:6-214     on-shell ep, kinematics rebuilt from (ein, eout) assuming ELASTIC ep
//   inter        brem.f:218-247   interference integral
//   srad         brem.f:251-283   ultra-relativistic limit (diagnostic print only)
//   spence       brem.f:287-320   the approximation actually used (0 for |x|<=1, -ln^2|x|/2 otherwise)
//   bremos       brem.f:344-574   off-shell version (p_i%e and ami forced to mp)
//   interPrime   brem.f:581-596   derivative of inter w.r.t. de
// Not ported (never called in SIMC): ssum (brem.f:324-340), vac/del (brem.f:600-653).
//
// Numerics (D-010, P5.F1 Q9): every literal is a double; x**2 of a variable is x*x; x**0.5 is libm pow (NOT
// sqrt); PARAMETER powers such as ame**4 are folded by gfortran with correct rounding (kAme4, kAm4).
// produce_output prints to `log` (if non-null) in the Fortran order; the values are unchanged.
#ifndef SIMC_RADIATION_BREM_H
#define SIMC_RADIATION_BREM_H

#include <iosfwd>

#include "simc/radiation/kernel/RadcState.h"

namespace simc {
namespace radiation {

/// Fortran: brem.f:brem (l.6-214). Returns brem; bsoft, bhard, dbsoft (per MeV) are output arguments.
double brem(const BremFlags& f, double ein, double eout, double egamma, bool radiate_proton, double& bsoft,
            double& bhard, double& dbsoft, std::ostream* log = nullptr);

/// Fortran: brem.f:inter (l.218-247).
double inter(bool calculate_spence, double alpha, double ar1, double ar2, double e1, double e2, double de);

/// Values srad computes (it has no output arguments; they are only printed).
struct SradValues {
    double dsoft{0}, dsoftz{0}, dsoftzz{0};
};
/// Fortran: brem.f:srad (l.251-283). Prints if produce_output; returns the three computed values.
SradValues srad(double ak, double akp, double eang, double q2, double ame, double am, double ap, double de,
                bool produce_output, std::ostream* log = nullptr);

/// Fortran: brem.f:spence (l.287-320).
double spence(double ax);

/// Fortran: brem.f:bremos (l.344-574). Momenta in MeV/c, energies in MeV.
double bremos(const BremFlags& f, double egamma, double k_ix, double k_iy, double k_iz, double k_fx, double k_fy,
              double k_fz, double p_ix, double p_iy, double p_iz, double p_fx, double p_fy, double p_fz, double p_fe,
              bool radiate_proton, double& bsoft, double& bhard, double& dbsoft, std::ostream* log = nullptr);

/// Fortran: brem.f:inter_prime (l.581-596).
double interPrime(double alpha, double ar1, double ar2, double de);

}  // namespace radiation
}  // namespace simc

#endif  // SIMC_RADIATION_BREM_H
