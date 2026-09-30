// include/simc/core/random/Gauss1.h
//
// The gauss1.f algorithm on top of any uniform source (IRandom::uniform()).
// Fortran: gauss1.f:gauss1 (l.1-30). Polar Box-Muller; only v1 is used; a pair is rejected
// (and a NEW pair drawn) if s > 1 or s == 0, and again if |g| > nsigmax.
// FORTRAN-QUIRK(D-012 F6): s == 1 exactly is accepted and gives g = +-0.
// FORTRAN-QUIRK(D-012 F2): u1 == 0.5 exactly gives g = 0 (RANLUX 2^-24 granularity; with the
// 53-bit Mt64Random stream this practically never happens).
#ifndef SIMC_CORE_RANDOM_GAUSS1_H
#define SIMC_CORE_RANDOM_GAUSS1_H

namespace simc {

class IRandom;

/// Fortran: gauss1.f:gauss1(nsigmax). Consumes 2 uniforms per attempt from rng.uniform().
/// If the uniform source returns NaN (exhausted ScriptedRandom queue) the result is NaN,
/// exactly as the Fortran harness stub behaves (no comparison is true for NaN).
double gauss1Polar(IRandom& rng, double nsigmax);

}  // namespace simc

#endif  // SIMC_CORE_RANDOM_GAUSS1_H
