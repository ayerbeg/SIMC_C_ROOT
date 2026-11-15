/**
 * @file
 * @ingroup core
 * @brief simc::gauss1Polar: the truncated-Gaussian algorithm of SIMC on top of any uniform source.
 *
 * Polar Box-Muller: only v1 of each accepted pair is used; a pair is rejected (and a NEW pair drawn) if s > 1 or
 * s == 0, and again if |g| > nsigmax. FORTRAN-QUIRK(D-012 F6): s == 1 exactly is accepted and gives g = +-0.
 * FORTRAN-QUIRK(D-012 F2): u1 == 0.5 exactly gives g = 0 (RANLUX 2^-24 granularity; with the 53-bit simc::Mt64Random
 * stream this practically never happens).
 * @note Origin: Fortran SIMC gauss1.f:gauss1 (l.1-30).
 */
#ifndef SIMC_CORE_RANDOM_GAUSS1_H
#define SIMC_CORE_RANDOM_GAUSS1_H

namespace simc {

/** @addtogroup core
 *  @{ */

class IRandom;

/**
 * @brief A standard normal deviate truncated at |g| <= nsigmax, by the polar method of the original program.
 *
 * Consumes 2 uniforms per attempt from rng.uniform(). If the uniform source returns NaN (an exhausted
 * simc::ScriptedRandom queue) the result is NaN, exactly as the reference test harness behaves (no comparison is true
 * for NaN).
 * @param rng The uniform source (only IRandom::uniform() is called).
 * @param nsigmax The truncation in standard deviations.
 * @return The deviate.
 * @note Origin: Fortran SIMC gauss1.f:gauss1(nsigmax).
 */
double gauss1Polar(IRandom& rng, double nsigmax);

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_RANDOM_GAUSS1_H
