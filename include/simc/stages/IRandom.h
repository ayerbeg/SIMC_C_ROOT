/**
 * @file
 * @ingroup stages
 * @brief simc::IRandom: the random-number source of a run.
 *
 * Decision D-002: the random sequence need not match the original program, the distributions must. The default source
 * is simc::Mt64Random ("mt64"); simc::RanluxRandom ("ranlux") reproduces the original RANLUX stream bit for bit;
 * simc::ScriptedRandom replays fixed numbers for tests.
 * @note Origin: Fortran SIMC grnd (call_ranlux.f / RANLUX) and gauss1.f.
 */
#ifndef SIMC_STAGES_IRANDOM_H
#define SIMC_STAGES_IRANDOM_H

#include <cstdint>

#include "simc/stages/IStage.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief Source of uniform and truncated-Gaussian random numbers.
 *
 * Every random number of a run comes from one IRandom (one per worker in the parallel loop). The number and order of
 * calls is part of the defined behaviour of each kernel. Not thread-safe: one object per thread.
 */
class IRandom : public IStage {
public:
    /**
     * @brief A uniform random number.
     * @return A number in [0, 1) (the built-in sources never return 0 either).
     * @note Origin: Fortran SIMC grnd().
     */
    virtual double uniform() = 0;
    /**
     * @brief A standard normal random number truncated at |x| <= nsig_max by rejection.
     * @param nsig_max The truncation in standard deviations (> 0).
     * @return The Gaussian deviate (mean 0, sigma 1 before truncation).
     * @note Origin: Fortran SIMC gauss1.f:gauss1(nsigmax).
     */
    virtual double gauss1(double nsig_max) = 0;
    /**
     * @brief Re-seeds the source.
     * @param s The seed; the same seed gives the same sequence.
     * @note Origin: Fortran SIMC sgrnd(random_seed) / RLUXGO.
     */
    virtual void seed(std::uint64_t s) = 0;
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_IRANDOM_H
