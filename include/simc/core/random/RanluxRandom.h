/**
 * @file
 * @ingroup core
 * @brief simc::RanluxRandom ("ranlux"): the random-number generator of the original SIMC, ported bit for bit.
 *
 * With it a SIMC++ run can reproduce the random stream of an original run exactly (bitwise replays of whole runs;
 * decision D-002 keeps "mt64" the default). Task P13.2.
 *  - RANLUX: Luescher's luxury subtract-with-borrow generator (F. James 1993), RLUXGO (initialisation from one
 *    integer), RLUXIN / RLUXUT (the 25-integer state); every value is a multiple of 2^-48, so the arithmetic is exact in
 *    double precision.
 *  - grnd: uniform() hands out RVEC(1..999) of each block of 1000 numbers and refills when latest reaches 1000 (the
 *    1000th number of every block is never used); seed(s) = RLUXGO(3, s, 0, 0).
 *  - gauss1(): simc::gauss1Polar on this uniform stream.
 *
 * JSON configuration (optional): `{ "seed": <integer> }` (default: the `.inp` random_seed), or
 * `{ "state": [25 integers] }` (a saved RANLUX state, the `*_start_random_state.dat` content).
 * @note Origin: Fortran SIMC cern/ranlux.f (compiled with -fdefault-real-8), call_ranlux.f (sgrnd, grnd), gauss1.f;
 *       simc.f:107-109.
 */
#ifndef SIMC_CORE_RANDOM_RANLUXRANDOM_H
#define SIMC_CORE_RANDOM_RANLUXRANDOM_H

#include <array>
#include <cstdint>
#include <string>

#include "simc/stages/IRandom.h"

namespace simc {

/** @addtogroup core
 *  @{ */

/**
 * @brief Bitwise port of the RANLUX generator and of the grnd/sgrnd wrappers of the original program.
 *
 * Holds the full generator state (no globals); not thread-safe.
 */
class RanluxRandom : public IRandom {
public:
    /// Creates an unseeded generator (as a fresh process: the first call initialises it like RANLUX does).
    RanluxRandom();

    std::string name() const override { return "ranlux"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext&) override {}

    /**
     * @brief The next number of the grnd stream.
     * @return A uniform number in (0, 1), a multiple of 2^-48.
     * @note Origin: Fortran SIMC call_ranlux.f:grnd.
     */
    double uniform() override;
    /**
     * @brief gauss1 on the grnd stream (simc::gauss1Polar).
     * @param nsig_max The truncation in standard deviations.
     * @return The truncated Gaussian deviate.
     * @note Origin: Fortran SIMC gauss1.f.
     */
    double gauss1(double nsig_max) override;
    /**
     * @brief Seeds as sgrnd(seed) = RLUXGO(3, seed, 0, 0) and resets grnd's block pointer as a fresh process.
     * @param s The seed (used as a 32-bit integer, as the original).
     * @note Origin: Fortran SIMC call_ranlux.f:sgrnd.
     */
    void seed(std::uint64_t s) override;

    /**
     * @brief Initialises the generator (RLUXGO).
     * @param lux Luxury level 0..4 (SIMC uses 3).
     * @param ins The seed integer.
     * @param k1 Number of numbers already generated, modulo 10^9 (skip-ahead), normally 0.
     * @param k2 Number of numbers already generated, multiple of 10^9, normally 0.
     * @note Origin: Fortran SIMC ranlux.f:RLUXGO(lux, ins, k1, k2).
     */
    void rluxgo(int lux, int ins, int k1, int k2);
    /**
     * @brief Restores a saved state (RLUXIN).
     * @param isdext The 25 integers written by rluxut().
     * @note Origin: Fortran SIMC ranlux.f:RLUXIN(isdext).
     */
    void rluxin(const std::array<int, 25>& isdext);
    /**
     * @brief Saves the state (RLUXUT).
     * @return The 25 integers that rluxin() accepts.
     * @note Origin: Fortran SIMC ranlux.f:RLUXUT(isdext).
     */
    std::array<int, 25> rluxut() const;
    /**
     * @brief Generates a vector of numbers (RANLUX).
     * @param rvec Output array of at least @p lenv numbers in (0, 1).
     * @param lenv Number of values to generate.
     * @note Origin: Fortran SIMC ranlux.f:RANLUX(rvec, lenv).
     */
    void ranlux(double* rvec, int lenv);

private:
    double seeds_[25]{};   // SEEDS(1..24)
    int next_[25]{};       // NEXT(1..24)
    int i24_{24}, j24_{10}, in24_{0}, kount_{0}, mkount_{0}, luxlev_{3}, nskip_{199}, inseed_{0};
    double carry_{0}, twom24_{0}, twom12_{0};
    bool notyet_{true};
    // call_ranlux.f:grnd SAVEd locals (static zero under -fno-automatic)
    double rvec_[1001]{};
    int latest_{0};
};

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_RANDOM_RANLUXRANDOM_H
