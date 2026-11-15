/**
 * @file
 * @ingroup core
 * @brief simc::Mt64Random: the default random source ("mt64"), std::mt19937_64.
 *
 * Decision D-002: the random sequence need not match the original program (RANLUX), only the distributions must.
 *  - uniform(): 53-bit uniform strictly inside (0,1): ((k >> 11) + 0.5) * 2^-53. Like RANLUX it never returns 0 (nor 1).
 *  - gauss1():  the gauss1.f algorithm (simc::gauss1Polar) on this uniform stream.
 *
 * Registered as IRandom "mt64". JSON configuration (optional): `{ "seed": <unsigned integer> }`.
 *
 * Parallel event loop (`--threads N > 1`, task P15.4): worker k gets its own engine, seeded by seedStream(base, k) =
 * std::seed_seq{base & 0xffffffff, base >> 32, k, 0x53494d43 ("SIMC")}, where base = BaseSeed(config, params) is the
 * seed the single-thread run would use. A stream depends only on (base, k) (not on N); every stream differs from the
 * single-thread one (seed(base) initialises the engine with the mt19937_64 integer recursion, not with a seed_seq).
 */
#ifndef SIMC_CORE_RANDOM_MT64RANDOM_H
#define SIMC_CORE_RANDOM_MT64RANDOM_H

#include <cstdint>
#include <random>
#include <string>

#include "simc/stages/IRandom.h"

namespace simc {

/** @addtogroup core
 *  @{ */

/**
 * @brief The default random source: std::mt19937_64 with a 53-bit uniform in (0, 1) (see the file description).
 *
 * Not thread-safe; the parallel loop gives every worker its own instance (seedStream()).
 */
class Mt64Random : public IRandom {
public:
    /// The seed used when neither the configuration nor the input gives one (the std::mt19937_64 default).
    static constexpr std::uint64_t kDefaultSeed = 5489u;   // std::mt19937_64 default

    /**
     * @brief Creates a source seeded with @p s.
     * @param s The seed (configure() may re-seed).
     */
    explicit Mt64Random(std::uint64_t s = kDefaultSeed) : eng_(s) {}

    std::string name() const override { return "mt64"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext&) override {}

    double uniform() override { return (static_cast<double>(eng_() >> 11) + 0.5) * 0x1.0p-53; }
    double gauss1(double nsig_max) override;
    void seed(std::uint64_t s) override { eng_.seed(s); }

    /**
     * @brief The seed configure() applies.
     * @param config The stage configuration (key "seed", optional).
     * @param params The input parameters (random_seed).
     * @return config "seed" if present, else params.random_seed if > 0, else #kDefaultSeed.
     */
    static std::uint64_t BaseSeed(const Json& config, const SimcParameters& params);
    /**
     * @brief Seeds the independent stream of one worker of the parallel event loop (see the file description).
     * @param base The single-thread seed (BaseSeed()).
     * @param stream The worker index k.
     */
    void seedStream(std::uint64_t base, std::uint64_t stream);

private:
    std::mt19937_64 eng_;
};

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_RANDOM_MT64RANDOM_H
