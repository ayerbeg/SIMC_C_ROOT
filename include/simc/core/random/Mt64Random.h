// include/simc/core/random/Mt64Random.h
//
// Default IRandom: std::mt19937_64. D-002: the random sequence need not match the Fortran
// (grnd = RANLUX, call_ranlux.f); only the distributions must.
//   uniform(): 53-bit uniform strictly inside (0,1): ((k >> 11) + 0.5) * 2^-53. Like RANLUX it
//              never returns 0 (nor 1).
//   gauss1():  the gauss1.f algorithm (gauss1Polar) on this uniform stream.
// Registered as IRandom "mt64". JSON config (optional): { "seed": <unsigned integer> }.
// Parallel event loop (task P15.4, --threads N > 1): worker k gets its own engine, seeded by seedStream(base, k) =
// std::seed_seq{base & 0xffffffff, base >> 32, k, 0x53494d43 ("SIMC")}, where base = BaseSeed(config, params) is the seed
// the single-thread run would use. A stream depends only on (base, k) (not on N); every stream differs from the
// single-thread one (seed(base) initialises the engine with the mt19937_64 integer recursion, not with a seed_seq).
#ifndef SIMC_CORE_RANDOM_MT64RANDOM_H
#define SIMC_CORE_RANDOM_MT64RANDOM_H

#include <cstdint>
#include <random>
#include <string>

#include "simc/stages/IRandom.h"

namespace simc {

class Mt64Random : public IRandom {
public:
    static constexpr std::uint64_t kDefaultSeed = 5489u;   // std::mt19937_64 default

    explicit Mt64Random(std::uint64_t s = kDefaultSeed) : eng_(s) {}

    std::string name() const override { return "mt64"; }
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext&) override {}

    double uniform() override { return (static_cast<double>(eng_() >> 11) + 0.5) * 0x1.0p-53; }
    double gauss1(double nsig_max) override;
    void seed(std::uint64_t s) override { eng_.seed(s); }

    /// The seed configure() applies: config "seed", else params.random_seed if > 0, else kDefaultSeed.
    static std::uint64_t BaseSeed(const Json& config, const SimcParameters& params);
    /// Independent stream `stream` of the parallel event loop (see the header comment).
    void seedStream(std::uint64_t base, std::uint64_t stream);

private:
    std::mt19937_64 eng_;
};

}  // namespace simc

#endif  // SIMC_CORE_RANDOM_MT64RANDOM_H
