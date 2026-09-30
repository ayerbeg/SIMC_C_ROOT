// include/simc/core/random/Mt64Random.h
//
// Default IRandom: std::mt19937_64. D-002: the random sequence need not match the Fortran
// (grnd = RANLUX, call_ranlux.f); only the distributions must.
//   uniform(): 53-bit uniform strictly inside (0,1): ((k >> 11) + 0.5) * 2^-53. Like RANLUX it
//              never returns 0 (nor 1).
//   gauss1():  the gauss1.f algorithm (gauss1Polar) on this uniform stream.
// Registered as IRandom "mt64". JSON config (optional): { "seed": <unsigned integer> }.
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
    void init(RunSetup&) override {}

    double uniform() override { return (static_cast<double>(eng_() >> 11) + 0.5) * 0x1.0p-53; }
    double gauss1(double nsig_max) override;
    void seed(std::uint64_t s) override { eng_.seed(s); }

private:
    std::mt19937_64 eng_;
};

}  // namespace simc

#endif  // SIMC_CORE_RANDOM_MT64RANDOM_H
