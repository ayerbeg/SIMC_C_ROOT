// src/core/random/Mt64Random.cpp — see include/simc/core/random/Mt64Random.h.
#include "simc/core/random/Mt64Random.h"

#include "simc/core/Registry.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/Gauss1.h"

namespace simc {

void Mt64Random::configure(const Json& config, const SimcParameters& params) {
    // Fortran: simc.f:107 sgrnd(random_seed). The Fortran random_seed = -1 (time) is the driver's job.
    if (config.is_object() && config.contains("seed")) {
        seed(config.at("seed").get<std::uint64_t>());
    } else if (params.random_seed > 0) {
        seed(static_cast<std::uint64_t>(params.random_seed));
    }
}

std::uint64_t Mt64Random::BaseSeed(const Json& config, const SimcParameters& params) {
    if (config.is_object() && config.contains("seed")) return config.at("seed").get<std::uint64_t>();
    if (params.random_seed > 0) return static_cast<std::uint64_t>(params.random_seed);
    return kDefaultSeed;
}

// P15.4: worker streams of the parallel event loop.
void Mt64Random::seedStream(std::uint64_t base, std::uint64_t stream) {
    std::seed_seq seq{static_cast<std::uint32_t>(base & 0xffffffffu), static_cast<std::uint32_t>(base >> 32),
                      static_cast<std::uint32_t>(stream & 0xffffffffu), static_cast<std::uint32_t>(0x53494d43u)};
    eng_.seed(seq);
}

// Fortran: gauss1.f:gauss1
double Mt64Random::gauss1(double nsig_max) { return gauss1Polar(*this, nsig_max); }

}  // namespace simc

SIMC_REGISTER(simc::IRandom, "mt64", simc::Mt64Random)
