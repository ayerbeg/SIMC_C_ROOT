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

// Fortran: gauss1.f:gauss1
double Mt64Random::gauss1(double nsig_max) { return gauss1Polar(*this, nsig_max); }

}  // namespace simc

SIMC_REGISTER(simc::IRandom, "mt64", simc::Mt64Random)
