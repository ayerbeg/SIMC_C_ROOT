// include/simc/stages/IRandom.h
// Random-number source. D-002: the sequence need not match the Fortran, the
// distributions must (Fortran: grnd = mt19937.f/RANLUX, gauss1.f).
#ifndef SIMC_STAGES_IRANDOM_H
#define SIMC_STAGES_IRANDOM_H

#include <cstdint>

#include "simc/stages/IStage.h"

namespace simc {

class IRandom : public IStage {
public:
    /// Uniform in [0,1). Fortran: grnd().
    virtual double uniform() = 0;
    /// Standard normal truncated at |x| <= nsig_max by rejection. Fortran: gauss1.f:gauss1(nsigmax).
    virtual double gauss1(double nsig_max) = 0;
    /// Re-seed. Fortran: sgrnd(random_seed) / RLUXGO.
    virtual void seed(std::uint64_t s) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_IRANDOM_H
