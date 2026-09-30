// src/core/random/Gauss1.cpp — see include/simc/core/random/Gauss1.h.
#include "simc/core/random/Gauss1.h"

#include <cmath>

#include "simc/stages/IRandom.h"

namespace simc {

// Fortran: gauss1.f:gauss1 (l.1-30)
double gauss1Polar(IRandom& rng, double nsigmax) {
    for (;;) {
        const double u1 = rng.uniform();   // 1  u1 = grnd()
        const double u2 = rng.uniform();   //    u2 = grnd()
        const double v1 = 2.0 * u1 - 1.0;
        const double v2 = 2.0 * u2 - 1.0;
        const double s = v1 * v1 + v2 * v2;          // v1**2+v2**2
        if (s > 1. || s == 0) continue;              // goto 1
        const double g = v1 * std::sqrt(-2. * std::log(s) / s);
        if (std::fabs(g) > nsigmax) continue;        // truncate at nsigmax (goto 1)
        return g;
    }
}

}  // namespace simc
