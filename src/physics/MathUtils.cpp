// src/physics/MathUtils.cpp
// Implementation of mathematical special functions

#include "simc/MathUtils.h"
#include <cmath>
#include <stdexcept>

namespace simc {

// ============================================================================
// Spence Function (Dilogarithm)
// ============================================================================

double MathUtils::Spence(double x) {
    // Series expansion: Li₂(x) = Σ(n=1 to ∞) x^n/n²
    // Converges when |term/sum| < tolerance
    
    double sum = 0.0;
    double term = 1.0;  // x^n, starts at n=1
    
    for (int n = 1; n <= kMaxTerms; ++n) {
        term *= x;  // x^n
        double val = term / (n * n);
        sum += val;
        
        // Check convergence
        if (std::abs(sum) > 0.0 && std::abs(val) < std::abs(sum) * kTolerance) {
            break;
        }
    }
    
    return sum;
}

// ============================================================================
// Approximate Spence Function
// ============================================================================

double MathUtils::SpenceApprox(double x) {
    // Fast approximation used in brem.f when calculate_spence=false
    // This saves significant CPU time
    
    if (std::abs(x) <= 1.0) {
        return 0.0;
    } else {
        double lnx = std::log(std::abs(x));
        return -0.5 * lnx * lnx;
    }
}

} // namespace simc
