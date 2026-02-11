// include/simc/MathUtils.h
// Mathematical special functions for physics calculations

#ifndef SIMC_MATH_UTILS_H
#define SIMC_MATH_UTILS_H

#include <cmath>

namespace simc {

/**
 * @class MathUtils
 * @brief Mathematical special functions for radiative corrections
 * 
 * Provides mathematical functions needed for QED radiative correction
 * calculations, including the Spence function (dilogarithm).
 */
class MathUtils {
public:
    /**
     * @brief Gamma function Γ(x)
     * @param x Argument
     * @return Γ(x)
     * 
     * Uses C++11 std::tgamma() from <cmath>
     */
    static double Gamma(double x) {
        return std::tgamma(x);
    }
    
    /**
     * @brief Spence function (dilogarithm) Li₂(x)
     * @param x Argument
     * @return Li₂(x) = Σ(n=1 to ∞) x^n/n²
     * 
     * Series expansion converges when |term/sum| < tolerance.
     * Maximum 100 terms.
     * 
     * Used in exact bremsstrahlung calculations and Schwinger corrections.
     */
    static double Spence(double x);
    
    /**
     * @brief Approximate Spence function (fast version)
     * @param x Argument
     * @return Approximation: 0 if |x|≤1, -0.5·ln²(|x|) otherwise
     * 
     * This is the approximation used in brem.f when calculate_spence=false.
     * Saves significant CPU time at the cost of some accuracy.
     */
    static double SpenceApprox(double x);

private:
    static constexpr int kMaxTerms = 100;          ///< Maximum series terms
    static constexpr double kTolerance = 1.0e-4;   ///< Convergence criterion
};

} // namespace simc

#endif // SIMC_MATH_UTILS_H
