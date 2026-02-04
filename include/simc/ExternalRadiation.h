// include/simc/ExternalRadiation.h
// External bremsstrahlung corrections

#ifndef SIMC_EXTERNAL_RADIATION_H
#define SIMC_EXTERNAL_RADIATION_H

#include "simc/SimcConstants.h"
#include "simc/MathUtils.h"
#include <array>

namespace simc {

/**
 * @struct ExtRadResult
 * @brief Result of external radiation calculation
 */
struct ExtRadResult {
    double dbrem = 0.0;        ///< External correction δ_brem
    double dbrem_prime = 0.0;  ///< Derivative dδ_brem/dE
};

/**
 * @class ExternalRadiation
 * @brief External bremsstrahlung corrections
 * 
 * Handles external bremsstrahlung (radiation in target material)
 * using various approximations:
 * - Friedrich approximation
 * - Multiplicative φ corrections
 * - Dave's formulas
 * 
 * References:
 * - Friedrich & Lenz, Nucl. Phys. B 2, 441 (1968)
 * - Stein et al., Phys. Rev. D 12, 1884 (1975)
 */
class ExternalRadiation {
public:
    /**
     * @brief Friedrich approximation for external radiation
     * @param Ei Initial particle energy (MeV)
     * @param Ecutoff Cutoff energy (MeV)
     * @param trad Radiation thickness t_rad = bt/etatzai
     * @param etatzai Target thickness parameter η
     * @return ExtRadResult structure
     * 
     * Formula:
     * δ_brem = t_rad × [-(η-0.5) - η·ln(x) + η·x - 0.5·x²]
     * where x = Ecutoff/Ei
     */
    static ExtRadResult Friedrich(double Ei, double Ecutoff, double trad,
                                  double etatzai);
    
    /**
     * @brief Multiplicative φ correction factor
     * @param itail Which tail (0=both, 1=incoming e, 2=outgoing e)
     * @param E1 Energy 1 (MeV) - incoming electron
     * @param E2 Energy 2 (MeV) - outgoing electron
     * @param Egamma Photon energy (MeV)
     * @param extrad_flag External correction mode
     * @param bt Radiation thickness array [e, p]
     * @param g Exponent array g[0..4]
     * @param etatzai Target thickness parameter
     * @return φ correction factor
     * 
     * Three modes:
     * - extrad_flag = 1: φ = 1 (no correction)
     * - extrad_flag = 2: Dave's formulas
     * - extrad_flag = 3: Friedrich prescription
     */
    static double Phi(int itail, double E1, double E2, double Egamma,
                     int extrad_flag,
                     const std::array<double, 2>& bt,
                     const std::array<double, 5>& g,
                     double etatzai);
};

} // namespace simc

#endif // SIMC_EXTERNAL_RADIATION_H
