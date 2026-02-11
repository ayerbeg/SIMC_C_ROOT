// src/physics/ExternalRadiation.cpp
// Implementation of external bremsstrahlung corrections

#include "simc/physics/ExternalRadiation.h"
#include <cmath>
#include <stdexcept>

namespace simc {

// ============================================================================
// Friedrich Approximation
// ============================================================================

ExtRadResult ExternalRadiation::Friedrich(double Ei, double Ecutoff, 
                                          double trad, double etatzai) {
    ExtRadResult result;
    
    // x = Ecutoff / Ei
    double x = Ecutoff / Ei;
    
    // δ_brem = t_rad × [-(η-0.5) - η·ln(x) + η·x - 0.5·x²]
    result.dbrem = trad * (-(etatzai - 0.5) - etatzai * std::log(x) 
                          + etatzai * x - 0.5 * x * x);
    
    // dδ_brem/dE = -t_rad/Ei × [η/x - η + x]
    result.dbrem_prime = -trad / Ei * (etatzai / x - etatzai + x);
    
    return result;
}

// ============================================================================
// Multiplicative Phi Correction
// ============================================================================

double ExternalRadiation::Phi(int itail, double E1, double E2, double Egamma,
                              int extrad_flag,
                              const std::array<double, 2>& bt,
                              const std::array<double, 5>& g,
                              double etatzai) {
    // Store energies in array for easier indexing
    std::array<double, 2> E = {E1, E2};
    
    double phi = 1.0;  // Default: no correction
    
    if (extrad_flag == 1) {
        // Mode 1: No correction
        phi = 1.0;
        
    } else if (extrad_flag == 2) {
        // Mode 2: Dave's formulas
        
        if (itail == 0) {
            // Both tails
            phi = 1.0 - (bt[0]/E[0] + bt[1]/E[1]) / (g[1] + g[2]) * Egamma;
            
        } else if (itail == 1 || itail == 2) {
            // Single tail (1=incoming e, 2=outgoing e)
            int idx = itail - 1;  // Convert to array index (0 or 1)
            phi = 1.0 - bt[idx]/E[idx] / g[itail] * Egamma;
            
        } else {
            throw std::invalid_argument("Phi: Invalid itail for extrad_flag=2");
        }
        
    } else if (extrad_flag == 3) {
        // Mode 3: Friedrich prescription
        
        if (itail == 0) {
            throw std::invalid_argument(
                "Phi not defined for peaking approx (itail=0) with extrad_flag=3");
        }
        
        if (itail == 1 || itail == 2) {
            int idx = itail - 1;
            double x = Egamma / E[idx];
            double t = bt[idx] / etatzai;
            
            // φ = (1 - x + x²/η) × exp(t×[(η-0.5) - η·x + x²/2]) × Γ(1+bt)
            phi = (1.0 - x + x*x/etatzai) 
                * std::exp(t * ((etatzai - 0.5) - etatzai * x + x*x/2.0))
                * MathUtils::Gamma(1.0 + bt[idx]);
                
        } else {
            throw std::invalid_argument("Phi: Invalid itail for extrad_flag=3");
        }
        
    } else {
        throw std::invalid_argument("Phi: Invalid extrad_flag");
    }
    
    return phi;
}

} // namespace simc
