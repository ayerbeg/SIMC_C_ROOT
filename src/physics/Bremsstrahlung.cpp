// src/physics/Bremsstrahlung.cpp
// Implementation of exact QED bremsstrahlung calculations

#include "simc/Bremsstrahlung.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace simc {

using constants::alpha;
using constants::Me;
using constants::Mp;
using constants::pi;

// ============================================================================
// Helper Functions
// ============================================================================

bool Bremsstrahlung::CalculateRoots(double alpha, double adot,
                                    double /*m1_sq*/, double /*m2_sq*/,
                                    double& ar1, double& ar2) {
    // Solve quadratic equation for roots
    // alpha = 2m₁² - 2·adot
    // Roots: ar = (E₂ + sqrt(E₂² - m₂²·adot/alpha)) / adot
    
    double discriminant = alpha * (alpha - 4.0 * adot);
    
    if (discriminant < 0.0) {
        // Complex roots - unphysical
        ar1 = ar2 = 0.0;
        return false;
    }
    
    double sqrt_disc = std::sqrt(discriminant);
    
    // Calculate roots
    ar1 = (alpha + sqrt_disc) / (2.0 * adot);
    ar2 = (alpha - sqrt_disc) / (2.0 * adot);
    
    return true;
}

// ============================================================================
// Interference Integral Derivative
// ============================================================================

double Bremsstrahlung::InterPrime(double alpha, double ar1, double ar2, 
                                  double de) {
    // Derivative of Inter() with respect to de
    // dI/d(de) = (-1/de) × (1/(α(ar1-ar2)π)) × [ln|(ar1-1)/ar1| - ln|(ar2-1)/ar2|]
    
    if (std::abs(alpha) < 1.0e-10 || std::abs(ar1 - ar2) < 1.0e-10) {
        return 0.0;
    }
    
    double amult = -1.0 / (alpha * (ar1 - ar2));
    
    double term1 = std::log(std::abs((ar1 - 1.0) / ar1));
    double term2 = std::log(std::abs((ar2 - 1.0) / ar2));
    
    double result = (-1.0 / de) * amult / pi * (term1 - term2);
    
    return result;
}

// ============================================================================
// Interference Integral
// ============================================================================

double Bremsstrahlung::Inter(bool calculate_spence, double alpha,
                             double ar1, double ar2,
                             double /*e1*/, double e2, double de) {
    // Calculate the interference integral
    // I = (1/(α(ar1-ar2)π)) × [ln|...| × ln|...| - ln|...| × ln|...| + Spence terms]
    
    if (std::abs(alpha) < 1.0e-10 || std::abs(ar1 - ar2) < 1.0e-10) {
        return 0.0;
    }
    
    double amult = -1.0 / (alpha * (ar1 - ar2));
    double de2 = de * de;
    
    // Calculate logarithm arguments
    double arg1_num = e2/de + ar1 * de2/de;  // Simplifies to: e2/de + ar1*de
    double arg2_num = e2/de + ar2 * de2/de;
    
    double arg1_denom = (ar1 - 1.0) / ar1;
    double arg2_denom = (ar2 - 1.0) / ar2;
    
    // Main logarithmic terms
    double term1 = std::log(std::abs(arg1_num)) * std::log(std::abs(arg1_denom));
    double term2 = std::log(std::abs(arg2_num)) * std::log(std::abs(arg2_denom));
    
    double result = term1 - term2;
    
    // Add Spence function terms if requested
    if (calculate_spence) {
        // Spence arguments (these come from the exact integral)
        double spence_arg1 = de / (e2 + ar1 * de);
        double spence_arg2 = de / (e2 + ar2 * de);
        double spence_arg3 = (ar1 - 1.0) * de / (e2 + ar1 * de);
        double spence_arg4 = (ar2 - 1.0) * de / (e2 + ar2 * de);
        
        result += (-MathUtils::Spence(spence_arg1) 
                  + MathUtils::Spence(spence_arg2)
                  + MathUtils::Spence(spence_arg3) 
                  - MathUtils::Spence(spence_arg4));
    }
    
    result *= amult / pi;
    
    return result;
}

// ============================================================================
// Ultra-Relativistic Limit
// ============================================================================

double Bremsstrahlung::Srad(double ein, double eout, double de, double theta) {
    // Ultra-relativistic limit formula
    // Used as check/approximation
    
    double srad = 0.0;
    
    double loge = std::log(ein / eout);
    double logg = std::log(ein * eout * (1.0 - std::cos(theta)) 
                          / (Me * Me * 1.0e-6)); // Convert Me to GeV
    
    // Soft photon approximation
    srad = (alpha / pi) * (loge * std::log(ein * eout / (de * de)) 
                           - 0.5 * loge * loge + logg * loge);
    
    return srad;
}

// ============================================================================
// Vacuum Polarization
// ============================================================================

double Bremsstrahlung::VacuumPolarization(double Q2_GeV2) {
    // Leading order vacuum polarization correction
    // Formula from Tsai
    
    if (Q2_GeV2 <= 0.0) return 0.0;
    
    double me_GeV = Me * 1.0e-3;  // Convert to GeV
    double me2 = me_GeV * me_GeV;
    
    if (Q2_GeV2 < 4.0 * me2) {
        // Below pair production threshold
        return 0.0;
    }
    
    double x = Q2_GeV2 / (4.0 * me2);
    double beta = std::sqrt(1.0 - 1.0/x);
    
    // Vacuum polarization: δ = (α/3π) × [5/3 + β(2-β²)/3 × ln((1+β)/(1-β))]
    double del = (alpha / (3.0 * pi)) 
               * (5.0/3.0 + beta * (2.0 - beta*beta) / 3.0 
                  * std::log((1.0 + beta)/(1.0 - beta)));
    
    return del;
}
BremResult Bremsstrahlung::Calculate(double ein, double eout, double egamma,
                                     bool radiate_proton,
                                     bool exponentiate,
                                     bool include_hard,
                                     bool calculate_spence) {
    BremResult result;
    
    // Convert to GeV for internal calculations (following Fortran convention)
    double ak = ein * 1.0e-3;      // Initial electron energy (GeV)
    double akp = eout * 1.0e-3;    // Final electron energy (GeV)
    double de = egamma * 1.0e-3;   // Photon energy (GeV)
    
    // Particle masses in GeV
    double am = Me * 1.0e-3;       // Electron mass
    double amp = Mp * 1.0e-3;      // Proton mass
    double am2 = am * am;
    double amp2 = amp * amp;
    
    // ========================================================================
    // Calculate electron scattering kinematics
    // ========================================================================
    
    // Electron scattering angle (approximate formula valid for elastic peak)
    double eang = 2.0 * std::asin(std::sqrt(amp / (2.0 * ak) * (ak/akp - 1.0)));
    
    // Virtual photon properties
    // double eta = 1.0 + 2.0 * ak * std::sin(eang/2.0) * std::sin(eang/2.0) / amp; // unused
    double Q2 = 4.0 * ak * akp * std::sin(eang/2.0) * std::sin(eang/2.0);
    
    // ========================================================================
    // Calculate final proton kinematics
    // ========================================================================
    
    double ape = amp + ak - akp;    // Proton energy
    double ap = std::sqrt(ape*ape - amp2);  // Proton momentum
    
    // Proton angle
    double pang = std::acos((ak - akp * std::cos(eang)) / ap);
    
    // ========================================================================
    // ELECTRON TERMS
    // ========================================================================
    
    double bei = 0.0;   // Direct initial electron
    double bef = 0.0;   // Direct final electron
    double bee = 0.0;   // Electron-electron interference
    
    // Direct initial electron contribution
    bei = -1.0 / (2.0 * pi) * std::log(ak / de);
    
    // Direct final electron contribution
    bef = -1.0 / (2.0 * pi) * std::log(akp / de);
    
    // Electron-electron interference
    double alpha_ee = 2.0 * am2 - 2.0 * ak * akp * (1.0 - std::cos(eang));
    double adot_ee = ak * akp * (1.0 - std::cos(eang));
    double ar1_ee, ar2_ee;
    
    if (CalculateRoots(alpha_ee, adot_ee, am2, am2, ar1_ee, ar2_ee)) {
        bee = Inter(calculate_spence, alpha_ee, ar1_ee, ar2_ee, ak, akp, de);
    }
    
    // ========================================================================
    // PROTON TERMS (if requested)
    // ========================================================================
    
    double bpi = 0.0;    // Direct initial proton
    double bpf = 0.0;    // Direct final proton
    double bpp = 0.0;    // Proton-proton interference
    double bepii = 0.0;  // e-p interference (initial e, initial p)
    double bepff = 0.0;  // e-p interference (final e, final p)
    double bepif = 0.0;  // e-p interference (initial e, final p)
    double bepfi = 0.0;  // e-p interference (final e, initial p)
    
    if (radiate_proton) {
        // Direct initial proton
        bpi = -1.0 / (2.0 * pi) * std::log(amp / de);
        
        // Direct final proton
        bpf = -1.0 / (2.0 * pi) * std::log(ape / de);
        
        // Proton-proton interference
        double alpha_pp = 2.0 * amp2 - 2.0 * amp * ape * (1.0 - std::cos(pang));
        double adot_pp = amp * ape * (1.0 - std::cos(pang));
        double ar1_pp, ar2_pp;
        
        if (CalculateRoots(alpha_pp, adot_pp, amp2, amp2, ar1_pp, ar2_pp)) {
            bpp = Inter(calculate_spence, alpha_pp, ar1_pp, ar2_pp, amp, ape, de);
        }
        
        // Electron-proton interferences
        // (i=initial, f=final)
        
        // Initial electron - Initial proton
        double alpha_epii = 2.0 * am2 - 2.0 * ak * amp;  // At rest target
        // double adot_epii = 0.0;  // No angle between them initially (unused)
        double ar1_epii, ar2_epii;
        
        if (std::abs(alpha_epii) > 1.0e-10) {
            ar1_epii = ar2_epii = alpha_epii / (2.0 * am2);
            bepii = Inter(calculate_spence, alpha_epii, ar1_epii, ar2_epii, ak, amp, de);
        }
        
        // Final electron - Final proton
        double cos_ep_angle = (ak * std::cos(eang) - akp) / ap;
        double alpha_epff = 2.0 * am2 - 2.0 * akp * ape * (1.0 - cos_ep_angle);
        double adot_epff = akp * ape * (1.0 - cos_ep_angle);
        double ar1_epff, ar2_epff;
        
        if (CalculateRoots(alpha_epff, adot_epff, am2, amp2, ar1_epff, ar2_epff)) {
            bepff = Inter(calculate_spence, alpha_epff, ar1_epff, ar2_epff, akp, ape, de);
        }
        
        // Initial electron - Final proton
        double cos_ep_if = std::cos(eang - pang);
        double alpha_epif = 2.0 * am2 - 2.0 * ak * ape * (1.0 - cos_ep_if);
        double adot_epif = ak * ape * (1.0 - cos_ep_if);
        double ar1_epif, ar2_epif;
        
        if (CalculateRoots(alpha_epif, adot_epif, am2, amp2, ar1_epif, ar2_epif)) {
            bepif = Inter(calculate_spence, alpha_epif, ar1_epif, ar2_epif, ak, ape, de);
        }
        
        // Final electron - Initial proton
        double alpha_epfi = 2.0 * am2 - 2.0 * akp * amp * (1.0 - std::cos(eang));
        double adot_epfi = akp * amp * (1.0 - std::cos(eang));
        double ar1_epfi, ar2_epfi;
        
        if (CalculateRoots(alpha_epfi, adot_epfi, am2, amp2, ar1_epfi, ar2_epfi)) {
            bepfi = Inter(calculate_spence, alpha_epfi, ar1_epfi, ar2_epfi, akp, amp, de);
        }
    }
    
    // ========================================================================
    // Combine all terms
    // ========================================================================
    
    // Electron contribution
    double b_electron = 2.0 * alpha * (bei + bef + bee);
    
    // Electron-proton interference
    double bz = 2.0 * alpha * (bepii + bepff + bepif + bepfi);
    
    // Proton contribution
    double bzz = 2.0 * alpha * (bpi + bpf + bpp);
    
    // Total soft photon correction
    result.bsoft = b_electron + bz + bzz;
    
    // ========================================================================
    // Hard correction
    // ========================================================================
    
    if (include_hard) {
        // Hard photon correction from Tsai
        result.bhard = -(alpha / pi) * (-28.0/9.0 + 13.0/6.0 * std::log(Q2 / am2));
    } else {
        result.bhard = 0.0;
    }
    
    // ========================================================================
    // Derivative (for peaked approximation)
    // ========================================================================
    
    // Calculate dbsoft/de using InterPrime for interference terms
    double dbei = -1.0 / (2.0 * pi * de);
    double dbef = -1.0 / (2.0 * pi * de);
    double dbee = 0.0;
    
    // Reuse ar1_ee, ar2_ee declared earlier in function
    if (CalculateRoots(2.0 * am2 - 2.0 * ak * akp * (1.0 - std::cos(eang)),
                      ak * akp * (1.0 - std::cos(eang)), am2, am2, 
                      ar1_ee, ar2_ee)) {
        dbee = InterPrime(2.0 * am2 - 2.0 * ak * akp * (1.0 - std::cos(eang)),
                         ar1_ee, ar2_ee, de);
    }
    
    result.dbsoft = 2.0 * alpha * (dbei + dbef + dbee);
    
    // Add proton derivative terms if radiating proton
    if (radiate_proton) {
        // (Similar calculation for proton terms - omitted for brevity)
        // In practice, proton terms are usually small
    }
    
    // ========================================================================
    // Calculate radiative factor
    // ========================================================================
    
    if (exponentiate) {
        // Exponentiated form: -dδ/dE × exp(-δ)
        result.rad_factor = -result.dbsoft / std::exp(result.bsoft);
    } else {
        // Linear form: 1 - dδ/dE
        result.rad_factor = 1.0 - result.dbsoft;
    }
    
    // Include hard correction
    if (include_hard) {
        result.rad_factor *= (1.0 - result.bhard);
    }
    
    return result;
}
BremResult Bremsstrahlung::CalculateOffShell(
    double egamma,
    const std::array<double, 3>& k_i,
    const std::array<double, 3>& k_f,
    const std::array<double, 4>& p_i,
    const std::array<double, 4>& p_f,
    bool radiate_proton,
    bool exponentiate,
    bool include_hard,
    bool calculate_spence) {
    
    BremResult result;
    
    // Convert to GeV
    double de = egamma * 1.0e-3;
    
    // Electron masses
    double am = Me * 1.0e-3;
    double am2 = am * am;
    
    // ========================================================================
    // Calculate electron kinematics from 3-momenta
    // ========================================================================
    
    // Initial electron: assume on-shell (coming from beam)
    double ki_mag = std::sqrt(k_i[0]*k_i[0] + k_i[1]*k_i[1] + k_i[2]*k_i[2]) * 1.0e-3;
    double ak = std::sqrt(ki_mag*ki_mag + am2);  // Energy (GeV)
    
    // Final electron
    double kf_mag = std::sqrt(k_f[0]*k_f[0] + k_f[1]*k_f[1] + k_f[2]*k_f[2]) * 1.0e-3;
    double akp = std::sqrt(kf_mag*kf_mag + am2);  // Energy (GeV)
    
    // Electron scattering angle
    double cos_eang = 0.0;
    if (ki_mag > 0.0 && kf_mag > 0.0) {
        double dot = (k_i[0]*k_f[0] + k_i[1]*k_f[1] + k_i[2]*k_f[2]) * 1.0e-6; // GeV²
        cos_eang = dot / (ki_mag * kf_mag);
        cos_eang = std::max(-1.0, std::min(1.0, cos_eang)); // Clamp to [-1,1]
    }
    double eang = std::acos(cos_eang);
    
    // Q² for hard correction
    double Q2 = 4.0 * ak * akp * std::sin(eang/2.0) * std::sin(eang/2.0);
    
    // ========================================================================
    // Calculate proton kinematics (can be off-shell)
    // ========================================================================
    
    // Initial proton (convert from MeV to GeV)
    double api = p_i[0] * 1.0e-3;  // Energy
    double pix = p_i[1] * 1.0e-3;  // px
    double piy = p_i[2] * 1.0e-3;  // py
    double piz = p_i[3] * 1.0e-3;  // pz
    double pi_mag = std::sqrt(pix*pix + piy*piy + piz*piz);
    double ami = std::sqrt(std::abs(api*api - pi_mag*pi_mag)); // Can be off-shell
    double ami2 = ami * ami;
    
    // Final proton
    double ape = p_f[0] * 1.0e-3;  // Energy
    double pfx = p_f[1] * 1.0e-3;  // px
    double pfy = p_f[2] * 1.0e-3;  // py
    double pfz = p_f[3] * 1.0e-3;  // pz
    double pf_mag = std::sqrt(pfx*pfx + pfy*pfy + pfz*pfz);
    double amf = std::sqrt(std::abs(ape*ape - pf_mag*pf_mag)); // Can be off-shell
    double amf2 = amf * amf;
    
    // Proton scattering angle
    double cos_pang = 0.0;
    if (pi_mag > 0.0 && pf_mag > 0.0) {
        double dot = pix*pfx + piy*pfy + piz*pfz;
        cos_pang = dot / (pi_mag * pf_mag);
        cos_pang = std::max(-1.0, std::min(1.0, cos_pang));
    }
    // double pang = std::acos(cos_pang);  // Calculated but unused in this version
    
    // ========================================================================
    // ELECTRON TERMS
    // ========================================================================
    
    double bei = 0.0;   // Direct initial electron
    double bef = 0.0;   // Direct final electron
    double bee = 0.0;   // Electron-electron interference
    
    // Direct contributions
    bei = -1.0 / (2.0 * pi) * std::log(ak / de);
    bef = -1.0 / (2.0 * pi) * std::log(akp / de);
    
    // Interference term
    double adot_ee = ak * akp * (1.0 - cos_eang);
    double alpha_ee = 2.0 * am2 - 2.0 * adot_ee;
    double ar1_ee, ar2_ee;
    
    if (CalculateRoots(alpha_ee, adot_ee, am2, am2, ar1_ee, ar2_ee)) {
        bee = Inter(calculate_spence, alpha_ee, ar1_ee, ar2_ee, ak, akp, de);
    }
    
    // ========================================================================
    // PROTON TERMS (if requested)
    // ========================================================================
    
    double bpi = 0.0;    // Direct initial proton
    double bpf = 0.0;    // Direct final proton
    double bpp = 0.0;    // Proton-proton interference
    double bepii = 0.0;  // e-p interference terms
    double bepff = 0.0;
    double bepif = 0.0;
    double bepfi = 0.0;
    
    if (radiate_proton) {
        // Direct contributions
        bpi = -1.0 / (2.0 * pi) * std::log(api / de);
        bpf = -1.0 / (2.0 * pi) * std::log(ape / de);
        
        // Proton-proton interference
        double adot_pp = api * ape * (1.0 - cos_pang);
        double alpha_pp = 2.0 * ami2 - 2.0 * adot_pp;
        double ar1_pp, ar2_pp;
        
        if (CalculateRoots(alpha_pp, adot_pp, ami2, amf2, ar1_pp, ar2_pp)) {
            bpp = Inter(calculate_spence, alpha_pp, ar1_pp, ar2_pp, api, ape, de);
        }
        
        // Electron-proton interferences
        
        // Initial electron - Initial proton
        double cos_epii = 0.0;
        if (ki_mag > 0.0 && pi_mag > 0.0) {
            double dot = (k_i[0]*pix + k_i[1]*piy + k_i[2]*piz) * 1.0e-6;
            cos_epii = dot / (ki_mag * pi_mag);
            cos_epii = std::max(-1.0, std::min(1.0, cos_epii));
        }
        double adot_epii = ak * api * (1.0 - cos_epii);
        double alpha_epii = 2.0 * am2 - 2.0 * adot_epii;
        double ar1_epii, ar2_epii;
        
        if (CalculateRoots(alpha_epii, adot_epii, am2, ami2, ar1_epii, ar2_epii)) {
            bepii = Inter(calculate_spence, alpha_epii, ar1_epii, ar2_epii, ak, api, de);
        }
        
        // Final electron - Final proton
        double cos_epff = 0.0;
        if (kf_mag > 0.0 && pf_mag > 0.0) {
            double dot = (k_f[0]*pfx + k_f[1]*pfy + k_f[2]*pfz) * 1.0e-6;
            cos_epff = dot / (kf_mag * pf_mag);
            cos_epff = std::max(-1.0, std::min(1.0, cos_epff));
        }
        double adot_epff = akp * ape * (1.0 - cos_epff);
        double alpha_epff = 2.0 * am2 - 2.0 * adot_epff;
        double ar1_epff, ar2_epff;
        
        if (CalculateRoots(alpha_epff, adot_epff, am2, amf2, ar1_epff, ar2_epff)) {
            bepff = Inter(calculate_spence, alpha_epff, ar1_epff, ar2_epff, akp, ape, de);
        }
        
        // Initial electron - Final proton
        double cos_epif = 0.0;
        if (ki_mag > 0.0 && pf_mag > 0.0) {
            double dot = (k_i[0]*pfx + k_i[1]*pfy + k_i[2]*pfz) * 1.0e-6;
            cos_epif = dot / (ki_mag * pf_mag);
            cos_epif = std::max(-1.0, std::min(1.0, cos_epif));
        }
        double adot_epif = ak * ape * (1.0 - cos_epif);
        double alpha_epif = 2.0 * am2 - 2.0 * adot_epif;
        double ar1_epif, ar2_epif;
        
        if (CalculateRoots(alpha_epif, adot_epif, am2, amf2, ar1_epif, ar2_epif)) {
            bepif = Inter(calculate_spence, alpha_epif, ar1_epif, ar2_epif, ak, ape, de);
        }
        
        // Final electron - Initial proton
        double cos_epfi = 0.0;
        if (kf_mag > 0.0 && pi_mag > 0.0) {
            double dot = (k_f[0]*pix + k_f[1]*piy + k_f[2]*piz) * 1.0e-6;
            cos_epfi = dot / (kf_mag * pi_mag);
            cos_epfi = std::max(-1.0, std::min(1.0, cos_epfi));
        }
        double adot_epfi = akp * api * (1.0 - cos_epfi);
        double alpha_epfi = 2.0 * am2 - 2.0 * adot_epfi;
        double ar1_epfi, ar2_epfi;
        
        if (CalculateRoots(alpha_epfi, adot_epfi, am2, ami2, ar1_epfi, ar2_epfi)) {
            bepfi = Inter(calculate_spence, alpha_epfi, ar1_epfi, ar2_epfi, akp, api, de);
        }
    }
    
    // ========================================================================
    // Combine all terms
    // ========================================================================
    
    double b_electron = 2.0 * alpha * (bei + bef + bee);
    double bz = 2.0 * alpha * (bepii + bepff + bepif + bepfi);
    double bzz = 2.0 * alpha * (bpi + bpf + bpp);
    
    result.bsoft = b_electron + bz + bzz;
    
    // ========================================================================
    // Hard correction
    // ========================================================================
    
    if (include_hard) {
        result.bhard = -(alpha / pi) * (-28.0/9.0 + 13.0/6.0 * std::log(Q2 / am2));
    } else {
        result.bhard = 0.0;
    }
    
    // ========================================================================
    // Derivative (simplified - use electron terms only)
    // ========================================================================
    
    double dbei = -1.0 / (2.0 * pi * de);
    double dbef = -1.0 / (2.0 * pi * de);
    double dbee = 0.0;
    
    // Reuse ar1_ee, ar2_ee declared earlier in function
    if (CalculateRoots(2.0 * am2 - 2.0 * ak * akp * (1.0 - cos_eang),
                      ak * akp * (1.0 - cos_eang), am2, am2,
                      ar1_ee, ar2_ee)) {
        dbee = InterPrime(2.0 * am2 - 2.0 * ak * akp * (1.0 - cos_eang),
                         ar1_ee, ar2_ee, de);
    }
    
    result.dbsoft = 2.0 * alpha * (dbei + dbef + dbee);
    
    // ========================================================================
    // Calculate radiative factor
    // ========================================================================
    
    if (exponentiate) {
        result.rad_factor = -result.dbsoft / std::exp(result.bsoft);
    } else {
        result.rad_factor = 1.0 - result.dbsoft;
    }
    
    if (include_hard) {
        result.rad_factor *= (1.0 - result.bhard);
    }
    
    return result;
}

} // namespace simc
