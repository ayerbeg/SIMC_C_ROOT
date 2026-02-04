// src/physics/RadiativeCorrections.cpp
// Main radiative corrections implementation - Part 1: Constructor and helpers

#include "simc/RadiativeCorrections.h"
#include "simc/Bremsstrahlung.h"
#include "simc/ExternalRadiation.h"
#include "simc/MathUtils.h"
#include <cmath>
#include <stdexcept>

namespace simc {

using constants::alpha;
using constants::Me;
using constants::Mp;
using constants::Mn;
using constants::pi;
using constants::euler;

// ============================================================================
// Constructor / Destructor
// ============================================================================

RadiativeCorrections::RadiativeCorrections(RandomGenerator& rng)
    : rng_(rng), brem_(nullptr), extrad_(nullptr) {
    // Initialize with default configuration
    state_.Reset();
}

RadiativeCorrections::~RadiativeCorrections() {
    // Cleanup (brem_ and extrad_ are just for organization, 
    // actual functions are static)
}

// ============================================================================
// Initialize
// ============================================================================

void RadiativeCorrections::Initialize(const RadConfig& config) {
    config_ = config;
    state_.Reset();
}

// ============================================================================
// BasicRad - Generate photon energy using power law
// ============================================================================

bool RadiativeCorrections::BasicRad(int itail, double Egamma_lo, double Egamma_hi,
                                    double& Egamma, double& weight, 
                                    double& val_reciprocal) {
    // Map itail=0 to itail=4 (for both tails case)
    int idx = (itail == 0) ? 4 : itail;
    
    // Check if radiation is enabled for this tail
    if (config_.g[idx] <= 0.0) {
        Egamma = 0.0;
        weight = 1.0;
        val_reciprocal = 0.0;
        return true;
    }
    
    // Check valid energy range
    if (Egamma_hi <= Egamma_lo || Egamma_hi <= 0.0) {
        Egamma = 0.0;
        weight = 0.0;
        val_reciprocal = 0.0;
        return false;
    }
    
    // Calculate powers: E^g
    double g = config_.g[idx];
    double power_hi = std::pow(Egamma_hi, g);
    double power_lo = (Egamma_lo > 0.0) ? std::pow(Egamma_lo, g) : 0.0;
    
    // Generate random y in [ymin, 1] where ymin = power_lo/power_hi
    double ymin = power_lo / power_hi;
    double y = ymin + rng_.Uniform() * (1.0 - ymin);
    
    // Invert to get x, then Egamma
    double x = std::pow(y, 1.0 / g);
    Egamma = x * Egamma_hi;
    
    // Calculate val_reciprocal = 1 / [generating function value]
    // Generating function: f(E) = g * E^(g-1) / (E_max^g - E_min^g)
    if (Egamma > 0.0) {
        val_reciprocal = std::pow(Egamma, 1.0 - g) * (power_hi - power_lo) / g;
    } else {
        val_reciprocal = 0.0;
    }
    
    // Calculate weight = normalization constant × (E_max^g - E_min^g)
    weight = config_.c[idx] / g * (power_hi - power_lo);
    
    return true;
}

// ============================================================================
// LambdaDave - Extended peaking approximation lambda
// ============================================================================

double RadiativeCorrections::LambdaDave(int itail, bool plus_flag, 
                                        bool doing_proton,
                                        double e1, double e2, double e3, 
                                        double p3, double theta) {
    double lambda = 0.0;
    
    if (itail == 1) {
        // Incoming electron
        lambda = (alpha / pi) * (2.0 * std::log(2.0 * e1 / Me) - 1.0);
        
    } else if (itail == 2) {
        // Outgoing electron  
        lambda = (alpha / pi) * (2.0 * std::log(2.0 * e2 / Me) - 1.0);
        
    } else if (itail == 3) {
        // Outgoing hadron
        if (p3 > 0.0 && e3 > 0.0) {
            lambda = (alpha / pi) * ((e3 / p3) * std::log((e3 + p3) / (e3 - p3)) - 2.0);
            
            // Safety check: if lambda comes out negative, set to zero
            if (lambda < 0.0) {
                lambda = 0.0;
            }
        }
    }
    
    // Add plus term if requested (for itail 1 or 2)
    if (plus_flag && itail < 3) {
        double plus_term = std::log((1.0 - std::cos(theta)) / 2.0);
        
        if (doing_proton) {
            plus_term += 2.0 * std::log(e1 / e2);
        }
        
        lambda += (alpha / pi) * plus_term;
    }
    
    return lambda;
}

// ============================================================================
// Schwinger - Schwinger correction (leading order QED)
// ============================================================================

double RadiativeCorrections::Schwinger(double Ecutoff, const SimcEvent& vertex,
                                       double target_mass, bool include_hard,
                                       double& dsoft, double& dhard) {
    // Calculate kinematics
    double lq = std::log(vertex.Q2 / (Me * Me)) - 1.0;
    double s2 = std::sin(vertex.e_theta / 2.0);
    s2 = s2 * s2;
    
    double b = 1.0 + 2.0 * vertex.nu * s2 / (target_mass);
    
    // Spence function term
    double spence = MathUtils::Spence(1.0 - s2) - 2.5893784;
    
    // Soft correction
    double eta_sq = config_.etta * config_.etta;
    dsoft = (alpha / pi) * lq * std::log(vertex.Ein / (eta_sq * vertex.e_E * b / (Ecutoff * Ecutoff)));
    
    // Hard correction
    dhard = 0.0;
    if (include_hard) {
        double log_ratio = std::log(vertex.Ein / vertex.e_E);
        dhard = -(alpha / pi) * (2.166666 * lq + spence - log_ratio * log_ratio / 2.0);
    }
    
    // Calculate Schwinger correction factor
    double schwinger = 1.0;
    
    if (config_.use_expon) {
        schwinger = std::exp(dsoft);
    } else {
        schwinger = 1.0 + dsoft;
    }
    
    if (include_hard) {
        schwinger = schwinger / (1.0 - dhard);
    }
    
    return schwinger;
}
using constants::alpha;
using constants::euler;
using constants::pi;

double RadiativeCorrections::PeakedRadWeight(const SimcEvent& vertex, 
                                             double Egamma,
                                             double emin, double emax,
                                             double /*basicrad_val_reciprocal*/,
                                             double basicrad_weight) {
    double weight = 1.0;
    
    // ========================================================================
    // Calculate internal corrections (exact bremsstrahlung)
    // ========================================================================
    
    double bsoft_int_max = 0.0;
    double bsoft_int_min = 0.0;
    // double dbsoft_int = 0.0;  // Not currently used in weight calculation
    
    if (config_.intcor_mode == 0) {
        // Use exact bremsstrahlung calculation
        
        BremResult brem_max, brem_min;
        
        if (config_.use_offshell_rad) {
            // Off-shell calculation - would need to construct momentum arrays
            // For now, fall back to on-shell
            // TODO: Implement off-shell path when Fermi momentum is available
            brem_max = Bremsstrahlung::Calculate(vertex.Ein, emax, Egamma,
                                                config_.rad_proton_this_ev,
                                                config_.exponentiate,
                                                config_.include_hard,
                                                config_.calculate_spence);
            brem_min = Bremsstrahlung::Calculate(vertex.Ein, emin, Egamma,
                                                config_.rad_proton_this_ev,
                                                config_.exponentiate,
                                                config_.include_hard,
                                                config_.calculate_spence);
        } else {
            // On-shell calculation
            brem_max = Bremsstrahlung::Calculate(vertex.Ein, emax, Egamma,
                                                config_.rad_proton_this_ev,
                                                config_.exponentiate,
                                                config_.include_hard,
                                                config_.calculate_spence);
            brem_min = Bremsstrahlung::Calculate(vertex.Ein, emin, Egamma,
                                                config_.rad_proton_this_ev,
                                                config_.exponentiate,
                                                config_.include_hard,
                                                config_.calculate_spence);
        }
        
        bsoft_int_max = brem_max.bsoft;
        bsoft_int_min = brem_min.bsoft;
        // dbsoft_int = brem_max.dbsoft;  // Derivative not used in current weight formula
        
    } else {
        // Skip exact calculation (intcor_mode = 1)
        bsoft_int_max = 0.0;
        bsoft_int_min = 0.0;
        // dbsoft_int = 0.0;
    }
    
    // ========================================================================
    // Calculate external corrections
    // ========================================================================
    
    double phi_ext_max = 1.0;
    double phi_ext_min = 1.0;
    
    if (config_.extrad_flag <= 2) {
        // Use ExternalRadiation::Phi for modes 1 and 2
        phi_ext_max = ExternalRadiation::Phi(0, vertex.Ein, emax, Egamma,
                                            config_.extrad_flag,
                                            config_.bt, config_.g,
                                            config_.etatzai);
        phi_ext_min = ExternalRadiation::Phi(0, vertex.Ein, emin, Egamma,
                                            config_.extrad_flag,
                                            config_.bt, config_.g,
                                            config_.etatzai);
                                            
    } else if (config_.extrad_flag == 3) {
        // Friedrich approximation
        double trad1 = config_.bt[0] / config_.etatzai;
        // double trad2 = config_.bt[1] / config_.etatzai;  // unused in current implementation
        
        ExtRadResult ext_max = ExternalRadiation::Friedrich(vertex.Ein, emax, 
                                                           trad1, config_.etatzai);
        ExtRadResult ext_min = ExternalRadiation::Friedrich(vertex.Ein, emin,
                                                           trad1, config_.etatzai);
        
        // Convert to phi factors (exp form)
        phi_ext_max = std::exp(ext_max.dbrem);
        phi_ext_min = std::exp(ext_min.dbrem);
    }
    
    // ========================================================================
    // Combine corrections using Mo & Tsai prescription
    // ========================================================================
    
    // Get g values for the combined case
    double g_tot = config_.g_int_val + config_.g_ext_val;
    // double c_tot = config_.c_int[0] + config_.c_ext[0];  // unused in current formula
    
    // Peaked approximation weight formula:
    // w = (c_ext/g_ext) × [exp(-δ_int_max) × E_max^g_ext - exp(-δ_int_min) × E_min^g_ext]
    //     × exp(-γ_E × g) / Γ(1+g) × [Gamma function products]
    //     × φ_ext
    
    double term1 = std::exp(-bsoft_int_max) * std::pow(emax, config_.g_ext_val);
    double term2 = std::exp(-bsoft_int_min) * std::pow(emin, config_.g_ext_val);
    
    weight = (config_.c_ext[0] / config_.g_ext_val) * (term1 - term2);
    
    // Add Euler constant and Gamma function terms
    double euler_term = std::exp(-euler * g_tot);
    double gamma_denom = MathUtils::Gamma(1.0 + g_tot);
    
    weight *= euler_term / gamma_denom;
    
    // Multiply by radiation thickness Gamma functions
    if (config_.bt[0] > 0.0 || config_.bt[1] > 0.0) {
        double bt_tot = config_.bt[0] + config_.bt[1];
        double gamma_num = MathUtils::Gamma(1.0 + g_tot - bt_tot)
                         * MathUtils::Gamma(1.0 + config_.bt[0])
                         * MathUtils::Gamma(1.0 + config_.bt[1]);
        double gamma_den = MathUtils::Gamma(1.0 + g_tot);
        
        weight *= gamma_num / gamma_den;
    }
    
    // Multiply by phi corrections
    weight *= (phi_ext_max + phi_ext_min) / 2.0;  // Average phi
    
    // Correct for basicrad generation
    weight /= basicrad_weight;
    
    return weight;
}
using constants::Me;

bool RadiativeCorrections::CompleteEvent(SimcEvent& vertex,
                                         double target_mass,
                                         double hadron_mass,
                                         bool doing_hyd_elast,
                                         bool doing_deuterium,
                                         bool doing_pion,
                                         bool doing_kaon,
                                         bool doing_delta) {
    
    // ========================================================================
    // Calculate unit vectors for electron
    // ========================================================================
    
    vertex.ue_x = std::sin(vertex.e_theta) * std::cos(vertex.e_phi);
    vertex.ue_y = std::sin(vertex.e_theta) * std::sin(vertex.e_phi);
    vertex.ue_z = std::cos(vertex.e_theta);
    
    // ========================================================================
    // Calculate basic kinematics
    // ========================================================================
    
    vertex.nu = vertex.Ein - vertex.e_E;
    vertex.Q2 = 2.0 * vertex.Ein * vertex.e_E * (1.0 - vertex.ue_z);
    vertex.q = std::sqrt(vertex.Q2 + vertex.nu * vertex.nu);
    
    // Virtual photon unit vector
    vertex.uq_x = -vertex.e_E * vertex.ue_x / vertex.q;
    vertex.uq_y = -vertex.e_E * vertex.ue_y / vertex.q;
    vertex.uq_z = (vertex.Ein - vertex.e_E * vertex.ue_z) / vertex.q;
    
    // ========================================================================
    // H(e,e'p) Elastic - Calculate proton from momentum conservation
    // ========================================================================
    
    if (doing_hyd_elast) {
        // For elastic scattering: p = q (momentum conservation)
        vertex.up_x = vertex.uq_x;
        vertex.up_y = vertex.uq_y;
        vertex.up_z = vertex.uq_z;
        
        vertex.p_P = vertex.q;
        vertex.p_E = std::sqrt(vertex.p_P * vertex.p_P + hadron_mass * hadron_mass);
        
        vertex.p_theta = std::acos(vertex.up_z);
        vertex.p_phi = std::atan2(vertex.up_y, vertex.up_x);
        if (vertex.p_phi < 0.0) vertex.p_phi += 2.0 * constants::pi;
        
        vertex.Em = 0.0;
        vertex.Pm = 0.0;
        
        return true;
    }
    
    // ========================================================================
    // D(e,e'p) - Solve for proton momentum
    // ========================================================================
    
    if (doing_deuterium) {
        // Binding energy
        vertex.Em = Mp + Mn - target_mass;  // ~2.22 MeV for deuteron
        
        // Recoil mass (neutron for D(e,e'p))
        double Mrec = target_mass - Mp + vertex.Em;  // = Mn
        
        // Solve quadratic equation for proton energy
        // Energy conservation: Ein + M = Ee + Ep + sqrt(Pm² + Mrec²)
        // Momentum conservation: q = pp + pm
        // This gives a quadratic in Ep
        
        double a = -vertex.q * (vertex.uq_x * vertex.up_x 
                              + vertex.uq_y * vertex.up_y 
                              + vertex.uq_z * vertex.up_z);
        double b = vertex.q * vertex.q;
        double c = vertex.nu + target_mass;
        double t = c*c - b + hadron_mass*hadron_mass - Mrec*Mrec;
        
        double QA = 4.0 * (a*a - c*c);
        double QB = 4.0 * c * t;
        double QC = -4.0 * a*a * hadron_mass*hadron_mass - t*t;
        
        double radical = QB*QB - 4.0*QA*QC;
        if (radical < 0.0) return false;  // Unphysical
        
        // Take the positive (forward-going) solution
        vertex.p_E = (-QB - std::sqrt(radical)) / (2.0 * QA);
        
        if (vertex.p_E <= hadron_mass) return false;  // Below mass threshold
        
        vertex.p_P = std::sqrt(vertex.p_E*vertex.p_E - hadron_mass*hadron_mass);
        
        return true;
    }
    
    // ========================================================================
    // Pion/Kaon Production - Similar to deuterium but with pion/kaon mass
    // ========================================================================
    
    if (doing_pion || doing_kaon || doing_delta) {
        // Similar quadratic solve as deuterium
        // For now, use hadron mass and target mass provided
        
        double Mrec_struck = Mp;  // Struck nucleon mass (simplification)
        
        double a = -vertex.q * (vertex.uq_x * vertex.up_x 
                              + vertex.uq_y * vertex.up_y 
                              + vertex.uq_z * vertex.up_z);
        double b = vertex.q * vertex.q;
        double c = vertex.nu + target_mass;
        double t = c*c - b + hadron_mass*hadron_mass - Mrec_struck*Mrec_struck;
        
        double QA = 4.0 * (a*a - c*c);
        double QB = 4.0 * c * t;
        double QC = -4.0 * a*a * hadron_mass*hadron_mass - t*t;
        
        double radical = QB*QB - 4.0*QA*QC;
        if (radical < 0.0) return false;
        
        vertex.p_E = (-QB - std::sqrt(radical)) / (2.0 * QA);
        
        double E_rec = c - vertex.p_E;
        if (E_rec <= Mrec_struck) return false;  // Unphysical
        
        if (vertex.p_E <= hadron_mass) return false;
        
        vertex.p_P = std::sqrt(vertex.p_E*vertex.p_E - hadron_mass*hadron_mass);
        
        return true;
    }
    
    // ========================================================================
    // Calculate missing momentum for all cases
    // ========================================================================
    
    vertex.Pmx = vertex.p_P * vertex.up_x - vertex.q * vertex.uq_x;
    vertex.Pmy = vertex.p_P * vertex.up_y - vertex.q * vertex.uq_y;
    vertex.Pmz = vertex.p_P * vertex.up_z - vertex.q * vertex.uq_z;
    vertex.Pm = std::sqrt(vertex.Pmx*vertex.Pmx + vertex.Pmy*vertex.Pmy + vertex.Pmz*vertex.Pmz);
    
    vertex.Em = vertex.nu + target_mass - vertex.p_E;
    
    // ========================================================================
    // Calculate W (invariant mass)
    // ========================================================================
    
    double W2 = target_mass*target_mass + 2.0*target_mass*vertex.nu - vertex.Q2;
    vertex.W = (W2 > 0.0) ? std::sqrt(W2) : -std::sqrt(-W2);
    
    return true;
}
using constants::pi;

bool RadiativeCorrections::Generate(SimcEvent& vertex, SimcEvent& orig,
                                    double target_mass, double hadron_mass,
                                    double& gen_weight,
                                    bool doing_heavy, bool doing_deuterium,
                                    bool doing_hyd_elast,
                                    bool doing_pion, bool doing_kaon,
                                    bool doing_delta) {
    
    // Reset state
    state_.Reset();
    state_.Egamma_used.fill(0.0);
    double rad_weight = 1.0;
    
    // ========================================================================
    // Decide which tail(s) to radiate
    // ========================================================================
    
    int ntail = config_.ntail;  // 0=decide randomly, 1/2/3=specific tail
    
    if (ntail == 0 && config_.one_tail == 0) {
        // Decide randomly which tails radiate based on frac probabilities
        double rn = rng_.Uniform();
        
        // Check probabilities for each tail
        if (config_.doing_tail[0] && rn < config_.frac[0]) {
            ntail = 1;  // Incoming electron only
        } else if (config_.doing_tail[1] && rn < config_.frac[0] + config_.frac[1]) {
            ntail = 2;  // Outgoing electron only
        } else if (config_.doing_tail[2] && rn < config_.frac[0] + config_.frac[1] + config_.frac[2]) {
            ntail = 3;  // Outgoing hadron only
        } else {
            ntail = 0;  // All tails simultaneously
        }
    }
    
    // ========================================================================
    // TAIL 1: Incoming Electron Radiation (before scattering)
    // ========================================================================
    
    if (config_.doing_tail[0] && (ntail == 0 || ntail == 1)) {
        
        double Egamma_min = 0.0;
        double Egamma_max = config_.Egamma1_max;
        
        // Calculate limits based on physics type
        if (doing_heavy) {
            // A(e,e'p): Limit based on Em range
            // Not fully implemented - would need VERTEXedge limits
            Egamma_min = 0.0;  // Placeholder
            Egamma_max = std::min(Egamma_max, 100.0);  // Placeholder
            
        } else if (doing_hyd_elast) {
            // H(e,e'p): Limit based on scattered electron energy range
            // Solve: Ee' = Ein*Mp / (Mp + Ein*(1-cos(theta)))
            // This gives limits on Ein that produce Ee' in acceptance
            
            // For now, use simple limits
            Egamma_max = std::min(Egamma_max, vertex.Ein * 0.5);
            
        } else if (doing_deuterium) {
            // D(e,e'p): Similar to hydrogen but with deuteron
            Egamma_max = std::min(Egamma_max, vertex.Ein * 0.5);
            
        } else if (doing_pion || doing_kaon || doing_delta) {
            // Pion/kaon production
            Egamma_max = std::min(Egamma_max, vertex.Ein - vertex.e_E);
        }
        
        // Enforce resolution limit
        Egamma_max = std::min(Egamma_max, config_.Egamma_res_limit);
        
        // Check valid range
        if (Egamma_max <= Egamma_min) {
            state_.success = false;
            return false;
        }
        
        // Generate photon energy
        double Egamma1, weight1, val_recip1;
        if (!BasicRad(1, Egamma_min, Egamma_max, Egamma1, weight1, val_recip1)) {
            state_.success = false;
            return false;
        }
        
        state_.Egamma_used[0] = Egamma1;
        state_.basicrad_weight[0] = weight1;
        state_.basicrad_val_reciprocal[0] = val_recip1;
        
        // Subtract radiation from incoming electron
        vertex.Ein -= Egamma1;
        
        // Recalculate kinematics
        if (!CompleteEvent(vertex, target_mass, hadron_mass,
                          doing_hyd_elast, doing_deuterium,
                          doing_pion, doing_kaon, doing_delta)) {
            state_.success = false;
            return false;
        }
        
        // Calculate radiative weight
        if (config_.rad_flag <= 1) {
            // Use peaked approximation with exact bremsstrahlung
            double peaked_weight = PeakedRadWeight(vertex, Egamma1,
                                                  Egamma_min, Egamma_max,
                                                  val_recip1, weight1);
            rad_weight *= peaked_weight;
        } else {
            // Use basicrad weight with phi correction
            rad_weight *= weight1;
            
            double phi = ExternalRadiation::Phi(1, vertex.Ein + Egamma1, 
                                               vertex.e_E, Egamma1,
                                               config_.extrad_flag,
                                               config_.bt, config_.g,
                                               config_.etatzai);
            rad_weight *= phi;
        }
    }
    
    // ========================================================================
    // TAIL 2: Outgoing Electron Radiation (after scattering)
    // ========================================================================
    
    if (config_.doing_tail[1] && (ntail == 0 || ntail == 2)) {
        
        double Egamma_min = 0.0;
        double Egamma_max = config_.Egamma2_max;
        
        // Limits from electron arm acceptance
        // The electron loses energy Egamma after scattering
        // So detected energy is vertex.e_E - Egamma
        // This must be within spectrometer acceptance
        
        // For now, use simple limit
        Egamma_max = std::min(Egamma_max, vertex.e_E * 0.5);
        Egamma_max = std::min(Egamma_max, config_.Egamma_res_limit);
        
        if (Egamma_max <= Egamma_min) {
            state_.success = false;
            return false;
        }
        
        // Generate photon energy
        double Egamma2, weight2, val_recip2;
        if (!BasicRad(2, Egamma_min, Egamma_max, Egamma2, weight2, val_recip2)) {
            state_.success = false;
            return false;
        }
        
        state_.Egamma_used[1] = Egamma2;
        state_.basicrad_weight[1] = weight2;
        state_.basicrad_val_reciprocal[1] = val_recip2;
        
        // For tail 2, we don't need to recalculate kinematics
        // (the scattering already happened, just photon emitted afterward)
        
        // Calculate radiative weight
        if (config_.rad_flag <= 1) {
            double peaked_weight = PeakedRadWeight(vertex, Egamma2,
                                                  Egamma_min, Egamma_max,
                                                  val_recip2, weight2);
            rad_weight *= peaked_weight;
        } else {
            rad_weight *= weight2;
            
            double phi = ExternalRadiation::Phi(2, vertex.Ein,
                                               vertex.e_E, Egamma2,
                                               config_.extrad_flag,
                                               config_.bt, config_.g,
                                               config_.etatzai);
            rad_weight *= phi;
        }
    }
    
    // ========================================================================
    // TAIL 3: Outgoing Hadron Radiation (after scattering)
    // ========================================================================
    
    if (config_.rad_proton_this_ev && config_.doing_tail[2] && 
        (ntail == 0 || ntail == 3)) {
        
        double Egamma_min = 0.0;
        double Egamma_max = config_.Egamma3_max;
        
        // Similar to tail 2 but for hadron
        Egamma_max = std::min(Egamma_max, vertex.p_E * 0.5);
        Egamma_max = std::min(Egamma_max, config_.Egamma_res_limit);
        
        if (Egamma_max <= Egamma_min) {
            state_.success = false;
            return false;
        }
        
        // Generate photon energy
        double Egamma3, weight3, val_recip3;
        if (!BasicRad(3, Egamma_min, Egamma_max, Egamma3, weight3, val_recip3)) {
            state_.success = false;
            return false;
        }
        
        state_.Egamma_used[2] = Egamma3;
        state_.basicrad_weight[2] = weight3;
        state_.basicrad_val_reciprocal[2] = val_recip3;
        
        // Calculate radiative weight
        if (config_.rad_flag <= 1) {
            double peaked_weight = PeakedRadWeight(vertex, Egamma3,
                                                  Egamma_min, Egamma_max,
                                                  val_recip3, weight3);
            rad_weight *= peaked_weight;
        } else {
            rad_weight *= weight3;
        }
    }
    
    // ========================================================================
    // Create 'orig' event (vertex + radiation added back)
    // ========================================================================
    
    orig = vertex;  // Copy all fields
    
    // Add radiation back
    orig.Ein = vertex.Ein + state_.Egamma_used[0];
    orig.e_E = vertex.e_E + state_.Egamma_used[1];
    orig.e_P = orig.e_E;  // Electron is relativistic
    
    orig.p_E = vertex.p_E + state_.Egamma_used[2];
    orig.p_P = std::sqrt(orig.p_E * orig.p_E - hadron_mass * hadron_mass);
    
    // Recalculate orig kinematics
    orig.nu = orig.Ein - orig.e_E;
    orig.Q2 = 2.0 * orig.Ein * orig.e_E * (1.0 - orig.ue_z);
    orig.q = std::sqrt(orig.Q2 + orig.nu * orig.nu);
    
    // Missing momentum
    orig.Pmx = orig.p_P * orig.up_x - orig.q * orig.uq_x;
    orig.Pmy = orig.p_P * orig.up_y - orig.q * orig.uq_y;
    orig.Pmz = orig.p_P * orig.up_z - orig.q * orig.uq_z;
    orig.Pm = std::sqrt(orig.Pmx*orig.Pmx + orig.Pmy*orig.Pmy + orig.Pmz*orig.Pmz);
    
    orig.Em = orig.nu + target_mass - orig.p_E;
    
    // W invariant mass
    double W2 = target_mass*target_mass + 2.0*target_mass*orig.nu - orig.Q2;
    orig.W = (W2 > 0.0) ? std::sqrt(W2) : -std::sqrt(-W2);
    
    // ========================================================================
    // Apply final weights
    // ========================================================================
    
    state_.rad_weight = rad_weight;
    
    // Apply hard correction factor
    gen_weight *= rad_weight / config_.hardcorfac;
    
    // Apply tail fraction if forcing single tail mode
    if (config_.one_tail != 0 && ntail == 2) {
        gen_weight *= config_.frac[1];  // Tail 2 fraction
    }
    
    state_.success = true;
    return true;
}

} // namespace simc
