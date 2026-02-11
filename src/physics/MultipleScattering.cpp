// src/physics/MultipleScattering.cpp
// Implementation of multiple scattering calculations
// Ported from musc.f and musc_ext.f

#include "simc/physics/MultipleScattering.h"
#include "simc/core/SimcConstants.h"
#include <cmath>
#include <array>
#include <algorithm>
#include <iostream>

namespace simc {
using namespace constants;

// ============================================================================
// Main Multiple Scattering Calculation
// ============================================================================

MultipleScattering::ScatterAngles MultipleScattering::Calculate(
    double length,
    const Material& material,
    double momentum,
    double charge,
    double mass,
    RandomGenerator& rng,
    bool use_moliere) {
    
    ScatterAngles angles;
    
    if (length <= 0.0) {
        return angles;
    }
    
    double rho = material.density;
    double X0_gcm2 = material.X0;
    double X0_cm = X0_gcm2 / rho;
    double thickness = length / X0_cm;
    
    angles.theta_rms = CalculateRMS(thickness, momentum, charge, mass);
    
    if (angles.theta_rms <= 0.0) {
        return angles;
    }
    
    if (use_moliere || thickness < 0.01) {
        angles = SampleMoliere(thickness, momentum, charge, mass, rng);
        angles.theta_rms = CalculateRMS(thickness, momentum, charge, mass);
    } else {
        angles = SampleGaussian(angles.theta_rms, rng);
    }
    
    return angles;
}

// ============================================================================
// Highland Formula for RMS Angle
// ============================================================================

double MultipleScattering::CalculateRMS(double thickness,
                                        double momentum,
                                        double charge,
                                        double mass) {
    if (thickness <= 0.0) {
        return 0.0;
    }
    
    double beta = CalculateBeta(momentum, mass);
    
    if (beta <= 0.0 || beta >= 1.0) {
        return 0.0;
    }
    
    double beta_p = beta * momentum;
    const double highland_const = 13.6;
    
    double theta_rms = (highland_const / beta_p) * std::abs(charge) * std::sqrt(thickness);
    
    if (thickness > 0.0) {
        double correction = 1.0 + 0.038 * std::log(thickness);
        theta_rms *= correction;
    }
    
    return theta_rms;
}

// ============================================================================
// Gaussian Approximation Sampling
// ============================================================================

MultipleScattering::ScatterAngles MultipleScattering::SampleGaussian(
    double theta_rms, 
    RandomGenerator& rng) {
    
    ScatterAngles angles;
    angles.theta_rms = theta_rms;
    
    if (theta_rms <= 0.0) {
        return angles;
    }
    
    angles.theta_x = rng.Gaussian(0.0, theta_rms);
    angles.theta_y = rng.Gaussian(0.0, theta_rms);
    
    return angles;
}

// ============================================================================
// Molière Theory Sampling
// ============================================================================

MultipleScattering::ScatterAngles MultipleScattering::SampleMoliere(
    double thickness,
    double momentum,
    double charge,
    double mass,
    RandomGenerator& rng) {
    
    ScatterAngles angles;
    
    double beta = CalculateBeta(momentum, mass);
    
    if (beta <= 0.0 || beta >= 1.0) {
        return angles;
    }
    
    double Z = std::round(0.5 * 1.0);
    
    double chi_alpha = CalculateChiAlpha(Z, beta);
    double B = CalculateMoliereB(chi_alpha, Z);
    double theta_space = SampleMoliereAngle(B, chi_alpha, rng);
    
    double phi = rng.Uniform(0.0, twopi);
    
    angles.theta_x = theta_space * std::cos(phi);
    angles.theta_y = theta_space * std::sin(phi);
    angles.theta_rms = CalculateRMS(thickness, momentum, charge, mass);
    
    return angles;
}

// ============================================================================
// Molière Theory Parameters
// ============================================================================

double MultipleScattering::CalculateChiAlpha(double Z, double beta) {
    if (beta <= 0.0) {
        return 0.0;
    }
    
    double Z_2_3 = std::pow(Z, 2.0/3.0);
    double beta2 = beta * beta;
    double Z_alpha = Z * alpha;
    double Z_alpha_beta = Z_alpha / beta;
    double correction = 1.0 + 3.34 * Z_alpha_beta * Z_alpha_beta;
    double chi_alpha = 2.007e-5 * Z_2_3 / (beta2 * correction);
    
    return chi_alpha;
}

double MultipleScattering::CalculateMoliereB(double chi_alpha, double Z) {
    if (chi_alpha <= 0.0) {
        return 0.0;
    }
    
    double Z_alpha = Z * alpha;
    double chi_c = 1.13 + 3.76 * Z_alpha * Z_alpha;
    double B = std::log(chi_c / chi_alpha);
    
    if (B < 0.0) {
        B = 0.0;
    }
    
    return B;
}

// ============================================================================
// Molière Angle Sampling (Internal)
// ============================================================================

double MultipleScattering::SampleMoliereAngle(double B, 
                                              double chi_alpha, 
                                              RandomGenerator& rng) {
    double theta_rms_sq = chi_alpha * B;
    double theta_rms = std::sqrt(std::max(theta_rms_sq, 0.0));
    
    double u = rng.Uniform();
    
    if (u < 0.95) {
        return std::abs(rng.Gaussian(0.0, theta_rms));
    } else {
        double tail_factor = 2.0 + rng.Exponential(1.0);
        return theta_rms * tail_factor;
    }
}

double MultipleScattering::MoliereCDF(double theta, double B, double chi_alpha) {
    double theta_rms = std::sqrt(chi_alpha * B);
    if (theta_rms <= 0.0) {
        return 0.0;
    }
    
    double u = theta / (theta_rms * std::sqrt(2.0));
    return 0.5 * (1.0 + std::erf(u));
}

std::array<double, 3> MultipleScattering::MoliereFunctions(double u, double B) {
    std::array<double, 3> f = {0.0, 0.0, 0.0};
    f[0] = std::exp(-u * u);
    f[1] = 0.0;
    f[2] = 0.0;
    return f;
}

// ============================================================================
// Utility Functions
// ============================================================================

double MultipleScattering::GetRadiationLength(const Material& material) {
    double X0_gcm2 = material.X0;
    double rho = material.density;
    return X0_gcm2 / rho;
}

double MultipleScattering::CalculateBeta(double momentum, double mass) {
    if (mass <= 0.0) {
        return 1.0;
    }
    
    double E = std::sqrt(momentum * momentum + mass * mass);
    return momentum / E;
}

double MultipleScattering::CalculateBetaGamma(double momentum, double mass) {
    if (mass <= 0.0) {
        return 1e10;
    }
    
    return momentum / mass;
}

// ============================================================================
// MultipleScatteringTable Implementation - CORRECTED
// ============================================================================

MultipleScatteringTable::MultipleScatteringTable(const Material& material,
                                                 double momentum_min,
                                                 double momentum_max,
                                                 int n_points)
    : material_(material),
      momentum_min_(momentum_min),
      momentum_max_(momentum_max) {
    
    momentum_points_.resize(n_points);
    rms_points_.resize(n_points);
    
    double log_pmin = std::log(momentum_min);
    double log_pmax = std::log(momentum_max);
    double dlog_p = (log_pmax - log_pmin) / (n_points - 1);
    
    // CORRECTED: Store theta/sqrt(x/X0) to make scaling work properly
    for (int i = 0; i < n_points; ++i) {
        double log_p = log_pmin + i * dlog_p;
        momentum_points_[i] = std::exp(log_p);
        
        // Calculate Highland formula for unit thickness (x/X0 = 1)
        // This gives us the "per sqrt(rad_length)" value
        double theta_per_sqrt_x = MultipleScattering::CalculateRMS(
            1.0, momentum_points_[i], 1.0, Me);
        
        rms_points_[i] = theta_per_sqrt_x;
    }
}

double MultipleScatteringTable::GetRMS(double momentum, 
                                       double length,
                                       double charge,
                                       double mass) const {
    
    // Get theta/sqrt(x/X0) from table
    double theta_per_sqrt_x = Interpolate(momentum);
    
    // Calculate actual thickness in radiation lengths
    double X0_cm = MultipleScattering::GetRadiationLength(material_);
    double thickness = length / X0_cm;
    
    // Scale by sqrt(thickness): theta = (theta/sqrt(x/X0)) * sqrt(x/X0)
    double theta_rms = theta_per_sqrt_x * std::sqrt(thickness);
    
    // Apply Lynch-Dahl correction for actual thickness
    if (thickness > 0.0) {
        double correction = 1.0 + 0.038 * std::log(thickness);
        theta_rms *= correction;
    }
    
    // Apply particle-specific correction (table built for electrons)
    double beta = MultipleScattering::CalculateBeta(momentum, mass);
    double beta_ref = MultipleScattering::CalculateBeta(momentum, Me);
    double particle_correction = (beta_ref / beta) * std::abs(charge);
    
    return theta_rms * particle_correction;
}

double MultipleScatteringTable::Interpolate(double momentum) const {
    if (momentum <= momentum_min_) {
        return rms_points_.front();
    }
    if (momentum >= momentum_max_) {
        return rms_points_.back();
    }
    
    double log_p = std::log(momentum);
    double log_pmin = std::log(momentum_min_);
    double log_pmax = std::log(momentum_max_);
    
    int n = momentum_points_.size();
    double dlog_p = (log_pmax - log_pmin) / (n - 1);
    
    int i = static_cast<int>((log_p - log_pmin) / dlog_p);
    i = std::max(0, std::min(i, n - 2));
    
    double log_p1 = std::log(momentum_points_[i]);
    double log_p2 = std::log(momentum_points_[i + 1]);
    
    double f = (log_p - log_p1) / (log_p2 - log_p1);
    
    return rms_points_[i] * (1.0 - f) + rms_points_[i + 1] * f;
}

} // namespace simc
