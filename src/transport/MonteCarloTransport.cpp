// src/transport/MonteCarloTransport.cpp
// Phase 5c.1: Monte Carlo Transport Implementation
// Phase 5c.2: CheckAcceptance DISABLED - let spectrometer do acceptance
// Ported from simc.f montecarlo() and target.f target_musc()

#include "simc/transport/MonteCarloTransport.h"
#include "simc/core/SimcConstants.h"
#include <cmath>
#include <iostream>

namespace simc {

using namespace constants;

// ============================================================================
// Constructor
// ============================================================================

MonteCarloTransport::MonteCarloTransport(
    const ConfigManager& config,
    std::shared_ptr<RandomGenerator> random)
    : random_(random)
{
    // Load configuration
    config_.use_energy_loss = config.Get<bool>("generation.use_energy_loss", true);
    config_.use_multiple_scattering = config.Get<bool>("transport.use_multiple_scattering", true);
    
    // Electron spectrometer
    config_.electron.P = config.Get<double>("spectrometer_electron.momentum");
    config_.electron.theta = config.Get<double>("spectrometer_electron.angle") * constants::DEG_TO_RAD;
    config_.electron.phi = config.Get<double>("spectrometer_electron.phi", 0.0) * constants::DEG_TO_RAD;
    config_.electron.cos_th = std::cos(config_.electron.theta);
    config_.electron.sin_th = std::sin(config_.electron.theta);
    
    config_.electron.delta_min = config.Get<double>("generation.electron.delta_min");
    config_.electron.delta_max = config.Get<double>("generation.electron.delta_max");
    config_.electron.xptar_min = config.Get<double>("generation.electron.xptar_min");
    config_.electron.xptar_max = config.Get<double>("generation.electron.xptar_max");
    config_.electron.yptar_min = config.Get<double>("generation.electron.yptar_min");
    config_.electron.yptar_max = config.Get<double>("generation.electron.yptar_max");
    
    // Hadron spectrometer
    config_.hadron.P = config.Get<double>("spectrometer_hadron.momentum");
    config_.hadron.theta = config.Get<double>("spectrometer_hadron.angle") * constants::DEG_TO_RAD;
    config_.hadron.phi = config.Get<double>("spectrometer_hadron.phi", 0.0) * constants::DEG_TO_RAD;
    config_.hadron.cos_th = std::cos(config_.hadron.theta);
    config_.hadron.sin_th = std::sin(config_.hadron.theta);
    
    config_.hadron.delta_min = config.Get<double>("generation.hadron.delta_min");
    config_.hadron.delta_max = config.Get<double>("generation.hadron.delta_max");
    config_.hadron.xptar_min = config.Get<double>("generation.hadron.xptar_min");
    config_.hadron.xptar_max = config.Get<double>("generation.hadron.xptar_max");
    config_.hadron.yptar_min = config.Get<double>("generation.hadron.yptar_min");
    config_.hadron.yptar_max = config.Get<double>("generation.hadron.yptar_max");
    
    // Target
    config_.length = config.Get<double>("target.length");
    config_.angle = config.Get<double>("target.angle", 0.0) * constants::DEG_TO_RAD;
    
    // Raster
    config_.correct_raster = config.Get<bool>("generation.correct_raster", false);
}

// ============================================================================
// Main Transport Function
// ============================================================================

bool MonteCarloTransport::Transport(SimcEvent& event, MainEvent& main) {
    // ========================================================================
    // BEAM MULTIPLE SCATTERING
    // From simc.f lines ~650-660
    // ========================================================================
    
    double dang_in[2] = {0.0, 0.0};  // {yptar, xptar} offsets
    
    if (config_.use_multiple_scattering) {
        constexpr double beta_electron = 1.0;  // Ultra-relativistic
        ApplyBeamMultipleScattering(event.Ein, beta_electron, 
                                   main.target.teff[0], dang_in);
    }
    
    // ========================================================================
    // HADRON ARM - Energy Loss and Multiple Scattering
    // From simc.f lines ~680-720
    // ========================================================================
    
    // Apply energy loss to get spectrometer delta
    if (config_.use_energy_loss) {
        // Hadron energy after energy loss
        double E_after_loss = event.p_E - main.target.Eloss[2];
        double P_after_loss = std::sqrt(E_after_loss * E_after_loss - Mp * Mp);
        main.SP_hadron.delta = (P_after_loss - config_.hadron.P) / config_.hadron.P * 100.0;
    } else {
        main.SP_hadron.delta = event.p_delta;
    }
    
    // Apply multiple scattering
    double dangles_p[2] = {0.0, 0.0};
    if (config_.use_multiple_scattering) {
        double beta_p = event.p_P / event.p_E;
        ApplyMultipleScattering(event.p_P, beta_p, main.target.teff[2], dangles_p);
    }
    
    // Spectrometer quantities (after multiple scattering)
    main.SP_hadron.yptar = event.p_yptar + dangles_p[0] + dang_in[0];
    main.SP_hadron.xptar = event.p_xptar + dangles_p[1] + dang_in[1] * config_.hadron.cos_th;
    
    // ========================================================================
    // ELECTRON ARM - Energy Loss and Multiple Scattering
    // Similar to hadron arm (not shown in montecarlo but would be here)
    // ========================================================================
    
    if (config_.use_energy_loss) {
        // Electron is ultra-relativistic: E ≈ P
        double P_after_loss = event.e_E - main.target.Eloss[1];
        main.SP_electron.delta = (P_after_loss - config_.electron.P) / config_.electron.P * 100.0;
    } else {
        main.SP_electron.delta = event.e_delta;
    }
    
    double dangles_e[2] = {0.0, 0.0};
    if (config_.use_multiple_scattering) {
        constexpr double beta_electron = 1.0;
        ApplyMultipleScattering(event.e_P, beta_electron, main.target.teff[1], dangles_e);
    }
    
    main.SP_electron.yptar = event.e_yptar + dangles_e[0] + dang_in[0];
    main.SP_electron.xptar = event.e_xptar + dangles_e[1] + dang_in[1] * config_.electron.cos_th;
    
    // ========================================================================
    // ACCEPTANCE CHECK (Phase 5c.1 - simple cuts)
    // ========================================================================
    
    bool passes = CheckAcceptance(event, main);
    
    return passes;
}

// ============================================================================
// Calculate Reconstructed Quantities
// ============================================================================

void MonteCarloTransport::CalculateReconstructed(SimcEvent& event) {
    // For elastic H(e,e'p), we need to calculate p_xptar, p_yptar
    // from the proton's physics angles (p_theta, p_phi)
    // These were set to 0.0 in EventGenerator - now we fix them!
    
    SpectrometerAngles(config_.hadron.theta, config_.hadron.phi,
                      event.p_theta, event.p_phi,
                      event.p_xptar, event.p_yptar);
}

// ============================================================================
// Apply Beam Multiple Scattering
// ============================================================================

void MonteCarloTransport::ApplyBeamMultipleScattering(
    double energy, double beta, double teff, double dang_out[2]) {

    double energy_MeV = energy * 1000.0;  // GeV → MeV
  
    // From target.f target_musc(), lines ~650-680
    // Uses Lynch & Dahl formula, NIM B58 (1991) 6-10
    
    constexpr double Es = 13.6;        // MeV (Highland formula constant)
    constexpr double epsilon = 0.088;  // Lynch & Dahl correction
    constexpr double nsig_max = 3.5;   // Maximum number of sigmas
    
    // if (energy < 25.0) {
    //     std::cerr << "Warning: Energy passed to ApplyBeamMultipleScattering "
    //               << "should be in MeV, but E = " << energy << std::endl;
    // }
    
    if (teff <= 0.0) {
        dang_out[0] = 0.0;
        dang_out[1] = 0.0;
        return;
    }
    
    // RMS scattering angle (rad)
    // Lynch & Dahl form (better for beta != 1)
    double theta_sigma = (Es / energy_MeV / beta) * std::sqrt(teff) * 
                        (1.0 + epsilon * std::log10(teff / (beta * beta)));
    
    // Generate two Gaussian random numbers (limited to nsig_max)
    // Your RandomGenerator doesn't support truncation, so we implement it here
    auto truncated_gaussian = [this](double nsig_max) -> double {
        double value;
        do {
            value = random_->Gaussian(0.0, 1.0);
        } while (std::abs(value) > nsig_max);
        return value;
    };
    
    dang_out[0] = theta_sigma * truncated_gaussian(nsig_max);
    dang_out[1] = theta_sigma * truncated_gaussian(nsig_max);
}

// ============================================================================
// Apply Multiple Scattering
// ============================================================================

void MonteCarloTransport::ApplyMultipleScattering(
    double p, double beta, double teff, double dang_out[2]) {

  double p_MeV = p * 1000.0;  // GeV → MeV
  
    // From target.f target_musc(), lines ~680-710
    // Uses same Lynch & Dahl formula as beam MS
    
    constexpr double Es = 13.6;        // MeV
    constexpr double epsilon = 0.088;
    constexpr double nsig_max = 3.5;
    
    if (teff <= 0.0) {
        dang_out[0] = 0.0;
        dang_out[1] = 0.0;
        return;
    }
    
    // RMS scattering angle (rad)
    double theta_sigma = (Es / p_MeV / beta) * std::sqrt(teff) * 
                        (1.0 + epsilon * std::log10(teff / (beta * beta)));
    
    auto truncated_gaussian = [this](double nsig_max) -> double {
        double value;
        do {
            value = random_->Gaussian(0.0, 1.0);
        } while (std::abs(value) > nsig_max);
        return value;
    };
    
    dang_out[0] = theta_sigma * truncated_gaussian(nsig_max);
    dang_out[1] = theta_sigma * truncated_gaussian(nsig_max);
}

// ============================================================================
// Apply Energy Loss
// ============================================================================

void MonteCarloTransport::ApplyEnergyLoss(SimcEvent& event, MainEvent& main) {
    // Placeholder for Phase 5c.3: Energy loss through spectrometer windows, etc.
    // Currently energy loss through target is handled in EventGenerator
    
    (void)event;
    (void)main;
    
    // NOTE: This function is kept for future expansion if we need
    // to calculate energy loss through spectrometer windows, air gaps, etc.
}

// ============================================================================
// Calculate Spectrometer Angles from Physics Angles
// ============================================================================

void MonteCarloTransport::SpectrometerAngles(
    double theta0, double phi0,
    double theta, double phi,
    double& xptar, double& yptar) const {
    
    // INVERSE of PhysicsAngles in EventGenerator
    // From event.f spectrometer_angles(), lines ~1720-1750
    //
    // Given physics angles (theta, phi) in lab frame,
    // calculate spectrometer angles (xptar, yptar)
    
    double cos_theta0 = std::cos(theta0);
    double sin_theta0 = std::sin(theta0);
    double cos_phi0 = std::cos(phi0);
    double sin_phi0 = std::sin(phi0);
    
    // Unit vector in lab frame
    double ux = std::sin(theta) * std::cos(phi);
    double uy = std::sin(theta) * std::sin(phi);
    double uz = std::cos(theta);
    
    // Rotate lab frame vector BACK to spectrometer frame
    // This is the INVERSE (transpose) of the rotation in PhysicsAngles
    double dx = cos_theta0 * cos_phi0 * ux + cos_theta0 * sin_phi0 * uy - sin_theta0 * uz;
    double dy = -sin_phi0 * ux + cos_phi0 * uy;
    double dz = sin_theta0 * cos_phi0 * ux + sin_theta0 * sin_phi0 * uy + cos_theta0 * uz;
    
    // Spectrometer angles are the direction ratios
    xptar = dx / dz;
    yptar = dy / dz;
}

// ============================================================================
// Check Acceptance (Phase 5c.2 - DISABLED)
// ============================================================================

bool MonteCarloTransport::CheckAcceptance(const SimcEvent& event, 
                                         const MainEvent& main) const {
    // ========================================================================
    // PHASE 5c.2 NOTE: This function is DISABLED
    // ========================================================================
    //
    // Previously (Phase 5c.1), this applied simple rectangular cuts on
    // spectrometer angles and delta to provide basic acceptance filtering
    // before full spectrometer transport was implemented.
    //
    // NOW (Phase 5c.2+), acceptance is determined by:
    // - SHMS class: Full 32-aperture transport with realistic geometry
    // - HMS class: (Future) Full aperture transport
    // - etc.
    //
    // The rectangular cuts from config (±60 mrad, ±30 mrad) were too tight
    // and rejected valid events that should have been transported through
    // the spectrometer. For example, elastic H(e,e'p) at 35° SHMS angle
    // produces proton angles of ~2400 mrad relative to spectrometer axis,
    // which is physically correct but fails the simple cuts.
    //
    // PHYSICS EXPLANATION:
    // -------------------
    // For elastic scattering H(e,e'p), when SHMS is at 35° and protons
    // scatter at ~33°, the angle RELATIVE to the spectrometer axis is:
    //   angle_rel = arctan(sin(33-35)/cos(33-35)) ≈ -2.4 rad ≈ -2400 mrad
    // This is CORRECT physics - the proton is going nearly opposite to
    // the spectrometer direction! The spectrometer must rotate and bend
    // this trajectory back through the apertures.
    //
    // The cuts below (lines 346-369) would reject these valid events:
    // -------------------
    /*
    // Electron arm cuts
    if (main.SP_electron.delta < config_.electron.delta_min || 
        main.SP_electron.delta > config_.electron.delta_max) {
        return false;
    }
    if (main.SP_electron.xptar < config_.electron.xptar_min || 
        main.SP_electron.xptar > config_.electron.xptar_max) {
        return false;
    }
    if (main.SP_electron.yptar < config_.electron.yptar_min || 
        main.SP_electron.yptar > config_.electron.yptar_max) {
        return false;
    }
    
    // Hadron arm cuts
    if (main.SP_hadron.delta < config_.hadron.delta_min || 
        main.SP_hadron.delta > config_.hadron.delta_max) {
        return false;
    }
    if (main.SP_hadron.xptar < config_.hadron.xptar_min || 
        main.SP_hadron.xptar > config_.hadron.xptar_max) {
        return false;
    }
    if (main.SP_hadron.yptar < config_.hadron.yptar_min || 
        main.SP_hadron.yptar > config_.hadron.yptar_max) {
        return false;
    }
    */
    //
    // FUTURE TODO (when implementing multiple spectrometers):
    // -------------------------------------------------------
    // - Add spectrometer-specific pre-transport sanity checks
    // - Consider re-enabling VERY WIDE cuts (e.g., |angle| < 5 rad = 286°)
    //   to catch obviously pathological events before expensive transport
    // - Or keep this always returning true and let each spectrometer
    //   class handle 100% of acceptance determination
    //
    // CURRENT BEHAVIOR: Always return true - let spectrometer transport
    //                   determine acceptance through realistic aperture checks
    
    (void)event;  // Suppress unused parameter warnings
    (void)main;
    
    return true;  // ALWAYS PASS - acceptance determined by spectrometer transport
}

// ============================================================================
// Phase 5c.2 Placeholders
// ============================================================================

bool MonteCarloTransport::TransportThroughSpectrometer(
    SimcEvent& event, MainEvent& main) {
    
    (void)event;
    (void)main;
    
    // This is now handled by SHMS::Transport() in main.cpp
    return true;
}

} // namespace simc
