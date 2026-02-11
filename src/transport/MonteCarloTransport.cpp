// src/transport/MonteCarloTransport.cpp
// Phase 5c.1: Monte Carlo Transport Implementation
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
    
    // From target.f target_musc(), lines ~650-680
    // Uses Lynch & Dahl formula, NIM B58 (1991) 6-10
    
    constexpr double Es = 13.6;        // MeV (Highland formula constant)
    constexpr double epsilon = 0.088;  // Lynch & Dahl correction
    constexpr double nsig_max = 3.5;   // Maximum number of sigmas
    
    if (energy < 25.0) {
        std::cerr << "Warning: Energy passed to ApplyBeamMultipleScattering "
                  << "should be in MeV, but E = " << energy << std::endl;
    }
    
    if (teff <= 0.0) {
        dang_out[0] = 0.0;
        dang_out[1] = 0.0;
        return;
    }
    
    // RMS scattering angle (rad)
    // Lynch & Dahl form (better for beta != 1)
    double theta_sigma = (Es / energy / beta) * std::sqrt(teff) * 
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
    double momentum, double beta, double teff, double dangles[2]) {
    
    // Exact same formula as beam scattering
    // From shared/musc.f and target.f target_musc()
    
    constexpr double Es = 13.6;        // MeV
    constexpr double epsilon = 0.088;
    constexpr double nsig_max = 3.5;
    
    if (momentum < 25.0) {
        std::cerr << "Warning: Momentum passed to ApplyMultipleScattering "
                  << "should be in MeV/c, but p = " << momentum << std::endl;
    }
    
    if (teff <= 0.0) {
        dangles[0] = 0.0;
        dangles[1] = 0.0;
        return;
    }
    
    // RMS scattering angle (rad)
    double theta_sigma = (Es / momentum / beta) * std::sqrt(teff) * 
                        (1.0 + epsilon * std::log10(teff / (beta * beta)));
    
    // Generate two independent Gaussian scattering angles (truncated)
    auto truncated_gaussian = [this](double nsig_max) -> double {
        double value;
        do {
            value = random_->Gaussian(0.0, 1.0);
        } while (std::abs(value) > nsig_max);
        return value;
    };
    
    dangles[0] = theta_sigma * truncated_gaussian(nsig_max);
    dangles[1] = theta_sigma * truncated_gaussian(nsig_max);
}

// ============================================================================
// Apply Energy Loss
// ============================================================================

void MonteCarloTransport::ApplyEnergyLoss(SimcEvent& event, MainEvent& main) {
    // Energy loss is already calculated in EventGenerator::TripThruTarget()
    // and stored in main.target.Eloss[0,1,2]
    //
    // Here we just apply those corrections to get spectrometer quantities
    // This is already done in Transport() above
    
    // Suppress unused parameter warnings
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
// Check Acceptance (Phase 5c.1 - Simple Cuts)
// ============================================================================

bool MonteCarloTransport::CheckAcceptance(const SimcEvent& event, 
                                         const MainEvent& main) const {
    // Phase 5c.1: Simple delta and angle cuts
    // Phase 5c.2: Full spectrometer MC with apertures
    
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
    
    // Phase 5c.1: Simple cuts only
    // Phase 5c.2 will add:
    // - Aperture checking at each optical element
    // - Focal plane position checks
    // - Detector acceptance
    
    // Suppress unused parameter warning
    (void)event;
    
    return true;
}

// ============================================================================
// Phase 5c.2 Placeholders
// ============================================================================

bool MonteCarloTransport::TransportThroughSpectrometer(
    SimcEvent& event, MainEvent& main) {
    
    // Suppress unused parameter warnings
    (void)event;
    (void)main;
    
    // NOT IMPLEMENTED IN PHASE 5c.1
    //
    // Phase 5c.2 will implement:
    // - Call mc_hms, mc_sos, mc_shms, mc_hrsr, mc_hrsl
    // - Transport through magnetic elements
    // - Check apertures at each element
    // - Calculate focal plane quantities
    
    std::cerr << "ERROR: TransportThroughSpectrometer not implemented in Phase 5c.1\n";
    std::cerr << "       Full spectrometer MCs will be added in Phase 5c.2\n";
    
    return false;  // Always fail if called
}

void MonteCarloTransport::CalculateFocalPlane(const SimcEvent& event, 
                                             MainEvent& main) {
    // Suppress unused parameter warnings
    (void)event;
    (void)main;
    
    // NOT IMPLEMENTED IN PHASE 5c.1
    //
    // Phase 5c.2 will calculate:
    // - Focal plane position (xfp, yfp)
    // - Focal plane angles (dxfp, dyfp)
    // - Reconstruction back to target
    
    std::cerr << "ERROR: CalculateFocalPlane not implemented in Phase 5c.1\n";
    std::cerr << "       Full reconstruction will be added in Phase 5c.2\n";
}

} // namespace simc
