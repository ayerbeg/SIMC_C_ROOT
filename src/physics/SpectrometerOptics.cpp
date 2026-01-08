// src/physics/SpectrometerOptics.cpp
// Implementation of spectrometer optics and transport
// Ported from transp.f, mc_hms.f, mc_shms.f, mc_sos.f

#include "simc/SpectrometerOptics.h"
#include "simc/SimcConstants.h"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace simc {
using namespace constants;

// ============================================================================
// SpectrometerOptics Base Class Implementation
// ============================================================================

SpectrometerOptics::SpectrometerOptics(const std::string& name,
                                       const std::string& forward_matrix_file,
                                       const std::string& recon_matrix_file)
    : name_(name) {
    
    // Load COSY matrices
    forward_matrix_ = std::make_unique<CosyMatrix>();
    recon_matrix_ = std::make_unique<CosyMatrix>();
    
    if (!forward_matrix_file.empty()) {
        if (!forward_matrix_->LoadFromFile(forward_matrix_file)) {
            std::cerr << "Warning: Failed to load forward matrix: " 
                      << forward_matrix_file << std::endl;
        }
    }
    
    if (!recon_matrix_file.empty()) {
        if (!recon_matrix_->LoadFromFile(recon_matrix_file)) {
            std::cerr << "Warning: Failed to load reconstruction matrix: " 
                      << recon_matrix_file << std::endl;
        }
    }
}

// ============================================================================
// Forward Transport - Port from transp.f lines 200-500
// ============================================================================

SpectrometerOptics::TransportResult 
SpectrometerOptics::TransportForward(const ArmState& state,
                                     double momentum,
                                     double mass,
                                     double charge,
                                     RandomGenerator& rng,
                                     bool do_energy_loss,
                                     bool do_multiple_scattering) {
    TransportResult result;
    result.success = false;
    
    // Suppress unused parameter warning for charge (used in condition checks)
    (void)charge;  // Will be used when full transport implemented
    
    // Check if forward matrix is loaded
    if (!forward_matrix_ || !forward_matrix_->IsLoaded()) {
        result.stop_reason = "Forward matrix not loaded";
        return result;
    }
    
    // CRITICAL: Convert to COSY input vector
    // From transp.f lines 210-220:
    // Input vector is [x, xp, y, yp, delta] in cm, rad, rad, rad, PERCENT
    double input[5];
    input[0] = 0.0;           // x at target (cm) - usually 0
    input[1] = state.xptar;   // xp = dx/dz (rad)
    input[2] = 0.0;           // y at target (cm) - usually 0
    input[3] = state.yptar;   // yp = dy/dz (rad)
    input[4] = state.delta;   // delta in PERCENT
    
    // Apply forward COSY matrix
    // From transp.f lines 240-260
    double output[5];
    forward_matrix_->Apply(input, output);
    
    // Extract output: [x_fp, xp_fp, y_fp, yp_fp, delta_fp]
    double x_fp = output[0];    // cm
    double xp_fp = output[1];   // rad
    double y_fp = output[2];    // cm
    double yp_fp = output[3];   // rad
    double delta_fp = output[4]; // percent
    
    // Check apertures during transport
    // From mc_hms.f lines 300-450, mc_shms.f lines 250-400
    int num_planes = GetNumAperturePlanes();
    
    for (int iplane = 0; iplane < num_planes; ++iplane) {
        // For simplicity, we check at focal plane coordinates
        // In full version, would transport to each plane and check
        // This is a simplified version for initial implementation
        
        auto aperture_check = CheckAperture(x_fp, y_fp, iplane);
        
        if (!aperture_check.passed) {
            result.stop_plane = iplane;
            result.stop_reason = "Hit aperture: " + aperture_check.name;
            return result;
        }
    }
    
    // Apply multiple scattering if requested
    // From transp.f lines 400-450
    if (do_multiple_scattering && charge != 0.0) {
        // Approximate: apply MS at focal plane
        // Full version would apply at each drift section
        
        // Typical drift through air to focal plane ~ 5 meters
        double drift_length = 500.0; // cm
        Material air;
        air.name = "Air";
        air.Z = 7.5;  // Approximate for air (N2 + O2)
        air.A = 14.5;
        air.density = 0.001225; // g/cm^3
        
        auto ms_angles = MultipleScattering::Calculate(
            drift_length, air, momentum, mass, charge, rng);
        
        xp_fp += ms_angles.theta_x;
        yp_fp += ms_angles.theta_y;
    }
    
    // Apply energy loss if requested
    // From transp.f lines 470-500
    if (do_energy_loss && charge != 0.0) {
        // Energy loss through air to focal plane
        double drift_length = 500.0; // cm
        Material air;
        air.name = "Air";
        air.Z = 7.5;
        air.A = 14.5;
        air.density = 0.001225;
        
        double E_particle = std::sqrt(momentum*momentum + mass*mass);
        
        double Eloss = EnergyLoss::Calculate(
            drift_length, air.density, air.Z, air.A,
            E_particle, mass, rng,
            EnergyLoss::LossType::MOST_PROBABLE);
        
        // Update momentum and delta
        double p_new = std::sqrt((E_particle - Eloss)*(E_particle - Eloss) - mass*mass);
        delta_fp = (p_new - p_central_) / p_central_ * 100.0; // PERCENT!
    }
    
    // Calculate path length (approximate)
    // From transp.f lines 520-540
    // Path length ~ sqrt(1 + xp^2 + yp^2) * z_distance
    double path_factor = std::sqrt(1.0 + xp_fp*xp_fp + yp_fp*yp_fp);
    result.path_length = path_factor * 500.0; // Approximate focal plane distance
    
    // Fill focal plane state
    result.focal_plane.x = x_fp;
    result.focal_plane.y = y_fp;
    result.focal_plane.dx = xp_fp;
    result.focal_plane.dy = yp_fp;
    result.focal_plane.path = result.path_length;
    
    // Fill final arm state
    result.final_state.delta = delta_fp;
    result.final_state.xptar = state.xptar; // At target
    result.final_state.yptar = state.yptar;
    result.final_state.P = p_central_ * (1.0 + delta_fp / 100.0); // MeV/c
    result.final_state.E = std::sqrt(result.final_state.P * result.final_state.P + mass*mass);
    
    result.success = true;
    return result;
}

// ============================================================================
// Reconstruction - Port from transp.f lines 600-700
// ============================================================================

ArmState SpectrometerOptics::Reconstruct(const FocalPlaneState& focal_plane,
                                         double momentum) {
    ArmState reconstructed;
    
    // Suppress unused parameter warning
    (void)momentum;  // Used in future enhancement for momentum calculation
    
    // Check if reconstruction matrix is loaded
    if (!recon_matrix_ || !recon_matrix_->IsLoaded()) {
        std::cerr << "Warning: Reconstruction matrix not loaded" << std::endl;
        return reconstructed;
    }
    
    // CRITICAL: Convert focal plane to COSY input vector
    // From transp.f lines 610-630:
    // Input: [x_fp, xp_fp, y_fp, yp_fp, delta_fp]
    // Units: cm, rad, cm, rad, PERCENT
    
    double input[5];
    input[0] = focal_plane.x;   // cm
    input[1] = focal_plane.dx;  // rad (xp at focal plane)
    input[2] = focal_plane.y;   // cm
    input[3] = focal_plane.dy;  // rad (yp at focal plane)
    
    // Delta at focal plane - need to calculate from measured momentum
    // This would come from tracking or TOF in real data
    // For now, use a placeholder
    input[4] = 0.0; // Will be updated if momentum is known
    
    // Apply reconstruction matrix
    // From transp.f lines 650-670
    double output[5];
    recon_matrix_->Apply(input, output);
    
    // Extract target coordinates
    // Output: [x_tar, xp_tar, y_tar, yp_tar, delta_tar]
    reconstructed.xptar = output[1]; // xp at target (rad)
    reconstructed.yptar = output[3]; // yp at target (rad)
    reconstructed.delta = output[4]; // delta (PERCENT)
    
    // Calculate momentum and energy
    reconstructed.P = p_central_ * (1.0 + reconstructed.delta / 100.0);
    // Energy depends on particle type - need mass as input
    // For now, leave as momentum only
    
    return reconstructed;
}

// ============================================================================
// Drift Transport with MS and Energy Loss
// Port from transp.f lines 750-850
// ============================================================================

void SpectrometerOptics::TransportDrift(ArmState& state,
                                        double length,
                                        const Material* material,
                                        double momentum,
                                        double mass,
                                        double charge,
                                        RandomGenerator& rng,
                                        bool do_energy_loss,
                                        bool do_multiple_scattering) {
    // From transp.f lines 760-780
    // Simple drift: project forward by length
    // In full version, this would be broken into steps
    
    if (material == nullptr) {
        return; // Vacuum drift, no MS or dE/dx
    }
    
    // Apply multiple scattering
    // From transp.f lines 800-820
    if (do_multiple_scattering && charge != 0.0) {
        auto ms_angles = MultipleScattering::Calculate(
            length, *material, momentum, mass, charge, rng);
        
        // Add MS angles to particle trajectory
        state.xptar += ms_angles.theta_x;
        state.yptar += ms_angles.theta_y;
    }
    
    // Apply energy loss
    // From transp.f lines 830-850
    if (do_energy_loss && charge != 0.0) {
        double E_particle = std::sqrt(momentum*momentum + mass*mass);
        
        double Eloss = EnergyLoss::Calculate(
            length, material->density, material->Z, material->A,
            E_particle, mass, rng,
            EnergyLoss::LossType::MOST_PROBABLE);
        
        // Update momentum
        double E_new = E_particle - Eloss;
        double p_new = std::sqrt(E_new*E_new - mass*mass);
        
        // Update delta (PERCENT!)
        state.delta = (p_new - p_central_) / p_central_ * 100.0;
        state.P = p_new;
        state.E = E_new;
    }
}

// ============================================================================
// COSY Vector Conversions
// Port from transp.f lines 900-950
// ============================================================================

std::vector<double> SpectrometerOptics::StateToCosyVector(const ArmState& state) const {
    // From transp.f lines 910-920
    // CRITICAL: delta is in PERCENT!
    std::vector<double> vec(5);
    vec[0] = 0.0;          // x at target (cm) - usually 0
    vec[1] = state.xptar;  // xp (rad)
    vec[2] = 0.0;          // y at target (cm) - usually 0
    vec[3] = state.yptar;  // yp (rad)
    vec[4] = state.delta;  // delta in PERCENT
    return vec;
}

ArmState SpectrometerOptics::CosyVectorToState(const std::vector<double>& vec, 
                                                double momentum) const {
    // From transp.f lines 930-950
    ArmState state;
    state.xptar = vec[1];  // rad
    state.yptar = vec[3];  // rad
    state.delta = vec[4];  // PERCENT
    state.P = momentum * (1.0 + state.delta / 100.0); // MeV/c
    return state;
}

void SpectrometerOptics::ApplyCosyMatrix(const CosyMatrix& matrix, ArmState& state) {
    auto input = StateToCosyVector(state);
    double input_arr[5] = {input[0], input[1], input[2], input[3], input[4]};
    double output[5];
    
    matrix.Apply(input_arr, output);
    
    state = CosyVectorToState(
        std::vector<double>(output, output + 5), 
        p_central_);
}

// ============================================================================
// HMS Optics Implementation
// Port from mc_hms.f lines 100-600
// ============================================================================

HMSOptics::HMSOptics(const std::string& forward_file,
                     const std::string& recon_file)
    : SpectrometerOptics("HMS", forward_file, recon_file) {
}

SpectrometerOptics::ApertureCheck 
HMSOptics::CheckAperture(double x, double y, int plane) const {
    // From mc_hms.f lines 200-450 and apertures_hms.inc
    
    ApertureCheck check;
    check.passed = true;
    check.plane = plane;
    check.x = x;
    check.y = y;
    
    // HMS has 8 aperture planes (from apertures_hms.inc)
    switch (plane) {
        case 0: // Octagon entrance
            check.name = "Octagon entrance";
            check.passed = CheckOctagon(x, y);
            break;
            
        case 1: // Dipole entrance
            check.name = "Dipole entrance";
            check.passed = CheckDipole(x, y);
            break;
            
        case 2: // Dipole exit
            check.name = "Dipole exit";
            check.passed = CheckDipole(x, y);
            break;
            
        case 3: // Q1 entrance
            check.name = "Q1 entrance";
            check.passed = CheckQuad1(x, y);
            break;
            
        case 4: // Q1 exit
            check.name = "Q1 exit";
            check.passed = CheckQuad1(x, y);
            break;
            
        case 5: // Q2 aperture
            check.name = "Q2";
            check.passed = CheckQuad2(x, y);
            break;
            
        case 6: // Q3 aperture
            check.name = "Q3";
            check.passed = CheckQuad3(x, y);
            break;
            
        case 7: // Detector hut
            check.name = "Detector hut";
            // Detector hut is very large, rarely a constraint
            check.passed = (std::abs(x) < 100.0 && std::abs(y) < 100.0);
            break;
            
        default:
            check.passed = true;
    }
    
    return check;
}

bool HMSOptics::CheckOctagon(double x, double y) const {
    // From apertures_hms.inc lines 20-50
    // Octagonal aperture with 17.145 cm radius
    // From mc_hms.f: octagon has 8 sides
    
    const double r_oct = 17.145; // cm
    const double sqrt2 = 1.41421356;
    
    // Check rectangular bounds first (fast check)
    if (std::abs(x) > r_oct || std::abs(y) > r_oct) {
        return false;
    }
    
    // Check 45-degree corners
    // Octagon condition: |x| + |y| <= r * sqrt(2)
    if (std::abs(x) + std::abs(y) > r_oct * sqrt2) {
        return false;
    }
    
    return true;
}

bool HMSOptics::CheckDipole(double x, double y) const {
    // From apertures_hms.inc lines 60-80
    // Dipole: ±30 cm horizontal, ±12.5 cm vertical
    // From mc_hms.f lines 280-290
    
    const double x_max = 30.0;  // cm
    const double y_max = 12.5;  // cm
    
    return (std::abs(x) <= x_max && std::abs(y) <= y_max);
}

bool HMSOptics::CheckQuad1(double x, double y) const {
    // From apertures_hms.inc lines 100-110
    // Q1: circular aperture, r = 12.5 cm
    // From mc_hms.f lines 320-330
    
    const double r_q1 = 12.5; // cm
    double r2 = x*x + y*y;
    
    return (r2 <= r_q1 * r_q1);
}

bool HMSOptics::CheckQuad2(double x, double y) const {
    // From apertures_hms.inc lines 130-140
    // Q2: circular aperture, r = 30 cm
    // From mc_hms.f lines 350-360
    
    const double r_q2 = 30.0; // cm
    double r2 = x*x + y*y;
    
    return (r2 <= r_q2 * r_q2);
}

bool HMSOptics::CheckQuad3(double x, double y) const {
    // From apertures_hms.inc lines 160-170
    // Q3: circular aperture, r = 30 cm
    // From mc_hms.f lines 380-390
    
    const double r_q3 = 30.0; // cm
    double r2 = x*x + y*y;
    
    return (r2 <= r_q3 * r_q3);
}

// ============================================================================
// SHMS Optics Implementation
// Port from mc_shms.f lines 100-700
// ============================================================================

SHMSOptics::SHMSOptics(const std::string& forward_file,
                       const std::string& recon_file)
    : SpectrometerOptics("SHMS", forward_file, recon_file) {
}

SpectrometerOptics::ApertureCheck 
SHMSOptics::CheckAperture(double x, double y, int plane) const {
    // From mc_shms.f lines 250-500 and apertures_shms.inc
    
    ApertureCheck check;
    check.passed = true;
    check.plane = plane;
    check.x = x;
    check.y = y;
    
    // SHMS has 9 aperture planes
    switch (plane) {
        case 0: // Entrance aperture
            check.name = "Entrance";
            check.passed = CheckEntrance(x, y);
            break;
            
        case 1: // HB dipole entrance
            check.name = "HB dipole entrance";
            check.passed = CheckHBDipole(x, y);
            break;
            
        case 2: // HB dipole exit
            check.name = "HB dipole exit";
            check.passed = CheckHBDipole(x, y);
            break;
            
        case 3: // Q1 aperture
            check.name = "Q1";
            check.passed = (x*x + y*y <= 14.92*14.92); // r = 14.92 cm
            break;
            
        case 4: // Q2 aperture
            check.name = "Q2";
            check.passed = (x*x + y*y <= 30.0*30.0); // r = 30 cm
            break;
            
        case 5: // Q3 aperture
            check.name = "Q3";
            check.passed = (x*x + y*y <= 30.0*30.0); // r = 30 cm
            break;
            
        case 6: // Collimator
            check.name = "Collimator";
            check.passed = CheckCollimator(x, y);
            break;
            
        case 7: // Detector hut entrance
            check.name = "Hut entrance";
            check.passed = (std::abs(x) < 50.0 && std::abs(y) < 50.0);
            break;
            
        case 8: // Detector stack
            check.name = "Detector stack";
            check.passed = (std::abs(x) < 60.0 && std::abs(y) < 60.0);
            break;
            
        default:
            check.passed = true;
    }
    
    return check;
}

bool SHMSOptics::CheckEntrance(double x, double y) const {
    // From apertures_shms.inc lines 30-50
    // Entrance: ±60 cm horizontal, ±25 cm vertical
    // From mc_shms.f lines 280-290
    
    return (std::abs(x) <= 60.0 && std::abs(y) <= 25.0);
}

bool SHMSOptics::CheckHBDipole(double x, double y) const {
    // From apertures_shms.inc lines 70-90
    // HB Dipole: ±50 cm horizontal, ±25 cm vertical
    // From mc_shms.f lines 310-320
    
    return (std::abs(x) <= 50.0 && std::abs(y) <= 25.0);
}

bool SHMSOptics::CheckCollimator(double x, double y) const {
    // From apertures_shms.inc lines 150-180
    // Two collimator options: LARGE and SMALL
    // From mc_shms.f lines 380-410
    
    if (collimator_type_ == "LARGE") {
        // Large collimator: ±6 cm both directions
        return (std::abs(x) <= 6.0 && std::abs(y) <= 6.0);
    } else { // SMALL
        // Small collimator: ±3 cm both directions
        return (std::abs(x) <= 3.0 && std::abs(y) <= 3.0);
    }
}

void SHMSOptics::SetCollimator(const std::string& type) {
    // From mc_shms.f lines 150-160
    if (type == "LARGE" || type == "large") {
        collimator_type_ = "LARGE";
    } else if (type == "SMALL" || type == "small") {
        collimator_type_ = "SMALL";
    } else {
        std::cerr << "Warning: Unknown collimator type: " << type 
                  << ". Using LARGE." << std::endl;
        collimator_type_ = "LARGE";
    }
}

// ============================================================================
// SOS Optics Implementation
// Port from mc_sos.f lines 100-500
// ============================================================================

SOSOptics::SOSOptics(const std::string& forward_file,
                     const std::string& recon_file)
    : SpectrometerOptics("SOS", forward_file, recon_file) {
}

SpectrometerOptics::ApertureCheck 
SOSOptics::CheckAperture(double x, double y, int plane) const {
    // From mc_sos.f lines 200-400 and apertures_sos.inc
    
    ApertureCheck check;
    check.passed = true;
    check.plane = plane;
    check.x = x;
    check.y = y;
    
    // SOS has 6 aperture planes (simpler than HMS/SHMS)
    switch (plane) {
        case 0: // Entrance
            check.name = "Entrance";
            check.passed = (std::abs(x) < 30.0 && std::abs(y) < 15.0);
            break;
            
        case 1: // Dipole 1
            check.name = "Dipole 1";
            check.passed = (std::abs(x) < 25.0 && std::abs(y) < 12.0);
            break;
            
        case 2: // Q1
            check.name = "Q1";
            check.passed = (x*x + y*y <= 15.0*15.0);
            break;
            
        case 3: // Q2
            check.name = "Q2";
            check.passed = (x*x + y*y <= 20.0*20.0);
            break;
            
        case 4: // Dipole 2
            check.name = "Dipole 2";
            check.passed = (std::abs(x) < 25.0 && std::abs(y) < 12.0);
            break;
            
        case 5: // Detector hut
            check.name = "Detector hut";
            check.passed = (std::abs(x) < 40.0 && std::abs(y) < 30.0);
            break;
            
        default:
            check.passed = true;
    }
    
    return check;
}

// ============================================================================
// HRS Optics Implementation
// Port from mc_hrsl.f and mc_hrsr.f
// ============================================================================

HRSOptics::HRSOptics(const std::string& name,
                     const std::string& forward_file,
                     const std::string& recon_file)
    : SpectrometerOptics(name, forward_file, recon_file) {
}

SpectrometerOptics::ApertureCheck 
HRSOptics::CheckAperture(double x, double y, int plane) const {
    // From mc_hrsl.f / mc_hrsr.f
    // HRS apertures (both L and R have similar geometry)
    
    ApertureCheck check;
    check.passed = true;
    check.plane = plane;
    check.x = x;
    check.y = y;
    
    // HRS has 7 aperture planes
    switch (plane) {
        case 0: // Septum (if used)
            check.name = "Septum";
            check.passed = true; // Usually not a constraint
            break;
            
        case 1: // Q1 entrance
            check.name = "Q1 entrance";
            check.passed = (x*x + y*y <= 12.5*12.5);
            break;
            
        case 2: // Dipole entrance
            check.name = "Dipole entrance";
            check.passed = (std::abs(x) < 40.0 && std::abs(y) < 12.0);
            break;
            
        case 3: // Dipole exit
            check.name = "Dipole exit";
            check.passed = (std::abs(x) < 40.0 && std::abs(y) < 12.0);
            break;
            
        case 4: // Q2
            check.name = "Q2";
            check.passed = (x*x + y*y <= 30.0*30.0);
            break;
            
        case 5: // Q3
            check.name = "Q3";
            check.passed = (x*x + y*y <= 30.0*30.0);
            break;
            
        case 6: // Detector hut
            check.name = "Detector hut";
            check.passed = (std::abs(x) < 100.0 && std::abs(y) < 100.0);
            break;
            
        default:
            check.passed = true;
    }
    
    return check;
}

// ============================================================================
// Factory Functions
// ============================================================================

std::unique_ptr<SpectrometerOptics> CreateSpectrometerOptics(
    const std::string& type,
    const std::string& forward_file,
    const std::string& recon_file) {

    // Get default files if not specified
    auto [fwd, rec] = GetDefaultMatrixFiles(type);
    std::string fwd_file = forward_file.empty() ? fwd : forward_file;
    std::string rec_file = recon_file.empty() ? rec : recon_file;
    
    // Create appropriate spectrometer
    if (type == "HMS" || type == "hms") {
        return std::make_unique<HMSOptics>(fwd_file, rec_file);
    } else if (type == "SHMS" || type == "shms") {
        return std::make_unique<SHMSOptics>(fwd_file, rec_file);
    } else if (type == "SOS" || type == "sos") {
        return std::make_unique<SOSOptics>(fwd_file, rec_file);
    } else if (type == "HRSL" || type == "hrsl") {
        return std::make_unique<HRSOptics>("HRSL", fwd_file, rec_file);
    } else if (type == "HRSR" || type == "hrsr") {
        return std::make_unique<HRSOptics>("HRSR", fwd_file, rec_file);
    } else {
        std::cerr << "Error: Unknown spectrometer type: " << type << std::endl;
        return nullptr;
    }
}

std::pair<std::string, std::string> GetDefaultMatrixFiles(const std::string& type) {
    // Default matrix file paths (relative to data directory)
    
    if (type == "HMS" || type == "hms") {
        return {"hms/forward_cosy.dat", "hms/recon_cosy.dat"};
    } else if (type == "SHMS" || type == "shms") {
        return {"shms/shms_forward.dat", "shms/shms_recon.dat"};
    } else if (type == "SOS" || type == "sos") {
        return {"sos/forward_cosy.dat", "sos/recon_cosy.dat"};
    } else if (type == "HRSL" || type == "hrsl") {
        return {"hrsl/hrs_forward_cosy.dat", "hrsl/hrs_recon_cosy.dat"};
    } else if (type == "HRSR" || type == "hrsr") {
        return {"hrsr/hrs_forward_cosy.dat", "hrsr/hrs_recon_cosy.dat"};
    } else {
        return {"", ""};
    }
}

} // namespace simc
