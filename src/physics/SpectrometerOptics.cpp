// src/physics/SpectrometerOptics.cpp
// Implementation of spectrometer optics and transport
// Ported from transp.f, mc_hms.f, mc_shms.f, mc_sos.f
// CORRECTED VERSION - Fixed units (mrad vs rad) and sequential transport

#include "simc/physics/SpectrometerOptics.h"
#include "simc/core/SimcConstants.h"
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
    
    // Suppress unused parameter warnings
    (void)charge;
    (void)mass;
    (void)do_energy_loss;
    (void)do_multiple_scattering;
    (void)rng;
    
    // Check if forward matrix is loaded
    if (!forward_matrix_ || !forward_matrix_->IsLoaded()) {
        result.stop_reason = "Forward matrix not loaded";
        return result;
    }
    
    // CRITICAL UNIT CONVERSION (from transp.f lines 224-228):
    // Fortran uses COSY-7 units:
    //   Positions: cm
    //   Angles: mrad (milliradians), NOT radians!
    //   Delta: percent
    //
    // ArmState has angles in RADIANS, so we must convert:
    //   1 rad = 1000 mrad
    
    double input[5];
    input[0] = 0.0;                       // x at target (cm) - usually 0
    input[1] = state.xptar * 1000.0;      // Convert rad to mrad
    input[2] = 0.0;                       // y at target (cm) - usually 0
    input[3] = state.yptar * 1000.0;      // Convert rad to mrad
    input[4] = state.delta;               // delta in PERCENT (no conversion)
    
    // Apply forward COSY matrix
    // Output: [x_fp(cm), xp_fp(mrad), y_fp(cm), yp_fp(mrad), dL(cm)]
    double output[5];
    forward_matrix_->Apply(input, output);
    
    // Extract output
    double x_fp = output[0];              // cm
    double xp_fp = output[1];             // mrad
    double y_fp = output[2];              // cm
    double yp_fp = output[3];             // mrad
    double delta_L = output[4];           // cm (path length correction)
    
    // Check apertures during transport
    // NOTE: This is simplified - in full version would check at each plane
    // For now, check at focal plane
    int num_planes = GetNumAperturePlanes();
    
    for (int iplane = 0; iplane < num_planes; ++iplane) {
        auto aperture_check = CheckAperture(x_fp, y_fp, iplane);
        
        if (!aperture_check.passed) {
            result.stop_plane = iplane;
            result.stop_reason = "Hit aperture: " + aperture_check.name;
            return result;
        }
    }
    
    // Calculate path length
    // From transp.f line 420: pathlen = pathlen + (zd + delta_z)
    // where delta_z = -output[4]
    double nominal_path = 500.0; // Approximate focal plane distance (cm)
    result.path_length = nominal_path - delta_L;
    
    // Fill focal plane state
    // CONVERT back to radians for FocalPlaneState
    result.focal_plane.x = x_fp;
    result.focal_plane.y = y_fp;
    result.focal_plane.dx = xp_fp / 1000.0;  // Convert mrad to rad
    result.focal_plane.dy = yp_fp / 1000.0;  // Convert mrad to rad
    result.focal_plane.path = result.path_length;
    
    // Fill final arm state
    result.final_state.delta = state.delta;  // Delta unchanged by transport
    result.final_state.xptar = state.xptar;  // At target (input)
    result.final_state.yptar = state.yptar;  // At target (input)
    result.final_state.P = momentum * (1.0 + state.delta / 100.0); // MeV/c
    result.final_state.E = std::sqrt(result.final_state.P * result.final_state.P + mass*mass);
    
    result.success = true;
    return result;
}

// ============================================================================
// Reconstruction - Port from transp.f reconstruction subroutine
// ============================================================================

ArmState SpectrometerOptics::Reconstruct(const FocalPlaneState& focal_plane,
                                         double momentum) {
    ArmState reconstructed;
    
    // Check if reconstruction matrix is loaded
    if (!recon_matrix_ || !recon_matrix_->IsLoaded()) {
        std::cerr << "Warning: Reconstruction matrix not loaded" << std::endl;
        return reconstructed;
    }
    
    // CRITICAL UNIT CONVERSION:
    // Input focal plane coordinates are in radians
    // Must convert to COSY units (mrad) before applying matrix
    
    double input[5];
    input[0] = focal_plane.x;             // cm (no conversion)
    input[1] = focal_plane.dx * 1000.0;   // Convert rad to mrad
    input[2] = focal_plane.y;             // cm (no conversion)
    input[3] = focal_plane.dy * 1000.0;   // Convert rad to mrad
    input[4] = 0.0;  // Delta at focal plane (would come from tracking/TOF)
    
    // Apply reconstruction matrix
    double output[5];
    recon_matrix_->Apply(input, output);
    
    // Extract target coordinates
    // Output: [x_tar(cm), xp_tar(mrad), y_tar(cm), yp_tar(mrad), delta_tar(%)]
    // CONVERT angles back to radians
    reconstructed.xptar = output[1] / 1000.0;  // Convert mrad to rad
    reconstructed.yptar = output[3] / 1000.0;  // Convert mrad to rad
    reconstructed.delta = output[4];           // percent (no conversion)
    
    // Calculate momentum
    reconstructed.P = momentum * (1.0 + reconstructed.delta / 100.0);
    
    return reconstructed;
}

// ============================================================================
// Helper Functions
// ============================================================================

std::vector<double> SpectrometerOptics::StateToCosyVector(const ArmState& state) const {
    // Convert ArmState to COSY input vector
    // CRITICAL: Convert angles from radians to milliradians
    std::vector<double> vec(5);
    vec[0] = 0.0;                    // x at target (cm)
    vec[1] = state.xptar * 1000.0;   // xp (mrad)
    vec[2] = 0.0;                    // y at target (cm)
    vec[3] = state.yptar * 1000.0;   // yp (mrad)
    vec[4] = state.delta;            // delta (PERCENT)
    return vec;
}

ArmState SpectrometerOptics::CosyVectorToState(const std::vector<double>& vec, 
                                                double momentum) const {
    // Convert COSY output vector to ArmState
    // CRITICAL: Convert angles from milliradians to radians
    ArmState state;
    state.xptar = vec[1] / 1000.0;   // Convert mrad to rad
    state.yptar = vec[3] / 1000.0;   // Convert mrad to rad
    state.delta = vec[4];            // PERCENT (no conversion)
    state.P = momentum * (1.0 + state.delta / 100.0);
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

void SpectrometerOptics::TransportDrift(ArmState& state,
                                        double length,
                                        const Material* material,
                                        double momentum,
                                        double mass,
                                        double charge,
                                        RandomGenerator& rng,
                                        bool do_energy_loss,
                                        bool do_multiple_scattering) {
    // Simplified drift transport
    // Full version would be in sequential transport
    
    if (material == nullptr) {
        return; // Vacuum drift
    }
    
    // Apply multiple scattering
    if (do_multiple_scattering && charge != 0.0) {
        auto ms_angles = MultipleScattering::Calculate(
            length, *material, momentum, mass, charge, rng);
        
        state.xptar += ms_angles.theta_x;
        state.yptar += ms_angles.theta_y;
    }
    
    // Apply energy loss
    if (do_energy_loss && charge != 0.0) {
        double E_particle = std::sqrt(momentum*momentum + mass*mass);
        
        double Eloss = EnergyLoss::Calculate(
            length, material->density, material->Z, material->A,
            E_particle, mass, rng,
            EnergyLoss::LossType::MOST_PROBABLE);
        
        double E_new = E_particle - Eloss;
        double p_new = std::sqrt(E_new*E_new - mass*mass);
        
        state.delta = (p_new - p_central_) / p_central_ * 100.0;
        state.P = p_new;
        state.E = E_new;
    }
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
    
    const double r_oct = 17.145; // cm
    const double sqrt2 = 1.41421356;
    
    // Check rectangular bounds first (fast check)
    if (std::abs(x) > r_oct || std::abs(y) > r_oct) {
        return false;
    }
    
    // Check 45-degree corners
    if (std::abs(x) + std::abs(y) > r_oct * sqrt2) {
        return false;
    }
    
    return true;
}

bool HMSOptics::CheckDipole(double x, double y) const {
    // From apertures_hms.inc
    // Dipole: ±30 cm horizontal, ±12.5 cm vertical
    
    const double x_max = 30.0;  // cm
    const double y_max = 12.5;  // cm
    
    return (std::abs(x) <= x_max && std::abs(y) <= y_max);
}

bool HMSOptics::CheckQuad1(double x, double y) const {
    // Q1: circular aperture, r = 12.5 cm
    const double r_q1 = 12.5;
    return (x*x + y*y <= r_q1*r_q1);
}

bool HMSOptics::CheckQuad2(double x, double y) const {
    // Q2: circular aperture, r = 30 cm
    const double r_q2 = 30.0;
    return (x*x + y*y <= r_q2*r_q2);
}

bool HMSOptics::CheckQuad3(double x, double y) const {
    // Q3: circular aperture, r = 30 cm
    const double r_q3 = 30.0;
    return (x*x + y*y <= r_q3*r_q3);
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
    // From mc_shms.f and apertures_shms.inc
    
    ApertureCheck check;
    check.passed = true;
    check.plane = plane;
    check.x = x;
    check.y = y;
    
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
            check.passed = (x*x + y*y <= 14.92*14.92);
            break;
            
        case 4: // Q2 aperture
            check.name = "Q2";
            check.passed = (x*x + y*y <= 30.0*30.0);
            break;
            
        case 5: // Q3 aperture
            check.name = "Q3";
            check.passed = (x*x + y*y <= 30.0*30.0);
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
    // Entrance: ±60 cm horizontal, ±25 cm vertical
    return (std::abs(x) <= 60.0 && std::abs(y) <= 25.0);
}

bool SHMSOptics::CheckHBDipole(double x, double y) const {
    // HB Dipole: ±50 cm horizontal, ±25 cm vertical
    return (std::abs(x) <= 50.0 && std::abs(y) <= 25.0);
}

bool SHMSOptics::CheckCollimator(double x, double y) const {
    // Two collimator options: LARGE and SMALL
    if (collimator_type_ == "LARGE") {
        // Large collimator: ±6 cm both directions
        return (std::abs(x) <= 6.0 && std::abs(y) <= 6.0);
    } else { // SMALL
        // Small collimator: ±3 cm both directions
        return (std::abs(x) <= 3.0 && std::abs(y) <= 3.0);
    }
}

void SHMSOptics::SetCollimator(const std::string& type) {
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
// Port from mc_sos.f
// ============================================================================

SOSOptics::SOSOptics(const std::string& forward_file,
                     const std::string& recon_file)
    : SpectrometerOptics("SOS", forward_file, recon_file) {
}

SpectrometerOptics::ApertureCheck 
SOSOptics::CheckAperture(double x, double y, int plane) const {
    // From mc_sos.f and apertures_sos.inc
    
    ApertureCheck check;
    check.passed = true;
    check.plane = plane;
    check.x = x;
    check.y = y;
    
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
    // HRS apertures (both L and R have similar geometry)
    
    ApertureCheck check;
    check.passed = true;
    check.plane = plane;
    check.x = x;
    check.y = y;
    
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
    // Default matrix file paths (relative to data/matrices directory)
    
    if (type == "HMS" || type == "hms") {
        return {"data/matrices/hms/forward_cosy.dat", 
                "data/matrices/hms/recon_cosy.dat"};
    } else if (type == "SHMS" || type == "shms") {
        return {"data/matrices/shms/shms_forward.dat", 
                "data/matrices/shms/shms_recon.dat"};
    } else if (type == "SOS" || type == "sos") {
        return {"data/matrices/sos/forward_cosy.dat", 
                "data/matrices/sos/recon_cosy.dat"};
    } else if (type == "HRSL" || type == "hrsl") {
        return {"data/matrices/hrsl/hrs_forward_cosy.dat", 
                "data/matrices/hrsl/hrs_recon_cosy.dat"};
    } else if (type == "HRSR" || type == "hrsr") {
        return {"data/matrices/hrsr/hrs_forward_cosy.dat", 
                "data/matrices/hrsr/hrs_recon_cosy.dat"};
    } else {
        return {"", ""};
    }
}

} // namespace simc
