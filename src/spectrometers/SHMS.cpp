#include "simc/spectrometers/SHMS.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace simc {

namespace {
    constexpr double DEG_TO_RAD = M_PI / 180.0;
    constexpr double RAD_TO_DEG = 180.0 / M_PI;
}

SHMS::SHMS() {
    // Initialize with default values (already done via member initializers)
}

bool SHMS::Transport(TrackState& track) {
    stats_.total_events++;
    
    // Transport through each element in sequence
    if (!TransportHoldingBox(track)) return false;
    if (!TransportQ1(track)) return false;
    if (!TransportQ2(track)) return false;
    if (!TransportQ3(track)) return false;
    if (!TransportDipole(track)) return false;
    if (!TransportHut(track)) return false;
    
    stats_.accepted++;
    return true;
}

bool SHMS::Reconstruct(const FocalPlaneState& fp, TargetState& target) {
    // TODO: Implement reconstruction in Session 4
    // For now, just copy values
    target.x = fp.x;
    target.xp = fp.xp;
    target.y = fp.y;
    target.yp = fp.yp;
    target.delta = fp.delta;
    return true;
}
// ============================================================================
// PHASE 5e: Add this implementation to src/spectrometers/SHMS.cpp
// Insert after the Reconstruct() method (around line 44)
// ============================================================================

bool SHMS::GetFocalPlane(const TrackState& track, FocalPlaneState& fp) const {
    /**
     * Extract focal plane coordinates from transported track state.
     * 
     * After Transport() completes successfully, the TrackState contains
     * the particle coordinates at the focal plane detector. This method
     * extracts those values into a FocalPlaneState structure.
     * 
     * IMPORTANT: This should ONLY be called AFTER Transport() has succeeded!
     * 
     * UNIT CONVERSIONS:
     * - track.x, track.y are in cm → fp.x, fp.y in cm (no conversion)
     * - track.dx, track.dy are slopes (dimensionless) → fp.xp, fp.yp in mrad
     * - track.delta is in % → fp.delta in % (no conversion)
     * 
     * COSY Convention:
     * - Slopes are defined as dx/dz and dy/dz
     * - For small angles: angle (mrad) ≈ 1000 * slope
     * 
     * Phase 5e - Week 1, Day 1
     */
    
    // Copy position coordinates (already in cm)
    fp.x = track.x;
    fp.y = track.y;
    
    // Convert slopes to mrad
    // track.dx and track.dy are dimensionless slopes (dx/dz, dy/dz)
    // fp.xp and fp.yp are in mrad
    fp.xp = track.dx * 1000.0;  // slope → mrad
    fp.yp = track.dy * 1000.0;  // slope → mrad
    
    // Copy momentum deviation (already in %)
    fp.delta = track.delta;
    
    return true;
}




  
bool SHMS::LoadMatrices(const std::string& forward_file, const std::string& recon_file) {
    if (!ParseMatrixFile(forward_file, forward_matrix_)) {
        std::cerr << "Failed to load forward matrix: " << forward_file << std::endl;
        return false;
    }
    
    if (!ParseMatrixFile(recon_file, recon_matrix_)) {
        std::cerr << "Failed to load reconstruction matrix: " << recon_file << std::endl;
        return false;
    }
    
    return true;
}

// ============================================================================
// Helper Functions
// ============================================================================

void SHMS::Project(TrackState& track, double zdrift) {
    /**
     * Calculate new track transverse coordinates after drifting in a field-free
     * region for distance zdrift from current track location.
     * 
     * Based on shared/project.f
     */
    
    // Calculate path length including slope correction
    double path_correction = std::sqrt(1.0 + track.dx*track.dx + track.dy*track.dy);
    track.pathlen += zdrift * path_correction;
    
    // Update transverse positions
    track.x += track.dx * zdrift;
    track.y += track.dy * zdrift;
    track.z += zdrift;
    
    // Note: Decay handling would go here if needed (not implemented yet)
}

void SHMS::TranspMatrix(TrackState& track, int class_id, double zdrift) {
    /**
     * Transport particle through spectrometer element using matrix elements.
     * 
     * Based on shared/transp.f
     * 
     * The matrix transformation uses COSY polynomial format:
     * For each output variable (X, XP, Y, YP, dL), sum over all terms:
     *   output = sum(coeff * x^ex * xp^exp * y^ey * yp^eyp * delta^ed)
     */
    
    if (class_id < 0 || class_id >= MatrixElements::MAX_CLASSES) {
        std::cerr << "Invalid matrix class ID: " << class_id << std::endl;
        return;
    }
    
    const auto& transformation = forward_matrix_.classes[class_id];
    
    if (transformation.n_terms == 0) {
        std::cerr << "Matrix class " << class_id << " has no terms!" << std::endl;
        return;
    }
    
    // Pack input coordinates (convert to COSY units)
    std::array<double, 5> ray;
    ray[0] = track.x;              // cm
    ray[1] = track.dx * 1000.0;    // mrad
    ray[2] = track.y;              // cm
    ray[3] = track.dy * 1000.0;    // mrad
    ray[4] = track.delta;          // %
    
    // Compute COSY sums
    std::array<double, 5> sum{0.0};
    
    for (const auto& term : transformation.terms) {
        double product = 1.0;
        
        // Calculate term = x^ex * xp^exp * y^ey * yp^eyp * delta^ed
        for (int j = 0; j < 5; ++j) {
            if (term.expon[j] != 0) {
                product *= std::pow(ray[j], term.expon[j]);
            }
        }
        
        // Add contribution to each output variable
        for (int j = 0; j < 5; ++j) {
            sum[j] += product * term.coeff[j];
        }
    }
    
    // Unpack output coordinates
    track.x = sum[0];              // cm
    track.dx = sum[1] / 1000.0;    // mrad -> slope
    track.y = sum[2];              // cm
    track.dy = sum[3] / 1000.0;    // mrad -> slope
    double delta_z = -sum[4];      // Path length correction (cm)
    
    // Update path length
    track.pathlen += (zdrift + delta_z);
    track.z += (zdrift + delta_z);
    
    // Note: Decay handling would go here if needed
}

void SHMS::RotateVAxis(double angle_deg, double& x, double& y, double dx, double dy) {
    /**
     * Calculate new trajectory coordinates in reference frame rotated about
     * vertical axis by angle_deg relative to central ray.
     * 
     * Based on shared/rotate_vaxis.f
     * 
     * Right-handed TRANSPORT coordinates:
     * - Rotation about negative X-axis
     */
    
    double angle_rad = angle_deg * DEG_TO_RAD;
    double tan_th = std::tan(angle_rad);
    double sin_th = std::sin(angle_rad);
    double cos_th = std::cos(angle_rad);
    
    double alpha = dy;    // dydzs
    double beta = dx;     // dxdzs
    
    double alpha_p = (alpha + tan_th) / (1.0 - alpha * tan_th);
    double beta_p = beta / (cos_th - alpha * sin_th);
    
    double y_initial = y;
    y = y_initial * (cos_th + alpha_p * sin_th);
    x = x + y_initial * beta_p * sin_th;
}

void SHMS::RotateHAxis(double angle_deg, double& x, double& y, double dx, double dy) {
    /**
     * Calculate new trajectory coordinates in reference frame rotated about
     * horizontal axis by angle_deg relative to central ray.
     * 
     * Based on shared/rotate_haxis.f
     * 
     * Right-handed TRANSPORT coordinates:
     * - Rotation about negative Y-axis
     */
    
    double angle_rad = angle_deg * DEG_TO_RAD;
    double tan_th = std::tan(angle_rad);
    double sin_th = std::sin(angle_rad);
    double cos_th = std::cos(angle_rad);
    
    double alpha = dx;    // dxdzs
    double beta = dy;     // dydzs
    
    double alpha_p = (alpha + tan_th) / (1.0 - alpha * tan_th);
    double beta_p = beta / (cos_th - alpha * sin_th);
    
    double x_initial = x;
    x = x_initial * (cos_th + alpha_p * sin_th);
    y = y + x_initial * beta_p * sin_th;
}

// ============================================================================
// Aperture Checking
// ============================================================================

bool SHMS::CheckHBAperture(double x, double y, int aperture_id) {
    /**
     * Check if particle passes through Holding Box aperture.
     * HB apertures are asymmetric in Y and tilted.
     */
    
    if (aperture_id < 0 || aperture_id >= 4) {
        return false;
    }
    
    // Check vertical (X) limit
    if (x*x > apertures_.r_HBx * apertures_.r_HBx) {
        return false;
    }
    
    // Check horizontal (Y) limits (asymmetric)
    if (y > apertures_.r_HBfyp[aperture_id] || y < apertures_.r_HBfym[aperture_id]) {
        return false;
    }
    
    return true;
}

bool SHMS::CheckQuadAperture(double x, double y, double radius) {
    /**
     * Check if particle passes through circular quadrupole aperture.
     */
    return (x*x + y*y) <= (radius * radius);
}

bool SHMS::CheckDipoleAperture(double x, double y, int section) {
    /**
     * Check if particle passes through dipole aperture.
     * Dipole apertures are circular but may be tilted and offset.
     */
    
    if (section < 0 || section >= 12) {
        return false;
    }
    
    // Check circular aperture
    return (x*x + y*y) <= (apertures_.r_D1 * apertures_.r_D1);
}

bool SHMS::CheckCollimatorAperture(double x, double y, bool entrance) {
    /**
     * Check if particle passes through collimator aperture (rectangular).
     */
    
    double x_limit = entrance ? apertures_.coll_x_in : apertures_.coll_x_out;
    double y_limit = entrance ? apertures_.coll_y_in : apertures_.coll_y_out;
    
    return (std::abs(x) <= x_limit) && (std::abs(y) <= y_limit);
}

// ============================================================================
// Element Transport Functions
// ============================================================================

bool SHMS::TransportHoldingBox(TrackState& track) {
    /**
     * Transport through Holding Box (Bender) with 4 aperture checks.
     * 
     * Based on mc_shms.f lines ~180-250
     * 
     * Sequence:
     * 1. Drift to HB entrance, check aperture (tilted, offset)
     * 2. Drift to HB mag entrance, check aperture (tilted, offset)
     * 3. MATRIX transport through HB magnet, check aperture (tilted, offset)
     * 4. Drift to HB exit, check aperture (tilted, offset)
     */
    
    // 1. HB Entrance
    Project(track, drifts_.hb_in);
    
    double xt = track.x;
    double yt = track.y;
    RotateVAxis(apertures_.hb_tilt[0], xt, yt, track.dx, track.dy);
    yt += apertures_.hb_offset[0];
    
    if (!CheckHBAperture(xt, yt, 0)) {
        stats_.stop_hb_in++;
        return false;
    }
    
    // 2. HB Magnetic Entrance
    Project(track, drifts_.hb_men);
    
    xt = track.x;
    yt = track.y;
    RotateVAxis(apertures_.hb_tilt[1], xt, yt, track.dx, track.dy);
    yt += apertures_.hb_offset[1];
    
    if (!CheckHBAperture(xt, yt, 1)) {
        stats_.stop_hb_men++;
        return false;
    }
    
    // 3. HB Magnetic Exit (uses MATRIX transport - class 2)
    TranspMatrix(track, 2, drifts_.hb_mex);
    
    xt = track.x;
    yt = track.y;
    RotateVAxis(apertures_.hb_tilt[2], xt, yt, track.dx, track.dy);
    yt += apertures_.hb_offset[2];
    
    if (!CheckHBAperture(xt, yt, 2)) {
        stats_.stop_hb_mex++;
        return false;
    }
    
    // 4. HB Exit
    Project(track, drifts_.hb_out);
    
    xt = track.x;
    yt = track.y;
    RotateVAxis(apertures_.hb_tilt[3], xt, yt, track.dx, track.dy);
    yt += apertures_.hb_offset[3];
    
    if (!CheckHBAperture(xt, yt, 3)) {
        stats_.stop_hb_out++;
        return false;
    }
    
    return true;
}

bool SHMS::TransportQ1(TrackState& track) {
    /**
     * Transport through Q1 (Quadrupole 1) with 5 aperture checks.
     * 
     * Based on mc_shms.f lines ~250-320
     */
    
    // 1. Q1 Entrance
    Project(track, drifts_.q1_in);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q1)) {
        stats_.stop_q1_in++;
        return false;
    }
    
    // 2. Q1 Magnetic Entrance
    Project(track, drifts_.q1_men);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q1)) {
        stats_.stop_q1_men++;
        return false;
    }
    
    // 3. Q1 Mid (uses MATRIX transport - class 6)
    TranspMatrix(track, 6, drifts_.q1_mid);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q1)) {
        stats_.stop_q1_mid++;
        return false;
    }
    
    // 4. Q1 Magnetic Exit (uses MATRIX transport - class 7)
    TranspMatrix(track, 7, drifts_.q1_mex);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q1)) {
        stats_.stop_q1_mex++;
        return false;
    }
    
    // 5. Q1 Exit
    Project(track, drifts_.q1_out);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q1)) {
        stats_.stop_q1_out++;
        return false;
    }
    
    return true;
}

bool SHMS::TransportQ2(TrackState& track) {
    /**
     * Transport through Q2 (Quadrupole 2) with 5 aperture checks.
     */
    
    // 1. Q2 Entrance
    Project(track, drifts_.q2_in);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q2)) {
        stats_.stop_q2_in++;
        return false;
    }
    
    // 2. Q2 Magnetic Entrance
    Project(track, drifts_.q2_men);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q2)) {
        stats_.stop_q2_men++;
        return false;
    }
    
    // 3. Q2 Mid (uses MATRIX transport - class 10)
    TranspMatrix(track, 10, drifts_.q2_mid);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q2)) {
        stats_.stop_q2_mid++;
        return false;
    }
    
    // 4. Q2 Magnetic Exit (uses MATRIX transport - class 11)
    TranspMatrix(track, 11, drifts_.q2_mex);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q2)) {
        stats_.stop_q2_mex++;
        return false;
    }
    
    // 5. Q2 Exit
    Project(track, drifts_.q2_out);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q2)) {
        stats_.stop_q2_out++;
        return false;
    }
    
    return true;
}

bool SHMS::TransportQ3(TrackState& track) {
    /**
     * Transport through Q3 (Quadrupole 3) with 5 aperture checks.
     */
    
    // 1. Q3 Entrance
    Project(track, drifts_.q3_in);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q3)) {
        stats_.stop_q3_in++;
        return false;
    }
    
    // 2. Q3 Magnetic Entrance
    Project(track, drifts_.q3_men);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q3)) {
        stats_.stop_q3_men++;
        return false;
    }
    
    // 3. Q3 Mid (uses MATRIX transport - class 14)
    TranspMatrix(track, 14, drifts_.q3_mid);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q3)) {
        stats_.stop_q3_mid++;
        return false;
    }
    
    // 4. Q3 Magnetic Exit (uses MATRIX transport - class 15)
    TranspMatrix(track, 15, drifts_.q3_mex);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q3)) {
        stats_.stop_q3_mex++;
        return false;
    }
    
    // 5. Q3 Exit
    Project(track, drifts_.q3_out);
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q3)) {
        stats_.stop_q3_out++;
        return false;
    }
    
    return true;
}

bool SHMS::TransportDipole(TrackState& track) {
    /**
     * Transport through Dipole (D1) with 12 aperture checks.
     * 
     * Based on mc_shms.f lines ~630-830
     * 
     * Aperture notes from Dave Potterveld (in mc_shms.f):
     * Each aperture is a tilted offset circle of 30 cm radius.
     * Positive tilt angles mean "stub toe before hitting head".
     * 
     * Matrix classes: 23-30 for mid sections 1-7 and exit
     */
    
    // 1. Q3 to D1 transition (mechanical entrance)
    Project(track, drifts_.q3_d1_trans);
    
    // Simple circular check (no tilt, no offset)
    if (!CheckDipoleAperture(track.x, track.y, 0)) {
        stats_.stop_d1_in++;  // Note: using stop_d1_in for this check
        return false;
    }
    
    // 2. D1 Flare (entrance cone)
    Project(track, drifts_.d1_flare);
    
    // Tilted aperture: angle = 9.2°, X offset = -3.5 cm
    double xt = track.x;
    double yt = track.y;
    RotateHAxis(9.200, xt, yt, track.dx, track.dy);
    xt -= 3.5;
    
    if (!CheckDipoleAperture(xt, yt, 1)) {
        stats_.stop_d1_flare++;
        return false;
    }
    
    // 3. D1 Magnetic Entrance (drift, not matrix!)
    Project(track, drifts_.d1_men);
    
    // Tilted aperture: angle = 9.2°, X offset = +2.82 cm
    xt = track.x;
    yt = track.y;
    RotateHAxis(9.200, xt, yt, track.dx, track.dy);
    xt += 2.82;
    
    if (!CheckDipoleAperture(xt, yt, 2)) {
        stats_.stop_d1_men++;
        return false;
    }
    
    // 4-10. D1 Mid sections (7 sections with MATRIX transport)
    // From mc_shms.f: classes 23, 24, 25, 26, 27, 28, 29
    
    struct DipoleSection {
        int matrix_class;
        double angle;
        double offset;
    };
    
    const DipoleSection sections[7] = {
        {23,  6.9,  +8.05},   // Mid 1
        {24,  4.6, +11.75},   // Mid 2
        {25,  2.3, +13.96},   // Mid 3
        {26,  0.0, +14.70},   // Mid 4
        {27, -2.3, +13.96},   // Mid 5
        {28, -4.6, +11.75},   // Mid 6
        {29, -6.9,  +8.05}    // Mid 7
    };
    
    for (int i = 0; i < 7; ++i) {
        TranspMatrix(track, sections[i].matrix_class, drifts_.d1_mid[i]);
        
        xt = track.x;
        yt = track.y;
        RotateHAxis(sections[i].angle, xt, yt, track.dx, track.dy);
        xt += sections[i].offset;
        
        if (!CheckDipoleAperture(xt, yt, 3 + i)) {
            // Increment appropriate counter
            switch(i) {
                case 0: stats_.stop_d1_mid1++; break;
                case 1: stats_.stop_d1_mid2++; break;
                case 2: stats_.stop_d1_mid3++; break;
                case 3: stats_.stop_d1_mid4++; break;
                case 4: stats_.stop_d1_mid5++; break;
                case 5: stats_.stop_d1_mid6++; break;
                case 6: stats_.stop_d1_mid7++; break;
            }
            return false;
        }
    }
    
    // 11. D1 Magnetic Exit (MATRIX transport - class 30)
    TranspMatrix(track, 30, drifts_.d1_mex);
    
    // Tilted aperture: angle = -9.2°, X offset = +2.82 cm
    xt = track.x;
    yt = track.y;
    RotateHAxis(-9.2, xt, yt, track.dx, track.dy);
    xt += 2.82;
    
    if (!CheckDipoleAperture(xt, yt, 10)) {
        stats_.stop_d1_mex++;
        return false;
    }
    
    // 12. D1 Mechanical Exit (drift)
    Project(track, drifts_.d1_out);
    
    // Tilted aperture: angle = -9.2°, X offset = -6.88 cm
    xt = track.x;
    yt = track.y;
    RotateHAxis(-9.20, xt, yt, track.dx, track.dy);
    xt -= 6.88;
    
    if (!CheckDipoleAperture(xt, yt, 11)) {
        stats_.stop_d1_out++;
        return false;
    }
    
    return true;
}

bool SHMS::TransportHut(TrackState& track) {
    /**
     * Transport through detector hut to focal plane.
     * 
     * Based on mc_shms.f lines ~850-870 and mc_shms_hut.f
     * 
     * In Fortran, this drifts to a position before the focal plane
     * (accounting for Cherenkov entrance or vacuum exit), then
     * mc_shms_hut handles the rest including collimator checks.
     * 
     * For now, we simplify by just drifting to focal plane with
     * optional collimator.
     */
    
    if (apertures_.use_collimator) {
        // Drift partway to focal plane
        Project(track, drifts_.fp * 0.5);
        
        // Collimator entrance check
        if (!CheckCollimatorAperture(track.x, track.y, true)) {
            stats_.stop_coll_in++;
            return false;
        }
        
        // Drift through collimator (assume ~20cm length)
        Project(track, 20.0);
        
        // Collimator exit check
        if (!CheckCollimatorAperture(track.x, track.y, false)) {
            stats_.stop_coll_out++;
            return false;
        }
        
        // Drift remaining distance to focal plane
        Project(track, drifts_.fp * 0.5 - 20.0);
    } else {
        // Direct drift to focal plane (no collimator)
        Project(track, drifts_.fp);
    }
    
    return true;
}

// ============================================================================
// Matrix File Parsing
// ============================================================================

bool SHMS::ParseMatrixFile(const std::string& filepath, MatrixElements& matrices) {
    /**
     * Parse COSY matrix file format.
     * 
     * Format:
     * !NAME: <name>
     * !REGION: <region>
     * !OFFSET: <offset>
     * !LENGTH: <length>
     * <5 coeffs> <6 exponents>
     * ...
     * ------------------------------------------------------------------------------
     */
    
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Cannot open matrix file: " << filepath << std::endl;
        return false;
    }
    
    std::string line;
    int current_class = -1;
    double current_length = 0.0;
    
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;
        
        // Check for separator (new transformation class)
        if (line.find("----") != std::string::npos) {
            current_class++;
            if (current_class >= MatrixElements::MAX_CLASSES) {
                std::cerr << "Too many transformation classes in matrix file!" << std::endl;
                return false;
            }
            
            matrices.classes[current_class].length = current_length;
            current_length = 0.0;
            continue;
        }
        
        // Parse length comment
        if (line.find("!LENGTH:") != std::string::npos) {
            size_t pos = line.find("!LENGTH:");
            std::string length_str = line.substr(pos + 8);
            current_length = std::stod(length_str) * 100.0;  // Convert m to cm
            continue;
        }
        
        // Skip other comments
        if (line[0] == '!') {
            continue;
        }
        
        // Parse matrix element line
        if (current_class >= 0) {
            std::istringstream iss(line);
            MatrixTerm term;
            
            // Read 5 coefficients
            for (int i = 0; i < 5; ++i) {
                if (!(iss >> term.coeff[i])) {
                    std::cerr << "Error parsing coefficient " << i << " in line: " << line << std::endl;
                    return false;
                }
            }
            
            // Read 6 exponent digits (as a single integer, then split)
            std::string expon_str;
            if (!(iss >> expon_str)) {
                std::cerr << "Error parsing exponents in line: " << line << std::endl;
                return false;
            }
            
            // Convert 6-digit string to individual exponents
            if (expon_str.length() >= 6) {
                for (int i = 0; i < 5; ++i) {
                    term.expon[i] = expon_str[i] - '0';
                }
                // Note: 6th digit is time-of-flight (ignored)
            }
            
            matrices.classes[current_class].terms.push_back(term);
            matrices.classes[current_class].n_terms++;
        }
    }
    
    file.close();
    
    std::cout << "Loaded matrix file: " << filepath << std::endl;
    std::cout << "  Found " << (current_class + 1) << " transformation classes" << std::endl;
    
    return true;
}

} // namespace simc
