#include "simc/spectrometers/HMS.h"
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

HMS::HMS() {
    // Initialize with default values (already done via member initializers)
}

bool HMS::Transport(TrackState& track) {
    stats_.hSTOP_trials++;
    
    // Transport through each element in sequence
    if (!TransportCollimator(track)) return false;
    if (!TransportQ1(track)) return false;
    if (!TransportQ2(track)) return false;
    if (!TransportQ3(track)) return false;
    if (!TransportDipole(track)) return false;
    if (!TransportPipes(track)) return false;
    if (!TransportHut(track)) return false;
    
    stats_.hSTOP_successes++;
    return true;
}

bool HMS::Reconstruct(const FocalPlaneState& fp, TargetState& target) {
    // TODO: Implement reconstruction
    // For now, just copy values
    target.x = fp.x;
    target.xp = fp.xp;
    target.y = fp.y;
    target.yp = fp.yp;
    target.delta = fp.delta;
    return true;
}

bool HMS::LoadMatrices(const std::string& forward_file, const std::string& recon_file) {
    if (!ParseMatrixFile(forward_file, forward_matrix_)) {
        std::cerr << "Failed to load forward matrix: " << forward_file << std::endl;
        return false;
    }
    
    if (!ParseMatrixFile(recon_file, recon_matrix_)) {
        std::cerr << "Failed to load reconstruction matrix: " << recon_file << std::endl;
        return false;
    }
    
    // Extract drift lengths from transformation classes
    // HMS has 12 classes: 1,4,7,10,12 are drifts; 2-3,5-6,8-9,11 are matrix
    for (int i = 0; i < MatrixElements::MAX_CLASSES; ++i) {
        if (forward_matrix_.classes[i].is_drift || 
            forward_matrix_.classes[i].length > 0.0) {
            
            switch(i) {
                case 1:
                    drifts_.target_to_q1 = forward_matrix_.classes[i].length;
                    break;
                case 4:
                    drifts_.q1_to_q2 = forward_matrix_.classes[i].length;
                    break;
                case 7:
                    drifts_.q2_to_q3 = forward_matrix_.classes[i].length;
                    break;
                case 10:
                    drifts_.q3_to_d1 = forward_matrix_.classes[i].length;
                    break;
                case 12:
                    drifts_.d1_to_fp = forward_matrix_.classes[i].length;
                    break;
            }
        }
    }
    
    return true;
}

// ============================================================================
// Helper Functions
// ============================================================================

void HMS::Project(TrackState& track, double zdrift) {
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

void HMS::TranspMatrix(TrackState& track, int class_id, double zdrift) {
    /**
     * Transport particle through spectrometer element using matrix elements.
     * 
     * Based on shared/transp.f
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

void HMS::RotateHAxis(double angle_deg, double& x, double& y) {
    /**
     * Calculate new trajectory coordinates in reference frame rotated about
     * horizontal axis by angle_deg relative to central ray.
     * 
     * Based on shared/rotate_haxis.f
     * 
     * For HMS dipole aperture checking, we use simplified rotation
     * (without dx, dy arguments) since we're just rotating position.
     */
    
    double angle_rad = angle_deg * DEG_TO_RAD;
    double cos_a = std::cos(angle_rad);
    double sin_a = std::sin(angle_rad);
    
    double x_new = x * cos_a + y * sin_a;
    double y_new = -x * sin_a + y * cos_a;
    
    x = x_new;
    y = y_new;
}

// ============================================================================
// Aperture Checking
// ============================================================================

bool HMS::CheckCollimatorEntrance(double x, double y) {
    /**
     * Check collimator entrance (octagonal).
     * 
     * Based on mc_hms.f lines ~125-145
     */
    
    // Apply offsets
    double x_check = x - collimator_.x_off;
    double y_check = y - collimator_.y_off;
    
    // Horizontal slit check
    if (std::abs(y_check) > collimator_.h_entr) {
        return false;
    }
    
    // Vertical slit check
    if (std::abs(x_check) > collimator_.v_entr) {
        return false;
    }
    
    // Octagonal corner cuts
    // if (abs(x - x_off) > (-v_entr/h_entr * abs(y - y_off) + 3*v_entr/2))
    double threshold = -collimator_.v_entr / collimator_.h_entr * 
                       std::abs(y_check) + 3.0 * collimator_.v_entr / 2.0;
    
    if (std::abs(x_check) > threshold) {
        return false;
    }
    
    return true;
}

bool HMS::CheckCollimatorExit(double x, double y) {
    /**
     * Check collimator exit (octagonal, slightly larger).
     */
    
    // Apply offsets
    double x_check = x - collimator_.x_off;
    double y_check = y - collimator_.y_off;
    
    // Horizontal slit check
    if (std::abs(y_check) > collimator_.h_exit) {
        return false;
    }
    
    // Vertical slit check
    if (std::abs(x_check) > collimator_.v_exit) {
        return false;
    }
    
    // Octagonal corner cuts
    double threshold = -collimator_.v_exit / collimator_.h_exit * 
                       std::abs(y_check) + 3.0 * collimator_.v_exit / 2.0;
    
    if (std::abs(x_check) > threshold) {
        return false;
    }
    
    return true;
}

bool HMS::CheckQuadAperture(double x, double y, double radius) {
    /**
     * Check if particle passes through circular quadrupole aperture.
     */
    return (x*x + y*y) <= (radius * radius);
}

bool HMS::HitDipole(double x, double y) {
    /**
     * Check if particle hits dipole aperture (complex 6-region shape).
     * 
     * Based on hit_dipole() function in mc_hms.f lines ~506-550
     * 
     * Returns true if particle is OUTSIDE aperture (hit = bad)
     * Returns false if particle is inside aperture (miss = good)
     * 
     * Uses symmetry to check only first quadrant.
     */
    
    // Use absolute values (exploit symmetry)
    double x_local = std::abs(x);
    double y_local = std::abs(y);
    
    // Check 6 regions defining the aperture
    bool check1 = (x_local <= apertures_.x_d1 && y_local <= apertures_.y_d1);
    bool check2 = (x_local <= apertures_.x_d2 && y_local <= apertures_.y_d2);
    bool check3 = (x_local <= apertures_.x_d3 && y_local <= apertures_.y_d3);
    bool check4 = (x_local <= apertures_.x_d4 && y_local <= apertures_.y_d4);
    
    // Rounded corner (circle centered at x_d5, y_d5)
    double dx = x_local - apertures_.x_d5;
    double dy = y_local - apertures_.y_d5;
    bool check5 = (dx*dx + dy*dy) <= (apertures_.r_d5 * apertures_.r_d5);
    
    // Slanted piece (line: y = a*x + b)
    bool check6 = (x_local >= apertures_.x_d4 && 
                   x_local <= apertures_.x_d3 &&
                   (y_local - apertures_.a_d6 * x_local - apertures_.b_d6) <= 0.0);
    
    // Inside aperture if in ANY region
    bool inside_aperture = check1 || check2 || check3 || check4 || check5 || check6;
    
    // Fortran returns "hit_dipole" (outside), we return same convention
    return !inside_aperture;
}

bool HMS::CheckPipeAperture(double x, double y, double x_offset, double y_offset, double r_sq) {
    /**
     * Check if particle passes through circular pipe.
     * Uses precomputed radius squared for efficiency.
     */
    
    double dx = x - x_offset;
    double dy = y - y_offset;
    
    return (dx*dx + dy*dy) <= r_sq;
}

// ============================================================================
// Element Transport Functions
// ============================================================================

bool HMS::TransportCollimator(TrackState& track) {
    /**
     * Transport through HMS-100 collimator (octagonal).
     * 
     * Based on mc_hms.f lines ~115-180
     * 
     * Collimator can be disabled for testing (wide-open mode).
     * Detailed pion absorption model (mc_hms_coll) not implemented yet.
     */
    
    if (!collimator_.use_collimator) {
        // Wide open - skip collimator entirely
        return true;
    }
    
    // If using HMS collimator simulation (for pions/muons)
    if (use_hms_coll_sim_ && track.m2 > 100.0*100.0 && track.m2 < 200.0*200.0) {
        // TODO: Implement detailed pion absorption (mc_hms_coll.f)
        // For now, just use aperture checks
        std::cerr << "HMS collimator simulation not implemented yet!" << std::endl;
        stats_.hSTOP_coll++;
        return false;
    }
    
    // Standard aperture checks (entrance and exit)
    
    // Drift to collimator entrance
    Project(track, collimator_.z_entr - track.z);
    
    // Check entrance aperture
    if (!CheckCollimatorEntrance(track.x, track.y)) {
        // Determine which part of aperture was hit
        double x_check = track.x - collimator_.x_off;
        double y_check = track.y - collimator_.y_off;
        
        if (std::abs(y_check) > collimator_.h_entr) {
            stats_.hSTOP_slit_hor++;
        } else if (std::abs(x_check) > collimator_.v_entr) {
            stats_.hSTOP_slit_vert++;
        } else {
            stats_.hSTOP_slit_oct++;
        }
        return false;
    }
    
    // Drift to collimator exit
    Project(track, collimator_.z_exit - collimator_.z_entr);
    
    // Check exit aperture
    if (!CheckCollimatorExit(track.x, track.y)) {
        double x_check = track.x - collimator_.x_off;
        double y_check = track.y - collimator_.y_off;
        
        if (std::abs(y_check) > collimator_.h_exit) {
            stats_.hSTOP_slit_hor++;
        } else if (std::abs(x_check) > collimator_.v_exit) {
            stats_.hSTOP_slit_vert++;
        } else {
            stats_.hSTOP_slit_oct++;
        }
        return false;
    }
    
    return true;
}

bool HMS::TransportQ1(TrackState& track) {
    /**
     * Transport through Q1 (Quadrupole 1) with 3 aperture checks.
     * 
     * Based on mc_hms.f lines ~182-202
     * 
     * Transformations:
     * 1. Drift to Q1 entrance (class 1)
     * 2. Matrix through Q1 to 2/3 point (class 2)
     * 3. Matrix from 2/3 to Q1 exit (class 3)
     */
    
    // Calculate drift from current position to Q1 entrance
    // Current z position depends on whether collimator was used
    double z_start = collimator_.use_collimator ? collimator_.z_exit : track.z;
    double zdrift = drifts_.target_to_q1 - z_start;
    
    // NOTE: Fortran checks if transformation 1 is a drift
    // We assume it is and use Project()
    Project(track, zdrift);
    
    // Check Q1 entrance
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q1)) {
        stats_.hSTOP_Q1_in++;
        return false;
    }
    
    // Transport to Q1 2/3 point (matrix transformation class 2)
    TranspMatrix(track, 2, drifts_.q1_mid);
    
    // Check Q1 mid aperture
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q1)) {
        stats_.hSTOP_Q1_mid++;
        return false;
    }
    
    // Transport to Q1 exit (matrix transformation class 3)
    TranspMatrix(track, 3, drifts_.q1_exit);
    
    // Check Q1 exit aperture
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q1)) {
        stats_.hSTOP_Q1_out++;
        return false;
    }
    
    return true;
}

bool HMS::TransportQ2(TrackState& track) {
    /**
     * Transport through Q2 (Quadrupole 2) with 3 aperture checks.
     * 
     * Based on mc_hms.f lines ~204-224
     */
    
    // Drift to Q2 entrance (transformation 4)
    Project(track, drifts_.q1_to_q2);
    
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q2)) {
        stats_.hSTOP_Q2_in++;
        return false;
    }
    
    // Transport to Q2 2/3 point (matrix transformation class 5)
    TranspMatrix(track, 5, drifts_.q2_mid);
    
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q2)) {
        stats_.hSTOP_Q2_mid++;
        return false;
    }
    
    // Transport to Q2 exit (matrix transformation class 6)
    TranspMatrix(track, 6, drifts_.q2_exit);
    
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q2)) {
        stats_.hSTOP_Q2_out++;
        return false;
    }
    
    return true;
}

bool HMS::TransportQ3(TrackState& track) {
    /**
     * Transport through Q3 (Quadrupole 3) with 3 aperture checks.
     * 
     * Based on mc_hms.f lines ~226-246
     */
    
    // Drift to Q3 entrance (transformation 7)
    Project(track, drifts_.q2_to_q3);
    
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q3)) {
        stats_.hSTOP_Q3_in++;
        return false;
    }
    
    // Transport to Q3 2/3 point (matrix transformation class 8)
    TranspMatrix(track, 8, drifts_.q3_mid);
    
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q3)) {
        stats_.hSTOP_Q3_mid++;
        return false;
    }
    
    // Transport to Q3 exit (matrix transformation class 9)
    TranspMatrix(track, 9, drifts_.q3_exit);
    
    if (!CheckQuadAperture(track.x, track.y, apertures_.r_Q3)) {
        stats_.hSTOP_Q3_out++;
        return false;
    }
    
    return true;
}

bool HMS::TransportDipole(TrackState& track) {
    /**
     * Transport through Dipole (D1) with complex aperture shape.
     * 
     * Based on mc_hms.f lines ~248-290
     * 
     * Dipole aperture has complex 6-region shape (Niculescu model).
     * Entrance and exit are tilted by ±6 degrees.
     */
    
    // Drift to D1 entrance (transformation 10)
    Project(track, drifts_.q3_to_d1);
    
    // Rotate coordinates for tilted entrance aperture (-6 degrees)
    double xt = track.x;
    double yt = track.y;
    RotateHAxis(-6.0, xt, yt);
    
    if (HitDipole(xt, yt)) {
        stats_.hSTOP_D1_in++;
        return false;
    }
    
    // Transport through dipole (matrix transformation class 11)
    TranspMatrix(track, 11, drifts_.d1_length);
    
    // Rotate coordinates for tilted exit aperture (+6 degrees)
    xt = track.x;
    yt = track.y;
    RotateHAxis(+6.0, xt, yt);
    
    if (HitDipole(xt, yt)) {
        stats_.hSTOP_D1_out++;
        return false;
    }
    
    // Check odd-shaped interface piece (30.48 cm radius pipe with rectangular cut)
    // From mc_hms.f lines ~285-290
    if ((std::pow(xt - apertures_.pipe_x_offset, 2) + 
         std::pow(yt - apertures_.pipe_y_offset, 2) > 30.48*30.48) ||
        (std::abs(yt - apertures_.pipe_y_offset) > 20.5232)) {
        stats_.hSTOP_D1_out++;
        return false;
    }
    
    return true;
}

bool HMS::TransportPipes(TrackState& track) {
    /**
     * Transport through post-dipole vacuum pipes (3 pipes).
     * 
     * Based on mc_hms.f lines ~292-320
     * 
     * Three pipes:
     * 1. 26.65" pipe (exit at z_pipe1)
     * 2. 117" pipe (exit at z_pipe2)
     * 3. 45.5" pipe (exit at z_pipe3)
     */
    
    // Pipe 1: Drift and check exit of 26.65" pipe
    Project(track, apertures_.z_pipe1);
    
    if (!CheckPipeAperture(track.x, track.y, 
                          apertures_.pipe_x_offset, 
                          apertures_.pipe_y_offset,
                          apertures_.r_pipe1_sq)) {
        stats_.hSTOP_D1_out++;  // Note: Fortran uses D1_out counter for pipes
        return false;
    }
    
    // Pipe 2: Drift and check exit of 117" pipe
    Project(track, apertures_.z_pipe2 - apertures_.z_pipe1);
    
    if (!CheckPipeAperture(track.x, track.y,
                          apertures_.pipe_x_offset,
                          apertures_.pipe_y_offset,
                          apertures_.r_pipe2_sq)) {
        stats_.hSTOP_D1_out++;
        return false;
    }
    
    // Pipe 3: Drift and check exit of 45.5" pipe
    Project(track, apertures_.z_pipe3 - apertures_.z_pipe2);
    
    if (!CheckPipeAperture(track.x, track.y,
                          apertures_.pipe_x_offset,
                          apertures_.pipe_y_offset,
                          apertures_.r_pipe3_sq)) {
        stats_.hSTOP_D1_out++;
        return false;
    }
    
    return true;
}

bool HMS::TransportHut(TrackState& track) {
    /**
     * Transport through detector hut to focal plane.
     * 
     * Based on mc_hms.f lines ~322-340 and mc_hms_hut.f
     * 
     * From mc_hms.f comments:
     * "Note that we do NOT transport (project) to focal plane. We will do this
     *  in mc_hms_hut.f so that it can take care of all of the decay, mult. scatt,
     *  and apertures."
     * 
     * Initial z position for mc_hms_hut is -147.48 cm so that the sum of four
     * drift lengths between pipe and focal plane is 625.0 cm:
     * (64.77 + 297.18 + 115.57 + 147.48 = 625)
     * 
     * For simplified implementation, we just drift to focal plane.
     */
    
    // Calculate remaining drift to focal plane
    // From mc_hms.f: zdrift = driftdist(spectr,12) - z_dip3
    // where driftdist(spectr,12) is the drift from D1 exit to focal plane
    double zdrift = drifts_.d1_to_fp - apertures_.z_pipe3;
    
    // Drift to focal plane
    Project(track, zdrift);
    
    // TODO: Implement detector hut checks (DC1, DC2, scintillators, calorimeter)
    // For now, count as success if we got this far
    stats_.hSTOP_hut++;
    
    return true;
}

// ============================================================================
// Matrix File Parsing
// ============================================================================

bool HMS::ParseMatrixFile(const std::string& filepath, MatrixElements& matrices) {
    /**
     * Parse COSY matrix file format.
     * 
     * Same format as SHMS matrices.
     */
    
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Cannot open matrix file: " << filepath << std::endl;
        return false;
    }
    
    std::string line;
    int current_class = -1;
    double current_length = 0.0;
    bool is_drift = false;
    
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;
        
        // Check for separator (new transformation class)
        if (line.find("----") != std::string::npos) {
            // Finalize previous class
            if (current_class >= 0) {
                matrices.classes[current_class].length = current_length;
                matrices.classes[current_class].is_drift = is_drift;
            }
            
            // Start new class
            current_class++;
            if (current_class >= MatrixElements::MAX_CLASSES) {
                std::cerr << "Too many transformation classes in matrix file!" << std::endl;
                return false;
            }
            
            current_length = 0.0;
            is_drift = false;
            continue;
        }
        
        // Parse length comment
        if (line.find("!LENGTH:") != std::string::npos) {
            size_t pos = line.find("!LENGTH:");
            std::string length_str = line.substr(pos + 8);
            current_length = std::stod(length_str) * 100.0;  // Convert m to cm
            continue;
        }
        
        // Check for drift indicator
        if (line.find("!DRIFT") != std::string::npos || 
            line.find("drift") != std::string::npos) {
            is_drift = true;
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
            
            // Try to read coefficients - be flexible about how many
            std::vector<double> coeffs;
            double val;
            while (iss >> val) {
                coeffs.push_back(val);
            }
            
                      
            // Check if we have at least 6 values (5 coeffs + exponents OR 4 coeffs + exponents)
            if (coeffs.size() < 5) {
                std::cerr << "Error: Expected at least 5 values in line: " << line << std::endl;
                std::cerr << "  Got " << coeffs.size() << " values" << std::endl;
                return false;
            }
            
            // Last value should be exponents (integer)
            int expon_int = static_cast<int>(coeffs.back());
            coeffs.pop_back();  // Remove exponents from coeffs
            
            // Now coeffs contains only the matrix coefficients (should be 4 or 5)
            // Fill in the term structure
            for (size_t i = 0; i < coeffs.size() && i < 5; ++i) {
                term.coeff[i] = coeffs[i];
            }
            // Fill remaining with zeros if needed
            for (size_t i = coeffs.size(); i < 5; ++i) {
                term.coeff[i] = 0.0;
            }
            
            // Parse exponents from integer
            // Format: ABCDEF where each digit is an exponent
            std::string expon_str = std::to_string(expon_int);
            while (expon_str.length() < 6) {
                expon_str = "0" + expon_str;  // Pad with leading zeros
            }
            
            for (int i = 0; i < 5; ++i) {
                term.expon[i] = expon_str[i] - '0';
            }
            
            matrices.classes[current_class].terms.push_back(term);
            matrices.classes[current_class].n_terms++;
        }
    }
    
    // Finalize last class
    if (current_class >= 0) {
        matrices.classes[current_class].length = current_length;
        matrices.classes[current_class].is_drift = is_drift;
    }
    
    file.close();
    
    std::cout << "Loaded HMS matrix file: " << filepath << std::endl;
    std::cout << "  Found " << (current_class + 1) << " transformation classes" << std::endl;
    
    return true;
}

} // namespace simc
