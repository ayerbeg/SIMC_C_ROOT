// src/spectrometers/SOS.cpp
// SOS (Short Orbit Spectrometer) Implementation - Phase 5c.5

#include "simc/spectrometers/SOS.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>

namespace simc {

SOS::SOS() {
    /**
     * Initialize SOS geometry
     * Based on sos/mc_sos.f and sos/apertures_sos.inc
     */
    
    // Quadrupole aperture (from apertures_sos.inc: r2_quad = 163.84 cm^2)
    quad_.radius = 12.8;  // cm (sqrt(163.84))
    
    // BM1 (First Bending Magnet) aperture
    // From apertures_sos.inc - NOTE: Fortran uses rotated coordinates
    // Using entrance values (IN) for conservative estimate
    bm1_.x_min = -8.0;   // cm (w_bm01 = 8.0, horizontal half-gap)
    bm1_.x_max = 8.0;
    bm1_.y_min = -30.52; // cm (b_bm01_in, bottom limit at entrance)
    bm1_.y_max = 41.32;  // cm (t_bm01_in, top limit at entrance)
    
    // BM2 (Second Bending Magnet) aperture
    // Using entrance values (IN) for conservative estimate
    bm2_.x_min = -8.0;   // cm (w_bm02 = 8.0, horizontal half-gap)
    bm2_.x_max = 8.0;
    bm2_.y_min = -66.38; // cm (b_bm02_in, bottom limit at entrance)
    bm2_.y_max = 51.73;  // cm (t_bm02_in, top limit at entrance)
    
    // Drift lengths - will be extracted from matrix files
    drifts_.target_to_quad = 0.0;
    drifts_.quad_entrance = 0.0;
    drifts_.quad_mid = 0.0;
    drifts_.quad_exit = 0.0;
    drifts_.quad_to_bm1 = 0.0;
    drifts_.bm1_length = 0.0;
    drifts_.bm1_to_bm2 = 0.0;
    drifts_.bm2_length = 0.0;
    drifts_.bm2_to_fp = 0.0;
}

bool SOS::Transport(TrackState& track) {
    /**
     * Transport particle through SOS spectrometer
     * Returns true if particle makes it through, false if stopped by aperture
     * 
     * Based on sos/mc_sos.f sequence:
     * 1. Collimator (octagonal slit at z=126.3 cm)
     * 2. Quadrupole (entrance, mid, exit)
     * 3. BM1 (first dipole)
     * 4. BM2 (second dipole)
     * 5. Detector hut
     */
    
    ++stats_.sSTOP_trials;
    
    // Check collimator (octagonal slit) - FROM FORTRAN mc_sos.f
    // Entrance: h_entr = 7.201 cm, v_entr = 4.696 cm
    // Exit: h_exit = 7.567 cm, v_exit = 4.935 cm
    // Position: z_entr = 126.3 cm, thickness = 6.3 cm
    
    // Project to collimator entrance (z = 126.3 cm from target)
    double z_coll_entr = 126.3;  // cm
    double x_coll = track.x + z_coll_entr * track.dx;
    double y_coll = track.y + z_coll_entr * track.dy;
    
    // Check octagonal aperture at entrance
    // Fortran checks: |ys| < h_entr, |xs| < v_entr, and octagonal constraint
    const double h_entr = 7.201;  // cm (horizontal half-width)
    const double v_entr = 4.696;  // cm (vertical half-height)
    
    if (std::abs(y_coll) > h_entr || std::abs(x_coll) > v_entr ||
        std::abs(x_coll) > (-v_entr/h_entr * std::abs(y_coll) + 3.0*v_entr/2.0)) {
        ++stats_.sSTOP_slit;
        return false;
    }
    
    // Check collimator exit (z = 132.6 cm)
    double z_coll_exit = 132.6;  // cm (126.3 + 6.3)
    x_coll = track.x + z_coll_exit * track.dx;
    y_coll = track.y + z_coll_exit * track.dy;
    
    const double h_exit = 7.567;  // cm
    const double v_exit = 4.935;  // cm
    
    if (std::abs(y_coll) > h_exit || std::abs(x_coll) > v_exit ||
        std::abs(x_coll) > (-v_exit/h_exit * std::abs(y_coll) + 3.0*v_exit/2.0)) {
        ++stats_.sSTOP_slit;
        return false;
    }
    
    // Continue with spectrometer elements
    if (!TransportQuadrupole(track)) return false;
    if (!TransportBM1(track)) return false;
    if (!TransportBM2(track)) return false;
    if (!TransportHut(track)) return false;
    
    ++stats_.sSTOP_successes;
    return true;
}

bool SOS::TransportQuadrupole(TrackState& track) {
    /**
     * Transport through quadrupole (3 regions)
     * Based on mc_sos.f
     */
    
    // Transformation class 0: QUAD_ENTRANCE
    // Drift to entrance + check aperture
    TranspMatrix(track, 0, drifts_.quad_entrance);
    
    double r_sq = quad_.radius * quad_.radius;
    if (!CheckQuadAperture(track.x, track.y, r_sq)) {
        ++stats_.sSTOP_quad_in;
        return false;
    }
    
    // Transformation class 1: QUAD_MID
    TranspMatrix(track, 1, drifts_.quad_mid);
    
    if (!CheckQuadAperture(track.x, track.y, r_sq)) {
        ++stats_.sSTOP_quad_mid;
        return false;
    }
    
    // Transformation class 2: QUAD_EXIT
    TranspMatrix(track, 2, drifts_.quad_exit);
    
    if (!CheckQuadAperture(track.x, track.y, r_sq)) {
        ++stats_.sSTOP_quad_out;
        return false;
    }
    
    return true;
}

bool SOS::TransportBM1(TrackState& track) {
    /**
     * Transport through first bending magnet (BM01)
     */
    
    // Drift to BM1 entrance
    Project(track, drifts_.quad_to_bm1);
    
    // Check entrance aperture
    if (!CheckDipoleAperture(track.x, track.y, bm1_)) {
        ++stats_.sSTOP_bm1_in;
        return false;
    }
    
    // Transformation class 3: BM01_ENTRANCE or similar
    // (class index depends on matrix file structure)
    TranspMatrix(track, 3, drifts_.bm1_length);
    
    // Check exit aperture
    if (!CheckDipoleAperture(track.x, track.y, bm1_)) {
        ++stats_.sSTOP_bm1_out;
        return false;
    }
    
    return true;
}

bool SOS::TransportBM2(TrackState& track) {
    /**
     * Transport through second bending magnet (BM02)
     */
    
    // Drift between dipoles
    Project(track, drifts_.bm1_to_bm2);
    
    // Check entrance aperture
    if (!CheckDipoleAperture(track.x, track.y, bm2_)) {
        ++stats_.sSTOP_bm2_in;
        return false;
    }
    
    // Transport through BM2 (class index depends on matrix structure)
    TranspMatrix(track, 4, drifts_.bm2_length);
    
    // Check exit aperture
    if (!CheckDipoleAperture(track.x, track.y, bm2_)) {
        ++stats_.sSTOP_bm2_out;
        return false;
    }
    
    return true;
}

bool SOS::TransportHut(TrackState& track) {
    /**
     * Transport through detector hut to focal plane
     */
    
    // Drift to focal plane
    Project(track, drifts_.bm2_to_fp);
    
    // SOS detector hut is fairly open - most particles make it through
    // Could add detector acceptance cuts here if needed
    
    ++stats_.sSTOP_hut;
    return true;
}

bool SOS::CheckQuadAperture(double x, double y, double r_sq) {
    /**
     * Check if particle hits circular quadrupole aperture
     * Returns false if particle hits (outside aperture)
     */
    return (x*x + y*y) <= r_sq;
}

bool SOS::CheckDipoleAperture(double x, double y, const DipoleAperture& aperture) {
    /**
     * Check rectangular dipole aperture
     * Returns false if particle hits
     */
    if (x < aperture.x_min || x > aperture.x_max) return false;
    if (y < aperture.y_min || y > aperture.y_max) return false;
    return true;
}

void SOS::Project(TrackState& track, double drift_length) {
    /**
     * Drift particle through field-free region
     * Simple linear projection
     */
    track.x += track.dx * drift_length;
    track.y += track.dy * drift_length;
    track.z += drift_length;
    track.pathlen += drift_length * std::sqrt(1.0 + track.dx*track.dx + track.dy*track.dy);
}

void SOS::TranspMatrix(TrackState& track, int class_index, double drift_length) {
    /**
     * Apply matrix transformation for a given class
     * 5th/6th order COSY polynomial transport
     */
    
    if (class_index >= forward_matrix_.n_classes) {
        std::cerr << "Warning: Invalid transformation class " << class_index << std::endl;
        return;
    }
    
    const auto& trans = forward_matrix_.classes[class_index];
    
    // If it's a drift, just project
    if (trans.is_drift) {
        Project(track, drift_length);
        return;
    }
    
    // Save initial state
    double x0 = track.x;
    double dx0 = track.dx;
    double y0 = track.y;
    double dy0 = track.dy;
    double delta0 = track.delta / 100.0;  // Convert % to fraction
    
    // Initialize new state
    double x_new = 0.0;
    double dx_new = 0.0;
    double y_new = 0.0;
    double dy_new = 0.0;
    
    // Evaluate polynomial
    for (const auto& term : trans.terms) {
        // Calculate monomial: x^i * dx^j * y^k * dy^l * delta^m
        double monomial = 1.0;
        
        if (term.expon[0] > 0) monomial *= std::pow(x0, term.expon[0]);
        if (term.expon[1] > 0) monomial *= std::pow(dx0, term.expon[1]);
        if (term.expon[2] > 0) monomial *= std::pow(y0, term.expon[2]);
        if (term.expon[3] > 0) monomial *= std::pow(dy0, term.expon[3]);
        if (term.expon[4] > 0) monomial *= std::pow(delta0, term.expon[4]);
        
        // Add to appropriate output variable
        x_new += term.coeff[0] * monomial;
        dx_new += term.coeff[1] * monomial;
        y_new += term.coeff[2] * monomial;
        dy_new += term.coeff[3] * monomial;
        // term.coeff[4] would be path length change (often zero)
    }
    
    // Update track state (convert back to cm)
    track.x = x_new * 100.0;   // COSY uses meters
    track.dx = dx_new;
    track.y = y_new * 100.0;
    track.dy = dy_new;
    
    // Update z and path length
    track.z += drift_length;
    track.pathlen += drift_length;
}

void SOS::Reconstruct(const FocalPlaneState& fp, TargetState& target) {
    /**
     * Reconstruct target quantities from focal plane
     * Uses reconstruction matrix
     */
    
    // Convert to meters for matrix
    double x = fp.x / 100.0;
    double dx = fp.dx;
    double y = fp.y / 100.0;
    double dy = fp.dy;
    
    // Initialize target quantities
    target.xptar = 0.0;
    target.yptar = 0.0;
    target.delta = 0.0;
    target.ytar = 0.0;
    
    // Apply reconstruction matrix (similar to forward but different coefficients)
    for (int iclass = 0; iclass < recon_matrix_.n_classes; ++iclass) {
        const auto& trans = recon_matrix_.classes[iclass];
        
        for (const auto& term : trans.terms) {
            double monomial = 1.0;
            
            if (term.expon[0] > 0) monomial *= std::pow(x, term.expon[0]);
            if (term.expon[1] > 0) monomial *= std::pow(dx, term.expon[1]);
            if (term.expon[2] > 0) monomial *= std::pow(y, term.expon[2]);
            if (term.expon[3] > 0) monomial *= std::pow(dy, term.expon[3]);
            
            target.xptar += term.coeff[0] * monomial;
            target.yptar += term.coeff[1] * monomial;
            target.delta += term.coeff[2] * monomial;
            target.ytar += term.coeff[3] * monomial;
        }
    }
    
    // Convert delta from fraction to %
    target.delta *= 100.0;
    // Convert ytar back to cm
    target.ytar *= 100.0;
}


bool SOS::GetFocalPlane(const TrackState& track, FocalPlaneState& fp) const {
    // Extract focal plane coordinates
    // NOTE: SOS uses dx/dy (slopes), NOT xp/yp (mrad)
    // NO unit conversion needed
    
    fp.x = track.x;
    fp.y = track.y;
    fp.dx = track.dx;  // NO conversion - stays as slope
    fp.dy = track.dy;  // NO conversion - stays as slope
    
    return true;
}
  
bool SOS::LoadMatrices(const std::string& forward_file, 
                       const std::string& recon_file) {
    /**
     * Load forward and reconstruction matrix files
     */
    
    if (!ParseMatrixFile(forward_file, forward_matrix_)) {
        std::cerr << "Failed to load forward matrix: " << forward_file << std::endl;
        return false;
    }
    
    if (!ParseMatrixFile(recon_file, recon_matrix_)) {
        std::cerr << "Failed to load reconstruction matrix: " << recon_file << std::endl;
        return false;
    }
    
    // Extract drift lengths from transformation classes
    if (forward_matrix_.n_classes >= 3) {
        drifts_.quad_entrance = forward_matrix_.classes[0].length;
        drifts_.quad_mid = forward_matrix_.classes[1].length;
        drifts_.quad_exit = forward_matrix_.classes[2].length;
    }
    
    return true;
}

bool SOS::ParseMatrixFile(const std::string& filepath, MatrixElements& matrices) {
    /**
     * Parse COSY matrix file format
     * Same format as HMS/SHMS
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
        if (line.empty()) continue;
        
        // Check for separator (new transformation class)
        if (line.find("----") != std::string::npos) {
            if (current_class >= 0) {
                matrices.classes[current_class].length = current_length;
                matrices.classes[current_class].is_drift = is_drift;
            }
            
            current_class++;
            if (current_class >= MatrixElements::MAX_CLASSES) {
                std::cerr << "Too many transformation classes!" << std::endl;
                return false;
            }
            
            current_length = 0.0;
            is_drift = false;
            continue;
        }
        
        // Parse length
        if (line.find("!LENGTH:") != std::string::npos) {
            size_t pos = line.find("!LENGTH:");
            std::string length_str = line.substr(pos + 8);
            current_length = std::stod(length_str) * 100.0;  // m to cm
            continue;
        }
        
        // Check for drift
        if (line.find("!DRIFT") != std::string::npos || 
            line.find("drift") != std::string::npos) {
            is_drift = true;
            continue;
        }
        
        // Skip other comments
        if (line[0] == '!') continue;
        
        // Parse matrix element line
        if (current_class >= 0) {
            std::istringstream iss(line);
            MatrixTerm term;
            
            // Read all values flexibly
            std::vector<double> coeffs;
            double val;
            while (iss >> val) {
                coeffs.push_back(val);
            }
            
            if (coeffs.size() < 5) continue;  // Need at least 5 values
            
            // Last value is exponents (integer)
            int expon_int = static_cast<int>(coeffs.back());
            coeffs.pop_back();
            
            // Fill coefficients
            for (size_t i = 0; i < coeffs.size() && i < 5; ++i) {
                term.coeff[i] = coeffs[i];
            }
            for (size_t i = coeffs.size(); i < 5; ++i) {
                term.coeff[i] = 0.0;
            }
            
            // Parse exponents
            std::string expon_str = std::to_string(expon_int);
            while (expon_str.length() < 6) {
                expon_str = "0" + expon_str;
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
        matrices.n_classes = current_class + 1;
    }
    
    file.close();
    
    std::cout << "Loaded SOS matrix file: " << filepath << std::endl;
    std::cout << "  Found " << matrices.n_classes << " transformation classes" << std::endl;
    
    return true;
}

} // namespace simc
