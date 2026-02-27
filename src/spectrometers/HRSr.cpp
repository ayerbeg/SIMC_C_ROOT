/**
 * HRSr.cpp - Complete implementation with all 32 aperture checks
 * Based on mc_hrsr.f from simc_gfortran
 */

#include "simc/spectrometers/HRSr.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

namespace simc {

// Helper function to rotate coordinates (like Fortran rotate_haxis)
static void RotateHAxis(double angle_deg, double& x, double& y) {
    double angle_rad = angle_deg * M_PI / 180.0;
    double cos_a = std::cos(angle_rad);
    double sin_a = std::sin(angle_rad);
    double x_rot = x * cos_a - y * sin_a;
    double y_rot = x * sin_a + y * cos_a;
    x = x_rot;
    y = y_rot;
}

HRSr::HRSr() {
    // Initialize apertures (from apertures_hrsr.inc)
    q1_.radius = 15.0;      // Q1 radius (cm)
    q2_.radius = 30.22;     // Q2 radius (cm)
    q3_.radius = 30.22;     // Q3 radius (cm)
}

bool HRSr::Transport(TrackState& track) {
    /**
     * Transport through HRSr spectrometer
     * Following EXACT sequence from mc_hrsr.f with ALL 32 aperture checks
     * 
     * Key differences from simplified version:
     * - Uses Project() for all field-free drifts (matches Fortran)
     * - Tracks cumulative z position with ztmp (matches Fortran ztmp logic)
     * - Uses driftdist to know where each transp() call starts
     * - Checks apertures at exact z positions from Fortran
     */
    
    ++stats_.rSTOP_trials;
    
    // Collimator dimensions (from mc_hrsr.f lines 53-65)
    const double h_entr = 3.145;   // horizontal half-width at entrance (cm)
    const double v_entr = 6.090;   // vertical half-height at entrance (cm)
    const double h_exit = 3.340;   // horizontal half-width at exit (cm)
    const double v_exit = 6.485;   // vertical half-height at exit (cm)
    const double y_off = 0.0;      // horizontal offset
    const double x_off = 0.0;      // vertical offset
    const double z_entr = 110.0;   // entrance z position (cm)
    const double z_exit = z_entr + 8.0;  // exit (8 cm thick)
    
    double ztmp = 0.0;  // Tracks cumulative distance from reference point
    double zdrift = 0.0;
    double xt, yt;  // Temporary coordinates
    
    // APERTURE 1: Circular aperture before collimator (z=65.686 cm)
    // Fortran lines 151-159
    zdrift = 65.686;
    ztmp = zdrift;
    Project(track, zdrift);
    if (std::sqrt(track.x*track.x + track.y*track.y) > 7.3787) {
        ++stats_.rSTOP_slit;
        return false;
    }
    
    // APERTURE 2: Circular aperture (z=80.436 cm)
    // Fortran lines 161-169
    zdrift = 80.436 - ztmp;
    ztmp = 80.436;
    Project(track, zdrift);
    if (std::sqrt(track.x*track.x + track.y*track.y) > 7.4092) {
        ++stats_.rSTOP_slit;
        return false;
    }
    
    // APERTURE 3: Collimator entrance (z=110 cm)
    // Fortran lines 173-181
    zdrift = z_entr - ztmp;
    Project(track, zdrift);
    if (std::abs(track.y - y_off) > h_entr) {
        ++stats_.rSTOP_slit;
        return false;
    }
    if (std::abs(track.x - x_off) > v_entr) {
        ++stats_.rSTOP_slit;
        return false;
    }
    
    // APERTURE 4: Collimator exit (z=118 cm)
    // Fortran lines 185-197
    zdrift = z_exit - z_entr;
    Project(track, zdrift);
    if (std::abs(track.y - y_off) > h_exit) {
        ++stats_.rSTOP_slit;
        return false;
    }
    if (std::abs(track.x - x_off) > v_exit) {
        ++stats_.rSTOP_slit;
        return false;
    }
    
    // APERTURE 5: Before Q1 (z=135.064 cm)
    // Fortran lines 201-209
    ztmp = 135.064;
    zdrift = ztmp - z_exit;
    Project(track, zdrift);
    if (std::sqrt(track.x*track.x + track.y*track.y) > 12.5222) {
        ++stats_.rSTOP_Q1_in;
        return false;
    }
    
    // Drift to Q1 entrance
    // Fortran lines 213-215: "zdrift = driftdist(spectr,1) - ztmp"
    // driftdist(spectr,1) = 159.03 cm (from hrs_driftlengths.dat)
    zdrift = 159.03 - ztmp;
    Project(track, zdrift);
    
    // APERTURE 6: Q1 entrance
    // Fortran lines 216-222
    static int q1_check_count = 0;
    q1_check_count++;
    if (q1_check_count == 1) {
        std::cout << "DEBUG: First particle at Q1 entrance check: x=" << track.x 
                  << " y=" << track.y << " r=" << std::sqrt(track.x*track.x + track.y*track.y)
                  << " (limit r=" << q1_.radius << ")" << std::endl;
    }
    
    if (track.x*track.x + track.y*track.y > q1_.radius*q1_.radius) {
        ++stats_.rSTOP_Q1_in;
        return false;
    }
    
    // APERTURE 7: Q1 mid (transformation class 2)
    // Fortran lines 226-232
    TranspMatrix(track, 1, 62.75333333);  // Class index 1 = Fortran class 2
    if (track.x*track.x + track.y*track.y > q1_.radius*q1_.radius) {
        ++stats_.rSTOP_Q1_mid;
        return false;
    }
    
    // APERTURE 8: Q1 exit (transformation class 3)
    // Fortran lines 236-242
    TranspMatrix(track, 2, 31.37666667);  // Class index 2 = Fortran class 3
    if (track.x*track.x + track.y*track.y > q1_.radius*q1_.radius) {
        ++stats_.rSTOP_Q1_out;
        return false;
    }
    
    // APERTURE 9-10: Between Q1 and Q2
    // Fortran lines 246-262
    // Q1 exit is at z=253.16 cm
    zdrift = 300.464 - 253.16;
    ztmp = zdrift;  // Distance from Q1 exit
    Project(track, zdrift);
    if (std::sqrt(track.x*track.x + track.y*track.y) > 14.9225) {
        ++stats_.rSTOP_Q1_out;
        return false;
    }
    
    zdrift = 314.464 - 300.464;
    ztmp = ztmp + zdrift;  // Cumulative distance from Q1 exit
    Project(track, zdrift);
    if (std::sqrt(track.x*track.x + track.y*track.y) > 20.9550) {
        ++stats_.rSTOP_Q2_in;
        return false;
    }
    
    // Drift to Q2 entrance
    // Fortran lines 266-268: "zdrift = driftdist(spectr,4) - ztmp"
    // driftdist(spectr,4) = 370.36 - 253.16 = 117.20 cm from Q1 exit
    zdrift = 117.20 - ztmp;
    Project(track, zdrift);
    
    // APERTURE 11: Q2 entrance
    // Fortran lines 269-275
    if (track.x*track.x + track.y*track.y > q2_.radius*q2_.radius) {
        ++stats_.rSTOP_Q2_in;
        return false;
    }
    
    // APERTURE 12: Q2 mid (transformation class 5)
    // Fortran lines 279-285
    TranspMatrix(track, 4, 121.77333333);  // Class index 4 = Fortran class 5
    if (track.x*track.x + track.y*track.y > q2_.radius*q2_.radius) {
        ++stats_.rSTOP_Q2_mid;
        return false;
    }
    
    // APERTURE 13: Q2 exit (transformation class 6)
    // Fortran lines 289-295
    TranspMatrix(track, 5, 60.88666667);  // Class index 5 = Fortran class 6
    if (track.x*track.x + track.y*track.y > q2_.radius*q2_.radius) {
        ++stats_.rSTOP_Q2_out;
        return false;
    }
    
    // APERTURE 14-16: Between Q2 and D1
    // Fortran lines 299-325
    // Q2 exit is at z=553.020 cm
    zdrift = 609.664 - 553.020;
    ztmp = zdrift;  // Distance from Q2 exit
    Project(track, zdrift);
    if (std::sqrt(track.x*track.x + track.y*track.y) > 30.0073) {
        ++stats_.rSTOP_Q2_out;
        return false;
    }
    
    zdrift = 641.800 - 609.664;
    ztmp = ztmp + zdrift;
    Project(track, zdrift);
    if (std::sqrt(track.x*track.x + track.y*track.y) > 30.0073) {
        ++stats_.rSTOP_Q2_out;
        return false;
    }
    
    zdrift = 819.489 - 641.800;
    ztmp = ztmp + zdrift;
    Project(track, zdrift);
    if (std::abs(track.x) > 50.0 || std::abs(track.y) > 15.0) {
        ++stats_.rSTOP_D1_in;
        return false;
    }
    
    // Drift to D1 entrance
    // Fortran lines 329-331: "zdrift = driftdist(spectr,7) - ztmp"
    // driftdist(spectr,7) = 996.1 - 553.02 = 443.08 cm from Q2 exit
    zdrift = 443.08 - ztmp;
    Project(track, zdrift);
    
    // APERTURE 17: D1 entrance (rotated -30°)
    // Fortran lines 332-346
    xt = track.x;
    yt = track.y;
    RotateHAxis(-30.0, xt, yt);
    if (std::abs(xt - 2.500) > 52.5) {  // -50 < x < +55
        ++stats_.rSTOP_D1_in;
        return false;
    }
    if (std::abs(yt) + 0.01861 * xt > 12.5) {  // tan(1.066°) ≈ 0.01861
        ++stats_.rSTOP_D1_in;
        return false;
    }
    
    // APERTURE 18: D1 exit (transformation class 8, rotated +30°)
    // Fortran lines 350-364
    TranspMatrix(track, 7, 659.73445725);  // Class index 7 = Fortran class 8
    xt = track.x;
    yt = track.y;
    RotateHAxis(30.0, xt, yt);
    if (std::abs(xt - 2.500) > 52.5) {
        ++stats_.rSTOP_D1_out;
        return false;
    }
    if (std::abs(yt) + 0.01861 * xt > 12.5) {
        ++stats_.rSTOP_D1_out;
        return false;
    }
    
    // APERTURE 19-21: Between D1 and Q3
    // Fortran lines 369-393
    // D1 exit is at z=1655.83446 cm
    zdrift = 1745.33546 - 1655.83446;
    ztmp = zdrift;  // Distance from D1 exit
    Project(track, zdrift);
    if (std::sqrt(track.x*track.x + track.y*track.y) > 30.3276) {
        ++stats_.rSTOP_D1_out;
        return false;
    }
    if (std::abs(track.x) > 50.0 || std::abs(track.y) > 15.0) {
        ++stats_.rSTOP_D1_out;
        return false;
    }
    
    zdrift = 1759.00946 - 1745.33546;
    ztmp = ztmp + zdrift;
    Project(track, zdrift);
    if (std::sqrt(track.x*track.x + track.y*track.y) > 30.3276) {
        ++stats_.rSTOP_Q3_in;
        return false;
    }
    
    // Drift to Q3 entrance
    // Fortran lines 397-399: "zdrift = driftdist(spectr,9) - ztmp"
    // driftdist(spectr,9) = 1815.08446 - 1655.83446 = 159.25 cm from D1 exit
    zdrift = 159.25 - ztmp;
    Project(track, zdrift);
    
    // APERTURE 22: Q3 entrance
    // Fortran lines 400-406
    if (track.x*track.x + track.y*track.y > q3_.radius*q3_.radius) {
        ++stats_.rSTOP_Q3_in;
        return false;
    }
    
    // APERTURE 23: Q3 mid (transformation class 10)
    // Fortran lines 410-416
    TranspMatrix(track, 9, 121.7866667);  // Class index 9 = Fortran class 10
    if (track.x*track.x + track.y*track.y > q3_.radius*q3_.radius) {
        ++stats_.rSTOP_Q3_mid;
        return false;
    }
    
    // APERTURE 24: Q3 exit (transformation class 11)
    // Fortran lines 420-426
    TranspMatrix(track, 10, 60.89333333);  // Class index 10 = Fortran class 11
    if (track.x*track.x + track.y*track.y > q3_.radius*q3_.radius) {
        ++stats_.rSTOP_Q3_out;
        return false;
    }
    
    // APERTURE 25-26: After Q3
    // Fortran lines 430-449
    // Q3 exit is at z=1909.214457 cm (1997.76446 in Fortran, but that's cumulative)
    // Actually: Q3 OUT transformation ends at 1909.214457 based on driftlengths.dat
    zdrift = 2080.38746 - 1909.214457;
    ztmp = zdrift;  // Distance from Q3 exit
    Project(track, zdrift);
    if (std::abs(track.x) > 35.56 || std::abs(track.y) > 17.145) {
        ++stats_.rSTOP_Q3_out;
        return false;
    }
    
    zdrift = 2327.47246 - 2080.38746;  // Vacuum window
    ztmp = ztmp + zdrift;
    Project(track, zdrift);
    if (std::abs(track.x) > 99.76635 || std::abs(track.y) > 17.145) {
        ++stats_.rSTOP_Q3_out;
        return false;
    }
    
    // Final drift to focal plane
    // Fortran lines 455-456: "zdrift = driftdist(spectr,12) - ztmp"
    // driftdist(spectr,12) = 345.23 cm from Q3 exit (total = 2254.444457)
    zdrift = 345.23 - ztmp;
    Project(track, zdrift);
    
    // Made it to the hut!
    ++stats_.rSTOP_hut;
    ++stats_.rSTOP_successes;
    return true;
}

bool HRSr::Reconstruct(const FocalPlaneState& fp, TargetState& target) {
    /**
     * Reconstruct target quantities from focal plane
     * Uses reconstruction matrix
     */
    
    // COSY 7 matrices use: cm for positions, mrad for angles
    double x = fp.x;                  // Keep in cm
    double dx = fp.dx * 1000.0;       // slope → mrad
    double y = fp.y;                  // Keep in cm
    double dy = fp.dy * 1000.0;       // slope → mrad
    
    // Initialize target quantities
    target.xptar = 0.0;
    target.yptar = 0.0;
    target.delta = 0.0;
    target.ytar = 0.0;
    
    // Apply reconstruction matrix
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
    
    return true;
}

bool HRSr::GetFocalPlane(const TrackState& track, FocalPlaneState& fp) const {
    // Extract focal plane coordinates
    // NOTE: HRSr uses dx/dy (slopes), NOT xp/yp (mrad)
    
    fp.x = track.x;
    fp.y = track.y;
    fp.dx = track.dx;  // NO conversion
    fp.dy = track.dy;  // NO conversion
    
    return true;
}


  
bool HRSr::LoadMatrices(const std::string& forward_file, const std::string& recon_file) {
    if (!ParseMatrixFile(forward_file, forward_matrix_)) {
        return false;
    }
    if (!ParseMatrixFile(recon_file, recon_matrix_)) {
        return false;
    }
    return true;
}

bool HRSr::ParseMatrixFile(const std::string& filepath, MatrixElements& matrices) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open matrix file: " << filepath << std::endl;
        return false;
    }
    
    matrices.n_classes = 0;
    int current_class = -1;
    double current_length = 0.0;
    bool is_drift = false;
    
    std::string line;
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty() || line.find_first_not_of(" \t\n\r") == std::string::npos) {
            continue;
        }
        
        // Check for new transformation class (starts with !NAME:)
        if (line.find("!NAME:") != std::string::npos) {
            // Finalize previous class
            if (current_class >= 0) {
                matrices.classes[current_class].length = current_length;
                matrices.classes[current_class].is_drift = is_drift;
            }
            
            // Start new class
            ++current_class;
            if (current_class >= MatrixElements::MAX_CLASSES) {
                std::cerr << "Too many transformation classes" << std::endl;
                return false;
            }
            matrices.classes[current_class].n_terms = 0;
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
        if (line[0] == '!') {
            continue;
        }
        
        // Skip separator lines
        if (line.find("---") != std::string::npos) {
            continue;
        }
        
        // Parse matrix term
        if (current_class >= 0) {
            std::istringstream iss(line);
            MatrixTerm term;
            
            // Read all values flexibly (handles missing spaces in Fortran files)
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
            
            // Parse exponents (pad to 6 digits if needed)
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
    
    std::cout << "Loaded HRSr matrix file: " << filepath << std::endl;
    std::cout << "  Found " << matrices.n_classes << " transformation classes" << std::endl;
    
    return true;
}

void HRSr::TranspMatrix(TrackState& track, int class_num, double pathlen) {
    /**
     * Apply COSY transformation matrix to track state
     * Similar to Fortran transp() subroutine
     * 
     * CRITICAL: COSY 7 uses cm for positions, mrad for angles
     */
    
    if (class_num >= forward_matrix_.n_classes) {
        std::cerr << "ERROR: Invalid transformation class " << class_num << std::endl;
        return;
    }
    
    const auto& trans = forward_matrix_.classes[class_num];
    
    // Use drift length from matrix if pathlen is 0
    double drift_length = (pathlen > 0.0) ? pathlen : trans.length;
    
    // COSY 7 matrices use: cm for positions, mrad for angles
    // Do NOT convert cm to meters!
    double x0 = track.x;              // cm (not meters!)
    double dx0 = track.dx * 1000.0;   // slope → mrad
    double y0 = track.y;              // cm (not meters!)
    double dy0 = track.dy * 1000.0;   // slope → mrad
    double delta0 = track.delta / 100.0;  // Convert % to fraction
    
    // Initialize new values
    double x_new = 0.0, dx_new = 0.0, y_new = 0.0, dy_new = 0.0;
    
    // Apply matrix transformation
    for (const auto& term : trans.terms) {
        double monomial = 1.0;
        
        // Compute monomial: x^i * dx^j * y^k * dy^l * delta^m
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
    }
    
    // Update track state (already in cm, just convert angles from mrad)
    track.x = x_new;                  // Already in cm
    track.dx = dx_new / 1000.0;       // mrad → slope
    track.y = y_new;                  // Already in cm
    track.dy = dy_new / 1000.0;       // mrad → slope
    track.z += drift_length;
    track.pathlen += drift_length;
}

void HRSr::Project(TrackState& track, double distance) {
    /**
     * Project particle through field-free region
     * Matches Fortran project() subroutine
     */
    track.x += distance * track.dx;
    track.y += distance * track.dy;
    track.z += distance;
    track.pathlen += distance;
}

} // namespace simc
