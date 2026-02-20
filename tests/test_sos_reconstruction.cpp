// tests/test_sos_reconstruction.cpp
// Phase 5e: SOS Reconstruction Validation Test
// Tests round-trip accuracy: target → focal plane → target

#include "simc/spectrometers/SOS.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>

using namespace simc;

int main() {
    std::cout << "=== SOS RECONSTRUCTION VALIDATION ===" << std::endl;
    std::cout << std::endl;
    
    // Create SOS instance
    SOS sos;
    
    // Load matrices
    if (!sos.LoadMatrices("../data/matrices/sos/forward_cosy.dat",
                          "../data/matrices/sos/recon_cosy.dat")) {
        std::cerr << "Failed to load SOS matrices!" << std::endl;
        return 1;
    }
    
    std::cout << "Matrices loaded successfully." << std::endl;
    std::cout << std::endl;
    
    // Test parameters
    const int n_events = 500;
    const double xp_range = 0.030;    // ±30 mrad
    const double yp_range = 0.030;    // ±30 mrad
    const double delta_range = 15.0;  // ±15% (SOS has large acceptance)
    const double y_range = 0.5;       // ±0.5 cm
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_angle(-1.0, 1.0);
    std::uniform_real_distribution<> dis_y(-y_range, y_range);
    
    // Statistics
    double sum_xp_error_sq = 0.0;
    double sum_yp_error_sq = 0.0;
    double sum_delta_error_sq = 0.0;
    int n_success = 0;
    
    std::cout << "Testing " << n_events << " events:" << std::endl;
    std::cout << "  XP range: ±" << xp_range*1000 << " mrad" << std::endl;
    std::cout << "  YP range: ±" << yp_range*1000 << " mrad" << std::endl;
    std::cout << "  Delta range: ±" << delta_range << "%" << std::endl;
    std::cout << std::endl;
    
    for (int i = 0; i < n_events; ++i) {
        // Generate random target coordinates
        double xp_true = dis_angle(gen) * xp_range;
        double yp_true = dis_angle(gen) * yp_range;
        double delta_true = dis_angle(gen) * delta_range;
        double y_true = dis_y(gen);
        
        // Create track state
        SOS::TrackState track;
        track.x = 0.0;  // Assume beam at x=0
        track.y = y_true;
        track.z = 0.0;
        track.dx = xp_true;  // slopes
        track.dy = yp_true;
        track.delta = delta_true;
        track.p = 1.0 * (1.0 + delta_true/100.0);  // 1 GeV/c central momentum
        track.m2 = 0.000511 * 0.000511;  // electron mass²
        
        // Forward transport through SOS
        if (!sos.Transport(track)) {
            continue;  // Particle didn't make it through
        }
        
        // Get focal plane coordinates
        SOS::FocalPlaneState fp;
        fp.x = sos.GetFocalPlaneX();
        fp.dx = sos.GetFocalPlaneDX();
        fp.y = sos.GetFocalPlaneY();
        fp.dy = sos.GetFocalPlaneDY();
        fp.delta = sos.GetFocalPlaneDelta();
        
        // Reconstruct target coordinates
        SOS::TargetState target;
        target.y = y_true;  // Known from raster
        
        if (!sos.Reconstruct(fp, target)) {
            std::cerr << "Reconstruction failed for event " << i << std::endl;
            continue;
        }
        
        // Calculate errors
        double xp_error = (target.xp - xp_true) * 1000.0;  // mrad
        double yp_error = (target.yp - yp_true) * 1000.0;  // mrad
        double delta_error = target.delta - delta_true;     // %
        
        sum_xp_error_sq += xp_error * xp_error;
        sum_yp_error_sq += yp_error * yp_error;
        sum_delta_error_sq += delta_error * delta_error;
        n_success++;
        
        // Print first few events
        if (i < 5) {
            std::cout << "Event " << i+1 << ":" << std::endl;
            std::cout << "  Target (true):  xp=" << std::setw(8) << std::fixed << std::setprecision(6) << xp_true
                      << ", yp=" << std::setw(8) << yp_true
                      << ", δ=" << std::setw(7) << std::setprecision(2) << delta_true << "%" << std::endl;
            std::cout << "  Focal plane:    x=" << std::setw(8) << std::setprecision(3) << fp.x << " cm"
                      << ", dx=" << std::setw(10) << std::setprecision(6) << fp.dx
                      << ", y=" << std::setw(8) << std::setprecision(3) << fp.y << " cm"
                      << ", dy=" << std::setw(10) << std::setprecision(6) << fp.dy << std::endl;
            std::cout << "  Target (recon): xp=" << std::setw(8) << target.xp
                      << ", yp=" << std::setw(8) << target.yp
                      << ", δ=" << std::setw(7) << std::setprecision(2) << target.delta << "%" << std::endl;
            std::cout << "  Errors:         xp=" << std::setw(7) << std::setprecision(2) << xp_error << " mrad"
                      << ", yp=" << std::setw(7) << yp_error << " mrad"
                      << ", δ=" << std::setw(7) << delta_error << "%" << std::endl;
            std::cout << std::endl;
        }
    }
    
    if (n_success == 0) {
        std::cerr << "No successful events!" << std::endl;
        return 1;
    }
    
    // Calculate RMS errors
    double xp_rms = std::sqrt(sum_xp_error_sq / n_success);
    double yp_rms = std::sqrt(sum_yp_error_sq / n_success);
    double delta_rms = std::sqrt(sum_delta_error_sq / n_success);
    
    std::cout << "=== RESULTS (" << n_success << " successful events) ===" << std::endl;
    std::cout << "RMS Errors:" << std::endl;
    std::cout << "  XP:    " << std::setw(6) << std::fixed << std::setprecision(2) << xp_rms << " mrad";
    if (xp_rms < 5.0) std::cout << " ✓ GOOD";
    std::cout << std::endl;
    
    std::cout << "  YP:    " << std::setw(6) << yp_rms << " mrad";
    if (yp_rms < 1.0) std::cout << " ✓ EXCELLENT";
    else if (yp_rms < 5.0) std::cout << " ✓ GOOD";
    std::cout << std::endl;
    
    std::cout << "  Delta: " << std::setw(6) << delta_rms << " %";
    if (delta_rms < 0.5) std::cout << " ✓ EXCELLENT";
    else if (delta_rms < 1.0) std::cout << " ✓ GOOD";
    std::cout << std::endl;
    
    std::cout << std::endl;
    std::cout << "Acceptance: " << (100.0 * n_success / n_events) << "%" << std::endl;
    std::cout << std::endl;
    
    std::cout << "=== INTERPRETATION ===" << std::endl;
    std::cout << "SOS has large acceptance (±40% delta, ±30 mrad angles)." << std::endl;
    std::cout << "Expected reconstruction accuracy:" << std::endl;
    std::cout << "  - XP, YP: 3-5 mrad RMS (acceptable for large acceptance)" << std::endl;
    std::cout << "  - Delta: < 1% RMS (good)" << std::endl;
    std::cout << std::endl;
    
    return 0;
}
