/**
 * @file test_hms_reconstruction.cpp
 * @brief Test HMS Phase 5e reconstruction implementation
 * 
 * PHASE 5e TEST PROGRAM
 * 
 * Tests the complete forward → focal plane → reconstruction → target chain:
 * 1. Generate random target coordinates
 * 2. Transport through HMS (forward)
 * 3. Extract focal plane coordinates
 * 4. Reconstruct back to target
 * 5. Compare reconstructed vs. original target coordinates
 * 
 * Expected results:
 * - Reconstruction should reproduce target angles to ~1 mrad precision
 * - Reconstruction should reproduce target delta to ~0.1% precision
 * - Any deviations indicate bugs in reconstruction matrix or unit handling
 */

#include "simc/spectrometers/HMS.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <vector>

using namespace simc;

struct ReconstructionTest {
    // Original target coordinates
    double tgt_xp_true;
    double tgt_yp_true;
    double tgt_delta_true;
    
    // Focal plane coordinates (from transport)
    double fp_x;
    double fp_dx;
    double fp_y;
    double fp_dy;
    double fp_delta;
    
    // Reconstructed target coordinates
    double tgt_xp_recon;
    double tgt_yp_recon;
    double tgt_delta_recon;
    
    // Differences
    double xp_diff;
    double yp_diff;
    double delta_diff;
};

void PrintTestHeader() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "HMS PHASE 5e RECONSTRUCTION TEST" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

void PrintStats(const std::vector<ReconstructionTest>& tests) {
    if (tests.empty()) return;
    
    // Calculate statistics
    double xp_mean = 0.0, xp_rms = 0.0;
    double yp_mean = 0.0, yp_rms = 0.0;
    double delta_mean = 0.0, delta_rms = 0.0;
    
    for (const auto& t : tests) {
        xp_mean += t.xp_diff;
        yp_mean += t.yp_diff;
        delta_mean += t.delta_diff;
        
        xp_rms += t.xp_diff * t.xp_diff;
        yp_rms += t.yp_diff * t.yp_diff;
        delta_rms += t.delta_diff * t.delta_diff;
    }
    
    int n = tests.size();
    xp_mean /= n;
    yp_mean /= n;
    delta_mean /= n;
    
    xp_rms = std::sqrt(xp_rms / n);
    yp_rms = std::sqrt(yp_rms / n);
    delta_rms = std::sqrt(delta_rms / n);
    
    std::cout << "\n=== RECONSTRUCTION ACCURACY ===" << std::endl;
    std::cout << "Events tested: " << n << std::endl;
    std::cout << "\nTarget XP (slopes):" << std::endl;
    std::cout << "  Mean difference: " << std::scientific << std::setprecision(3) 
              << xp_mean << " (dimensionless)" << std::endl;
    std::cout << "  RMS difference:  " << xp_rms << " (dimensionless)" << std::endl;
    std::cout << "  RMS in mrad:     " << std::fixed << std::setprecision(3) 
              << xp_rms * 1000.0 << " mrad" << std::endl;
    
    std::cout << "\nTarget YP (slopes):" << std::endl;
    std::cout << "  Mean difference: " << std::scientific << std::setprecision(3) 
              << yp_mean << " (dimensionless)" << std::endl;
    std::cout << "  RMS difference:  " << yp_rms << " (dimensionless)" << std::endl;
    std::cout << "  RMS in mrad:     " << std::fixed << std::setprecision(3) 
              << yp_rms * 1000.0 << " mrad" << std::endl;
    
    std::cout << "\nTarget Delta:" << std::endl;
    std::cout << "  Mean difference: " << std::scientific << std::setprecision(3) 
              << delta_mean << " %" << std::endl;
    std::cout << "  RMS difference:  " << delta_rms << " %" << std::endl;
    
    // Quality assessment
    std::cout << "\n=== QUALITY ASSESSMENT ===" << std::endl;
    bool xp_good = (xp_rms * 1000.0) < 1.0;  // < 1 mrad
    bool yp_good = (yp_rms * 1000.0) < 1.0;  // < 1 mrad
    bool delta_good = delta_rms < 0.1;        // < 0.1 %
    
    std::cout << "XP reconstruction:    " << (xp_good ? "✓ GOOD" : "✗ POOR") 
              << " (RMS = " << std::fixed << std::setprecision(3) 
              << xp_rms * 1000.0 << " mrad)" << std::endl;
    std::cout << "YP reconstruction:    " << (yp_good ? "✓ GOOD" : "✗ POOR")
              << " (RMS = " << yp_rms * 1000.0 << " mrad)" << std::endl;
    std::cout << "Delta reconstruction: " << (delta_good ? "✓ GOOD" : "✗ POOR")
              << " (RMS = " << delta_rms << " %)" << std::endl;
    
    if (xp_good && yp_good && delta_good) {
        std::cout << "\n✓ ALL TESTS PASSED - Reconstruction is working correctly!" << std::endl;
    } else {
        std::cout << "\n✗ TESTS FAILED - Check unit conversions and matrix application!" << std::endl;
    }
}

int main() {
    PrintTestHeader();
    
    // Initialize HMS
    HMS hms;
    
    // Load matrix files
    std::string forward_file = "../data/matrices/hms/forward_cosy.dat";
    std::string recon_file = "../data/matrices/hms/recon_cosy.dat";
    
    std::cout << "Loading HMS matrices..." << std::endl;
    if (!hms.LoadMatrices(forward_file, recon_file)) {
        std::cerr << "ERROR: Failed to load HMS matrices!" << std::endl;
        std::cerr << "  Forward: " << forward_file << std::endl;
        std::cerr << "  Recon:   " << recon_file << std::endl;
        return 1;
    }
    std::cout << "✓ Matrices loaded successfully\n" << std::endl;
    
    // Configuration
    const int num_events = 1000;
    const double central_momentum = 8.8;  // GeV/c
    
    // Use wide-open mode (no collimator) for maximum statistics
    hms.SetUseCollimator(false);
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(12345);  // Fixed seed for reproducibility
    
    // Realistic HMS acceptance limits (±30 mrad angles, ±8% delta)
    std::uniform_real_distribution<> xptar_dist(-0.030, +0.030);  // slopes
    std::uniform_real_distribution<> yptar_dist(-0.030, +0.030);  // slopes
    std::uniform_real_distribution<> delta_dist(-8.0, +8.0);      // %
    std::uniform_real_distribution<> ytar_dist(-0.5, +0.5);       // cm (raster)
    
    std::vector<ReconstructionTest> results;
    
    std::cout << "Running " << num_events << " reconstruction tests..." << std::endl;
    std::cout << "  Angle range: ±30 mrad" << std::endl;
    std::cout << "  Delta range: ±8 %" << std::endl;
    std::cout << "  Target Y:    ±0.5 cm (raster)\n" << std::endl;
    
    int accepted = 0;
    
    for (int i = 0; i < num_events; ++i) {
        HMS::TrackState track;
        
        // Generate random target coordinates
        track.x = 0.0;  // Always at beam center
        track.y = ytar_dist(gen);  // Raster position (cm)
        track.z = 0.0;  // Target position
        
        track.dx = xptar_dist(gen);  // Target slope (dimensionless)
        track.dy = yptar_dist(gen);  // Target slope (dimensionless)
        track.delta = delta_dist(gen);  // Momentum deviation (%)
        
        track.p = central_momentum * (1.0 + track.delta / 100.0);
        track.m2 = 0.000511 * 0.000511;  // Electron mass squared
        
        // Store true target values
        double tgt_xp_true = track.dx;
        double tgt_yp_true = track.dy;
        double tgt_y_true = track.y;
        double tgt_delta_true = track.delta;
        
        // Forward transport through HMS
        bool pass = hms.Transport(track);
        
        if (!pass) {
            continue;  // Skip events that don't make it through
        }
        
        accepted++;
        
        // Get focal plane coordinates (stored by Transport())
        HMS::FocalPlaneState fp;
        fp.x = hms.GetFocalPlaneX();      // cm
        fp.dx = hms.GetFocalPlaneDX();    // slope (dimensionless)
        fp.y = hms.GetFocalPlaneY();      // cm
        fp.dy = hms.GetFocalPlaneDY();    // slope (dimensionless)
        fp.delta = hms.GetFocalPlaneDelta();  // %
        
        // Reconstruct back to target
        HMS::TargetState target;
        target.y = tgt_y_true;  // Known from raster (required input)
        
        if (!hms.Reconstruct(fp, target)) {
            std::cerr << "ERROR: Reconstruction failed for event " << i << std::endl;
            continue;
        }
        
        // Calculate differences
        ReconstructionTest test;
        test.tgt_xp_true = tgt_xp_true;
        test.tgt_yp_true = tgt_yp_true;
        test.tgt_delta_true = tgt_delta_true;
        
        test.fp_x = fp.x;
        test.fp_dx = fp.dx;
        test.fp_y = fp.y;
        test.fp_dy = fp.dy;
        test.fp_delta = fp.delta;
        
        test.tgt_xp_recon = target.xp;
        test.tgt_yp_recon = target.yp;
        test.tgt_delta_recon = target.delta;
        
        test.xp_diff = target.xp - tgt_xp_true;
        test.yp_diff = target.yp - tgt_yp_true;
        test.delta_diff = target.delta - tgt_delta_true;
        
        results.push_back(test);
        
        // Print first few events in detail
        if (accepted <= 5) {
            std::cout << "\n--- Event " << accepted << " ---" << std::endl;
            std::cout << "Target (true):  xp=" << std::fixed << std::setprecision(6) 
                      << tgt_xp_true << ", yp=" << tgt_yp_true 
                      << ", δ=" << std::setprecision(2) << tgt_delta_true << "%" << std::endl;
            std::cout << "Focal plane:    x=" << std::setprecision(2) << fp.x 
                      << " cm, dx=" << std::setprecision(6) << fp.dx
                      << ", y=" << std::setprecision(2) << fp.y 
                      << " cm, dy=" << std::setprecision(6) << fp.dy << std::endl;
            std::cout << "Target (recon): xp=" << std::setprecision(6) << target.xp 
                      << ", yp=" << target.yp 
                      << ", δ=" << std::setprecision(2) << target.delta << "%" << std::endl;
            std::cout << "Differences:    Δxp=" << std::scientific << std::setprecision(2) 
                      << test.xp_diff << " (" << std::fixed << std::setprecision(3) 
                      << test.xp_diff*1000.0 << " mrad)"
                      << ", Δyp=" << std::scientific << std::setprecision(2) << test.yp_diff 
                      << " (" << std::fixed << std::setprecision(3) << test.yp_diff*1000.0 << " mrad)"
                      << ", Δδ=" << std::setprecision(3) << test.delta_diff << "%" << std::endl;
        }
    }
    
    std::cout << "\n=== TRANSPORT SUMMARY ===" << std::endl;
    std::cout << "Events generated:  " << num_events << std::endl;
    std::cout << "Events accepted:   " << accepted << std::endl;
    std::cout << "Acceptance:        " << std::fixed << std::setprecision(1) 
              << (100.0 * accepted / num_events) << "%" << std::endl;
    std::cout << "Events tested:     " << results.size() << std::endl;
    
    // Print statistics
    PrintStats(results);
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "TEST COMPLETE" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    return 0;
}
