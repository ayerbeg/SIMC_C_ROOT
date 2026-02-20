/**
 * @file test_hms_coordinate_check.cpp
 * @brief Diagnostic test for coordinate system issues
 * 
 * Tests specific simple cases to identify sign/coordinate problems
 */

#include "simc/spectrometers/HMS.h"
#include <iostream>
#include <iomanip>

using namespace simc;

int main() {
    std::cout << "\n=== HMS COORDINATE SYSTEM DIAGNOSTIC ===\n" << std::endl;
    
    HMS hms;
    
    if (!hms.LoadMatrices("../data/matrices/hms/forward_cosy.dat",
                          "../data/matrices/hms/recon_cosy.dat")) {
        std::cerr << "Failed to load matrices!" << std::endl;
        return 1;
    }
    
    hms.SetUseCollimator(false);
    
    std::cout << "Testing simple cases to check coordinate conventions:\n" << std::endl;
    
    // Test 1: Zero angles, zero delta
    {
        std::cout << "TEST 1: Zero angles, zero delta (central ray)" << std::endl;
        HMS::TrackState track;
        track.x = 0.0;
        track.y = 0.0;
        track.z = 0.0;
        track.dx = 0.0;  // Zero slope
        track.dy = 0.0;  // Zero slope
        track.delta = 0.0;
        track.p = 8.8;
        track.m2 = 0.000511 * 0.000511;
        
        double true_xp = track.dx;
        double true_yp = track.dy;
        double true_delta = track.delta;
        
        if (hms.Transport(track)) {
            HMS::FocalPlaneState fp;
            fp.x = hms.GetFocalPlaneX();
            fp.dx = hms.GetFocalPlaneDX();
            fp.y = hms.GetFocalPlaneY();
            fp.dy = hms.GetFocalPlaneDY();
            fp.delta = hms.GetFocalPlaneDelta();
            
            HMS::TargetState target;
            target.y = 0.0;
            hms.Reconstruct(fp, target);
            
            std::cout << "  Target (true):  xp=" << std::fixed << std::setprecision(6) << true_xp
                      << ", yp=" << true_yp << ", δ=" << std::setprecision(2) << true_delta << "%" << std::endl;
            std::cout << "  Focal plane:    x=" << std::setprecision(3) << fp.x << " cm, dx=" 
                      << std::setprecision(6) << fp.dx << ", y=" << std::setprecision(3) << fp.y 
                      << " cm, dy=" << std::setprecision(6) << fp.dy << std::endl;
            std::cout << "  Target (recon): xp=" << std::setprecision(6) << target.xp
                      << ", yp=" << target.yp << ", δ=" << std::setprecision(2) << target.delta << "%" << std::endl;
            std::cout << "  EXPECT: All reconstructed values ≈ 0" << std::endl;
        } else {
            std::cout << "  FAILED transport!" << std::endl;
        }
        std::cout << std::endl;
    }
    
    // Test 2: Positive X angle only
    {
        std::cout << "TEST 2: Positive X angle (+10 mrad), zero Y angle, zero delta" << std::endl;
        HMS::TrackState track;
        track.x = 0.0;
        track.y = 0.0;
        track.z = 0.0;
        track.dx = 0.010;  // +10 mrad as slope
        track.dy = 0.0;
        track.delta = 0.0;
        track.p = 8.8;
        track.m2 = 0.000511 * 0.000511;
        
        double true_xp = track.dx;
        double true_yp = track.dy;
        double true_delta = track.delta;
        
        if (hms.Transport(track)) {
            HMS::FocalPlaneState fp;
            fp.x = hms.GetFocalPlaneX();
            fp.dx = hms.GetFocalPlaneDX();
            fp.y = hms.GetFocalPlaneY();
            fp.dy = hms.GetFocalPlaneDY();
            fp.delta = hms.GetFocalPlaneDelta();
            
            HMS::TargetState target;
            target.y = 0.0;
            hms.Reconstruct(fp, target);
            
            std::cout << "  Target (true):  xp=" << std::fixed << std::setprecision(6) << true_xp
                      << ", yp=" << true_yp << ", δ=" << std::setprecision(2) << true_delta << "%" << std::endl;
            std::cout << "  Focal plane:    x=" << std::setprecision(3) << fp.x << " cm, dx=" 
                      << std::setprecision(6) << fp.dx << ", y=" << std::setprecision(3) << fp.y 
                      << " cm, dy=" << std::setprecision(6) << fp.dy << std::endl;
            std::cout << "  Target (recon): xp=" << std::setprecision(6) << target.xp
                      << ", yp=" << target.yp << ", δ=" << std::setprecision(2) << target.delta << "%" << std::endl;
            std::cout << "  EXPECT: Reconstructed xp ≈ +0.010, yp ≈ 0, δ ≈ 0" << std::endl;
            std::cout << "  ERROR: xp_diff = " << (target.xp - true_xp)*1000.0 << " mrad" << std::endl;
        } else {
            std::cout << "  FAILED transport!" << std::endl;
        }
        std::cout << std::endl;
    }
    
    // Test 3: Negative X angle only
    {
        std::cout << "TEST 3: Negative X angle (-10 mrad), zero Y angle, zero delta" << std::endl;
        HMS::TrackState track;
        track.x = 0.0;
        track.y = 0.0;
        track.z = 0.0;
        track.dx = -0.010;  // -10 mrad as slope
        track.dy = 0.0;
        track.delta = 0.0;
        track.p = 8.8;
        track.m2 = 0.000511 * 0.000511;
        
        double true_xp = track.dx;
        double true_yp = track.dy;
        double true_delta = track.delta;
        
        if (hms.Transport(track)) {
            HMS::FocalPlaneState fp;
            fp.x = hms.GetFocalPlaneX();
            fp.dx = hms.GetFocalPlaneDX();
            fp.y = hms.GetFocalPlaneY();
            fp.dy = hms.GetFocalPlaneDY();
            fp.delta = hms.GetFocalPlaneDelta();
            
            HMS::TargetState target;
            target.y = 0.0;
            hms.Reconstruct(fp, target);
            
            std::cout << "  Target (true):  xp=" << std::fixed << std::setprecision(6) << true_xp
                      << ", yp=" << true_yp << ", δ=" << std::setprecision(2) << true_delta << "%" << std::endl;
            std::cout << "  Focal plane:    x=" << std::setprecision(3) << fp.x << " cm, dx=" 
                      << std::setprecision(6) << fp.dx << ", y=" << std::setprecision(3) << fp.y 
                      << " cm, dy=" << std::setprecision(6) << fp.dy << std::endl;
            std::cout << "  Target (recon): xp=" << std::setprecision(6) << target.xp
                      << ", yp=" << target.yp << ", δ=" << std::setprecision(2) << target.delta << "%" << std::endl;
            std::cout << "  EXPECT: Reconstructed xp ≈ -0.010, yp ≈ 0, δ ≈ 0" << std::endl;
            std::cout << "  ERROR: xp_diff = " << (target.xp - true_xp)*1000.0 << " mrad" << std::endl;
        } else {
            std::cout << "  FAILED transport!" << std::endl;
        }
        std::cout << std::endl;
    }
    
    // Test 4: Positive delta only
    {
        std::cout << "TEST 4: Zero angles, positive delta (+5%)" << std::endl;
        HMS::TrackState track;
        track.x = 0.0;
        track.y = 0.0;
        track.z = 0.0;
        track.dx = 0.0;
        track.dy = 0.0;
        track.delta = 5.0;
        track.p = 8.8 * 1.05;
        track.m2 = 0.000511 * 0.000511;
        
        double true_xp = track.dx;
        double true_yp = track.dy;
        double true_delta = track.delta;
        
        if (hms.Transport(track)) {
            HMS::FocalPlaneState fp;
            fp.x = hms.GetFocalPlaneX();
            fp.dx = hms.GetFocalPlaneDX();
            fp.y = hms.GetFocalPlaneY();
            fp.dy = hms.GetFocalPlaneDY();
            fp.delta = hms.GetFocalPlaneDelta();
            
            HMS::TargetState target;
            target.y = 0.0;
            hms.Reconstruct(fp, target);
            
            std::cout << "  Target (true):  xp=" << std::fixed << std::setprecision(6) << true_xp
                      << ", yp=" << true_yp << ", δ=" << std::setprecision(2) << true_delta << "%" << std::endl;
            std::cout << "  Focal plane:    x=" << std::setprecision(3) << fp.x << " cm, dx=" 
                      << std::setprecision(6) << fp.dx << ", y=" << std::setprecision(3) << fp.y 
                      << " cm, dy=" << std::setprecision(6) << fp.dy << std::endl;
            std::cout << "  Target (recon): xp=" << std::setprecision(6) << target.xp
                      << ", yp=" << target.yp << ", δ=" << std::setprecision(2) << target.delta << "%" << std::endl;
            std::cout << "  EXPECT: Reconstructed xp ≈ 0, yp ≈ 0, δ ≈ +5.0%" << std::endl;
            std::cout << "  ERROR: delta_diff = " << (target.delta - true_delta) << "%" << std::endl;
        } else {
            std::cout << "  FAILED transport!" << std::endl;
        }
        std::cout << std::endl;
    }
    
    std::cout << "=== DIAGNOSTIC COMPLETE ===\n" << std::endl;
    std::cout << "If reconstruction is working correctly:" << std::endl;
    std::cout << "  - All tests should show errors < 0.1 mrad for angles" << std::endl;
    std::cout << "  - Test 4 might show larger delta error (energy loss effects)" << std::endl;
    std::cout << "  - Signs should match (positive input → positive output)" << std::endl;
    
    return 0;
}
