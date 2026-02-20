/**
 * @file test_xp_delta_coupling.cpp
 * @brief Test XP reconstruction with different deltas
 */

#include "simc/spectrometers/HMS.h"
#include <iostream>
#include <iomanip>

using namespace simc;

int main() {
    std::cout << "\n=== XP-DELTA COUPLING TEST ===\n" << std::endl;
    
    HMS hms;
    
    if (!hms.LoadMatrices("../data/matrices/hms/forward_cosy.dat",
                          "../data/matrices/hms/recon_cosy.dat")) {
        return 1;
    }
    
    hms.SetUseCollimator(false);
    
    struct TestCase {
        double xp;
        double delta;
        std::string name;
    };
    
    std::vector<TestCase> tests = {
        {0.010, 0.0, "xp=+10mrad, δ=0%"},
        {0.010, 5.0, "xp=+10mrad, δ=+5%"},
        {0.010, -5.0, "xp=+10mrad, δ=-5%"},
        {0.0, 5.0, "xp=0, δ=+5%"},
        {-0.010, 0.0, "xp=-10mrad, δ=0%"},
        {-0.010, 5.0, "xp=-10mrad, δ=+5%"},
        {-0.010, -5.0, "xp=-10mrad, δ=-5%"}
    };
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Testing if XP reconstruction is independent of delta:\n" << std::endl;
    
    for (const auto& test : tests) {
        HMS::TrackState track;
        track.x = 0.0;
        track.y = 0.0;
        track.z = 0.0;
        track.dx = test.xp;
        track.dy = 0.0;
        track.delta = test.delta;
        track.p = 8.8 * (1.0 + test.delta / 100.0);
        track.m2 = 0.000511 * 0.000511;
        
        if (!hms.Transport(track)) {
            std::cout << test.name << ": FAILED transport" << std::endl;
            continue;
        }
        
        HMS::FocalPlaneState fp;
        fp.x = hms.GetFocalPlaneX();
        fp.dx = hms.GetFocalPlaneDX();
        fp.y = hms.GetFocalPlaneY();
        fp.dy = hms.GetFocalPlaneDY();
        fp.delta = hms.GetFocalPlaneDelta();
        
        HMS::TargetState target;
        target.y = 0.0;
        
        if (!hms.Reconstruct(fp, target)) {
            std::cout << test.name << ": FAILED reconstruction" << std::endl;
            continue;
        }
        
        double xp_error = (target.xp - test.xp) * 1000.0;  // mrad
        double delta_error = target.delta - test.delta;      // %
        
        std::cout << std::setw(25) << std::left << test.name << " | ";
        std::cout << "FP: x=" << std::setw(8) << std::right << fp.x 
                  << " cm, dx=" << std::setw(9) << fp.dx << " | ";
        std::cout << "Recon: xp=" << std::setw(9) << target.xp 
                  << ", δ=" << std::setw(6) << std::setprecision(2) << target.delta << "% | ";
        std::cout << "Error: xp=" << std::setw(6) << std::setprecision(2) << xp_error 
                  << " mrad, δ=" << std::setw(5) << delta_error << "%" << std::endl;
    }
    
    std::cout << "\n=== EXPECTED BEHAVIOR ===" << std::endl;
    std::cout << "For fixed XP, changing delta should:" << std::endl;
    std::cout << "  1. Change focal plane X (dispersion)" << std::endl;
    std::cout << "  2. Change focal plane dx slightly" << std::endl;
    std::cout << "  3. Reconstruct back to SAME XP (within ~0.5 mrad)" << std::endl;
    std::cout << "  4. Reconstruct correct delta" << std::endl;
    std::cout << "\nIf XP errors are large (>1 mrad) when delta changes," << std::endl;
    std::cout << "there's a problem with XP-delta coupling in reconstruction." << std::endl;
    
    return 0;
}
