/**
 * @file test_hms_recon_regions.cpp  
 * @brief Test reconstruction in different phase space regions
 */

#include "simc/spectrometers/HMS.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <vector>

using namespace simc;

struct RegionTest {
    std::string name;
    double xp_min, xp_max;
    double yp_min, yp_max;
    double delta_min, delta_max;
    
    double xp_rms;
    double yp_rms;
    double delta_rms;
    int n_events;
};

int main() {
    std::cout << "\n=== HMS RECONSTRUCTION BY PHASE SPACE REGION ===\n" << std::endl;
    
    HMS hms;
    
    if (!hms.LoadMatrices("../data/matrices/hms/forward_cosy.dat",
                          "../data/matrices/hms/recon_cosy.dat")) {
        std::cerr << "Failed to load matrices!" << std::endl;
        return 1;
    }
    
    hms.SetUseCollimator(false);
    
    std::random_device rd;
    std::mt19937 gen(12345);
    
    std::uniform_real_distribution<> ytar_dist(-0.5, +0.5);
    
    const double central_p = 8.8;
    const int events_per_region = 200;
    
    // Define test regions
    std::vector<RegionTest> regions = {
        {"Central (±5 mrad, ±2%)",   -0.005, +0.005, -0.005, +0.005, -2.0, +2.0, 0,0,0,0},
        {"Medium angles (±15 mrad, ±5%)", -0.015, +0.015, -0.015, +0.015, -5.0, +5.0, 0,0,0,0},
        {"Large angles (±25 mrad, ±8%)",  -0.025, +0.025, -0.025, +0.025, -8.0, +8.0, 0,0,0,0},
        {"Full acceptance (±30 mrad, ±8%)", -0.030, +0.030, -0.030, +0.030, -8.0, +8.0, 0,0,0,0}
    };
    
    for (auto& region : regions) {
        std::cout << "Testing: " << region.name << std::endl;
        
        std::uniform_real_distribution<> xp_dist(region.xp_min, region.xp_max);
        std::uniform_real_distribution<> yp_dist(region.yp_min, region.yp_max);
        std::uniform_real_distribution<> delta_dist(region.delta_min, region.delta_max);
        
        double xp_sum2 = 0, yp_sum2 = 0, delta_sum2 = 0;
        int accepted = 0;
        
        for (int i = 0; i < events_per_region; ++i) {
            HMS::TrackState track;
            track.x = 0.0;
            track.y = ytar_dist(gen);
            track.z = 0.0;
            track.dx = xp_dist(gen);
            track.dy = yp_dist(gen);
            track.delta = delta_dist(gen);
            track.p = central_p * (1.0 + track.delta / 100.0);
            track.m2 = 0.000511 * 0.000511;
            
            double true_xp = track.dx;
            double true_yp = track.dy;
            double true_delta = track.delta;
            double true_y = track.y;
            
            if (!hms.Transport(track)) continue;
            
            HMS::FocalPlaneState fp;
            fp.x = hms.GetFocalPlaneX();
            fp.dx = hms.GetFocalPlaneDX();
            fp.y = hms.GetFocalPlaneY();
            fp.dy = hms.GetFocalPlaneDY();
            fp.delta = hms.GetFocalPlaneDelta();
            
            HMS::TargetState target;
            target.y = true_y;
            
            if (!hms.Reconstruct(fp, target)) continue;
            
            double xp_diff = target.xp - true_xp;
            double yp_diff = target.yp - true_yp;
            double delta_diff = target.delta - true_delta;
            
            xp_sum2 += xp_diff * xp_diff;
            yp_sum2 += yp_diff * yp_diff;
            delta_sum2 += delta_diff * delta_diff;
            accepted++;
        }
        
        if (accepted > 0) {
            region.xp_rms = std::sqrt(xp_sum2 / accepted);
            region.yp_rms = std::sqrt(yp_sum2 / accepted);
            region.delta_rms = std::sqrt(delta_sum2 / accepted);
            region.n_events = accepted;
            
            std::cout << "  Events: " << accepted << "/" << events_per_region << std::endl;
            std::cout << "  XP RMS:    " << std::fixed << std::setprecision(3) 
                      << region.xp_rms * 1000.0 << " mrad";
            if (region.xp_rms * 1000.0 < 1.0) std::cout << " ✓";
            std::cout << std::endl;
            
            std::cout << "  YP RMS:    " << region.yp_rms * 1000.0 << " mrad";
            if (region.yp_rms * 1000.0 < 1.0) std::cout << " ✓";
            std::cout << std::endl;
            
            std::cout << "  Delta RMS: " << region.delta_rms << " %";
            if (region.delta_rms < 0.5) std::cout << " ✓";
            std::cout << "\n" << std::endl;
        }
    }
    
    std::cout << "=== SUMMARY ===" << std::endl;
    std::cout << "Region                           XP (mrad)  YP (mrad)  Delta (%)  Status" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    
    for (const auto& r : regions) {
        std::cout << std::left << std::setw(32) << r.name;
        std::cout << std::right << std::fixed << std::setprecision(2);
        std::cout << std::setw(10) << r.xp_rms * 1000.0;
        std::cout << std::setw(11) << r.yp_rms * 1000.0;
        std::cout << std::setw(11) << r.delta_rms;
        
        bool xp_good = (r.xp_rms * 1000.0) < 1.0;
        bool yp_good = (r.yp_rms * 1000.0) < 1.0;
        bool delta_good = r.delta_rms < 0.5;
        
        if (xp_good && yp_good && delta_good) {
            std::cout << "   ✓ GOOD";
        } else if ((r.xp_rms * 1000.0) < 3.0 && (r.yp_rms * 1000.0) < 1.0 && r.delta_rms < 1.0) {
            std::cout << "   ~ OK";
        } else {
            std::cout << "   ✗ POOR";
        }
        std::cout << std::endl;
    }
    
    std::cout << "\n=== INTERPRETATION ===" << std::endl;
    std::cout << "Good:  XP<1mrad, YP<1mrad, Delta<0.5%" << std::endl;
    std::cout << "OK:    XP<3mrad, YP<1mrad, Delta<1.0%" << std::endl;
    std::cout << "Poor:  Worse than OK criteria" << std::endl;
    std::cout << "\nNote: Errors naturally increase at larger angles/deltas" << std::endl;
    std::cout << "due to higher-order aberrations." << std::endl;
    
    return 0;
}
