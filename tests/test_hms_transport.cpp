/**
 * @file test_hms_transport.cpp
 * @brief Standalone test for HMS spectrometer transport
 * 
 * Tests HMS transport with random angles and momentum deviation.
 * Should give ~80-90% acceptance with wide-open mode (no collimator).
 */

#include "simc/spectrometers/HMS.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

using namespace simc;

int main() {
    std::cout << "=== HMS Transport Standalone Test ===" << std::endl;
    std::cout << std::endl;
    
    // Initialize HMS
    HMS hms;
    
    // Load matrix files
    std::string forward_file = "../data/matrices/hms/forward_cosy.dat";
    std::string recon_file = "../data/matrices/hms/recon_cosy.dat";
    
    if (!hms.LoadMatrices(forward_file, recon_file)) {
        std::cerr << "Failed to load HMS matrices!" << std::endl;
        return 1;
    }
    
    std::cout << "Loaded HMS matrices successfully." << std::endl;
    std::cout << std::endl;
    
    // Configuration
    const int num_events = 10000;
    const double central_momentum = 8.8;  // GeV/c (example electron momentum)
    
    // Test configurations
    struct TestConfig {
        std::string name;
        bool use_collimator;
        double delta_min;
        double delta_max;
        double xptar_min;
        double xptar_max;
        double yptar_min;
        double yptar_max;
    };
    
    std::vector<TestConfig> tests = {
        {
            "Wide Open (No Collimator)",
            false,  // no collimator
            -10.0,  // delta min (%)
            +10.0,  // delta max (%)
            -60.0,  // xptar min (mrad)
            +60.0,  // xptar max (mrad)
            -60.0,  // yptar min (mrad)
            +60.0   // yptar max (mrad)
        },
        {
            "With HMS-100 Collimator",
            true,   // use collimator
            -10.0,
            +10.0,
            -60.0,
            +60.0,
            -60.0,
            +60.0
        },
        {
            "Tight Angles (±30 mrad)",
            false,
            -10.0,
            +10.0,
            -30.0,
            +30.0,
            -30.0,
            +30.0
        }
    };
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Run tests
    for (const auto& test : tests) {
        std::cout << "===========================================" << std::endl;
        std::cout << "Test: " << test.name << std::endl;
        std::cout << "===========================================" << std::endl;
        std::cout << "  Collimator: " << (test.use_collimator ? "YES" : "NO") << std::endl;
        std::cout << "  Delta range: " << test.delta_min << " to " << test.delta_max << " %" << std::endl;
        std::cout << "  Angle range (X): " << test.xptar_min << " to " << test.xptar_max << " mrad" << std::endl;
        std::cout << "  Angle range (Y): " << test.yptar_min << " to " << test.yptar_max << " mrad" << std::endl;
        std::cout << std::endl;
        
        // Configure HMS
        hms.SetUseCollimator(test.use_collimator);
        hms.ResetStatistics();
        
        // Uniform distributions
        std::uniform_real_distribution<> delta_dist(test.delta_min, test.delta_max);
        std::uniform_real_distribution<> xptar_dist(test.xptar_min / 1000.0, test.xptar_max / 1000.0);  // mrad -> rad
        std::uniform_real_distribution<> yptar_dist(test.yptar_min / 1000.0, test.yptar_max / 1000.0);
        
        // Transport events
        int accepted = 0;
        
        for (int i = 0; i < num_events; ++i) {
            HMS::TrackState track;
            
            // Initial position (at target)
            track.x = 0.0;  // cm
            track.y = 0.0;  // cm
            track.z = 0.0;  // cm (target position)
            
            // Random angles (target coordinates)
            track.dx = xptar_dist(gen);  // dX/dZ (slope)
            track.dy = yptar_dist(gen);  // dY/dZ (slope)
            
            // Random momentum deviation
            track.delta = delta_dist(gen);  // %
            
            // Calculate momentum
            track.p = central_momentum * (1.0 + track.delta / 100.0);  // GeV/c
            
            // Assume electron mass
            track.m2 = 0.000511 * 0.000511;  // (GeV)^2
            
            // Transport through HMS
            bool pass = hms.Transport(track);
            
            if (pass) {
                accepted++;
            }
            
            // Debug first few events
            if (i < 5) {
                std::cout << "Event " << i << ": ";
                std::cout << "dx=" << std::setw(6) << std::fixed << std::setprecision(1) 
                         << track.dx*1000.0 << " mrad, ";
                std::cout << "dy=" << std::setw(6) << track.dy*1000.0 << " mrad, ";
                std::cout << "delta=" << std::setw(6) << std::setprecision(2) << track.delta << "% ";
                std::cout << (pass ? "PASS" : "FAIL") << std::endl;
            }
        }
        
        // Print results
        std::cout << std::endl;
        std::cout << "Results:" << std::endl;
        std::cout << "  Events tried:    " << num_events << std::endl;
        std::cout << "  Events accepted: " << accepted << std::endl;
        std::cout << "  Acceptance:      " << std::fixed << std::setprecision(1) 
                  << (100.0 * accepted / num_events) << "%" << std::endl;
        std::cout << std::endl;
        
        // Print detailed statistics
        const auto& stats = hms.GetStatistics();
        std::cout << "Rejection Statistics:" << std::endl;
        
        if (test.use_collimator) {
            int coll_total = stats.hSTOP_slit_hor + stats.hSTOP_slit_vert + 
                           stats.hSTOP_slit_oct + stats.hSTOP_coll;
            std::cout << "  Collimator:   " << coll_total << std::endl;
            std::cout << "    Horizontal: " << stats.hSTOP_slit_hor << std::endl;
            std::cout << "    Vertical:   " << stats.hSTOP_slit_vert << std::endl;
            std::cout << "    Octagonal:  " << stats.hSTOP_slit_oct << std::endl;
        }
        
        int q1_total = stats.hSTOP_Q1_in + stats.hSTOP_Q1_mid + stats.hSTOP_Q1_out;
        int q2_total = stats.hSTOP_Q2_in + stats.hSTOP_Q2_mid + stats.hSTOP_Q2_out;
        int q3_total = stats.hSTOP_Q3_in + stats.hSTOP_Q3_mid + stats.hSTOP_Q3_out;
        int d1_total = stats.hSTOP_D1_in + stats.hSTOP_D1_out;
        
        std::cout << "  Q1:           " << q1_total << std::endl;
        std::cout << "    Entrance:   " << stats.hSTOP_Q1_in << std::endl;
        std::cout << "    Mid:        " << stats.hSTOP_Q1_mid << std::endl;
        std::cout << "    Exit:       " << stats.hSTOP_Q1_out << std::endl;
        
        std::cout << "  Q2:           " << q2_total << std::endl;
        std::cout << "    Entrance:   " << stats.hSTOP_Q2_in << std::endl;
        std::cout << "    Mid:        " << stats.hSTOP_Q2_mid << std::endl;
        std::cout << "    Exit:       " << stats.hSTOP_Q2_out << std::endl;
        
        std::cout << "  Q3:           " << q3_total << std::endl;
        std::cout << "    Entrance:   " << stats.hSTOP_Q3_in << std::endl;
        std::cout << "    Mid:        " << stats.hSTOP_Q3_mid << std::endl;
        std::cout << "    Exit:       " << stats.hSTOP_Q3_out << std::endl;
        
        std::cout << "  Dipole:       " << d1_total << std::endl;
        std::cout << "    Entrance:   " << stats.hSTOP_D1_in << std::endl;
        std::cout << "    Exit+Pipes: " << stats.hSTOP_D1_out << std::endl;
        
        std::cout << "  Hut:          " << stats.hSTOP_hut << std::endl;
        
        std::cout << std::endl;
    }
    
    std::cout << "=== HMS Transport Test Complete ===" << std::endl;
    
    return 0;
}
