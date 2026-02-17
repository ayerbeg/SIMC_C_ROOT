// tests/test_hrsr_transport.cpp
// HRSr Standalone Transport Test

#include "simc/spectrometers/HRSl.h"
#include <iostream>
#include <iomanip>
#include <random>

using namespace simc;

struct TestConfig {
    std::string name;
    double theta_max;   // mrad
    double delta_max;   // percent
    int nevents;
};

int main() {
    std::cout << "==========================================\n";
    std::cout << "  HRSl Standalone Transport Test\n";
    std::cout << "==========================================\n\n";
    
    // Initialize HRSl
    HRSl hrsl;
    
    std::cout << "Loading HRSl matrix files...\n";
    if (!hrsl.LoadMatrices("../data/matrices/hrsl/hrs_forward_cosy.dat",
                           "../data/matrices/hrsl/hrs_recon_cosy.dat")) {
        std::cerr << "ERROR: Failed to load matrix files!\n";
        return 1;
    }
    std::cout << "HRSl initialized successfully.\n\n";
    
    // Test configurations
    std::vector<TestConfig> configs = {
        {"Wide Open (±80 mrad, ±4.5%)", 80.0, 4.5, 10000},
        {"Medium (±60 mrad, ±3.0%)", 60.0, 3.0, 10000},
        {"Tight (±40 mrad, ±2.0%)", 40.0, 2.0, 10000}
    };
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Run tests
    for (const auto& config : configs) {
        std::cout << "==========================================\n";
        std::cout << "Test: " << config.name << "\n";
        std::cout << "==========================================\n";
        
        hrsl.ResetStats();
        
        std::uniform_real_distribution<> angle_dist(-config.theta_max, config.theta_max);
        std::uniform_real_distribution<> delta_dist(-config.delta_max, config.delta_max);
        std::uniform_real_distribution<> target_dist(-0.1, 0.1);  // ±1mm target offset
        
        int accepted = 0;
        
        for (int i = 0; i < config.nevents; ++i) {
            HRSl::TrackState track;
            
            // Generate random particle at target
            track.x = 0.0;  // cm
            track.y = target_dist(gen);  // cm
            track.z = 0.0;
            track.dx = angle_dist(gen) / 1000.0;  // mrad → dimensionless
            track.dy = angle_dist(gen) / 1000.0;
            track.delta = delta_dist(gen);  // percent
            track.p = 1.0 * (1.0 + track.delta/100.0);  // 1 GeV/c central momentum (typical for HRS)
            track.m2 = 0.13957 * 0.13957;  // pion mass squared
            
            // Transport through HRSr
            if (hrsl.Transport(track)) {
                ++accepted;
            }
        }
        
        auto stats = hrsl.GetStats();
        double acceptance = 100.0 * accepted / config.nevents;
        
        std::cout << "\nResults:\n";
        std::cout << "  Events:      " << config.nevents << "\n";
        std::cout << "  Accepted:    " << accepted << "\n";
        std::cout << "  Acceptance:  " << std::fixed << std::setprecision(1) 
                  << acceptance << "%\n";
        std::cout << "\nRejection Breakdown:\n";
        std::cout << "  Collimator:    " << stats.lSTOP_slit << "\n";
        std::cout << "  Q1 entrance:   " << stats.lSTOP_Q1_in << "\n";
        std::cout << "  Q1 mid:        " << stats.lSTOP_Q1_mid << "\n";
        std::cout << "  Q1 exit:       " << stats.lSTOP_Q1_out << "\n";
        std::cout << "  Q2 entrance:   " << stats.lSTOP_Q2_in << "\n";
        std::cout << "  Q2 mid:        " << stats.lSTOP_Q2_mid << "\n";
        std::cout << "  Q2 exit:       " << stats.lSTOP_Q2_out << "\n";
        std::cout << "  D1 entrance:   " << stats.lSTOP_D1_in << "\n";
        std::cout << "  D1 exit:       " << stats.lSTOP_D1_out << "\n";
        std::cout << "  Q3 entrance:   " << stats.lSTOP_Q3_in << "\n";
        std::cout << "  Q3 mid:        " << stats.lSTOP_Q3_mid << "\n";
        std::cout << "  Q3 exit:       " << stats.lSTOP_Q3_out << "\n";
        std::cout << "  Made to hut:   " << stats.lSTOP_hut << "\n";
        std::cout << "\n";
    }
    
    std::cout << "==========================================\n";
    std::cout << "HRSl Standalone Test Complete!\n";
    std::cout << "==========================================\n\n";
    
    return 0;
}
