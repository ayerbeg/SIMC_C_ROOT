// tests/test_sos_transport.cpp
// Standalone test for SOS spectrometer - Phase 5c.5

#include "simc/spectrometers/SOS.h"
#include <iostream>
#include <iomanip>
#include <random>

using namespace simc;

int main() {
    std::cout << "\n";
    std::cout << "==========================================\n";
    std::cout << "   SOS Standalone Transport Test\n";
    std::cout << "==========================================\n";
    std::cout << "\n";
    
    // Initialize SOS
    SOS sos;
    
    std::cout << "Loading SOS matrix files..." << std::endl;
    if (!sos.LoadMatrices("../data/matrices/sos/forward_cosy.dat",
                          "../data/matrices/sos/recon_cosy.dat")) {
        std::cerr << "Failed to load SOS matrices!" << std::endl;
        return 1;
    }
    std::cout << "SOS initialized successfully.\n\n";
    
    // Random number generator
    std::mt19937 rng(42);  // Fixed seed for reproducibility
    std::uniform_real_distribution<double> angle_dist(-0.08, 0.08);  // ±80 mrad
    std::uniform_real_distribution<double> delta_dist(-40.0, 40.0);  // ±40% (SOS large acceptance!)
    
    // Test configurations
    struct TestConfig {
        std::string name;
        double angle_range;  // mrad
        double delta_min;
        double delta_max;
        int nevents;
    };
    
    std::vector<TestConfig> configs = {
        {"Wide Open (±80 mrad, ±40%)", 80.0, -40.0, 40.0, 10000},
        {"Medium (±60 mrad, ±30%)", 60.0, -30.0, 30.0, 10000},
        {"Tight (±40 mrad, ±20%)", 40.0, -20.0, 20.0, 10000}
    };
    
    for (const auto& config : configs) {
        std::cout << "==========================================\n";
        std::cout << "Test: " << config.name << "\n";
        std::cout << "==========================================\n";
        
        sos.ResetStats();
        
        int accepted = 0;
        
        for (int i = 0; i < config.nevents; ++i) {
            SOS::TrackState track;
            
            // Initial position at target
            track.x = 0.0;
            track.y = 0.0;
            track.z = 0.0;
            
            // Random angles within config range
            double angle_rad = config.angle_range / 1000.0;  // mrad to rad
            std::uniform_real_distribution<double> angle_gen(-angle_rad, angle_rad);
            track.dx = angle_gen(rng);
            track.dy = angle_gen(rng);
            
            // Random momentum deviation
            std::uniform_real_distribution<double> delta_gen(config.delta_min, config.delta_max);
            track.delta = delta_gen(rng);
            
            // Typical SOS momentum and particle (pion)
            track.p = 0.8;  // GeV/c (mid-range for SOS)
            track.m2 = 0.13957 * 0.13957;  // Pion mass squared
            
            // Transport through SOS
            if (sos.Transport(track)) {
                ++accepted;
            }
        }
        
        auto stats = sos.GetStats();
        double acceptance = 100.0 * accepted / config.nevents;
        
        std::cout << "\nResults:\n";
        std::cout << "  Events:      " << config.nevents << "\n";
        std::cout << "  Accepted:    " << accepted << "\n";
        std::cout << "  Acceptance:  " << std::fixed << std::setprecision(1) 
                  << acceptance << "%\n";
        std::cout << "\nRejection Breakdown:\n";
        std::cout << "  Collimator:    " << stats.sSTOP_slit << "\n";
        std::cout << "  Quad entrance: " << stats.sSTOP_quad_in << "\n";
        std::cout << "  Quad mid:      " << stats.sSTOP_quad_mid << "\n";
        std::cout << "  Quad exit:     " << stats.sSTOP_quad_out << "\n";
        std::cout << "  BM1 entrance:  " << stats.sSTOP_bm1_in << "\n";
        std::cout << "  BM1 exit:      " << stats.sSTOP_bm1_out << "\n";
        std::cout << "  BM2 entrance:  " << stats.sSTOP_bm2_in << "\n";
        std::cout << "  BM2 exit:      " << stats.sSTOP_bm2_out << "\n";
        std::cout << "  Made to hut:   " << stats.sSTOP_hut << "\n";
        std::cout << "\n";
    }
    
    std::cout << "==========================================\n";
    std::cout << "SOS Standalone Test Complete!\n";
    std::cout << "==========================================\n\n";
    
    return 0;
}
