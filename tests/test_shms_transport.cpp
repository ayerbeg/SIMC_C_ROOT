#include "simc/spectrometers/SHMS.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace simc;

int main() {
    std::cout << "=== SHMS Transport Test ===" << std::endl << std::endl;
    
    // Create SHMS instance
    SHMS shms;
    
    // Load matrix files
    std::cout << "Loading matrix files..." << std::endl;
    bool loaded = shms.LoadMatrices(
        "shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_forward.dat",
        "shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_recon_fit.dat"
    );
    
    if (!loaded) {
        std::cerr << "Failed to load matrix files!" << std::endl;
        std::cerr << "Make sure matrix files are in the current directory." << std::endl;
        return 1;
    }
    
    std::cout << "Matrix files loaded successfully!" << std::endl << std::endl;
    
    // Test 1: Central trajectory (should pass)
    std::cout << "Test 1: Central trajectory" << std::endl;
    {
        SHMS::TrackState track;
        track.x = 0.0;      // cm
        track.y = 0.0;      // cm
        track.z = 0.0;      // cm
        track.dx = 0.0;     // slope
        track.dy = 0.0;     // slope
        track.delta = 0.0;  // %
        track.p = 3.0;      // GeV/c
        track.m2 = 0.139*0.139;  // Pion mass squared
        
        bool accepted = shms.Transport(track);
        
        std::cout << "  Accepted: " << (accepted ? "YES" : "NO") << std::endl;
        std::cout << "  Final position: (" << track.x << ", " << track.y << ", " << track.z << ") cm" << std::endl;
        std::cout << "  Path length: " << track.pathlen << " cm" << std::endl << std::endl;
    }
    
    // Test 2: Off-axis trajectory (moderate angles)
    std::cout << "Test 2: Off-axis trajectory (+10 mrad X, +5 mrad Y)" << std::endl;
    {
        SHMS::TrackState track;
        track.x = 0.5;           // cm
        track.y = 0.3;           // cm
        track.z = 0.0;           // cm
        track.dx = 0.010;        // 10 mrad
        track.dy = 0.005;        // 5 mrad
        track.delta = 2.0;       // +2%
        track.p = 3.0;           // GeV/c
        track.m2 = 0.139*0.139;
        
        bool accepted = shms.Transport(track);
        
        std::cout << "  Accepted: " << (accepted ? "YES" : "NO") << std::endl;
        std::cout << "  Final position: (" << track.x << ", " << track.y << ", " << track.z << ") cm" << std::endl;
        std::cout << "  Path length: " << track.pathlen << " cm" << std::endl << std::endl;
    }
    
    // Test 3: Large angle (likely to be rejected)
    std::cout << "Test 3: Large angle trajectory (+50 mrad X)" << std::endl;
    {
        SHMS::TrackState track;
        track.x = 2.0;           // cm
        track.y = 1.0;           // cm
        track.z = 0.0;           // cm
        track.dx = 0.050;        // 50 mrad
        track.dy = 0.020;        // 20 mrad
        track.delta = 5.0;       // +5%
        track.p = 3.0;           // GeV/c
        track.m2 = 0.139*0.139;
        
        bool accepted = shms.Transport(track);
        
        std::cout << "  Accepted: " << (accepted ? "YES" : "NO") << std::endl;
        if (accepted) {
            std::cout << "  Final position: (" << track.x << ", " << track.y << ", " << track.z << ") cm" << std::endl;
            std::cout << "  Path length: " << track.pathlen << " cm" << std::endl;
        }
        std::cout << std::endl;
    }
    
    // Test 4: Monte Carlo acceptance test
    std::cout << "Test 4: Monte Carlo acceptance test (1000 events)" << std::endl;
    {
        shms.ResetStatistics();
        
        const int n_events = 1000;
        const double x_rms = 0.5;     // cm
        const double y_rms = 0.5;     // cm
        const double dx_rms = 0.01;   // 10 mrad
        const double dy_rms = 0.01;   // 10 mrad
        const double delta_rms = 3.0; // %
        
        // Simple uniform random generator
        auto uniform = []() { return (double)rand() / RAND_MAX; };
        auto gaussian = [&uniform]() {
            // Box-Muller transform
            double u1 = uniform();
            double u2 = uniform();
            return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
        };
        
        srand(12345);  // Fixed seed for reproducibility
        
        for (int i = 0; i < n_events; ++i) {
            SHMS::TrackState track;
            track.x = gaussian() * x_rms;
            track.y = gaussian() * y_rms;
            track.z = 0.0;
            track.dx = gaussian() * dx_rms;
            track.dy = gaussian() * dy_rms;
            track.delta = gaussian() * delta_rms;
            track.p = 3.0;
            track.m2 = 0.139*0.139;
            
            shms.Transport(track);
        }
        
        const auto& stats = shms.GetStatistics();
        
        std::cout << "  Total events: " << stats.total_events << std::endl;
        std::cout << "  Accepted: " << stats.accepted << std::endl;
        std::cout << "  Acceptance: " << (stats.GetAcceptance() * 100.0) << "%" << std::endl;
        std::cout << std::endl;
        
        std::cout << "  Rejection breakdown:" << std::endl;
        std::cout << "    HB: in=" << stats.stop_hb_in 
                  << " men=" << stats.stop_hb_men 
                  << " mex=" << stats.stop_hb_mex 
                  << " out=" << stats.stop_hb_out << std::endl;
        
        std::cout << "    Q1: in=" << stats.stop_q1_in 
                  << " men=" << stats.stop_q1_men 
                  << " mid=" << stats.stop_q1_mid 
                  << " mex=" << stats.stop_q1_mex 
                  << " out=" << stats.stop_q1_out << std::endl;
        
        std::cout << "    Q2: in=" << stats.stop_q2_in 
                  << " men=" << stats.stop_q2_men 
                  << " mid=" << stats.stop_q2_mid 
                  << " mex=" << stats.stop_q2_mex 
                  << " out=" << stats.stop_q2_out << std::endl;
        
        std::cout << "    Q3: in=" << stats.stop_q3_in 
                  << " men=" << stats.stop_q3_men 
                  << " mid=" << stats.stop_q3_mid 
                  << " mex=" << stats.stop_q3_mex 
                  << " out=" << stats.stop_q3_out << std::endl;
        
        std::cout << "    D1: in=" << stats.stop_d1_in 
                  << " flr=" << stats.stop_d1_flare 
                  << " men=" << stats.stop_d1_men << std::endl;
        std::cout << "        mid1-7=" << stats.stop_d1_mid1 
                  << "," << stats.stop_d1_mid2 
                  << "," << stats.stop_d1_mid3 
                  << "," << stats.stop_d1_mid4 
                  << "," << stats.stop_d1_mid5 
                  << "," << stats.stop_d1_mid6 
                  << "," << stats.stop_d1_mid7 << std::endl;
        std::cout << "        mex=" << stats.stop_d1_mex 
                  << " out=" << stats.stop_d1_out << std::endl;
    }
    
    std::cout << std::endl << "=== SHMS Transport Test Complete ===" << std::endl;
    
    return 0;
}
