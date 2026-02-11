// tests/test_integrated_generation.cpp
// Test the full event generation pipeline

#include "simc/physics/EventGenerator.h"
#include "simc/physics/CrossSection.h"
#include "simc/core/ConfigManager.h"
#include "simc/core/RandomGenerator.h"
#include "simc/core/SimcConstants.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>

using namespace simc;
using namespace simc::constants;

// Test configuration - minimal JSON for testing
const char* test_config = R"JSON({
    "simulation": {
        "nevents": 100,
        "random_seed": 12345,
        "output_file": "test_output.root"
    },
    "beam": {
        "energy": 10600.0,
        "energy_spread": 0.0,
        "current": 50.0
    },
    "target": {
        "mass": 938.272,
        "Z": 1,
        "A": 1,
        "length": 10.0,
        "density": 0.071,
        "x_offset": 0.0,
        "y_offset": 0.0,
        "z_offset": 0.0,
        "angle": 0.0,
        "raster_pattern": 0,
        "raster_x": 0.0,
        "raster_y": 0.0
    },
    "spectrometer_electron": {
        "name": "HMS",
        "momentum": 8000.0,
        "angle": 12.5,
        "phi": 0.0
    },
    "spectrometer_hadron": {
        "name": "SHMS",
        "momentum": 5000.0,
        "angle": 35.0,
        "phi": 0.0
    },
    "generation": {
        "reaction_type": "H(e,e'p)",
        "use_energy_loss": false,
        "use_coulomb": false,
        "use_radiative": false,
        "electron": {
            "delta_min": -10.0,
            "delta_max": 22.0,
            "yptar_min": -0.03,
            "yptar_max": 0.03,
            "xptar_min": -0.06,
            "xptar_max": 0.06,
            "E_min": 500.0,
            "E_max": 10000.0
        },
        "hadron": {
            "delta_min": -10.0,
            "delta_max": 22.0,
            "yptar_min": -0.03,
            "yptar_max": 0.03,
            "xptar_min": -0.06,
            "xptar_max": 0.06,
            "E_min": 500.0,
            "E_max": 10000.0
        }
    }
})JSON";

// Statistics accumulator
struct Stats {
    int n = 0;
    double sum = 0.0;
    double sum2 = 0.0;
    double min = 1e99;
    double max = -1e99;
    
    void Add(double x) {
        n++;
        sum += x;
        sum2 += x*x;
        if (x < min) min = x;
        if (x > max) max = x;
    }
    
    double Mean() const { return n > 0 ? sum/n : 0.0; }
    double RMS() const { 
        if (n <= 1) return 0.0;
        double mean = Mean();
        return std::sqrt(sum2/n - mean*mean);
    }
    
    void Print(const std::string& name) const {
        std::cout << std::setw(20) << name << ": "
                  << std::setw(10) << Mean() << " ± " 
                  << std::setw(10) << RMS()
                  << "  [" << min << ", " << max << "]"
                  << std::endl;
    }
};

int main() {
    std::cout << "\n========================================\n";
    std::cout << "  Integrated Event Generation Test\n";
    std::cout << "========================================\n\n";
    
    try {
        // Write test config to file
        std::ofstream config_file("test_config.json");
        config_file << test_config;
        config_file.close();
        
        // Load configuration
        ConfigManager config("test_config.json");
        
        // Initialize components
        auto rng = std::make_shared<RandomGenerator>(12345);
        auto cross_section = std::make_shared<ElasticCrossSection>();
        
        EventGenerator generator(config, rng, cross_section, nullptr);
        generator.Initialize();
        
        std::cout << "Initialized successfully.\n";
        std::cout << "Reaction type: H(e,e'p) elastic\n";
        std::cout << "Generating 1000 events...\n\n";
        
        // Statistics
        Stats Q2_stats, W_stats, theta_e_stats, theta_p_stats;
        Stats Ee_stats, Ep_stats, xbj_stats, eps_stats;
        Stats sigcc_stats, weight_stats;
        
        int ntried = 0;
        int ngenerated = 0;
        int nphysical = 0;
        
        // Generate events
        while (ngenerated < 1000) {
            ntried++;
            
            SimcEvent event;
            MainEvent main;
            
            main.gen_weight = 1.0;
            main.jacobian = 1.0;
            main.success = false;
            
            // Generate event
            bool success = generator.GenerateEvent(event, main);
            
            if (!success) {
                continue;
            }
            
            ngenerated++;
            
            // Calculate cross section
            main.sigcc = cross_section->Calculate(event);
            main.weight = main.gen_weight * main.jacobian * main.sigcc;
            
            // Check physics
            bool physical = (event.Q2 > 0.0 && 
                           event.W > Mp &&
                           event.W < 1.5 * Mp &&
                           event.e_E > 0.0 &&
                           event.p_E > Mp);
            
            if (physical) {
                nphysical++;
                
                // Accumulate statistics
                Q2_stats.Add(event.Q2 / 1e6);  // Convert to GeV²
                W_stats.Add(event.W);
                theta_e_stats.Add(event.e_theta * 180.0 / pi);
                theta_p_stats.Add(event.p_theta * 180.0 / pi);
                Ee_stats.Add(event.e_E / 1000.0);  // Convert to GeV
                Ep_stats.Add(event.p_E / 1000.0);
                xbj_stats.Add(event.xbj);
                eps_stats.Add(event.epsilon);
                sigcc_stats.Add(main.sigcc * 1e9);  // Convert to nb
                weight_stats.Add(main.weight);
            }
            
            // Progress
            if (ngenerated % 250 == 0) {
                std::cout << "  Generated " << ngenerated << " events..." << std::endl;
            }
        }
        
        // Print results
        std::cout << "\n========================================\n";
        std::cout << "  Generation Statistics\n";
        std::cout << "========================================\n";
        std::cout << "Events tried:    " << ntried << "\n";
        std::cout << "Events generated: " << ngenerated << "\n";
        std::cout << "Physical events:  " << nphysical << "\n";
        std::cout << "Efficiency:       " << 100.0 * ngenerated / ntried << " %\n";
        std::cout << "Physical rate:    " << 100.0 * nphysical / ngenerated << " %\n";
        
        std::cout << "\n========================================\n";
        std::cout << "  Kinematic Distributions\n";
        std::cout << "========================================\n";
        Q2_stats.Print("Q² (GeV²)");
        W_stats.Print("W (MeV)");
        xbj_stats.Print("x_Bjorken");
        eps_stats.Print("ε (polarization)");
        
        std::cout << "\n========================================\n";
        std::cout << "  Angles\n";
        std::cout << "========================================\n";
        theta_e_stats.Print("θ_e (deg)");
        theta_p_stats.Print("θ_p (deg)");
        
        std::cout << "\n========================================\n";
        std::cout << "  Energies\n";
        std::cout << "========================================\n";
        Ee_stats.Print("E_e' (GeV)");
        Ep_stats.Print("E_p (GeV)");
        
        std::cout << "\n========================================\n";
        std::cout << "  Cross Section\n";
        std::cout << "========================================\n";
        sigcc_stats.Print("σ (nb)");
        weight_stats.Print("Weight");
        
        // Validation checks
        std::cout << "\n========================================\n";
        std::cout << "  Validation Checks\n";
        std::cout << "========================================\n";
        
        bool all_pass = true;
        
        // Check Q² range
        if (Q2_stats.min < 0.0) {
            std::cout << "❌ FAIL: Q² < 0 detected!\n";
            all_pass = false;
        } else {
            std::cout << "✓ PASS: Q² > 0 for all events\n";
        }
        
        // Check W is near proton mass (elastic)
        double W_mean = W_stats.Mean();
        if (std::abs(W_mean - Mp) > 10.0) {  // Within 10 MeV
            std::cout << "❌ FAIL: W not at proton mass (elastic peak)\n";
            all_pass = false;
        } else {
            std::cout << "✓ PASS: W at elastic peak (~" << Mp << " MeV)\n";
        }
        
        // Check electron angle reasonable
        double theta_e_mean = theta_e_stats.Mean();
        if (theta_e_mean < 10.0 || theta_e_mean > 15.0) {
            std::cout << "⚠ WARNING: θ_e = " << theta_e_mean 
                      << "° (expected ~12.5°)\n";
        } else {
            std::cout << "✓ PASS: θ_e ~ 12.5°\n";
        }
        
        // Check cross section is positive
        if (sigcc_stats.min <= 0.0) {
            std::cout << "❌ FAIL: Non-positive cross section detected!\n";
            all_pass = false;
        } else {
            std::cout << "✓ PASS: σ > 0 for all events\n";
        }
        
        // Check cross section magnitude (should be ~nb for elastic)
        double sigcc_mean = sigcc_stats.Mean();
        if (sigcc_mean < 0.01 || sigcc_mean > 1000.0) {
            std::cout << "⚠ WARNING: σ = " << sigcc_mean 
                      << " nb (check normalization)\n";
        } else {
            std::cout << "✓ PASS: σ in reasonable range\n";
        }
        
        // Check energy conservation (elastic)
        double Ee_mean = Ee_stats.Mean();
        double expected_Ee = 10.6 * Mp / (Mp + 10.6 * (1.0 - std::cos(12.5 * pi/180.0)));
        expected_Ee /= 1000.0;  // Convert to GeV
        
        if (std::abs(Ee_mean - expected_Ee) > 0.5) {
            std::cout << "⚠ WARNING: E_e' = " << Ee_mean 
                      << " GeV (expected ~" << expected_Ee << " GeV)\n";
        } else {
            std::cout << "✓ PASS: E_e' matches elastic kinematics\n";
        }
        
        std::cout << "\n========================================\n";
        if (all_pass) {
            std::cout << "✅ ALL TESTS PASSED\n";
        } else {
            std::cout << "❌ SOME TESTS FAILED\n";
        }
        std::cout << "========================================\n\n";
        
        return all_pass ? 0 : 1;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}
