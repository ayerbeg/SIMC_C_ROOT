// src/main.cpp
// Main program for SIMC Monte Carlo

#include "simc/SimcEvent.h"
#include "simc/ConfigManager.h"
#include "simc/RandomGenerator.h"
#include "simc/OutputManager.h"
#include "simc/SimcConstants.h"

#include <iostream>
#include <cmath>

using namespace simc;
using namespace simc::constants;

/**
 * @brief Simple test event generator
 * 
 * This is a placeholder that generates dummy events for testing.
 * In Phase 3, this will be replaced with the full physics event generator.
 */
void GenerateTestEvent(SimcEvent& evt, MainEvent& main, RandomGenerator& rng) {
    // Generate random kinematic quantities
    evt.Ein = 10600.0;  // 10.6 GeV
    
    // Electron arm
    evt.e_delta = rng.Gaussian(0.0, 2.0);  // ~2% resolution
    evt.e_xptar = rng.Gaussian(0.0, 0.002); // ~2 mrad
    evt.e_yptar = rng.Gaussian(0.0, 0.002);
    evt.e_theta = 12.5 * DEG_TO_RAD + evt.e_yptar;
    evt.e_phi = 0.0 + evt.e_xptar;
    evt.e_P = 8000.0 * (1.0 + evt.e_delta/100.0);
    evt.e_E = evt.e_P;  // Electron is essentially massless
    
    // Hadron arm
    evt.p_delta = rng.Gaussian(0.0, 2.0);
    evt.p_xptar = rng.Gaussian(0.0, 0.003);
    evt.p_yptar = rng.Gaussian(0.0, 0.003);
    evt.p_theta = 35.0 * DEG_TO_RAD + evt.p_yptar;
    evt.p_phi = 0.0 + evt.p_xptar;
    evt.p_P = 5000.0 * (1.0 + evt.p_delta/100.0);
    evt.p_E = std::sqrt(evt.p_P*evt.p_P + Mp2);
    
    // Calculate kinematic quantities
    evt.nu = evt.Ein - evt.e_E;
    evt.Q2 = 2.0 * evt.Ein * evt.e_E * (1.0 - std::cos(evt.e_theta));
    evt.q = std::sqrt(evt.Q2 + evt.nu*evt.nu);
    evt.W = std::sqrt(Mp2 + 2.0*Mp*evt.nu - evt.Q2);
    evt.xbj = evt.Q2 / (2.0 * Mp * evt.nu);
    evt.epsilon = 1.0 / (1.0 + 2.0*(1.0 + evt.nu*evt.nu/evt.Q2)*
                         std::tan(evt.e_theta/2.0)*std::tan(evt.e_theta/2.0));
    
    // Missing quantities (simple calculation)
    evt.Em = evt.nu - evt.p_E + Mp;
    evt.Pm = std::abs(rng.Gaussian(100.0, 50.0));  // Dummy value
    
    // Set weights (dummy values for testing)
    main.weight = 1.0;
    main.jacobian = 1.0;
    main.sigcc = 1.0e-3;  // 1 nb
    main.success = (std::abs(evt.e_delta) < 10.0 && 
                    std::abs(evt.p_delta) < 15.0);
}

int main(int argc, char** argv) {
    std::cout << "\n";
    std::cout << "==========================================\n";
    std::cout << "   SIMC C++/ROOT Monte Carlo - Phase 2   \n";
    std::cout << "==========================================\n";
    std::cout << "\n";
    
    // Default configuration file
    std::string config_file = "data/config/default.json";
    if (argc > 1) {
        config_file = argv[1];
    }
    
    try {
        // Load configuration
        std::cout << "Loading configuration from: " << config_file << std::endl;
        ConfigManager config(config_file);
        config.Print();
        
        // Get parameters
        int nevents = config.GetNEvents();
        unsigned int seed = config.GetRandomSeed();
        std::string output_file = config.GetOutputFile();
        
        std::cout << "\nSimulation parameters:" << std::endl;
        std::cout << "  Number of events: " << nevents << std::endl;
        std::cout << "  Random seed:      " << seed << std::endl;
        std::cout << "  Output file:      " << output_file << std::endl;
        std::cout << "\n";
        
        // Initialize random generator
        RandomGenerator rng(seed);
        
        // Initialize output
        OutputManager output(output_file);
        output.Initialize();
        
        // Event generation loop
        std::cout << "Generating events..." << std::endl;
        
        int ngenerated = 0;
        int ncontribute = 0;
        
        for (int i = 0; i < nevents; ++i) {
            // Progress indicator
            if ((i+1) % 1000 == 0) {
                std::cout << "  Generated " << (i+1) << " events..." << std::endl;
            }
            
            // Generate event
            SimcEvent evt;
            MainEvent main;
            GenerateTestEvent(evt, main, rng);
            
            ++ngenerated;
            
            // Fill generated tree
            output.FillGenerated(evt);
            
            // Check if event passes cuts
            if (main.success) {
                ++ncontribute;
                output.FillEvent(evt, main);
                output.FillHistograms(evt, main.weight, false);
            }
            
            output.IncrementTried();
        }
        
        // Set final statistics
        output.SetStatistics(nevents, ngenerated, ncontribute);
        
        // Finalize and write output
        std::cout << "\nFinalizing output..." << std::endl;
        output.Finalize();
        
        // Print summary
        double efficiency = 100.0 * ncontribute / static_cast<double>(nevents);
        std::cout << "\n==========================================\n";
        std::cout << "Simulation complete!" << std::endl;
        std::cout << "  Events generated:   " << ngenerated << std::endl;
        std::cout << "  Events passed cuts: " << ncontribute << std::endl;
        std::cout << "  Efficiency:         " << efficiency << " %" << std::endl;
        std::cout << "==========================================\n";
        std::cout << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
