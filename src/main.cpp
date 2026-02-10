// src/main.cpp - Phase 5c.1 COMPLETE
// Integrates MC transport: multiple scattering, energy loss, acceptance

#include "simc/SimcEvent.h"
#include "simc/ConfigManager.h"
#include "simc/RandomGenerator.h"
#include "simc/OutputManager.h"
#include "simc/EventGenerator.h"
#include "simc/MonteCarloTransport.h"
#include "simc/CrossSection.h"
#include "simc/SimcConstants.h"

#include <iostream>
#include <memory>
#include <cmath>

using namespace simc;
using namespace simc::constants;

int main(int argc, char** argv) {
    std::cout << "\n";
    std::cout << "==========================================\n";
    std::cout << "   SIMC C++/ROOT Monte Carlo - Phase 5c.1\n";
    std::cout << "   MC Transport (Core) Complete           \n";
    std::cout << "==========================================\n";
    std::cout << "\n";
    
    std::string config_file = "data/config/default.json";
    if (argc > 1) {
        config_file = argv[1];
    }
    
    try {
        // ====================================================================
        // INITIALIZATION
        // ====================================================================
        
        std::cout << "Loading configuration: " << config_file << std::endl;
        ConfigManager config(config_file);
        
        int nevents = config.GetNEvents();
        unsigned int seed = config.GetRandomSeed();
        std::string output_filename = config.GetOutputFile();
        
        std::cout << "  Events: " << nevents << std::endl;
        std::cout << "  Seed:   " << seed << std::endl;
        std::cout << "  Output: " << output_filename << std::endl;
        std::cout << "\n";
        
        // Initialize random generator
        auto rng = std::make_shared<RandomGenerator>(seed);
        
        // Initialize cross section - use ElasticCrossSection (concrete class)
        auto cross_section = std::make_shared<ElasticCrossSection>();
        
        // Spectrometer optics not needed until Phase 5c
        std::shared_ptr<SpectrometerOptics> optics = nullptr;
        
        // Initialize event generator
        std::cout << "Initializing event generator..." << std::endl;
        EventGenerator generator(config, rng, cross_section, optics);
        if (!generator.Initialize()) {
            throw std::runtime_error("Failed to initialize event generator");
        }
        std::cout << "Event generator initialized." << std::endl;
        
        // Initialize output
        OutputManager output(output_filename);
        output.Initialize();
        
        // ====================================================================
        // EVENT GENERATION LOOP
        // ====================================================================
        
        std::cout << "\nGenerating events..." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        
        int ntried = 0;
        int ngenerated = 0;
        int ncontribute = 0;
        
        // Initialize MC transport (Phase 5c.1)
        MonteCarloTransport transport(config, rng);
        
        generator.ResetStatistics();
        
        while (ngenerated < nevents) {
            ++ntried;
            
            // Progress indicator
            if (ntried % 1000 == 1) {
                std::cout << "  Event " << ntried 
                          << " (" << ngenerated << " generated)" << std::endl;
            }
            
            // Generate event
            SimcEvent event;
            MainEvent main;
            
            main.gen_weight = 1.0;
            main.jacobian = 1.0;
            main.weight = 1.0;
            main.sigcc = 1.0;
            main.success = false;
            
            // ================================================================
            // PHASE 5b: Generate physics event
            // ================================================================
            bool success = generator.GenerateEvent(event, main);
            if (!success) {
                continue;
            }
            
            // Calculate cross section
            main.sigcc = cross_section->Calculate(event);
            
            // Total weight
            main.weight = main.gen_weight * main.jacobian * main.sigcc;
            
            // Count as generated
            ++ngenerated;
            main.success = true;
            
            // Fill GENERATED quantities
            output.FillGenerated(event);
            output.FillHistograms(event, main.weight, true);  // generated=true
            
            // ================================================================
            // PHASE 5c.1: Monte Carlo Transport
            // - Apply multiple scattering (beam, electron, hadron)
            // - Apply energy loss corrections
            // - Check acceptance (simple delta/angle cuts)
            // ================================================================
            bool passes_acceptance = transport.Transport(event, main);
            
            if (passes_acceptance) {
                // ============================================================
                // PHASE 5c.1: Calculate Reconstructed Quantities
                // For H(e,e'p) elastic, calculate p_xptar, p_yptar
                // that were set to 0.0 in EventGenerator
                // ============================================================
                transport.CalculateReconstructed(event);
                
                // ============================================================
                // Copy spectrometer quantities back to event for storage
                // These include energy loss corrections and multiple scattering
                // ============================================================
                event.e_delta = main.SP_electron.delta;
                event.e_xptar = main.SP_electron.xptar;
                event.e_yptar = main.SP_electron.yptar;
                event.p_delta = main.SP_hadron.delta;
                // p_xptar and p_yptar already set by CalculateReconstructed()
                
                // Fill event tree and reconstructed histograms
                ++ncontribute;
                output.FillEvent(event, main);
                output.FillHistograms(event, main.weight, false);  // reconstructed
            }
            
            output.IncrementTried();
        }
        
        // ====================================================================
        // FINALIZATION
        // ====================================================================
        
        output.SetStatistics(ntried, ngenerated, ncontribute);
        
        std::cout << "\nFinalizing output..." << std::endl;
        output.Finalize();
        
        // Print summary
        double efficiency = 100.0 * ncontribute / static_cast<double>(ntried);
        double acceptance = 100.0 * ngenerated / static_cast<double>(ntried);
        
        std::cout << "\n==========================================\n";
        std::cout << "Phase 5c.1 Complete!\n";
        std::cout << "--------------------------------------------\n";
        std::cout << "  Events tried:       " << ntried << "\n";
        std::cout << "  Events generated:   " << ngenerated << "\n";
        std::cout << "  Events passed cuts: " << ncontribute << "\n";
        std::cout << "--------------------------------------------\n";
        std::cout << "  Generation eff:     " << acceptance << " %\n";
        std::cout << "  Acceptance eff:     " << efficiency << " %\n";
        std::cout << "==========================================\n";
        std::cout << "\n";
        
        std::cout << "Output written to: " << output_filename << "\n";
        std::cout << "\nPhase 5c.1 Features Implemented:\n";
        std::cout << "  ✓ Multiple scattering (beam, e, p)\n";
        std::cout << "  ✓ Energy loss corrections\n";
        std::cout << "  ✓ Hadron spectrometer angles (elastic)\n";
        std::cout << "  ✓ Simple acceptance (delta/angle cuts)\n";
        std::cout << "  ✓ Reconstructed histograms filled\n";
        std::cout << "\nNext: Phase 5c.2 - Full Spectrometer MCs\n";
        std::cout << "  - mc_hms, mc_sos, mc_shms transport\n";
        std::cout << "  - Focal plane calculations\n";
        std::cout << "  - Detailed aperture checking\n";
        std::cout << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Simulation terminated.\n" << std::endl;
        return 1;
    }
    
    return 0;
}
