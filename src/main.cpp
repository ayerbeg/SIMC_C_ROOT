// src/main.cpp - Phase 5c.2 IN PROGRESS
// Integrating SHMS spectrometer transport

#include "simc/core/SimcEvent.h"
#include "simc/core/ConfigManager.h"
#include "simc/core/RandomGenerator.h"
#include "simc/io/OutputManager.h"
#include "simc/physics/EventGenerator.h"
#include "simc/transport/MonteCarloTransport.h"
#include "simc/spectrometers/SHMS.h"
#include "simc/physics/CrossSection.h"
#include "simc/core/SimcConstants.h"

#include <iostream>
#include <memory>
#include <cmath>

using namespace simc;
using namespace simc::constants;

int main(int argc, char** argv) {
    std::cout << "\n";
    std::cout << "==========================================\n";
    std::cout << "   SIMC C++/ROOT Monte Carlo - Phase 5c.2\n";
    std::cout << "   SHMS Spectrometer Transport            \n";
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
        
        // ====================================================================
        // PHASE 5c.2: Initialize SHMS Spectrometer
        // ====================================================================
        std::cout << "Initializing SHMS spectrometer..." << std::endl;
        
        SHMS shms;
        
        // Load matrix files (using Fortran-style symlink names)
        std::cout << "  Loading matrix files..." << std::endl;
	if (!shms.LoadMatrices("../data/matrices/shms/shms_forward.dat", "../data/matrices/shms/shms_recon.dat")) {
            throw std::runtime_error("Failed to load SHMS matrix files. "
                                   "Run: cd data/shms && ls -l shms_*.dat");
        }
        std::cout << "  Matrix files loaded successfully" << std::endl;
        
        // Get initial statistics (all zeros)
        auto initial_stats = shms.GetStats();
        std::cout << "SHMS initialized." << std::endl;
        std::cout << "\n";
        
        // Spectrometer optics (Phase 5c.1 style) - kept for compatibility
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
        
        // Initialize MC transport (Phase 5c.1 - core transport)
        MonteCarloTransport transport(config, rng);
        
        generator.ResetStatistics();
        
        while (ngenerated < nevents) {
            ++ntried;
            
            // Progress indicator
            if (ntried % 1000 == 1) {
                std::cout << "  Event " << ntried 
                          << " (" << ngenerated << " generated, "
                          << ncontribute << " accepted)" << std::endl;
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
            // PHASE 5c.1: Core Monte Carlo Transport
            // - Apply multiple scattering (beam, electron, hadron)
            // - Apply energy loss corrections
            // - Calculate spectrometer angles
            // ================================================================
            bool passes_core_transport = transport.Transport(event, main);
            
            if (!passes_core_transport) {
                // Event failed core transport (beam/target effects)
                continue;
            }
            
            // ================================================================
            // PHASE 5c.2: SHMS Spectrometer Transport
            // - Full aperture checking (32 apertures)
            // - Matrix element transport
            // - Focal plane calculations
            // ================================================================
            
            // Create track state from spectrometer coordinates
            SHMS::TrackState track;
            track.x = 0.0;  // Start at target center
            track.y = 0.0;
            track.dx = main.SP_hadron.xptar;  // Spectrometer angles (rad)
            track.dy = main.SP_hadron.yptar;
            track.delta = main.SP_hadron.delta;  // Momentum percent
            track.z = 0.0;
            track.pathlen = 0.0;
            
            // Transport through SHMS
            bool accepted_by_shms = shms.Transport(track);
            
            if (accepted_by_shms) {
                // ============================================================
                // PHASE 5c.2: Store focal plane coordinates
                // ============================================================
                // TODO: Add focal plane variables to SimcEvent
                // For now, just count as accepted
                
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
        
        // Get final SHMS statistics
        auto final_stats = shms.GetStats();
        
        // Print summary
        double generation_eff = 100.0 * ngenerated / static_cast<double>(ntried);
        double acceptance = 100.0 * ncontribute / static_cast<double>(ngenerated);
        
        std::cout << "\n==========================================\n";
        std::cout << "Phase 5c.2 Session 1 Complete!\n";
        std::cout << "--------------------------------------------\n";
        std::cout << "Event Statistics:\n";
        std::cout << "  Events tried:       " << ntried << "\n";
        std::cout << "  Events generated:   " << ngenerated 
                  << " (" << generation_eff << "%)\n";
        std::cout << "  Events accepted:    " << ncontribute 
                  << " (" << acceptance << "%)\n";
        std::cout << "--------------------------------------------\n";
        std::cout << "SHMS Transport Statistics:\n";
        std::cout << "  Total transported:  " << final_stats.total_events << "\n";
        std::cout << "  Accepted:           " << final_stats.accepted << "\n";
        
        // Show rejection breakdown
        int total_rejected = final_stats.total_events - final_stats.accepted;
        if (total_rejected > 0) {
            std::cout << "  Rejected:           " << total_rejected << "\n";
            std::cout << "\n  Rejection Breakdown:\n";
            
            if (final_stats.shmsSTOP_HB > 0)
                std::cout << "    Holding Box:      " << final_stats.shmsSTOP_HB << "\n";
            if (final_stats.shmsSTOP_Q1 > 0)
                std::cout << "    Q1 apertures:     " << final_stats.shmsSTOP_Q1 << "\n";
            if (final_stats.shmsSTOP_Q2 > 0)
                std::cout << "    Q2 apertures:     " << final_stats.shmsSTOP_Q2 << "\n";
            if (final_stats.shmsSTOP_Q3 > 0)
                std::cout << "    Q3 apertures:     " << final_stats.shmsSTOP_Q3 << "\n";
            if (final_stats.shmsSTOP_D1 > 0)
                std::cout << "    Dipole:           " << final_stats.shmsSTOP_D1 << "\n";
            if (final_stats.shmsSTOP_COLL > 0)
                std::cout << "    Collimator:       " << final_stats.shmsSTOP_COLL << "\n";
            if (final_stats.shmsSTOP_HUT > 0)
                std::cout << "    Hut/Detectors:    " << final_stats.shmsSTOP_HUT << "\n";
        }
        
        std::cout << "==========================================\n";
        std::cout << "\n";
        
        std::cout << "Output written to: " << output_filename << "\n";
        std::cout << "\nPhase 5c.2 Session 1 Implemented:\n";
        std::cout << "  ✓ SHMS class with 32 aperture checks\n";
        std::cout << "  ✓ Matrix element transport (forward)\n";
        std::cout << "  ✓ Holding Box (4 tilted apertures)\n";
        std::cout << "  ✓ Quadrupoles Q1, Q2, Q3 (15 apertures)\n";
        std::cout << "  ✓ Dipole D1 (12 tilted apertures)\n";
        std::cout << "  ✓ Hut transport to focal plane\n";
        std::cout << "  ✓ Detailed rejection statistics\n";
        std::cout << "\nPhase 5c.1 Features (Still Active):\n";
        std::cout << "  ✓ Multiple scattering (beam, e, p)\n";
        std::cout << "  ✓ Energy loss corrections\n";
        std::cout << "  ✓ Hadron spectrometer angles (elastic)\n";
        std::cout << "\nNext Steps:\n";
        std::cout << "  - Session 2: Implement Reconstruct() (inverse transport)\n";
        std::cout << "  - Session 3: Add focal plane variables to SimcEvent\n";
        std::cout << "  - Session 4: Validate against Fortran SIMC output\n";
        std::cout << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Simulation terminated.\n" << std::endl;
        return 1;
    }
    
    return 0;
}
