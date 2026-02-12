// src/main.cpp - Phase 5c.2 DEBUG VERSION
// Shows exactly what's failing

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
#include <iomanip>
#include <memory>
#include <cmath>

using namespace simc;
using namespace simc::constants;

int main(int argc, char** argv) {
    std::cout << "\n";
    std::cout << "==========================================\n";
    std::cout << "   SIMC C++/ROOT Monte Carlo - Phase 5c.2\n";
    std::cout << "   SHMS Transport - DEBUG VERSION    \n";
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
        
        // Initialize cross section
        auto cross_section = std::make_shared<ElasticCrossSection>();
        
        // ====================================================================
        // PHASE 5c.2: Initialize SHMS Spectrometer
        // ====================================================================
        std::cout << "Initializing SHMS spectrometer..." << std::endl;
        
        SHMS shms;
        
        std::cout << "  Loading matrix files..." << std::endl;
        if (!shms.LoadMatrices("data/matrices/shms/shms_forward.dat", 
                                "data/matrices/shms/shms_recon.dat")) {
            throw std::runtime_error("Failed to load SHMS matrix files.");
        }
        std::cout << "  Matrix files loaded successfully" << std::endl;
        std::cout << "SHMS initialized." << std::endl;
        std::cout << "\n";
        
        // Spectrometer optics
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
        
        MonteCarloTransport transport(config, rng);
        generator.ResetStatistics();
        
        while (ngenerated < nevents) {
            ++ntried;
            
            if (ntried % 1000 == 1) {
                std::cout << "  Event " << ntried 
                          << " (" << ngenerated << " generated, "
                          << ncontribute << " accepted)" << std::endl;
            }
            
            SimcEvent event;
            MainEvent main;
            
            main.gen_weight = 1.0;
            main.jacobian = 1.0;
            main.weight = 1.0;
            main.sigcc = 1.0;
            main.success = false;
            
            // ================================================================
            // Generate physics event
            // ================================================================
            bool success = generator.GenerateEvent(event, main);
            if (!success) {
                continue;
            }
            
            main.sigcc = cross_section->Calculate(event);
            main.weight = main.gen_weight * main.jacobian * main.sigcc;
            
            ++ngenerated;
            main.success = true;
            
            output.FillGenerated(event);
            output.FillHistograms(event, main.weight, true);
            
            // ================================================================
            // KEY FIX: Calculate hadron angles from physics BEFORE transport
            // ================================================================
            transport.CalculateReconstructed(event);
            
            // ================================================================
            // Core Monte Carlo Transport
            // ================================================================
            bool passes_core_transport = transport.Transport(event, main);
            
            // ================================================================
            // DEBUG: Print details for first 10 events BEFORE checking pass/fail
            // ================================================================
            if (ngenerated <= 10) {
                std::cout << "\nEvent " << ngenerated << " details:" << std::endl;
                std::cout << "  Physics angles:" << std::endl;
                std::cout << "    p_theta = " << (event.p_theta * 180.0 / M_PI) << " deg" << std::endl;
                std::cout << "    p_phi   = " << (event.p_phi * 180.0 / M_PI) << " deg" << std::endl;
                std::cout << "  Event xptar/yptar (from physics):" << std::endl;
                std::cout << "    p_xptar = " << (event.p_xptar * 1000.0) << " mrad" << std::endl;
                std::cout << "    p_yptar = " << (event.p_yptar * 1000.0) << " mrad" << std::endl;
                std::cout << "  After transport (with MS):" << std::endl;
                std::cout << "    SP_hadron.xptar = " << (main.SP_hadron.xptar * 1000.0) << " mrad" << std::endl;
                std::cout << "    SP_hadron.yptar = " << (main.SP_hadron.yptar * 1000.0) << " mrad" << std::endl;
                std::cout << "    SP_hadron.delta = " << main.SP_hadron.delta << " %" << std::endl;
                std::cout << "    SP_electron.xptar = " << (main.SP_electron.xptar * 1000.0) << " mrad" << std::endl;
                std::cout << "    SP_electron.yptar = " << (main.SP_electron.yptar * 1000.0) << " mrad" << std::endl;
                std::cout << "    SP_electron.delta = " << main.SP_electron.delta << " %" << std::endl;
                std::cout << "  Core transport: " << (passes_core_transport ? "PASS" : "FAIL") << std::endl;
                
                if (!passes_core_transport) {
                    std::cout << "  REJECTION REASONS:" << std::endl;
                    if (main.SP_hadron.delta < -10.0 || main.SP_hadron.delta > 22.0) {
                        std::cout << "    - Hadron delta out of [-10, 22]%" << std::endl;
                    }
                    if (main.SP_hadron.xptar < -0.06 || main.SP_hadron.xptar > 0.06) {
                        std::cout << "    - Hadron xptar out of ±60 mrad" << std::endl;
                    }
                    if (main.SP_hadron.yptar < -0.03 || main.SP_hadron.yptar > 0.03) {
                        std::cout << "    - Hadron yptar out of ±30 mrad" << std::endl;
                    }
                    if (main.SP_electron.delta < -10.0 || main.SP_electron.delta > 22.0) {
                        std::cout << "    - Electron delta out of [-10, 22]%" << std::endl;
                    }
                    if (main.SP_electron.xptar < -0.06 || main.SP_electron.xptar > 0.06) {
                        std::cout << "    - Electron xptar out of ±60 mrad" << std::endl;
                    }
                    if (main.SP_electron.yptar < -0.03 || main.SP_electron.yptar > 0.03) {
                        std::cout << "    - Electron yptar out of ±30 mrad" << std::endl;
                    }
                }
            }
            
            if (!passes_core_transport) {
                continue;
            }
            
            // ================================================================
            // SHMS Spectrometer Transport
            // ================================================================
            
            SHMS::TrackState track;
            track.x = 0.0;
            track.y = 0.0;
            track.dx = main.SP_hadron.xptar;
            track.dy = main.SP_hadron.yptar;
            track.delta = main.SP_hadron.delta;
            track.z = 0.0;
            track.pathlen = 0.0;
            track.p = event.p_P / 1000.0;  // MeV/c → GeV/c
            track.m2 = 0.938272 * 0.938272;  // Proton mass^2
            
            if (ngenerated <= 10) {
                std::cout << "  SHMS input: dx=" << (track.dx*1000.0) << " mrad, "
                          << "dy=" << (track.dy*1000.0) << " mrad, "
                          << "delta=" << track.delta << "%, "
                          << "p=" << track.p << " GeV" << std::endl;
            }
            
            bool accepted_by_shms = shms.Transport(track);
            
            if (ngenerated <= 10) {
                std::cout << "  SHMS result: " << (accepted_by_shms ? "ACCEPT" : "REJECT") << std::endl;
            }
            
            if (accepted_by_shms) {
                event.e_delta = main.SP_electron.delta;
                event.e_xptar = main.SP_electron.xptar;
                event.e_yptar = main.SP_electron.yptar;
                event.p_delta = main.SP_hadron.delta;
                
                ++ncontribute;
                output.FillEvent(event, main);
                output.FillHistograms(event, main.weight, false);
            }
            
            output.IncrementTried();
        }
        
        // ====================================================================
        // FINALIZATION
        // ====================================================================
        
        output.SetStatistics(ntried, ngenerated, ncontribute);
        
        std::cout << "\nFinalizing output..." << std::endl;
        output.Finalize();
        
        auto final_stats = shms.GetStats();
        
        double generation_eff = 100.0 * ngenerated / static_cast<double>(ntried);
        double acceptance = 100.0 * ncontribute / static_cast<double>(ngenerated);
        
        std::cout << "\n==========================================\n";
        std::cout << "Phase 5c.2 Debug Complete!\n";
        std::cout << "--------------------------------------------\n";
        std::cout << "Event Statistics:\n";
        std::cout << "  Events tried:       " << ntried << "\n";
        std::cout << "  Events generated:   " << ngenerated 
                  << " (" << std::fixed << std::setprecision(1) << generation_eff << "%)\n";
        std::cout << "  Events accepted:    " << ncontribute 
                  << " (" << acceptance << "%)\n";
        std::cout << "--------------------------------------------\n";
        std::cout << "SHMS Transport Statistics:\n";
        std::cout << "  Total transported:  " << final_stats.total_events << "\n";
        std::cout << "  Accepted:           " << final_stats.accepted << "\n";
        
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
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Simulation terminated.\n" << std::endl;
        return 1;
    }
    
    return 0;
}
