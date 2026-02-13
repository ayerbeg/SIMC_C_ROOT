// src/main.cpp - Phase 5c.4 VERSION
// Two-arm coincidence: HMS (electron) + SHMS (hadron)

#include "simc/core/SimcEvent.h"
#include "simc/core/ConfigManager.h"
#include "simc/core/RandomGenerator.h"
#include "simc/io/OutputManager.h"
#include "simc/physics/EventGenerator.h"
#include "simc/transport/MonteCarloTransport.h"
#include "simc/spectrometers/HMS.h"    // PHASE 5c.4: Added HMS for electron arm
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
    std::cout << "   SIMC C++/ROOT Monte Carlo - Phase 5c.4\n";
    std::cout << "   Two-Arm Coincidence (HMS + SHMS)    \n";
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
        // PHASE 5c.4: Initialize BOTH Spectrometers (HMS + SHMS)
        // ====================================================================
        std::cout << "Initializing spectrometers..." << std::endl;
        
        // HMS (electron arm) - PHASE 5c.4: NEW!
        HMS hms;
        std::cout << "  Loading HMS matrix files..." << std::endl;
        if (!hms.LoadMatrices("data/matrices/hms/forward_cosy.dat", 
                               "data/matrices/hms/recon_cosy.dat")) {
            throw std::runtime_error("Failed to load HMS matrix files.");
        }
        std::cout << "  HMS initialized." << std::endl;
        
        // SHMS (hadron arm) - from Phase 5c.2
        SHMS shms;
        std::cout << "  Loading SHMS matrix files..." << std::endl;
        if (!shms.LoadMatrices("data/matrices/shms/shms_forward.dat", 
                                "data/matrices/shms/shms_recon.dat")) {
            throw std::runtime_error("Failed to load SHMS matrix files.");
        }
        std::cout << "  SHMS initialized." << std::endl;
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
            // Calculate hadron angles from physics BEFORE transport
            // ================================================================
            transport.CalculateReconstructed(event);
            
            // ================================================================
            // Core Monte Carlo Transport
            // ================================================================
            bool passes_core_transport = transport.Transport(event, main);
            
            // ================================================================
            // DEBUG: Print details for first 10 events
            // ================================================================
            if (ngenerated <= 10) {
                std::cout << "\nEvent " << ngenerated << " details:" << std::endl;
                std::cout << "  Electron (HMS):" << std::endl;
                std::cout << "    xptar = " << (main.SP_electron.xptar * 1000.0) << " mrad" << std::endl;
                std::cout << "    yptar = " << (main.SP_electron.yptar * 1000.0) << " mrad" << std::endl;
                std::cout << "    delta = " << main.SP_electron.delta << " %" << std::endl;
                std::cout << "  Hadron (SHMS):" << std::endl;
                std::cout << "    xptar = " << (main.SP_hadron.xptar * 1000.0) << " mrad" << std::endl;
                std::cout << "    yptar = " << (main.SP_hadron.yptar * 1000.0) << " mrad" << std::endl;
                std::cout << "    delta = " << main.SP_hadron.delta << " %" << std::endl;
                std::cout << "  Core transport: " << (passes_core_transport ? "PASS" : "FAIL") << std::endl;
            }
            
            if (!passes_core_transport) {
                continue;
            }
            
            // ================================================================
            // HMS Spectrometer Transport (Electron Arm) - PHASE 5c.4: NEW!
            // ================================================================
            
            HMS::TrackState e_track;
            e_track.x = 0.0;
            e_track.y = 0.0;
            e_track.dx = main.SP_electron.xptar;  // Electron angles (NOT hadron!)
            e_track.dy = main.SP_electron.yptar;
            e_track.delta = main.SP_electron.delta;
            e_track.z = 0.0;
            e_track.pathlen = 0.0;
            e_track.p = event.e_P / 1000.0;  // MeV/c → GeV/c (electron momentum)
            e_track.m2 = 0.000511 * 0.000511;  // Electron mass^2 (NOT proton!)
            
            if (ngenerated <= 10) {
                std::cout << "  HMS input: dx=" << (e_track.dx*1000.0) << " mrad, "
                          << "dy=" << (e_track.dy*1000.0) << " mrad, "
                          << "delta=" << e_track.delta << "%, "
                          << "p=" << e_track.p << " GeV" << std::endl;
            }
            
            bool accepted_by_hms = hms.Transport(e_track);
            
            if (ngenerated <= 10) {
                std::cout << "  HMS result: " << (accepted_by_hms ? "ACCEPT" : "REJECT") << std::endl;
            }
            
            // PHASE 5c.4: Skip event if electron rejected by HMS
            // This is an optimization - we don't need to check SHMS if HMS already rejected
            if (!accepted_by_hms) {
                continue;
            }
            
            // ================================================================
            // SHMS Spectrometer Transport (Hadron Arm) - from Phase 5c.2
            // ================================================================
            
            SHMS::TrackState h_track;
            h_track.x = 0.0;
            h_track.y = 0.0;
            h_track.dx = main.SP_hadron.xptar;  // Hadron angles
            h_track.dy = main.SP_hadron.yptar;
            h_track.delta = main.SP_hadron.delta;
            h_track.z = 0.0;
            h_track.pathlen = 0.0;
            h_track.p = event.p_P / 1000.0;  // MeV/c → GeV/c (hadron momentum)
            h_track.m2 = 0.938272 * 0.938272;  // Proton mass^2
            
            if (ngenerated <= 10) {
                std::cout << "  SHMS input: dx=" << (h_track.dx*1000.0) << " mrad, "
                          << "dy=" << (h_track.dy*1000.0) << " mrad, "
                          << "delta=" << h_track.delta << "%, "
                          << "p=" << h_track.p << " GeV" << std::endl;
            }
            
            bool accepted_by_shms = shms.Transport(h_track);
            
            if (ngenerated <= 10) {
                std::cout << "  SHMS result: " << (accepted_by_shms ? "ACCEPT" : "REJECT") << std::endl;
            }
            
            // ================================================================
            // PHASE 5c.4: Require BOTH spectrometers to accept (coincidence)
            // 
            // This matches Fortran SIMC logic:
            //   success = ok_E_arm .and. ok_P_arm
            // 
            // For two-arm coincidence, both spectrometers must accept the event.
            // This is the standard mode for (e,e'p) reactions.
            // ================================================================
            if (accepted_by_hms && accepted_by_shms) {
                // Store final reconstructed quantities
                event.e_delta = main.SP_electron.delta;
                event.e_xptar = main.SP_electron.xptar;
                event.e_yptar = main.SP_electron.yptar;
                event.p_delta = main.SP_hadron.delta;
                
                // Increment coincidence counter
                ++ncontribute;
                
                // Write to output file
                output.FillEvent(event, main);
                output.FillHistograms(event, main.weight, false);
                
                if (ngenerated <= 10) {
                    std::cout << "  COINCIDENCE: Event accepted by BOTH spectrometers!" << std::endl;
                }
            }
            
            output.IncrementTried();
        }
        
        // ====================================================================
        // FINALIZATION
        // ====================================================================
        
        output.SetStatistics(ntried, ngenerated, ncontribute);
        
        std::cout << "\nFinalizing output..." << std::endl;
        output.Finalize();
        
        // Get statistics from both spectrometers
        auto hms_stats = hms.GetStats();
        auto shms_stats = shms.GetStats();
        
        double generation_eff = 100.0 * ngenerated / static_cast<double>(ntried);
        double acceptance = 100.0 * ncontribute / static_cast<double>(ngenerated);
        
        std::cout << "\n==========================================\n";
        std::cout << "Phase 5c.4 Complete - Two-Arm Coincidence!\n";
        std::cout << "--------------------------------------------\n";
        std::cout << "Event Statistics:\n";
        std::cout << "  Events tried:       " << ntried << "\n";
        std::cout << "  Events generated:   " << ngenerated 
                  << " (" << std::fixed << std::setprecision(1) << generation_eff << "%)\n";
        std::cout << "  Events accepted:    " << ncontribute 
                  << " (" << acceptance << "%)\n";
        std::cout << "--------------------------------------------\n";
        
        // ================================================================
        // PHASE 5c.4: Print HMS statistics
        // ================================================================
        std::cout << "HMS Transport Statistics (Electron Arm):\n";
        std::cout << "  Total transported:  " << hms_stats.hSTOP_trials << "\n";
        std::cout << "  Accepted:           " << hms_stats.hSTOP_successes << "\n";
        std::cout << "  Acceptance:         " << std::fixed << std::setprecision(1)
                  << (100.0 * hms_stats.hSTOP_successes / hms_stats.hSTOP_trials) << "%\n";
        
        int hms_rejected = hms_stats.hSTOP_trials - hms_stats.hSTOP_successes;
        if (hms_rejected > 0) {
            std::cout << "\n  HMS Rejection Breakdown:\n";
            
            int coll_total = hms_stats.hSTOP_slit_hor + hms_stats.hSTOP_slit_vert + 
                           hms_stats.hSTOP_slit_oct + hms_stats.hSTOP_coll;
            int q1_total = hms_stats.hSTOP_Q1_in + hms_stats.hSTOP_Q1_mid + hms_stats.hSTOP_Q1_out;
            int q2_total = hms_stats.hSTOP_Q2_in + hms_stats.hSTOP_Q2_mid + hms_stats.hSTOP_Q2_out;
            int q3_total = hms_stats.hSTOP_Q3_in + hms_stats.hSTOP_Q3_mid + hms_stats.hSTOP_Q3_out;
            int d1_total = hms_stats.hSTOP_D1_in + hms_stats.hSTOP_D1_out;
            
            if (coll_total > 0)
                std::cout << "    Collimator:       " << coll_total << "\n";
            if (q1_total > 0)
                std::cout << "    Q1 apertures:     " << q1_total << "\n";
            if (q2_total > 0)
                std::cout << "    Q2 apertures:     " << q2_total << "\n";
            if (q3_total > 0)
                std::cout << "    Q3 apertures:     " << q3_total << "\n";
            if (d1_total > 0)
                std::cout << "    Dipole+Pipes:     " << d1_total << "\n";
            if (hms_stats.hSTOP_hut > 0)
                std::cout << "    Hut/Detectors:    " << hms_stats.hSTOP_hut << "\n";
        }
        
        std::cout << "\n";
        
        // ================================================================
        // Print SHMS statistics (from Phase 5c.2)
        // ================================================================
        std::cout << "SHMS Transport Statistics (Hadron Arm):\n";
        std::cout << "  Total transported:  " << shms_stats.total_events << "\n";
        std::cout << "  Accepted:           " << shms_stats.accepted << "\n";
        std::cout << "  Acceptance:         " << std::fixed << std::setprecision(1)
                  << (100.0 * shms_stats.accepted / shms_stats.total_events) << "%\n";
        
        int shms_rejected = shms_stats.total_events - shms_stats.accepted;
        if (shms_rejected > 0) {
            std::cout << "\n  SHMS Rejection Breakdown:\n";
            
            if (shms_stats.shmsSTOP_HB > 0)
                std::cout << "    Holding Box:      " << shms_stats.shmsSTOP_HB << "\n";
            if (shms_stats.shmsSTOP_Q1 > 0)
                std::cout << "    Q1 apertures:     " << shms_stats.shmsSTOP_Q1 << "\n";
            if (shms_stats.shmsSTOP_Q2 > 0)
                std::cout << "    Q2 apertures:     " << shms_stats.shmsSTOP_Q2 << "\n";
            if (shms_stats.shmsSTOP_Q3 > 0)
                std::cout << "    Q3 apertures:     " << shms_stats.shmsSTOP_Q3 << "\n";
            if (shms_stats.shmsSTOP_D1 > 0)
                std::cout << "    Dipole:           " << shms_stats.shmsSTOP_D1 << "\n";
            if (shms_stats.shmsSTOP_COLL > 0)
                std::cout << "    Collimator:       " << shms_stats.shmsSTOP_COLL << "\n";
            if (shms_stats.shmsSTOP_HUT > 0)
                std::cout << "    Hut/Detectors:    " << shms_stats.shmsSTOP_HUT << "\n";
        }
        
        std::cout << "\n";
        
        // ================================================================
        // PHASE 5c.4: Coincidence statistics
        // ================================================================
        std::cout << "Coincidence Statistics:\n";
        std::cout << "  Both accepted:      " << ncontribute << "\n";
        std::cout << "  Coincidence rate:   " << std::fixed << std::setprecision(1)
                  << (100.0 * ncontribute / ngenerated) << "%\n";
        std::cout << "\n";
        std::cout << "Note: Coincidence rate = HMS acceptance × SHMS acceptance\n";
        std::cout << "      (for events that pass core transport)\n";
        
        std::cout << "==========================================\n";
        std::cout << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Simulation terminated.\n" << std::endl;
        return 1;
    }
    
    return 0;
}
