// src/main.cpp - Phase 5e: Focal Plane Integration
// COMPLETE MODIFIED VERSION - Week 1, Days 2-5
// Adds focal plane coordinate extraction after spectrometer transport

#include "simc/core/SimcEvent.h"
#include "simc/core/ConfigManager.h"
#include "simc/core/RandomGenerator.h"
#include "simc/io/OutputManager.h"
#include "simc/physics/EventGenerator.h"
#include "simc/transport/MonteCarloTransport.h"
#include "simc/physics/CrossSection.h"
#include "simc/core/SimcConstants.h"

// All 5 spectrometers
#include "simc/spectrometers/HMS.h"
#include "simc/spectrometers/SOS.h"
#include "simc/spectrometers/SHMS.h"
#include "simc/spectrometers/HRSr.h"
#include "simc/spectrometers/HRSl.h"

#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <string>
#include <stdexcept>

using namespace simc;
using namespace simc::constants;


// ============================================================================
// Spectrometer Transport Wrapper - PHASE 5e MODIFIED
// ============================================================================
struct SpectrometerTransport {
    int electron_arm_id{5};  // Default: SHMS (id=5)
    int hadron_arm_id{5};
    
    // All 5 spectrometer instances
    std::unique_ptr<HMS> hms;
    std::unique_ptr<SOS> sos;
    std::unique_ptr<SHMS> shms;
    std::unique_ptr<HRSr> hrsr;
    std::unique_ptr<HRSl> hrsl;
    
    bool Initialize(const ConfigManager& config) {
        // Get spectrometer configuration from config
        electron_arm_id = config.Get<int>("spectrometers.electron_arm", 5);  // Default: SHMS
        hadron_arm_id = config.Get<int>("spectrometers.hadron_arm", 5);      // Default: SHMS
        
        std::cout << "Initializing spectrometers..." << std::endl;
        std::cout << "  Electron arm: " << GetSpecName(electron_arm_id) << " (id=" << electron_arm_id << ")" << std::endl;
        std::cout << "  Hadron arm:   " << GetSpecName(hadron_arm_id) << " (id=" << hadron_arm_id << ")" << std::endl;
        
        // Initialize required spectrometers
        bool success = true;
        
        if (NeedsSpec(1)) {  // HMS
            hms = std::make_unique<HMS>();
            success &= LoadMatrices(hms.get(), "hms", config);
        }
        
        if (NeedsSpec(2)) {  // SOS
            sos = std::make_unique<SOS>();
            success &= LoadMatrices(sos.get(), "sos", config);
        }
        
        if (NeedsSpec(5) || NeedsSpec(6)) {  // SHMS
            shms = std::make_unique<SHMS>();
            success &= LoadMatrices(shms.get(), "shms", config);
        }
        
        if (NeedsSpec(3)) {  // HRSr
            hrsr = std::make_unique<HRSr>();
            success &= LoadMatrices(hrsr.get(), "hrsr", config);
        }
        
        if (NeedsSpec(4)) {  // HRSl
            hrsl = std::make_unique<HRSl>();
            success &= LoadMatrices(hrsl.get(), "hrsl", config);
        }
        
        return success;
    }
    
    // Transport electron through its spectrometer
    bool TransportElectron(const SimcEvent& event, const MainEvent& main) {
        return TransportParticle(event, main, electron_arm_id, true);
    }
    
    // Transport hadron through its spectrometer
    bool TransportHadron(const SimcEvent& event, const MainEvent& main) {
        return TransportParticle(event, main, hadron_arm_id, false);
    }
    
    // ========================================================================
    // PHASE 5e: Focal Plane Extraction Methods
    // ========================================================================
    
    /**
     * @brief Extract electron focal plane coordinates after successful transport
     * @param main MainEvent structure to fill with focal plane data
     * 
     * Converts spectrometer-specific focal plane states to global FocalPlaneState:
     * - HMS/SHMS: xp,yp in mrad → convert to dx,dy slopes (*1/1000)
     * - SOS/HRS: dx,dy already slopes → copy directly
     */
    void ExtractElectronFocalPlane(MainEvent& main) {
        switch (electron_arm_id) {
            case 1:  // HMS
                if (hms) {
                    HMS::FocalPlaneState fp;
                    hms->GetFocalPlane(last_hms_track_, fp);
                    main.FP_electron.x = fp.x;
                    main.FP_electron.y = fp.y;
                    main.FP_electron.dx = fp.xp / 1000.0;  // mrad → slope
                    main.FP_electron.dy = fp.yp / 1000.0;  // mrad → slope
                    main.FP_electron.path = last_hms_track_.pathlen;
                }
                break;
            
            case 2:  // SOS
                if (sos) {
                    SOS::FocalPlaneState fp;
                    sos->GetFocalPlane(last_sos_track_, fp);
                    main.FP_electron.x = fp.x;
                    main.FP_electron.y = fp.y;
                    main.FP_electron.dx = fp.dx;  // already slope
                    main.FP_electron.dy = fp.dy;  // already slope
                    main.FP_electron.path = last_sos_track_.pathlen;
                }
                break;
            
            case 5:  // SHMS
            case 6:  // SHMS_P5
                if (shms) {
                    SHMS::FocalPlaneState fp;
                    shms->GetFocalPlane(last_shms_track_, fp);
                    main.FP_electron.x = fp.x;
                    main.FP_electron.y = fp.y;
                    main.FP_electron.dx = fp.xp / 1000.0;  // mrad → slope
                    main.FP_electron.dy = fp.yp / 1000.0;  // mrad → slope
                    main.FP_electron.path = last_shms_track_.pathlen;
                }
                break;
            
            case 3:  // HRSr
                if (hrsr) {
                    HRSr::FocalPlaneState fp;
                    hrsr->GetFocalPlane(last_hrsr_track_, fp);
                    main.FP_electron.x = fp.x;
                    main.FP_electron.y = fp.y;
                    main.FP_electron.dx = fp.dx;
                    main.FP_electron.dy = fp.dy;
                    main.FP_electron.path = last_hrsr_track_.pathlen;
                }
                break;
            
            case 4:  // HRSl
                if (hrsl) {
                    HRSl::FocalPlaneState fp;
                    hrsl->GetFocalPlane(last_hrsl_track_, fp);
                    main.FP_electron.x = fp.x;
                    main.FP_electron.y = fp.y;
                    main.FP_electron.dx = fp.dx;
                    main.FP_electron.dy = fp.dy;
                    main.FP_electron.path = last_hrsl_track_.pathlen;
                }
                break;
        }
    }
    
    /**
     * @brief Extract hadron focal plane coordinates after successful transport
     */
    void ExtractHadronFocalPlane(MainEvent& main) {
        switch (hadron_arm_id) {
            case 1:  // HMS
                if (hms) {
                    HMS::FocalPlaneState fp;
                    hms->GetFocalPlane(last_hms_track_, fp);
                    main.FP_hadron.x = fp.x;
                    main.FP_hadron.y = fp.y;
                    main.FP_hadron.dx = fp.xp / 1000.0;
                    main.FP_hadron.dy = fp.yp / 1000.0;
                    main.FP_hadron.path = last_hms_track_.pathlen;
                }
                break;
            
            case 2:  // SOS
                if (sos) {
                    SOS::FocalPlaneState fp;
                    sos->GetFocalPlane(last_sos_track_, fp);
                    main.FP_hadron.x = fp.x;
                    main.FP_hadron.y = fp.y;
                    main.FP_hadron.dx = fp.dx;
                    main.FP_hadron.dy = fp.dy;
                    main.FP_hadron.path = last_sos_track_.pathlen;
                }
                break;
            
            case 5:  // SHMS
            case 6:  // SHMS_P5
                if (shms) {
                    SHMS::FocalPlaneState fp;
                    shms->GetFocalPlane(last_shms_track_, fp);
                    main.FP_hadron.x = fp.x;
                    main.FP_hadron.y = fp.y;
                    main.FP_hadron.dx = fp.xp / 1000.0;
                    main.FP_hadron.dy = fp.yp / 1000.0;
                    main.FP_hadron.path = last_shms_track_.pathlen;
                }
                break;
            
            case 3:  // HRSr
                if (hrsr) {
                    HRSr::FocalPlaneState fp;
                    hrsr->GetFocalPlane(last_hrsr_track_, fp);
                    main.FP_hadron.x = fp.x;
                    main.FP_hadron.y = fp.y;
                    main.FP_hadron.dx = fp.dx;
                    main.FP_hadron.dy = fp.dy;
                    main.FP_hadron.path = last_hrsr_track_.pathlen;
                }
                break;
            
            case 4:  // HRSl
                if (hrsl) {
                    HRSl::FocalPlaneState fp;
                    hrsl->GetFocalPlane(last_hrsl_track_, fp);
                    main.FP_hadron.x = fp.x;
                    main.FP_hadron.y = fp.y;
                    main.FP_hadron.dx = fp.dx;
                    main.FP_hadron.dy = fp.dy;
                    main.FP_hadron.path = last_hrsl_track_.pathlen;
                }
                break;
        }
    }
    
    // Print statistics for all active spectrometers
    void PrintStatistics(std::ostream& os) const {
        os << "\n==========================================\n";
        os << "Spectrometer Transport Statistics:\n";
        os << "==========================================\n";
        
        if (hms) PrintHMSStats(os);
        if (sos) PrintSOSStats(os);
        if (shms) PrintSHMSStats(os);
        if (hrsr) PrintHRSrStats(os);
        if (hrsl) PrintHRSlStats(os);
    }
    
private:
    // ========================================================================
    // PHASE 5e: Track storage for focal plane extraction
    // ========================================================================
    HMS::TrackState last_hms_track_;
    SOS::TrackState last_sos_track_;
    SHMS::TrackState last_shms_track_;
    HRSr::TrackState last_hrsr_track_;
    HRSl::TrackState last_hrsl_track_;
    
    bool NeedsSpec(int spec_id) const {
        return electron_arm_id == spec_id || hadron_arm_id == spec_id;
    }
    
    std::string GetSpecName(int spec_id) const {
        switch (spec_id) {
            case 1: return "HMS";
            case 2: return "SOS";
            case 3: return "HRSr";
            case 4: return "HRSl";
            case 5: return "SHMS";
            case 6: return "SHMS (P5)";
            default: return "NONE";
        }
    }
    
    template<typename T>
    bool LoadMatrices(T* spec, const std::string& name, const ConfigManager& config) {
        std::string data_dir = config.Get<std::string>("paths.data_dir", "data");
        std::string forward, recon;
        
        // Determine correct filenames based on spectrometer
        if (name == "hms" || name == "sos") {
            // HMS and SOS use: forward_cosy.dat, recon_cosy.dat
            forward = data_dir + "/matrices/" + name + "/forward_cosy.dat";
            recon = data_dir + "/matrices/" + name + "/recon_cosy.dat";
        } else if (name == "hrsr" || name == "hrsl") {
            // HRSr and HRSl use: hrs_forward_cosy.dat, hrs_recon_cosy.dat
            forward = data_dir + "/matrices/" + name + "/hrs_forward_cosy.dat";
            recon = data_dir + "/matrices/" + name + "/hrs_recon_cosy.dat";
        } else if (name == "shms") {
            // SHMS uses: shms_forward.dat, shms_recon.dat
            forward = data_dir + "/matrices/" + name + "/shms_forward.dat";
            recon = data_dir + "/matrices/" + name + "/shms_recon.dat";
        } else {
            std::cerr << "  ERROR: Unknown spectrometer name: " << name << std::endl;
            return false;
        }
        
        std::cout << "  Loading " << name << " matrices..." << std::endl;
        if (!spec->LoadMatrices(forward, recon)) {
            std::cerr << "  ERROR: Failed to load " << name << " matrix files" << std::endl;
            return false;
        }
        std::cout << "  " << name << " matrices loaded successfully" << std::endl;
        return true;
    }
    
    bool TransportParticle(const SimcEvent& event, const MainEvent& main, 
                          int spec_id, bool is_electron) {
        // Get spectrometer quantities from MainEvent
        const auto& SP = is_electron ? main.SP_electron : main.SP_hadron;
        
        switch (spec_id) {
            case 1:  // HMS
                return TransportHMS(event, SP);
            case 2:  // SOS
                return TransportSOS(event, SP);
            case 5:  // SHMS
            case 6:  // SHMS_P5
                return TransportSHMS(event, SP);
            case 3:  // HRSr
                return TransportHRSr(event, SP);
            case 4:  // HRSl
                return TransportHRSl(event, SP);
            default:
                return true;  // No spectrometer = accept
        }
    }
    
    // ========================================================================
    // PHASE 5e MODIFIED: Save track after successful transport
    // ========================================================================
    
    bool TransportHMS(const SimcEvent& event, const ArmState& SP) {
        HMS::TrackState track;
        track.x = 0.0;      // At target
        track.y = 0.0;
        track.z = 0.0;
        track.dx = SP.xptar;  // Slopes (dimensionless)
        track.dy = SP.yptar;
        track.delta = SP.delta;  // %
        track.p = event.p_P / 1000.0;  // MeV/c → GeV/c
        track.m2 = 0.938272 * 0.938272;  // Proton mass² (GeV²)
        
        bool success = hms->Transport(track);
        if (success) {
            last_hms_track_ = track;  // Phase 5e: Save for focal plane extraction
        }
        return success;
    }
    
    bool TransportSOS(const SimcEvent& event, const ArmState& SP) {
        SOS::TrackState track;
        track.x = 0.0;
        track.y = 0.0;
        track.z = 0.0;
        track.dx = SP.xptar;
        track.dy = SP.yptar;
        track.delta = SP.delta;
        track.p = event.p_P / 1000.0;  // MeV/c → GeV/c
        track.m2 = 0.938272 * 0.938272;
        
        bool success = sos->Transport(track);
        if (success) {
            last_sos_track_ = track;  // Phase 5e
        }
        return success;
    }
    
    bool TransportSHMS(const SimcEvent& event, const ArmState& SP) {
        SHMS::TrackState track;
        track.x = 0.0;
        track.y = 0.0;
        track.z = 0.0;
        track.dx = SP.xptar;
        track.dy = SP.yptar;
        track.delta = SP.delta;
        track.p = event.p_P / 1000.0;  // MeV/c → GeV/c
        track.m2 = 0.938272 * 0.938272;
        
        bool success = shms->Transport(track);
        if (success) {
            last_shms_track_ = track;  // Phase 5e
        }
        return success;
    }
    
    bool TransportHRSr(const SimcEvent& event, const ArmState& SP) {
        HRSr::TrackState track;
        track.x = 0.0;
        track.y = 0.0;
        track.z = 0.0;
        track.dx = SP.xptar;
        track.dy = SP.yptar;
        track.delta = SP.delta;
        track.p = event.p_P / 1000.0;  // MeV/c → GeV/c
        track.m2 = 0.938272 * 0.938272;
        
        bool success = hrsr->Transport(track);
        if (success) {
            last_hrsr_track_ = track;  // Phase 5e
        }
        return success;
    }
    
    bool TransportHRSl(const SimcEvent& event, const ArmState& SP) {
        HRSl::TrackState track;
        track.x = 0.0;
        track.y = 0.0;
        track.z = 0.0;
        track.dx = SP.xptar;
        track.dy = SP.yptar;
        track.delta = SP.delta;
        track.p = event.p_P / 1000.0;  // MeV/c → GeV/c
        track.m2 = 0.938272 * 0.938272;
        
        bool success = hrsl->Transport(track);
        if (success) {
            last_hrsl_track_ = track;  // Phase 5e
        }
        return success;
    }
    
    void PrintHMSStats(std::ostream& os) const {
        const auto& stats = hms->GetStats();
        os << "\nHMS Statistics:\n";
        os << "  Trials:       " << stats.hSTOP_trials << "\n";
        os << "  Successes:    " << stats.hSTOP_successes << "\n";
        os << "  Acceptance:   " << std::fixed << std::setprecision(2) 
           << (100.0 * stats.hSTOP_successes / std::max(1, stats.hSTOP_trials)) << "%\n";
        os << "  Rejections:\n";
        if (stats.hSTOP_slit_hor + stats.hSTOP_slit_vert + stats.hSTOP_slit_oct > 0)
            os << "    Collimator:   " << (stats.hSTOP_slit_hor + stats.hSTOP_slit_vert + stats.hSTOP_slit_oct) << "\n";
        if (stats.hSTOP_Q1_in + stats.hSTOP_Q1_mid + stats.hSTOP_Q1_out > 0)
            os << "    Q1:           " << (stats.hSTOP_Q1_in + stats.hSTOP_Q1_mid + stats.hSTOP_Q1_out) << "\n";
        if (stats.hSTOP_Q2_in + stats.hSTOP_Q2_mid + stats.hSTOP_Q2_out > 0)
            os << "    Q2:           " << (stats.hSTOP_Q2_in + stats.hSTOP_Q2_mid + stats.hSTOP_Q2_out) << "\n";
        if (stats.hSTOP_Q3_in + stats.hSTOP_Q3_mid + stats.hSTOP_Q3_out > 0)
            os << "    Q3:           " << (stats.hSTOP_Q3_in + stats.hSTOP_Q3_mid + stats.hSTOP_Q3_out) << "\n";
        if (stats.hSTOP_D1_in + stats.hSTOP_D1_out > 0)
            os << "    Dipole:       " << (stats.hSTOP_D1_in + stats.hSTOP_D1_out) << "\n";
    }
    
    void PrintSOSStats(std::ostream& os) const {
        const auto& stats = sos->GetStats();
        os << "\nSOS Statistics:\n";
        os << "  Trials:       " << stats.sSTOP_trials << "\n";
        os << "  Successes:    " << stats.sSTOP_successes << "\n";
        os << "  Acceptance:   " << std::fixed << std::setprecision(2)
           << (100.0 * stats.sSTOP_successes / std::max(1, stats.sSTOP_trials)) << "%\n";
    }
    
    void PrintSHMSStats(std::ostream& os) const {
        const auto& stats = shms->GetStats();
        os << "\nSHMS Statistics:\n";
        os << "  Trials:       " << stats.total_events << "\n";
        os << "  Successes:    " << stats.accepted << "\n";
        os << "  Acceptance:   " << std::fixed << std::setprecision(2)
           << (100.0 * stats.accepted / std::max(1, stats.total_events)) << "%\n";
        int total_rej = stats.total_events - stats.accepted;
        if (total_rej > 0) {
            os << "  Rejections:\n";
            if (stats.shmsSTOP_HB > 0) os << "    Holding Box:  " << stats.shmsSTOP_HB << "\n";
            if (stats.shmsSTOP_Q1 > 0) os << "    Q1:           " << stats.shmsSTOP_Q1 << "\n";
            if (stats.shmsSTOP_Q2 > 0) os << "    Q2:           " << stats.shmsSTOP_Q2 << "\n";
            if (stats.shmsSTOP_Q3 > 0) os << "    Q3:           " << stats.shmsSTOP_Q3 << "\n";
            if (stats.shmsSTOP_D1 > 0) os << "    Dipole:       " << stats.shmsSTOP_D1 << "\n";
        }
    }
    
    void PrintHRSrStats(std::ostream& os) const {
        const auto& stats = hrsr->GetStats();
        os << "\nHRSr Statistics:\n";
        os << "  Trials:       " << stats.rSTOP_trials << "\n";
        os << "  Successes:    " << stats.rSTOP_successes << "\n";
        os << "  Acceptance:   " << std::fixed << std::setprecision(2)
           << (100.0 * stats.rSTOP_successes / std::max(1, stats.rSTOP_trials)) << "%\n";
    }
    
    void PrintHRSlStats(std::ostream& os) const {
        const auto& stats = hrsl->GetStats();
        os << "\nHRSl Statistics:\n";
        os << "  Trials:       " << stats.lSTOP_trials << "\n";
        os << "  Successes:    " << stats.lSTOP_successes << "\n";
        os << "  Acceptance:   " << std::fixed << std::setprecision(2)
           << (100.0 * stats.lSTOP_successes / std::max(1, stats.lSTOP_trials)) << "%\n";
    }
};


// ============================================================================
// Main Program
// ============================================================================
int main(int argc, char** argv) {
    std::cout << "\n";
    std::cout << "==========================================\n";
    std::cout << "   SIMC C++/ROOT Monte Carlo - Phase 5e\n";
    std::cout << "   Focal Plane Integration              \n";
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
        
        // Initialize subsystems
        auto rng = std::make_shared<RandomGenerator>(seed);
        auto cross_section = std::make_shared<ElasticCrossSection>();
        std::shared_ptr<SpectrometerOptics> optics = nullptr;
        
        // Initialize spectrometers
        SpectrometerTransport spectrometers;
        if (!spectrometers.Initialize(config)) {
            throw std::runtime_error("Failed to initialize spectrometers");
        }
        std::cout << "\n";
        
        // Initialize event generator
        std::cout << "Initializing event generator..." << std::endl;
        EventGenerator generator(config, rng, cross_section, optics);
        if (!generator.Initialize()) {
            throw std::runtime_error("Failed to initialize event generator");
        }
        std::cout << "Event generator initialized.\n" << std::endl;
        
        // Initialize transport (multiple scattering, energy loss)
        MonteCarloTransport transport(config, rng);
        
        // Initialize output
        OutputManager output(output_filename);
        output.Initialize();
        
        // ====================================================================
        // EVENT GENERATION LOOP
        // ====================================================================
        
        std::cout << "Generating events..." << std::endl;
        std::cout << "--------------------------------------------" << std::endl;
        
        int ntried = 0;
        int ngenerated = 0;
        int ncontribute = 0;
        
        generator.ResetStatistics();
        
        while (ngenerated < nevents) {
            ++ntried;
            
            // Progress update
            if (ntried % 5000 == 1) {
                std::cout << "  Event " << ntried 
                          << " (" << ngenerated << " generated, "
                          << ncontribute << " accepted)" << std::endl;
            }
            
            SimcEvent event;
            MainEvent main;
            
            // Initialize main event
            main.gen_weight = 1.0;
            main.jacobian = 1.0;
            main.weight = 1.0;
            main.sigcc = 1.0;
            main.success = false;
            
            // ================================================================
            // 1. Generate physics event
            // ================================================================
            bool success = generator.GenerateEvent(event, main);
            if (!success) continue;
            
            // Calculate cross section
            main.sigcc = cross_section->Calculate(event);
            main.weight = main.gen_weight * main.jacobian * main.sigcc;
            
            ++ngenerated;
            main.success = true;
            
            // Fill generated distributions
            output.FillGenerated(event);
            output.FillHistograms(event, main.weight, true);
            
            // ================================================================
            // 2. Calculate hadron spectrometer angles from physics
            // ================================================================
            transport.CalculateReconstructed(event);
            
            // ================================================================
            // 3. Core Monte Carlo Transport (energy loss, multiple scattering)
            // ================================================================
            bool passes_core_transport = transport.Transport(event, main);
            if (!passes_core_transport) continue;
            
            // ================================================================
            // 4. Spectrometer Transport - BOTH ARMS
            // ================================================================
            
            // Transport electron
            bool electron_ok = spectrometers.TransportElectron(event, main);
            if (!electron_ok) continue;
            
            // ================================================================
            // PHASE 5e: Extract electron focal plane coordinates
            // ================================================================
            spectrometers.ExtractElectronFocalPlane(main);
            
            // Transport hadron
            bool hadron_ok = spectrometers.TransportHadron(event, main);
            if (!hadron_ok) continue;
            
            // ================================================================
            // PHASE 5e: Extract hadron focal plane coordinates
            // ================================================================
            spectrometers.ExtractHadronFocalPlane(main);
            
            // ================================================================
            // 5. Event passed all cuts - fill accepted distributions
            // ================================================================
            
            // Update event with final spectrometer quantities
            event.e_delta = main.SP_electron.delta;
            event.e_xptar = main.SP_electron.xptar;
            event.e_yptar = main.SP_electron.yptar;
            event.p_delta = main.SP_hadron.delta;
            event.p_xptar = main.SP_hadron.xptar;
            event.p_yptar = main.SP_hadron.yptar;
            
            ++ncontribute;
            output.FillEvent(event, main);
            output.FillHistograms(event, main.weight, false);
            output.IncrementTried();
        }
        
        // ====================================================================
        // FINALIZATION
        // ====================================================================
        
        output.SetStatistics(ntried, ngenerated, ncontribute);
        
        std::cout << "\nFinalizing output..." << std::endl;
        output.Finalize();
        
        // ====================================================================
        // SUMMARY STATISTICS
        // ====================================================================
        
        double generation_eff = 100.0 * ngenerated / static_cast<double>(ntried);
        double acceptance = 100.0 * ncontribute / static_cast<double>(ngenerated);
        
        std::cout << "\n==========================================\n";
        std::cout << "Phase 5e Integration Complete!\n";
        std::cout << "--------------------------------------------\n";
        std::cout << "Event Statistics:\n";
        std::cout << "  Events tried:       " << ntried << "\n";
        std::cout << "  Events generated:   " << ngenerated 
                  << " (" << std::fixed << std::setprecision(1) 
                  << generation_eff << "%)\n";
        std::cout << "  Events accepted:    " << ncontribute 
                  << " (" << acceptance << "%)\n";
        
        // Print spectrometer statistics
        spectrometers.PrintStatistics(std::cout);
        
        std::cout << "==========================================\n";
        std::cout << "Output written to: " << output_filename << "\n";
        std::cout << "Focal plane data saved in FP_e_* and FP_p_* branches\n";
        std::cout << "==========================================\n";
        std::cout << "\n";
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        std::cerr << "Simulation terminated.\n" << std::endl;
        return 1;
    }
    
    return 0;
}
