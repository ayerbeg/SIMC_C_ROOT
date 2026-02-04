// tests/test_event_generator_diagnostic.cpp
// DIAGNOSTIC version to trace exactly where event generation fails

#include "simc/EventGenerator.h"
#include "simc/ConfigManager.h"
#include "simc/CrossSection.h"
#include "simc/RandomGenerator.h"
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace simc;

// Create test configuration
nlohmann::json create_test_config() {
    return nlohmann::json{
        {"beam", {
            {"energy", 2200.0},
            {"energy_spread", 0.0}
        }},
        {"target", {
            {"mass", 938.272},
            {"Z", 1},
            {"A", 1},
            {"length", 10.0},
            {"density", 0.071},
            {"x_offset", 0.0},
            {"y_offset", 0.0},
            {"z_offset", 0.0},
            {"angle", 0.0},
            {"raster_pattern", 0},
            {"raster_x", 0.0},
            {"raster_y", 0.0}
        }},
        {"generation", {
            {"reaction_type", "elastic"},
            {"use_energy_loss", false},
            {"use_coulomb", false},
            {"use_radiative", false},
            {"beam_xwidth", 0.01},
            {"beam_ywidth", 0.01},
            {"electron", {
                {"delta_min", -10.0},
                {"delta_max", 22.0},
                {"yptar_min", -0.030},
                {"yptar_max", 0.030},
                {"xptar_min", -0.060},
                {"xptar_max", 0.060},
                {"E_min", 500.0},
                {"E_max", 2000.0}
            }},
            {"hadron", {
                {"delta_min", -10.0},
                {"delta_max", 22.0},
                {"yptar_min", -0.030},
                {"yptar_max", 0.030},
                {"xptar_min", -0.060},
                {"xptar_max", 0.060},
                {"E_min", 500.0},
                {"E_max", 2000.0}
            }}
        }},
        {"spectrometer_electron", {
            {"type", "HMS"},
            {"momentum", 1000.0},
            {"angle", 12.5},
            {"phi", 0.0}
        }},
        {"spectrometer_hadron", {
            {"type", "SHMS"},
            {"momentum", 1000.0},
            {"angle", 50.0},
            {"phi", 0.0}
        }}
    };
}

// Manual step-by-step event generation with diagnostics
void test_step_by_step() {
    std::cout << "========================================" << std::endl;
    std::cout << "Step-by-Step Event Generation Diagnostic" << std::endl;
    std::cout << "========================================" << std::endl;
    
    nlohmann::json config_json = create_test_config();
    
    ConfigManager cfg;
    std::ostringstream oss;
    oss << config_json;
    cfg.LoadFromString(oss.str());
    
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<ElasticCrossSection>();
    
    EventGenerator gen(cfg, random, xs, nullptr);
    gen.Initialize();
    
    std::cout << "\nAttempting single event generation:" << std::endl;
    std::cout << "====================================" << std::endl;
    
    // Try to generate one event with detailed output
    SimcEvent event;
    MainEvent main;
    
    // Step 1: Generate vertex
    std::cout << "\nStep 1: Generate Vertex" << std::endl;
    gen.GenerateVertex(main.target);
    std::cout << "  Vertex: x=" << main.target.x 
              << ", y=" << main.target.y 
              << ", z=" << main.target.z << " cm" << std::endl;
    
    // Step 2: Generate phase space
    std::cout << "\nStep 2: Generate Phase Space" << std::endl;
    bool ps_ok = gen.GeneratePhaseSpace(event, main);
    std::cout << "  Phase space OK: " << (ps_ok ? "YES" : "NO") << std::endl;
    
    if (!ps_ok) {
        std::cout << "  ERROR: Phase space generation failed!" << std::endl;
        return;
    }
    
    std::cout << "  Generated electron arm:" << std::endl;
    std::cout << "    delta  = " << event.e_delta << " %" << std::endl;
    std::cout << "    xptar  = " << event.e_xptar << " rad" << std::endl;
    std::cout << "    yptar  = " << event.e_yptar << " rad" << std::endl;
    std::cout << "    theta  = " << event.e_theta * 180.0/3.14159 << " deg" << std::endl;
    std::cout << "    phi    = " << event.e_phi * 180.0/3.14159 << " deg" << std::endl;
    std::cout << "    P      = " << event.e_P << " MeV/c" << std::endl;
    std::cout << "    E      = " << event.e_E << " MeV" << std::endl;
    std::cout << "  Generated hadron arm:" << std::endl;
    std::cout << "    delta  = " << event.p_delta << " %" << std::endl;
    std::cout << "    xptar  = " << event.p_xptar << " rad" << std::endl;
    std::cout << "    yptar  = " << event.p_yptar << " rad" << std::endl;
    std::cout << "    P      = " << event.p_P << " MeV/c" << std::endl;
    std::cout << "  Beam energy: " << event.Ein << " MeV" << std::endl;
    
    // Step 3: Complete event (solves kinematics)
    std::cout << "\nStep 3: Complete Event (Solve Kinematics)" << std::endl;
    bool complete_ok = gen.CompleteEvent(event, main);
    std::cout << "  Complete OK: " << (complete_ok ? "YES" : "NO") << std::endl;
    
    if (!complete_ok) {
        std::cout << "  ERROR: CompleteEvent failed!" << std::endl;
        std::cout << "  This could be due to SolveHydrogenElastic returning false" << std::endl;
        return;
    }
    
    std::cout << "  After solving elastic kinematics:" << std::endl;
    std::cout << "    e_E      = " << event.e_E << " MeV (RECALCULATED)" << std::endl;
    std::cout << "    e_P      = " << event.e_P << " MeV/c (RECALCULATED)" << std::endl;
    std::cout << "    e_delta  = " << event.e_delta << " % (RECALCULATED)" << std::endl;
    std::cout << "    p_E      = " << event.p_E << " MeV" << std::endl;
    std::cout << "    p_P      = " << event.p_P << " MeV/c" << std::endl;
    std::cout << "    Q2       = " << event.Q2 << " MeV^2" << std::endl;
    std::cout << "    nu       = " << event.nu << " MeV" << std::endl;
    
    // Step 4: Calculate kinematics
    std::cout << "\nStep 4: Calculate Basic Kinematics" << std::endl;
    gen.CalculateBasicKinematics(event);
    std::cout << "  W        = " << event.W << " MeV" << std::endl;
    std::cout << "  xbj      = " << event.xbj << std::endl;
    std::cout << "  epsilon  = " << event.epsilon << std::endl;
    
    // Step 5: Calculate missing momentum
    std::cout << "\nStep 5: Calculate Missing Momentum" << std::endl;
    gen.CalculateMissingMomentum(event);
    std::cout << "  Pmiss    = " << event.Pmiss << " MeV/c" << std::endl;
    std::cout << "  Em       = " << event.Em << " MeV" << std::endl;
    
    // Step 6: Calculate jacobian
    std::cout << "\nStep 6: Calculate Jacobian" << std::endl;
    main.jacobian = gen.CalculateJacobian(event, main);
    std::cout << "  Jacobian = " << main.jacobian << std::endl;
    
    // Step 7: Validate physics
    std::cout << "\nStep 7: Validate Physics" << std::endl;
    
    // Manual validation checks
    bool energy_ok = event.e_E < event.Ein;
    bool electron_mass_ok = event.e_E >= 0.511;
    bool proton_mass_ok = event.p_E >= 938.272;
    bool q2_ok = event.Q2 > 0.0;
    
    std::cout << "  Energy conservation: " << (energy_ok ? "PASS" : "FAIL") << std::endl;
    std::cout << "    e_E (" << event.e_E << ") < Ein (" << event.Ein << "): " 
              << (energy_ok ? "YES" : "NO") << std::endl;
    std::cout << "  Electron mass: " << (electron_mass_ok ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Proton mass: " << (proton_mass_ok ? "PASS" : "FAIL") << std::endl;
    std::cout << "  Q2 positive: " << (q2_ok ? "PASS" : "FAIL") << std::endl;
    
    // Step 8: Check generation cuts
    std::cout << "\nStep 8: Check Generation Cuts" << std::endl;
    std::cout << "  Generation limits:" << std::endl;
    std::cout << "    electron.delta: [-10.0, 22.0]" << std::endl;
    std::cout << "    electron.E:     [500.0, 2000.0]" << std::endl;
    std::cout << "  Event values:" << std::endl;
    std::cout << "    electron.delta: " << event.e_delta << std::endl;
    std::cout << "    electron.E:     " << event.e_E << std::endl;
    
    bool delta_ok = (event.e_delta >= -10.0 && event.e_delta <= 22.0);
    bool E_ok = (event.e_E >= 500.0 && event.e_E <= 2000.0);
    
    std::cout << "  Delta in range: " << (delta_ok ? "PASS" : "FAIL") << std::endl;
    std::cout << "  E in range: " << (E_ok ? "PASS" : "FAIL") << std::endl;
    
    if (!delta_ok) {
        std::cout << "\n*** PROBLEM IDENTIFIED ***" << std::endl;
        std::cout << "The electron delta is OUTSIDE generation limits after" << std::endl;
        std::cout << "SolveHydrogenElastic recalculated it!" << std::endl;
        std::cout << "\nThis happens because:" << std::endl;
        std::cout << "1. Phase space generates random delta in [-10, 22]%" << std::endl;
        std::cout << "2. This gives initial e_E and e_P" << std::endl;
        std::cout << "3. SolveHydrogenElastic recalculates e_E from elastic constraint" << std::endl;
        std::cout << "4. Then recalculates delta = 100*(e_P - P_central)/P_central" << std::endl;
        std::cout << "5. New delta is often outside [-10, 22]% range!" << std::endl;
        std::cout << "\nSOLUTION: For elastic scattering, need to generate in" << std::endl;
        std::cout << "electron angles (theta_e, phi_e) instead of delta!" << std::endl;
    }
    
    if (!E_ok) {
        std::cout << "\n*** PROBLEM IDENTIFIED ***" << std::endl;
        std::cout << "The electron energy is OUTSIDE generation limits!" << std::endl;
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Diagnostic complete" << std::endl;
    std::cout << "========================================" << std::endl;
}

int main() {
    test_step_by_step();
    return 0;
}
