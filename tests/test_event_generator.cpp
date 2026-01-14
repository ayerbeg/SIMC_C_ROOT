// tests/test_event_generator.cpp
// Unit tests for EventGenerator
// CORRECTED to use EXACT config paths from EventGenerator.cpp

#include "simc/EventGenerator.h"
#include "simc/ConfigManager.h"
#include "simc/CrossSection.h"
#include "simc/SpectrometerOptics.h"
#include "simc/RandomGenerator.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <sstream>

using namespace simc;

// Helper function for floating point comparison
bool approx_equal(double a, double b, double eps = 1e-6) {
    return std::abs(a - b) < eps;
}

// Create standard test configuration
nlohmann::json create_test_config() {
    return nlohmann::json{
        {"beam", {
            {"energy", 2200.0},  // MeV
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
            {"reaction_type", "elastic"},  // lowercase!
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
            {"momentum", 1000.0},  // MeV/c
            {"angle", 12.5},       // degrees
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

// Test 1: Constructor and initialization
void test_constructor() {
    std::cout << "Test 1: Constructor and Initialization..." << std::endl;
    
    nlohmann::json config_json = create_test_config();
    
    ConfigManager cfg;
    std::ostringstream oss;
    oss << config_json;
    cfg.LoadFromString(oss.str());
    
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<ElasticCrossSection>();
    
    EventGenerator gen(cfg, random, xs, nullptr);
    assert(gen.GetReactionType() == ReactionType::ELASTIC);
    
    bool init_success = gen.Initialize();
    assert(init_success);
    
    std::cout << "  PASSED" << std::endl;
}

// Test 2: Vertex generation
void test_vertex_generation() {
    std::cout << "Test 2: Vertex Generation..." << std::endl;
    
    nlohmann::json config_json = create_test_config();
    
    ConfigManager cfg;
    std::ostringstream oss;
    oss << config_json;
    cfg.LoadFromString(oss.str());
    
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<ElasticCrossSection>();
    
    EventGenerator gen(cfg, random, xs, nullptr);
    gen.Initialize();
    
    // Generate many vertices and check distribution
    const int N = 10000;
    double sum_x = 0.0, sum_x2 = 0.0;
    double sum_y = 0.0, sum_y2 = 0.0;
    double sum_z = 0.0, sum_z2 = 0.0;
    
    for (int i = 0; i < N; ++i) {
        TargetInfo target;
        gen.GenerateVertex(target);
        
        sum_x += target.x;
        sum_x2 += target.x * target.x;
        sum_y += target.y;
        sum_y2 += target.y * target.y;
        sum_z += target.z;
        sum_z2 += target.z * target.z;
        
        // Check z is within target length
        assert(std::abs(target.z) <= 5.0); // length/2
    }
    
    double mean_x = sum_x / N;
    double mean_y = sum_y / N;
    double mean_z = sum_z / N;
    
    double var_x = sum_x2 / N - mean_x * mean_x;
    double var_y = sum_y2 / N - mean_y * mean_y;
    
    double sigma_x = std::sqrt(var_x);
    double sigma_y = std::sqrt(var_y);
    
    std::cout << "  Mean x: " << mean_x << " cm (expect ~0)" << std::endl;
    std::cout << "  Mean y: " << mean_y << " cm (expect ~0)" << std::endl;
    std::cout << "  Mean z: " << mean_z << " cm (expect ~0)" << std::endl;
    std::cout << "  Sigma x: " << sigma_x << " cm (expect ~0.01)" << std::endl;
    std::cout << "  Sigma y: " << sigma_y << " cm (expect ~0.01)" << std::endl;
    
    // Check means are near zero
    assert(std::abs(mean_x) < 0.001);
    assert(std::abs(mean_y) < 0.001);
    assert(std::abs(mean_z) < 0.1);
    
    // Check sigmas are approximately correct (within 20% due to sampling)
    assert(std::abs(sigma_x - 0.01) < 0.002);
    assert(std::abs(sigma_y - 0.01) < 0.002);
    
    std::cout << "  PASSED" << std::endl;
}

// Test 3: Hydrogen elastic kinematics
void test_hydrogen_elastic() {
    std::cout << "Test 3: Hydrogen Elastic Kinematics..." << std::endl;
    
    nlohmann::json config_json = create_test_config();
    
    ConfigManager cfg;
    std::ostringstream oss;
    oss << config_json;
    cfg.LoadFromString(oss.str());
    
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<ElasticCrossSection>();
    
    EventGenerator gen(cfg, random, xs, nullptr);
    gen.Initialize();
    
    // Generate events
    int n_success = 0;
    const int N = 1000;
    
    std::cout << "  Attempting to generate " << N << " events..." << std::endl;
    
    for (int i = 0; i < N; ++i) {
        SimcEvent event;
        MainEvent main;
        
        if (gen.GenerateEvent(event, main)) {
            n_success++;
            
            // Check energy conservation
            double E_initial = event.Ein + kProtonMass;
            double E_final = event.e_E + event.p_E;
            
            if (n_success == 1) {
                // Print first event details
                std::cout << "  First event details:" << std::endl;
                std::cout << "    E_initial = " << E_initial << " MeV" << std::endl;
                std::cout << "    E_final   = " << E_final << " MeV" << std::endl;
                std::cout << "    E_beam    = " << event.Ein << " MeV" << std::endl;
                std::cout << "    E_e'      = " << event.e_E << " MeV" << std::endl;
                std::cout << "    E_p       = " << event.p_E << " MeV" << std::endl;
                std::cout << "    theta_e   = " << event.e_theta * kRadToDeg << " deg" << std::endl;
                std::cout << "    Q2        = " << event.Q2 << " MeV^2" << std::endl;
                std::cout << "    nu        = " << event.nu << " MeV" << std::endl;
                std::cout << "    P_miss    = " << event.Pmiss << " MeV/c" << std::endl;
            }
            
            // Energy conservation (1 MeV tolerance)
            assert(approx_equal(E_initial, E_final, 1.0));
            
            // Check Q² is positive
            assert(event.Q2 > 0.0);
            
            // Check missing momentum is small for elastic
            assert(event.Pmiss < 10.0); // Should be ~0
        }
    }
    
    std::cout << "  Generated " << n_success << " / " << N << " events" << std::endl;
    
    if (n_success == 0) {
        std::cerr << "  ERROR: No events generated! Check SolveHydrogenElastic implementation." << std::endl;
        std::cerr << "  Generation limits:" << std::endl;
        std::cerr << "    electron.delta: [" << config_json["generation"]["electron"]["delta_min"] 
                  << ", " << config_json["generation"]["electron"]["delta_max"] << "]" << std::endl;
        std::cerr << "    electron.E: [" << config_json["generation"]["electron"]["E_min"] 
                  << ", " << config_json["generation"]["electron"]["E_max"] << "]" << std::endl;
    }
    
    assert(n_success > 0);
    
    std::cout << "  PASSED" << std::endl;
}

// Test 4: Angle conversions
void test_angle_conversions() {
    std::cout << "Test 4: Angle Conversions..." << std::endl;
    
    nlohmann::json config_json = create_test_config();
    config_json["spectrometer_hadron"]["phi"] = 90.0;  // SOS position
    
    ConfigManager cfg;
    std::ostringstream oss;
    oss << config_json;
    cfg.LoadFromString(oss.str());
    
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<ElasticCrossSection>();
    
    EventGenerator gen(cfg, random, xs, nullptr);
    gen.Initialize();
    
    // Test round-trip conversion
    double theta0 = 12.5 * kDegToRad;
    double phi0 = 90.0 * kDegToRad;
    
    double xptar_orig = 0.02; // rad
    double yptar_orig = 0.01; // rad
    
    // Convert to physics angles
    double theta, phi;
    gen.PhysicsAngles(theta0, phi0, xptar_orig, yptar_orig, theta, phi);
    
    // Convert back to spectrometer angles
    double xptar_back, yptar_back;
    gen.SpectrometerAngles(theta0, phi0, theta, phi, xptar_back, yptar_back);
    
    std::cout << "  Original: xptar=" << xptar_orig << ", yptar=" << yptar_orig << std::endl;
    std::cout << "  Physics: theta=" << theta * kRadToDeg << "°, phi=" << phi * kRadToDeg << "°" << std::endl;
    std::cout << "  Back: xptar=" << xptar_back << ", yptar=" << yptar_back << std::endl;
    
    // Check round-trip (should recover original values within tolerance)
    assert(approx_equal(xptar_orig, xptar_back, 1e-4));
    assert(approx_equal(yptar_orig, yptar_back, 1e-4));
    
    std::cout << "  PASSED" << std::endl;
}

// Test 5: Phase space generation
void test_phase_space() {
    std::cout << "Test 5: Phase Space Generation..." << std::endl;
    
    nlohmann::json config_json = create_test_config();
    
    ConfigManager cfg;
    std::ostringstream oss;
    oss << config_json;
    cfg.LoadFromString(oss.str());
    
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<ElasticCrossSection>();
    
    EventGenerator gen(cfg, random, xs, nullptr);
    gen.Initialize();
    
    // Generate phase space and check ranges
    const int N = 1000;
    int n_in_range = 0;
    
    for (int i = 0; i < N; ++i) {
        SimcEvent event;
        MainEvent main;
        
        if (gen.GeneratePhaseSpace(event, main)) {
            n_in_range++;
            
            // Check electron arm is in range
            assert(event.e_delta >= -10.0 && event.e_delta <= 22.0);
            assert(event.e_xptar >= -0.060 && event.e_xptar <= 0.060);
            assert(event.e_yptar >= -0.030 && event.e_yptar <= 0.030);
            
            // Check hadron arm is in range
            assert(event.p_delta >= -10.0 && event.p_delta <= 22.0);
            assert(event.p_xptar >= -0.060 && event.p_xptar <= 0.060);
            assert(event.p_yptar >= -0.030 && event.p_yptar <= 0.030);
        }
    }
    
    std::cout << "  Generated " << n_in_range << " / " << N << " in acceptance" << std::endl;
    assert(n_in_range > 0);
    
    std::cout << "  PASSED" << std::endl;
}

// Test 6: Statistics tracking
void test_statistics() {
    std::cout << "Test 6: Statistics Tracking..." << std::endl;
    
    nlohmann::json config_json = create_test_config();
    
    ConfigManager cfg;
    std::ostringstream oss;
    oss << config_json;
    cfg.LoadFromString(oss.str());
    
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<ElasticCrossSection>();
    
    EventGenerator gen(cfg, random, xs, nullptr);
    gen.Initialize();
    
    // Reset statistics
    gen.ResetStatistics();
    assert(gen.GetNGenerated() == 0);
    assert(gen.GetNAccepted() == 0);
    
    // Generate some events
    const int N = 100;
    for (int i = 0; i < N; ++i) {
        SimcEvent event;
        MainEvent main;
        gen.GenerateEvent(event, main);
    }
    
    // Check statistics are tracked
    assert(gen.GetNGenerated() == N);
    assert(gen.GetNAccepted() > 0);
    
    std::cout << "  Generated: " << gen.GetNGenerated() << std::endl;
    std::cout << "  Accepted: " << gen.GetNAccepted() << std::endl;
    std::cout << "  Efficiency: " << 100.0 * gen.GetNAccepted() / gen.GetNGenerated() << "%" << std::endl;
    
    std::cout << "  PASSED" << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "EventGenerator Unit Tests" << std::endl;
    std::cout << "========================================" << std::endl;
    
    try {
        test_constructor();
        test_vertex_generation();
        test_hydrogen_elastic();
        test_angle_conversions();
        test_phase_space();
        test_statistics();
        
        std::cout << "========================================" << std::endl;
        std::cout << "All tests PASSED!" << std::endl;
        std::cout << "========================================" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "========================================" << std::endl;
        std::cerr << "Test FAILED with exception: " << e.what() << std::endl;
        std::cerr << "========================================" << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "========================================" << std::endl;
        std::cerr << "Test FAILED with unknown exception" << std::endl;
        std::cerr << "========================================" << std::endl;
        return 1;
    }
}
