// tests/test_event_generator.cpp
// Unit tests for EventGenerator

#include "simc/EventGenerator.h"
#include "simc/ConfigManager.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

using namespace simc;

// Helper function for floating point comparison
bool approx_equal(double a, double b, double eps = 1e-6) {
    return std::abs(a - b) < eps;
}

// Test 1: Constructor and initialization
void test_constructor() {
    std::cout << "Test 1: Constructor and Initialization..." << std::endl;
    
    // Create minimal config
    nlohmann::json config = {
        {"generation", {
            {"reaction_type", "H(e,e'p)"},
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
        {"target", {
            {"mass", 938.272},
            {"Z", 1},
            {"A", 1},
            {"length", 10.0},
            {"x_offset", 0.0},
            {"y_offset", 0.0},
            {"z_offset", 0.0},
            {"angle", 0.0},
            {"raster_pattern", 0},
            {"raster_x", 0.0},
            {"raster_y", 0.0}
        }},
        {"beam", {
            {"energy", 2200.0},
            {"energy_spread", 0.0}
        }},
        {"spectrometer_electron", {
            {"momentum", 1000.0},
            {"angle", 12.5},
            {"phi", 0.0}
        }},
        {"spectrometer_hadron", {
            {"momentum", 1000.0},
            {"angle", 50.0},
            {"phi", 0.0}
        }}
    };
    
    ConfigManager cfg(config);
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<CrossSection>();
    auto optics = std::make_shared<SpectrometerOptics>();
    
    EventGenerator gen(cfg, random, xs, optics);
    assert(gen.GetReactionType() == ReactionType::kHydElastic);
    
    bool init_success = gen.Initialize();
    assert(init_success);
    
    std::cout << "  PASSED" << std::endl;
}

// Test 2: Vertex generation
void test_vertex_generation() {
    std::cout << "Test 2: Vertex Generation..." << std::endl;
    
    nlohmann::json config = {
        {"generation", {
            {"reaction_type", "H(e,e'p)"},
            {"use_energy_loss", false},
            {"use_coulomb", false},
            {"beam_xwidth", 0.01},
            {"beam_ywidth", 0.01},
            {"electron", {
                {"delta_min", -10.0}, {"delta_max", 22.0},
                {"yptar_min", -0.030}, {"yptar_max", 0.030},
                {"xptar_min", -0.060}, {"xptar_max", 0.060},
                {"E_min", 500.0}, {"E_max", 2000.0}
            }},
            {"hadron", {
                {"delta_min", -10.0}, {"delta_max", 22.0},
                {"yptar_min", -0.030}, {"yptar_max", 0.030},
                {"xptar_min", -0.060}, {"xptar_max", 0.060},
                {"E_min", 500.0}, {"E_max", 2000.0}
            }}
        }},
        {"target", {
            {"mass", 938.272}, {"Z", 1}, {"A", 1},
            {"length", 10.0},
            {"x_offset", 0.0}, {"y_offset", 0.0}, {"z_offset", 0.0},
            {"angle", 0.0},
            {"raster_pattern", 0},
            {"raster_x", 0.0}, {"raster_y", 0.0}
        }},
        {"beam", {{"energy", 2200.0}, {"energy_spread", 0.0}}},
        {"spectrometer_electron", {{"momentum", 1000.0}, {"angle", 12.5}, {"phi", 0.0}}},
        {"spectrometer_hadron", {{"momentum", 1000.0}, {"angle", 50.0}, {"phi", 0.0}}}
    };
    
    ConfigManager cfg(config);
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<CrossSection>();
    auto optics = std::make_shared<SpectrometerOptics>();
    
    EventGenerator gen(cfg, random, xs, optics);
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
    
    // Check sigmas are approximately correct (within 10%)
    assert(std::abs(sigma_x - 0.01) < 0.002);
    assert(std::abs(sigma_y - 0.01) < 0.002);
    
    std::cout << "  PASSED" << std::endl;
}

// Test 3: Hydrogen elastic kinematics
void test_hydrogen_elastic() {
    std::cout << "Test 3: Hydrogen Elastic Kinematics..." << std::endl;
    
    nlohmann::json config = {
        {"generation", {
            {"reaction_type", "H(e,e'p)"},
            {"use_energy_loss", false},
            {"use_coulomb", false},
            {"beam_xwidth", 0.01},
            {"beam_ywidth", 0.01},
            {"electron", {
                {"delta_min", -10.0}, {"delta_max", 22.0},
                {"yptar_min", -0.030}, {"yptar_max", 0.030},
                {"xptar_min", -0.060}, {"xptar_max", 0.060},
                {"E_min", 500.0}, {"E_max", 2000.0}
            }},
            {"hadron", {
                {"delta_min", -10.0}, {"delta_max", 22.0},
                {"yptar_min", -0.030}, {"yptar_max", 0.030},
                {"xptar_min", -0.060}, {"xptar_max", 0.060},
                {"E_min", 500.0}, {"E_max", 2000.0}
            }}
        }},
        {"target", {
            {"mass", 938.272}, {"Z", 1}, {"A", 1},
            {"length", 10.0},
            {"x_offset", 0.0}, {"y_offset", 0.0}, {"z_offset", 0.0},
            {"angle", 0.0},
            {"raster_pattern", 0},
            {"raster_x", 0.0}, {"raster_y", 0.0}
        }},
        {"beam", {{"energy", 2200.0}, {"energy_spread", 0.0}}},
        {"spectrometer_electron", {{"momentum", 1000.0}, {"angle", 12.5}, {"phi", 0.0}}},
        {"spectrometer_hadron", {{"momentum", 1000.0}, {"angle", 50.0}, {"phi", 0.0}}}
    };
    
    ConfigManager cfg(config);
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<CrossSection>();
    auto optics = std::make_shared<SpectrometerOptics>();
    
    EventGenerator gen(cfg, random, xs, optics);
    gen.Initialize();
    
    // Generate events
    int n_success = 0;
    const int N = 1000;
    
    for (int i = 0; i < N; ++i) {
        SimcEvent vertex;
        MainEvent main;
        
        if (gen.GenerateEvent(vertex, main)) {
            n_success++;
            
            // Check energy conservation
            double E_initial = vertex.Ein + kProtonMass;
            double E_final = vertex.e_E + vertex.p_E;
            assert(approx_equal(E_initial, E_final, 1.0)); // 1 MeV tolerance
            
            // Check momentum conservation (proton = q for elastic)
            assert(approx_equal(vertex.p_P, vertex.q, 1.0));
            
            // Check Q² calculation
            double Q2_calc = 2.0 * vertex.Ein * vertex.e_E * (1.0 - vertex.ue_z);
            assert(approx_equal(vertex.Q2, Q2_calc, 1.0));
            
            // Check missing momentum is zero
            assert(vertex.Pmiss < 1.0); // Should be ~0
            assert(std::abs(vertex.Em) < 1.0); // Should be ~0
        }
    }
    
    std::cout << "  Generated " << n_success << " / " << N << " events" << std::endl;
    assert(n_success > 0);
    
    std::cout << "  PASSED" << std::endl;
}

// Test 4: Angle conversions
void test_angle_conversions() {
    std::cout << "Test 4: Angle Conversions..." << std::endl;
    
    nlohmann::json config = {
        {"generation", {
            {"reaction_type", "H(e,e'p)"},
            {"use_energy_loss", false},
            {"use_coulomb", false},
            {"beam_xwidth", 0.01},
            {"beam_ywidth", 0.01},
            {"electron", {
                {"delta_min", -10.0}, {"delta_max", 22.0},
                {"yptar_min", -0.030}, {"yptar_max", 0.030},
                {"xptar_min", -0.060}, {"xptar_max", 0.060},
                {"E_min", 500.0}, {"E_max", 2000.0}
            }},
            {"hadron", {
                {"delta_min", -10.0}, {"delta_max", 22.0},
                {"yptar_min", -0.030}, {"yptar_max", 0.030},
                {"xptar_min", -0.060}, {"xptar_max", 0.060},
                {"E_min", 500.0}, {"E_max", 2000.0}
            }}
        }},
        {"target", {
            {"mass", 938.272}, {"Z", 1}, {"A", 1},
            {"length", 10.0},
            {"x_offset", 0.0}, {"y_offset", 0.0}, {"z_offset", 0.0},
            {"angle", 0.0},
            {"raster_pattern", 0},
            {"raster_x", 0.0}, {"raster_y", 0.0}
        }},
        {"beam", {{"energy", 2200.0}, {"energy_spread", 0.0}}},
        {"spectrometer_electron", {{"momentum", 1000.0}, {"angle", 12.5}, {"phi", 0.0}}},
        {"spectrometer_hadron", {{"momentum", 1000.0}, {"angle", 50.0}, {"phi", 90.0}}}
    };
    
    ConfigManager cfg(config);
    auto random = std::make_shared<RandomGenerator>(12345);
    auto xs = std::make_shared<CrossSection>();
    auto optics = std::make_shared<SpectrometerOptics>();
    
    EventGenerator gen(cfg, random, xs, optics);
    
    // Test round-trip conversion
    double theta0 = 12.5 * kDegToRad;
    double phi0 = 90.0 * kDegToRad; // SOS position
    
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
    
    // Check round-trip (should recover original values)
    assert(approx_equal(xptar_orig, xptar_back, 1e-4));
    assert(approx_equal(yptar_orig, yptar_back, 1e-4));
    
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
        
        std::cout << "========================================" << std::endl;
        std::cout << "All tests PASSED!" << std::endl;
        std::cout << "========================================" << std::endl;
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Test FAILED with exception: " << e.what() << std::endl;
        return 1;
    }
}
