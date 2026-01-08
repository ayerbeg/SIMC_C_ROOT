// tests/test_config.cpp
// Unit tests for SimcEvent (standalone version)

#include "simc/ConfigManager.h"
#include <iostream>

// Forward declarations from test_main.cpp
void ASSERT_TRUE(bool condition, const std::string& msg);
void ASSERT_FALSE(bool condition, const std::string& msg);
void ASSERT_EQ(double a, double b, const std::string& msg);
template<typename T> void ASSERT_EQ_T(T a, T b, const std::string& msg);

using namespace simc;

void TestConfig() {
    // Create test configuration
    std::string test_json = R"({
        "beam": {
            "energy": 10.6,
            "energy_spread": 0.05
        },
        "target": {
            "type": "LH2",
            "thickness": 0.1,
            "A": 1,
            "Z": 1
        },
        "reaction": {
            "type": "elastic"
        },
        "monte_carlo": {
            "nevents": 10000,
            "random_seed": 12345
        },
        "output": {
            "filename": "test_output.root"
        }
    })";
    
    ConfigManager config;
    
    // Test 1: Load configuration
    ASSERT_TRUE(config.LoadFromString(test_json), "Failed to load config");
    ASSERT_TRUE(config.IsLoaded(), "Config not marked as loaded");
    
    // Test 2: Get basic values
    double energy = config.Get<double>("beam.energy");
    ASSERT_EQ(energy, 10.6, "beam.energy wrong");
    
    std::string target = config.Get<std::string>("target.type");
    ASSERT_EQ_T(target, std::string("LH2"), "target.type wrong");
    
    int nevents = config.Get<int>("monte_carlo.nevents");
    ASSERT_EQ_T(nevents, 10000, "nevents wrong");
    
    // Test 3: Get with default
    double existing = config.Get<double>("beam.energy", 5.0);
    ASSERT_EQ(existing, 10.6, "existing value with default wrong");
    
    double missing = config.Get<double>("beam.missing_param", 99.0);
    ASSERT_EQ(missing, 99.0, "default value not returned");
    
    // Test 4: Has method
    ASSERT_TRUE(config.Has("beam.energy"), "Has() failed for existing");
    ASSERT_FALSE(config.Has("beam.nonexistent"), "Has() true for nonexistent");
    
    // Test 5: Get structured data
    int n = config.GetNEvents();
    ASSERT_EQ_T(n, 10000, "GetNEvents wrong");
    
    unsigned int seed = config.GetRandomSeed();
    ASSERT_EQ_T(seed, 12345u, "GetRandomSeed wrong");
    
    std::cout << "  All config tests passed\n";
}

