// tests/test_main.cpp
// Simple standalone test framework (no Google Test dependency)

#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <sstream>

// ============================================================================
// Simple Test Framework
// ============================================================================

struct TestResult {
    std::string name;
    bool passed;
    std::string message;
};

class TestRunner {
public:
    void AddTest(const std::string& name, std::function<void()> test) {
        tests_.push_back({name, test});
    }
    
    int Run() {
        std::cout << "\n========================================\n";
        std::cout << "Running SIMC Tests\n";
        std::cout << "========================================\n\n";
        
        int passed = 0;
        int failed = 0;
        
        for (const auto& test : tests_) {
            std::cout << "Running: " << test.name << "... ";
            try {
                test.func();
                std::cout << "PASSED\n";
                ++passed;
            } catch (const std::exception& e) {
                std::cout << "FAILED\n";
                std::cout << "  Error: " << e.what() << "\n";
                ++failed;
            }
        }
        
        std::cout << "\n========================================\n";
        std::cout << "Results: " << passed << " passed, " << failed << " failed\n";
        std::cout << "========================================\n\n";
        
        return (failed == 0) ? 0 : 1;
    }
    
private:
    struct Test {
        std::string name;
        std::function<void()> func;
    };
    std::vector<Test> tests_;
};

// ============================================================================
// Assertion Helpers
// ============================================================================

void ASSERT_TRUE(bool condition, const std::string& msg = "") {
    if (!condition) {
        throw std::runtime_error("Assertion failed: " + msg);
    }
}

void ASSERT_FALSE(bool condition, const std::string& msg = "") {
    ASSERT_TRUE(!condition, msg);
}

void ASSERT_EQ(double a, double b, const std::string& msg = "") {
    if (std::abs(a - b) > 1e-10) {
        throw std::runtime_error("Values not equal: " + std::to_string(a) + 
                               " != " + std::to_string(b) + " " + msg);
    }
}

void ASSERT_NEAR(double a, double b, double tol, const std::string& msg = "") {
    if (std::abs(a - b) > tol) {
        throw std::runtime_error("Values not near: " + std::to_string(a) + 
                               " vs " + std::to_string(b) + " (tol=" + 
                               std::to_string(tol) + ") " + msg);
    }
}

// Template function - must be in header or defined here
template<typename T>
void ASSERT_EQ_T(T a, T b, const std::string& msg = "") {
    if (a != b) {
        std::ostringstream oss;
        oss << "Values not equal";
        if (!msg.empty()) {
            oss << ": " << msg;
        }
        oss << " (got " << a << ", expected " << b << ")";
        throw std::runtime_error(oss.str());
    }
}

// Explicit instantiations for types we use
template void ASSERT_EQ_T<int>(int, int, const std::string&);
template void ASSERT_EQ_T<unsigned int>(unsigned int, unsigned int, const std::string&);
template void ASSERT_EQ_T<unsigned long long>(unsigned long long, unsigned long long, const std::string&);
template void ASSERT_EQ_T<std::string>(std::string, std::string, const std::string&);

// ============================================================================
// Forward Declarations of Test Functions
// ============================================================================
void TestEvent();
void TestConfig();
void TestRandom();
void TestCrossSection();
void TestEnergyLoss();
void TestMultipleScattering();
void TestSpectrometerOptics();

// ============================================================================
// Main
// ============================================================================
int main() {
    TestRunner runner;
    
    // Register tests
    runner.AddTest("EventTests", TestEvent);
    runner.AddTest("ConfigTests", TestConfig);
    runner.AddTest("RandomTests", TestRandom);
    runner.AddTest("CrossSectionTests", TestCrossSection);
    runner.AddTest("EnergyLossTests", TestEnergyLoss);
    runner.AddTest("MultipleScatteringTests", TestMultipleScattering);
    runner.AddTest("SpectrometerOpticsTests", TestSpectrometerOptics);
    
    return runner.Run();
}
