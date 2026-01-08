// tests/test_random.cpp
// Unit tests for RandomGenerator (standalone version)

#include "simc/RandomGenerator.h"
#include <iostream>
#include <cmath>

// Forward declarations from test_main.cpp
void ASSERT_TRUE(bool condition, const std::string& msg = "");
void ASSERT_EQ(double a, double b, const std::string& msg = "");
void ASSERT_NEAR(double a, double b, double tol, const std::string& msg = "");
template<typename T> void ASSERT_EQ_T(T a, T b, const std::string& msg = "");

using namespace simc;

void TestRandom() {
    RandomGenerator rng(12345);  // Fixed seed for reproducibility
    
    // Test 1: Uniform range
    for (int i = 0; i < 1000; ++i) {
        double r = rng.Uniform();
        ASSERT_TRUE(r >= 0.0 && r < 1.0, "Uniform not in range");
    }
    
    // Test 2: Uniform min/max
    double min = 5.0, max = 10.0;
    for (int i = 0; i < 100; ++i) {
        double r = rng.Uniform(min, max);
        ASSERT_TRUE(r >= min && r < max, "Uniform(min,max) not in range");
    }
    
    // Test 3: Gaussian distribution
    double mean = 100.0;
    double sigma = 10.0;
    double sum = 0.0;
    double sum2 = 0.0;
    int n = 10000;
    
    for (int i = 0; i < n; ++i) {
        double r = rng.Gaussian(mean, sigma);
        sum += r;
        sum2 += r * r;
    }
    
    double calc_mean = sum / n;
    double calc_sigma = std::sqrt(sum2/n - calc_mean*calc_mean);
    
    // Should be close to expected (within 5%)
    ASSERT_NEAR(calc_mean, mean, mean * 0.05, "Gaussian mean off");
    ASSERT_NEAR(calc_sigma, sigma, sigma * 0.05, "Gaussian sigma off");
    
    // Test 4: UniformInt
    int imin = 1, imax = 10;
    for (int i = 0; i < 100; ++i) {
        int r = rng.UniformInt(imin, imax);
        ASSERT_TRUE(r >= imin && r <= imax, "UniformInt not in range");
    }
    
    // Test 5: Boolean
    int true_count = 0;
    n = 10000;
    for (int i = 0; i < n; ++i) {
        if (rng.Boolean(0.5)) ++true_count;
    }
    double fraction = static_cast<double>(true_count) / n;
    ASSERT_NEAR(fraction, 0.5, 0.05, "Boolean not ~50%");
    
    // Test 6: Reproducibility
    RandomGenerator rng1(12345);
    RandomGenerator rng2(12345);
    for (int i = 0; i < 100; ++i) {
        ASSERT_EQ(rng1.Uniform(), rng2.Uniform(), "Not reproducible");
    }
    
    // Test 7: Counter
    rng.ResetCount();
    ASSERT_EQ_T(rng.GetCount(), 0ull, "Counter not reset");
    for (int i = 0; i < 10; ++i) rng.Uniform();
    ASSERT_EQ_T(rng.GetCount(), 10ull, "Counter not incrementing");
    
    std::cout << "  All random tests passed\n";
}
