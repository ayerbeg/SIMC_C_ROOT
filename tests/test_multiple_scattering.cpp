// tests/test_multiple_scattering.cpp
// Unit tests for multiple scattering calculations

#include "simc/physics/MultipleScattering.h"
#include "simc/core/SimcConstants.h"
#include <iostream>
#include <cmath>

// Forward declarations from test_main.cpp
void ASSERT_TRUE(bool condition, const std::string& msg = "");
void ASSERT_NEAR(double a, double b, double tol, const std::string& msg = "");

using namespace simc;
using namespace simc::constants;

void TestMultipleScattering() {
    std::cout << "\n  Testing Multiple Scattering...\n";
    
    RandomGenerator rng(12345);  // Fixed seed
    
    // ========================================================================
    // Test 1: Beta Calculation
    // ========================================================================
    {
        // Test for electron at 8 GeV
        double p = 8000.0;  // MeV/c
        double m = Me;      // electron mass
        double beta = MultipleScattering::CalculateBeta(p, m);
        
        // At 8 GeV, electron is highly relativistic, beta ~ 1
        ASSERT_TRUE(beta > 0.999, "High energy electron should have beta ~1");
        ASSERT_TRUE(beta < 1.0, "Beta should be less than 1");
        
        // Test for proton at 1 GeV
        p = 1000.0;  // MeV/c
        m = Mp;      // proton mass
        beta = MultipleScattering::CalculateBeta(p, m);
        
        // At 1 GeV, proton beta ~ 0.729
        ASSERT_NEAR(beta, 0.729, 0.05, "1 GeV proton beta");
        
        std::cout << "    Beta calculation test: PASSED\n";
    }
    
    // ========================================================================
    // Test 2: Beta*Gamma Calculation
    // ========================================================================
    {
        double p = 8000.0;
        double m = Me;
        double bg = MultipleScattering::CalculateBetaGamma(p, m);
        
        // For electron at 8 GeV: β*γ = p/m ≈ 8000/0.511 ≈ 15655
        ASSERT_NEAR(bg, 15655.0, 500.0, "Beta*gamma for 8 GeV electron");
        
        std::cout << "    Beta*gamma calculation test: PASSED\n";
    }
    
    // ========================================================================
    // Test 3: Highland RMS Formula
    // ========================================================================
    {
        // Test Highland formula for known case
        // 8 GeV electron through 10 cm of LH2
        
        Material lh2 = Material::LiquidHydrogen();
        double X0_cm = MultipleScattering::GetRadiationLength(lh2);
        double thickness = 10.0 / X0_cm;  // thickness in radiation lengths
        
        double p = 8000.0;  // MeV/c
        double charge = 1.0;
        double mass = Me;
        
        double theta_rms = MultipleScattering::CalculateRMS(thickness, p, charge, mass);
        
        // Should be positive and small (< 1 mrad for thin target)
        ASSERT_TRUE(theta_rms > 0.0, "RMS angle should be positive");
        ASSERT_TRUE(theta_rms < 0.01, "RMS angle should be small (< 10 mrad)");
        
        std::cout << "    Highland RMS for 8 GeV e- through 10 cm LH2: " 
                  << theta_rms * 1000.0 << " mrad\n";
        std::cout << "    Highland formula test: PASSED\n";
    }
    
    // ========================================================================
    // Test 4: RMS Scaling with Thickness
    // ========================================================================
    {
        // Multiple scattering should scale as sqrt(thickness)
        Material lh2 = Material::LiquidHydrogen();
        double X0_cm = MultipleScattering::GetRadiationLength(lh2);
        
        double p = 5000.0;
        double charge = 1.0;
        double mass = Me;
        
        double thickness1 = 1.0 / X0_cm;   // 1 cm
        double thickness2 = 4.0 / X0_cm;   // 4 cm
        
        double rms1 = MultipleScattering::CalculateRMS(thickness1, p, charge, mass);
        double rms2 = MultipleScattering::CalculateRMS(thickness2, p, charge, mass);
        
        // rms2 should be approximately sqrt(4) = 2 times rms1
        double ratio = rms2 / rms1;
        ASSERT_NEAR(ratio, 2.0, 0.2, "RMS should scale as sqrt(thickness)");
        
        std::cout << "    Thickness scaling test: PASSED\n";
    }
    
    // ========================================================================
    // Test 5: RMS Scaling with Momentum
    // ========================================================================
    {
        // At high energies, theta_rms ~ 1/p
        Material lh2 = Material::LiquidHydrogen();
        double X0_cm = MultipleScattering::GetRadiationLength(lh2);
        double thickness = 10.0 / X0_cm;
        
        double p1 = 4000.0;  // 4 GeV
        double p2 = 8000.0;  // 8 GeV
        double charge = 1.0;
        double mass = Me;
        
        double rms1 = MultipleScattering::CalculateRMS(thickness, p1, charge, mass);
        double rms2 = MultipleScattering::CalculateRMS(thickness, p2, charge, mass);
        
        // rms1 should be approximately 2 times rms2 (inverse momentum)
        double ratio = rms1 / rms2;
        ASSERT_NEAR(ratio, 2.0, 0.2, "RMS should scale as 1/p");
        
        std::cout << "    Momentum scaling test: PASSED\n";
    }
    
    // ========================================================================
    // Test 6: Gaussian Sampling
    // ========================================================================
    {
        // Test that Gaussian sampling produces correct distribution
        double theta_rms = 0.001;  // 1 mrad
        int N = 10000;
        
        double sum_x = 0.0;
        double sum_x2 = 0.0;
        double sum_y = 0.0;
        double sum_y2 = 0.0;
        
        for (int i = 0; i < N; ++i) {
            auto angles = MultipleScattering::SampleGaussian(theta_rms, rng);
            sum_x += angles.theta_x;
            sum_x2 += angles.theta_x * angles.theta_x;
            sum_y += angles.theta_y;
            sum_y2 += angles.theta_y * angles.theta_y;
        }
        
        double mean_x = sum_x / N;
        double mean_y = sum_y / N;
        double rms_x = std::sqrt(sum_x2 / N - mean_x * mean_x);
        double rms_y = std::sqrt(sum_y2 / N - mean_y * mean_y);
        
        // Mean should be ~0
        ASSERT_NEAR(mean_x, 0.0, 0.0001, "Gaussian mean_x should be 0");
        ASSERT_NEAR(mean_y, 0.0, 0.0001, "Gaussian mean_y should be 0");
        
        // RMS should match input
        ASSERT_NEAR(rms_x, theta_rms, 0.0001, "Gaussian RMS_x should match");
        ASSERT_NEAR(rms_y, theta_rms, 0.0001, "Gaussian RMS_y should match");
        
        std::cout << "    Gaussian sampling test: PASSED\n";
    }
    
    // ========================================================================
    // Test 7: Full Scattering Calculation (Gaussian)
    // ========================================================================
    {
        Material lh2 = Material::LiquidHydrogen();
        double length = 10.0;  // cm
        double p = 8000.0;     // MeV/c
        double charge = 1.0;
        double mass = Me;
        
        auto angles = MultipleScattering::Calculate(
            length, lh2, p, charge, mass, rng, false);  // Gaussian mode
        
        ASSERT_TRUE(angles.theta_rms > 0.0, "RMS should be positive");
        ASSERT_TRUE(std::abs(angles.theta_x) < 0.1, "theta_x should be reasonable");
        ASSERT_TRUE(std::abs(angles.theta_y) < 0.1, "theta_y should be reasonable");
        
        std::cout << "    Full calculation (Gaussian): PASSED\n";
    }
    
    // ========================================================================
    // Test 8: Full Scattering Calculation (Molière)
    // ========================================================================
    {
        Material lh2 = Material::LiquidHydrogen();
        double length = 1.0;   // cm (thin absorber for Molière)
        double p = 8000.0;     // MeV/c
        double charge = 1.0;
        double mass = Me;
        
        auto angles = MultipleScattering::Calculate(
            length, lh2, p, charge, mass, rng, true);  // Molière mode
        
        ASSERT_TRUE(angles.theta_rms > 0.0, "RMS should be positive");
        ASSERT_TRUE(std::abs(angles.theta_x) < 0.1, "theta_x should be reasonable");
        ASSERT_TRUE(std::abs(angles.theta_y) < 0.1, "theta_y should be reasonable");
        
        std::cout << "    Full calculation (Molière): PASSED\n";
    }
    
    // ========================================================================
    // Test 9: Chi Alpha Calculation
    // ========================================================================
    {
        double Z = 1.0;      // Hydrogen
        double beta = 0.99;  // Relativistic
        
        double chi_alpha = MultipleScattering::CalculateChiAlpha(Z, beta);
        
        ASSERT_TRUE(chi_alpha > 0.0, "Chi_alpha should be positive");
        ASSERT_TRUE(chi_alpha < 1.0, "Chi_alpha should be small");
        
        std::cout << "    Chi_alpha calculation test: PASSED\n";
    }
    
    // ========================================================================
    // Test 10: Molière B Parameter
    // ========================================================================
    {
        double chi_alpha = 1e-5;
        double Z = 1.0;
        
        double B = MultipleScattering::CalculateMoliereB(chi_alpha, Z);
        
        ASSERT_TRUE(B > 0.0, "Molière B should be positive");
        ASSERT_TRUE(B < 20.0, "Molière B should be reasonable");
        
        std::cout << "    Molière B parameter test: PASSED\n";
    }
    
    // ========================================================================
    // Test 11: Different Materials
    // ========================================================================
    {
        Material lh2 = Material::LiquidHydrogen();
        Material fe = Material::Iron();
        
        double length = 1.0;  // 1 cm
        double p = 5000.0;
        double charge = 1.0;
        double mass = Me;
        
        auto angles_lh2 = MultipleScattering::Calculate(
            length, lh2, p, charge, mass, rng, false);
        
        auto angles_fe = MultipleScattering::Calculate(
            length, fe, p, charge, mass, rng, false);
        
        // Iron should scatter more than LH2 (shorter radiation length)
        ASSERT_TRUE(angles_fe.theta_rms > angles_lh2.theta_rms,
                   "Iron should scatter more than LH2");
        
        std::cout << "    Material comparison test: PASSED\n";
    }
    
    // ========================================================================
    // Test 12: Zero Length Edge Case
    // ========================================================================
    {
        Material lh2 = Material::LiquidHydrogen();
        double length = 0.0;
        double p = 5000.0;
        double charge = 1.0;
        double mass = Me;
        
        auto angles = MultipleScattering::Calculate(
            length, lh2, p, charge, mass, rng, false);
        
        ASSERT_TRUE(angles.theta_rms == 0.0, "Zero length should give zero scattering");
        ASSERT_TRUE(angles.theta_x == 0.0, "Zero length should give zero theta_x");
        ASSERT_TRUE(angles.theta_y == 0.0, "Zero length should give zero theta_y");
        
        std::cout << "    Zero length test: PASSED\n";
    }
    
    // ========================================================================
    // Test 13: Charge Dependence
    // ========================================================================
    {
        Material lh2 = Material::LiquidHydrogen();
        double length = 10.0;
        double p = 5000.0;
        double mass = Me;
        
        double X0_cm = MultipleScattering::GetRadiationLength(lh2);
        double thickness = length / X0_cm;
        
        double rms_z1 = MultipleScattering::CalculateRMS(thickness, p, 1.0, mass);
        double rms_z2 = MultipleScattering::CalculateRMS(thickness, p, 2.0, mass);
        
        // RMS should scale linearly with charge
        double ratio = rms_z2 / rms_z1;
        ASSERT_NEAR(ratio, 2.0, 0.1, "RMS should scale linearly with charge");
        
        std::cout << "    Charge dependence test: PASSED\n";
    }
    
    // ========================================================================
    // Test 14: Multiple Scattering Table
    // ========================================================================
    {
        Material lh2 = Material::LiquidHydrogen();
        
        MultipleScatteringTable table(lh2, 1000.0, 10000.0, 50);
        
        double p = 5000.0;
        double length = 10.0;
        double charge = 1.0;
        double mass = Me;  // Use electron mass for comparison
        
        double rms_table = table.GetRMS(p, length, charge, mass);
        
        // Should give reasonable value
        ASSERT_TRUE(rms_table > 0.0, "Table RMS should be positive");
        ASSERT_TRUE(rms_table < 0.1, "Table RMS should be reasonable");
        
        // Compare with direct calculation
        double X0_cm = MultipleScattering::GetRadiationLength(lh2);
        double thickness = length / X0_cm;
        double rms_direct = MultipleScattering::CalculateRMS(thickness, p, charge, mass);
        
        // Should be very similar for electrons (within 10%)
        // The table is built for electrons, so should match closely
        double rel_diff = std::abs(rms_table - rms_direct) / rms_direct;
        
        std::cout << "      Table RMS:  " << rms_table * 1000.0 << " mrad\n";
        std::cout << "      Direct RMS: " << rms_direct * 1000.0 << " mrad\n";
        std::cout << "      Rel diff:   " << rel_diff * 100.0 << " %\n";
        
        ASSERT_TRUE(rel_diff < 0.15, "Table should match direct calculation within 15%");
        
        std::cout << "    Multiple scattering table test: PASSED\n";
    }
    
    std::cout << "\n  All multiple scattering tests PASSED!\n";
}
