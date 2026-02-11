// tests/test_energy_loss.cpp
// Unit tests for energy loss calculations

#include "simc/physics/EnergyLoss.h"
#include "simc/core/SimcConstants.h"
#include <iostream>
#include <cmath>

// Forward declarations from test_main.cpp
void ASSERT_TRUE(bool condition, const std::string& msg = "");
void ASSERT_NEAR(double a, double b, double tol, const std::string& msg = "");

using namespace simc;
using namespace simc::constants;

void TestEnergyLoss() {
    std::cout << "\n  Testing Energy Loss...\n";
    
    RandomGenerator rng(12345);  // Fixed seed
    
    // ========================================================================
    // Test 1: Ionization Potential
    // ========================================================================
    {
        double I_H = EnergyLoss::GetIonizationPotential(1.0);
        ASSERT_NEAR(I_H, 21.8e-06, 1e-8, "Hydrogen I should be 21.8 eV");
        
        double I_C = EnergyLoss::GetIonizationPotential(6.0);
        ASSERT_TRUE(I_C > I_H, "Carbon I should be larger than hydrogen");
        ASSERT_TRUE(I_C < 100e-06, "Carbon I should be reasonable");
        
        std::cout << "    Ionization potential test: PASSED\n";
    }
    
    // ========================================================================
    // Test 2: Plasma Frequency
    // ========================================================================
    {
        // Liquid hydrogen
        double hnup = EnergyLoss::PlasmaFrequency(0.0708, 1.0, 1.008);
        ASSERT_TRUE(hnup > 0.0, "Plasma frequency should be positive");
        ASSERT_TRUE(hnup < 1e-3, "Plasma frequency should be small (keV range)");
        
        std::cout << "    Plasma frequency test: PASSED\n";
    }
    
    // ========================================================================
    // Test 3: Density Correction
    // ========================================================================
    {
        // Test at low energy (should be ~0)
        double beta_low = 0.5;
        double gamma_low = 1.0 / std::sqrt(1.0 - beta_low*beta_low);
        double delta_low = EnergyLoss::DensityCorrection(beta_low, gamma_low, 1.0, 1.008, 0.0708);
        ASSERT_TRUE(delta_low >= 0.0, "Density correction should be non-negative");
        ASSERT_TRUE(delta_low < 1.0, "Density correction should be < 1 at low energy");
        
        // Test at high energy (should be larger)
        double beta_high = 0.99;
        double gamma_high = 1.0 / std::sqrt(1.0 - beta_high*beta_high);
        double delta_high = EnergyLoss::DensityCorrection(beta_high, gamma_high, 1.0, 1.008, 0.0708);
        ASSERT_TRUE(delta_high > delta_low, "Density correction should increase with energy");
        
        std::cout << "    Density correction test: PASSED\n";
    }
    
    // ========================================================================
    // Test 4: Xi Parameter
    // ========================================================================
    {
        double thickness = 0.708;  // 0.708 g/cm² (10 cm of LH2)
        double xi = EnergyLoss::CalcXi(thickness, 1.0, 1.008, 0.9);
        
        ASSERT_TRUE(xi > 0.0, "Xi should be positive");
        ASSERT_TRUE(xi < 10.0, "Xi should be reasonable");
        
        std::cout << "    Xi parameter test: PASSED\n";
    }
    
    // ========================================================================
    // Test 5: Most Probable Energy Loss
    // ========================================================================
    {
        double thickness = 0.708;  // g/cm²
        double beta = 0.9;
        double gamma = 1.0 / std::sqrt(1.0 - beta*beta);
        
        double Eloss_mp = EnergyLoss::MostProbableEnergyLoss(
            thickness, 1.0, 1.008, beta, gamma, 0.0708);
        
        ASSERT_TRUE(Eloss_mp >= 0.0, "Most probable loss should be non-negative");
        ASSERT_TRUE(Eloss_mp < 100.0, "Most probable loss should be reasonable (< 100 MeV)");
        
        std::cout << "    Most probable energy loss: " << Eloss_mp << " MeV\n";
        std::cout << "    Most probable loss test: PASSED\n";
    }
    
    // ========================================================================
    // Test 6: Full Energy Loss Calculation
    // ========================================================================
    {
        // 8 GeV electron through 10 cm of LH2
        double length = 10.0;      // cm
        double density = 0.0708;   // g/cm³
        double E_particle = 8000.0; // MeV
        double M_particle = Me;     // electron mass
        
        double Eloss = EnergyLoss::Calculate(
            length, density, 1.0, 1.008, E_particle, M_particle, rng,
            EnergyLoss::LossType::MOST_PROBABLE);
        
        ASSERT_TRUE(Eloss >= 0.0, "Energy loss should be non-negative");
        ASSERT_TRUE(Eloss < E_particle, "Energy loss should be less than particle energy");
        ASSERT_TRUE(Eloss < 100.0, "Energy loss through 10 cm LH2 should be < 100 MeV");
        
        std::cout << "    8 GeV electron through 10 cm LH2: " << Eloss << " MeV\n";
        std::cout << "    Energy loss calculation test: PASSED\n";
    }
    
    // ========================================================================
    // Test 7: Sampled Energy Loss
    // ========================================================================
    {
        // Test that sampled loss varies but stays reasonable
        double length = 10.0;
        double density = 0.0708;
        double E_particle = 8000.0;
        double M_particle = Me;
        
        double sum = 0.0;
        int N = 100;
        
        for (int i = 0; i < N; ++i) {
            double Eloss = EnergyLoss::Calculate(
                length, density, 1.0, 1.008, E_particle, M_particle, rng,
                EnergyLoss::LossType::SAMPLED);
            
            ASSERT_TRUE(Eloss >= 0.0 && Eloss < E_particle, 
                       "Sampled loss should be in physical range");
            sum += Eloss;
        }
        
        double mean_loss = sum / N;
        std::cout << "    Mean sampled loss (100 events): " << mean_loss << " MeV\n";
        std::cout << "    Sampled energy loss test: PASSED\n";
    }
    
    // ========================================================================
    // Test 8: Different Loss Types
    // Physics: Landau distribution with different lambda values
    // From enerloss_new.f:
    //   MINIMUM (lambda=3.0): Lower tail → HIGHER energy loss
    //   MOST_PROBABLE (lambda=1.0): Peak of distribution
    //   MAXIMUM (lambda=0.0067): Upper tail → LOWER energy loss
    // Therefore: MINIMUM > MOST_PROBABLE > MAXIMUM
    // ========================================================================
    {
        double length = 10.0;
        double density = 0.0708;
        double E_particle = 8000.0;
        double M_particle = Me;
        
        double loss_min = EnergyLoss::Calculate(
            length, density, 1.0, 1.008, E_particle, M_particle, rng,
            EnergyLoss::LossType::MINIMUM);
        
        double loss_mp = EnergyLoss::Calculate(
            length, density, 1.0, 1.008, E_particle, M_particle, rng,
            EnergyLoss::LossType::MOST_PROBABLE);
        
        double loss_max = EnergyLoss::Calculate(
            length, density, 1.0, 1.008, E_particle, M_particle, rng,
            EnergyLoss::LossType::MAXIMUM);
        
        // Correct physics from Fortran: MINIMUM > MOST_PROBABLE > MAXIMUM
        ASSERT_TRUE(loss_min > loss_mp, "Minimum ionizing (lambda=3) > most probable (lambda=1)");
        ASSERT_TRUE(loss_mp > loss_max, "Most probable (lambda=1) > maximum (lambda=0.0067)");
        
        std::cout << "    Minimum loss:       " << loss_min << " MeV (lambda=3.0)\n";
        std::cout << "    Most probable loss: " << loss_mp << " MeV (lambda=1.0)\n";
        std::cout << "    Maximum loss:       " << loss_max << " MeV (lambda=0.0067)\n";
        std::cout << "    Loss types test: PASSED\n";
    }
    
    // ========================================================================
    // Test 9: Material Database
    // ========================================================================
    {
        Material lh2 = Material::LiquidHydrogen();
        ASSERT_TRUE(lh2.Z == 1.0, "LH2 should have Z=1");
        ASSERT_NEAR(lh2.density, 0.0708, 0.001, "LH2 density");
        
        Material fe = Material::Iron();
        ASSERT_TRUE(fe.Z == 26.0, "Iron should have Z=26");
        ASSERT_TRUE(fe.density > lh2.density, "Iron denser than LH2");
        
        // Test GetByName
        Material lh2_by_name = Material::GetByName("LH2");
        ASSERT_TRUE(lh2_by_name.Z == 1.0, "GetByName should work");
        
        std::cout << "    Material database test: PASSED\n";
    }
    
    // ========================================================================
    // Test 10: Zero Thickness
    // ========================================================================
    {
        double Eloss = EnergyLoss::Calculate(
            0.0, 0.0708, 1.0, 1.008, 8000.0, Me, rng);
        
        ASSERT_TRUE(Eloss == 0.0, "Zero thickness should give zero loss");
        
        std::cout << "    Zero thickness test: PASSED\n";
    }
    
    // ========================================================================
    // Test 11: Energy Dependence
    // ========================================================================
    {
        double length = 10.0;
        double density = 0.0708;
        
        // Test that higher energy particles lose less (relatively)
        double E1 = 1000.0;  // 1 GeV
        double E2 = 10000.0; // 10 GeV
        
        double loss1 = EnergyLoss::Calculate(
            length, density, 1.0, 1.008, E1, Me, rng,
            EnergyLoss::LossType::MOST_PROBABLE);
        
        double loss2 = EnergyLoss::Calculate(
            length, density, 1.0, 1.008, E2, Me, rng,
            EnergyLoss::LossType::MOST_PROBABLE);
        
        // dE/dx decreases with energy in the relativistic region
        // but absolute loss might be similar - just check it's reasonable
        ASSERT_TRUE(loss1 >= 0.0 && loss2 >= 0.0, "Both losses should be non-negative");
        ASSERT_TRUE(std::abs(loss1 - loss2) < 50.0, "Losses should be similar in relativistic region");
        
        std::cout << "    Energy dependence test: PASSED\n";
    }
    
    std::cout << "\n  All energy loss tests PASSED!\n";
}
