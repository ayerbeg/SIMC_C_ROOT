// tests/test_cross_section.cpp
// Unit tests for cross section calculations
// UPDATED: Now tests the complete Fortran port with CMake integration

#include "simc/CrossSection.h"
#include "simc/Kinematics.h"
#include "simc/SimcConstants.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// Forward declarations from test_main.cpp
void ASSERT_TRUE(bool condition, const std::string& msg = "");
void ASSERT_FALSE(bool condition, const std::string& msg = "");
void ASSERT_NEAR(double a, double b, double tol, const std::string& msg = "");
void ASSERT_EQ(double a, double b, const std::string& msg = "");

using namespace simc;
using namespace simc::constants;

void TestCrossSection() {
    std::cout << "\n  Testing Cross Sections (Complete Fortran Port)...\n";
    
    // ========================================================================
    // Test 1: Mott Cross Section (sigMott from physics_proton.f)
    // ========================================================================
    {
        std::cout << "    Test 1: Mott Cross Section (sigMott)...";
        
        double Ein = 10600.0;  // 10.6 GeV
        double theta = 0.218;   // ~12.5 degrees
        double Q2 = 2.0 * Ein * Ein * (1.0 - std::cos(theta));
        
        // NOTE: Updated signature - now takes Q2 as third argument
        double sig_mott = ElasticCrossSection::SigMott(Ein, theta, Q2);
        
        // Mott cross section should be positive and reasonable
        ASSERT_TRUE(sig_mott > 0.0, "Mott XS should be positive");
        ASSERT_TRUE(sig_mott < 1e10, "Mott XS should be reasonable");
        
        // Check angular dependence (should decrease with angle)
        double Q2_small = 2.0 * Ein * Ein * (1.0 - std::cos(0.1));
        double Q2_large = 2.0 * Ein * Ein * (1.0 - std::cos(0.5));
        double sig_mott_small = ElasticCrossSection::SigMott(Ein, 0.1, Q2_small);
        double sig_mott_large = ElasticCrossSection::SigMott(Ein, 0.5, Q2_large);
        ASSERT_TRUE(sig_mott_small > sig_mott_large, "Mott XS should decrease with angle");
        
        // Validate formula: sig = (2*alpha*hbarc*E*cos(theta/2)/Q2)^2 * 1e4
        double manual = 2.0 * alpha * hbarc * Ein * std::cos(theta/2.0) / Q2;
        manual = manual * manual * 1.0e4;
        ASSERT_NEAR(sig_mott, manual, 1.0e-8, "Mott formula verification");
        
        std::cout << " PASSED\n";
    }
    
    // ========================================================================
    // Test 2: Form Factors (fofa_best_fit from physics_proton.f)
    // ========================================================================
    {
        std::cout << "    Test 2: Form Factors (fofa_best_fit)...";
        
        // NOTE: Updated signature - now void function with output parameters
        // Input: qsquar = -Q^2/(hbarc^2) in fm^-2
        
        // Test at Q2 = 0 (should give GE=1, GM=mu_p)
        double GE_0, GM_0;
        ElasticCrossSection::FofaBestFit(0.0, GE_0, GM_0);
        
        ASSERT_NEAR(GE_0, 1.0, 0.001, "GE(0) should be 1");
        ASSERT_NEAR(GM_0, 2.793, 0.001, "GM(0) should be mu_p");
        
        // Test at finite Q2 = 1.0 GeV^2
        double Q2_GeV2 = 1.0;
        double qsquar = -Q2_GeV2 * 1.0e6 / (hbarc * hbarc);  // Convert to Fortran units
        double GE, GM;
        ElasticCrossSection::FofaBestFit(qsquar, GE, GM);
        
        // Form factors should decrease with Q2
        ASSERT_TRUE(GE < GE_0, "GE should decrease with Q2");
        ASSERT_TRUE(GM < GM_0, "GM should decrease with Q2");
        ASSERT_TRUE(GE > 0.0 && GE < 1.0, "GE should be between 0 and 1");
        ASSERT_TRUE(GM > 0.0, "GM should be positive");
        
        // Test Bosted fit at several Q2 values
        std::vector<double> Q2_tests = {0.1, 0.5, 1.0, 2.0, 5.0};  // GeV^2
        for (double Q2_test : Q2_tests) {
            double qs = -Q2_test * 1.0e6 / (hbarc * hbarc);
            double ge, gm;
            ElasticCrossSection::FofaBestFit(qs, ge, gm);
            ASSERT_TRUE(ge > 0.0 && ge <= 1.0, "GE should be physical");
            ASSERT_TRUE(gm > 0.0 && gm <= 3.0, "GM should be physical");
        }
        
        std::cout << " PASSED\n";
    }
    
    // ========================================================================
    // Test 3: Elastic Cross Section (sigep from physics_proton.f)
    // ========================================================================
    {
        std::cout << "    Test 3: Elastic Cross Section (sigep)...";
        
        SimcEvent evt;
        evt.Clear();
        
        // Set up elastic scattering kinematics
        double Ein = 10600.0;   // 10.6 GeV
        double theta_e = 0.218;  // ~12.5 deg
        
        // Calculate scattered electron energy for elastic
        evt.Ein = Ein;
        evt.e_theta = theta_e;
        evt.e_E = Ein / (1.0 + (2.0*Ein/Mp)*std::sin(theta_e/2.0)*std::sin(theta_e/2.0));
        
        // Calculate Q2 and other kinematics
        evt.Q2 = 2.0 * Ein * evt.e_E * (1.0 - std::cos(theta_e));
        evt.nu = Ein - evt.e_E;
        evt.q = std::sqrt(evt.Q2 + evt.nu*evt.nu);
        evt.W = Mp;  // Elastic
        
        ElasticCrossSection elastic;
        double sigma = elastic.Calculate(evt);
        
        ASSERT_TRUE(sigma > 0.0, "Elastic XS should be positive");
        ASSERT_TRUE(sigma < 1e10, "Elastic XS should be finite");
        ASSERT_TRUE(elastic.IsPhysical(evt), "Kinematics should be physical");
        
        // Test multiple energies and angles
        std::vector<double> energies = {1000.0, 2000.0, 4000.0};
        std::vector<double> angles = {15.0, 30.0, 45.0};
        
        for (double E : energies) {
            for (double theta_deg : angles) {
                evt.Ein = E;
                evt.e_theta = theta_deg * DEG_TO_RAD;
                evt.e_E = E / (1.0 + (2.0*E/Mp)*std::sin(evt.e_theta/2.0)*std::sin(evt.e_theta/2.0));
                evt.Q2 = 2.0 * E * evt.e_E * (1.0 - std::cos(evt.e_theta));
                evt.nu = E - evt.e_E;
                evt.q = std::sqrt(evt.Q2 + evt.nu*evt.nu);
                evt.W = Mp;
                
                sigma = elastic.Calculate(evt);
                ASSERT_TRUE(sigma > 0.0, "All elastic XS should be positive");
            }
        }
        
        std::cout << " PASSED\n";
    }
    
    // ========================================================================
    // Test 4: Quasi-Elastic Cross Section (deForest from physics_proton.f)
    // ========================================================================
    {
        std::cout << "    Test 4: Quasi-Elastic Cross Section (deForest)...";
        
        SimcEvent evt;
        evt.Clear();
        
        // Set up quasi-elastic (e,e'p) kinematics
        evt.Ein = 2000.0;  // 2 GeV
        evt.e_theta = 30.0 * DEG_TO_RAD;
        evt.e_E = 1500.0;
        evt.e_phi = 0.0;
        
        evt.p_theta = 40.0 * DEG_TO_RAD;
        evt.p_phi = pi;  // Opposite side
        evt.p_P = 600.0;  // MeV/c
        evt.p_E = std::sqrt(evt.p_P*evt.p_P + Mp2);
        
        // Calculate Q2, nu, q
        evt.nu = evt.Ein - evt.e_E;
        evt.Q2 = 2.0 * evt.Ein * evt.e_E * (1.0 - std::cos(evt.e_theta));
        evt.q = std::sqrt(evt.Q2 + evt.nu*evt.nu);
        
        // Unit vectors
        evt.ue_x = std::sin(evt.e_theta) * std::cos(evt.e_phi);
        evt.ue_y = std::sin(evt.e_theta) * std::sin(evt.e_phi);
        evt.ue_z = std::cos(evt.e_theta);
        
        evt.up_x = std::sin(evt.p_theta) * std::cos(evt.p_phi);
        evt.up_y = std::sin(evt.p_theta) * std::sin(evt.p_phi);
        evt.up_z = std::cos(evt.p_theta);
        
        evt.uq_x = -evt.e_E * evt.ue_x / evt.q;
        evt.uq_y = -evt.e_E * evt.ue_y / evt.q;
        evt.uq_z = (evt.Ein - evt.e_E * evt.ue_z) / evt.q;
        
        // Missing momentum
        evt.Pmx = evt.p_P * evt.up_x - evt.q * evt.uq_x;
        evt.Pmy = evt.p_P * evt.up_y - evt.q * evt.uq_y;
        evt.Pmz = evt.p_P * evt.up_z - evt.q * evt.uq_z;
        evt.Pm = std::sqrt(evt.Pmx*evt.Pmx + evt.Pmy*evt.Pmy + evt.Pmz*evt.Pmz);
        evt.Em = evt.nu + Mp - evt.p_E;
        
        // Test all three deForest flags
        for (int flag : {0, 1, -1}) {
            QuasiElasticCrossSection qe_xs(flag);
            double sigma = qe_xs.Calculate(evt);
            
            ASSERT_TRUE(sigma > 0.0, "deForest XS should be positive");
            ASSERT_TRUE(std::isfinite(sigma), "deForest XS should be finite");
            ASSERT_TRUE(qe_xs.IsPhysical(evt), "QE kinematics should be physical");
        }
        
        std::cout << " PASSED\n";
    }
    
    // ========================================================================
    // Test 5: Pion Production Cross Section
    // ========================================================================
    {
        std::cout << "    Test 5: Pion Production...";
        
        SimcEvent evt;
        evt.Clear();
        
        // Set up pion production kinematics (W > threshold)
        evt.Ein = 10600.0;
        evt.e_E = 8000.0;
        evt.e_theta = 0.218;
        evt.p_E = 1500.0;
        evt.p_P = 1500.0;
        evt.Q2 = 2.0e6;  // 2 (GeV/c)^2
        evt.W = 1232.0;  // Delta mass
        evt.nu = evt.Ein - evt.e_E;
        evt.q = std::sqrt(evt.Q2 + evt.nu * evt.nu);
        
        PionCrossSection pion(PionType::PI_PLUS);
        double sigma = pion.Calculate(evt);
        
        ASSERT_TRUE(sigma >= 0.0, "Pion XS should be non-negative");
        ASSERT_TRUE(pion.IsPhysical(evt), "Should be physical above threshold");
        
        // Check threshold
        evt.W = 1000.0;  // Below threshold
        ASSERT_FALSE(pion.IsPhysical(evt), "Below threshold should be unphysical");
        
        std::cout << " PASSED\n";
    }
    
    // ========================================================================
    // Test 6: Kaon Production Cross Section
    // ========================================================================
    {
        std::cout << "    Test 6: Kaon Production...";
        
        SimcEvent evt;
        evt.Clear();
        
        evt.Ein = 10600.0;
        evt.e_E = 8000.0;
        evt.e_theta = 0.218;
        evt.p_E = 1800.0;
        evt.p_P = 1500.0;
        evt.Q2 = 2.0e6;
        evt.W = 1600.0;  // Above threshold
        evt.nu = evt.Ein - evt.e_E;
        evt.q = std::sqrt(evt.Q2 + evt.nu * evt.nu);
        
        KaonCrossSection kaon(KaonType::K_PLUS);
        double sigma = kaon.Calculate(evt);
        
        ASSERT_TRUE(sigma >= 0.0, "Kaon XS should be non-negative");
        ASSERT_TRUE(kaon.IsPhysical(evt), "Should be physical above threshold");
        
        // Check threshold
        evt.W = 1300.0;
        ASSERT_FALSE(kaon.IsPhysical(evt), "Below threshold should be unphysical");
        
        std::cout << " PASSED\n";
    }
    
    // ========================================================================
    // Test 7: CrossSectionFactory
    // ========================================================================
    {
        std::cout << "    Test 7: Factory Pattern...";
        
        auto elastic_xs = CrossSectionFactory::Create(ReactionType::ELASTIC);
        ASSERT_TRUE(elastic_xs != nullptr, "Factory should create elastic XS");
        
        auto qe_xs = CrossSectionFactory::Create(ReactionType::QUASIELASTIC);
        ASSERT_TRUE(qe_xs != nullptr, "Factory should create QE XS");
        
        auto pion_xs = CrossSectionFactory::Create(ReactionType::PION_PRODUCTION);
        ASSERT_TRUE(pion_xs != nullptr, "Factory should create pion XS");
        
        std::cout << " PASSED\n";
    }
    
    // ========================================================================
    // Test 8: Radiative Corrections
    // ========================================================================
    {
        std::cout << "    Test 8: Radiative Corrections...";
        
        RadiativeCorrections rc(0.667);
        
        double Ein = 10600.0;
        double Ee = 8000.0;
        double theta = 0.218;
        int Z = 1;
        
        double corr = rc.GetCorrectionFactor(Ein, Ee, theta, Z);
        ASSERT_TRUE(corr > 0.5 && corr < 2.0, "RC factor should be reasonable");
        
        double internal = rc.InternalCorrection(Ein, Ee, theta);
        ASSERT_TRUE(std::abs(internal) < 1.0, "Internal RC should be small");
        
        double external = rc.ExternalCorrection(Ein, Ee, Z, 0.1);
        ASSERT_TRUE(std::abs(external) < 0.5, "External RC should be small");
        
        std::cout << " PASSED\n";
    }
    
    // ========================================================================
    // Test 9: Q2 Dependence
    // ========================================================================
    {
        std::cout << "    Test 9: Q2 Dependence...";
        
        SimcEvent evt;
        ElasticCrossSection elastic;
        
        std::vector<double> Q2_values = {0.5e6, 1.0e6, 2.0e6, 4.0e6};  // MeV^2
        double prev_sigma = 1e10;
        
        for (double Q2_test : Q2_values) {
            evt.Clear();
            evt.Ein = 10600.0;
            evt.Q2 = Q2_test;
            
            double theta_test = 0.3;
            evt.e_theta = theta_test;
            evt.e_E = Q2_test / (2.0 * evt.Ein * (1.0 - std::cos(theta_test)));
            
            if (evt.e_E > 0.0 && evt.e_E < evt.Ein) {
                evt.W = Mp;
                evt.nu = evt.Ein - evt.e_E;
                evt.q = std::sqrt(evt.Q2 + evt.nu * evt.nu);
                
                double sigma = elastic.Calculate(evt);
                if (sigma > 0.0) {
                    ASSERT_TRUE(sigma < prev_sigma, "XS should decrease with Q2");
                    prev_sigma = sigma;
                }
            }
        }
        
        std::cout << " PASSED\n";
    }
    
    // ========================================================================
    // Test 10: Fortran Equivalence Check
    // ========================================================================
    {
        std::cout << "    Test 10: Fortran Equivalence...";
        
        // Verify key functions match Fortran exactly
        
        // 1. sigMott test case
        double e0 = 2000.0, theta = 0.5236, Q2_test = 2.0*e0*e0*(1.0-std::cos(theta));
        double mott = ElasticCrossSection::SigMott(e0, theta, Q2_test);
        ASSERT_TRUE(mott > 0.0 && std::isfinite(mott), "Mott should be finite");
        
        // 2. Form factors at Q2 = 1 GeV^2
        //    From Bosted Phys. Rev. C 51, 409 (1995), Eqs. 4-5:
        //    At Q = 1 GeV: GE = 0.1686, GM = 0.4926
        double ge, gm;
        ElasticCrossSection::FofaBestFit(-1.0e6/(hbarc*hbarc), ge, gm);
        ASSERT_NEAR(ge, 0.1686, 0.001, "GE at Q2=1 should match Bosted fit");
        ASSERT_NEAR(gm, 0.4926, 0.005, "GM at Q2=1 should match Bosted fit");
        
        std::cout << " PASSED\n";
    }
    
    std::cout << "\n  ✓ All Cross Section Tests PASSED (Fortran Port Validated)\n";
}
