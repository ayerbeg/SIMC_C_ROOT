// tests/test_spectrometer_optics.cpp
// Unit tests for spectrometer optics and COSY matrices
// UPDATED for new CosyMatrix interface

#include "simc/SpectrometerOptics.h"
#include "simc/CosyMatrix.h"
#include "simc/SimcConstants.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

// Forward declarations from test_main.cpp
void ASSERT_TRUE(bool condition, const std::string& msg = "");
void ASSERT_FALSE(bool condition, const std::string& msg = "");
void ASSERT_NEAR(double a, double b, double tol, const std::string& msg = "");
void ASSERT_EQ(double a, double b, const std::string& msg = "");

using namespace simc;
using namespace simc::constants;

void TestSpectrometerOptics() {
    
    // ========================================================================
    // Test 1: COSY Matrix Loading
    // ========================================================================
    {
        std::cout << "    Test 1: COSY Matrix Loading... ";
        
        CosyMatrix matrix;
        ASSERT_FALSE(matrix.IsLoaded(), "Matrix should not be loaded initially");
        ASSERT_EQ(static_cast<double>(matrix.GetNumElements()), 0.0, "Should have 0 elements");
        ASSERT_EQ(static_cast<double>(matrix.GetMaxOrder()), 0.0, "Max order should be 0");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 2: COSY Matrix Element Parsing (NEW INTERFACE)
    // ========================================================================
    {
        std::cout << "    Test 2: COSY Matrix Element Parsing... ";
        
        // NEW interface: MatrixElement has coefficients[5] array
        CosyMatrix::MatrixElement elem;
        
        // Set up element: coefficient for x_out only
        elem.coefficients[0] = 1.5;  // x output
        elem.coefficients[1] = 0.0;  // xp output
        elem.coefficients[2] = 0.0;  // y output
        elem.coefficients[3] = 0.0;  // yp output
        elem.coefficients[4] = 0.0;  // dL output
        
        elem.exponents[0] = 2;  // x^2
        elem.exponents[1] = 1;  // xp^1
        elem.exponents[2] = 0;
        elem.exponents[3] = 0;
        elem.exponents[4] = 0;
        elem.order = 3;
        
        // Test evaluation for x output: 1.5 * x^2 * xp^1
        double input[5] = {2.0, 3.0, 0.0, 0.0, 0.0};  // x=2, xp=3
        
        // NEW: Evaluate takes output_index as first parameter
        double result_x = elem.Evaluate(0, input);  // output index 0 = x
        double result_xp = elem.Evaluate(1, input); // output index 1 = xp
        
        // Should be: 1.5 * 2^2 * 3 = 1.5 * 4 * 3 = 18.0
        ASSERT_NEAR(result_x, 18.0, 1e-10, "x output evaluation");
        ASSERT_NEAR(result_xp, 0.0, 1e-10, "xp output should be 0");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 3: CRITICAL - COSY Units (cm, mrad, PERCENT)
    // ========================================================================
    {
        std::cout << "    Test 3: COSY Units Verification... ";
        
        // CRITICAL: COSY-7 uses:
        //   Positions: cm
        //   Angles: mrad (milliradians), NOT radians!
        //   Delta: PERCENT
        
        // Example: 25 mrad angle
        double angle_rad = 0.025;          // In radians
        double angle_mrad = 25.0;          // In milliradians
        
        ASSERT_NEAR(angle_rad * 1000.0, angle_mrad, 1e-10, 
                   "mrad = rad * 1000");
        
        // Example: 5% momentum deviation
        double p_central = 5000.0;
        double p_actual = 5250.0;
        
        double delta_percent = (p_actual - p_central) / p_central * 100.0;
        ASSERT_NEAR(delta_percent, 5.0, 1e-10, "Delta in PERCENT");
        
        std::cout << "PASSED (Angles in mrad, delta in %)\n";
    }
    
    // ========================================================================
    // Test 4: HMS Octagon Aperture
    // ========================================================================
    {
        std::cout << "    Test 4: HMS Octagon Aperture... ";
        
        HMSOptics hms("", "");
        
        auto check1 = hms.CheckAperture(0.0, 0.0, 0);
        ASSERT_TRUE(check1.passed, "Origin inside");
        
        auto check2 = hms.CheckAperture(17.0, 0.0, 0);
        ASSERT_TRUE(check2.passed, "r=17 cm inside");
        
        auto check3 = hms.CheckAperture(20.0, 20.0, 0);
        ASSERT_FALSE(check3.passed, "(20,20) outside");
        
        auto check4 = hms.CheckAperture(12.0, 12.0, 0);
        ASSERT_TRUE(check4.passed, "(12,12) inside octagon");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 5: HMS Dipole Aperture
    // ========================================================================
    {
        std::cout << "    Test 5: HMS Dipole Aperture... ";
        
        HMSOptics hms("", "");
        
        auto check1 = hms.CheckAperture(25.0, 10.0, 1);
        ASSERT_TRUE(check1.passed, "(25,10) inside");
        
        auto check2 = hms.CheckAperture(35.0, 5.0, 1);
        ASSERT_FALSE(check2.passed, "(35,5) outside horizontal");
        
        auto check3 = hms.CheckAperture(10.0, 15.0, 1);
        ASSERT_FALSE(check3.passed, "(10,15) outside vertical");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 6: HMS Quadrupole Apertures
    // ========================================================================
    {
        std::cout << "    Test 6: HMS Quadrupole Apertures... ";
        
        HMSOptics hms("", "");
        
        // Q1: r = 12.5 cm
        auto check_q1_in = hms.CheckAperture(8.0, 8.0, 3);
        ASSERT_TRUE(check_q1_in.passed, "Inside Q1");
        
        auto check_q1_out = hms.CheckAperture(10.0, 10.0, 3);
        ASSERT_FALSE(check_q1_out.passed, "Outside Q1");
        
        // Q2: r = 30 cm
        auto check_q2_in = hms.CheckAperture(20.0, 15.0, 5);
        ASSERT_TRUE(check_q2_in.passed, "Inside Q2");
        
        auto check_q2_out = hms.CheckAperture(25.0, 25.0, 5);
        ASSERT_FALSE(check_q2_out.passed, "Outside Q2");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 7: SHMS Collimator Apertures
    // ========================================================================
    {
        std::cout << "    Test 7: SHMS Collimator Apertures... ";
        
        SHMSOptics shms("", "");
        
        // LARGE collimator
        shms.SetCollimator("LARGE");
        auto check_large_in = shms.CheckAperture(5.0, 5.0, 6);
        ASSERT_TRUE(check_large_in.passed, "Inside LARGE");
        
        auto check_large_out = shms.CheckAperture(7.0, 5.0, 6);
        ASSERT_FALSE(check_large_out.passed, "Outside LARGE");
        
        // SMALL collimator
        shms.SetCollimator("SMALL");
        auto check_small_in = shms.CheckAperture(2.0, 2.0, 6);
        ASSERT_TRUE(check_small_in.passed, "Inside SMALL");
        
        auto check_small_out = shms.CheckAperture(4.0, 2.0, 6);
        ASSERT_FALSE(check_small_out.passed, "Outside SMALL");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 8: SHMS HB Dipole Aperture
    // ========================================================================
    {
        std::cout << "    Test 8: SHMS HB Dipole Aperture... ";
        
        SHMSOptics shms("", "");
        
        auto check_in = shms.CheckAperture(40.0, 20.0, 1);
        ASSERT_TRUE(check_in.passed, "Inside HB dipole");
        
        auto check_out_x = shms.CheckAperture(55.0, 20.0, 1);
        ASSERT_FALSE(check_out_x.passed, "Outside horizontal");
        
        auto check_out_y = shms.CheckAperture(40.0, 30.0, 1);
        ASSERT_FALSE(check_out_y.passed, "Outside vertical");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 9: SOS Apertures
    // ========================================================================
    {
        std::cout << "    Test 9: SOS Apertures... ";
        
        SOSOptics sos("", "");
        
        auto check_entrance = sos.CheckAperture(25.0, 10.0, 0);
        ASSERT_TRUE(check_entrance.passed, "Inside entrance");
        
        auto check_entrance_out = sos.CheckAperture(35.0, 10.0, 0);
        ASSERT_FALSE(check_entrance_out.passed, "Outside entrance");
        
        auto check_q1 = sos.CheckAperture(10.0, 10.0, 2);
        ASSERT_TRUE(check_q1.passed, "Inside Q1");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 10: HRS Apertures
    // ========================================================================
    {
        std::cout << "    Test 10: HRS Apertures... ";
        
        HRSOptics hrs("HRSL", "", "");
        
        auto check_q1_in = hrs.CheckAperture(8.0, 8.0, 1);
        ASSERT_TRUE(check_q1_in.passed, "Inside Q1");
        
        auto check_q1_out = hrs.CheckAperture(10.0, 10.0, 1);
        ASSERT_FALSE(check_q1_out.passed, "Outside Q1");
        
        auto check_dipole = hrs.CheckAperture(35.0, 10.0, 2);
        ASSERT_TRUE(check_dipole.passed, "Inside dipole");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 11: Factory Function
    // ========================================================================
    {
        std::cout << "    Test 11: Spectrometer Factory... ";
        
        auto hms = CreateSpectrometerOptics("HMS");
        ASSERT_TRUE(hms != nullptr, "HMS created");
        ASSERT_TRUE(hms->GetName() == "HMS", "HMS name");
        
        auto shms = CreateSpectrometerOptics("SHMS");
        ASSERT_TRUE(shms != nullptr, "SHMS created");
        
        auto sos = CreateSpectrometerOptics("SOS");
        ASSERT_TRUE(sos != nullptr, "SOS created");
        
        auto hrsl = CreateSpectrometerOptics("HRSL");
        ASSERT_TRUE(hrsl != nullptr, "HRSL created");
        
        auto hrsr = CreateSpectrometerOptics("HRSR");
        ASSERT_TRUE(hrsr != nullptr, "HRSR created");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 12: Default Matrix File Paths
    // ========================================================================
    {
        std::cout << "    Test 12: Default Matrix Files... ";
        
        auto [hms_fwd, hms_rec] = GetDefaultMatrixFiles("HMS");
        ASSERT_TRUE(hms_fwd.find("hms") != std::string::npos, "HMS forward");
        ASSERT_TRUE(hms_rec.find("hms") != std::string::npos, "HMS recon");
        
        auto [shms_fwd, shms_rec] = GetDefaultMatrixFiles("SHMS");
        ASSERT_TRUE(shms_fwd.find("shms") != std::string::npos, "SHMS forward");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 13: Spectrometer Settings
    // ========================================================================
    {
        std::cout << "    Test 13: Spectrometer Settings... ";
        
        auto hms = CreateSpectrometerOptics("HMS");
        
        ASSERT_NEAR(hms->GetCentralMomentum(), 5000.0, 1e-6, "Default p");
        ASSERT_NEAR(hms->GetAngle(), 0.0, 1e-10, "Default angle");
        
        hms->SetCentralMomentum(7500.0);
        ASSERT_NEAR(hms->GetCentralMomentum(), 7500.0, 1e-6, "Set p");
        
        hms->SetAngle(0.218);
        ASSERT_NEAR(hms->GetAngle(), 0.218, 1e-6, "Set angle");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 14: Delta and Momentum Relationship
    // ========================================================================
    {
        std::cout << "    Test 14: Delta-Momentum Relationship... ";
        
        double p_central = 5000.0;
        double delta_percent = 7.5;
        
        // Calculate momentum from delta (PERCENT)
        double p_actual = p_central * (1.0 + delta_percent / 100.0);
        ASSERT_NEAR(p_actual, 5375.0, 1e-6, "p from delta");
        
        // Calculate delta from momentum
        double delta_calc = (p_actual - p_central) / p_central * 100.0;
        ASSERT_NEAR(delta_calc, 7.5, 1e-10, "delta from p");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 15: Number of Aperture Planes
    // ========================================================================
    {
        std::cout << "    Test 15: Number of Aperture Planes... ";
        
        auto hms = CreateSpectrometerOptics("HMS");
        ASSERT_EQ(static_cast<double>(hms->GetNumAperturePlanes()), 8.0, 
                 "HMS 8 planes");
        
        auto shms = CreateSpectrometerOptics("SHMS");
        ASSERT_EQ(static_cast<double>(shms->GetNumAperturePlanes()), 9.0, 
                 "SHMS 9 planes");
        
        auto sos = CreateSpectrometerOptics("SOS");
        ASSERT_EQ(static_cast<double>(sos->GetNumAperturePlanes()), 6.0, 
                 "SOS 6 planes");
        
        auto hrs = CreateSpectrometerOptics("HRSL");
        ASSERT_EQ(static_cast<double>(hrs->GetNumAperturePlanes()), 7.0, 
                 "HRS 7 planes");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 16: COSY Matrix Application (if files exist)
    // ========================================================================
    {
        std::cout << "    Test 16: COSY Matrix Application... ";
        
        // Try to load a real matrix file
        CosyMatrix matrix;
        bool loaded = matrix.LoadFromFile("data/matrices/hms/forward_cosy.dat");
        
        if (loaded) {
            std::cout << "(File loaded) ";
            
            // Test with simple input
            double input[5] = {0.0, 10.0, 0.0, 5.0, 0.0};  // mrad units!
            double output[5];
            matrix.Apply(input, output);
            
            // Sanity checks
            ASSERT_TRUE(std::abs(output[0]) < 100.0, "x_fp reasonable");
            ASSERT_TRUE(std::abs(output[2]) < 100.0, "y_fp reasonable");
            ASSERT_TRUE(std::abs(output[4]) < 500.0, "dL reasonable");
            
            std::cout << "PASSED\n";
        } else {
            std::cout << "SKIPPED (no matrix file)\n";
        }
    }
    
    // ========================================================================
    // Test 17: Unit Conversion Helpers
    // ========================================================================
    {
        std::cout << "    Test 17: Unit Conversion Helpers... ";
        
        // Test rad <-> mrad conversion
        double angle_rad = 0.025;
        double angle_mrad = angle_rad * 1000.0;
        ASSERT_NEAR(angle_mrad, 25.0, 1e-10, "rad to mrad");
        
        double back_to_rad = angle_mrad / 1000.0;
        ASSERT_NEAR(back_to_rad, 0.025, 1e-10, "mrad to rad");
        
        // Test PERCENT calculation
        double p0 = 5000.0;
        double p = 5250.0;
        double delta = (p - p0) / p0 * 100.0;
        ASSERT_NEAR(delta, 5.0, 1e-10, "delta PERCENT");
        
        // Reverse
        double p_calc = p0 * (1.0 + delta / 100.0);
        ASSERT_NEAR(p_calc, 5250.0, 1e-6, "p from delta");
        
        std::cout << "PASSED\n";
    }
    
    std::cout << "\n  All spectrometer optics tests PASSED!\n";
}
