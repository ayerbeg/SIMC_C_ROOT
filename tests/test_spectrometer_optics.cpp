// tests/test_spectrometer_optics.cpp
// Unit tests for spectrometer optics and COSY matrices

#include "simc/SpectrometerOptics.h"
#include "simc/CosyMatrix.h"
#include "simc/SimcConstants.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "simc/SimcConfig.h"

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
        
        // Note: This test requires actual COSY matrix files
        // For now, we test the structure
        
        CosyMatrix matrix;
        ASSERT_FALSE(matrix.IsLoaded(), "Matrix should not be loaded initially");
        ASSERT_EQ(static_cast<double>(matrix.GetNumElements()), 0.0, "Should have 0 elements");
        ASSERT_EQ(static_cast<double>(matrix.GetMaxOrder()), 0.0, "Max order should be 0");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 2: COSY Matrix Parsing
    // ========================================================================
    {
        std::cout << "    Test 2: COSY Matrix Element Parsing... ";
        
        CosyMatrix::MatrixElement elem;
        elem.output_index = 0;  // x output
        elem.coefficient = 1.5;
        elem.exponents[0] = 2;  // x^2
        elem.exponents[1] = 1;  // xp^1
        elem.exponents[2] = 0;
        elem.exponents[3] = 0;
        elem.exponents[4] = 0;
        elem.order = 3;
        
        // Test evaluation: 1.5 * x^2 * xp^1
        double input[5] = {2.0, 3.0, 0.0, 0.0, 0.0};  // x=2, xp=3
        double result = elem.Evaluate(input);
        
        // Should be: 1.5 * 2^2 * 3 = 1.5 * 4 * 3 = 18.0
        ASSERT_NEAR(result, 18.0, 1e-10, "Element evaluation");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 3: CRITICAL - Delta Units (PERCENT vs Fraction)
    // ========================================================================
    {
        std::cout << "    Test 3: Delta Units Verification... ";
        
        // This is CRITICAL: Fortran SIMC uses PERCENT for delta!
        // From transp.f line 45: "delta is in PERCENT, NOT fractional"
        
        double p_central = 5000.0;  // 5 GeV/c
        double p_actual = 5250.0;   // 5.25 GeV/c
        
        // CORRECT calculation (PERCENT):
        double delta_percent = (p_actual - p_central) / p_central * 100.0;
        ASSERT_NEAR(delta_percent, 5.0, 1e-10, "Delta should be 5%");
        
        // WRONG calculation (fraction) - DO NOT USE:
        double delta_fraction = (p_actual - p_central) / p_central;
        ASSERT_NEAR(delta_fraction, 0.05, 1e-10, "Fraction is 0.05");
        
        // Verify they are different by factor of 100
        ASSERT_NEAR(delta_percent, delta_fraction * 100.0, 1e-10, 
                   "Percent = Fraction * 100");
        
        std::cout << "PASSED (Delta in PERCENT: " << delta_percent << "%)\n";
    }
    
    // ========================================================================
    // Test 4: HMS Octagon Aperture
    // ========================================================================
    {
        std::cout << "    Test 4: HMS Octagon Aperture... ";
        
        // From apertures_hms.inc: octagon radius = 17.145 cm
        HMSOptics hms("", "");  // Empty files for geometry-only test
        
        // Point inside octagon (0, 0) - should pass
        auto check1 = hms.CheckAperture(0.0, 0.0, 0);
        ASSERT_TRUE(check1.passed, "Origin should be inside");
        
        // Point on edge (17, 0) - should pass
        auto check2 = hms.CheckAperture(17.0, 0.0, 0);
        ASSERT_TRUE(check2.passed, "Point at r=17 cm should be inside");
        
        // Point outside octagon (20, 20) - should fail
        auto check3 = hms.CheckAperture(20.0, 20.0, 0);
        ASSERT_FALSE(check3.passed, "Point at (20,20) should be outside");
        
        // Point on 45° corner - test octagon vs circle
        // Octagon allows |x| + |y| <= r*sqrt(2) ≈ 24.25 cm
        auto check4 = hms.CheckAperture(12.0, 12.0, 0);
        ASSERT_TRUE(check4.passed, "Point (12,12) should be inside octagon");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 5: HMS Dipole Aperture
    // ========================================================================
    {
        std::cout << "    Test 5: HMS Dipole Aperture... ";
        
        // From apertures_hms.inc: dipole ±30 cm horizontal, ±12.5 cm vertical
        HMSOptics hms("", "");
        
        // Inside dipole
        auto check1 = hms.CheckAperture(25.0, 10.0, 1);
        ASSERT_TRUE(check1.passed, "(25,10) inside dipole");
        
        // Outside horizontal
        auto check2 = hms.CheckAperture(35.0, 5.0, 1);
        ASSERT_FALSE(check2.passed, "(35,5) outside dipole horizontal");
        
        // Outside vertical
        auto check3 = hms.CheckAperture(10.0, 15.0, 1);
        ASSERT_FALSE(check3.passed, "(10,15) outside dipole vertical");
        
        // On edge
        auto check4 = hms.CheckAperture(30.0, 12.5, 1);
        ASSERT_TRUE(check4.passed, "(30,12.5) on dipole edge");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 6: HMS Quadrupole Apertures
    // ========================================================================
    {
        std::cout << "    Test 6: HMS Quadrupole Apertures... ";
        
        HMSOptics hms("", "");
        
        // Q1: r = 12.5 cm (plane 3)
        auto check_q1_in = hms.CheckAperture(8.0, 8.0, 3);
        ASSERT_TRUE(check_q1_in.passed, "Inside Q1");
        
        auto check_q1_out = hms.CheckAperture(10.0, 10.0, 3);
        // r = sqrt(10^2 + 10^2) = 14.14 cm > 12.5 cm, should be outside
        ASSERT_FALSE(check_q1_out.passed, "Outside Q1");
        
        // Q2: r = 30 cm (plane 5)
        auto check_q2_in = hms.CheckAperture(20.0, 15.0, 5);
        ASSERT_TRUE(check_q2_in.passed, "Inside Q2");
        
        auto check_q2_out = hms.CheckAperture(25.0, 25.0, 5);
        // r = sqrt(25^2 + 25^2) = 35.35 cm > 30 cm, should be outside
        ASSERT_FALSE(check_q2_out.passed, "Outside Q2");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 7: SHMS Collimator Apertures
    // ========================================================================
    {
        std::cout << "    Test 7: SHMS Collimator Apertures... ";
        
        SHMSOptics shms("", "");
        
        // Test LARGE collimator (±6 cm)
        shms.SetCollimator("LARGE");
        auto check_large_in = shms.CheckAperture(5.0, 5.0, 6);
        ASSERT_TRUE(check_large_in.passed, "Inside LARGE collimator");
        
        auto check_large_out = shms.CheckAperture(7.0, 5.0, 6);
        ASSERT_FALSE(check_large_out.passed, "Outside LARGE collimator");
        
        // Test SMALL collimator (±3 cm)
        shms.SetCollimator("SMALL");
        auto check_small_in = shms.CheckAperture(2.0, 2.0, 6);
        ASSERT_TRUE(check_small_in.passed, "Inside SMALL collimator");
        
        auto check_small_out = shms.CheckAperture(4.0, 2.0, 6);
        ASSERT_FALSE(check_small_out.passed, "Outside SMALL collimator");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 8: SHMS HB Dipole Aperture
    // ========================================================================
    {
        std::cout << "    Test 8: SHMS HB Dipole Aperture... ";
        
        // From apertures_shms.inc: HB dipole ±50 cm horizontal, ±25 cm vertical
        SHMSOptics shms("", "");
        
        auto check_in = shms.CheckAperture(40.0, 20.0, 1);
        ASSERT_TRUE(check_in.passed, "Inside HB dipole");
        
        auto check_out_x = shms.CheckAperture(55.0, 20.0, 1);
        ASSERT_FALSE(check_out_x.passed, "Outside HB dipole horizontal");
        
        auto check_out_y = shms.CheckAperture(40.0, 30.0, 1);
        ASSERT_FALSE(check_out_y.passed, "Outside HB dipole vertical");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 9: Coordinate System Convention
    // ========================================================================
    {
        std::cout << "    Test 9: TRANSPORT Coordinate Convention... ";
        
        // Verify TRANSPORT convention:
        // x = horizontal (dispersive), positive = beam right
        // y = vertical, positive = up
        // xp = dx/dz, yp = dy/dz in radians
        // delta in PERCENT
        
        ArmState state;
        state.xptar = 0.02;   // 20 mrad horizontal
        state.yptar = 0.03;   // 30 mrad vertical
        state.delta = 5.0;    // 5% momentum deviation
        
        // Verify the state is set correctly
        ASSERT_NEAR(state.xptar, 0.02, 1e-10, "xptar set");
        ASSERT_NEAR(state.yptar, 0.03, 1e-10, "yptar set");
        ASSERT_NEAR(state.delta, 5.0, 1e-10, "delta in PERCENT");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 10: SOS Apertures
    // ========================================================================
    {
        std::cout << "    Test 10: SOS Apertures... ";
        
        SOSOptics sos("", "");
        
        // Entrance (±30 cm x, ±15 cm y)
        auto check_entrance = sos.CheckAperture(25.0, 10.0, 0);
        ASSERT_TRUE(check_entrance.passed, "Inside SOS entrance");
        
        auto check_entrance_out = sos.CheckAperture(35.0, 10.0, 0);
        ASSERT_FALSE(check_entrance_out.passed, "Outside SOS entrance");
        
        // Q1 (r = 15 cm)
        auto check_q1 = sos.CheckAperture(10.0, 10.0, 2);
        ASSERT_TRUE(check_q1.passed, "Inside SOS Q1");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 11: HRS Apertures
    // ========================================================================
    {
        std::cout << "    Test 11: HRS Apertures... ";
        
        HRSOptics hrs("HRSL", "", "");
        
        // Q1 (r = 12.5 cm)
        auto check_q1_in = hrs.CheckAperture(8.0, 8.0, 1);
        ASSERT_TRUE(check_q1_in.passed, "Inside HRS Q1");
        
        auto check_q1_out = hrs.CheckAperture(10.0, 10.0, 1);
        ASSERT_FALSE(check_q1_out.passed, "Outside HRS Q1");
        
        // Dipole (±40 cm x, ±12 cm y)
        auto check_dipole = hrs.CheckAperture(35.0, 10.0, 2);
        ASSERT_TRUE(check_dipole.passed, "Inside HRS dipole");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 12: Factory Function
    // ========================================================================
    {
        std::cout << "    Test 12: Spectrometer Factory... ";
	
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
    // Test 13: Default Matrix File Paths
    // ========================================================================
    {
        std::cout << "    Test 13: Default Matrix Files... ";
        
        auto [hms_fwd, hms_rec] = GetDefaultMatrixFiles("HMS");
        ASSERT_TRUE(hms_fwd.find("hms") != std::string::npos, "HMS forward path");
        ASSERT_TRUE(hms_rec.find("hms") != std::string::npos, "HMS recon path");
        
        auto [shms_fwd, shms_rec] = GetDefaultMatrixFiles("SHMS");
        ASSERT_TRUE(shms_fwd.find("shms") != std::string::npos, "SHMS forward path");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 14: Central Momentum and Angle
    // ========================================================================
    {
        std::cout << "    Test 14: Spectrometer Settings... ";
        
        auto hms = CreateSpectrometerOptics("HMS");
        
        // Test default values
        ASSERT_NEAR(hms->GetCentralMomentum(), 5000.0, 1e-6, "Default p_central");
        ASSERT_NEAR(hms->GetAngle(), 0.0, 1e-10, "Default angle");
        
        // Test setters
        hms->SetCentralMomentum(7500.0);
        ASSERT_NEAR(hms->GetCentralMomentum(), 7500.0, 1e-6, "Set p_central");
        
        hms->SetAngle(0.218);  // ~12.5 degrees
        ASSERT_NEAR(hms->GetAngle(), 0.218, 1e-6, "Set angle");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 15: Delta and Momentum Relationship
    // ========================================================================
    {
        std::cout << "    Test 15: Delta and Momentum Relationship... ";
        
        HMSOptics hms("", "");
        hms.SetCentralMomentum(5000.0);
        
        // Create test state
        ArmState state_in;
        state_in.xptar = 0.025;  // 25 mrad
        state_in.yptar = 0.035;  // 35 mrad
        state_in.delta = 7.5;    // 7.5%
        state_in.P = 5000.0 * (1.0 + 7.5 / 100.0); // 5375 MeV/c
        
        // Verify momentum calculation from delta
        double p_expected = 5000.0 * (1.0 + 7.5 / 100.0);
        ASSERT_NEAR(state_in.P, p_expected, 1e-6, "Momentum from delta");
        ASSERT_NEAR(state_in.P, 5375.0, 1e-6, "Momentum value");
        
        // Verify reverse calculation
        double delta_calc = (state_in.P - 5000.0) / 5000.0 * 100.0;
        ASSERT_NEAR(delta_calc, 7.5, 1e-10, "Delta from momentum");
        
        std::cout << "PASSED\n";
    }
    
    // ========================================================================
    // Test 16: Multiple Aperture Planes
    // ========================================================================
    {
        std::cout << "    Test 16: Multiple Aperture Planes... ";
        
        auto hms = CreateSpectrometerOptics("HMS");
        ASSERT_EQ(static_cast<double>(hms->GetNumAperturePlanes()), 8.0, 
                 "HMS has 8 planes");
        
        auto shms = CreateSpectrometerOptics("SHMS");
        ASSERT_EQ(static_cast<double>(shms->GetNumAperturePlanes()), 9.0, 
                 "SHMS has 9 planes");
        
        auto sos = CreateSpectrometerOptics("SOS");
        ASSERT_EQ(static_cast<double>(sos->GetNumAperturePlanes()), 6.0, 
                 "SOS has 6 planes");
        
        auto hrs = CreateSpectrometerOptics("HRSL");
        ASSERT_EQ(static_cast<double>(hrs->GetNumAperturePlanes()), 7.0, 
                 "HRS has 7 planes");
        
        std::cout << "PASSED\n";
    }
    
    std::cout << "\n  All spectrometer optics tests PASSED!\n";
}
