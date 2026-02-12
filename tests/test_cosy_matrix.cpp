// tests/test_cosy_matrix.cpp
// Test program to verify COSY matrix parser and transport
// Tests the CORRECTED parser that properly handles Fortran format

#include "simc/CosyMatrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace simc;

void test_basic_load() {
    std::cout << "\n=== Test 1: Basic Matrix Loading ===" << std::endl;
    
    CosyMatrix matrix;
    
    // Try to load HMS forward matrix
    bool success = matrix.LoadFromFile("../data/matrices/hms/forward_cosy.dat");
    
    if (success) {
        std::cout << "✓ Matrix loaded successfully" << std::endl;
        std::cout << "  Elements: " << matrix.GetNumElements() << std::endl;
        std::cout << "  Max order: " << matrix.GetMaxOrder() << std::endl;
        std::cout << "  Is loaded: " << (matrix.IsLoaded() ? "Yes" : "No") << std::endl;
        
        if (matrix.IsDrift()) {
            std::cout << "  Type: Drift (" << matrix.GetDriftDistance() << " cm)" << std::endl;
        }
        
        // Sanity check: HMS matrices typically have 300-600 elements
        if (matrix.GetNumElements() < 100) {
            std::cout << "⚠ WARNING: Element count seems too low!" << std::endl;
            std::cout << "  This suggests the parser may not be working correctly." << std::endl;
        }
        
    } else {
        std::cout << "✗ Failed to load matrix" << std::endl;
        std::cout << "  Make sure you're running from the build directory!" << std::endl;
    }
}

void test_simple_transport() {
    std::cout << "\n=== Test 2: Simple Transport ===" << std::endl;
    
    CosyMatrix matrix;
    if (!matrix.LoadFromFile("../data/matrices/hms/forward_cosy.dat")) {
        std::cout << "✗ Cannot load matrix - skipping test" << std::endl;
        return;
    }
    
    std::cout << "Loaded " << matrix.GetNumElements() << " matrix elements" << std::endl;
    
    // Test with simple on-axis particle
    // CRITICAL UNITS (from transp.f):
    //   x, y: cm
    //   xp, yp: mrad (milliradians) - NOT radians!
    //   delta: percent
    
    double input[5] = {
        0.0,    // x = 0 cm (on axis)
        10.0,   // xp = 10 mrad = 0.010 rad
        0.0,    // y = 0 cm (on axis)
        5.0,    // yp = 5 mrad = 0.005 rad
        0.0     // delta = 0% (central momentum)
    };
    
    double output[5];
    matrix.Apply(input, output);
    
    std::cout << "\nInput (target coordinates):" << std::endl;
    std::cout << "  x  = " << std::fixed << std::setprecision(4) << input[0] << " cm" << std::endl;
    std::cout << "  xp = " << input[1] << " mrad (" << input[1]/1000.0 << " rad)" << std::endl;
    std::cout << "  y  = " << input[2] << " cm" << std::endl;
    std::cout << "  yp = " << input[3] << " mrad (" << input[3]/1000.0 << " rad)" << std::endl;
    std::cout << "  δ  = " << input[4] << " %" << std::endl;
    
    std::cout << "\nOutput (focal plane):" << std::endl;
    std::cout << "  x  = " << output[0] << " cm" << std::endl;
    std::cout << "  xp = " << output[1] << " mrad (" << output[1]/1000.0 << " rad)" << std::endl;
    std::cout << "  y  = " << output[2] << " cm" << std::endl;
    std::cout << "  yp = " << output[3] << " mrad (" << output[3]/1000.0 << " rad)" << std::endl;
    std::cout << "  dL = " << output[4] << " cm (path length correction)" << std::endl;
    
    // Sanity checks for HMS
    // Expected: focal plane positions should be reasonable
    bool sane = true;
    if (std::abs(output[0]) > 200.0) {
        std::cout << "⚠ Warning: x_fp = " << output[0] << " cm seems too large" << std::endl;
        sane = false;
    }
    if (std::abs(output[2]) > 200.0) {
        std::cout << "⚠ Warning: y_fp = " << output[2] << " cm seems too large" << std::endl;
        sane = false;
    }
    if (std::abs(output[4]) > 1000.0) {
        std::cout << "⚠ Warning: dL = " << output[4] << " cm seems too large" << std::endl;
        sane = false;
    }
    
    if (sane) {
        std::cout << "✓ Output values look reasonable" << std::endl;
    }
}

void test_reconstruction() {
    std::cout << "\n=== Test 3: Reconstruction ===" << std::endl;
    
    CosyMatrix recon_matrix;
    if (!recon_matrix.LoadFromFile("../data/matrices/hms/recon_cosy.dat")) {
        std::cout << "✗ Cannot load reconstruction matrix - skipping test" << std::endl;
        return;
    }
    
    std::cout << "Loaded " << recon_matrix.GetNumElements() << " reconstruction elements" << std::endl;
    
    // Focal plane coordinates (example)
    // Units: cm, mrad, cm, mrad, %
    double fp_input[5] = {
        5.2,    // x_fp (cm)
        -8.5,   // xp_fp (mrad)
        2.1,    // y_fp (cm)
        3.7,    // yp_fp (mrad)
        1.5     // delta (%) - need to provide this from tracking/TOF
    };
    
    double target_output[5];
    recon_matrix.Apply(fp_input, target_output);
    
    std::cout << "\nFocal plane input:" << std::endl;
    std::cout << "  x  = " << std::fixed << std::setprecision(4) << fp_input[0] << " cm" << std::endl;
    std::cout << "  xp = " << fp_input[1] << " mrad" << std::endl;
    std::cout << "  y  = " << fp_input[2] << " cm" << std::endl;
    std::cout << "  yp = " << fp_input[3] << " mrad" << std::endl;
    std::cout << "  δ  = " << fp_input[4] << " %" << std::endl;
    
    std::cout << "\nReconstructed target:" << std::endl;
    std::cout << "  x_tar  = " << target_output[0] << " cm (should be ~0)" << std::endl;
    std::cout << "  xp_tar = " << target_output[1] << " mrad (" 
              << target_output[1]/1000.0 << " rad)" << std::endl;
    std::cout << "  y_tar  = " << target_output[2] << " cm (should be ~0)" << std::endl;
    std::cout << "  yp_tar = " << target_output[3] << " mrad (" 
              << target_output[3]/1000.0 << " rad)" << std::endl;
    std::cout << "  δ_tar  = " << target_output[4] << " %" << std::endl;
    
    // Sanity check: target x and y should be near zero
    if (std::abs(target_output[0]) < 1.0 && std::abs(target_output[2]) < 1.0) {
        std::cout << "✓ Reconstructed target positions near zero (as expected)" << std::endl;
    } else {
        std::cout << "⚠ Warning: Target positions seem off-axis" << std::endl;
    }
}

void test_round_trip() {
    std::cout << "\n=== Test 4: Round-Trip (Forward + Reconstruction) ===" << std::endl;
    
    CosyMatrix forward, recon;
    if (!forward.LoadFromFile("../data/matrices/hms/forward_cosy.dat")) {
        std::cout << "✗ Cannot load forward matrix - skipping test" << std::endl;
        return;
    }
    if (!recon.LoadFromFile("../data/matrices/hms/recon_cosy.dat")) {
        std::cout << "✗ Cannot load recon matrix - skipping test" << std::endl;
        return;
    }
    
    // Original target angles (in mrad!)
    double original_xptar = 15.0;  // mrad
    double original_yptar = 8.0;   // mrad
    double original_delta = 2.5;   // %
    
    std::cout << "Original target coordinates:" << std::endl;
    std::cout << "  x_tar  = 0.0 cm" << std::endl;
    std::cout << "  xp_tar = " << std::fixed << std::setprecision(4) 
              << original_xptar << " mrad (" << original_xptar/1000.0 << " rad)" << std::endl;
    std::cout << "  y_tar  = 0.0 cm" << std::endl;
    std::cout << "  yp_tar = " << original_yptar << " mrad (" 
              << original_yptar/1000.0 << " rad)" << std::endl;
    std::cout << "  δ      = " << original_delta << " %" << std::endl;
    
    // Forward transport
    double target_in[5] = {0.0, original_xptar, 0.0, original_yptar, original_delta};
    double fp_out[5];
    forward.Apply(target_in, fp_out);
    
    std::cout << "\nFocal plane (after forward transport):" << std::endl;
    std::cout << "  x_fp  = " << fp_out[0] << " cm" << std::endl;
    std::cout << "  xp_fp = " << fp_out[1] << " mrad" << std::endl;
    std::cout << "  y_fp  = " << fp_out[2] << " cm" << std::endl;
    std::cout << "  yp_fp = " << fp_out[3] << " mrad" << std::endl;
    std::cout << "  dL    = " << fp_out[4] << " cm" << std::endl;
    
    // Reconstruction
    // Note: In real data analysis, delta would come from tracking/TOF
    // Here we use the known delta for testing
    double recon_in[5] = {fp_out[0], fp_out[1], fp_out[2], fp_out[3], original_delta};
    double recon_out[5];
    recon.Apply(recon_in, recon_out);
    
    std::cout << "\nReconstructed target coordinates:" << std::endl;
    std::cout << "  x_tar  = " << recon_out[0] << " cm" << std::endl;
    std::cout << "  xp_tar = " << recon_out[1] << " mrad (" 
              << recon_out[1]/1000.0 << " rad)" << std::endl;
    std::cout << "  y_tar  = " << recon_out[2] << " cm" << std::endl;
    std::cout << "  yp_tar = " << recon_out[3] << " mrad (" 
              << recon_out[3]/1000.0 << " rad)" << std::endl;
    std::cout << "  δ      = " << recon_out[4] << " %" << std::endl;
    
    // Check accuracy
    double xp_diff = std::abs(recon_out[1] - original_xptar);
    double yp_diff = std::abs(recon_out[3] - original_yptar);
    double delta_diff = std::abs(recon_out[4] - original_delta);
    
    std::cout << "\nRound-trip reconstruction errors:" << std::endl;
    std::cout << "  Δxp = " << xp_diff << " mrad" << std::endl;
    std::cout << "  Δyp = " << yp_diff << " mrad" << std::endl;
    std::cout << "  Δδ  = " << delta_diff << " %" << std::endl;
    
    // Typical reconstruction errors should be small
    // < 1 mrad for angles, < 0.1% for delta
    if (xp_diff < 1.0 && yp_diff < 1.0 && delta_diff < 0.1) {
        std::cout << "✓ Round-trip successful! Errors are acceptable." << std::endl;
    } else {
        std::cout << "⚠ Warning: Round-trip errors are larger than expected" << std::endl;
        std::cout << "  This might indicate a problem with the matrices or transport." << std::endl;
    }
}

void test_zero_angle() {
    std::cout << "\n=== Test 5: Zero Angle Particle (Central Ray) ===" << std::endl;
    
    CosyMatrix matrix;
    if (!matrix.LoadFromFile("../data/matrices/hms/forward_cosy.dat")) {
        std::cout << "✗ Cannot load matrix - skipping test" << std::endl;
        return;
    }
    
    // Particle on central ray (all zeros)
    double input[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double output[5];
    matrix.Apply(input, output);
    
    std::cout << "Central ray transport:" << std::endl;
    std::cout << "  Input:  all zeros" << std::endl;
    std::cout << "  Output: x=" << std::fixed << std::setprecision(6) << output[0] 
              << ", xp=" << output[1] << ", y=" << output[2] 
              << ", yp=" << output[3] << ", dL=" << output[4] << std::endl;
    
    // For central ray, output should also be near zero (within numerical precision)
    double max_val = std::max({std::abs(output[0]), std::abs(output[1]), 
                                std::abs(output[2]), std::abs(output[3])});
    
    if (max_val < 1e-10) {
        std::cout << "✓ Central ray transported correctly (all outputs ~0)" << std::endl;
    } else {
        std::cout << "⚠ Warning: Central ray has non-zero output" << std::endl;
    }
}

void test_print_matrix() {
    std::cout << "\n=== Test 6: Matrix Details ===" << std::endl;
    
    CosyMatrix matrix;
    if (!matrix.LoadFromFile("../data/matrices/hms/forward_cosy.dat")) {
        std::cout << "✗ Cannot load matrix - skipping test" << std::endl;
        return;
    }
    
    matrix.Print();
}

int main(int argc, char** argv) {
    std::cout << "========================================" << std::endl;
    std::cout << "SIMC++ COSY Matrix Test Suite" << std::endl;
    std::cout << "Testing CORRECTED parser" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Check if we should run verbose mode
    bool verbose = false;
    if (argc > 1 && std::string(argv[1]) == "-v") {
        verbose = true;
    }
    
    // Run all tests
    test_basic_load();
    test_simple_transport();
    test_zero_angle();
    test_reconstruction();
    test_round_trip();
    
    if (verbose) {
        test_print_matrix();
    } else {
        std::cout << "\nRun with -v flag to see detailed matrix information" << std::endl;
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "All tests complete!" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    return 0;
}
