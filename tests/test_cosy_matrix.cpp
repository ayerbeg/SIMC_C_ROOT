// test_cosy_matrix.cpp
// Simple test program to verify COSY matrix fixes
// Compile: g++ -std=c++17 -I../include test_cosy_matrix.cpp ../src/physics/CosyMatrix.cpp -o test_cosy

#include "simc/CosyMatrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace simc;

void test_basic_load() {
    std::cout << "\n=== Test 1: Basic Matrix Loading ===" << std::endl;
    
    CosyMatrix matrix;
    bool success = matrix.LoadFromFile("../data/matrices/hms/forward_cosy.dat");
    
    if (success) {
        std::cout << "✓ Matrix loaded successfully" << std::endl;
        std::cout << "  Elements: " << matrix.GetNumElements() << std::endl;
        std::cout << "  Max order: " << matrix.GetMaxOrder() << std::endl;
        if (matrix.IsDrift()) {
            std::cout << "  Type: Drift (" << matrix.GetDriftDistance() << " cm)" << std::endl;
        }
    } else {
        std::cout << "✗ Failed to load matrix" << std::endl;
    }
}

void test_simple_transport() {
    std::cout << "\n=== Test 2: Simple Transport ===" << std::endl;
    
    CosyMatrix matrix;
    if (!matrix.LoadFromFile("../data/matrices/hms/forward_cosy.dat")) {
        std::cout << "✗ Cannot load matrix" << std::endl;
        return;
    }
    
    // Test with simple on-axis particle
    // Units: x(cm), xp(mrad), y(cm), yp(mrad), delta(%)
    double input[5] = {
        0.0,    // x = 0 cm (on axis)
        10.0,   // xp = 10 mrad = 0.010 rad
        0.0,    // y = 0 cm (on axis)
        5.0,    // yp = 5 mrad = 0.005 rad
        0.0     // delta = 0% (central momentum)
    };
    
    double output[5];
    matrix.Apply(input, output);
    
    std::cout << "Input (target):" << std::endl;
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
    std::cout << "  dL = " << output[4] << " cm" << std::endl;
    
    // Sanity checks
    bool sane = true;
    if (std::abs(output[0]) > 100.0) {
        std::cout << "⚠ Warning: x_fp seems too large" << std::endl;
        sane = false;
    }
    if (std::abs(output[2]) > 100.0) {
        std::cout << "⚠ Warning: y_fp seems too large" << std::endl;
        sane = false;
    }
    if (std::abs(output[4]) > 500.0) {
        std::cout << "⚠ Warning: path length correction seems too large" << std::endl;
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
        std::cout << "✗ Cannot load reconstruction matrix" << std::endl;
        return;
    }
    
    // Focal plane coordinates (example)
    double fp_input[5] = {
        5.2,    // x_fp (cm)
        -8.5,   // xp_fp (mrad)
        2.1,    // y_fp (cm)
        3.7,    // yp_fp (mrad)
        0.0     // delta (%)
    };
    
    double target_output[5];
    recon_matrix.Apply(fp_input, target_output);
    
    std::cout << "Focal plane input:" << std::endl;
    std::cout << "  x  = " << std::fixed << std::setprecision(4) << fp_input[0] << " cm" << std::endl;
    std::cout << "  xp = " << fp_input[1] << " mrad" << std::endl;
    std::cout << "  y  = " << fp_input[2] << " cm" << std::endl;
    std::cout << "  yp = " << fp_input[3] << " mrad" << std::endl;
    
    std::cout << "\nReconstructed target:" << std::endl;
    std::cout << "  x_tar  = " << target_output[0] << " cm (should be ~0)" << std::endl;
    std::cout << "  xp_tar = " << target_output[1] << " mrad (" 
              << target_output[1]/1000.0 << " rad)" << std::endl;
    std::cout << "  y_tar  = " << target_output[2] << " cm (should be ~0)" << std::endl;
    std::cout << "  yp_tar = " << target_output[3] << " mrad (" 
              << target_output[3]/1000.0 << " rad)" << std::endl;
    std::cout << "  δ_tar  = " << target_output[4] << " %" << std::endl;
    
    std::cout << "✓ Reconstruction complete" << std::endl;
}

void test_round_trip() {
    std::cout << "\n=== Test 4: Round-Trip (Forward + Reconstruction) ===" << std::endl;
    
    CosyMatrix forward, recon;
    if (!forward.LoadFromFile("../data/matrices/hms/forward_cosy.dat")) {
        std::cout << "✗ Cannot load forward matrix" << std::endl;
        return;
    }
    if (!recon.LoadFromFile("../data/matrices/hms/recon_cosy.dat")) {
        std::cout << "✗ Cannot load recon matrix" << std::endl;
        return;
    }
    
    // Original target angles
    double original_xptar = 15.0;  // mrad
    double original_yptar = 8.0;   // mrad
    double original_delta = 2.5;   // %
    
    // Forward transport
    double target_in[5] = {0.0, original_xptar, 0.0, original_yptar, original_delta};
    double fp_out[5];
    forward.Apply(target_in, fp_out);
    
    // Reconstruction
    double recon_in[5] = {fp_out[0], fp_out[1], fp_out[2], fp_out[3], original_delta};
    double recon_out[5];
    recon.Apply(recon_in, recon_out);
    
    std::cout << "Original target angles:" << std::endl;
    std::cout << "  xp = " << std::fixed << std::setprecision(4) << original_xptar << " mrad" << std::endl;
    std::cout << "  yp = " << original_yptar << " mrad" << std::endl;
    std::cout << "  δ  = " << original_delta << " %" << std::endl;
    
    std::cout << "\nFocal plane:" << std::endl;
    std::cout << "  x  = " << fp_out[0] << " cm" << std::endl;
    std::cout << "  y  = " << fp_out[2] << " cm" << std::endl;
    
    std::cout << "\nReconstructed target angles:" << std::endl;
    std::cout << "  xp = " << recon_out[1] << " mrad" << std::endl;
    std::cout << "  yp = " << recon_out[3] << " mrad" << std::endl;
    std::cout << "  δ  = " << recon_out[4] << " %" << std::endl;
    
    // Check accuracy
    double xp_diff = std::abs(recon_out[1] - original_xptar);
    double yp_diff = std::abs(recon_out[3] - original_yptar);
    double delta_diff = std::abs(recon_out[4] - original_delta);
    
    std::cout << "\nReconstruction errors:" << std::endl;
    std::cout << "  Δxp = " << xp_diff << " mrad" << std::endl;
    std::cout << "  Δyp = " << yp_diff << " mrad" << std::endl;
    std::cout << "  Δδ  = " << delta_diff << " %" << std::endl;
    
    // Should be very small (< 0.1 mrad for angles, < 0.01% for delta)
    if (xp_diff < 0.1 && yp_diff < 0.1 && delta_diff < 0.01) {
        std::cout << "✓ Round-trip successful! Errors are acceptable." << std::endl;
    } else {
        std::cout << "⚠ Warning: Round-trip errors seem large" << std::endl;
    }
}

void test_print_matrix() {
    std::cout << "\n=== Test 5: Matrix Details ===" << std::endl;
    
    CosyMatrix matrix;
    if (!matrix.LoadFromFile("../data/matrices/hms/forward_cosy.dat")) {
        std::cout << "✗ Cannot load matrix" << std::endl;
        return;
    }
    
    matrix.Print();
}

int main(int argc, char** argv) {
    std::cout << "========================================" << std::endl;
    std::cout << "SIMC++ COSY Matrix Test Suite" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Run all tests
    test_basic_load();
    test_simple_transport();
    test_reconstruction();
    test_round_trip();
    
    if (argc > 1 && std::string(argv[1]) == "-v") {
        test_print_matrix();
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "All tests complete!" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    return 0;
}
