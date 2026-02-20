/**
 * @file test_matrix_delta_terms.cpp
 * @brief Check if forward matrices contain delta-dependent terms
 */

#include "simc/physics/CosyMatrix.h"
#include <iostream>
#include <iomanip>

using namespace simc;

int main() {
    std::cout << "\n=== FORWARD MATRIX DELTA DEPENDENCE CHECK ===\n" << std::endl;
    
    CosyMatrix matrix;
    
    if (!matrix.LoadFromFile("../data/matrices/hms/forward_cosy.dat")) {
        std::cerr << "Failed to load forward matrix!" << std::endl;
        return 1;
    }
    
    std::cout << "Loaded forward matrix with " << matrix.GetNumElements() << " elements\n" << std::endl;
    
    // Check how many elements have delta dependence
    const auto& elements = matrix.GetElements();
    
    int count_delta_dependent = 0;
    int count_first_order_delta = 0;
    int count_higher_order_delta = 0;
    
    for (const auto& elem : elements) {
        int delta_exp = elem.exponents[4];  // Delta is 5th variable (index 4)
        
        if (delta_exp > 0) {
            count_delta_dependent++;
            
            if (elem.order == 1) {
                count_first_order_delta++;
            } else {
                count_higher_order_delta++;
            }
            
            // Print first few delta-dependent terms
            if (count_delta_dependent <= 10) {
                std::cout << "Term " << count_delta_dependent << ": ";
                std::cout << "exponents=[" << elem.exponents[0] << "," << elem.exponents[1] 
                          << "," << elem.exponents[2] << "," << elem.exponents[3] 
                          << "," << elem.exponents[4] << "], ";
                std::cout << "order=" << elem.order << ", ";
                std::cout << "coeffs=[";
                for (int i = 0; i < 5; ++i) {
                    std::cout << std::scientific << std::setprecision(2) << elem.coefficients[i];
                    if (i < 4) std::cout << ", ";
                }
                std::cout << "]" << std::endl;
            }
        }
    }
    
    std::cout << "\n=== SUMMARY ===" << std::endl;
    std::cout << "Total elements: " << elements.size() << std::endl;
    std::cout << "Delta-dependent terms: " << count_delta_dependent << std::endl;
    std::cout << "  First-order: " << count_first_order_delta << std::endl;
    std::cout << "  Higher-order: " << count_higher_order_delta << std::endl;
    
    if (count_delta_dependent == 0) {
        std::cout << "\n⚠️  WARNING: Forward matrix has NO delta-dependent terms!" << std::endl;
        std::cout << "This explains why delta doesn't affect focal plane position." << std::endl;
    } else {
        std::cout << "\n✓ Forward matrix has delta-dependent terms." << std::endl;
        std::cout << "The problem must be in how delta is being applied." << std::endl;
    }
    
    // Test a simple case manually
    std::cout << "\n=== MANUAL TEST ===" << std::endl;
    std::cout << "Testing if matrix output changes with different delta values:\n" << std::endl;
    
    double input_zero_delta[5] = {0.0, 0.0, 0.0, 0.0, 0.0};  // All zeros
    double output_zero_delta[5];
    matrix.Apply(input_zero_delta, output_zero_delta);
    
    double input_plus_delta[5] = {0.0, 0.0, 0.0, 0.0, 5.0};  // delta = 5%
    double output_plus_delta[5];
    matrix.Apply(input_plus_delta, output_plus_delta);
    
    std::cout << "Input: x=0, xp=0, y=0, yp=0, delta=0%" << std::endl;
    std::cout << "Output: x=" << std::fixed << std::setprecision(6) << output_zero_delta[0] 
              << ", xp=" << output_zero_delta[1] 
              << ", y=" << output_zero_delta[2]
              << ", yp=" << output_zero_delta[3]
              << ", dL=" << output_zero_delta[4] << std::endl;
    
    std::cout << "\nInput: x=0, xp=0, y=0, yp=0, delta=5%" << std::endl;
    std::cout << "Output: x=" << output_plus_delta[0] 
              << ", xp=" << output_plus_delta[1]
              << ", y=" << output_plus_delta[2]
              << ", yp=" << output_plus_delta[3]
              << ", dL=" << output_plus_delta[4] << std::endl;
    
    std::cout << "\nDifferences:" << std::endl;
    std::cout << "  Δx = " << (output_plus_delta[0] - output_zero_delta[0]) << " (should be non-zero)" << std::endl;
    std::cout << "  Δxp = " << (output_plus_delta[1] - output_zero_delta[1]) << " (should be non-zero)" << std::endl;
    
    return 0;
}
