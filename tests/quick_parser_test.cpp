// tests/quick_parser_test.cpp
// Quick test to verify COSY parser handles all format variations
// This test uses the actual CosyMatrix class

#include "simc/CosyMatrix.h"
#include <iostream>
#include <iomanip>

void TestSingleLine(const std::string& test_line, const std::string& description) {
    simc::CosyMatrix matrix;
    std::vector<double> coeffs;
    std::vector<int> exponents;
    
    std::cout << "\n=== Testing: " << description << " ===" << std::endl;
    std::cout << "Input: " << test_line << std::endl;
    
    if (matrix.ParseMatrixLine(test_line, coeffs, exponents)) {
        std::cout << "✓ Parsed successfully!" << std::endl;
        std::cout << "  Coefficients (" << coeffs.size() << "): ";
        for (size_t i = 0; i < coeffs.size(); ++i) {
            std::cout << std::setw(12) << coeffs[i];
            if (i < coeffs.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
        std::cout << "  Exponents: ";
        for (int e : exponents) {
            std::cout << e;
        }
        std::cout << std::endl;
    } else {
        std::cout << "❌ Parse failed!" << std::endl;
    }
}

void TestFile(const std::string& filename) {
    std::cout << "\n=== Testing File: " << filename << " ===" << std::endl;
    
    simc::CosyMatrix matrix;
    if (matrix.LoadFromFile(filename)) {
        std::cout << "✓ File loaded successfully!" << std::endl;
        
        const auto& elements = matrix.GetElements();
        if (elements.size() > 0) {
            std::cout << "\nFirst 3 elements:" << std::endl;
            for (size_t i = 0; i < std::min(size_t(3), elements.size()); ++i) {
                std::cout << "  [" << i << "] coeffs=[";
                for (int j = 0; j < 5; ++j) {
                    std::cout << std::setw(10) << elements[i].coefficients[j];
                    if (j < 4) std::cout << ", ";
                }
                std::cout << "] exp=";
                for (int j = 0; j < 5; ++j) {
                    std::cout << elements[i].exponents[j];
                }
                std::cout << " order=" << elements[i].order << std::endl;
            }
        }
    } else {
        std::cout << "❌ Failed to load file!" << std::endl;
    }
}

int main(int argc, char** argv) {
    std::cout << "========================================" << std::endl;
    std::cout << "COSY Matrix Parser Validation" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Test various problematic line formats found in actual files
    std::cout << "\n=== Testing Different Line Formats ===" << std::endl;
    
    // Format 1: HMS/SHMS format with 5 coefficients + 6 digits
    TestSingleLine("  0.0000000E+00 0.0000000E+00 0.0000000E+00 0.0000000E+00  1.000000     000010",
                   "HMS format with spaces before exponents (6 digits)");
    
    // Format 2: Standard 6-digit format
    TestSingleLine("  0.7779354     -3.321846     0.0000000E+00 0.0000000E+00 0.0000000E+00 100000",
                   "Standard 6-digit exponent format");
    
    // Format 3: HRS format with 4 coefficients + 5 digits
    TestSingleLine("  0.777935    -3.32185     0.000000000E+00 0.000000000E+00 10000",
                   "HRS format with 4 coefficients (5 digits)");
    
    // Format 4: Concatenated negative numbers (NO SPACE)
    TestSingleLine("  0.348558022     0.273541506E-01-0.671779120E-02 0.258721514     10000",
                   "Concatenated numbers (E-01-0.671)");
    
    // Format 5: More concatenated negatives
    TestSingleLine("  -3.15053452    -0.249126126     0.611924309E-01 0.918777414E-01 01000",
                   "Multiple negative numbers");
    
    // Format 6: Mixed positive/negative
    TestSingleLine("  0.177375605    -0.859578863     -2.10043200     0.155654165     00010",
                   "Mixed positive and negative");
    
    // Format 7: Very small numbers
    TestSingleLine("  0.0000000E+00 0.0000000E+00 0.0000000E+00 0.132568347E-15 000000",
                   "Very small scientific notation");
    
    // Test actual files from different spectrometers
    std::cout << "\n\n=== Testing Actual Matrix Files ===" << std::endl;
    
    // Allow command-line file specification
    if (argc > 1) {
        std::cout << "\nTesting user-specified file:" << std::endl;
        TestFile(argv[1]);
    } else {
        // Test a representative sample
        std::cout << "\n--- HMS Matrices ---" << std::endl;
        TestFile("../data/matrices/hms/forward_cosy.dat");
        TestFile("../data/matrices/hms/recon_cosy.dat");
        
        std::cout << "\n--- HRS Left Matrices ---" << std::endl;
        TestFile("../data/matrices/hrsl/hrs_forward_cosy.dat");
        TestFile("../data/matrices/hrsl/hrs_recon_cosy.dat");
        
        std::cout << "\n--- SHMS Matrices ---" << std::endl;
        TestFile("../data/matrices/shms/shms_forward.dat");
        TestFile("../data/matrices/shms/shms_recon.dat");
        
        std::cout << "\n--- SOS Matrices ---" << std::endl;
        TestFile("../data/matrices/sos/forward_cosy.dat");
        TestFile("../data/matrices/sos/recon_cosy.dat");
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "✓ Parser validation complete!" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}
