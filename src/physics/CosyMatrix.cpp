// src/physics/CosyMatrix.cpp
// Implementation of COSY matrix operations
// Ported from transp.f (Fortran SIMC)

#include "simc/CosyMatrix.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>

namespace simc {

// ============================================================================
// CosyMatrix Implementation
// ============================================================================

bool CosyMatrix::LoadFromFile(const std::string& filename) {
    // Port from transp.f lines 50-150
    // Read COSY-generated matrix elements from file
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open COSY matrix file: " << filename << std::endl;
        return false;
    }
    
    elements_.clear();
    max_order_ = 0;
    
    std::string line;
    int line_num = 0;
    
    while (std::getline(file, line)) {
        ++line_num;
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '!' || line[0] == '#') {
            continue;
        }
        
        // Try to parse matrix element
        MatrixElement element;
        if (ParseLine(line, element)) {
            elements_.push_back(element);
            max_order_ = std::max(max_order_, element.order);
        } else {
            // Not an error - could be header line
            // Only warn if line looks like data but failed to parse
            if (line.find_first_of("0123456789") != std::string::npos) {
                std::cerr << "Warning: Could not parse line " << line_num 
                         << " in " << filename << std::endl;
            }
        }
    }
    
    if (elements_.empty()) {
        std::cerr << "Error: No matrix elements loaded from " << filename << std::endl;
        return false;
    }
    
    is_loaded_ = true;
    
    std::cout << "Loaded COSY matrix: " << filename << std::endl;
    std::cout << "  Elements: " << elements_.size() << std::endl;
    std::cout << "  Max order: " << max_order_ << std::endl;
    
    return true;
}

bool CosyMatrix::ParseLine(const std::string& line, MatrixElement& element) {
    // Port from transp.f lines 100-130
    // Parse format: output_index coefficient order ex0 ex1 ex2 ex3 ex4
    
    std::istringstream iss(line);
    
    int output_idx_fortran;  // 1-indexed from file
    double coeff;
    int order;
    int ex[5];
    
    // Read values
    if (!(iss >> output_idx_fortran >> coeff >> order >> 
          ex[0] >> ex[1] >> ex[2] >> ex[3] >> ex[4])) {
        return false;
    }
    
    // Validate
    if (output_idx_fortran < 1 || output_idx_fortran > 5) {
        return false;
    }
    
    // Convert to 0-indexed
    element.output_index = output_idx_fortran - 1;
    element.coefficient = coeff;
    element.order = order;
    
    // Copy exponents
    for (int i = 0; i < 5; ++i) {
        element.exponents[i] = ex[i];
    }
    
    // Verify order matches sum of exponents
    int sum = 0;
    for (int i = 0; i < 5; ++i) {
        sum += ex[i];
    }
    
    if (sum != order) {
        std::cerr << "Warning: Order mismatch (stated=" << order 
                  << ", computed=" << sum << ")" << std::endl;
    }
    
    return true;
}

void CosyMatrix::Apply(const double input[5], double output[5]) const {
    // Port from transp.f lines 200-250
    // Apply transfer map: output = M(input)
    
    if (!is_loaded_) {
        std::cerr << "Error: COSY matrix not loaded" << std::endl;
        return;
    }
    
    // Initialize output to zero
    for (int i = 0; i < 5; ++i) {
        output[i] = 0.0;
    }
    
    // Sum contributions from all matrix elements
    for (const auto& element : elements_) {
        int idx = element.output_index;
        output[idx] += element.Evaluate(input);
    }
}

std::vector<CosyMatrix::MatrixElement> CosyMatrix::GetElementsForOutput(int output_index) const {
    std::vector<MatrixElement> result;
    
    for (const auto& elem : elements_) {
        if (elem.output_index == output_index) {
            result.push_back(elem);
        }
    }
    
    return result;
}

void CosyMatrix::Print() const {
    std::cout << "\n=== COSY Matrix ===" << std::endl;
    std::cout << "Loaded: " << (is_loaded_ ? "Yes" : "No") << std::endl;
    std::cout << "Elements: " << elements_.size() << std::endl;
    std::cout << "Max order: " << max_order_ << std::endl;
    
    if (!is_loaded_ || elements_.empty()) {
        return;
    }
    
    // Print summary by output variable
    const char* var_names[5] = {"x", "xp", "y", "yp", "delta"};
    
    for (int i = 0; i < 5; ++i) {
        auto elems = GetElementsForOutput(i);
        std::cout << "\n" << var_names[i] << "_out: " << elems.size() << " terms" << std::endl;
        
        // Print first few terms
        int n_print = std::min(5, (int)elems.size());
        for (int j = 0; j < n_print; ++j) {
            const auto& e = elems[j];
            std::cout << "  " << e.coefficient << " * ";
            for (int k = 0; k < 5; ++k) {
                if (e.exponents[k] > 0) {
                    std::cout << var_names[k];
                    if (e.exponents[k] > 1) {
                        std::cout << "^" << e.exponents[k];
                    }
                    std::cout << " ";
                }
            }
            std::cout << std::endl;
        }
        if (elems.size() > 5) {
            std::cout << "  ... (" << (elems.size() - 5) << " more terms)" << std::endl;
        }
    }
    std::cout << "===================" << std::endl;
}

// ============================================================================
// CosyMatrixSequential Implementation
// ============================================================================

void CosyMatrixSequential::AddElement(const CosyMatrix& matrix,
                                      const std::string& type,
                                      double length) {
    Element elem;
    elem.matrix = matrix;
    elem.type = type;
    elem.length = length;
    elements_.push_back(elem);
}

void CosyMatrixSequential::Apply(const double input[5], double output[5]) const {
    // Copy input to temporary
    double temp_in[5], temp_out[5];
    std::memcpy(temp_in, input, 5 * sizeof(double));
    
    // Apply each element sequentially
    for (const auto& elem : elements_) {
        elem.matrix.Apply(temp_in, temp_out);
        std::memcpy(temp_in, temp_out, 5 * sizeof(double));
    }
    
    // Copy final result to output
    std::memcpy(output, temp_out, 5 * sizeof(double));
}

} // namespace simc
