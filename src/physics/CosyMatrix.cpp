// src/physics/CosyMatrix.cpp
// Implementation of COSY matrix operations
// Ported from transp.f (Fortran SIMC)
// CORRECTED VERSION - Proper Fortran format parsing

#include "simc/CosyMatrix.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cctype>

namespace simc {

// ============================================================================
// CosyMatrix Implementation
// ============================================================================

bool CosyMatrix::LoadFromFile(const std::string& filename) {
    // Port from transp.f lines 340-530
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
        
        // Skip empty lines
        if (line.empty()) {
            continue;
        }
        
        // Skip comment lines (starting with '!')
        if (line[0] == '!') {
            continue;
        }
        
        // Skip separator lines (lines starting with space and dashes)
        // Format: " --------------..."
        if (line.size() > 1 && line[0] == ' ' && line[1] == '-') {
            continue;
        }
        
        // Try to parse matrix element
        MatrixElement element;
        if (ParseLine(line, element)) {
            elements_.push_back(element);
            max_order_ = std::max(max_order_, element.order);
        } else {
            // Only warn if line contains digits (likely data line that failed)
            if (line.find_first_of("0123456789") != std::string::npos &&
                line[0] != '!') {
                std::cerr << "Warning: Could not parse line " << line_num 
                         << " in " << filename << std::endl;
                std::cerr << "Line: " << line << std::endl;
            }
        }
    }
    
    if (elements_.empty()) {
        std::cerr << "Error: No matrix elements loaded from " << filename << std::endl;
        return false;
    }
    
    is_loaded_ = true;
    
    // Check if this is a drift transformation
    CheckIfDrift();
    
    std::cout << "Loaded COSY matrix: " << filename << std::endl;
    std::cout << "  Elements: " << elements_.size() << std::endl;
    std::cout << "  Max order: " << max_order_ << std::endl;
    if (is_drift_) {
        std::cout << "  Type: Drift (length = " << drift_distance_ << " cm)" << std::endl;
    }
    
    return true;
}

bool CosyMatrix::ParseLine(const std::string& line, MatrixElement& element) {
    // Port from transp.f lines 419-426
    // 
    // Fortran format: 1200 format(1x,5g14.7,1x,6i1)
    // This means:
    //   1x        = skip 1 character
    //   5g14.7    = read 5 floating-point numbers (general format, 14 chars each)
    //   1x        = skip 1 character
    //   6i1       = read 6 ONE-DIGIT integers (NOT a 6-digit string!)
    //
    // Example line:
    //    1.000000     0.0000000E+00 0.0000000E+00 0.0000000E+00 0.0000000E+00 100000
    //    ^          ^                                                          ^
    //    skip 1     5 coefficients                                             6 digits
    //
    // The key insight: Fortran's format "6i1" reads 6 consecutive single digits
    // without whitespace, so "100000" is read as [1,0,0,0,0,0]
    
    std::istringstream iss(line);
    
    // Read 5 coefficients (they're space-separated in the file)
    for (int i = 0; i < 5; ++i) {
        if (!(iss >> element.coefficients[i])) {
            return false;
        }
    }
    
    // Read the 6-character exponent field
    // In the actual files, this appears as a 6-digit number like "100000"
    // But we need to parse it as 6 individual digits
    std::string exp_str;
    if (!(iss >> exp_str)) {
        return false;
    }
    
    // Exponents should be exactly 6 characters
    if (exp_str.length() != 6) {
        return false;
    }
    
    // Verify all characters are digits
    for (int i = 0; i < 6; ++i) {
        if (!std::isdigit(static_cast<unsigned char>(exp_str[i]))) {
            std::cerr << "Error: Non-digit in exponent string: " << exp_str << std::endl;
            return false;
        }
    }
    
    // Parse each digit as an exponent
    // File format: [x_exp, xp_exp, y_exp, yp_exp, TOF_exp, delta_exp]
    // Our storage: [x, xp, y, yp, delta]  (we ignore TOF at position 4)
    
    element.exponents[0] = exp_str[0] - '0';  // x
    element.exponents[1] = exp_str[1] - '0';  // xp
    element.exponents[2] = exp_str[2] - '0';  // y
    element.exponents[3] = exp_str[3] - '0';  // yp
    // Skip exp_str[4] - that's TOF, should always be 0
    element.exponents[4] = exp_str[5] - '0';  // delta
    
    // Calculate total order (sum of exponents, excluding TOF)
    element.order = 0;
    for (int i = 0; i < 5; ++i) {
        element.order += element.exponents[i];
    }
    
    // Verify TOF exponent is zero (as Fortran code requires)
    int tof_exp = exp_str[4] - '0';
    if (tof_exp != 0) {
        // Check if any coefficient is non-zero
        bool has_nonzero = false;
        for (int i = 0; i < 4; ++i) {  // Check first 4 coefficients (x, xp, y, yp)
            if (element.coefficients[i] != 0.0) {
                has_nonzero = true;
                break;
            }
        }
        
        if (has_nonzero) {
            std::cerr << "Error: Non-zero TOF term with non-zero coefficients!" << std::endl;
            return false;
        }
        
        // Fortran code ignores this term - we'll return false to skip it
        return false;
    }
    
    return true;
}

void CosyMatrix::Apply(const double input[5], double output[5]) const {
    // Port from transp.f lines 248-265
    // Apply transfer map: output = M(input)
    
    if (!is_loaded_) {
        std::cerr << "Error: COSY matrix not loaded" << std::endl;
        return;
    }
    
    // Initialize output to zero
    // From transp.f lines 248-250:
    //   do i = 1,5
    //      sum(i) = 0.
    //   enddo
    for (int i = 0; i < 5; ++i) {
        output[i] = 0.0;
    }
    
    // Sum contributions from all matrix elements
    // From transp.f lines 254-264:
    //   do i = 1,n_terms(spectr,k)
    //     term = 1.0e0
    //     do j = 1,5
    //       temp = 1.0e0
    //       if (expon(spectr,j,i,k).ne.0.) temp = ray(j)**expon(spectr,j,i,k)
    //       term = term*temp
    //     enddo
    //     sum(1) = sum(1) + term*coeff(spectr,1,i,k)
    //     ...
    //   enddo
    
    for (const auto& element : elements_) {
        // Calculate the product of input variables raised to their exponents
        double term = 1.0;
        for (int j = 0; j < 5; ++j) {
            if (element.exponents[j] != 0) {
                term *= std::pow(input[j], element.exponents[j]);
            }
        }
        
        // Add contribution to each output component
        for (int i = 0; i < 5; ++i) {
            output[i] += term * element.coefficients[i];
        }
    }
}

void CosyMatrix::CheckIfDrift() {
    // Port from transp.f lines 463-493
    // Check if this matrix represents a field-free drift
    //
    // A pure drift transformation has only first-order terms with:
    //   x_out  = x_in  + distance * xp_in    (coeff[x][x]=1, coeff[x][xp]=dist)
    //   xp_out = xp_in                        (coeff[xp][xp]=1)
    //   y_out  = y_in  + distance * yp_in    (coeff[y][y]=1, coeff[y][yp]=dist)
    //   yp_out = yp_in                        (coeff[yp][yp]=1)
    //   All other coefficients = 0
    
    const double coeff_min = 1.0e-14;
    
    is_drift_ = true;
    drift_distance_ = 0.0;
    
    for (const auto& elem : elements_) {
        int order = elem.order;
        
        // From transp.f line 469: if (order.eq.1) then
        if (order == 1) {
            // Get exponents
            int e_x  = elem.exponents[0];  // x
            int e_xp = elem.exponents[1];  // xp (theta)
            int e_y  = elem.exponents[2];  // y
            int e_yp = elem.exponents[3];  // yp (phi)
            int e_d  = elem.exponents[4];  // delta
            
            // Get coefficients
            double c_x  = elem.coefficients[0];  // x output
            double c_xp = elem.coefficients[1];  // xp output
            double c_y  = elem.coefficients[2];  // y output
            double c_yp = elem.coefficients[3];  // yp output
            
            // Check each first-order term
            if (e_x == 1 && e_xp == 0 && e_y == 0 && e_yp == 0 && e_d == 0) {
                // x term: x_out coeff should be 1
                if (std::abs(c_x - 1.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_xp - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_y - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp - 0.0) > coeff_min) is_drift_ = false;
            }
            else if (e_x == 0 && e_xp == 1 && e_y == 0 && e_yp == 0 && e_d == 0) {
                // xp term
                drift_distance_ = c_x * 1000.0;  // Convert from m to cm
                if (std::abs(c_xp - 1.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_y - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp - 0.0) > coeff_min) is_drift_ = false;
            }
            else if (e_x == 0 && e_xp == 0 && e_y == 1 && e_yp == 0 && e_d == 0) {
                // y term
                if (std::abs(c_x - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_xp - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_y - 1.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp - 0.0) > coeff_min) is_drift_ = false;
            }
            else if (e_x == 0 && e_xp == 0 && e_y == 0 && e_yp == 1 && e_d == 0) {
                // yp term
                if (std::abs(c_x - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_xp - 0.0) > coeff_min) is_drift_ = false;
                // c_y should equal drift_distance (in meters/1000)
                if (std::abs(drift_distance_ - c_y * 1000.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp - 1.0) > coeff_min) is_drift_ = false;
            }
        } else {
            // Higher-order terms: all coefficients should be zero
            double csum = std::abs(elem.coefficients[0]) + 
                         std::abs(elem.coefficients[1]) +
                         std::abs(elem.coefficients[2]) + 
                         std::abs(elem.coefficients[3]);
            if (csum > coeff_min) {
                is_drift_ = false;
            }
        }
    }
}

void CosyMatrix::Print() const {
    std::cout << "\n=== COSY Matrix ===" << std::endl;
    std::cout << "Loaded: " << (is_loaded_ ? "Yes" : "No") << std::endl;
    std::cout << "Elements: " << elements_.size() << std::endl;
    std::cout << "Max order: " << max_order_ << std::endl;
    
    if (is_drift_) {
        std::cout << "Type: Drift (distance = " << drift_distance_ << " cm)" << std::endl;
    }
    
    if (!is_loaded_ || elements_.empty()) {
        return;
    }
    
    // Print summary by output variable
    const char* var_names[5] = {"x", "xp", "y", "yp", "dL"};
    
    // Count non-zero coefficients for each output
    int count[5] = {0, 0, 0, 0, 0};
    for (const auto& elem : elements_) {
        for (int i = 0; i < 5; ++i) {
            if (elem.coefficients[i] != 0.0) {
                count[i]++;
            }
        }
    }
    
    for (int i = 0; i < 5; ++i) {
        std::cout << "\n" << var_names[i] << "_out: " << count[i] << " terms" << std::endl;
        
        // Print first few non-zero terms
        int n_printed = 0;
        for (const auto& elem : elements_) {
            if (elem.coefficients[i] != 0.0 && n_printed < 5) {
                std::cout << "  " << elem.coefficients[i] << " * ";
                bool first = true;
                for (int k = 0; k < 5; ++k) {
                    if (elem.exponents[k] > 0) {
                        if (!first) std::cout << "*";
                        std::cout << var_names[k];
                        if (elem.exponents[k] > 1) {
                            std::cout << "^" << elem.exponents[k];
                        }
                        first = false;
                    }
                }
                if (first) {
                    std::cout << "1";  // Constant term
                }
                std::cout << std::endl;
                n_printed++;
            }
        }
        if (count[i] > 5) {
            std::cout << "  ... (" << (count[i] - 5) << " more terms)" << std::endl;
        }
    }
    std::cout << "===================" << std::endl;
}

// ============================================================================
// CosyMatrixSequential Implementation
// ============================================================================

void CosyMatrixSequential::AddElement(CosyMatrix matrix,
                                      const std::string& name,
                                      double length) {
    Element elem;
    elem.matrix = matrix;
    elem.name = name;
    elem.length = length;
    elements_.push_back(elem);
}

void CosyMatrixSequential::Apply(const double input[5], double output[5], 
                                 double& total_path_length) const {
    // Copy input to temporary
    double temp_in[5], temp_out[5];
    std::memcpy(temp_in, input, 5 * sizeof(double));
    
    total_path_length = 0.0;
    
    // Apply each element sequentially
    for (const auto& elem : elements_) {
        elem.matrix.Apply(temp_in, temp_out);
        
        // Accumulate path length from 5th component (dL)
        // From transp.f line 279: delta_z = -sum(5)
        // From transp.f line 331: pathlen = pathlen + (zd + delta_z)
        total_path_length += temp_out[4];
        
        // Copy output to input for next element
        std::memcpy(temp_in, temp_out, 5 * sizeof(double));
    }
    
    // Copy final result to output
    std::memcpy(output, temp_out, 5 * sizeof(double));
}

std::string CosyMatrixSequential::GetElementName(size_t index) const {
    if (index >= elements_.size()) {
        return "";
    }
    return elements_[index].name;
}

void CosyMatrixSequential::Print() const {
    std::cout << "\n=== COSY Sequential Matrix ===" << std::endl;
    std::cout << "Number of elements: " << elements_.size() << std::endl;
    
    for (size_t i = 0; i < elements_.size(); ++i) {
        const auto& elem = elements_[i];
        std::cout << "\nElement " << (i+1) << ": " << elem.name << std::endl;
        std::cout << "  Length: " << elem.length << " m" << std::endl;
        std::cout << "  Matrix elements: " << elem.matrix.GetNumElements() << std::endl;
        std::cout << "  Max order: " << elem.matrix.GetMaxOrder() << std::endl;
        if (elem.matrix.IsDrift()) {
            std::cout << "  Type: Drift (" << elem.matrix.GetDriftDistance() << " cm)" << std::endl;
        }
    }
    std::cout << "===============================" << std::endl;
}

} // namespace simc
