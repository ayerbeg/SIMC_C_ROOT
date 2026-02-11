// src/physics/CosyMatrix.cpp
// Implementation of COSY matrix operations
// Ported from transp.f (Fortran SIMC)
// ROBUST VERSION - Handles all Fortran formatting variations

#include "simc/physics/CosyMatrix.h"
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
    
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open COSY matrix file: " << filename << std::endl;
        return false;
    }
    
    elements_.clear();
    max_order_ = 0;
    
    // Read and parse the file
    std::string line;
    int line_num = 0;
    int parse_errors = 0;
    int comment_lines = 0;
    int separator_lines = 0;
    const int max_errors_to_show = 3;
    
    while (std::getline(infile, line)) {
        line_num++;
        
        // Check for comments first
        std::string trimmed = line;
        trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
        if (trimmed.empty() || trimmed[0] == '!' || trimmed[0] == '#') {
            comment_lines++;
            continue;
        }
        
        // Check for separator lines
        bool is_separator = true;
        for (char c : trimmed) {
            if (c != '-' && c != ' ' && c != '\t' && c != '\r' && c != '\n') {
                is_separator = false;
                break;
            }
        }
        if (is_separator) {
            separator_lines++;
            continue;
        }
        
        std::vector<double> coeffs;
        std::vector<int> exponents;
        
        if (ParseMatrixLine(line, coeffs, exponents)) {
            MatrixElement elem;
            for (size_t i = 0; i < coeffs.size() && i < 5; ++i) {
                elem.coefficients[i] = coeffs[i];
            }
            for (size_t i = 0; i < exponents.size() && i < 5; ++i) {
                elem.exponents[i] = exponents[i];
            }
            
            // Calculate order (sum of exponents)
            int order = 0;
            for (int exp : exponents) order += exp;
            elem.order = order;
            if (order > max_order_) max_order_ = order;
            
            elements_.push_back(elem);
        } else {
            parse_errors++;
            if (parse_errors <= max_errors_to_show) {
                std::cerr << "Warning: Could not parse line " << line_num 
                         << " in " << filename << std::endl;
                std::cerr << "Line: " << line << std::endl;
            } else if (parse_errors == max_errors_to_show + 1) {
                std::cerr << "... (suppressing further parse warnings)" << std::endl;
            }
        }
    }
    
    infile.close();
    
    if (elements_.empty()) {
        std::cerr << "Error: No valid matrix elements found in " << filename << std::endl;
        std::cerr << "  Lines read: " << line_num << std::endl;
        std::cerr << "  Comments: " << comment_lines << std::endl;
        std::cerr << "  Separators: " << separator_lines << std::endl;
        std::cerr << "  Parse errors: " << parse_errors << std::endl;
        return false;
    }
    
    is_loaded_ = true;
    
    // Check if this is a drift transformation
    CheckIfDrift();
    
    std::cout << "Loaded COSY matrix: " << filename << std::endl;
    std::cout << "  Total lines: " << line_num << std::endl;
    std::cout << "  Comments: " << comment_lines << std::endl;
    std::cout << "  Separators: " << separator_lines << std::endl;
    std::cout << "  Data elements: " << elements_.size() << std::endl;
    if (parse_errors > 0) {
        std::cout << "  Parse warnings: " << parse_errors << " (may be normal for some formats)" << std::endl;
    }
    std::cout << "  Max order: " << max_order_ << std::endl;
    
    return true;
}

bool CosyMatrix::ParseMatrixLine(const std::string& line, 
                                  std::vector<double>& coeffs,
                                  std::vector<int>& exponents) {
    coeffs.clear();
    exponents.clear();
    
    // Skip empty lines and comments
    std::string trimmed = line;
    trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
    trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);
    
    if (trimmed.empty() || trimmed[0] == '!' || trimmed[0] == '#') {
        return false;
    }
    
    // Check for separator lines
    bool is_separator = true;
    for (char c : trimmed) {
        if (c != '-' && c != ' ' && c != '\t') {
            is_separator = false;
            break;
        }
    }
    if (is_separator) return false;
    
    // Find rightmost group of 5-6 consecutive digits for exponents
    size_t exp_start = std::string::npos;
    size_t exp_end = std::string::npos;
    
    for (size_t i = trimmed.length(); i > 0; --i) {
        size_t pos = i - 1;
        if (std::isdigit(trimmed[pos])) {
            if (exp_end == std::string::npos) {
                exp_end = pos;
            }
            exp_start = pos;
        } else if (exp_end != std::string::npos) {
            size_t len = exp_end - exp_start + 1;
            if (len >= 5 && len <= 6) {
                break;
            } else {
                exp_start = std::string::npos;
                exp_end = std::string::npos;
            }
        }
    }
    
    if (exp_start == std::string::npos || exp_end == std::string::npos) {
        return false;
    }
    
    size_t exp_len = exp_end - exp_start + 1;
    if (exp_len < 5 || exp_len > 6) {
        return false;
    }
    
    // Extract exactly 5 exponent digits (ignore 6th if present - that's TOF)
    // Format: [x, xp, y, yp, delta] or [x, xp, y, yp, TOF, delta]
    if (exp_len == 5) {
        // 5 digits: direct mapping
        for (size_t i = 0; i < 5; ++i) {
            if (std::isdigit(trimmed[exp_start + i])) {
                exponents.push_back(trimmed[exp_start + i] - '0');
            } else {
                return false;
            }
        }
    } else {
        // 6 digits: positions 0,1,2,3,5 (skip position 4 which is TOF)
        for (size_t i = 0; i < 6; ++i) {
            if (i == 4) continue;  // Skip TOF exponent
            if (std::isdigit(trimmed[exp_start + i])) {
                exponents.push_back(trimmed[exp_start + i] - '0');
            } else {
                return false;
            }
        }
    }
    
    if (exponents.size() != 5) return false;
    
    // Parse coefficients from the part before exponents
    // Need special handling for Fortran scientific notation without spaces
    std::string coeff_part = trimmed.substr(0, exp_start);
    
    // Manual tokenization to handle concatenated numbers like "E-01-0.67"
    size_t pos = 0;
    while (pos < coeff_part.length()) {
        // Skip whitespace
        while (pos < coeff_part.length() && std::isspace(coeff_part[pos])) {
            pos++;
        }
        if (pos >= coeff_part.length()) break;
        
        // Extract a number
        size_t start = pos;
        bool has_decimal = false;
        bool has_exp = false;
        
        // Handle leading sign
        if (coeff_part[pos] == '+' || coeff_part[pos] == '-') {
            pos++;
        }
        
        // Read digits and decimal point
        while (pos < coeff_part.length()) {
            char c = coeff_part[pos];
            if (std::isdigit(c)) {
                pos++;
            } else if (c == '.' && !has_decimal && !has_exp) {
                has_decimal = true;
                pos++;
            } else if ((c == 'E' || c == 'e' || c == 'D' || c == 'd') && !has_exp) {
                has_exp = true;
                pos++;
                // Handle exponent sign
                if (pos < coeff_part.length() && 
                    (coeff_part[pos] == '+' || coeff_part[pos] == '-')) {
                    pos++;
                }
            } else {
                break;
            }
        }
        
        if (pos > start) {
            std::string num_str = coeff_part.substr(start, pos - start);
            // Replace D with E for stod
            for (char& c : num_str) {
                if (c == 'D' || c == 'd') c = 'E';
            }
            try {
                double val = std::stod(num_str);
                coeffs.push_back(val);
            } catch (...) {
                return false;
            }
        }
    }
    
    // Should have 4 or 5 coefficients
    if (coeffs.size() < 4 || coeffs.size() > 5) {
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
    for (int i = 0; i < 5; ++i) {
        output[i] = 0.0;
    }
    
    // Sum contributions from all matrix elements
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
    
    const double coeff_min = 1.0e-14;
    
    is_drift_ = true;
    drift_distance_ = 0.0;
    
    for (const auto& elem : elements_) {
        int order = elem.order;
        
        if (order == 1) {
            // Get exponents
            int e_x  = elem.exponents[0];
            int e_xp = elem.exponents[1];
            int e_y  = elem.exponents[2];
            int e_yp = elem.exponents[3];
            int e_d  = elem.exponents[4];
            
            // Get coefficients
            double c_x  = elem.coefficients[0];
            double c_xp = elem.coefficients[1];
            double c_y  = elem.coefficients[2];
            double c_yp = elem.coefficients[3];
            
            // Check each first-order term
            if (e_x == 1 && e_xp == 0 && e_y == 0 && e_yp == 0 && e_d == 0) {
                if (std::abs(c_x - 1.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_xp - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_y - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp - 0.0) > coeff_min) is_drift_ = false;
            }
            else if (e_x == 0 && e_xp == 1 && e_y == 0 && e_yp == 0 && e_d == 0) {
                drift_distance_ = c_x * 1000.0;  // Convert from m to cm
                if (std::abs(c_xp - 1.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_y - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp - 0.0) > coeff_min) is_drift_ = false;
            }
            else if (e_x == 0 && e_xp == 0 && e_y == 1 && e_yp == 0 && e_d == 0) {
                if (std::abs(c_x - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_xp - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_y - 1.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp - 0.0) > coeff_min) is_drift_ = false;
            }
            else if (e_x == 0 && e_xp == 0 && e_y == 0 && e_yp == 1 && e_d == 0) {
                if (std::abs(c_x - 0.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_xp - 0.0) > coeff_min) is_drift_ = false;
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
                    std::cout << "1";
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
        
        // Accumulate path length from 5th component
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

const std::vector<CosyMatrix::MatrixElement>& CosyMatrix::GetElements() const {
    return elements_;
}

} // namespace simc
