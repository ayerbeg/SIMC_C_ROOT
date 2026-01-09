// src/physics/CosyMatrix.cpp
// Implementation of COSY matrix operations
// Ported from transp.f (Fortran SIMC)
// CORRECTED VERSION - Fixed format parsing and units

#include "simc/CosyMatrix.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <iomanip>

namespace simc {

// ============================================================================
// CosyMatrix Implementation
// ============================================================================

bool CosyMatrix::LoadFromFile(const std::string& filename) {
    // Port from transp.f lines 50-420
    // Read COSY-generated matrix elements from file
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open COSY matrix file: " << filename << std::endl;
        return false;
    }
    
    elements_.clear();
    max_order_ = 0;
    is_drift_ = true;  // Assume drift until proven otherwise
    drift_distance_ = 0.0;
    
    std::string line;
    int line_num = 0;
    
    while (std::getline(file, line)) {
        ++line_num;
        
        // Skip empty lines
        if (line.empty()) {
            continue;
        }
        
        // Skip comment lines (start with '!')
        if (line[0] == '!' || line[0] == '#') {
            // Check for LENGTH comment (transp.f lines 170-175)
            if (line.find("!LENGTH:") != std::string::npos) {
                // Extract length (in meters, convert to cm)
                size_t pos = line.find("!LENGTH:");
                std::string length_str = line.substr(pos + 8);
                std::istringstream iss(length_str);
                double length_m;
                if (iss >> length_m) {
                    drift_distance_ = length_m * 100.0; // Convert m to cm
                }
            }
            continue;
        }
        
        // Skip separator lines
        if (line.find("---") != std::string::npos) {
            continue;
        }
        
        // Try to parse matrix element
        MatrixElement element;
        if (ParseLine(line, element)) {
            elements_.push_back(element);
            max_order_ = std::max(max_order_, element.order);
        } else {
            // Not necessarily an error - could be header line
            // Only warn if line looks like numerical data but failed to parse
            if (line.find_first_of("0123456789") != std::string::npos &&
                line.find("NAME:") == std::string::npos &&
                line.find("REGION:") == std::string::npos &&
                line.find("OFFSET:") == std::string::npos) {
                std::cerr << "Warning: Could not parse line " << line_num 
                         << " in " << filename << std::endl;
                std::cerr << "  Line: " << line << std::endl;
            }
        }
    }
    
    file.close();
    
    if (elements_.empty()) {
        std::cerr << "Error: No matrix elements loaded from " << filename << std::endl;
        return false;
    }
    
    // Check if this is a pure drift (transp.f lines 340-410)
    CheckIfDrift();
    
    is_loaded_ = true;
    
    std::cout << "Loaded COSY matrix: " << filename << std::endl;
    std::cout << "  Elements: " << elements_.size() << " lines" << std::endl;
    std::cout << "  Max order: " << max_order_ << std::endl;
    if (is_drift_) {
        std::cout << "  Type: Pure drift" << std::endl;
        std::cout << "  Distance: " << drift_distance_ << " cm" << std::endl;
    }
    
    return true;
}

bool CosyMatrix::ParseLine(const std::string& line, MatrixElement& element) {
    // Port from transp.f lines 1200 format and lines 280-320
    // Format: c1 c2 c3 c4 c5  ex0 ex1 ex2 ex3 exTOF exDelta
    // Where:
    //   c1-c5: coefficients for x, xp, y, yp, dL (5 doubles in g14.7 format)
    //   ex0-ex5: exponents for x, xp, y, yp, TOF, delta (6 single digits)
    
    std::istringstream iss(line);
    
    // Read 5 coefficients
    double c[5];
    for (int i = 0; i < 5; ++i) {
        if (!(iss >> c[i])) {
            return false;
        }
        element.coefficients[i] = c[i];
    }
    
    // Read 6 single-digit exponents
    int ex[6];
    for (int i = 0; i < 6; ++i) {
        if (!(iss >> ex[i])) {
            return false;
        }
    }
    
    // Store exponents (skip TOF exponent which is ex[4])
    element.exponents[0] = ex[0];  // x
    element.exponents[1] = ex[1];  // xp
    element.exponents[2] = ex[2];  // y
    element.exponents[3] = ex[3];  // yp
    element.exponents[4] = ex[5];  // delta (NOT ex[4] which is TOF!)
    
    // Calculate total order
    element.order = 0;
    for (int i = 0; i < 5; ++i) {
        element.order += element.exponents[i];
    }
    
    // Check if TOF exponent is used with non-zero coefficients
    // (transp.f lines 310-320)
    if (ex[4] != 0) {
        bool has_nonzero = false;
        for (int i = 0; i < 4; ++i) {  // Check x, xp, y, yp (not dL)
            if (element.coefficients[i] != 0.0) {
                has_nonzero = true;
                break;
            }
        }
        if (has_nonzero) {
            std::cerr << "Warning: Non-zero TOF term found (should not happen)" << std::endl;
        }
        // Ignore this element by returning false
        return false;
    }
    
    return true;
}

void CosyMatrix::CheckIfDrift() {
    // Port from transp.f lines 340-410
    // Check if this transformation is a simple field-free drift
    // 
    // For a drift, first-order terms must be:
    //   x_out  = x_in          (coefficient = 1.0)
    //   xp_out = xp_in         (coefficient = 1.0)
    //   y_out  = y_in          (coefficient = 1.0)
    //   yp_out = yp_in         (coefficient = 1.0)
    //   x_out  = distance * xp_in  (where distance is in meters in file, so coeff = distance*1000 in mm)
    //   y_out  = distance * yp_in
    // 
    // All higher-order terms must be zero
    
    const double coeff_min = 1.0e-14;  // Minimum coefficient to consider non-zero
    
    is_drift_ = true;
    double computed_drift = 0.0;
    
    for (const auto& elem : elements_) {
        int order = elem.order;
        int ex_x  = elem.exponents[0];
        int ex_xp = elem.exponents[1];
        int ex_y  = elem.exponents[2];
        int ex_yp = elem.exponents[3];
        
        double c_x  = elem.coefficients[0];
        double c_xp = elem.coefficients[1];
        double c_y  = elem.coefficients[2];
        double c_yp = elem.coefficients[3];
        double csum = std::abs(c_x) + std::abs(c_xp) + std::abs(c_y) + std::abs(c_yp);
        
        if (order == 1) {
            // First-order terms
            if (ex_x == 1 && ex_xp == 0 && ex_y == 0 && ex_yp == 0) {
                // x term: should be x_out = x_in (c_x = 1.0, others = 0)
                if (std::abs(c_x - 1.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_xp) > coeff_min) is_drift_ = false;
                if (std::abs(c_y) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp) > coeff_min) is_drift_ = false;
            }
            else if (ex_x == 0 && ex_xp == 1 && ex_y == 0 && ex_yp == 0) {
                // xp term: x_out = distance * xp_in, xp_out = xp_in
                computed_drift = c_x * 1000.0;  // Convert from m to cm (coeff is in m, mrad->mrad)
                if (std::abs(c_xp - 1.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_y) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp) > coeff_min) is_drift_ = false;
            }
            else if (ex_x == 0 && ex_xp == 0 && ex_y == 1 && ex_yp == 0) {
                // y term: should be y_out = y_in (c_y = 1.0, others = 0)
                if (std::abs(c_x) > coeff_min) is_drift_ = false;
                if (std::abs(c_xp) > coeff_min) is_drift_ = false;
                if (std::abs(c_y - 1.0) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp) > coeff_min) is_drift_ = false;
            }
            else if (ex_x == 0 && ex_xp == 0 && ex_y == 0 && ex_yp == 1) {
                // yp term: y_out = distance * yp_in, yp_out = yp_in
                if (std::abs(computed_drift - c_y * 1000.0) > 0.01) is_drift_ = false;
                if (std::abs(c_x) > coeff_min) is_drift_ = false;
                if (std::abs(c_xp) > coeff_min) is_drift_ = false;
                if (std::abs(c_yp - 1.0) > coeff_min) is_drift_ = false;
            }
        } else {
            // Higher-order terms: all should be zero for a drift
            if (csum > coeff_min) {
                is_drift_ = false;
            }
        }
    }
    
    if (is_drift_) {
        // Verify computed drift matches file comment
        if (drift_distance_ > 0.01 && 
            std::abs(drift_distance_ - computed_drift) > 0.01) {
            std::cerr << "Warning: Drift distance mismatch" << std::endl;
            std::cerr << "  From file comment: " << drift_distance_ << " cm" << std::endl;
            std::cerr << "  From matrix elements: " << computed_drift << " cm" << std::endl;
        }
        if (computed_drift > 0.01) {
            drift_distance_ = computed_drift;
        }
    }
}

void CosyMatrix::Apply(const double input[5], double output[5]) const {
    // Port from transp.f lines 230-260
    // Apply transfer map: output = M(input)
    // 
    // CRITICAL: Input and output are in COSY units:
    //   [x(cm), xp(mrad), y(cm), yp(mrad), delta(%)]
    
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
        // Each element contributes to all 5 outputs
        for (int out_idx = 0; out_idx < 5; ++out_idx) {
            output[out_idx] += element.Evaluate(out_idx, input);
        }
    }
}

void CosyMatrix::Print() const {
    std::cout << "\n=== COSY Matrix ===" << std::endl;
    std::cout << "Loaded: " << (is_loaded_ ? "Yes" : "No") << std::endl;
    std::cout << "Elements: " << elements_.size() << " lines" << std::endl;
    std::cout << "Max order: " << max_order_ << std::endl;
    
    if (is_drift_) {
        std::cout << "Type: Pure drift" << std::endl;
        std::cout << "Distance: " << drift_distance_ << " cm" << std::endl;
    } else {
        std::cout << "Type: Magnet transformation" << std::endl;
    }
    
    if (!is_loaded_ || elements_.empty()) {
        return;
    }
    
    // Print summary - count non-zero terms for each output
    const char* var_names[5] = {"x", "xp", "y", "yp", "dL"};
    int counts[5] = {0, 0, 0, 0, 0};
    
    for (const auto& elem : elements_) {
        for (int i = 0; i < 5; ++i) {
            if (elem.coefficients[i] != 0.0) {
                counts[i]++;
            }
        }
    }
    
    std::cout << "\nNon-zero terms per output:" << std::endl;
    for (int i = 0; i < 5; ++i) {
        std::cout << "  " << var_names[i] << ": " << counts[i] << " terms" << std::endl;
    }
    
    // Print first few elements (all 5 coefficients)
    std::cout << "\nFirst few matrix elements:" << std::endl;
    int n_print = std::min(5, (int)elements_.size());
    for (int i = 0; i < n_print; ++i) {
        const auto& e = elements_[i];
        
        // Print exponents
        std::cout << "  Exponents: [";
        for (int j = 0; j < 5; ++j) {
            std::cout << e.exponents[j];
            if (j < 4) std::cout << ",";
        }
        std::cout << "] (order " << e.order << ")" << std::endl;
        
        // Print coefficients
        std::cout << "    Coefficients: [" << std::scientific << std::setprecision(6);
        for (int j = 0; j < 5; ++j) {
            std::cout << e.coefficients[j];
            if (j < 4) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
    if (elements_.size() > 5) {
        std::cout << "  ... (" << (elements_.size() - 5) << " more elements)" << std::endl;
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
    elem.matrix = std::move(matrix);
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
        
        // temp_out[4] contains dL (path length correction in cm)
        // Real path length = nominal - dL (from transp.f line 260)
        double delta_z = -temp_out[4];
        
        // Update total path length
        // Nominal length is elem.length (in meters), convert to cm
        double nominal_length = elem.length * 100.0;  // m to cm
        total_path_length += (nominal_length + delta_z);
        
        // Copy output to input for next iteration
        // Note: delta (temp_out[4]) is unchanged by transformation
        temp_in[0] = temp_out[0];  // x
        temp_in[1] = temp_out[1];  // xp
        temp_in[2] = temp_out[2];  // y
        temp_in[3] = temp_out[3];  // yp
        // temp_in[4] = temp_in[4];  // delta unchanged
    }
    
    // Copy final result to output
    std::memcpy(output, temp_in, 5 * sizeof(double));
}

std::string CosyMatrixSequential::GetElementName(size_t index) const {
    if (index < elements_.size()) {
        return elements_[index].name;
    }
    return "";
}

void CosyMatrixSequential::Print() const {
    std::cout << "\n=== Sequential COSY Transport ===" << std::endl;
    std::cout << "Number of elements: " << elements_.size() << std::endl;
    
    for (size_t i = 0; i < elements_.size(); ++i) {
        const auto& elem = elements_[i];
        std::cout << "\n[" << i << "] " << elem.name << std::endl;
        std::cout << "    Length: " << elem.length << " m" << std::endl;
        std::cout << "    Matrix elements: " << elem.matrix.GetNumElements() << std::endl;
        std::cout << "    Max order: " << elem.matrix.GetMaxOrder() << std::endl;
        if (elem.matrix.IsDrift()) {
            std::cout << "    Type: Drift (" << elem.matrix.GetDriftDistance() << " cm)" << std::endl;
        }
    }
    std::cout << "====================================" << std::endl;
}

} // namespace simc
