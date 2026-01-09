// include/simc/CosyMatrix.h
// COSY matrix element representation and evaluation
// Ported from transp.f (COSY matrix parsing and evaluation)
// CORRECTED VERSION - Fixed format parsing and units

#ifndef SIMC_COSY_MATRIX_H
#define SIMC_COSY_MATRIX_H

#include <vector>
#include <string>
#include <map>
#include <cmath>

namespace simc {

/**
 * @class CosyMatrix
 * @brief Represents and evaluates COSY-generated transport matrices
 * 
 * COSY (COSmic raY) is a differential algebra code used to calculate
 * particle optics. It outputs transfer maps as polynomials.
 * 
 * CRITICAL FORMAT (from actual COSY files):
 * Each line contains:
 * - 5 coefficients: [c_x, c_xp, c_y, c_yp, c_dL] (scientific notation)
 * - 1 six-digit exponent string: "exponent_x exponent_xp exponent_y exponent_yp exponent_TOF exponent_delta"
 *   concatenated as a 6-digit number (e.g., "100000" means [1,0,0,0,0,0])
 * 
 * Example line from forward_cosy.dat:
 *   1.000000     0.0000000E+00 0.0000000E+00 0.0000000E+00 0.0000000E+00 100000
 * 
 * This means (parsing "100000" as [1,0,0,0,0,0]):
 *   x_out   += 1.000000  * x^1 * xp^0 * y^0 * yp^0 * delta^0
 *   xp_out  += 0.0       * x^1 (no contribution)
 *   y_out   += 0.0       * x^1 (no contribution)
 *   yp_out  += 0.0       * x^1 (no contribution)
 *   dL_out  += 0.0       * x^1 (no contribution)
 * 
 * Units (COSY-7):
 *   Input/Output positions: cm
 *   Input/Output angles: mrad (milliradians)
 *   Delta: percent (%)
 * 
 * Vector indices (0-indexed in code):
 *   0: x      - horizontal position (cm)
 *   1: xp     - dx/dz horizontal angle (mrad)
 *   2: y      - vertical position (cm)
 *   3: yp     - dy/dz vertical angle (mrad)
 *   4: delta  - momentum deviation (%)
 * 
 * Exponent order in 6-digit string:
 *   Digit 0: x exponent
 *   Digit 1: xp exponent
 *   Digit 2: y exponent
 *   Digit 3: yp exponent
 *   Digit 4: TOF exponent (should always be 0, ignored)
 *   Digit 5: delta exponent
 * 
 * Based on transp.f from Fortran SIMC
 */
class CosyMatrix {
public:
    /**
     * @struct MatrixElement
     * @brief Single term in the transfer map polynomial
     * 
     * CRITICAL: Unlike typical matrix elements, each COSY line
     * contains coefficients for ALL 5 outputs with a single set of input exponents
     */
    struct MatrixElement {
        double coefficients[5]{0,0,0,0,0}; ///< Coefficients for [x, xp, y, yp, dL] outputs
        int exponents[5]{0,0,0,0,0};       ///< Exponents for [x, xp, y, yp, delta] inputs
        int order{0};                      ///< Total order (sum of exponents)
        
        /**
         * @brief Evaluate this term for a specific output
         * @param output_index Which output (0=x, 1=xp, 2=y, 3=yp, 4=dL)
         * @param input Input vector [x, xp, y, yp, delta] in COSY units
         * @return Contribution to output
         */
        double Evaluate(int output_index, const double input[5]) const {
            double result = coefficients[output_index];
            
            // If coefficient is zero, skip calculation
            if (result == 0.0) return 0.0;
            
            // Compute product: coeff * x^ex[0] * xp^ex[1] * y^ex[2] * yp^ex[3] * delta^ex[4]
            for (int i = 0; i < 5; ++i) {
                if (exponents[i] > 0) {
                    result *= std::pow(input[i], exponents[i]);
                }
            }
            
            return result;
        }
    };
    
    /**
     * @brief Constructor
     */
    CosyMatrix() = default;
    
    /**
     * @brief Load COSY matrix from file
     * @param filename Path to COSY matrix file (.dat)
     * @return true if successful
     * 
     * File format (from actual files in data/matrices/):
     * - Comment lines start with '!'
     * - Separator lines: ' ---...'
     * - Data lines: c1 c2 c3 c4 c5  XXXXXX
     *   where c1-c5 are 5 coefficients (scientific notation)
     *   and XXXXXX is a 6-digit string of single-digit exponents
     * 
     * Example:
     *   0.7779354  -3.321846  0.0  0.0  0.0  100000
     *   ^^^^^^^^^  ^^^^^^^^^  ^^^  ^^^  ^^^  ^^^^^^
     *   c_x        c_xp       c_y  c_yp c_dL exponents
     * 
     * Exponents "100000" = [1,0,0,0,0,0] means x^1 * xp^0 * y^0 * yp^0 * TOF^0 * delta^0
     */
    bool LoadFromFile(const std::string& filename);
    
    /**
     * @brief Apply matrix to input vector
     * @param input Input vector [x(cm), xp(mrad), y(cm), yp(mrad), delta(%)]
     * @param output Output vector [x(cm), xp(mrad), y(cm), yp(mrad), dL(cm)]
     * 
     * From transp.f lines 230-250:
     * For each output variable, sum all matrix elements:
     *   output[i] = sum_over_terms( coefficients[i] * product(input[j]^exponent[j]) )
     * 
     * Note: dL (path length correction) is the 5th output
     */
    void Apply(const double input[5], double output[5]) const;
    
    /**
     * @brief Get maximum order in matrix
     */
    int GetMaxOrder() const { return max_order_; }
    
    /**
     * @brief Get number of matrix elements (lines)
     */
    size_t GetNumElements() const { return elements_.size(); }
    
    /**
     * @brief Check if matrix is loaded
     */
    bool IsLoaded() const { return is_loaded_; }
    
    /**
     * @brief Get drift distance (if this is a pure drift transformation)
     * @return Drift distance in cm, or 0 if not a drift
     */
    double GetDriftDistance() const { return drift_distance_; }
    
    /**
     * @brief Check if this is a pure drift transformation
     */
    bool IsDrift() const { return is_drift_; }
    
    /**
     * @brief Print matrix information
     */
    void Print() const;

private:
    std::vector<MatrixElement> elements_;   ///< All matrix elements
    int max_order_{0};                      ///< Maximum polynomial order
    bool is_loaded_{false};                 ///< Matrix loaded flag
    bool is_drift_{false};                  ///< Is this a field-free drift?
    double drift_distance_{0.0};            ///< Drift distance if is_drift_==true (cm)
    
    /**
     * @brief Parse a line from COSY file
     * @param line Line from file
     * @param element Output matrix element
     * @return true if line parsed successfully
     * 
     * Format: c1 c2 c3 c4 c5 XXXXXX
     * where XXXXXX is a 6-digit string (each digit 0-9)
     * representing exponents for [x, xp, y, yp, TOF, delta]
     */
    bool ParseLine(const std::string& line, MatrixElement& element);
    
    /**
     * @brief Check if element is consistent with drift transformation
     * Based on transp.f lines 340-410
     */
    void CheckIfDrift();
};

/**
 * @class CosyMatrixSequential
 * @brief Sequential COSY matrix (transport through multiple elements)
 * 
 * IMPORTANT: The COSY matrices in SIMC are SEQUENTIAL, meaning you apply
 * multiple transformations in order (one per spectrometer element) rather
 * than a single end-to-end transformation.
 * 
 * From transp.f: Each "class" represents transport from one z-plane to the next
 */
class CosyMatrixSequential {
public:
    /**
     * @brief Add transport element (transformation class)
     * @param matrix COSY matrix for this element
     * @param name Element name (e.g., "Q1_ENTRANCE", "DIPOLE")
     * @param length Element length in meters (from comments in .dat file)
     */
    void AddElement(CosyMatrix matrix, 
                   const std::string& name,
                   double length);
    
    /**
     * @brief Apply sequential transport through all elements
     * @param input Input state [x, xp, y, yp, delta]
     * @param output Output state [x, xp, y, yp, delta]
     * @param total_path_length Total path length accumulated (cm)
     */
    void Apply(const double input[5], double output[5], double& total_path_length) const;
    
    /**
     * @brief Get number of elements (transformation classes)
     */
    size_t GetNumElements() const { return elements_.size(); }
    
    /**
     * @brief Get element name by index
     */
    std::string GetElementName(size_t index) const;
    
    /**
     * @brief Print sequence information
     */
    void Print() const;

private:
    struct Element {
        CosyMatrix matrix;
        std::string name;
        double length{0.0};  // meters (from file comments)
    };
    
    std::vector<Element> elements_;
};

} // namespace simc

#endif // SIMC_COSY_MATRIX_H
