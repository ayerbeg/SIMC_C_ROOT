// include/simc/CosyMatrix.h
// COSY matrix element representation and evaluation
// Ported from transp.f (COSY matrix parsing and evaluation)
// ROBUST VERSION - Handles all Fortran formatting variations

#ifndef SIMC_COSY_MATRIX_H
#define SIMC_COSY_MATRIX_H

#include <vector>
#include <string>
#include <cmath>

namespace simc {

/**
 * @class CosyMatrix
 * @brief Represents and evaluates COSY-generated transport matrices
 * 
 * COSY (COSmic raY) is a differential algebra code used to calculate
 * particle optics. It outputs transfer maps as polynomials.
 * 
 * FILE FORMAT (from actual COSY files in data/matrices/):
 * 
 * Lines contain coefficients and exponents in various formats:
 * 
 * Format A (6-digit exponents):
 *   coeff1 coeff2 coeff3 coeff4 coeff5 XXXXXX
 *   where XXXXXX = [x_exp, xp_exp, y_exp, yp_exp, TOF_exp, delta_exp]
 *   Example: 1.000000 0.0 0.0 0.0 0.0 100000
 *            means x^1 * xp^0 * y^0 * yp^0 * TOF^0 * delta^0
 * 
 * Format B (5-digit exponents):
 *   coeff1 coeff2 coeff3 coeff4 coeff5 XXXXX
 *   where XXXXX = [x_exp, xp_exp, y_exp, yp_exp, delta_exp]
 *   Example: 0.777935 -3.32185 0.0 0.0 0.0 10000
 * 
 * The parser handles:
 * - Scientific notation (E, e, D, d)
 * - Concatenated numbers without spaces (e.g., "E-01-0.67")
 * - Variable whitespace
 * - 4 or 5 coefficients per line
 * - 5 or 6 digit exponent fields
 * 
 * Units (COSY-7):
 *   Positions: cm
 *   Angles: mrad (milliradians)
 *   Delta: percent (%)
 * 
 * Vector indices:
 *   0: x      - horizontal position (cm)
 *   1: xp     - dx/dz horizontal angle (mrad)
 *   2: y      - vertical position (cm)
 *   3: yp     - dy/dz vertical angle (mrad)
 *   4: delta  - momentum deviation (%)
 * 
 * Based on transp.f from Fortran SIMC
 */
class CosyMatrix {
public:
    /**
     * @struct MatrixElement
     * @brief Single term in the transfer map polynomial
     * 
     * Each COSY line contains coefficients for ALL 5 outputs
     * with a single set of input exponents
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
            
            if (result == 0.0) return 0.0;
            
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
     * Handles multiple Fortran format variations:
     * - Comment lines (starting with '!' or '#')
     * - Separator lines (' ---...')
     * - Data lines with 4-5 coefficients + 5-6 digit exponents
     * - Scientific notation (E, e, D, d)
     * - Numbers concatenated without spaces
     */
    bool LoadFromFile(const std::string& filename);
    
    /**
     * @brief Apply matrix to input vector
     * @param input Input vector [x(cm), xp(mrad), y(cm), yp(mrad), delta(%)]
     * @param output Output vector [x(cm), xp(mrad), y(cm), yp(mrad), dL(cm)]
     * 
     * For each output variable, sum all matrix elements:
     *   output[i] = sum( coeff[i] * product(input[j]^exp[j]) )
     */
    void Apply(const double input[5], double output[5]) const;
    
    /**
     * @brief Get maximum order in matrix
     */
    int GetMaxOrder() const { return max_order_; }
    
    /**
     * @brief Get number of matrix elements
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
     * @brief Get all matrix elements (for testing/inspection)
     */
    const std::vector<MatrixElement>& GetElements() const;
    
    /**
     * @brief Print matrix information
     */
    void Print() const;
    
    /**
     * @brief Parse a single line from COSY matrix file
     * @param line Line to parse
     * @param coeffs Output: parsed coefficients
     * @param exponents Output: parsed exponents
     * @return true if line parsed successfully
     * 
     * Made public for testing purposes
     */
    bool ParseMatrixLine(const std::string& line,
                        std::vector<double>& coeffs,
                        std::vector<int>& exponents);

private:
    std::vector<MatrixElement> elements_;   ///< All matrix elements
    int max_order_{0};                      ///< Maximum polynomial order
    bool is_loaded_{false};                 ///< Matrix loaded flag
    bool is_drift_{false};                  ///< Is this a field-free drift?
    double drift_distance_{0.0};            ///< Drift distance if drift (cm)
    
    /**
     * @brief Check if element is consistent with drift transformation
     */
    void CheckIfDrift();
};

/**
 * @class CosyMatrixSequential
 * @brief Sequential COSY matrix (transport through multiple elements)
 * 
 * The COSY matrices in SIMC are SEQUENTIAL - you apply multiple
 * transformations in order (one per spectrometer element) rather
 * than a single end-to-end transformation.
 */
class CosyMatrixSequential {
public:
    /**
     * @brief Add transport element
     * @param matrix COSY matrix for this element
     * @param name Element name (e.g., "Q1_ENTRANCE", "DIPOLE")
     * @param length Element length in meters
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
     * @brief Get number of elements
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
        double length{0.0};  // meters
    };
    
    std::vector<Element> elements_;
};

} // namespace simc

#endif // SIMC_COSY_MATRIX_H
