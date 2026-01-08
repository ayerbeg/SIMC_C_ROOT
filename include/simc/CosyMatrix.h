// include/simc/CosyMatrix.h
// COSY matrix element representation and evaluation
// Ported from transp.f (COSY matrix parsing and evaluation)

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
 * Format:
 * Each output variable (x_out, xp_out, y_out, yp_out, delta_out) is
 * expressed as a polynomial in input variables (x, xp, y, yp, delta).
 * 
 * Example COSY line:
 *   I  COEFFICIENT              ORDER EXPONENTS
 *   1  1.234567E-02              3   2 1 0 0 0
 * 
 * This means: x_out += 1.234567e-02 * x^2 * xp^1
 * 
 * Vector indices (Fortran SIMC convention, 1-indexed in file, 0-indexed in code):
 *   0: x      - horizontal position (cm)
 *   1: xp     - dx/dz horizontal angle (rad)
 *   2: y      - vertical position (cm)
 *   3: yp     - dy/dz vertical angle (rad)
 *   4: delta  - momentum deviation (%)
 * 
 * Based on transp.f from Fortran SIMC
 */
class CosyMatrix {
public:
    /**
     * @struct MatrixElement
     * @brief Single term in the transfer map polynomial
     */
    struct MatrixElement {
        int output_index{0};        ///< Which output variable (0-4)
        double coefficient{0.0};    ///< Coefficient of this term
        int exponents[5]{0,0,0,0,0}; ///< Exponents for [x, xp, y, yp, delta]
        int order{0};               ///< Total order (sum of exponents)
        
        /**
         * @brief Evaluate this term
         * @param input Input vector [x, xp, y, yp, delta]
         * @return Contribution to output
         */
        double Evaluate(const double input[5]) const {
            double result = coefficient;
            
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
     * File format (from transp.f lines 50-150):
     * - Comment lines start with '!'
     * - Data lines: output_index coefficient order ex0 ex1 ex2 ex3 ex4
     * - Output indices: 1=x, 2=xp, 3=y, 4=yp, 5=delta (convert to 0-indexed)
     */
    bool LoadFromFile(const std::string& filename);
    
    /**
     * @brief Apply matrix to input vector
     * @param input Input vector [x, xp, y, yp, delta]
     * @param output Output vector [x, xp, y, yp, delta]
     * 
     * From transp.f lines 200-250:
     * For each output variable, sum all matrix elements:
     *   output[i] = sum_over_terms( coeff * product(input[j]^exponent[j]) )
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
     * @brief Get matrix elements for specific output
     * @param output_index Output variable (0-4)
     * @return Vector of matrix elements for that output
     */
    std::vector<MatrixElement> GetElementsForOutput(int output_index) const;
    
    /**
     * @brief Print matrix information
     */
    void Print() const;

private:
    std::vector<MatrixElement> elements_;   ///< All matrix elements
    int max_order_{0};                      ///< Maximum polynomial order
    bool is_loaded_{false};                 ///< Matrix loaded flag
    
    /**
     * @brief Parse a line from COSY file
     * @param line Line from file
     * @param element Output matrix element
     * @return true if line parsed successfully
     */
    bool ParseLine(const std::string& line, MatrixElement& element);
};

/**
 * @class CosyMatrixSequential
 * @brief Sequential COSY matrix (transport through multiple elements)
 * 
 * Some spectrometers use sequential matrices where transport is
 * calculated element-by-element rather than as a single transfer map.
 * This is less common but supported in original SIMC.
 * 
 * From transp.f: sequential transport through dipoles, quads, drifts
 */
class CosyMatrixSequential {
public:
    /**
     * @brief Add transport element
     * @param matrix COSY matrix for this element
     * @param type Element type ("drift", "dipole", "quad")
     * @param length Element length (cm)
     */
    void AddElement(const CosyMatrix& matrix, 
                   const std::string& type,
                   double length);
    
    /**
     * @brief Apply sequential transport
     * @param input Input state
     * @param output Output state
     */
    void Apply(const double input[5], double output[5]) const;
    
    /**
     * @brief Get number of elements
     */
    size_t GetNumElements() const { return elements_.size(); }

private:
    struct Element {
        CosyMatrix matrix;
        std::string type;
        double length;
    };
    
    std::vector<Element> elements_;
};

} // namespace simc

#endif // SIMC_COSY_MATRIX_H
