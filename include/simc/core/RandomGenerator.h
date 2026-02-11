// include/simc/RandomGenerator.h
// Modern C++ random number generation for SIMC Monte Carlo

#ifndef SIMC_RANDOM_GENERATOR_H
#define SIMC_RANDOM_GENERATOR_H

#include <random>
#include <string>
#include <fstream>

namespace simc {

/**
 * @class RandomGenerator
 * @brief High-quality random number generator using Mersenne Twister
 * 
 * This class provides random number generation functionality to replace
 * the Fortran RANLUX routines. It uses std::mt19937_64 (64-bit Mersenne
 * Twister) which has excellent statistical properties and a period of 2^19937-1.
 * 
 * Features:
 * - Multiple distribution types (uniform, gaussian, exponential)
 * - State persistence (save/restore for reproducibility)
 * - Thread-safe if each thread has its own instance
 * 
 * Example usage:
 * @code
 * RandomGenerator rng(12345);  // Fixed seed
 * double x = rng.Uniform();    // [0, 1)
 * double y = rng.Gaussian(0.0, 1.0);  // mean=0, sigma=1
 * rng.SaveState("random_state.dat");
 * @endcode
 */
class RandomGenerator {
public:
    // ========================================================================
    // Constructors
    // ========================================================================
    
    /**
     * @brief Default constructor with random seed from system
     */
    RandomGenerator();
    
    /**
     * @brief Constructor with specific seed
     * @param seed Random seed (0 = use random_device)
     */
    explicit RandomGenerator(unsigned int seed);
    
    /**
     * @brief Destructor
     */
    ~RandomGenerator() = default;
    
    // ========================================================================
    // Random Number Generation
    // ========================================================================
    
    /**
     * @brief Generate uniform random number in [0, 1)
     * @return Random number in range [0.0, 1.0)
     */
    double Uniform();
    
    /**
     * @brief Generate uniform random number in [min, max)
     * @param min Minimum value (inclusive)
     * @param max Maximum value (exclusive)
     * @return Random number in range [min, max)
     */
    double Uniform(double min, double max);
    
    /**
     * @brief Generate Gaussian (normal) distributed random number
     * @param mean Mean of distribution
     * @param sigma Standard deviation
     * @return Gaussian distributed random number
     */
    double Gaussian(double mean = 0.0, double sigma = 1.0);
    
    /**
     * @brief Generate exponentially distributed random number
     * @param lambda Rate parameter (1/mean)
     * @return Exponentially distributed random number
     */
    double Exponential(double lambda = 1.0);
    
    /**
     * @brief Generate random integer in range [min, max]
     * @param min Minimum value (inclusive)
     * @param max Maximum value (inclusive)
     * @return Random integer in range [min, max]
     */
    int UniformInt(int min, int max);
    
    /**
     * @brief Generate random boolean with given probability
     * @param probability Probability of returning true (default 0.5)
     * @return Random boolean
     */
    bool Boolean(double probability = 0.5);
    
    // ========================================================================
    // Seed Management
    // ========================================================================
    
    /**
     * @brief Set random seed
     * @param seed New seed (0 = use random_device)
     */
    void SetSeed(unsigned int seed);
    
    /**
     * @brief Get current seed
     * @return Current seed value
     */
    unsigned int GetSeed() const { return current_seed_; }
    
    // ========================================================================
    // State Persistence
    // ========================================================================
    
    /**
     * @brief Save current generator state to file
     * @param filename Output filename
     * @return true if successful
     * 
     * This allows exact reproduction of random sequences by saving
     * and restoring the generator state.
     */
    bool SaveState(const std::string& filename) const;
    
    /**
     * @brief Restore generator state from file
     * @param filename Input filename
     * @return true if successful
     */
    bool RestoreState(const std::string& filename);
    
    // ========================================================================
    // Statistics and Testing
    // ========================================================================
    
    /**
     * @brief Get number of random numbers generated
     * @return Count of generated numbers
     */
    unsigned long long GetCount() const { return count_; }
    
    /**
     * @brief Reset counter
     */
    void ResetCount() { count_ = 0; }
    
    /**
     * @brief Discard n random numbers (for jumping ahead)
     * @param n Number of values to discard
     */
    void Discard(unsigned long long n);
    
private:
    // ========================================================================
    // Private Members
    // ========================================================================
    std::mt19937_64 engine_;                    ///< 64-bit Mersenne Twister
    std::uniform_real_distribution<double> uniform_dist_; ///< [0,1) distribution
    std::normal_distribution<double> gaussian_dist_;      ///< Gaussian distribution
    std::exponential_distribution<double> exp_dist_;      ///< Exponential distribution
    
    unsigned int current_seed_{0};              ///< Current seed value
    unsigned long long count_{0};               ///< Number of generated values
    
    // ========================================================================
    // Private Methods
    // ========================================================================
    
    /**
     * @brief Initialize with seed (0 = random)
     */
    void Initialize(unsigned int seed);
};

/**
 * @class GaussianGenerator
 * @brief Specialized Gaussian random number generator with truncation
 * 
 * This class provides Gaussian random numbers with optional truncation
 * at a specified number of standard deviations. Useful for multiple
 * scattering and resolution smearing.
 */
class GaussianGenerator {
public:
    /**
     * @brief Constructor
     * @param rng Reference to random generator
     * @param mean Mean of distribution
     * @param sigma Standard deviation
     * @param nsigma_max Maximum number of sigmas (0 = no truncation)
     */
    GaussianGenerator(RandomGenerator& rng, 
                     double mean = 0.0, 
                     double sigma = 1.0,
                     double nsigma_max = 0.0);
    
    /**
     * @brief Generate random number
     * @return Gaussian random number (possibly truncated)
     */
    double Generate();
    
    /**
     * @brief Set parameters
     */
    void SetParameters(double mean, double sigma, double nsigma_max = 0.0);
    
private:
    RandomGenerator& rng_;
    double mean_;
    double sigma_;
    double nsigma_max_;
};

/**
 * @brief Global random generator instance (optional convenience)
 * 
 * While generally each component should have its own generator,
 * this provides a global instance for simple use cases.
 */
extern RandomGenerator g_random;

} // namespace simc

#endif // SIMC_RANDOM_GENERATOR_H
