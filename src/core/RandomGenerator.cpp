// src/core/RandomGenerator.cpp
// Implementation of random number generation

#include "simc/RandomGenerator.h"
#include <iostream>
#include <stdexcept>

namespace simc {

// Global instance
RandomGenerator g_random;

// ============================================================================
// Constructors
// ============================================================================
RandomGenerator::RandomGenerator() {
    Initialize(0);  // Use random_device
}

RandomGenerator::RandomGenerator(unsigned int seed) {
    Initialize(seed);
}

// ============================================================================
// Initialization
// ============================================================================
void RandomGenerator::Initialize(unsigned int seed) {
    if (seed == 0) {
        // Use random_device to get a truly random seed
        std::random_device rd;
        current_seed_ = rd();
    } else {
        current_seed_ = seed;
    }
    
    engine_.seed(current_seed_);
    
    // Initialize distributions
    uniform_dist_ = std::uniform_real_distribution<double>(0.0, 1.0);
    gaussian_dist_ = std::normal_distribution<double>(0.0, 1.0);
    exp_dist_ = std::exponential_distribution<double>(1.0);
    
    count_ = 0;
}

// ============================================================================
// Random Number Generation
// ============================================================================
double RandomGenerator::Uniform() {
    ++count_;
    return uniform_dist_(engine_);
}

double RandomGenerator::Uniform(double min, double max) {
    ++count_;
    return min + (max - min) * uniform_dist_(engine_);
}

double RandomGenerator::Gaussian(double mean, double sigma) {
    ++count_;
    return mean + sigma * gaussian_dist_(engine_);
}

double RandomGenerator::Exponential(double lambda) {
    ++count_;
    std::exponential_distribution<double> dist(lambda);
    return dist(engine_);
}

int RandomGenerator::UniformInt(int min, int max) {
    ++count_;
    std::uniform_int_distribution<int> dist(min, max);
    return dist(engine_);
}

bool RandomGenerator::Boolean(double probability) {
    ++count_;
    return uniform_dist_(engine_) < probability;
}

// ============================================================================
// Seed Management
// ============================================================================
void RandomGenerator::SetSeed(unsigned int seed) {
    Initialize(seed);
}

// ============================================================================
// State Persistence
// ============================================================================
bool RandomGenerator::SaveState(const std::string& filename) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file for writing: " << filename << std::endl;
        return false;
    }
    
    try {
        // Save seed and count
        file.write(reinterpret_cast<const char*>(&current_seed_), sizeof(current_seed_));
        file.write(reinterpret_cast<const char*>(&count_), sizeof(count_));
        
        // Save engine state
        file << engine_;
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error saving random state: " << e.what() << std::endl;
        return false;
    }
}

bool RandomGenerator::RestoreState(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file for reading: " << filename << std::endl;
        return false;
    }
    
    try {
        // Restore seed and count
        file.read(reinterpret_cast<char*>(&current_seed_), sizeof(current_seed_));
        file.read(reinterpret_cast<char*>(&count_), sizeof(count_));
        
        // Restore engine state
        file >> engine_;
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error restoring random state: " << e.what() << std::endl;
        return false;
    }
}

// ============================================================================
// Utility Methods
// ============================================================================
void RandomGenerator::Discard(unsigned long long n) {
    engine_.discard(n);
    count_ += n;
}

// ============================================================================
// GaussianGenerator Implementation
// ============================================================================
GaussianGenerator::GaussianGenerator(RandomGenerator& rng, 
                                   double mean, 
                                   double sigma,
                                   double nsigma_max)
    : rng_(rng), mean_(mean), sigma_(sigma), nsigma_max_(nsigma_max) {
}

double GaussianGenerator::Generate() {
    if (nsigma_max_ <= 0.0) {
        // No truncation
        return rng_.Gaussian(mean_, sigma_);
    }
    
    // Generate with truncation
    double value;
    double max_dev = nsigma_max_ * sigma_;
    
    do {
        value = rng_.Gaussian(mean_, sigma_);
    } while (std::abs(value - mean_) > max_dev);
    
    return value;
}

void GaussianGenerator::SetParameters(double mean, double sigma, double nsigma_max) {
    mean_ = mean;
    sigma_ = sigma;
    nsigma_max_ = nsigma_max;
}

} // namespace simc
