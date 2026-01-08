// include/simc/ConfigManager.h
// Configuration management for SIMC Monte Carlo

#ifndef SIMC_CONFIG_MANAGER_H
#define SIMC_CONFIG_MANAGER_H

#include "SimcTypes.h"
#include <string>
#include <map>
#include <stdexcept>
#include <nlohmann/json.hpp>

namespace simc {

/**
 * @class ConfigManager
 * @brief Manages simulation configuration from JSON files
 * 
 * This class reads and validates configuration parameters from JSON files.
 * It provides type-safe access to hierarchical configuration data.
 * 
 * Example JSON structure:
 * @code{.json}
 * {
 *   "beam": {
 *     "energy": 10.6,
 *     "energy_spread": 0.05,
 *     "x_width": 0.2,
 *     "y_width": 0.2
 *   },
 *   "target": {
 *     "type": "LH2",
 *     "thickness": 0.1,
 *     "angle": 0.0
 *   },
 *   "reaction": {
 *     "type": "pion_production",
 *     "pion_type": "pi_plus"
 *   },
 *   "spectrometers": {
 *     "electron_arm": {
 *       "type": "HMS",
 *       "momentum": 8.0,
 *       "angle": 12.5
 *     }
 *   }
 * }
 * @endcode
 */
class ConfigManager {
public:
    // ========================================================================
    // Constructors
    // ========================================================================
    
    /**
     * @brief Default constructor
     */
    ConfigManager() = default;
    
    /**
     * @brief Constructor with config file
     * @param filename Path to JSON configuration file
     */
    explicit ConfigManager(const std::string& filename);
    
    /**
     * @brief Destructor
     */
    ~ConfigManager() = default;
    
    // ========================================================================
    // Configuration Loading
    // ========================================================================
    
    /**
     * @brief Load configuration from JSON file
     * @param filename Path to JSON configuration file
     * @return true if successful, false otherwise
     */
    bool LoadConfig(const std::string& filename);
    
    /**
     * @brief Load configuration from JSON string
     * @param json_string JSON configuration as string
     * @return true if successful, false otherwise
     */
    bool LoadFromString(const std::string& json_string);
    
    /**
     * @brief Check if configuration is loaded
     */
    bool IsLoaded() const { return is_loaded_; }
    
    // ========================================================================
    // Generic Access Methods
    // ========================================================================
    
    /**
     * @brief Get value from configuration using path notation
     * @tparam T Type to return
     * @param path Dot-separated path (e.g., "beam.energy")
     * @return Value at path
     * @throws std::runtime_error if path not found or type mismatch
     */
    template<typename T>
    T Get(const std::string& path) const;
    
    /**
     * @brief Get value with default fallback
     * @tparam T Type to return
     * @param path Dot-separated path
     * @param default_value Default value if path not found
     * @return Value at path or default
     */
    template<typename T>
    T Get(const std::string& path, const T& default_value) const;
    
    /**
     * @brief Check if path exists in configuration
     * @param path Dot-separated path
     * @return true if path exists
     */
    bool Has(const std::string& path) const;
    
    // ========================================================================
    // Structured Access Methods
    // ========================================================================
    
    /**
     * @brief Get beam configuration
     */
    BeamParameters GetBeamConfig() const;
    
    /**
     * @brief Get spectrometer configuration
     * @param arm_name Name of arm ("electron_arm" or "hadron_arm")
     */
    SpectrometerSettings GetSpecConfig(const std::string& arm_name) const;
    
    /**
     * @brief Get reaction type
     */
    ReactionType GetReactionType() const;
    
    /**
     * @brief Get number of events to generate
     */
    int GetNEvents() const;
    
    /**
     * @brief Get random seed
     */
    unsigned int GetRandomSeed() const;
    
    /**
     * @brief Get output filename
     */
    std::string GetOutputFile() const;
    
    // ========================================================================
    // Target Configuration
    // ========================================================================
    
    /**
     * @brief Get target type (e.g., "LH2", "LD2", "C12")
     */
    std::string GetTargetType() const;
    
    /**
     * @brief Get target thickness (g/cm^2)
     */
    double GetTargetThickness() const;
    
    /**
     * @brief Get target angle (degrees)
     */
    double GetTargetAngle() const;
    
    /**
     * @brief Get target mass number A
     */
    int GetTargetA() const;
    
    /**
     * @brief Get target atomic number Z
     */
    int GetTargetZ() const;
    
    // ========================================================================
    // Monte Carlo Settings
    // ========================================================================
    
    /**
     * @brief Check if using radiative corrections
     */
    bool UseRadiativeCorrections() const;
    
    /**
     * @brief Check if using energy loss corrections
     */
    bool UseEnergyLoss() const;
    
    /**
     * @brief Check if using multiple scattering
     */
    bool UseMultipleScattering() const;
    
    /**
     * @brief Check if using Coulomb corrections
     */
    bool UseCoulomb() const;
    
    /**
     * @brief Check if doing particle decay
     */
    bool DoDecay() const;
    
    // ========================================================================
    // Validation
    // ========================================================================
    
    /**
     * @brief Validate configuration
     * @return true if configuration is valid
     * @throws std::runtime_error with description if invalid
     */
    bool Validate() const;
    
    /**
     * @brief Print configuration to stdout
     */
    void Print() const;
    
    /**
     * @brief Save configuration to file
     * @param filename Output filename
     */
    void SaveToFile(const std::string& filename) const;
    
private:
    // ========================================================================
    // Private Members
    // ========================================================================
    nlohmann::json config_;     ///< JSON configuration object
    bool is_loaded_{false};     ///< Configuration loaded flag
    
    // ========================================================================
    // Private Helper Methods
    // ========================================================================
    
    /**
     * @brief Navigate to JSON node using path
     * @param path Dot-separated path
     * @return Pointer to JSON node or nullptr
     */
    const nlohmann::json* Navigate(const std::string& path) const;
    
    /**
     * @brief Convert string to ReactionType
     */
    ReactionType StringToReactionType(const std::string& str) const;
    
    /**
     * @brief Convert string to SpectrometerType
     */
    SpectrometerType StringToSpectrometerType(const std::string& str) const;
};

// ============================================================================
// Template Implementations
// ============================================================================

template<typename T>
T ConfigManager::Get(const std::string& path) const {
    if (!is_loaded_) {
        throw std::runtime_error("Configuration not loaded");
    }
    
    const nlohmann::json* node = Navigate(path);
    if (!node) {
        throw std::runtime_error("Configuration path not found: " + path);
    }
    
    try {
        return node->get<T>();
    } catch (const nlohmann::json::exception& e) {
        throw std::runtime_error("Type mismatch for path: " + path + 
                                 " - " + e.what());
    }
}

template<typename T>
T ConfigManager::Get(const std::string& path, const T& default_value) const {
    if (!is_loaded_) {
        return default_value;
    }
    
    const nlohmann::json* node = Navigate(path);
    if (!node) {
        return default_value;
    }
    
    try {
        return node->get<T>();
    } catch (const nlohmann::json::exception&) {
        return default_value;
    }
}

} // namespace simc

#endif // SIMC_CONFIG_MANAGER_H
