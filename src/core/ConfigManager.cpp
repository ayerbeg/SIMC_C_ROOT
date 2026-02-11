// src/core/ConfigManager.cpp
// Implementation of configuration management

#include "simc/core/ConfigManager.h"
#include "simc/core/SimcConstants.h"
#include <fstream>
#include <sstream>
#include <iostream>

using json = nlohmann::json;

namespace simc {

// ============================================================================
// Constructor
// ============================================================================
ConfigManager::ConfigManager(const std::string& filename) {
    LoadConfig(filename);
}

// ============================================================================
// Configuration Loading
// ============================================================================
bool ConfigManager::LoadConfig(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open configuration file: " << filename << std::endl;
        return false;
    }
    
    try {
        file >> config_;
        is_loaded_ = true;
        return true;
    } catch (const json::exception& e) {
        std::cerr << "Error parsing JSON: " << e.what() << std::endl;
        is_loaded_ = false;
        return false;
    }
}

bool ConfigManager::LoadFromString(const std::string& json_string) {
    try {
        config_ = json::parse(json_string);
        is_loaded_ = true;
        return true;
    } catch (const json::exception& e) {
        std::cerr << "Error parsing JSON string: " << e.what() << std::endl;
        is_loaded_ = false;
        return false;
    }
}

// ============================================================================
// Path Navigation
// ============================================================================
const json* ConfigManager::Navigate(const std::string& path) const {
    if (!is_loaded_) {
        return nullptr;
    }
    
    // Split path by '.'
    std::vector<std::string> keys;
    std::stringstream ss(path);
    std::string key;
    while (std::getline(ss, key, '.')) {
        keys.push_back(key);
    }
    
    // Navigate through JSON
    const json* current = &config_;
    for (const auto& k : keys) {
        if (!current->contains(k)) {
            return nullptr;
        }
        current = &(*current)[k];
    }
    
    return current;
}

bool ConfigManager::Has(const std::string& path) const {
    return Navigate(path) != nullptr;
}

// ============================================================================
// Structured Access Methods
// ============================================================================
BeamParameters ConfigManager::GetBeamConfig() const {
    BeamParameters beam;
    
    beam.energy = Get<double>("beam.energy", 10.6) * constants::GEV_TO_MEV;
    beam.energy_spread = Get<double>("beam.energy_spread", 0.0) * constants::GEV_TO_MEV;
    beam.x_width = Get<double>("beam.x_width", 0.2);
    beam.y_width = Get<double>("beam.y_width", 0.2);
    beam.x_offset = Get<double>("beam.x_offset", 0.0);
    beam.y_offset = Get<double>("beam.y_offset", 0.0);
    
    return beam;
}

SpectrometerSettings ConfigManager::GetSpecConfig(const std::string& arm_name) const {
    SpectrometerSettings spec;
    
    std::string base = "spectrometers." + arm_name;
    
    spec.P = Get<double>(base + ".momentum", 5.0) * constants::GEV_TO_MEV;
    spec.theta = Get<double>(base + ".angle", 15.0) * constants::DEG_TO_RAD;
    spec.phi = Get<double>(base + ".phi", 0.0) * constants::DEG_TO_RAD;
    
    // Offsets (optional)
    spec.offset.x = Get<double>(base + ".offset.x", 0.0);
    spec.offset.y = Get<double>(base + ".offset.y", 0.0);
    spec.offset.z = Get<double>(base + ".offset.z", 0.0);
    spec.offset.xptar = Get<double>(base + ".offset.xptar", 0.0);
    spec.offset.yptar = Get<double>(base + ".offset.yptar", 0.0);
    
    return spec;
}

ReactionType ConfigManager::GetReactionType() const {
    std::string type_str = Get<std::string>("reaction.type", "elastic");
    return StringToReactionType(type_str);
}

int ConfigManager::GetNEvents() const {
    return Get<int>("monte_carlo.nevents", 10000);
}

unsigned int ConfigManager::GetRandomSeed() const {
    return Get<unsigned int>("monte_carlo.random_seed", 0);
}

std::string ConfigManager::GetOutputFile() const {
    return Get<std::string>("output.filename", "simc_output.root");
}

// ============================================================================
// Target Configuration
// ============================================================================
std::string ConfigManager::GetTargetType() const {
    return Get<std::string>("target.type", "LH2");
}

double ConfigManager::GetTargetThickness() const {
    return Get<double>("target.thickness", 0.1);
}

double ConfigManager::GetTargetAngle() const {
    return Get<double>("target.angle", 0.0) * constants::DEG_TO_RAD;
}

int ConfigManager::GetTargetA() const {
    return Get<int>("target.A", 1);
}

int ConfigManager::GetTargetZ() const {
    return Get<int>("target.Z", 1);
}

// ============================================================================
// Monte Carlo Settings
// ============================================================================
bool ConfigManager::UseRadiativeCorrections() const {
    return Get<bool>("monte_carlo.radiative_corrections", true);
}

bool ConfigManager::UseEnergyLoss() const {
    return Get<bool>("monte_carlo.energy_loss", true);
}

bool ConfigManager::UseMultipleScattering() const {
    return Get<bool>("monte_carlo.multiple_scattering", true);
}

bool ConfigManager::UseCoulomb() const {
    return Get<bool>("monte_carlo.coulomb_corrections", true);
}

bool ConfigManager::DoDecay() const {
    return Get<bool>("monte_carlo.particle_decay", false);
}

// ============================================================================
// Validation
// ============================================================================
bool ConfigManager::Validate() const {
    if (!is_loaded_) {
        throw std::runtime_error("Configuration not loaded");
    }
    
    // Check required fields
    std::vector<std::string> required = {
        "beam.energy",
        "target.type",
        "reaction.type",
        "spectrometers.electron_arm.type",
        "spectrometers.electron_arm.momentum",
        "spectrometers.electron_arm.angle"
    };
    
    for (const auto& path : required) {
        if (!Has(path)) {
            throw std::runtime_error("Missing required configuration: " + path);
        }
    }
    
    // Validate ranges
    double beam_energy = Get<double>("beam.energy");
    if (beam_energy <= 0.0 || beam_energy > 100.0) {
        throw std::runtime_error("Invalid beam energy (must be 0-100 GeV)");
    }
    
    return true;
}

// ============================================================================
// Output Methods
// ============================================================================
void ConfigManager::Print() const {
    if (!is_loaded_) {
        std::cout << "Configuration not loaded" << std::endl;
        return;
    }
    
    std::cout << "========== SIMC Configuration ==========" << std::endl;
    std::cout << config_.dump(2) << std::endl;
    std::cout << "========================================" << std::endl;
}

void ConfigManager::SaveToFile(const std::string& filename) const {
    if (!is_loaded_) {
        throw std::runtime_error("Configuration not loaded");
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    file << config_.dump(2) << std::endl;
}

// ============================================================================
// String Conversion Methods
// ============================================================================
ReactionType ConfigManager::StringToReactionType(const std::string& str) const {
    if (str == "elastic") return ReactionType::ELASTIC;
    if (str == "quasielastic") return ReactionType::QUASIELASTIC;
    if (str == "pion_production") return ReactionType::PION_PRODUCTION;
    if (str == "kaon_production") return ReactionType::KAON_PRODUCTION;
    if (str == "delta_production") return ReactionType::DELTA_PRODUCTION;
    if (str == "rho_production") return ReactionType::RHO_PRODUCTION;
    if (str == "semi_inclusive") return ReactionType::SEMI_INCLUSIVE;
    
    throw std::runtime_error("Unknown reaction type: " + str);
}

SpectrometerType ConfigManager::StringToSpectrometerType(const std::string& str) const {
    if (str == "HMS") return SpectrometerType::HMS;
    if (str == "SHMS") return SpectrometerType::SHMS;
    if (str == "SOS") return SpectrometerType::SOS;
    if (str == "HRSL") return SpectrometerType::HRSL;
    if (str == "HRSR") return SpectrometerType::HRSR;
    
    throw std::runtime_error("Unknown spectrometer type: " + str);
}

} // namespace simc
