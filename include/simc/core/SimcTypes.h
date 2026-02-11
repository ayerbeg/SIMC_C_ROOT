// include/simc/SimcTypes.h
// Common type definitions for SIMC Monte Carlo

#ifndef SIMC_TYPES_H
#define SIMC_TYPES_H

#include <cmath>
#include <string>

namespace simc {

// ============================================================================
// 3D Vector Class
// ============================================================================
/**
 * @brief 3D vector for particle momenta, positions, and directions
 */
struct Vector3D {
    double x{0.0};  ///< X component
    double y{0.0};  ///< Y component
    double z{0.0};  ///< Z component
    
    /// Calculate magnitude of vector
    double Magnitude() const { 
        return std::sqrt(x*x + y*y + z*z); 
    }
    
    /// Calculate magnitude squared (more efficient when magnitude not needed)
    double MagnitudeSquared() const {
        return x*x + y*y + z*z;
    }
    
    /// Return normalized (unit) vector
    Vector3D Normalized() const {
        double mag = Magnitude();
        if (mag > 0.0) {
            return {x/mag, y/mag, z/mag};
        }
        return {0.0, 0.0, 0.0};
    }
    
    /// Dot product with another vector
    double Dot(const Vector3D& v) const {
        return x*v.x + y*v.y + z*v.z;
    }
    
    /// Cross product with another vector
    Vector3D Cross(const Vector3D& v) const {
        return {
            y*v.z - z*v.y,
            z*v.x - x*v.z,
            x*v.y - y*v.x
        };
    }
    
    /// Vector addition
    Vector3D operator+(const Vector3D& v) const {
        return {x + v.x, y + v.y, z + v.z};
    }
    
    /// Vector subtraction
    Vector3D operator-(const Vector3D& v) const {
        return {x - v.x, y - v.y, z - v.z};
    }
    
    /// Scalar multiplication
    Vector3D operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }
};

// ============================================================================
// Particle State in Spectrometer Coordinates
// ============================================================================
/**
 * @brief Particle state in both spectrometer and physics coordinates
 * 
 * Spectrometer coordinates (TRANSPORT system):
 * - delta: momentum deviation from central momentum (%)
 * - xptar: horizontal angle at target (radians)
 * - yptar: vertical angle at target (radians)
 * - z: position along beamline (cm)
 * 
 * Physics coordinates:
 * - theta: scattering angle from beamline (radians)
 * - phi: azimuthal angle (radians)
 * - E: particle energy (MeV)
 * - P: particle momentum (MeV/c)
 */
struct ArmState {
    // Spectrometer coordinates
    double delta{0.0};   ///< Momentum deviation (%)
    double xptar{0.0};   ///< Horizontal angle (rad)
    double yptar{0.0};   ///< Vertical angle (rad)
    double z{0.0};       ///< Position along beamline (cm)
    
    // Physics coordinates
    double theta{0.0};   ///< Scattering angle (rad)
    double phi{0.0};     ///< Azimuthal angle (rad)
    double E{0.0};       ///< Energy (MeV)
    double P{0.0};       ///< Momentum (MeV/c)
};

// ============================================================================
// Focal Plane State
// ============================================================================
/**
 * @brief Particle coordinates at spectrometer focal plane
 */
struct FocalPlaneState {
    double x{0.0};       ///< Horizontal position (cm)
    double y{0.0};       ///< Vertical position (cm)
    double dx{0.0};      ///< Horizontal angle (rad)
    double dy{0.0};      ///< Vertical angle (rad)
    double path{0.0};    ///< Path length from target (cm)
};

// ============================================================================
// Target Information
// ============================================================================
/**
 * @brief Target vertex and interaction properties
 */
struct TargetInfo {
    // Vertex position (cm)
    double x{0.0};
    double y{0.0};
    double z{0.0};
    
    // Raster position (cm)
    double rasterx{0.0};
    double rastery{0.0};
    
    // Energy loss corrections (MeV)
    // [0] = beam, [1] = electron, [2] = hadron
    double Eloss[3]{0.0, 0.0, 0.0};
    
    // Effective target thickness (radiation lengths)
    // [0] = beam, [1] = electron, [2] = hadron
    double teff[3]{0.0, 0.0, 0.0};
    
    // Coulomb correction (MeV)
    double Coulomb{0.0};
};

// ============================================================================
// Spectrometer Offsets
// ============================================================================
/**
 * @brief Spectrometer alignment offsets
 */
struct SpecOffset {
    double x{0.0};       ///< Horizontal offset (cm)
    double y{0.0};       ///< Vertical offset (cm)
    double z{0.0};       ///< Longitudinal offset (cm)
    double xptar{0.0};   ///< Horizontal angle offset (rad)
    double yptar{0.0};   ///< Vertical angle offset (rad)
};

// ============================================================================
// Beam Parameters
// ============================================================================
/**
 * @brief Beam properties
 */
struct BeamParameters {
    double energy{0.0};         ///< Central energy (MeV)
    double energy_spread{0.0};  ///< Energy spread (MeV, full width)
    double x_width{0.0};        ///< Horizontal width, 1 sigma (cm)
    double y_width{0.0};        ///< Vertical width, 1 sigma (cm)
    double x_offset{0.0};       ///< Horizontal offset (cm)
    double y_offset{0.0};       ///< Vertical offset (cm)
};

// ============================================================================
// Spectrometer Settings
// ============================================================================
/**
 * @brief Spectrometer configuration
 */
struct SpectrometerSettings {
    double P{0.0};              ///< Central momentum (MeV/c)
    double theta{0.0};          ///< Central angle (rad)
    double phi{0.0};            ///< Azimuthal angle (rad)
    SpecOffset offset;          ///< Alignment offsets
};

// ============================================================================
// Enumerations
// ============================================================================

/**
 * @brief Types of physics reactions
 */
enum class ReactionType {
    ELASTIC,              ///< Elastic scattering
    QUASIELASTIC,         ///< Quasi-elastic (e,e'p)
    PION_PRODUCTION,      ///< Pion electroproduction
    KAON_PRODUCTION,      ///< Kaon electroproduction
    DELTA_PRODUCTION,     ///< Delta resonance production
    RHO_PRODUCTION,       ///< Rho meson production
    SEMI_INCLUSIVE        ///< Semi-inclusive DIS
};

/**
 * @brief Spectrometer types
 */
enum class SpectrometerType {
    HMS,    ///< High Momentum Spectrometer
    SHMS,   ///< Super High Momentum Spectrometer
    SOS,    ///< Short Orbit Spectrometer
    HRSL,   ///< Hall A Left HRS
    HRSR    ///< Hall A Right HRS
};

/**
 * @brief Pion types
 */
enum class PionType {
    PI_PLUS,    ///< pi+
    PI_MINUS,   ///< pi-
    PI_ZERO     ///< pi0
};

/**
 * @brief Kaon types
 */
enum class KaonType {
    K_PLUS,     ///< K+
    K_MINUS     ///< K-
};

/**
 * @brief Final state particle types
 */
enum class FinalStateType {
    NUCLEON,    ///< Free nucleon
    DELTA,      ///< Delta resonance
    BOUND       ///< Bound state
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Convert reaction type to string
 */
inline std::string ReactionTypeToString(ReactionType type) {
    switch(type) {
        case ReactionType::ELASTIC: return "Elastic";
        case ReactionType::QUASIELASTIC: return "QuasiElastic";
        case ReactionType::PION_PRODUCTION: return "PionProduction";
        case ReactionType::KAON_PRODUCTION: return "KaonProduction";
        case ReactionType::DELTA_PRODUCTION: return "DeltaProduction";
        case ReactionType::RHO_PRODUCTION: return "RhoProduction";
        case ReactionType::SEMI_INCLUSIVE: return "SemiInclusive";
        default: return "Unknown";
    }
}

/**
 * @brief Convert spectrometer type to string
 */
inline std::string SpectrometerTypeToString(SpectrometerType type) {
    switch(type) {
        case SpectrometerType::HMS: return "HMS";
        case SpectrometerType::SHMS: return "SHMS";
        case SpectrometerType::SOS: return "SOS";
        case SpectrometerType::HRSL: return "HRSL";
        case SpectrometerType::HRSR: return "HRSR";
        default: return "Unknown";
    }
}

} // namespace simc

#endif // SIMC_TYPES_H
