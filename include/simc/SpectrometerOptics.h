// include/simc/SpectrometerOptics.h
// Spectrometer optics and transport using COSY matrices
// Ported from transp.f, mc_hms.f, mc_shms.f, mc_sos.f

#ifndef SIMC_SPECTROMETER_OPTICS_H
#define SIMC_SPECTROMETER_OPTICS_H

#include "SimcTypes.h"
#include "SimcEvent.h"
#include "CosyMatrix.h"
#include "EnergyLoss.h"
#include "MultipleScattering.h"
#include "RandomGenerator.h"
#include <string>
#include <vector>
#include <memory>

namespace simc {

/**
 * @class SpectrometerOptics
 * @brief Base class for spectrometer transport and optics
 * 
 * This class handles particle transport through magnetic spectrometers
 * using COSY-generated matrix elements. It includes:
 * - Forward transport (target to focal plane)
 * - Reconstruction (focal plane to target)
 * - Multiple scattering during transport
 * - Energy loss in materials
 * - Aperture checking at multiple planes
 * 
 * Coordinate System (TRANSPORT convention):
 * - x: horizontal (dispersive direction), positive = beam right
 * - y: vertical, positive = up
 * - z: along central ray
 * - xp = dx/dz (horizontal angle)
 * - yp = dy/dz (vertical angle)
 * - delta = (p - p0)/p0 * 100 (momentum deviation in %)
 * 
 * CRITICAL UNITS (must match COSY matrix format):
 * - Positions (x, y): cm
 * - Angles (xp, yp): MILLIRADIANS (mrad) [NOT radians!]
 * - Momentum (p): MeV/c
 * - Delta (δ): PERCENT (%) [NOT fraction!]
 * - Path lengths: cm
 * 
 * IMPORTANT: COSY matrices use milliradians for angles!
 * When interfacing with other code that uses radians, you must convert:
 *   angle_mrad = angle_rad * 1000.0
 *   angle_rad = angle_mrad / 1000.0
 * 
 * Based on Fortran SIMC transp.f and spectrometer-specific files
 */
class SpectrometerOptics {
public:
    /**
     * @struct TransportResult
     * @brief Result of transport calculation
     */
    struct TransportResult {
        bool success{false};          ///< Particle made it through
        ArmState final_state;         ///< Final coordinates
        FocalPlaneState focal_plane;  ///< Focal plane coordinates
        double path_length{0.0};      ///< Total path length (cm)
        int stop_plane{-1};           ///< Plane where particle stopped (-1 = success)
        std::string stop_reason;      ///< Why particle stopped
    };
    
    /**
     * @struct ApertureCheck
     * @brief Result of aperture checking
     */
    struct ApertureCheck {
        bool passed{true};            ///< Particle passed aperture
        int plane{-1};                ///< Plane number where stopped
        std::string name;             ///< Aperture name
        double x{0.0}, y{0.0};        ///< Position where stopped (cm)
    };
    
    /**
     * @brief Constructor
     * @param name Spectrometer name (HMS, SHMS, SOS, HRSL, HRSR)
     * @param forward_matrix_file Path to forward COSY matrix
     * @param recon_matrix_file Path to reconstruction COSY matrix
     */
    SpectrometerOptics(const std::string& name,
                      const std::string& forward_matrix_file,
                      const std::string& recon_matrix_file);
    
    virtual ~SpectrometerOptics() = default;
    
    /**
     * @brief Transport particle from target to focal plane (forward)
     * 
     * IMPORTANT: Input angles (state.xptar, state.yptar) should be in RADIANS
     * They will be converted to milliradians internally for COSY matrix application
     * 
     * @param state Initial state at target (spectrometer coordinates)
     *              Angles in RADIANS (will be converted to mrad internally)
     * @param momentum Particle momentum (MeV/c)
     * @param mass Particle mass (MeV/c²)
     * @param charge Particle charge (units of e)
     * @param rng Random number generator (for multiple scattering)
     * @param do_energy_loss Apply energy loss corrections
     * @param do_multiple_scattering Apply multiple scattering
     * @return TransportResult with final state and success flag
     */
    TransportResult TransportForward(const ArmState& state,
                                    double momentum,
                                    double mass,
                                    double charge,
                                    RandomGenerator& rng,
                                    bool do_energy_loss = true,
                                    bool do_multiple_scattering = true);
    
    /**
     * @brief Reconstruct target coordinates from focal plane
     * 
     * IMPORTANT: Output angles will be in RADIANS (converted from mrad internally)
     * 
     * @param focal_plane Measured focal plane coordinates
     *                    Angles in RADIANS (will be converted to mrad internally)
     * @param momentum Central momentum (MeV/c)
     * @return Reconstructed target coordinates with angles in RADIANS
     */
    ArmState Reconstruct(const FocalPlaneState& focal_plane,
                        double momentum);
    
    /**
     * @brief Check if particle passes apertures
     * @param x Horizontal position (cm)
     * @param y Vertical position (cm)
     * @param plane Aperture plane number
     * @return ApertureCheck result
     */
    virtual ApertureCheck CheckAperture(double x, double y, int plane) const = 0;
    
    /**
     * @brief Get number of aperture planes
     */
    virtual int GetNumAperturePlanes() const = 0;
    
    /**
     * @brief Get spectrometer name
     */
    std::string GetName() const { return name_; }
    
    /**
     * @brief Get central momentum (MeV/c)
     */
    double GetCentralMomentum() const { return p_central_; }
    
    /**
     * @brief Set central momentum (MeV/c)
     */
    void SetCentralMomentum(double p) { p_central_ = p; }
    
    /**
     * @brief Get spectrometer angle (radians)
     */
    double GetAngle() const { return theta_; }
    
    /**
     * @brief Set spectrometer angle (radians)
     */
    void SetAngle(double theta) { theta_ = theta; }

protected:
    std::string name_;                      ///< Spectrometer name
    double p_central_{5000.0};              ///< Central momentum (MeV/c)
    double theta_{0.0};                     ///< Spectrometer angle (radians)
    
    std::unique_ptr<CosyMatrix> forward_matrix_;   ///< Forward transport matrix
    std::unique_ptr<CosyMatrix> recon_matrix_;     ///< Reconstruction matrix
    
    /**
     * @brief Transport through a single drift section
     * @param state Current state (angles in RADIANS)
     * @param length Drift length (cm)
     * @param material Material in drift (nullptr = vacuum)
     * @param momentum Particle momentum (MeV/c)
     * @param mass Particle mass (MeV/c²)
     * @param charge Particle charge (units of e)
     * @param rng Random number generator
     * @param do_energy_loss Apply energy loss
     * @param do_multiple_scattering Apply multiple scattering
     */
    void TransportDrift(ArmState& state,
                       double length,
                       const Material* material,
                       double momentum,
                       double mass,
                       double charge,
                       RandomGenerator& rng,
                       bool do_energy_loss,
                       bool do_multiple_scattering);
    
    /**
     * @brief Apply COSY matrix to state vector
     * @param matrix COSY matrix to apply
     * @param state Input/output state
     * 
     * INTERNAL: Handles conversion between radians and milliradians
     */
    void ApplyCosyMatrix(const CosyMatrix& matrix, ArmState& state);
    
    /**
     * @brief Convert ArmState to COSY vector
     * 
     * CRITICAL: Converts angles from radians to milliradians!
     * Output: [x(cm), xp(mrad), y(cm), yp(mrad), delta(%)]
     * 
     * @param state Input state with angles in RADIANS
     * @return COSY vector with angles in MILLIRADIANS
     */
    std::vector<double> StateToCosyVector(const ArmState& state) const;
    
    /**
     * @brief Convert COSY vector to ArmState
     * 
     * CRITICAL: Converts angles from milliradians to radians!
     * Input: [x(cm), xp(mrad), y(cm), yp(mrad), delta(%)]
     * 
     * @param vec COSY vector with angles in MILLIRADIANS
     * @param momentum Central momentum (MeV/c)
     * @return ArmState with angles in RADIANS
     */
    ArmState CosyVectorToState(const std::vector<double>& vec, double momentum) const;
};

// ============================================================================
// HMS Spectrometer
// ============================================================================

/**
 * @class HMSOptics
 * @brief High Momentum Spectrometer (HMS) optics
 * 
 * HMS aperture geometry (from apertures_hms.inc and mc_hms.f):
 * - Octagon entrance aperture: 17.145 cm radius
 * - Dipole entrance/exit: ±30 cm horizontal, ±12.5 cm vertical
 * - Quadrupole apertures:
 *   * Q1: circular, r = 12.5 cm
 *   * Q2: circular, r = 30 cm
 *   * Q3: circular, r = 30 cm
 * - Detector hut: ±100 cm (large, rarely limiting)
 * 
 * Specifications:
 * - Central momentum range: 0.5 - 7.5 GeV/c
 * - Solid angle: ~6 msr (without collimator)
 * - Momentum acceptance: ±8-10%
 * - Angular acceptance: ±30-40 mrad (horizontal), ±65-80 mrad (vertical)
 */
class HMSOptics : public SpectrometerOptics {
public:
    explicit HMSOptics(const std::string& forward_file,
                      const std::string& recon_file);
    
    ApertureCheck CheckAperture(double x, double y, int plane) const override;
    int GetNumAperturePlanes() const override { return 8; }
    
private:
    /**
     * @brief HMS aperture checking functions
     * All dimensions in cm (from apertures_hms.inc)
     */
    bool CheckOctagon(double x, double y) const;  ///< Octagon: r=17.145 cm
    bool CheckDipole(double x, double y) const;   ///< Dipole: ±30 × ±12.5 cm
    bool CheckQuad1(double x, double y) const;    ///< Q1: r=12.5 cm
    bool CheckQuad2(double x, double y) const;    ///< Q2: r=30 cm
    bool CheckQuad3(double x, double y) const;    ///< Q3: r=30 cm
};

// ============================================================================
// SHMS Spectrometer
// ============================================================================

/**
 * @class SHMSOptics
 * @brief Super High Momentum Spectrometer (SHMS) optics
 * 
 * SHMS aperture geometry (from apertures_shms.inc and mc_shms.f):
 * - Entrance aperture: ±60 cm horizontal, ±25 cm vertical
 * - HB (Horizontal Bend) dipole: ±50 cm horizontal, ±25 cm vertical
 * - Quadrupole apertures:
 *   * Q1: circular, r = 14.92 cm
 *   * Q2: circular, r = 30 cm
 *   * Q3: circular, r = 30 cm
 * - Collimators:
 *   * LARGE: ±6 cm both directions
 *   * SMALL: ±3 cm both directions
 * - Detector hut: ±60 cm
 * 
 * Specifications:
 * - Central momentum range: 2.0 - 11.0 GeV/c
 * - Solid angle: ~4 msr (LARGE collimator)
 * - Momentum acceptance: ±20-25%
 * - Angular acceptance: ±50 mrad (horizontal), ±75-100 mrad (vertical)
 */
class SHMSOptics : public SpectrometerOptics {
public:
    explicit SHMSOptics(const std::string& forward_file,
                       const std::string& recon_file);
    
    ApertureCheck CheckAperture(double x, double y, int plane) const override;
    int GetNumAperturePlanes() const override { return 9; }
    
    /**
     * @brief Set collimator type
     * @param type "LARGE" (±6 cm) or "SMALL" (±3 cm)
     */
    void SetCollimator(const std::string& type);
    
private:
    std::string collimator_type_{"LARGE"};  ///< Collimator type
    
    bool CheckEntrance(double x, double y) const;    ///< Entrance: ±60 × ±25 cm
    bool CheckHBDipole(double x, double y) const;    ///< HB dipole: ±50 × ±25 cm
    bool CheckCollimator(double x, double y) const;  ///< LARGE ±6cm or SMALL ±3cm
};

// ============================================================================
// SOS Spectrometer
// ============================================================================

/**
 * @class SOSOptics
 * @brief Short Orbit Spectrometer (SOS) optics
 * 
 * SOS is a smaller spectrometer used in Hall C for lower momenta
 * 
 * Specifications:
 * - Central momentum range: 0.3 - 1.7 GeV/c
 * - Solid angle: ~7 msr
 * - Momentum acceptance: ±20%
 * 
 * Aperture geometry (approximate, from mc_sos.f):
 * - 6 aperture planes through dipoles and quadrupoles
 */
class SOSOptics : public SpectrometerOptics {
public:
    explicit SOSOptics(const std::string& forward_file,
                      const std::string& recon_file);
    
    ApertureCheck CheckAperture(double x, double y, int plane) const override;
    int GetNumAperturePlanes() const override { return 6; }
};

// ============================================================================
// HRS Spectrometers (Hall A)
// ============================================================================

/**
 * @class HRSOptics
 * @brief Hall A High Resolution Spectrometer
 * 
 * Both HRSL (Left) and HRSR (Right) use the same optics with
 * mirrored geometry
 * 
 * Specifications:
 * - Central momentum range: 0.3 - 4.3 GeV/c
 * - Momentum acceptance: ±4.5%
 * - Angular acceptance: ±28 mrad (horizontal), ±60 mrad (vertical)
 * - Solid angle: ~6 msr
 * 
 * Aperture geometry (from mc_hrsl.f / mc_hrsr.f):
 * - 7 aperture planes through quadrupoles, dipoles, and detector hut
 */
class HRSOptics : public SpectrometerOptics {
public:
    explicit HRSOptics(const std::string& name,  // "HRSL" or "HRSR"
                      const std::string& forward_file,
                      const std::string& recon_file);
    
    ApertureCheck CheckAperture(double x, double y, int plane) const override;
    int GetNumAperturePlanes() const override { return 7; }
};

// ============================================================================
// Factory Functions
// ============================================================================

/**
 * @brief Create spectrometer optics object
 * @param type Spectrometer type (HMS, SHMS, SOS, HRSL, HRSR)
 * @param forward_file Forward matrix file (empty = use default)
 * @param recon_file Reconstruction file (empty = use default)
 * @return Unique pointer to SpectrometerOptics
 */
std::unique_ptr<SpectrometerOptics> CreateSpectrometerOptics(
    const std::string& type,
    const std::string& forward_file = "",
    const std::string& recon_file = "");

/**
 * @brief Get default matrix file paths
 * @param type Spectrometer type
 * @return {forward_file, recon_file}
 */
std::pair<std::string, std::string> GetDefaultMatrixFiles(const std::string& type);

} // namespace simc

#endif // SIMC_SPECTROMETER_OPTICS_H
