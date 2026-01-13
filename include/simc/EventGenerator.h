// include/simc/EventGenerator.h
// Event generation for SIMC Monte Carlo
// Ported from event.f, jacobians.f

#ifndef SIMC_EVENT_GENERATOR_H
#define SIMC_EVENT_GENERATOR_H

#include "SimcTypes.h"
#include "SimcEvent.h"
#include "SimcConstants.h"
#include "CrossSection.h"
#include "SpectrometerOptics.h"
#include "RandomGenerator.h"
#include "Kinematics.h"
#include "ConfigManager.h"
#include <memory>
#include <string>
#include <vector>

namespace simc {

// Use constants from SimcConstants.h
using namespace simc::constants;

/**
 * @struct GenerationLimits
 * @brief Phase space generation limits
 * 
 * Defines the kinematic ranges for Monte Carlo generation.
 * All angles in radians, energies in MeV, delta in percent.
 */
struct GenerationLimits {
    // Electron arm limits
    struct {
        double delta_min, delta_max;      ///< Momentum deviation (%)
        double yptar_min, yptar_max;      ///< Vertical angle (rad)
        double xptar_min, xptar_max;      ///< Horizontal angle (rad)
        double E_min, E_max;              ///< Energy (MeV)
    } electron;
    
    // Hadron arm limits
    struct {
        double delta_min, delta_max;
        double yptar_min, yptar_max;
        double xptar_min, xptar_max;
        double E_min, E_max;
    } hadron;
    
    // Additional generation constraints
    double sumEgen_min{0.0}, sumEgen_max{0.0};  ///< Sum energy constraint (MeV)
    double Trec_min{0.0}, Trec_max{0.0};        ///< Recoil kinetic energy (MeV)
    
    // Beam parameters
    double xwid{0.01}, ywid{0.01};  ///< Beam widths, 1 sigma (cm)
};

/**
 * @struct TargetProperties
 * @brief Target material and geometry properties
 */
struct TargetProperties {
    double mass{Mp};              ///< Target mass (MeV)
    int Z{1};                     ///< Atomic number
    int A{1};                     ///< Mass number
    double length{10.0};          ///< Target length (cm)
    double density{0.071};        ///< Density (g/cm³)
    double x_offset{0.0};         ///< X offset (cm)
    double y_offset{0.0};         ///< Y offset (cm)
    double z_offset{0.0};         ///< Z offset (cm)
    double angle{0.0};            ///< Target angle (rad)
    int raster_pattern{0};        ///< Raster pattern: 0=none, 1=bedpost, 2=circular, 3=flat
    double raster_x{0.0};         ///< Raster X size (cm)
    double raster_y{0.0};         ///< Raster Y size (cm)
    double coulomb_constant{0.0}; ///< Coulomb correction constant
    double coulomb_ave{0.0};      ///< Average Coulomb correction
};

/**
 * @struct FermiMomentum
 * @brief Fermi momentum distribution for nuclear targets
 */
struct FermiMomentum {
    double magnitude{0.0};        ///< |p_F| (MeV/c)
    double x{0.0}, y{0.0}, z{0.0};///< Components (MeV/c)
    double Em{0.0};               ///< Missing energy (MeV)
    
    void Clear() {
        magnitude = x = y = z = Em = 0.0;
    }
};

/**
 * @class EventGenerator
 * @brief Main event generator for SIMC Monte Carlo
 * 
 * This class handles Monte Carlo event generation including:
 * - Vertex position generation with beam profile and raster
 * - Phase space sampling based on reaction type
 * - Kinematic calculations
 * - Event weights and jacobians
 * - Integration with cross sections and spectrometer acceptance
 * 
 * Reference: event.f, jacobians.f from Fortran SIMC
 */
class EventGenerator {
public:
    // ========================================================================
    // Constructor and Configuration
    // ========================================================================
    
    /**
     * @brief Constructor
     * @param config Configuration manager
     * @param random Random number generator
     * @param cross_section Cross section calculator
     * @param optics Spectrometer optics
     */
    EventGenerator(
        const ConfigManager& config,
        std::shared_ptr<RandomGenerator> random,
        std::shared_ptr<CrossSection> cross_section,
        std::shared_ptr<SpectrometerOptics> optics
    );
    
    /**
     * @brief Destructor
     */
    ~EventGenerator() = default;
    
    /**
     * @brief Initialize generator with configuration
     * @return true if successful
     */
    bool Initialize();
    
    // ========================================================================
    // Main Event Generation
    // ========================================================================
    
    /**
     * @brief Generate a single event
     * @param vertex Output: generated vertex event
     * @param main Output: main event structure with weights
     * @return true if event passed all cuts
     * 
     * From event.f generate() subroutine
     */
    bool GenerateEvent(SimcEvent& vertex, MainEvent& main);
    
    /**
     * @brief Complete event kinematics from generated quantities
     * @param vertex Input/Output: event to complete
     * @param main Input: main event info
     * @return true if successful
     * 
     * From event.f complete_ev() subroutine
     */
    bool CompleteEvent(SimcEvent& vertex, MainEvent& main);
    
    /**
     * @brief Generate vertex position
     * @param target Output: target vertex info
     * 
     * Includes beam profile (gaussian), raster pattern, and position along target
     * From event.f generate() lines 150-230
     */
    void GenerateVertex(TargetInfo& target);
    
    /**
     * @brief Generate phase space variables
     * @param vertex Output: generated event
     * @param main Output: generation weights
     * @return true if generation successful
     * 
     * From event.f generate() lines 240-400
     */
    bool GeneratePhaseSpace(SimcEvent& vertex, MainEvent& main);
    
    // ========================================================================
    // Angle Conversions
    // ========================================================================
    
    /**
     * @brief Convert spectrometer angles to physics angles
     * @param theta0 Central spectrometer angle (rad)
     * @param phi0 Central spectrometer phi (rad)
     * @param xptar Horizontal target angle (rad)
     * @param yptar Vertical target angle (rad)
     * @param theta Output: physics scattering angle (rad)
     * @param phi Output: physics azimuthal angle (rad)
     * 
     * From event.f physics_angles() subroutine
     */
    void PhysicsAngles(
        double theta0, double phi0,
        double xptar, double yptar,
        double& theta, double& phi
    ) const;
    
    /**
     * @brief Convert physics angles to spectrometer angles
     * @param theta0 Central spectrometer angle (rad)
     * @param phi0 Central spectrometer phi (rad)
     * @param theta Physics scattering angle (rad)
     * @param phi Physics azimuthal angle (rad)
     * @param xptar Output: horizontal target angle (rad)
     * @param yptar Output: vertical target angle (rad)
     * 
     * From event.f spectrometer_angles() subroutine
     */
    void SpectrometerAngles(
        double theta0, double phi0,
        double theta, double phi,
        double& xptar, double& yptar
    ) const;
    
    // ========================================================================
    // Kinematics Calculations
    // ========================================================================
    
    /**
     * @brief Calculate Q², nu, W, xbj from electron scattering
     * @param event Event to update
     */
    void CalculateBasicKinematics(SimcEvent& event) const;
    
    /**
     * @brief Calculate missing momentum components
     * @param event Event to update
     */
    void CalculateMissingMomentum(SimcEvent& event) const;
    
    /**
     * @brief Calculate angles between q and hadron
     * @param event Event to update
     */
    void CalculateHadronAngles(SimcEvent& event) const;
    
    // ========================================================================
    // Reaction-Specific Calculations
    // ========================================================================
    
    /**
     * @brief Solve for hydrogen elastic kinematics
     * @param event Event to complete
     * @return true if successful
     */
    bool SolveHydrogenElastic(SimcEvent& event) const;
    
    /**
     * @brief Solve for deuterium kinematics
     * @param event Event to complete
     * @param main Main event for jacobian
     * @return true if successful
     */
    bool SolveDeuterium(SimcEvent& event, MainEvent& main) const;
    
    /**
     * @brief Solve for pion/kaon production kinematics
     * @param event Event to complete
     * @param pfer Fermi momentum
     * @return true if successful
     */
    bool SolvePionKaon(SimcEvent& event, const FermiMomentum& pfer) const;
    
    // ========================================================================
    // Fermi Momentum and Spectral Functions
    // ========================================================================
    
    /**
     * @brief Generate Fermi momentum from distribution
     * @param pfer Output: Fermi momentum
     */
    void GenerateFermiMomentum(FermiMomentum& pfer);
    
    /**
     * @brief Generate missing energy for nuclear targets
     * @param pfer Fermi momentum magnitude
     * @param Em Output: missing energy
     */
    void GenerateMissingEnergy(double pfer, double& Em);
    
    // ========================================================================
    // Jacobians
    // ========================================================================
    
    /**
     * @brief Calculate phase space jacobian
     * @param vertex Event kinematics
     * @param main Main event info
     * @return Jacobian value
     */
    double CalculateJacobian(const SimcEvent& vertex, const MainEvent& main) const;
    
    /**
     * @brief Calculate jacobian for center-of-mass transformation
     * @param vertex Event kinematics
     * @param main Main event info
     * @return CM jacobian
     */
    double CalculateCMJacobian(const SimcEvent& vertex, const MainEvent& main) const;
    
    // ========================================================================
    // Particle Decays
    // ========================================================================
    
    /**
     * @brief Perform two-body decay in rest frame
     * @param parent Parent 4-vector (E, px, py, pz) in MeV
     * @param decay1 Output: first decay product 4-vector
     * @param decay2 Output: second decay product 4-vector
     * @param m1 First product mass (MeV)
     * @param m2 Second product mass (MeV)
     * @param costh Decay angle cosine in CM frame
     * @param phi Decay azimuthal angle (rad)
     */
    void DecayTwoBody(
        const std::array<double, 4>& parent,
        std::array<double, 4>& decay1,
        std::array<double, 4>& decay2,
        double m1, double m2,
        double costh, double phi
    ) const;
    
    // ========================================================================
    // Configuration and Access
    // ========================================================================
    
    /**
     * @brief Set reaction type
     */
    void SetReactionType(ReactionType type) { reaction_type_ = type; }
    
    /**
     * @brief Get reaction type
     */
    ReactionType GetReactionType() const { return reaction_type_; }
    
    /**
     * @brief Set generation limits
     */
    void SetGenerationLimits(const GenerationLimits& limits) { gen_limits_ = limits; }
    
    /**
     * @brief Get generation limits
     */
    const GenerationLimits& GetGenerationLimits() const { return gen_limits_; }
    
    /**
     * @brief Enable/disable energy loss corrections
     */
    void SetUseEnergyLoss(bool use) { use_energy_loss_ = use; }
    
    /**
     * @brief Enable/disable Coulomb corrections
     */
    void SetUseCoulomb(bool use) { use_coulomb_ = use; }
    
    /**
     * @brief Enable/disable radiative corrections
     */
    void SetUseRadiative(bool use) { use_radiative_ = use; }
    
    /**
     * @brief Get number of events generated
     */
    unsigned int GetNGenerated() const { return n_generated_; }
    
    /**
     * @brief Get number of events accepted
     */
    unsigned int GetNAccepted() const { return n_accepted_; }
    
    /**
     * @brief Reset statistics
     */
    void ResetStatistics();

private:
    // ========================================================================
    // Helper Methods
    // ========================================================================
    
    /**
     * @brief Calculate energy loss through target
     * @param particle_id Particle type (1=beam, 2=electron, 3=hadron)
     * @param z_position Position along target (cm)
     * @param energy Particle energy (MeV)
     * @param theta Scattering angle (rad)
     * @param eloss Output: energy loss (MeV)
     * @param teff Output: effective thickness (g/cm²)
     */
    void TripThruTarget(
        int particle_id,
        double z_position,
        double energy,
        double theta,
        double& eloss,
        double& teff
    ) const;
    
    /**
     * @brief Calculate unit vectors from angles
     * @param theta Polar angle (rad)
     * @param phi Azimuthal angle (rad)
     * @return Unit vector
     */
    Vector3D UnitVector(double theta, double phi) const;
    
    /**
     * @brief Check if event passes generation cuts
     */
    bool PassesGenerationCuts(const SimcEvent& event) const;
    
    /**
     * @brief Validate physical constraints
     */
    bool ValidatePhysics(const SimcEvent& event) const;
    
    // ========================================================================
    // Member Variables
    // ========================================================================
    
    // External dependencies
    std::shared_ptr<RandomGenerator> random_;
    std::shared_ptr<CrossSection> cross_section_;
    std::shared_ptr<SpectrometerOptics> optics_;
    
    // Configuration
    ReactionType reaction_type_{ReactionType::ELASTIC};
    GenerationLimits gen_limits_;
    TargetProperties target_props_;
    SpectrometerSettings spec_electron_;
    SpectrometerSettings spec_hadron_;
    
    // Beam parameters
    double beam_energy_{0.0};          ///< Nominal beam energy (MeV)
    double beam_energy_spread_{0.0};   ///< Energy spread (MeV)
    double beam_energy_vertex_{0.0};   ///< Average at vertex (MeV)
    
    // Flags
    bool use_energy_loss_{true};
    bool use_coulomb_{true};
    bool use_radiative_{true};
    bool initialized_{false};
    
    // Statistics
    unsigned int n_generated_{0};
    unsigned int n_accepted_{0};
    
    // Fermi momentum distribution (for nuclei)
    std::vector<double> pfer_values_;      ///< Momentum values (MeV/c)
    std::vector<double> pfer_probs_;       ///< Cumulative probabilities
    std::vector<double> nprot_theory_;     ///< Proton numbers
    
    // Constants
    static constexpr double kNSigmaMax = 3.0;     ///< Max sigma for gaussian
    static constexpr double kMinWeight = 1e-10;   ///< Minimum weight
};

} // namespace simc

#endif // SIMC_EVENT_GENERATOR_H
