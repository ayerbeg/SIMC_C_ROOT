// include/simc/EventGenerator.h
// Event generation for SIMC Monte Carlo
// Ported from event.f, jacobians.f

#ifndef SIMC_EVENT_GENERATOR_H
#define SIMC_EVENT_GENERATOR_H

#include "simc/core/SimcTypes.h"
#include "simc/core/SimcEvent.h"
#include "simc/core/SimcConstants.h"
#include "simc/physics/CrossSection.h"
#include "simc/physics/SpectrometerOptics.h"
#include "simc/core/RandomGenerator.h"
#include "simc/physics/Kinematics.h"
#include "simc/core/ConfigManager.h"
#include <memory>
#include <string>
#include <vector>
#include <Rtypes.h>

namespace simc {

// ============================================================================
// Import constants from SimcConstants.h and define k-prefixed versions
// ============================================================================
using constants::Mp;
using constants::Me;
using constants::Mn;
using constants::Mpi;
using constants::Mk;
using constants::pi;
using constants::degrad;

// Define k-prefixed constants for consistency with Fortran naming
constexpr double kProtonMass = Mp;
constexpr double kElectronMass = Me;
constexpr double kNeutronMass = Mn;
constexpr double kPionMass = Mpi;
constexpr double kKaonMass = Mk;
constexpr double kPi = pi;
constexpr double kRadToDeg = degrad;
constexpr double kDegToRad = pi / 180.0;

// ============================================================================
// Generation-Specific Structures
// ============================================================================

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
        double delta_min{-10.0}, delta_max{22.0};      ///< Momentum deviation (%)
        double yptar_min{-0.03}, yptar_max{0.03};      ///< Vertical angle (rad)
        double xptar_min{-0.06}, xptar_max{0.06};      ///< Horizontal angle (rad)
        double E_min{500.0}, E_max{2000.0};            ///< Energy (MeV)
    } electron;
    
    // Hadron arm limits
    struct {
        double delta_min{-10.0}, delta_max{22.0};
        double yptar_min{-0.03}, yptar_max{0.03};
        double xptar_min{-0.06}, xptar_max{0.06};
        double E_min{500.0}, E_max{2000.0};
    } hadron;
    
    // Additional generation constraints
    double sumEgen_min{0.0}, sumEgen_max{0.0};  ///< Sum energy constraint (MeV)
    double Trec_min{0.0}, Trec_max{0.0};        ///< Recoil kinetic energy (MeV)
    
    // Beam parameters
    double xwid{0.01}, ywid{0.01};  ///< Beam widths, 1 sigma (cm)
};

/**
 * @struct TargetProperties
 * @brief Extended target properties for generation
 * 
 * Extends TargetInfo with generation-specific parameters
 */
struct TargetProperties {
    double mass{kProtonMass};     ///< Target mass (MeV)
    int Z{1};                     ///< Atomic number
    int A{1};                     ///< Mass number
    double length{10.0};          ///< Target length (cm)
    double density{0.071};        ///< Density (g/cm³)
    double x_offset{0.0};         ///< X offset (cm)
    double y_offset{0.0};         ///< Y offset (cm)
    double z_offset{0.0};         ///< Z offset (cm)
    double angle{0.0};            ///< Target angle (rad)
    int raster_pattern{0};        ///< Raster: 0=none, 1=bedpost, 2=circular, 3=flat
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

// ============================================================================
// EventGenerator Class
// ============================================================================

/**
 * @class EventGenerator
 * @brief Main event generator for SIMC Monte Carlo
 * 
 * Reference: event.f, jacobians.f from Fortran SIMC
 */
class EventGenerator {
public:
    // ========================================================================
    // Constructor and Configuration
    // ========================================================================
    
    EventGenerator(
        const ConfigManager& config,
        std::shared_ptr<RandomGenerator> random,
        std::shared_ptr<CrossSectionBase> cross_section = nullptr,
        std::shared_ptr<SpectrometerOptics> optics = nullptr
    );
    
    ~EventGenerator() = default;
    
    bool Initialize();
    
    // ========================================================================
    // Main Event Generation
    // ========================================================================
    
    bool GenerateEvent(SimcEvent& event, MainEvent& main);
    bool CompleteEvent(SimcEvent& event, MainEvent& main);
    void GenerateVertex(TargetInfo& target);
    bool GeneratePhaseSpace(SimcEvent& event, MainEvent& main);
    
    // ========================================================================
    // Angle Conversions
    // ========================================================================
    
    void PhysicsAngles(
        Double_t theta0, Double_t phi0,
        Double_t xptar, Double_t yptar,
        Double_t& theta, Double_t& phi
    ) const;
    
    void SpectrometerAngles(
        Double_t theta0, Double_t phi0,
        Double_t theta, Double_t phi,
        Double_t& xptar, Double_t& yptar
    ) const;
    
    // ========================================================================
    // Kinematics Calculations
    // ========================================================================
    
    void CalculateBasicKinematics(SimcEvent& event) const;
    void CalculateMissingMomentum(SimcEvent& event) const;
    void CalculateHadronAngles(SimcEvent& event) const;
    
    // ========================================================================
    // Reaction-Specific Calculations
    // ========================================================================
    
    bool SolveHydrogenElastic(SimcEvent& event) const;
    bool SolveDeuterium(SimcEvent& event, MainEvent& main) const;
    bool SolvePionKaon(SimcEvent& event, const FermiMomentum& pfer) const;
    
    // ========================================================================
    // Fermi Momentum and Spectral Functions
    // ========================================================================
    
    void GenerateFermiMomentum(FermiMomentum& pfer);
    void GenerateMissingEnergy(double pfer, double& Em);
    
    // ========================================================================
    // Jacobians
    // ========================================================================
    
    Double_t CalculateJacobian(const SimcEvent& event, const MainEvent& main) const;
    Double_t CalculateCMJacobian(const SimcEvent& event, const MainEvent& main) const;
    
    // ========================================================================
    // Particle Decays
    // ========================================================================
    
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
    
    void SetReactionType(ReactionType type) { reaction_type_ = type; }
    ReactionType GetReactionType() const { return reaction_type_; }
    void SetGenerationLimits(const GenerationLimits& limits) { gen_limits_ = limits; }
    const GenerationLimits& GetGenerationLimits() const { return gen_limits_; }
    void SetUseEnergyLoss(bool use) { use_energy_loss_ = use; }
    void SetUseCoulomb(bool use) { use_coulomb_ = use; }
    void SetUseRadiative(bool use) { use_radiative_ = use; }
    unsigned int GetNGenerated() const { return n_generated_; }
    unsigned int GetNAccepted() const { return n_accepted_; }
    void ResetStatistics();

private:
    // ========================================================================
    // Helper Methods
    // ========================================================================
    
    void TripThruTarget(
        Int_t particle_id, Double_t z_position, 
        Double_t energy, Double_t theta,
        Double_t& eloss, Double_t& teff
    ) const;
    
    Vector3D UnitVector(double theta, double phi) const;
    bool PassesGenerationCuts(const SimcEvent& event) const;
    bool ValidatePhysics(const SimcEvent& event) const;
    
    // ========================================================================
    // Member Variables
    // ========================================================================
    
    std::shared_ptr<RandomGenerator> random_;
    std::shared_ptr<CrossSectionBase> cross_section_;
    std::shared_ptr<SpectrometerOptics> optics_;
    
    ReactionType reaction_type_{ReactionType::ELASTIC};
    GenerationLimits gen_limits_;
    TargetProperties target_props_;
    SpectrometerSettings spec_electron_;
    SpectrometerSettings spec_hadron_;
    
    double beam_energy_{0.0};
    double beam_energy_spread_{0.0};
    double beam_energy_vertex_{0.0};
    
    bool use_energy_loss_{true};
    bool use_coulomb_{true};
    bool use_radiative_{false};
    bool initialized_{false};
    
    unsigned int n_generated_{0};
    unsigned int n_accepted_{0};
    
    std::vector<double> pfer_values_;
    std::vector<double> pfer_probs_;
    std::vector<double> nprot_theory_;
    
    static constexpr double kNSigmaMax = 3.0;
    static constexpr double kMinWeight = 1e-10;
};

} // namespace simc

#endif // SIMC_EVENT_GENERATOR_H
