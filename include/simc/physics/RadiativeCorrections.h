// include/simc/RadiativeCorrections.h
// Radiative corrections for electron scattering

#ifndef SIMC_RADIATIVE_CORRECTIONS_H
#define SIMC_RADIATIVE_CORRECTIONS_H

#include "simc/core/SimcEvent.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/RandomGenerator.h"
#include "simc/physics/Kinematics.h"
#include <array>
#include <string>

namespace simc {

// Forward declarations
class Bremsstrahlung;
class ExternalRadiation;

/**
 * @struct RadConfig
 * @brief Configuration parameters for radiative corrections
 * 
 * Replaces Fortran COMMON blocks /radccom/ and /bremcom/.
 * Contains all parameters needed to configure the radiation generation.
 */
struct RadConfig {
    // ========================================================================
    // Radiation Modes
    // ========================================================================
    
    /// Radiation mode flag
    /// - 0: Best formulas, peaked basis (ntail, Egamma)
    /// - 1: BASICRAD only, peaked basis  
    /// - 2: BASICRAD only, (Egamma1,2,3) basis, no overlap
    /// - 3: BASICRAD only, all 3 tails simultaneously
    int rad_flag = 0;
    
    /// External radiation mode
    /// - 0: Use defaults
    /// - 1: Basic external (φ=1)
    /// - 2: Basic × φ correction
    /// - 3: Friedrich approximation
    int extrad_flag = 0;
    
    /// Internal correction mode
    /// - 0: Use brem() for exact calculation
    /// - 1: Skip exact bremsstrahlung (use basic formulas only)
    int intcor_mode = 0;
    
    /// Force single tail (0=no forcing, 1/2/3=force that tail)
    int one_tail = 0;
    
    /// Which tail to radiate this event (0=decide randomly, 1/2/3=specific tail)
    int ntail = 0;
    
    // ========================================================================
    // Flags
    // ========================================================================
    
    /// Use off-shell proton radiation (bremos instead of brem)
    bool use_offshell_rad = false;
    
    /// Radiate proton this event (tail 3)
    bool rad_proton_this_ev = false;
    
    /// Use hardwired parameters (bypass initialization)
    bool hardwired_rad = false;
    
    /// Which tails are active [tail 1, tail 2, tail 3]
    std::array<bool, 3> doing_tail = {false, false, false};
    
    /// Use exponentiation (exp(-δ)) vs linear (1-δ)
    bool use_expon = true;
    
    // ========================================================================
    // Generation Parameters
    // ========================================================================
    
    /// Exponents for power-law generation
    /// g[0] = both tails, g[1] = tail 1, g[2] = tail 2, g[3] = tail 3, g[4] = total
    std::array<double, 5> g = {0, 0, 0, 0, 0};
    
    /// Normalization constants
    std::array<double, 5> c = {0, 0, 0, 0, 0};
    
    /// Internal exponents
    std::array<double, 4> c_int = {0, 0, 0, 0};
    std::array<double, 4> g_int = {0, 0, 0, 0};
    
    /// External exponents  
    std::array<double, 4> c_ext = {0, 0, 0, 0};
    std::array<double, 4> g_ext = {0, 0, 0, 0};
    
    /// Single g_int, g_ext values (for some modes)
    double g_int_val = 0.0;
    double g_ext_val = 0.0;
    
    // ========================================================================
    // Tail Parameters
    // ========================================================================
    
    /// Probability each tail radiates [tail 1, tail 2, tail 3]
    std::array<double, 3> frac = {0, 0, 0};
    
    /// Lambda parameters (extended peaking approximation)
    std::array<double, 3> lambda = {0, 0, 0};
    
    /// Radiation thickness [0]=electron, [1]=proton
    std::array<double, 2> bt = {0, 0};
    
    // ========================================================================
    // Target Properties
    // ========================================================================
    
    /// Target thickness parameter (for external bremsstrahlung)
    double etatzai = 0.0;
    
    /// Related thickness parameter
    double etta = 0.0;
    
    /// Energy resolution FWHM (MeV)
    double fwhm = 0.0;
    
    // ========================================================================
    // Energy Limits
    // ========================================================================
    
    /// Minimum photon energies for each tail (MeV)
    std::array<double, 3> Egamma_min = {0, 0, 0};
    
    /// Maximum photon energies for each tail (MeV)
    std::array<double, 3> Egamma_max = {0, 0, 0};
    
    /// Resolution cutoff limit (MeV)
    double Egamma_res_limit = 0.0;
    
    /// Total maximum photon energy (MeV)
    double Egamma_tot_max = 0.0;
    
    /// Individual tail maxima (MeV)
    double Egamma1_max = 0.0;
    double Egamma2_max = 0.0;
    double Egamma3_max = 0.0;
    
    // ========================================================================
    // Correction Factors
    // ========================================================================
    
    /// Hard correction factor
    double hardcorfac = 1.0;
    
    // ========================================================================
    // Bremsstrahlung Control
    // ========================================================================
    
    /// Produce debug output
    bool produce_output = false;
    
    /// Use exp(-δ) instead of (1-δ)
    bool exponentiate = true;
    
    /// Include hard corrections
    bool include_hard = true;
    
    /// Calculate full Spence function (vs approximation)
    bool calculate_spence = false;
};

/**
 * @struct RadState
 * @brief Per-event radiation state
 * 
 * Tracks the state of radiation generation for a single event.
 */
struct RadState {
    /// Actual photon energies used for each tail (MeV)
    std::array<double, 3> Egamma_used = {0, 0, 0};
    
    /// Total radiative weight for this event
    double rad_weight = 1.0;
    
    /// Generation successful?
    bool success = false;
    
    /// Weight from basicrad() for each tail
    std::array<double, 3> basicrad_weight = {1, 1, 1};
    
    /// Reciprocal of generating function value
    std::array<double, 3> basicrad_val_reciprocal = {0, 0, 0};
    
    /// Reset state for new event
    void Reset() {
        Egamma_used.fill(0.0);
        rad_weight = 1.0;
        success = false;
        basicrad_weight.fill(1.0);
        basicrad_val_reciprocal.fill(0.0);
    }
};

/**
 * @class RadiativeCorrections
 * @brief Main class for radiative correction generation
 * 
 * Implements Mo & Tsai formalism for radiative corrections in 
 * electron scattering. Handles internal and external bremsstrahlung
 * for incoming electron, outgoing electron, and outgoing hadron.
 * 
 * References:
 * - Mo & Tsai, Rev. Mod. Phys. 41, 205 (1969)
 * - Stein et al., Phys. Rev. D 12, 1884 (1975)
 * - Maximon & Tjon, Phys. Rev. C 62, 054320 (2000)
 */
class RadiativeCorrections {
public:
    /**
     * @brief Constructor
     * @param rng Random number generator
     */
    explicit RadiativeCorrections(RandomGenerator& rng);
    
    /**
     * @brief Destructor
     */
    ~RadiativeCorrections();
    
    /**
     * @brief Initialize with configuration
     * @param config Radiation configuration
     */
    void Initialize(const RadConfig& config);
    
    /**
     * @brief Generate radiative corrections for an event
     * @param vertex Event at vertex (INPUT/OUTPUT: radiation subtracted)
     * @param orig Event before radiation (OUTPUT: filled with original values)
     * @param target_mass Target mass (MeV)
     * @param hadron_mass Detected hadron mass (MeV)
     * @param gen_weight Generation weight (INPUT/OUTPUT: multiplied by rad weight)
     * @param doing_heavy Is this A(e,e'p) with A>2?
     * @param doing_deuterium Is this D(e,e'p)?
     * @param doing_hyd_elast Is this H(e,e'p) elastic?
     * @param doing_pion Is this pion production?
     * @param doing_kaon Is this kaon production?
     * @param doing_delta Is this delta production?
     * @return True if radiation generation successful
     * 
     * MODIFIES:
     * - vertex: Energies reduced by radiated photons
     * - orig: Filled with original energies before radiation
     * - gen_weight: Multiplied by radiative weight
     */
    bool Generate(SimcEvent& vertex, SimcEvent& orig,
                  double target_mass, double hadron_mass,
                  double& gen_weight,
                  bool doing_heavy, bool doing_deuterium, bool doing_hyd_elast,
                  bool doing_pion, bool doing_kaon, bool doing_delta);
    
    /**
     * @brief Get current radiation state
     * @return Const reference to radiation state
     */
    const RadState& GetState() const { return state_; }
    
    /**
     * @brief Get configuration
     * @return Const reference to configuration
     */
    const RadConfig& GetConfig() const { return config_; }
    
    /**
     * @brief Set debug output
     * @param enable Enable debug output
     */
    void SetDebugOutput(bool enable) { config_.produce_output = enable; }

private:
    // Configuration and state
    RadConfig config_;
    RadState state_;
    RandomGenerator& rng_;
    
    // Helper classes (will be implemented)
    class Bremsstrahlung* brem_;
    class ExternalRadiation* extrad_;
    
    // Private helper methods
    
    /**
     * @brief Generate photon energy using power law
     * @param itail Which tail (1,2,3)
     * @param Egamma_lo Minimum energy (MeV)
     * @param Egamma_hi Maximum energy (MeV)
     * @param Egamma Generated energy (OUTPUT)
     * @param weight Generation weight (OUTPUT)
     * @param val_reciprocal 1/generating_function (OUTPUT)
     * @return True if successful
     */
    bool BasicRad(int itail, double Egamma_lo, double Egamma_hi,
                  double& Egamma, double& weight, double& val_reciprocal);
    
    /**
     * @brief Calculate precise radiative weight (peaking approximation)
     * @param vertex Event kinematics
     * @param Egamma Photon energy (MeV)
     * @param emin Minimum energy (MeV)
     * @param emax Maximum energy (MeV)
     * @param basicrad_val_reciprocal From BasicRad
     * @param basicrad_weight From BasicRad
     * @return Radiative weight
     */
    double PeakedRadWeight(const SimcEvent& vertex, double Egamma,
                          double emin, double emax,
                          double basicrad_val_reciprocal,
                          double basicrad_weight);
    
    /**
     * @brief Recalculate event kinematics after radiation
     * @param vertex Event to recalculate (modified)
     * @param target_mass Target mass (MeV)
     * @param hadron_mass Hadron mass (MeV)
     * @param doing_hyd_elast Is this H(e,e'p) elastic?
     * @param doing_deuterium Is this D(e,e'p)?
     * @param doing_pion Is this pion production?
     * @param doing_kaon Is this kaon production?
     * @param doing_delta Is this delta production?
     * @return True if kinematics are physical
     */
    bool CompleteEvent(SimcEvent& vertex,
                      double target_mass,
                      double hadron_mass,
                      bool doing_hyd_elast,
                      bool doing_deuterium,
                      bool doing_pion,
                      bool doing_kaon,
                      bool doing_delta);
    
    /**
     * @brief Calculate Schwinger correction
     * @param Ecutoff Cutoff energy (MeV)
     * @param vertex Event kinematics
     * @param target_mass Target mass (MeV)
     * @param include_hard Include hard corrections?
     * @param dsoft Soft correction (OUTPUT)
     * @param dhard Hard correction (OUTPUT)
     * @return Schwinger correction factor
     */
    double Schwinger(double Ecutoff, const SimcEvent& vertex,
                    double target_mass, bool include_hard,
                    double& dsoft, double& dhard);
    
    /**
     * @brief Calculate lambda parameter (extended peaking approx)
     * @param itail Which tail (1,2,3)
     * @param plus_flag Include plus term?
     * @param doing_proton Proton in final state?
     * @param e1 Incident electron energy (MeV)
     * @param e2 Scattered electron energy (MeV)
     * @param e3 Hadron energy (MeV)
     * @param p3 Hadron momentum (MeV/c)
     * @param theta Scattering angle (rad)
     * @return Lambda parameter
     */
    double LambdaDave(int itail, bool plus_flag, bool doing_proton,
                     double e1, double e2, double e3, double p3, double theta);
};

} // namespace simc

#endif // SIMC_RADIATIVE_CORRECTIONS_H
