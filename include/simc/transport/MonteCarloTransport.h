// include/simc/MonteCarloTransport.h
// Phase 5c.1: Monte Carlo Transport - Core functionality
// Multiple scattering, energy loss, and basic acceptance
// Full spectrometer MCs (Phase 5c.2) to be added later

#ifndef SIMC_MONTE_CARLO_TRANSPORT_H
#define SIMC_MONTE_CARLO_TRANSPORT_H

#include "simc/core/SimcEvent.h"
#include "simc/core/ConfigManager.h"
#include "simc/core/RandomGenerator.h"
#include <memory>

namespace simc {

/**
 * @brief Monte Carlo transport through target and spectrometer acceptance
 * 
 * Phase 5c.1 Implementation:
 * - Multiple scattering in target
 * - Energy loss corrections
 * - Calculate hadron spectrometer angles for elastic events
 * - Simple acceptance checking (delta/angle cuts)
 * 
 * Phase 5c.2 (Future):
 * - Full spectrometer magnetic transport (mc_hms, mc_sos, mc_shms)
 * - Focal plane calculations
 * - Detailed aperture checking
 */
class MonteCarloTransport {
public:
    /**
     * @brief Constructor
     * @param config Configuration manager
     * @param random Random number generator
     */
    MonteCarloTransport(const ConfigManager& config,
                       std::shared_ptr<RandomGenerator> random);
    
    /**
     * @brief Main transport function - applies all corrections
     * @param event The event to transport (modified in place)
     * @param main Main event information (modified in place)
     * @return true if event passes acceptance, false otherwise
     */
    bool Transport(SimcEvent& event, MainEvent& main);
    
    /**
     * @brief Calculate reconstructed quantities after transport
     * @param event Event with transported quantities
     * 
     * For elastic events, this calculates p_xptar, p_yptar from 
     * the proton's physics angles (p_theta, p_phi)
     */
    void CalculateReconstructed(SimcEvent& event);

private:
    // Configuration
    struct TransportConfig {
        bool use_energy_loss;
        bool use_multiple_scattering;
        
        // Spectrometer configuration
        struct SpectrometerConfig {
            double P;         // Central momentum (MeV/c)
            double theta;     // Central angle (rad)
            double phi;       // Azimuthal angle (rad)
            double cos_th;    // cos(theta)
            double sin_th;    // sin(theta)
            
            // Acceptance limits (for Phase 5c.1 simple cuts)
            double delta_min, delta_max;  // %
            double xptar_min, xptar_max;  // rad
            double yptar_min, yptar_max;  // rad
            
            // Offsets
            struct {
                double x, y, z;           // cm
                double xptar, yptar;      // rad
            } offset;
        } electron, hadron;
        
        // Target configuration
        double length;        // cm
        double angle;         // rad
        
        // Raster correction
        bool correct_raster;
        
    } config_;
    
    std::shared_ptr<RandomGenerator> random_;
    
    // === PHASE 5c.1 IMPLEMENTATIONS ===
    
    /**
     * @brief Apply multiple scattering to beam
     * @param energy Beam energy (MeV)
     * @param beta Particle velocity (v/c)
     * @param teff Target effective thickness (radiation lengths)
     * @param dang_out [out] Scattering angles (yptar, xptar) in rad
     */
    void ApplyBeamMultipleScattering(double energy, double beta, 
                                     double teff, double dang_out[2]);
    
    /**
     * @brief Apply multiple scattering to scattered particle
     * @param momentum Particle momentum (MeV/c)
     * @param beta Particle velocity (v/c)
     * @param teff Target effective thickness (radiation lengths)
     * @param dangles [out] Scattering angles (yptar, xptar) in rad
     */
    void ApplyMultipleScattering(double momentum, double beta, 
                                double teff, double dangles[2]);
    
    /**
     * @brief Apply energy loss corrections to get spectrometer quantities
     * @param event Event with true quantities
     * @param main Main event structure with target info
     */
    void ApplyEnergyLoss(SimcEvent& event, MainEvent& main);
    
    /**
     * @brief Calculate spectrometer angles from physics angles
     * @param theta0 Spectrometer central angle (rad)
     * @param phi0 Spectrometer azimuthal angle (rad)
     * @param theta Physics polar angle (rad)
     * @param phi Physics azimuthal angle (rad)
     * @param xptar [out] Spectrometer horizontal angle (rad)
     * @param yptar [out] Spectrometer vertical angle (rad)
     */
    void SpectrometerAngles(double theta0, double phi0,
                           double theta, double phi,
                           double& xptar, double& yptar) const;
    
    /**
     * @brief Simple acceptance check (Phase 5c.1)
     * @param event Event to check
     * @param main Main event structure
     * @return true if event passes acceptance cuts
     * 
     * Phase 5c.1: Simple delta/angle cuts
     * Phase 5c.2: Full spectrometer acceptance with apertures
     */
    bool CheckAcceptance(const SimcEvent& event, const MainEvent& main) const;
    
    // === PHASE 5c.2 PLACEHOLDERS (Future) ===
    
    /**
     * @brief Full spectrometer MC transport (Phase 5c.2)
     * @return true if event survives transport
     * 
     * NOT IMPLEMENTED IN PHASE 5c.1
     * Will call mc_hms, mc_sos, mc_shms, etc.
     */
    bool TransportThroughSpectrometer(SimcEvent& event, MainEvent& main);
    
    /**
     * @brief Calculate focal plane quantities (Phase 5c.2)
     * 
     * NOT IMPLEMENTED IN PHASE 5c.1
     */
    void CalculateFocalPlane(const SimcEvent& event, MainEvent& main);
};

} // namespace simc

#endif // SIMC_MONTE_CARLO_TRANSPORT_H
