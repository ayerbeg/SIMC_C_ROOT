// include/simc/SimcEvent.h
// Main event structure for SIMC Monte Carlo

#ifndef SIMC_EVENT_H
#define SIMC_EVENT_H

#include "SimcTypes.h"
#include "SimcConstants.h"
#include <Rtypes.h>  // ROOT types
#include <TObject.h>

namespace simc {

/**
 * @class SimcEvent
 * @brief Main event structure containing all kinematic information
 * 
 * This class stores both generated (vertex) and reconstructed quantities
 * for a single simulated event. It uses ROOT-compatible types to allow
 * direct storage in TTrees.
 * 
 * Convention:
 * - All energies in MeV
 * - All momenta in MeV/c  
 * - All angles in radians
 * - All positions in cm
 */
class SimcEvent : public TObject {
public:
    // ========================================================================
    // Constructors and Destructor
    // ========================================================================
    SimcEvent();
    virtual ~SimcEvent() = default;
    
    // ========================================================================
    // Beam and Kinematic Quantities
    // ========================================================================
    Double_t Ein;           ///< Incident beam energy (MeV)
    Double_t nu;            ///< Energy transfer (MeV)
    Double_t q;             ///< 3-momentum transfer magnitude (MeV/c)
    Double_t Q2;            ///< 4-momentum transfer squared (MeV^2)
    Double_t W;             ///< Invariant mass (MeV)
    Double_t xbj;           ///< Bjorken x
    Double_t epsilon;       ///< Virtual photon polarization parameter
    
    // ========================================================================
    // Missing Quantities
    // ========================================================================
    Double_t Em;            ///< Missing energy (MeV)
    Double_t Pm;            ///< Missing momentum magnitude (MeV/c)
    Double_t Emiss;         ///< Calculated missing energy (MeV)
    Double_t Pmiss;         ///< Calculated missing momentum (MeV/c)
    
    // Missing momentum components (lab frame)
    Double_t Pmx;           ///< Missing momentum x-component (MeV/c)
    Double_t Pmy;           ///< Missing momentum y-component (MeV/c)
    Double_t Pmz;           ///< Missing momentum z-component (MeV/c)
    
    // Missing momentum decomposition
    Double_t PmPar;         ///< Parallel component (along q) (MeV/c)
    Double_t PmPer;         ///< Perpendicular component (MeV/c)
    Double_t PmOop;         ///< Out-of-plane component (MeV/c)
    
    // ========================================================================
    // Recoil System
    // ========================================================================
    Double_t Mrec;          ///< Recoil system mass (MeV)
    Double_t Trec;          ///< Recoil kinetic energy (MeV)
    
    // ========================================================================
    // Angles and Kinematics
    // ========================================================================
    Double_t theta_pq;      ///< Angle between hadron and q (rad)
    Double_t phi_pq;        ///< Out-of-plane angle (rad)
    Double_t theta_tarq;    ///< Angle between target pol and q (rad)
    Double_t phi_targ;      ///< Target polarization angle (rad)
    Double_t beta;          ///< Beta angle for target polarization (rad)
    Double_t phi_s;         ///< Sivers angle (rad)
    Double_t phi_c;         ///< Collins angle (rad)
    
    // Semi-inclusive specific
    Double_t zhad;          ///< Hadron energy fraction z = E_h/nu
    Double_t pt2;           ///< Transverse momentum squared (GeV^2)
    
    // ========================================================================
    // Electron Arm
    // ========================================================================
    // Spectrometer coordinates
    Double_t e_delta;       ///< Momentum deviation (%)
    Double_t e_xptar;       ///< Horizontal angle (rad)
    Double_t e_yptar;       ///< Vertical angle (rad)
    Double_t e_z;           ///< Position along beamline (cm)
    
    // Physics coordinates
    Double_t e_theta;       ///< Scattering angle (rad)
    Double_t e_phi;         ///< Azimuthal angle (rad)
    Double_t e_E;           ///< Energy (MeV)
    Double_t e_P;           ///< Momentum (MeV/c)
    
    // ========================================================================
    // Hadron Arm
    // ========================================================================
    // Spectrometer coordinates
    Double_t p_delta;       ///< Momentum deviation (%)
    Double_t p_xptar;       ///< Horizontal angle (rad)
    Double_t p_yptar;       ///< Vertical angle (rad)
    Double_t p_z;           ///< Position along beamline (cm)
    
    // Physics coordinates
    Double_t p_theta;       ///< Scattering angle (rad)
    Double_t p_phi;         ///< Azimuthal angle (rad)
    Double_t p_E;           ///< Energy (MeV)
    Double_t p_P;           ///< Momentum (MeV/c)
    
    // ========================================================================
    // Unit Vectors
    // ========================================================================
    Double_t ue_x, ue_y, ue_z;  ///< Electron direction unit vector
    Double_t up_x, up_y, up_z;  ///< Hadron direction unit vector
    Double_t uq_x, uq_y, uq_z;  ///< Virtual photon direction unit vector
    
    // ========================================================================
    // Methods
    // ========================================================================
    
    /**
     * @brief Reset all values to zero
     */
    void Clear(Option_t* option = "") override;
    
    /**
     * @brief Print event information
     */
    void Print(Option_t* option = "") const override;
    
    /**
     * @brief Set electron state from ArmState
     */
    void SetElectron(const ArmState& arm);
    
    /**
     * @brief Set hadron state from ArmState
     */
    void SetHadron(const ArmState& arm);
    
    /**
     * @brief Set unit vectors from Vector3D
     */
    void SetUnitVectors(const Vector3D& ue, const Vector3D& up, const Vector3D& uq);
    
    /**
     * @brief Get electron state as ArmState
     */
    ArmState GetElectronState() const;
    
    /**
     * @brief Get hadron state as ArmState
     */
    ArmState GetHadronState() const;
    
    /**
     * @brief Get electron unit vector
     */
    Vector3D GetElectronVector() const;
    
    /**
     * @brief Get hadron unit vector
     */
    Vector3D GetHadronVector() const;
    
    /**
     * @brief Get photon unit vector
     */
    Vector3D GetPhotonVector() const;
    
    // ROOT dictionary generation
    ClassDef(SimcEvent, 1)
};

/**
 * @struct MainEvent
 * @brief Additional event information including weights and metadata
 * 
 * This structure contains quantities that are calculated during
 * the Monte Carlo process but are not part of the core event kinematics.
 */
struct MainEvent {
    // ========================================================================
    // Target Information
    // ========================================================================
    TargetInfo target;          ///< Target vertex and properties
    
    // ========================================================================
    // Spectrometer Quantities
    // ========================================================================
    ArmState SP_electron;       ///< Electron in spectrometer coordinates
    ArmState SP_hadron;         ///< Hadron in spectrometer coordinates
    ArmState RECON_electron;    ///< Reconstructed electron
    ArmState RECON_hadron;      ///< Reconstructed hadron
    
    // ========================================================================
    // Focal Plane
    // ========================================================================
    FocalPlaneState FP_electron;    ///< Electron focal plane state
    FocalPlaneState FP_hadron;      ///< Hadron focal plane state
    
    // ========================================================================
    // Weights and Cross Sections
    // ========================================================================
    Double_t weight{1.0};           ///< Total event weight
    Double_t jacobian{1.0};         ///< Phase space jacobian
    Double_t gen_weight{1.0};       ///< Generation weight
    Double_t SF_weight{1.0};        ///< Spectral function weight
    Double_t sigcc{1.0};            ///< Cross section (ub/MeV^n/sr^m)
    Double_t sigcc_recon{1.0};      ///< Reconstructed cross section
    Double_t sigcent{1.0};          ///< Central cross section
    
    // ========================================================================
    // Energy Shifts
    // ========================================================================
    Double_t Ein_shift{0.0};        ///< Beam energy shift
    Double_t Ee_shift{0.0};         ///< Electron energy shift
    
    // ========================================================================
    // Additional Kinematics
    // ========================================================================
    Double_t W{0.0};                ///< Invariant mass (MeV)
    Double_t epsilon{0.0};          ///< Virtual photon polarization
    Double_t theta_pq{0.0};         ///< Hadron-q angle (rad)
    Double_t phi_pq{0.0};           ///< Out-of-plane angle (rad)
    Double_t t{0.0};                ///< Mandelstam t (GeV^2)
    Double_t tmin{0.0};             ///< Minimum t (GeV^2)
    Double_t q2{0.0};               ///< Q^2 (GeV^2)
    
    // ========================================================================
    // Flags
    // ========================================================================
    Bool_t success{false};          ///< Event passed all cuts
    
    /**
     * @brief Calculate total weight
     */
    Double_t GetTotalWeight() const {
        return weight * jacobian * gen_weight * SF_weight * sigcc;
    }
    
    /**
     * @brief Reset all values
     */
    void Clear() {
        weight = jacobian = gen_weight = SF_weight = 1.0;
        sigcc = sigcc_recon = sigcent = 1.0;
        Ein_shift = Ee_shift = 0.0;
        W = epsilon = theta_pq = phi_pq = 0.0;
        t = tmin = q2 = 0.0;
        success = false;
    }
};

} // namespace simc

#endif // SIMC_EVENT_H
