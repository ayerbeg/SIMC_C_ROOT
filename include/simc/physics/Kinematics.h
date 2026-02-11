// include/simc/Kinematics.h
// Kinematic calculations for SIMC Monte Carlo

#ifndef SIMC_KINEMATICS_H
#define SIMC_KINEMATICS_H

#include "simc/core/SimcTypes.h"
#include "simc/core/SimcEvent.h"
#include "simc/core/SimcConstants.h"

namespace simc {

/**
 * @class Kinematics
 * @brief Calculate kinematic quantities for electron scattering
 * 
 * This class provides methods to calculate all relevant kinematic
 * quantities from basic measured quantities (beam energy, angles, momenta).
 * 
 * Convention:
 * - Lab frame: beam along z-axis
 * - All energies in MeV
 * - All momenta in MeV/c
 * - All angles in radians
 */
class Kinematics {
public:
    /**
     * @brief Calculate electron scattering kinematics
     * @param Ein Incident beam energy (MeV)
     * @param Ee Scattered electron energy (MeV)
     * @param theta_e Electron scattering angle (rad)
     * @param Ep Detected hadron energy (MeV)
     * @param theta_p Hadron angle (rad)
     * @param phi_pq Out-of-plane angle (rad)
     * @param M_target Target mass (MeV)
     * @param M_hadron Detected hadron mass (MeV)
     */
    static void Calculate(SimcEvent& evt, 
                         double Ein, double Ee, double theta_e,
                         double Ep, double theta_p, double phi_pq,
                         double M_target, double M_hadron);
    
    /**
     * @brief Calculate Q2 from electron kinematics
     * @param Ein Incident energy (MeV)
     * @param Ee Scattered energy (MeV)
     * @param theta Scattering angle (rad)
     * @return Q2 in MeV^2
     */
    static double CalcQ2(double Ein, double Ee, double theta);
    
    /**
     * @brief Calculate energy transfer nu
     * @param Ein Incident energy (MeV)
     * @param Ee Scattered energy (MeV)
     * @return nu in MeV
     */
    static double CalcNu(double Ein, double Ee);
    
    /**
     * @brief Calculate 3-momentum transfer
     * @param Ein Incident energy (MeV)
     * @param Ee Scattered energy (MeV)
     * @param theta Scattering angle (rad)
     * @return |q| in MeV/c
     */
    static double CalcQ(double Ein, double Ee, double theta);
    
    /**
     * @brief Calculate invariant mass W
     * @param M_target Target mass (MeV)
     * @param nu Energy transfer (MeV)
     * @param Q2 Four-momentum transfer squared (MeV^2)
     * @return W in MeV
     */
    static double CalcW(double M_target, double nu, double Q2);
    
    /**
     * @brief Calculate Bjorken x
     * @param Q2 Four-momentum transfer squared (MeV^2)
     * @param M_target Target mass (MeV)
     * @param nu Energy transfer (MeV)
     * @return x_Bjorken
     */
    static double CalcXbj(double Q2, double M_target, double nu);
    
    /**
     * @brief Calculate virtual photon polarization epsilon
     * @param Ein Incident energy (MeV)
     * @param Ee Scattered energy (MeV)
     * @param theta Scattering angle (rad)
     * @param Q2 Four-momentum transfer squared (MeV^2)
     * @return epsilon (0 to 1)
     */
    static double CalcEpsilon(double Ein, double Ee, double theta, double Q2);
    
    /**
     * @brief Calculate missing energy
     * @param M_target Target mass (MeV)
     * @param nu Energy transfer (MeV)
     * @param Ep Detected hadron energy (MeV)
     * @param M_recoil Recoil system mass (MeV)
     * @return E_miss in MeV
     */
    static double CalcEmiss(double M_target, double nu, double Ep, double M_recoil);
    
    /**
     * @brief Calculate missing momentum magnitude
     * @param q 3-momentum transfer magnitude (MeV/c)
     * @param pp Detected hadron momentum (MeV/c)
     * @param theta_pq Angle between p and q (rad)
     * @return P_miss in MeV/c
     */
    static double CalcPmiss(double q, double pp, double theta_pq);
    
    /**
     * @brief Calculate missing momentum vector
     * @param q_vec Virtual photon 3-momentum
     * @param p_vec Detected hadron 3-momentum
     * @return Missing momentum vector
     */
    static Vector3D CalcPmissVec(const Vector3D& q_vec, const Vector3D& p_vec);
    
    /**
     * @brief Calculate angle between hadron and virtual photon
     * @param theta_e Electron angle (rad)
     * @param phi_e Electron azimuth (rad)
     * @param theta_p Hadron angle (rad)
     * @param phi_p Hadron azimuth (rad)
     * @return theta_pq in radians
     */
    static double CalcThetaPQ(double theta_e, double phi_e, 
                              double theta_p, double phi_p);
    
    /**
     * @brief Calculate recoil mass
     * @param M_target Target mass (MeV)
     * @param nu Energy transfer (MeV)
     * @param Ep Detected hadron energy (MeV)
     * @param Pmiss Missing momentum (MeV/c)
     * @return M_recoil in MeV
     */
    static double CalcMrecoil(double M_target, double nu, double Ep, double Pmiss);
    
    /**
     * @brief Convert (theta, phi) to unit vector
     * @param theta Polar angle from z-axis (rad)
     * @param phi Azimuthal angle (rad)
     * @return Unit vector
     */
    static Vector3D AnglesToVector(double theta, double phi);
    
    /**
     * @brief Convert unit vector to (theta, phi)
     * @param vec Unit vector
     * @return {theta, phi} in radians
     */
    static std::pair<double, double> VectorToAngles(const Vector3D& vec);
    
    /**
     * @brief Calculate Mandelstam t
     * @param Q2 Four-momentum transfer squared (MeV^2)
     * @param M_target Target mass (MeV)
     * @param M_hadron Hadron mass (MeV)
     * @param theta_pq Angle between p and q (rad)
     * @param Ep Hadron energy (MeV)
     * @return t in MeV^2
     */
    static double CalcMandelstamT(double Q2, double M_target, double M_hadron,
                                  double theta_pq, double Ep);
    
    /**
     * @brief Calculate minimum t (kinematic limit)
     * @param Q2 Four-momentum transfer squared (MeV^2)
     * @param W Invariant mass (MeV)
     * @param M_target Target mass (MeV)
     * @param M_hadron Hadron mass (MeV)
     * @return t_min in MeV^2
     */
    static double CalcTmin(double Q2, double W, double M_target, double M_hadron);
};

/**
 * @class CoordinateTransform
 * @brief Transform between coordinate systems
 * 
 * Handles transformations between:
 * - Lab frame
 * - Spectrometer (TRANSPORT) coordinates
 * - Physics coordinates
 */
class CoordinateTransform {
public:
    /**
     * @brief Transform from lab frame to spectrometer coordinates
     * @param momentum Lab frame momentum (MeV/c)
     * @param central_momentum Central momentum (MeV/c)
     * @param central_angle Central angle (rad)
     * @return Spectrometer coordinates
     */
    static ArmState LabToSpec(const Vector3D& momentum, 
                              double central_momentum,
                              double central_angle);
    
    /**
     * @brief Transform from spectrometer to lab frame
     * @param spec Spectrometer coordinates
     * @param central_momentum Central momentum (MeV/c)
     * @param central_angle Central angle (rad)
     * @return Lab frame momentum vector
     */
    static Vector3D SpecToLab(const ArmState& spec,
                              double central_momentum,
                              double central_angle);
    
    /**
     * @brief Calculate momentum from spectrometer delta
     * @param delta Momentum deviation (%)
     * @param central_momentum Central momentum (MeV/c)
     * @return Actual momentum (MeV/c)
     */
    static double DeltaToMomentum(double delta, double central_momentum);
    
    /**
     * @brief Calculate delta from momentum
     * @param momentum Actual momentum (MeV/c)
     * @param central_momentum Central momentum (MeV/c)
     * @return Delta (%)
     */
    static double MomentumToDelta(double momentum, double central_momentum);
};

} // namespace simc

#endif // SIMC_KINEMATICS_H
