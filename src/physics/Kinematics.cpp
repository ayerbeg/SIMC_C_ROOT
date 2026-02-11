// src/physics/Kinematics.cpp
// Implementation of kinematic calculations

#include "simc/physics/Kinematics.h"
#include <cmath>

namespace simc {
using namespace constants;

// ============================================================================
// Main Calculation
// ============================================================================
void Kinematics::Calculate(SimcEvent& evt, 
                          double Ein, double Ee, double theta_e,
                          double Ep, double theta_p, double phi_pq,
                          double M_target, double M_hadron) {
    // Store beam energy
    evt.Ein = Ein;
    
    // Basic electron quantities
    double Pe = Ee;  // Electron is essentially massless
    evt.e_E = Ee;
    evt.e_P = Pe;
    evt.e_theta = theta_e;
    
    // Virtual photon kinematics
    evt.nu = CalcNu(Ein, Ee);
    evt.Q2 = CalcQ2(Ein, Ee, theta_e);
    evt.q = CalcQ(Ein, Ee, theta_e);
    evt.W = CalcW(M_target, evt.nu, evt.Q2);
    evt.xbj = CalcXbj(evt.Q2, M_target, evt.nu);
    evt.epsilon = CalcEpsilon(Ein, Ee, theta_e, evt.Q2);
    
    // Hadron quantities
    double Pp = std::sqrt(Ep*Ep - M_hadron*M_hadron);
    evt.p_E = Ep;
    evt.p_P = Pp;
    evt.p_theta = theta_p;
    evt.phi_pq = phi_pq;
    
    // Calculate angle between hadron and q
    evt.theta_pq = std::acos((evt.q*evt.q + Pp*Pp - evt.nu*evt.nu - Ep*Ep + 2*evt.nu*Ep) / 
                             (2.0 * evt.q * Pp));
    
    // Missing quantities
    evt.Pmiss = CalcPmiss(evt.q, Pp, evt.theta_pq);
    evt.Emiss = CalcEmiss(M_target, evt.nu, Ep, M_target);  // Assuming nucleon recoil
    
    // Missing momentum vector (in lab frame with beam along z)
    Vector3D q_vec = AnglesToVector(theta_e, 0.0);
    q_vec = q_vec * evt.q;
    
    Vector3D p_vec = AnglesToVector(theta_p, phi_pq);
    p_vec = p_vec * Pp;
    
    Vector3D Pm_vec = CalcPmissVec(q_vec, p_vec);
    evt.Pmx = Pm_vec.x;
    evt.Pmy = Pm_vec.y;
    evt.Pmz = Pm_vec.z;
    evt.Pm = Pm_vec.Magnitude();
    
    // Parallel and perpendicular components
    Vector3D q_unit = q_vec.Normalized();
    evt.PmPar = Pm_vec.Dot(q_unit);
    Vector3D Pm_perp = Pm_vec - q_unit * evt.PmPar;
    evt.PmPer = Pm_perp.Magnitude();
    
    // Recoil mass
    evt.Mrec = CalcMrecoil(M_target, evt.nu, Ep, evt.Pmiss);
    evt.Trec = std::sqrt(evt.Mrec*evt.Mrec + evt.Pmiss*evt.Pmiss) - evt.Mrec;
    
    // Unit vectors
    Vector3D ue = AnglesToVector(theta_e, 0.0);
    Vector3D up = AnglesToVector(theta_p, phi_pq);
    Vector3D uq = q_vec.Normalized();
    evt.SetUnitVectors(ue, up, uq);
}

// ============================================================================
// Q2 Calculation
// ============================================================================
double Kinematics::CalcQ2(double Ein, double Ee, double theta) {
    return 2.0 * Ein * Ee * (1.0 - std::cos(theta));
}

// ============================================================================
// Energy Transfer
// ============================================================================
double Kinematics::CalcNu(double Ein, double Ee) {
    return Ein - Ee;
}

// ============================================================================
// 3-Momentum Transfer
// ============================================================================
double Kinematics::CalcQ(double Ein, double Ee, double theta) {
    double Q2 = CalcQ2(Ein, Ee, theta);
    double nu = CalcNu(Ein, Ee);
    return std::sqrt(Q2 + nu*nu);
}

// ============================================================================
// Invariant Mass W
// ============================================================================
double Kinematics::CalcW(double M_target, double nu, double Q2) {
    double W2 = M_target*M_target + 2.0*M_target*nu - Q2;
    return (W2 > 0.0) ? std::sqrt(W2) : 0.0;
}

// ============================================================================
// Bjorken x
// ============================================================================
double Kinematics::CalcXbj(double Q2, double M_target, double nu) {
    return Q2 / (2.0 * M_target * nu);
}

// ============================================================================
// Virtual Photon Polarization
// ============================================================================
double Kinematics::CalcEpsilon(double Ein, double Ee, double theta, double Q2) {
    double nu = CalcNu(Ein, Ee);
    double tan2_half = std::tan(theta/2.0) * std::tan(theta/2.0);
    return 1.0 / (1.0 + 2.0 * (1.0 + nu*nu/Q2) * tan2_half);
}

// ============================================================================
// Missing Energy
// ============================================================================
double Kinematics::CalcEmiss(double M_target, double nu, double Ep, double M_recoil) {
    return M_target + nu - Ep - M_recoil;
}

// ============================================================================
// Missing Momentum
// ============================================================================
double Kinematics::CalcPmiss(double q, double pp, double theta_pq) {
    return std::sqrt(q*q + pp*pp - 2.0*q*pp*std::cos(theta_pq));
}

Vector3D Kinematics::CalcPmissVec(const Vector3D& q_vec, const Vector3D& p_vec) {
    return q_vec - p_vec;
}

// ============================================================================
// Angle Between Hadron and q
// ============================================================================
double Kinematics::CalcThetaPQ(double theta_e, double phi_e, 
                               double theta_p, double phi_p) {
    // Convert to unit vectors
    Vector3D ue = AnglesToVector(theta_e, phi_e);
    Vector3D up = AnglesToVector(theta_p, phi_p);
    
    // q is along incident - scattered electron
    Vector3D q = AnglesToVector(0.0, 0.0) - ue;  // Beam along z
    q = q.Normalized();
    
    // Angle between p and q
    double cos_theta = up.Dot(q);
    return std::acos(cos_theta);
}

// ============================================================================
// Recoil Mass
// ============================================================================
double Kinematics::CalcMrecoil(double M_target, double nu, double Ep, double Pmiss) {
    double Erec = M_target + nu - Ep;
    double M2 = Erec*Erec - Pmiss*Pmiss;
    return (M2 > 0.0) ? std::sqrt(M2) : 0.0;
}

// ============================================================================
// Angle/Vector Conversions
// ============================================================================
Vector3D Kinematics::AnglesToVector(double theta, double phi) {
    return {
        std::sin(theta) * std::cos(phi),
        std::sin(theta) * std::sin(phi),
        std::cos(theta)
    };
}

std::pair<double, double> Kinematics::VectorToAngles(const Vector3D& vec) {
    Vector3D unit = vec.Normalized();
    double theta = std::acos(unit.z);
    double phi = std::atan2(unit.y, unit.x);
    return {theta, phi};
}

// ============================================================================
// Mandelstam Variables
// ============================================================================
double Kinematics::CalcMandelstamT(double Q2, double M_target, double M_hadron,
                                   double theta_pq, double Ep) {
    double Pp = std::sqrt(Ep*Ep - M_hadron*M_hadron);
    double nu = Ep - M_target;  // Approximation
    double q = std::sqrt(Q2 + nu*nu);
    
    // t = (q - p)^2
    return M_hadron*M_hadron - Q2 - 2.0*q*Pp*std::cos(theta_pq) + 2.0*nu*Ep;
}

double Kinematics::CalcTmin(double Q2, double W, double M_target, double M_hadron) {
    double s = M_target*M_target;
    double term1 = (s - M_hadron*M_hadron - Q2) * (s - W*W);
    double term2 = 2.0 * s * (s + Q2);
    return (term1 - term2) / (2.0 * s);
}

// ============================================================================
// Coordinate Transformations
// ============================================================================
ArmState CoordinateTransform::LabToSpec(const Vector3D& momentum, 
                                        double central_momentum,
                                        double central_angle) {
    ArmState spec;
    
    // Calculate momentum magnitude
    double P = momentum.Magnitude();
    spec.P = P;
    
    // Calculate delta
    spec.delta = MomentumToDelta(P, central_momentum);
    
    // Rotate to spectrometer frame
    double cos_theta = std::cos(central_angle);
    double sin_theta = std::sin(central_angle);
    
    Vector3D p_spec;
    p_spec.x = momentum.x * cos_theta + momentum.z * sin_theta;
    p_spec.y = momentum.y;
    p_spec.z = -momentum.x * sin_theta + momentum.z * cos_theta;
    
    // Calculate angles in spectrometer frame
    spec.xptar = p_spec.x / p_spec.z;
    spec.yptar = p_spec.y / p_spec.z;
    
    // Calculate physics angles
    auto [theta, phi] = Kinematics::VectorToAngles(momentum);
    spec.theta = theta;
    spec.phi = phi;
    
    return spec;
}

Vector3D CoordinateTransform::SpecToLab(const ArmState& spec,
                                        double central_momentum,
                                        double central_angle) {
    // Get actual momentum
    double P = DeltaToMomentum(spec.delta, central_momentum);
    
    // Build momentum in spectrometer frame
    double pz_spec = P / std::sqrt(1.0 + spec.xptar*spec.xptar + spec.yptar*spec.yptar);
    Vector3D p_spec{spec.xptar * pz_spec, spec.yptar * pz_spec, pz_spec};
    
    // Rotate to lab frame
    double cos_theta = std::cos(central_angle);
    double sin_theta = std::sin(central_angle);
    
    Vector3D p_lab;
    p_lab.x = p_spec.x * cos_theta - p_spec.z * sin_theta;
    p_lab.y = p_spec.y;
    p_lab.z = p_spec.x * sin_theta + p_spec.z * cos_theta;
    
    return p_lab;
}

double CoordinateTransform::DeltaToMomentum(double delta, double central_momentum) {
    return central_momentum * (1.0 + delta / 100.0);
}

double CoordinateTransform::MomentumToDelta(double momentum, double central_momentum) {
    return 100.0 * (momentum - central_momentum) / central_momentum;
}

} // namespace simc
