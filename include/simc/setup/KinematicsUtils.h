// include/simc/setup/KinematicsUtils.h
//
// Deterministic kinematics utilities of the Fortran SIMC, ported bit-for-bit (task P2.1):
//   physics_angles       event.f:1572-1616   spectrometer slopes -> lab (theta, phi)
//   spectrometer_angles  event.f:1618-1647   lab (theta, phi) -> spectrometer slopes
//   loren                loren.f:1-26        Lorentz boost of a 4-vector
// Angles in rad; phi = 0 is down, pi/2 the SOS/SHMS/HRS-L side, 3pi/2 the HMS/HRS-R side.
//
// Known Fortran behaviour that is reproduced (D-012, P3.F1 findings F1/F7):
//  * physics_angles assumes phi0 = pi/2 or 3pi/2 (the cos(phi0) term is dropped). For other phi0 it
//    prints a warning (two lines) and returns the (wrong) formula result.
//  * spectrometer_angles is ill-conditioned for dy ~ 0: sqrt(1/c^2 - 1 - dx^2) of rounding noise gives
//    NaN (negative argument) or ~3e-9 instead of 0, and it is not the inverse of physics_angles at
//    theta0 = 90/120 deg or when the trajectory crosses the beam axis.
#ifndef SIMC_SETUP_KINEMATICSUTILS_H
#define SIMC_SETUP_KINEMATICSUTILS_H

#include <ostream>

namespace simc {

/// Fortran: event.f:physics_angles (l.1572-1616).
/// In: theta0, phi0 (central physics angles of the spectrometer), dx = xptar, dy = yptar.
/// Out: theta, phi of the event. Warnings go to `log` (if non-null) exactly as the Fortran prints them.
void PhysicsAngles(double theta0, double phi0, double dx, double dy, double& theta, double& phi,
                   std::ostream* log = nullptr);

/// Fortran: event.f:spectrometer_angles (l.1618-1647). NOTE: dx, dy are OUTPUTS.
/// In: theta0, phi0, theta, phi. Out: dx = xptar, dy = yptar.
void SpectrometerAngles(double theta0, double phi0, double& dx, double& dy, double theta, double phi);

/// Result of loren (energy, momentum components and |p| in the boosted frame).
struct LorenResult {
    double ef1{0}, pxf{0}, pyf{0}, pzf{0}, pf1{0};
};

/// Fortran: loren.f:loren (l.1-26). Boost (e,x,y,z) into the frame moving with velocity (bx,by,bz)
/// relative to the initial frame. gam is supplied by the caller (not recomputed from b).
LorenResult Loren(double gam, double bx, double by, double bz, double e, double x, double y, double z);

}  // namespace simc

#endif  // SIMC_SETUP_KINEMATICSUTILS_H
