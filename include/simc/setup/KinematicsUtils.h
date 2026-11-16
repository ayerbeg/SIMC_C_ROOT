/**
 * @file
 * @ingroup input
 * @brief Kinematics utilities ported bit-for-bit from the Fortran: physics_angles, spectrometer_angles, loren.
 *
 * Deterministic kinematics utilities of the Fortran SIMC, ported bit-for-bit (task P2.1):
 *   physics_angles       event.f:1572-1616   spectrometer slopes -> lab (theta, phi)
 *   spectrometer_angles  event.f:1618-1647   lab (theta, phi) -> spectrometer slopes
 *   loren                loren.f:1-26        Lorentz boost of a 4-vector
 * Angles in rad; phi = 0 is down, pi/2 the SOS/SHMS/HRS-L side, 3pi/2 the HMS/HRS-R side.
 *
 * Known Fortran behaviour that is reproduced (D-012, P3.F1 findings F1/F7):
 *  * physics_angles assumes phi0 = pi/2 or 3pi/2 (the cos(phi0) term is dropped). For other phi0 it
 *    prints a warning (two lines) and returns the (wrong) formula result.
 *  * spectrometer_angles is ill-conditioned for dy ~ 0: sqrt(1/c^2 - 1 - dx^2) of rounding noise gives
 *    NaN (negative argument) or ~3e-9 instead of 0, and it is not the inverse of physics_angles at
 *    theta0 = 90/120 deg or when the trajectory crosses the beam axis. The NaN is a FORTRAN-BUG (D-028 F2), corrected
 *    unless fortran_compat (the argument is clamped to 0); the ~3e-9 and the other limitations are kept.
 */
#ifndef SIMC_SETUP_KINEMATICSUTILS_H
#define SIMC_SETUP_KINEMATICSUTILS_H

#include <ostream>

namespace simc {

/// @addtogroup input
/// @{

/// Fortran: event.f:physics_angles (l.1572-1616).
/// In: theta0, phi0 (central physics angles of the spectrometer), dx = xptar, dy = yptar.
/// Out: theta, phi of the event. Warnings go to `log` (if non-null) exactly as the Fortran prints them.
/// @param theta0 central polar angle of the spectrometer [rad].
/// @param phi0 central azimuthal angle of the spectrometer [rad] (pi/2 or 3pi/2).
/// @param dx xptar [rad].
/// @param dy yptar [rad].
/// @param[out] theta lab polar angle of the event [rad].
/// @param[out] phi lab azimuthal angle of the event [rad].
/// @param log stream for the warnings (not owned), or nullptr.
void PhysicsAngles(double theta0, double phi0, double dx, double dy, double& theta, double& phi,
                   std::ostream* log = nullptr);

/// Fortran: event.f:spectrometer_angles (l.1618-1647). NOTE: dx, dy are OUTPUTS.
/// In: theta0, phi0, theta, phi. Out: dx = xptar, dy = yptar.
/// `fortran_compat` (D-031): true = sqrt of the (possibly rounding-negative) argument as in the Fortran (NaN dy);
/// false = the argument is clamped to >= 0 (FORTRAN-BUG D-028 F2 / D-012 F1, see the .cpp).
/// @param theta0 central polar angle of the spectrometer [rad].
/// @param phi0 central azimuthal angle of the spectrometer [rad].
/// @param[out] dx xptar [rad].
/// @param[out] dy yptar [rad].
/// @param theta lab polar angle of the event [rad].
/// @param phi lab azimuthal angle of the event [rad].
/// @param fortran_compat true: exact Fortran (NaN possible); false: argument clamped.
void SpectrometerAngles(double theta0, double phi0, double& dx, double& dy, double theta, double phi,
                        bool fortran_compat);

/// Result of loren (energy, momentum components and |p| in the boosted frame).
struct LorenResult {
    double ef1{0} /**< @copydoc pf1 */, pxf{0} /**< @copydoc pf1 */, pyf{0} /**< @copydoc pf1 */, pzf{0} /**< @copydoc pf1 */, pf1{0};  ///< energy, px, py, pz and |p| in the boosted frame (input units).
};

/// Fortran: loren.f:loren (l.1-26). Boost (e,x,y,z) into the frame moving with velocity (bx,by,bz)
/// relative to the initial frame. gam is supplied by the caller (not recomputed from b).
/// @param gam Lorentz factor of the boost.
/// @param bx,by,bz velocity of the new frame (units of c).
/// @param e energy of the 4-vector.
/// @param x,y,z momentum components of the 4-vector (same unit as @p e).
/// @return the boosted 4-vector and its |p|.
LorenResult Loren(double gam, double bx, double by, double bz, double e, double x, double y, double z);

/// @}

}  // namespace simc

#endif  // SIMC_SETUP_KINEMATICSUTILS_H
