/**
 * @file
 * @ingroup driver
 * @brief Angle between the polarized-target field and each spectrometer arm (simc.f:120-155, task P15.2b).
 *
 * simc.f:120-155 (task P15.2b): the angle between the polarized-target field and each spectrometer arm, passed to
 * trgInit (optics::TargetField::init) in degrees.
 */
#ifndef SIMC_DRIVER_TARGETFIELDANGLES_H
#define SIMC_DRIVER_TARGETFIELDANGLES_H

#include <string>
#include <vector>

namespace simc {
namespace driver {

/// @addtogroup driver
/// @{

/// @brief Result of TargetFieldArmAngles(): the field angle of each arm and the log of the computation.
///
/// simc.f:120-155: the angle between the target field and each arm (radians), from targ_Bangle/targ_Bphi and the
/// central angles (all in radians). `log` receives the Fortran write(*,*) lines. The angle of an arm stays 0 (the
/// SAVEd, zero-initialised local of the Fortran main program) when the error branch is taken (T3: no STOP).
/// compat = true keeps the Fortran bugs: T1 (simc.f:127/141 opposite-side test without the outer abs) and T2
/// (simc.f:139 p-arm same-side sign for targ_Bangle < theta_p). Default (corrected): abs(abs(dphi)-180) < .01 and the
/// e-arm formula for the p arm.
struct ArmFieldAngles {
    double earm{0.0}, parm{0.0};   ///< ang_targ_earm, ang_targ_parm [rad]
    std::vector<std::string> log;  ///< Fortran write(*,*) lines produced while computing the angles (one string per line).
};
/// @brief Compute the angle between the target field and each arm.
/// @param targ_Bangle field polar angle [rad].
/// @param targ_Bphi field azimuthal angle [rad].
/// @param e_theta central polar angle of the electron arm [rad].
/// @param e_phi central azimuthal angle of the electron arm [rad].
/// @param p_theta central polar angle of the hadron arm [rad].
/// @param p_phi central azimuthal angle of the hadron arm [rad].
/// @param compat true: keep the Fortran bugs T1/T2 described for ArmFieldAngles; false: corrected formulas.
/// @return the two angles [rad] and the log lines.
ArmFieldAngles TargetFieldArmAngles(double targ_Bangle, double targ_Bphi, double e_theta, double e_phi, double p_theta,
                         double p_phi, bool compat);

/// @}
}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_TARGETFIELDANGLES_H
