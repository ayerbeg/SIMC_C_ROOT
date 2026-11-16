/**
 * @file
 * @ingroup physics
 * @brief Proton (e,e'p) cross sections: sigMott, fofa_best_fit, sigep, deForest (physics_proton.f).
 *
 * Proton (e,e'p) cross sections, ported from simc_gfortran/physics_proton.f.
 * Free functions (deterministic, no state) + IPhysicsModel plug-ins:
 *
 *   "sigep"                 elastic H(e,e'p)  (event.f:1450-1452, doing_hyd_elast)
 *   "deforest"              deForest, variant from SimcParameters::deForest_flag (event.f:1454-1456)
 *   "deforest_cc1"          deForest with deForest_flag =  0 (sigcc1)
 *   "deforest_cc2"          deForest with deForest_flag = +1 (sigcc2)
 *   "deforest_cc1_onshell"  deForest with deForest_flag = -1 (sigcc1 on-shell)
 *
 * Numerics follow D-010: double literals (gfortran -fdefault-real-8), the Fortran evaluation
 * order (left-associative, no reassociation), std::pow for real exponents (Q**3. etc.),
 * dot products summed x+y+z. Validated at L1 against reference/functions/{sigmott,
 * fofa_best_fit,sigep,deforest} (tests/l1/test_l1_physics_proton.cpp).
 */
#ifndef SIMC_PHYSICS_PROTON_PROTONCROSSSECTIONS_H
#define SIMC_PHYSICS_PROTON_PROTONCROSSSECTIONS_H

#include <iosfwd>

#include "simc/core/Event.h"

namespace simc {
namespace proton {

/// @addtogroup physics
/// @{

/// Fortran: physics_proton.f:sigMott (l.176-191). Mott cross section, microbarn/sr.
/// e0 [MeV], theta [rad], Q2 [MeV^2].
/// @param e0 beam energy [MeV].
/// @param theta scattering angle [rad].
/// @param Q2 four-momentum transfer squared [MeV^2].
/// @return the Mott cross section [microbarn/sr].
double sigMott(double e0, double theta, double Q2);

/// Fortran: physics_proton.f:fofa_best_fit (l.137-171). Bosted fit (PRC 51, 409, Eqs. 4-5).
/// qsquar [fm^-2] is NEGATIVE for spacelike momentum transfer (callers pass -Q2/hbarc**2).
/// @param qsquar q^2 [fm^-2] (negative for spacelike).
/// @param[out] GE electric form factor.
/// @param[out] GM magnetic form factor.
void fofa_best_fit(double qsquar, double& GE, double& GM);

/// Fortran: physics_proton.f:sigep (l.1-21). Elastic e-p cross section, microbarn/sr.
/// Reads vertex.Q2, vertex.e.E, vertex.e.theta, vertex.Ein only.
/// @param vertex the vertex.
/// @return the cross section [microbarn/sr].
double sigep(const Event& vertex);

/// COMMON state read by deForest (simulate.inc): deForest_flag, Mh2, nevent (warnings only).
struct DeForestState {
    int deForest_flag{0};     ///< -1 sigcc1 on-shell, 0 sigcc1, +1 sigcc2
    double Mh2{0};            ///< hadron mass squared [MeV^2] (dbase.f:224: Mh*Mh, Mh = Mp for (e,e'p))
    long nevent{0};           ///< event counter printed in the warnings (simc.f nevent)
    std::ostream* log{nullptr};  ///< where WARNINGs go; nullptr = std::cout (Fortran unit 6)
    bool fortran_compat{false};  ///< D-031: true = the cos_phi warning can never print, as in the Fortran (D-014)
};

/// Fortran: physics_proton.f:deForest (l.25-133). Off-shell e-N cross section
/// K*sigma_eN in microbarn*MeV^2/sr^2. Does not modify `ev`.
/// @param ev the event.
/// @param state the COMMON state.
/// @return K*sigma_eN [microbarn*MeV^2/sr^2].
double deForest(const Event& ev, const DeForestState& state);

/// @}

}  // namespace proton
}  // namespace simc

#endif  // SIMC_PHYSICS_PROTON_PROTONCROSSSECTIONS_H
