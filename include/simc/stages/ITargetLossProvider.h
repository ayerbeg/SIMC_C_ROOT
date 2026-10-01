// include/simc/stages/ITargetLossProvider.h
//
// The start-up part of a target model: the hook between the run set-up (src/setup, BuildRunSetup) and the target
// plug-in. Every ITarget implements it (ITarget.h), so a new target plug-in also defines the losses that the set-up
// uses for the limits and the .hist (P7.4). The set-up calls, at exactly the point where the Fortran makes the call
// and with the Fortran arguments:
//
//   init.f:target_init (l.1-88)         whole routine: targ%L1, L2, X0, X0_cm, the most probable
//                                        losses targ%Eloss(i)%ave (zeroed if !using_Eloss) and
//                                        targ%teff(i)%ave, targ%Coulomb and Coulomb_constant
//   target.f:extreme_trip_thru_target   targ%Eloss(i)%min/max, targ%teff(i)%min/max,
//     (+ extreme_target_musc)            targ%musc_max(i), targ%musc_nsig_max  (called at init.f:245)
//
// BuildRunSetup applies the limits_init post-processing itself (zeroing Eloss min/max when
// !using_Eloss, musc_max when !mc_smear, adding the losses to the edges ...).
// driver::SimcRun passes a dedicated instance of the pipeline's target (created by name, not configured: the two
// calls below get everything through their arguments). BuildRunSetup without a provider uses material::HallCCryoTarget,
// the only target model of the Fortran.
#ifndef SIMC_STAGES_ITARGETLOSSPROVIDER_H
#define SIMC_STAGES_ITARGETLOSSPROVIDER_H

#include "simc/core/RunSetup.h"

namespace simc {

/// Arguments of init.f:target_init (dbase.f:551-552). Fortran units.
struct TargetInitArgs {
    int electron_arm{0}, hadron_arm{0};   ///< Fortran COMMON electron_arm/hadron_arm (window materials per arm)
    bool using_Eloss{false};
    bool using_Coulomb{false};
    double the_cent{0};    ///< [rad] spec%e%theta
    double thp_cent{0};    ///< [rad] spec%p%theta
    double Pp_cent{0};     ///< [MeV/c] spec%p%P
    double Mh2{0};         ///< [MeV^2]
    double Ebeam{0};       ///< [MeV]
    double Pe_cent{0};     ///< [MeV/c] spec%e%P
};

/// Arguments of target.f:extreme_trip_thru_target(ebeam, the, thp, pe, pp, z, m) as built by
/// init.f:limits_init l.231-246. Fortran units (MeV, rad, cm).
struct ExtremeTripArgs {
    int electron_arm{0}, hadron_arm{0};   ///< Fortran COMMON electron_arm/hadron_arm
    double Ebeam{0};
    Cut the;               ///< the_phys min/max [rad] (init.f:232-235)
    Cut thp;               ///< thp_phys min/max [rad] (init.f:236-239)
    Cut pe;                ///< edge%e%E before energy loss [MeV] (init.f:220-223)
    Cut pp;                ///< hadron momentum limits [MeV/c] (init.f:224-226)
    Cut z;                 ///< -/+ targ%length/2 [cm] (init.f:240-241)
    double m{0};           ///< Mh [MeV]
};

class ITargetLossProvider {
public:
    virtual ~ITargetLossProvider() = default;

    /// Fortran: init.f:target_init (l.1-88), called at dbase.f:551. `targ` holds A, Z, N, rho, thick,
    /// length, angle, can (dbase.f conversions). Must set L1, L2, X0, X0_cm, Eloss[i].ave (0 if
    /// !using_Eloss), teff[i].ave, Coulomb (ave/min/max) and Coulomb_constant (i = 0,1,2 for Fortran 1,2,3).
    /// Informational Fortran messages may be appended to `messages`.
    virtual void targetInit(const TargetInitArgs& args, TargetInfo& targ, std::vector<std::string>& messages) = 0;

    /// Fortran: target.f:extreme_trip_thru_target (called at init.f:245).
    /// Must set targ.Eloss[i].min/max, targ.teff[i].min/max, targ.musc_max[i], targ.musc_nsig_max.
    virtual void extremeLosses(const ExtremeTripArgs& args, TargetInfo& targ) = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_ITARGETLOSSPROVIDER_H
