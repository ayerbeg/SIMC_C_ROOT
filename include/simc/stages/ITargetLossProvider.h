/**
 * @file
 * @ingroup stages
 * @brief simc::ITargetLossProvider: the start-up part of a target model (losses used by the run set-up).
 *
 * The hook between the run set-up (simc::BuildRunSetup) and the target plug-in. Every simc::ITarget implements it, so
 * a new target plug-in also defines the losses that the set-up uses for the limits and the `.hist` (task P7.4). The
 * set-up calls, at exactly the point where the original program makes the call and with its arguments:
 *  - targetInit() = init.f:target_init (l.1-88): targ.L1, L2, X0, X0_cm, the most probable losses Eloss[i].ave (zeroed
 *    if !using_Eloss) and teff[i].ave, targ.Coulomb and Coulomb_constant;
 *  - extremeLosses() = target.f:extreme_trip_thru_target (+ extreme_target_musc), called at init.f:245:
 *    Eloss[i].min/max, teff[i].min/max, musc_max[i], musc_nsig_max.
 *
 * simc::BuildRunSetup applies the limits_init post-processing itself (zeroing Eloss min/max when !using_Eloss,
 * musc_max when !mc_smear, adding the losses to the edges, ...). simc::driver::SimcRun passes a dedicated instance of
 * the pipeline's target (created by name, not configured: the two calls get everything through their arguments).
 * Without a provider, simc::BuildRunSetup uses simc::material::HallCCryoTarget, the only target model of the original.
 */
#ifndef SIMC_STAGES_ITARGETLOSSPROVIDER_H
#define SIMC_STAGES_ITARGETLOSSPROVIDER_H

#include "simc/core/RunSetup.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief Arguments of the target initialisation (internal units).
 * @note Origin: Fortran SIMC init.f:target_init, called at dbase.f:551-552.
 */
struct TargetInitArgs {
    int electron_arm{0}, hadron_arm{0};   ///< Fortran COMMON electron_arm/hadron_arm (window materials per arm)
    bool using_Eloss{false};              ///< energy loss is simulated
    bool using_Coulomb{false};            ///< the Coulomb correction is applied
    double the_cent{0};    ///< [rad] spec%e%theta
    double thp_cent{0};    ///< [rad] spec%p%theta
    double Pp_cent{0};     ///< [MeV/c] spec%p%P
    double Mh2{0};         ///< [MeV^2]
    double Ebeam{0};       ///< [MeV]
    double Pe_cent{0};     ///< [MeV/c] spec%e%P
};

/**
 * @brief Arguments of the extreme-loss computation (internal units: MeV, rad, cm).
 * @note Origin: Fortran SIMC target.f:extreme_trip_thru_target(ebeam, the, thp, pe, pp, z, m) as built by
 *       init.f:limits_init l.231-246.
 */
struct ExtremeTripArgs {
    int electron_arm{0}, hadron_arm{0};   ///< Fortran COMMON electron_arm/hadron_arm
    double Ebeam{0};       ///< beam energy [MeV]
    Cut the;               ///< the_phys min/max [rad] (init.f:232-235)
    Cut thp;               ///< thp_phys min/max [rad] (init.f:236-239)
    Cut pe;                ///< edge%e%E before energy loss [MeV] (init.f:220-223)
    Cut pp;                ///< hadron momentum limits [MeV/c] (init.f:224-226)
    Cut z;                 ///< -/+ targ%length/2 [cm] (init.f:240-241)
    double m{0};           ///< Mh [MeV]
    bool fortran_compat{false};   ///< D-031 run option (RunSetup::fortran_compat): keep the Fortran bugs
};

/**
 * @brief The start-up part of a target model (see the file description).
 *
 * Called by the set-up only, before the event loop, from one thread. Implemented by every simc::ITarget.
 */
class ITargetLossProvider {
public:
    /// Virtual destructor (implementations are deleted through base pointers).
    virtual ~ITargetLossProvider() = default;

    /**
     * @brief The target initialisation: radiation lengths, average losses and the Coulomb correction.
     * @param args The arguments (internal units).
     * @param targ In/out: holds A, Z, N, rho, thick, length, angle, can (after the set-up conversions); must receive
     *        L1, L2, X0, X0_cm, Eloss[i].ave (0 if !using_Eloss), teff[i].ave, Coulomb (ave/min/max) and
     *        Coulomb_constant (i = 0, 1, 2 for beam, electron, hadron).
     * @param messages Informational messages of the original program may be appended.
     * @note Origin: Fortran SIMC init.f:target_init (l.1-88), called at dbase.f:551.
     */
    virtual void targetInit(const TargetInitArgs& args, TargetInfo& targ, std::vector<std::string>& messages) = 0;

    /**
     * @brief The extreme (minimum/maximum) losses over the acceptance, used for the generation limits.
     * @param args The arguments (internal units).
     * @param targ In/out: must receive Eloss[i].min/max, teff[i].min/max, musc_max[i] and musc_nsig_max.
     * @note Origin: Fortran SIMC target.f:extreme_trip_thru_target (called at init.f:245).
     */
    virtual void extremeLosses(const ExtremeTripArgs& args, TargetInfo& targ) = 0;
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_ITARGETLOSSPROVIDER_H
