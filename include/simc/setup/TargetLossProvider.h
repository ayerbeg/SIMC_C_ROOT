// include/simc/setup/TargetLossProvider.h
//
// Hook between the run setup (task P2.1) and the energy-loss / multiple-scattering code
// (task P3.1, src/material). The setup code needs, before the event loop, the pieces of the
// Fortran that call trip_thru_target / target_musc:
//
//   init.f:target_init l.44-58        most probable losses targ%Eloss(i)%ave, targ%teff(i)%ave
//                                     (trip_thru_target typeflag 4, i = 1 beam, 2 e', 3 hadron)
//   target.f:extreme_trip_thru_target targ%Eloss(i)%min/max, targ%teff(i)%min/max,
//     (+ extreme_target_musc)          targ%musc_max(i), targ%musc_nsig_max  (called at init.f:245)
//
// BuildRunSetup calls the provider at exactly the point where the Fortran makes the call, with the
// Fortran arguments, and then applies the Fortran post-processing itself (zeroing Eloss when
// !using_Eloss, musc_max when !mc_smear, adding the losses to the edges ...). Everything else of
// target_init (L1, L2, X0, X0_cm, Coulomb) needs no energy loss and is computed by the builder.
//
// P3.1 / P2.2 integration: implement this interface on top of the material/target code, e.g.
//   class MaterialTargetLosses : public ITargetLossProvider { ... uses TripThruTarget(...) ... };
// Without a provider the builder leaves those fields at 0 and sets
// RunSetup::target_losses_from_provider = false (fields derived from them are then incomplete
// unless using_Eloss = mc_smear = false, where the Fortran zeroes them anyway; teff always needs it).
#ifndef SIMC_SETUP_TARGETLOSSPROVIDER_H
#define SIMC_SETUP_TARGETLOSSPROVIDER_H

#include "simc/core/RunSetup.h"

namespace simc {

/// Arguments of init.f:target_init that the trip_thru_target calls use (l.44-58). Fortran units.
struct TargetInitArgs {
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

    /// Fortran: init.f:target_init l.44-55 (the three trip_thru_target(...,4) calls).
    /// Must set targ.Eloss[i].ave and targ.teff[i].ave (i = 0,1,2 for Fortran 1,2,3). `targ` already
    /// holds A, Z, N, rho, thick, length, angle, can, X0, X0_cm, L1, L2 (dbase.f + target_init l.19-43).
    /// The builder zeroes Eloss ave afterwards if !using_Eloss (init.f:56-60).
    virtual void mostProbableLosses(const TargetInitArgs& args, TargetInfo& targ) = 0;

    /// Fortran: target.f:extreme_trip_thru_target (called at init.f:245).
    /// Must set targ.Eloss[i].min/max, targ.teff[i].min/max, targ.musc_max[i], targ.musc_nsig_max.
    virtual void extremeLosses(const ExtremeTripArgs& args, TargetInfo& targ) = 0;
};

}  // namespace simc

#endif  // SIMC_SETUP_TARGETLOSSPROVIDER_H
