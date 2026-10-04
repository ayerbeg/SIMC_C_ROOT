// include/simc/material/HallCCryoTarget.h
//
// ITarget "hallc_cryo": the standard SIMC target (cryo cells can 1/2/3 and solid slabs) with the
// hard-coded Hall A/C window, air, kevlar and mylar thicknesses per spectrometer.
//   trip()                   Fortran: target.f:trip_thru_target (l.1-306)
//   extremeTripThruTarget()  Fortran: target.f:extreme_trip_thru_target (l.310-544)
//   targetMusc()             Fortran: target.f:target_musc (l.548-577), 3.5 sigma truncation
//   extremeTargetMusc()      Fortran: target.f:extreme_target_musc (entry, l.581-588)
//   targetInit()             Fortran: init.f:target_init (l.1-88)
//   configure()              Fortran: dbase.f:451-468 (targ conversions/checks), 239-244 (angles), 547
//   generateZ()              Fortran: event.f:185 main%target%z = (0.5-grnd())*targ%length+targ%zoffset
//
// State (Fortran COMMON / static locals):
//   * targ (COMMON /target_info/) -> targ();  electron_arm/hadron_arm -> setArms().
//   * The energy losses go through an IMaterialModel (default: own SimcMaterial "simc_material"),
//     which carries the enerloss_new static state (numerr, x). In SIMC only target.f calls
//     enerloss_new, so one target instance reproduces the Fortran process state.
//   Markers (docs/FORTRAN_DIFFERENCES.md): FORTRAN-QUIRK = documented Fortran behaviour, FORTRAN-DIFF-CANDIDATE = kept,
//   open question, FORTRAN-BUG = corrected unless fortran_compat.
//   * FORTRAN-DIFF-CANDIDATE(D-012 F4): th_can of the pudding can (can 2) is a static local. When
//     |z_can/R| > 1 + 1e-6 it is NOT updated (only a warning is printed) and the previous call's value
//     is used. Reproduced (member thCan_, initially 0).
//   * FORTRAN-BUG(D-012 F4), path not ported: electron_arm/hadron_arm 7/8 (calorimeter) leave s_air, s_kevlar,
//     s_mylar and forward_path stale. Out of scope until P14: trip() throws std::invalid_argument
//     for the calorimeter arms and for codes that are not spectrometers (core/SpectrometerId.h, D-032).
//   * FORTRAN-QUIRK: `goto (10,20,30) narm` falls through to label 10 for narm outside 1..3,
//     i.e. such a call is treated as the incoming beam. Reproduced.
//   * FORTRAN-BUG(D-017), corrected unless fortran_compat: extreme_trip_thru_target, hadron, can 2/3 (target.f:501) uses
//     the ELECTRON angle the%min to find zz (the hadron angle thp%min in the corrected mode).
//   * FORTRAN-DIFF-CANDIDATE(D-017): the beer-can corner in extreme_trip_thru_target uses a 1.25" radius
//     (target.f:368,370,376,456,458,464), trip_thru_target uses 1.325" (l.108,228). Kept.
#ifndef SIMC_MATERIAL_HALLCCRYOTARGET_H
#define SIMC_MATERIAL_HALLCCRYOTARGET_H

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "simc/material/TargetInfo.h"
#include "simc/stages/IMaterialModel.h"
#include "simc/stages/ITarget.h"

namespace simc {
namespace material {

/// Gaussian truncation of target_musc (target.f:555 parameter nsig_max = 3.5).
constexpr double kTargetMuscNsig = 3.5;

/// Arguments of init.f:target_init as dbase.f:551 passes them (filled by configure()).
struct TargetInitInputs {
    bool using_Eloss{false};
    bool using_Coulomb{false};   ///< already forced .false. for nint(Z) == 1 (dbase.f:547)
    double the_cent{0};          ///< spec%e%theta [rad] = abs(deg)/degrad (dbase.f:239)
    double thp_cent{0};          ///< spec%p%theta [rad] (dbase.f:242)
    double Pp_cent{0};           ///< spec%p%P [MeV/c]
    double Mh2{0};               ///< Mh*Mh [MeV^2]
    double Ebeam{0};             ///< [MeV]
    double Pe_cent{0};           ///< spec%e%P [MeV/c]
};

class HallCCryoTarget : public ITarget {
public:
    HallCCryoTarget();

    std::string name() const override { return "hallc_cryo"; }
    /// dbase.f conversions of the targ% input and the target_init arguments.
    /// JSON (all optional): { "Mh": MeV (default Mp; the reaction's hadron mass, dbase.f:215-226),
    ///                        "material": IMaterialModel name (default "simc_material"), "material_config": {...},
    ///                        "verbose": bool (default true; false = no Fortran write(6,*) messages) }
    void configure(const Json& config, const SimcParameters& params) override;
    /// Fortran: dbase.f:551 call target_init(...) with the configured inputs.
    void init(const StageContext& ctx) override;

    // ---- ITarget ---------------------------------------------------------------------------
    double generateZ(IRandom& rng) override;
    Traversal tripThruTarget(Arm arm, double zpos, double energy, double theta, double mass, LossMode mode,
                             IRandom& rng) override;
    /// target_musc (3.5 sigma).
    ScatteringAngles multipleScattering(double p, double beta, double teff, IRandom& rng) override;

    // ---- ITargetLossProvider (the set-up's start-up losses; arms taken from the arguments) ---------------------
    /// Fortran: init.f:target_init on the RunSetup target description (copied in, results copied back).
    void targetInit(const TargetInitArgs& args, simc::TargetInfo& targ, std::vector<std::string>& messages) override;
    /// Fortran: target.f:extreme_trip_thru_target (+ extreme_target_musc) on the RunSetup target description.
    void extremeLosses(const ExtremeTripArgs& args, simc::TargetInfo& targ) override;

    // ---- Fortran routines (same argument order) --------------------------------------------
    /// Fortran: target.f:trip_thru_target(narm,zpos,energy,theta,Eloss,radlen,mass,typeflag)
    void trip(int narm, double zpos, double energy, double theta, double& Eloss, double& radlen, double mass,
              int typeflag, IRandom& rng);
    /// Fortran: target.f:extreme_trip_thru_target(ebeam,the,thp,pe,pp,z,m). No RNG (typeflags 2/3).
    /// Writes targ().Eloss[i].min/max, teff[i].min/max, musc_max[i], musc_nsig_max.
    /// `fortran_compat` (D-031): true keeps the Fortran bug of target.f:501 (FORTRAN-BUG D-017, see the .cpp).
    void extremeTripThruTarget(double ebeam, const Limits& the, const Limits& thp, const Limits& pe,
                               const Limits& pp, const Limits& z, double m, bool fortran_compat);
    /// Fortran: target.f:target_musc(p,beta,teff,dangles)
    void targetMusc(double p, double beta, double teff, double dangles[2], IRandom& rng);
    /// Fortran: target.f:extreme_target_musc(p,beta,teff,dangle,r) (entry of target_musc)
    void extremeTargetMusc(double p, double beta, double teff, double& dangle, double& r) const;
    /// Fortran: init.f:target_init(using_Eloss,using_Coulomb,the_cent,thp_cent,Pp_cent,Mh2,ebeam,Pe_cent)
    void targetInit(bool using_Eloss, bool using_Coulomb, double the_cent, double thp_cent, double Pp_cent,
                    double Mh2, double ebeam, double Pe_cent);

    // ---- state -----------------------------------------------------------------------------
    TargetInfo& targ() { return targ_; }
    const TargetInfo& targ() const { return targ_; }
    void setArms(int electron_arm, int hadron_arm) {
        electron_arm_ = electron_arm;
        hadron_arm_ = hadron_arm;
    }
    int electronArm() const { return electron_arm_; }
    int hadronArm() const { return hadron_arm_; }
    const TargetInitInputs& initInputs() const { return init_; }
    TargetInitInputs& initInputs() { return init_; }

    void setMaterialModel(std::unique_ptr<IMaterialModel> m) { material_ = std::move(m); }
    IMaterialModel& materialModel() { return *material_; }

    /// Where the Fortran write(6,*) warnings go (default std::cout; nullptr = silent).
    /// Also forwarded to the material model if it is a SimcMaterial.
    void setLog(std::ostream* os);

    double thCan() const { return thCan_; }   ///< static th_can of the pudding can (F4)

private:
    TargetInfo targ_;
    int electron_arm_{0}, hadron_arm_{0};
    TargetInitInputs init_;
    std::unique_ptr<IMaterialModel> material_;   ///< owned; default SimcMaterial ("simc_material")
    std::ostream* log_;
    double thCan_{0};   ///< target.f: th_can (static local, -fno-automatic)
};

}  // namespace material
}  // namespace simc

#endif  // SIMC_MATERIAL_HALLCCRYOTARGET_H
