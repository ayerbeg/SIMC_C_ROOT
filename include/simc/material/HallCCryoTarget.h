/**
 * @file
 * @ingroup material
 * @brief simc::material::HallCCryoTarget: ITarget "hallc_cryo", the standard SIMC target (cryo cells and solid slabs).
 *
 * ITarget "hallc_cryo": the standard SIMC target (cryo cells can 1/2/3 and solid slabs) with the
 * hard-coded Hall A/C window, air, kevlar and mylar thicknesses per spectrometer.
 *   trip()                   Fortran: target.f:trip_thru_target (l.1-306)
 *   extremeTripThruTarget()  Fortran: target.f:extreme_trip_thru_target (l.310-544)
 *   targetMusc()             Fortran: target.f:target_musc (l.548-577), 3.5 sigma truncation
 *   extremeTargetMusc()      Fortran: target.f:extreme_target_musc (entry, l.581-588)
 *   targetInit()             Fortran: init.f:target_init (l.1-88)
 *   configure()              Fortran: dbase.f:451-468 (targ conversions/checks), 239-244 (angles), 547
 *   generateZ()              Fortran: event.f:185 main%target%z = (0.5-grnd())*targ%length+targ%zoffset
 *
 * State (Fortran COMMON / static locals):
 *   * targ (COMMON /target_info/) -> targ();  electron_arm/hadron_arm -> setArms().
 *   * The energy losses go through an IMaterialModel (default: own SimcMaterial "simc_material"),
 *     which carries the enerloss_new static state (numerr, x). In SIMC only target.f calls
 *     enerloss_new, so one target instance reproduces the Fortran process state.
 *   Markers (docs/translation/markdown/FORTRAN_DIFFERENCES.md): FORTRAN-QUIRK = documented Fortran behaviour, FORTRAN-DIFF-CANDIDATE = kept,
 *   open question, FORTRAN-BUG = corrected unless fortran_compat.
 *   * FORTRAN-DIFF-CANDIDATE(D-012 F4): th_can of the pudding can (can 2) is a static local. When
 *     |z_can/R| > 1 + 1e-6 it is NOT updated (only a warning is printed) and the previous call's value
 *     is used. Reproduced (member thCan_, initially 0).
 *   * FORTRAN-DIFF-CANDIDATE(D-012 F4, P14.3): target.f:73-102/189-222 have no window branch for the calorimeter
 *     arms 7/8, so s_air, s_kevlar, s_mylar and forward_path keep the values of the previous narm 2/3 call (static locals,
 *     -fno-automatic; init.f:53-56 target_init makes the first calls), s_Al is reset to 0 (l.20). Reproduced in both
 *     modes (members win_*): which windows (and which path in the target) the calorimeter arm should get needs the NPS
 *     geometry. Codes that are not spectrometers still throw std::invalid_argument (core/SpectrometerId.h, D-032).
 *   * FORTRAN-QUIRK: `goto (10,20,30) narm` falls through to label 10 for narm outside 1..3,
 *     i.e. such a call is treated as the incoming beam. Reproduced.
 *   * FORTRAN-BUG(D-017), corrected unless fortran_compat: extreme_trip_thru_target, hadron, can 2/3 (target.f:501) uses
 *     the ELECTRON angle the%min to find zz (the hadron angle thp%min in the corrected mode).
 *   * FORTRAN-DIFF-CANDIDATE(D-017): the beer-can corner in extreme_trip_thru_target uses a 1.25" radius
 *     (target.f:368,370,376,456,458,464), trip_thru_target uses 1.325" (l.108,228). Kept.
 */
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

/// @addtogroup material
/// @{

/// Gaussian truncation of target_musc (target.f:555 parameter nsig_max = 3.5).
constexpr double kTargetMuscNsig = 3.5;

/// Arguments of init.f:target_init as dbase.f:551 passes them (filled by configure()).
struct TargetInitInputs {
    bool using_Eloss{false};  ///< energy loss in the target is simulated.
    bool using_Coulomb{false};   ///< already forced .false. for nint(Z) == 1 (dbase.f:547)
    double the_cent{0};          ///< spec%e%theta [rad] = abs(deg)/degrad (dbase.f:239)
    double thp_cent{0};          ///< spec%p%theta [rad] (dbase.f:242)
    double Pp_cent{0};           ///< spec%p%P [MeV/c]
    double Mh2{0};               ///< Mh*Mh [MeV^2]
    double Ebeam{0};             ///< [MeV]
    double Pe_cent{0};           ///< spec%e%P [MeV/c]
};

/// @brief ITarget "hallc_cryo": the Hall C cryotarget (target.f, init.f:target_init) with energy loss, multiple
///        scattering and the start-up loss limits (ITargetLossProvider).
///
/// Holds Fortran SAVEd state (th_can, window path lengths). Owns its material model. Not thread-safe.
class HallCCryoTarget : public ITarget {
public:
    /// Creates the target with the default material model ("simc_material").
    HallCCryoTarget();

    std::string name() const override { return "hallc_cryo"; }
    /// dbase.f conversions of the targ% input and the target_init arguments.
    /// JSON (all optional): { "Mh": MeV (default Mp; the reaction's hadron mass, dbase.f:215-226),
    ///                        "material": IMaterialModel name (default "simc_material"), "material_config": {...},
    ///                        "verbose": bool (default true; false = no Fortran write(6,*) messages) }
    /// @param config the JSON config (see above).
    /// @param params the run parameters.
    void configure(const Json& config, const SimcParameters& params) override;
    /// Fortran: dbase.f:551 call target_init(...) with the configured inputs.
    /// @param ctx the stage context.
    void init(const StageContext& ctx) override;

    // ---- ITarget ---------------------------------------------------------------------------
    double generateZ(IRandom& rng) override;
    Traversal tripThruTarget(Arm arm, double zpos, double energy, double theta, double mass, LossMode mode,
                             IRandom& rng) override;
    /// target_musc (3.5 sigma).
    /// @param p momentum [MeV/c].
    /// @param beta velocity (v/c).
    /// @param teff effective thickness [radiation lengths].
    /// @param rng the random generator.
    /// @return the two scattering angles [rad].
    ScatteringAngles multipleScattering(double p, double beta, double teff, IRandom& rng) override;

    // ---- ITargetLossProvider (the set-up's start-up losses; arms taken from the arguments) ---------------------
    /// Fortran: init.f:target_init on the RunSetup target description (copied in, results copied back).
    /// @param args the target_init arguments.
    /// @param[in,out] targ the RunSetup target description.
    /// @param[out] messages the Fortran messages.
    void targetInit(const TargetInitArgs& args, simc::TargetInfo& targ, std::vector<std::string>& messages) override;
    /// Fortran: target.f:extreme_trip_thru_target (+ extreme_target_musc) on the RunSetup target description.
    /// @param args the extreme_trip_thru_target arguments.
    /// @param[in,out] targ the RunSetup target description.
    void extremeLosses(const ExtremeTripArgs& args, simc::TargetInfo& targ) override;

    // ---- Fortran routines (same argument order) --------------------------------------------
    /// Fortran: target.f:trip_thru_target(narm,zpos,energy,theta,Eloss,radlen,mass,typeflag)
    /// @param narm arm number (1 = beam, 2 = electron, 3 = hadron).
    /// @param zpos vertex z [cm].
    /// @param energy particle energy [MeV].
    /// @param theta particle angle [rad].
    /// @param[out] Eloss energy loss [MeV].
    /// @param[out] radlen thickness traversed [radiation lengths].
    /// @param mass particle mass [MeV].
    /// @param typeflag 1 = random loss, 2 = minimum, 3 = maximum, 4 = most probable.
    /// @param rng the random generator.
    void trip(int narm, double zpos, double energy, double theta, double& Eloss, double& radlen, double mass,
              int typeflag, IRandom& rng);
    /// Fortran: target.f:extreme_trip_thru_target(ebeam,the,thp,pe,pp,z,m). No RNG (typeflags 2/3).
    /// Writes targ().Eloss[i].min/max, teff[i].min/max, musc_max[i], musc_nsig_max.
    /// `fortran_compat` (D-031): true keeps the Fortran bug of target.f:501 (FORTRAN-BUG D-017, see the .cpp).
    /// @param ebeam beam energy [MeV].
    /// @param the electron angle range [rad].
    /// @param thp hadron angle range [rad].
    /// @param pe electron momentum range [MeV/c].
    /// @param pp hadron momentum range [MeV/c].
    /// @param z vertex z range [cm].
    /// @param m hadron mass [MeV].
    /// @param fortran_compat true: keep the Fortran bug D-017.
    void extremeTripThruTarget(double ebeam, const Limits& the, const Limits& thp, const Limits& pe,
                               const Limits& pp, const Limits& z, double m, bool fortran_compat);
    /// Fortran: target.f:target_musc(p,beta,teff,dangles)
    /// @param p momentum [MeV/c].
    /// @param beta velocity (v/c).
    /// @param teff effective thickness [radiation lengths].
    /// @param[out] dangles the two scattering angles [rad].
    /// @param rng the random generator.
    void targetMusc(double p, double beta, double teff, double dangles[2], IRandom& rng);
    /// Fortran: target.f:extreme_target_musc(p,beta,teff,dangle,r) (entry of target_musc)
    /// @param p momentum [MeV/c].
    /// @param beta velocity (v/c).
    /// @param teff effective thickness [radiation lengths].
    /// @param[out] dangle the maximum scattering angle [rad].
    /// @param[out] r the width of the angle distribution [rad].
    void extremeTargetMusc(double p, double beta, double teff, double& dangle, double& r) const;
    /// Fortran: init.f:target_init(using_Eloss,using_Coulomb,the_cent,thp_cent,Pp_cent,Mh2,ebeam,Pe_cent)
    /// @param using_Eloss energy loss is simulated.
    /// @param using_Coulomb Coulomb correction is applied.
    /// @param the_cent central electron angle [rad].
    /// @param thp_cent central hadron angle [rad].
    /// @param Pp_cent central hadron momentum [MeV/c].
    /// @param Mh2 hadron mass squared [MeV^2].
    /// @param ebeam beam energy [MeV].
    /// @param Pe_cent central electron momentum [MeV/c].
    void targetInit(bool using_Eloss, bool using_Coulomb, double the_cent, double thp_cent, double Pp_cent,
                    double Mh2, double ebeam, double Pe_cent);

    // ---- state -----------------------------------------------------------------------------
    /// @return the target description (modifiable).
    TargetInfo& targ() { return targ_; }
    /// @return the target description.
    const TargetInfo& targ() const { return targ_; }
    /// @param electron_arm electron arm code.
    /// @param hadron_arm hadron arm code.
    void setArms(int electron_arm, int hadron_arm) {
        electron_arm_ = electron_arm;
        hadron_arm_ = hadron_arm;
    }
    /// @return the electron arm code.
    int electronArm() const { return electron_arm_; }
    /// @return the hadron arm code.
    int hadronArm() const { return hadron_arm_; }
    /// @return the target_init inputs.
    const TargetInitInputs& initInputs() const { return init_; }
    /// @return the target_init inputs (modifiable).
    TargetInitInputs& initInputs() { return init_; }

    /// @param m the material model (ownership taken).
    void setMaterialModel(std::unique_ptr<IMaterialModel> m) { material_ = std::move(m); }
    /// @return the material model (owned by the target).
    IMaterialModel& materialModel() { return *material_; }

    /// Where the Fortran write(6,*) warnings go (default std::cout; nullptr = silent).
    /// Also forwarded to the material model if it is a SimcMaterial.
    /// @param os the stream (not owned), or nullptr.
    void setLog(std::ostream* os);

    /// @return th_can [rad].
    double thCan() const { return thCan_; }   ///< static th_can of the pudding can (F4)

private:
    TargetInfo targ_;
    int electron_arm_{0}, hadron_arm_{0};
    TargetInitInputs init_;
    std::unique_ptr<IMaterialModel> material_;   ///< owned; default SimcMaterial ("simc_material")
    std::ostream* log_;
    double thCan_{0};   ///< target.f: th_can (static local, -fno-automatic)
    /// target.f static locals s_air, s_kevlar, s_mylar, forward_path: reused by the calorimeter arms (D-012 F4)
    double win_air_{0}, win_kevlar_{0}, win_mylar_{0}, win_forward_{0};
};

/// @}

}  // namespace material
}  // namespace simc

#endif  // SIMC_MATERIAL_HALLCCRYOTARGET_H
