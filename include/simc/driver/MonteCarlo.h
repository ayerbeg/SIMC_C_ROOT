// include/simc/driver/MonteCarlo.h
//
// Per-event transport of both arms (task P7.1).
// Fortran: simc.f:montecarlo(orig,main,recon,success) (l.1310-1853).
//   * beam multiple scattering (target_musc with teff(1)) when mc_smear;
//   * hadron arm first (l.1374-1662): SP quantities (Eloss, target_musc with teff(3)), the conversion to
//     TRANSPORT coordinates (spectrometer phi, offsets, drift back to z = 0, xtar_init passed as fry to HMS/SHMS, the
//     raster fry to SOS/HRS), the single-arm MC (ISpectrometer "hms"/"sos"/"hrsr"/"hrsl"/"shms" by hadron_arm), RECON/FP, then the
//     recon arm quantities
//     (physics_angles, optional most-probable energy-loss correction);
//   * electron arm (l.1664-1849): the same with Eloss(2), Coulomb, teff(2), m2 = me2, decay off;
//   * success = .true. only if both arms pass.
//
// How the main loop calls it (simc.f:197-207): after a successful generate,
//     bool ok = mc.run(orig, main, recon, rng);
// with the SAME `main` and `recon` objects for every try: montecarlo writes main.SP/RECON/FP and the recon arm
// fields in place and returns early on a failed arm, so (as in SIMC) those structures keep partly stale values of
// earlier tries; main.ntup.resfac is reset here (simc.f:1352) and must also be reset by the loop (simc.f:178).
// main.ntup.decdist / Mh2_final are set from the hadron arm (the loop resets them per event, simc.f:174-177).
//
// Not ported (throw simc::NotPorted): the calorimeter arms 7/8
// (mc_calo), the polarized-target field (using_tgt_field: track_from_tgt / track_to_tgt). Arm code 6 (SHMS LSA)
// reaches mc_shms, which STOPs (simc::FortranStop).
//
// Fortran behaviour reproduced (docs/FORTRAN_DIFFERENCES.md):
//   * the electron-arm resmult goes through the SAVEd local tmpfact (simc.f:1334): if the electron-arm MC stops
//     before writing it (HMS before the hut), ntup%resfac gets the value of an EARLIER call; it is added before
//     the ok check, so it also happens for failed electron arms;
//   * the hadron-arm resmult is the argument ntup%resfac itself: left at 0 if the HMS/SOS hut is not reached;
//   * RECON%p%z / RECON%e%z are not written when the arm Monte Carlo is off (using_*_arm_montecarlo false): the
//     recon z keeps the previous value;
//   * recon%e%E is not corrected for Coulomb (l.1840-1842).
#ifndef SIMC_DRIVER_MONTECARLO_H
#define SIMC_DRIVER_MONTECARLO_H

#include <iosfwd>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include "simc/core/Event.h"
#include "simc/core/RunSetup.h"
#include "simc/stages/ISpectrometer.h"
#include "simc/stages/ITarget.h"
#include "simc/core/Errors.h"

namespace simc {

class IRandom;

namespace driver {


/// Diagnostics of one montecarlo call.
struct MonteCarloInfo {
    bool success{false};
    int fail_arm{0};                  ///< 0 none, 1 hadron arm, 2 electron arm
    bool p_called{false}, e_called{false};   ///< the arm MC was called
    std::string p_stop, e_stop;       ///< ISpectrometer stop reason ("" = passed / not called)
};

class MonteCarlo {
public:
    /// P7.2 (additive): an arm plug-in = ISpectrometer registry name + its JSON config.
    struct ArmChoice {
        std::string type;               ///< "" = from the Fortran arm code (DefaultArm)
        Json config = Json::object();   ///< merged over the default config of the arm code
    };
    /// The Fortran choice for an arm code (simc.f:montecarlo dispatch, core/SpectrometerId.h): HMS -> "hms",
    /// SOS -> "sos", HRS-R/L -> "hrsr"/"hrsl", SHMS -> "shms". Throws NotPorted for the calorimeters (mc_calo), std::invalid_argument for a
    /// code that is not a spectrometer.
    static ArmChoice DefaultArm(int code);

    /// Creates (through Registry<ISpectrometer>) and initialises the arm spectrometers of setup.params
    /// (electron_arm, hadron_arm). `target` provides target_musc and trip_thru_target (correct_Eloss).
    /// Throws NotPorted for the calorimeter arms, or for using_tgt_field.
    MonteCarlo(const RunSetup& setup, ITarget& target);
    /// P7.2 (additive): the same with the arm plug-ins chosen explicitly (e.g. from the JSON "pipeline").
    /// Arms with the same type and the same (merged) config share one object (one Fortran COMMON block of counters);
    /// otherwise each arm has its own object (P7.4).
    MonteCarlo(const RunSetup& setup, ITarget& target, const ArmChoice& hadron, const ArmChoice& electron);

    /// Fortran: simc.f:montecarlo(orig,main,recon,success). Returns success.
    bool run(const Event& orig, MainEvent& main, Event& recon, IRandom& rng);

    const MonteCarloInfo& info() const { return info_; }
    /// The spectrometer used for the hadron / electron arm (the same object if both arms are the same type).
    ISpectrometer& hadronArm() { return *pArm_; }
    ISpectrometer& electronArm() { return *eArm_; }

    /// Where the Fortran write(6,*) messages go (default std::cout; nullptr = silent).
    void setLog(std::ostream* os) { log_ = os; }

private:
    ISpectrometer* makeArm(int code, const ArmChoice& choice);

    const RunSetup& setup_;
    ITarget& target_;
    std::map<std::string, std::unique_ptr<ISpectrometer>> specs_;   ///< owned, one per (type, config): COMMON counters
    ISpectrometer* pArm_{nullptr};   ///< the hadron arm (points into specs_)
    ISpectrometer* eArm_{nullptr};   ///< the electron arm (may be the same object as pArm_)
    MonteCarloInfo info_;
    double tmpfact_{0};   ///< simc.f:1334 local tmpfact (static under -fno-automatic)
    std::ostream* log_;
};

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_MONTECARLO_H
