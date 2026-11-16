/**
 * @file
 * @ingroup driver
 * @brief simc::driver::MonteCarlo: per-event transport of both arms through the spectrometers (simc.f:montecarlo).
 *
 * Per-event transport of both arms (task P7.1).
 * Fortran: simc.f:montecarlo(orig,main,recon,success) (l.1310-1853).
 *   * beam multiple scattering (target_musc with teff(1)) when mc_smear;
 *   * hadron arm first (l.1374-1662): SP quantities (Eloss, target_musc with teff(3)), the conversion to
 *     TRANSPORT coordinates (spectrometer phi, offsets, drift back to z = 0, xtar_init passed as fry to HMS/SHMS, the
 *     raster fry to SOS/HRS), the single-arm MC (ISpectrometer "hms"/"sos"/"hrsr"/"hrsl"/"shms" by hadron_arm), RECON/FP, then the
 *     recon arm quantities
 *     (physics_angles, optional most-probable energy-loss correction);
 *   * electron arm (l.1664-1849): the same with Eloss(2), Coulomb, teff(2), m2 = me2, decay off;
 *   * success = .true. only if both arms pass.
 *
 * How the main loop calls it (simc.f:197-207): after a successful generate,
 *     bool ok = mc.run(orig, main, recon, rng);
 * with the SAME `main` and `recon` objects for every try: montecarlo writes main.SP/RECON/FP and the recon arm
 * fields in place and returns early on a failed arm, so (as in SIMC) those structures keep partly stale values of
 * earlier tries; main.ntup.resfac is reset here (simc.f:1352) and must also be reset by the loop (simc.f:178).
 * main.ntup.decdist / Mh2_final are set from the hadron arm (the loop resets them per event, simc.f:174-177).
 *
 * Calorimeter arms 7/8 (P14.3): "calo"; with doing_pizero the two photons go through driver::PizeroCalorimeter.
 * Polarized-target field (using_tgt_field, P15.2b): optics::TargetField (trg_track.f), owned here: track_from_tgt on
 * the vertex before the drift to z = 0 (simc.f:1425-1432, 1693-1701) and track_to_tgt after a successful arm
 * (simc.f:1573-1590, 1777-1794, ISpectrometer::reconstructLastTrack). A calorimeter arm with the field STOPs (T17).
 *
 * Fortran behaviour reproduced (docs/translation/markdown/FORTRAN_DIFFERENCES.md):
 *   * the electron-arm resmult goes through the SAVEd local tmpfact (simc.f:1334): if the electron-arm MC stops
 *     before writing it (HMS before the hut), ntup%resfac gets the value of an EARLIER call; it is added before
 *     the ok check, so it also happens for failed electron arms;
 *   * the hadron-arm resmult is the argument ntup%resfac itself: left at 0 if the HMS/SOS hut is not reached;
 *   * RECON%p%z / RECON%e%z are not written when the arm Monte Carlo is off (using_*_arm_montecarlo false): the
 *     recon z keeps the previous value;
 *   * recon%e%E is not corrected for Coulomb (l.1840-1842).
 */
#ifndef SIMC_DRIVER_MONTECARLO_H
#define SIMC_DRIVER_MONTECARLO_H

#include <iosfwd>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "simc/core/Event.h"
#include "simc/core/RunSetup.h"
#include "simc/optics/TargetField.h"
#include "simc/stages/ISpectrometer.h"
#include "simc/stages/ITarget.h"
#include "simc/core/Errors.h"

namespace simc {

class IRandom;

namespace driver {


/// @addtogroup driver
/// @{

/// Diagnostics of one montecarlo call.
struct MonteCarloInfo {
    bool success{false};  ///< both arms passed.
    int fail_arm{0};                  ///< 0 none, 1 hadron arm, 2 electron arm
    bool p_called{false} /**< @copydoc e_called */, e_called{false};   ///< the arm MC was called
    std::string p_stop /**< @copydoc e_stop */, e_stop;       ///< ISpectrometer stop reason ("" = passed / not called)
};

/// @brief Transport of one event through the hadron and electron arms (beam multiple scattering, energy loss,
///        TRANSPORT conversion, single-arm Monte Carlo, reconstruction); see the file description.
///
/// Owns the arm spectrometer plug-ins and the optional target-field tracking. Holds per-run Fortran state
/// (tmpfact, spectrometer counters): not thread-safe; use one object per worker.
class MonteCarlo {
public:
    /// P7.2 (additive): an arm plug-in = ISpectrometer registry name + its JSON config.
    struct ArmChoice {
        std::string type;               ///< "" = from the Fortran arm code (DefaultArm)
        Json config = Json::object();   ///< merged over the default config of the arm code
    };
    /// The Fortran choice for an arm code (simc.f:montecarlo dispatch, core/SpectrometerId.h): HMS -> "hms",
    /// SOS -> "sos", HRS-R/L -> "hrsr"/"hrsl", SHMS -> "shms", calorimeter 7/8 -> "calo". std::invalid_argument for a
    /// code that is not a spectrometer.
    /// @param code arm code (1-5, 7, 8; core/SpectrometerId.h).
    /// @return the registry name of the arm with an empty config.
    /// @throws std::invalid_argument if @p code is not a spectrometer.
    static ArmChoice DefaultArm(int code);

    /// Creates (through Registry<ISpectrometer>) and initialises the arm spectrometers of setup.params
    /// (electron_arm, hadron_arm). `target` provides target_musc and trip_thru_target (correct_Eloss).
    /// With using_tgt_field: computes the arm angles (simc.f:120-151) and loads the field map (trgInit).
    /// @param setup the run set-up; must outlive this object.
    /// @param target the target model; must outlive this object.
    MonteCarlo(const RunSetup& setup, ITarget& target);
    /// P7.2 (additive): the same with the arm plug-ins chosen explicitly (e.g. from the JSON "pipeline").
    /// Arms with the same type and the same (merged) config share one object (one Fortran COMMON block of counters);
    /// otherwise each arm has its own object (P7.4).
    /// @param setup the run set-up; must outlive this object.
    /// @param target the target model; must outlive this object.
    /// @param hadron plug-in of the hadron arm.
    /// @param electron plug-in of the electron arm.
    MonteCarlo(const RunSetup& setup, ITarget& target, const ArmChoice& hadron, const ArmChoice& electron);

    /// Fortran: simc.f:montecarlo(orig,main,recon,success). Returns success.
    /// @param orig the event after radiation (input).
    /// @param main per-event main record; SP/RECON/FP and ntup fields written in place (may keep stale values).
    /// @param recon reconstructed event, written in place (may keep stale values on a failed arm).
    /// @param rng the random generator of the run.
    /// @return true if both arms passed.
    bool run(const Event& orig, MainEvent& main, Event& recon, IRandom& rng);

    /// @return diagnostics of the last run() call.
    const MonteCarloInfo& info() const { return info_; }
    /// The spectrometer used for the hadron / electron arm (the same object if both arms are the same type).
    /// @return the hadron-arm spectrometer (owned by this object).
    ISpectrometer& hadronArm() { return *pArm_; }
    /// @return the electron-arm spectrometer (owned by this object; may be the hadron-arm object).
    ISpectrometer& electronArm() { return *eArm_; }

    /// Where the Fortran write(6,*) messages go (default std::cout; nullptr = silent).
    /// @param os output stream (not owned; must outlive this object), or nullptr for silence.
    void setLog(std::ostream* os) { log_ = os; }
    /// using_tgt_field: the field tracking (nullptr otherwise) and the simc.f:120-151 log lines.
    /// @return the field tracking (owned by this object), or nullptr without using_tgt_field.
    const optics::TargetField* targetField() const { return field_.get(); }
    /// @return the log lines of the field-angle computation (simc.f:120-151).
    const std::vector<std::string>& targetFieldLog() const { return fieldLog_; }

private:
    ISpectrometer* makeArm(int code, const ArmChoice& choice);
    void initTargetField();
    double fieldStheta(int arm, double stheta) const;

    const RunSetup& setup_;
    ITarget& target_;
    std::map<std::string, std::unique_ptr<ISpectrometer>> specs_;   ///< owned, one per (type, config): COMMON counters
    ISpectrometer* pArm_{nullptr};   ///< the hadron arm (points into specs_)
    ISpectrometer* eArm_{nullptr};   ///< the electron arm (may be the same object as pArm_)
    MonteCarloInfo info_;
    std::unique_ptr<optics::TargetField> field_;   ///< using_tgt_field (trg_track.f COMMON blocks)
    std::vector<std::string> fieldLog_;
    double tmpfact_{0};   ///< simc.f:1334 local tmpfact (static under -fno-automatic)
    std::ostream* log_;
};

/// @}

}  // namespace driver
}  // namespace simc

#endif  // SIMC_DRIVER_MONTECARLO_H
