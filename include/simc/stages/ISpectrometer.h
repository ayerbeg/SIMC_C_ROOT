/**
 * @file
 * @ingroup stages
 * @brief simc::ISpectrometer: a single-arm spectrometer Monte Carlo (transport through apertures and magnets,
 *        detector hut, reconstruction), with its input and output tracks.
 *
 * The physics content of the argument list of the original arm routines (mc_hms.f l.1-4) is carried by
 * simc::SpecInput / simc::SpecResult; the exact argument mirrors (McHmsArgs, McShmsArgs, ...) are private to the
 * spectrometer kernels.
 * @note Origin: Fortran SIMC hms/mc_hms.f, shms/mc_shms.f, sos/mc_sos.f, hrsr/, hrsl/, calo/.
 */
#ifndef SIMC_STAGES_ISPECTROMETER_H
#define SIMC_STAGES_ISPECTROMETER_H

#include <map>
#include <optional>
#include <string>

#include "simc/core/Errors.h"
#include "simc/stages/IStage.h"

namespace simc {

/** @addtogroup stages
 *  @{ */

/**
 * @brief Input track of a spectrometer arm, in the spectrometer TRANSPORT system at z = 0.
 *
 * Built by the driver's Monte Carlo of the try from the orig event (simc::driver::MonteCarlo).
 * @note Origin: Fortran SIMC simc.f:montecarlo l.1412-1441.
 */
struct SpecInput {
    double p_spec{0};      ///< central momentum (MeV/c)
    double th_spec{0};     ///< central angle (rad)
    double delta{0};       ///< dp/p (%)
    double x{0}, y{0}, z{0};        ///< track position in the TRANSPORT system [cm]
    double dxdz{0}, dydz{0};        ///< track slopes in the TRANSPORT system [dx/dz, dy/dz]
    double m2{0};          ///< particle mass squared (MeV^2)
    bool ms_flag{false};   ///< multiple scattering (mc_smear)
    bool wcs_flag{false};  ///< wire-chamber smearing (mc_smear)
    bool decay_flag{false};  ///< decay in flight of the particle (doing_decay)
    /// Recon argument `fry` of mc_hms/mc_shms: simc.f passes xtar_init (= x at z = 0, simc.f:1441/1462) (cm).
    double fry{0};
    /// simc.f:montecarlo l.1349-1355 fry / frx: the fast-raster position (-main%target%rastery / -rasterx if
    /// correct_raster, else 0; cm). mc_sos (simc.f:1472/1738) and mc_hrsr/mc_hrsl receive raster_fry as their recon
    /// argument `fry` (task P14.1).
    double raster_fry{0};
    double raster_frx{0};  ///< fast-raster x position (-main.target.rasterx if correct_raster, else 0) [cm]
    /// Path length already accumulated when the routine is called (simc.f passes `pathlen` in/out; 0 for a normal arm
    /// call, simc.f:1462/1727). Used by the pi0 two-photon calorimeter logic in fortran_compat (D-062 G3).
    double pathlen_start{0};
};

/**
 * @brief Output of a single-arm Monte Carlo: whether the track survived, its focal-plane track and the reconstructed
 *        target quantities.
 */
struct SpecResult {
    bool ok{false};                  ///< the track reached the detectors and was reconstructed
    std::string stop_reason;         ///< name of the aperture/detector that stopped the track ("" if ok)
    double x_fp{0}, dx_fp{0}, y_fp{0}, dy_fp{0};   ///< focal plane (cm, slopes); used when ok
    double delta_recon{0}, dxdz_recon{0}, dydz_recon{0}, y_recon{0};  ///< reconstructed target (%, slopes, cm); when ok
    double pathlen{0};               ///< cm
    /// Resolution factor of the track reconstruction (ntup%resfac; HMS/SHMS: 1 or 2 from the wire chambers), if the
    /// spectrometer determined one for this track (HMS/SHMS: the detector hut was reached). Empty otherwise: the
    /// driver then leaves ntup%resfac unchanged for the hadron arm and re-uses the previous value for the electron arm
    /// (FORTRAN-DIFF-CANDIDATE simc.f:1334, SAVEd tmpfact).
    std::optional<double> resolution_factor;
    double m2_final{0};              ///< particle mass squared after a possible decay (MeV^2)
    double decay_distance{0};        ///< decay distance (cm; Fortran COMMON decdist, 0 if no decay)
    double Mh2_final{0};             ///< hadron mass squared after a possible decay (Fortran COMMON Mh2_final; = m2 if none)
};

/**
 * @brief A spectrometer arm: transports a track from the target to the focal plane, applies the apertures, smears the
 *        detected track and reconstructs the target quantities.
 *
 * Implementations: simc::hms::HmsSpectrometer ("hms"), simc::shms::ShmsSpectrometer ("shms"),
 * simc::sos::SosSpectrometer ("sos"), simc::hrs::HrsrSpectrometer / simc::hrs::HrslSpectrometer ("hrsr", "hrsl"),
 * simc::calo::CaloSpectrometer ("calo"). Each arm of a run has its own instance, which counts the stopped tracks per
 * aperture; not thread-safe.
 */
class ISpectrometer : public IStage {
public:
    /**
     * @brief Transports one track through the spectrometer and reconstructs it.
     * @param in The track at the target, in the TRANSPORT system.
     * @param rng The random source (multiple scattering, decay, wire-chamber smearing).
     * @return The result; result.ok is false if an aperture stopped the track (result.stop_reason names it).
     */
    virtual SpecResult transport(const SpecInput& in, IRandom& rng) = 0;
    /**
     * @brief The per-aperture stop counters.
     * @return Counts keyed by aperture name, including "trials" and "successes" (valid while the stage lives).
     * @note Origin: Fortran SIMC hSTOP_*, shmsSTOP_*, ... counters.
     */
    virtual const std::map<std::string, long>& stopCounters() const = 0;
    /**
     * @brief Reconstructs again the focal-plane track of the LAST transport(), with another raster position.
     *
     * Used by the polarised-target field tracking (task P15.2b). Default: throws simc::NotPorted.
     * @param fry Vertical raster position given to the reconstruction [cm].
     * @param delta Out: reconstructed delta [%].
     * @param dydz Out: reconstructed dy/dz (delta_t).
     * @param dxdz Out: reconstructed dx/dz (delta_phi).
     * @param y Out: reconstructed y target [cm].
     * @throws simc::NotPorted if the spectrometer does not support it.
     * @note Origin: Fortran SIMC mc_<arm>_recon(delta_p, delta_t, delta_phi, y_tgt, fry) on COMMON /track/, as
     *       trg_track.f:track_to_tgt calls it.
     */
    virtual void reconstructLastTrack(double fry, double& delta, double& dydz, double& dxdz, double& y) {
        (void)fry, (void)delta, (void)dydz, (void)dxdz, (void)y;
        throw NotPorted("spectrometer '" + name() + "': no reconstruction for the target-field tracking (track_to_tgt)");
    }
};

/** @} */

}  // namespace simc

#endif  // SIMC_STAGES_ISPECTROMETER_H
