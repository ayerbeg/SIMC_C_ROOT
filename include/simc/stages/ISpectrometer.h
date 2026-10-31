// include/simc/stages/ISpectrometer.h
// Single-arm spectrometer Monte Carlo: forward transport through apertures and magnets,
// detector hut, reconstruction.
// Fortran: hms/mc_hms.f, shms/mc_shms.f, sos/mc_sos.f, hrsr/, hrsl/, calo/ — the physics content of the mc_hms
//          argument list (mc_hms.f l.1-4) is carried by SpecInput/SpecResult; the Fortran argument mirrors themselves
//          (McHmsArgs, McShmsArgs) are private to the spectrometer kernels.
#ifndef SIMC_STAGES_ISPECTROMETER_H
#define SIMC_STAGES_ISPECTROMETER_H

#include <map>
#include <optional>
#include <string>

#include "simc/stages/IStage.h"

namespace simc {

/// Input track in the spectrometer TRANSPORT system at z=0 (simc.f:montecarlo l.1412-1441).
struct SpecInput {
    double p_spec{0};      ///< central momentum (MeV/c)
    double th_spec{0};     ///< central angle (rad)
    double delta{0};       ///< dp/p (%)
    double x{0}, y{0}, z{0};        ///< cm
    double dxdz{0}, dydz{0};        ///< slopes
    double m2{0};          ///< particle mass squared (MeV^2)
    bool ms_flag{false};   ///< multiple scattering (mc_smear)
    bool wcs_flag{false};  ///< wire-chamber smearing (mc_smear)
    bool decay_flag{false};
    /// Recon argument `fry` of mc_hms/mc_shms: simc.f passes xtar_init (= x at z = 0, simc.f:1441/1462) (cm).
    double fry{0};
    /// simc.f:montecarlo l.1349-1355 fry / frx: the fast-raster position (-main%target%rastery / -rasterx if
    /// correct_raster, else 0; cm). mc_sos (simc.f:1472/1738) and mc_hrsr/mc_hrsl receive raster_fry as their recon
    /// argument `fry` (task P14.1).
    double raster_fry{0};
    double raster_frx{0};
};

/// Output of the single-arm MC.
struct SpecResult {
    bool ok{false};
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

class ISpectrometer : public IStage {
public:
    virtual SpecResult transport(const SpecInput& in, IRandom& rng) = 0;
    /// Per-aperture stop counters (Fortran hSTOP_*, shmsSTOP_*), including "trials" and "successes".
    virtual const std::map<std::string, long>& stopCounters() const = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_ISPECTROMETER_H
