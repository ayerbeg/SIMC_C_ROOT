// include/simc/stages/ISpectrometer.h
// Single-arm spectrometer Monte Carlo: forward transport through apertures and magnets,
// detector hut, reconstruction.
// Fortran: hms/mc_hms.f, shms/mc_shms.f, sos/mc_sos.f, hrsr/, hrsl/, calo/ — the argument
//          list of mc_hms (mc_hms.f l.1-4) is mirrored by SpecInput/SpecResult.
#ifndef SIMC_STAGES_ISPECTROMETER_H
#define SIMC_STAGES_ISPECTROMETER_H

#include <map>
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
    double fry{0};         ///< vertical beam position for raster reconstruction (cm)
};

/// Output of the single-arm MC.
struct SpecResult {
    bool ok{false};
    std::string stop_reason;         ///< name of the aperture that stopped the track ("" if ok)
    double x_fp{0}, dx_fp{0}, y_fp{0}, dy_fp{0};   ///< focal plane (cm, slopes)
    double delta_recon{0}, dxdz_recon{0}, dydz_recon{0}, y_recon{0};  ///< reconstructed target
    double pathlen{0};               ///< cm
    double resmult{0};               ///< DC resolution factor (ntup%resfac)
    double m2_final{0};              ///< mass squared after possible decay
};

class ISpectrometer : public IStage {
public:
    virtual SpecResult transport(const SpecInput& in, IRandom& rng) = 0;
    /// Per-aperture stop counters (Fortran hSTOP_*, shmsSTOP_*), including "trials" and "successes".
    virtual const std::map<std::string, long>& stopCounters() const = 0;
};

}  // namespace simc

#endif  // SIMC_STAGES_ISPECTROMETER_H
