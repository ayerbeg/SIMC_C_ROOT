/**
 * @file
 * @ingroup spectrometers
 * @brief simc::shms::ShmsSpectrometer: ISpectrometer "shms", the SHMS plug-in on the Fortran-faithful kernel McShms.
 *
 * ISpectrometer "shms": the SHMS (Hall C Super High Momentum Spectrometer) plug-in. It adapts the Fortran-faithful
 * kernel shms::McShms (kernel/McShms.h: mc_shms.f, mc_shms_coll.f, mc_shms_hut.f + hut.inc, mc_shms_recon.f on the
 * shared COSY engine) to the stage interface, exactly as simc.f:montecarlo calls mc_shms (l.1482-1489 hadron arm,
 * l.1747-1755 electron arm).
 *
 * JSON config (all optional):
 *   "matrix_dir"     : holds shms/shms_forward.dat, shms/shms_recon.dat (default: the repository copies)
 *   (the former "spectr" key is rejected: the SHMS is always spectrometer 5, D-032)
 *   "using_SHMScoll" : bool (default params.using_SHMScoll)
 *   "ctau"           : cm (default params.ctau)
 *   "ms_flag", "wcs_flag", "decay": bool, override the per-call SpecInput flags if present
 *   "use_sieve", "use_2sieve" : mock sieve options of mc_shms.f:58-66 (default false)
 *   "use_coll", "skip_hb"     : mc_shms.f compile-time options (defaults true / false, as in SIMC)
 *   "verbose"        : bool (default true)
 * Stop counters (stopCounters()): every shmsSTOP_* counter by name, plus "q2_mid_dup" (= q2_mid, simc.f:519-520 label
 * bug), without the COMMON variable stop_id.
 */
#ifndef SIMC_SPECTROMETERS_SHMS_SHMSSPECTROMETER_H
#define SIMC_SPECTROMETERS_SHMS_SHMSSPECTROMETER_H

#include <map>
#include <memory>
#include <optional>
#include <string>

#include "simc/stages/ISpectrometer.h"

namespace simc {
namespace shms {

class McShms;   // kernel (kernel/McShms.h)

/// @addtogroup spectrometers
/// @{

/// @brief ISpectrometer "shms": the SHMS single-arm Monte Carlo (apertures, COSY transport, hut, reconstruction).
///
/// Owns its kernel (the Fortran COMMON counters and SAVEd state). Not copyable, not thread-safe (one per worker).
class ShmsSpectrometer : public ISpectrometer {
public:
    /// Creates the plug-in with a default kernel.
    ShmsSpectrometer();
    /// Destroys the kernel.
    ~ShmsSpectrometer() override;
    /// Not copyable (owns Fortran state).
    ShmsSpectrometer(const ShmsSpectrometer&) = delete;
    /// Not copyable.
    /// @return (deleted)
    ShmsSpectrometer& operator=(const ShmsSpectrometer&) = delete;

    std::string name() const override { return "shms"; }
    void configure(const Json& config, const SimcParameters& params) override;
    /// Loads the COSY maps (spectrometer 5).
    /// @param ctx the stage context.
    void init(const StageContext& ctx) override;
    SpecResult transport(const SpecInput& in, IRandom& rng) override;
    const std::map<std::string, long>& stopCounters() const override;
    /// mc_shms_recon on the last hut track (trg_track.f:track_to_tgt, P15.2b).
    /// @param fry vertical beam position [cm].
    /// @param[out] delta reconstructed delta [%].
    /// @param[out] dydz reconstructed yptar.
    /// @param[out] dxdz reconstructed xptar.
    /// @param[out] y reconstructed ytar [cm].
    void reconstructLastTrack(double fry, double& delta, double& dydz, double& dxdz, double& y) override;

    /// The Fortran-faithful kernel (L1 tests, diagnostics).
    /// @return the kernel (owned by the plug-in).
    McShms& kernel() { return *k_; }
    /// @return the kernel (owned by the plug-in).
    const McShms& kernel() const { return *k_; }

private:
    std::unique_ptr<McShms> k_;
    bool usingColl_{false};
    double ctau_{0};
    std::optional<bool> msOverride_, wcsOverride_, decayOverride_;
    mutable std::map<std::string, long> cntMap_;
};

/// @}

}  // namespace shms
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_SHMS_SHMSSPECTROMETER_H
