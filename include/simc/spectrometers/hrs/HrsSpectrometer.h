/**
 * @file
 * @ingroup spectrometers
 * @brief simc::hrs::HrsSpectrometer: ISpectrometer "hrsr" and "hrsl", the Hall A HRS plug-ins on the kernel McHrs.
 *
 * ISpectrometer "hrsr" (Hall A HRS-R, arm code 3) and "hrsl" (HRS-L, arm code 4), task P14.2: one class with a side
 * parameter around the Fortran-faithful kernel hrs::McHrs (kernel/McHrs.h: mc_hrsr.f / mc_hrsl.f, *_hut.f, *_recon.f on
 * the shared COSY engine), called exactly as simc.f:montecarlo calls mc_hrsr / mc_hrsl (l.1473-1482 hadron arm,
 * l.1739-1748 electron arm: m2 = Mh2 / me2, decay = doing_decay / .false., ms = wcs = mc_smear, pathlen = 0, recon
 * fry = the RASTER fry, SpecInput::raster_fry).
 *
 * JSON config (all optional): "matrix_dir" (holds hrs/hrs_forward_cosy.dat, hrs/hrs_recon_cosy.dat), "ctau",
 * "ms_flag"/"wcs_flag"/"decay" (override the SpecInput flags), "use_sieve" (mc_hrs*.f:47), "verbose".
 * Stop counters: trials, slit_hor, slit_vert, Q1_in ... D1_out, hut, dc1, dc2, s1, s2, cal, successes (rSTOP_* / lSTOP_*).
 */
#ifndef SIMC_SPECTROMETERS_HRS_HRSSPECTROMETER_H
#define SIMC_SPECTROMETERS_HRS_HRSSPECTROMETER_H

#include <map>
#include <memory>
#include <optional>
#include <string>

#include "simc/stages/ISpectrometer.h"

namespace simc {
namespace hrs {

class McHrs;
enum class Side;

/// @addtogroup spectrometers
/// @{

/// @brief Common base of the HRS-R and HRS-L plug-ins: the HRS single-arm Monte Carlo of one side.
///
/// Owns its kernel. Not copyable, not thread-safe (one per worker).
class HrsSpectrometer : public ISpectrometer {
public:
    /// @param side the HRS side.
    explicit HrsSpectrometer(Side side);
    /// Destroys the kernel.
    ~HrsSpectrometer() override;
    /// Not copyable (owns Fortran state).
    HrsSpectrometer(const HrsSpectrometer&) = delete;
    /// Not copyable.
    /// @return (deleted)
    HrsSpectrometer& operator=(const HrsSpectrometer&) = delete;

    std::string name() const override;
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    SpecResult transport(const SpecInput& in, IRandom& rng) override;
    const std::map<std::string, long>& stopCounters() const override;
    /// mc_hrs_recon on the last hut track (trg_track.f:track_to_tgt, P15.2b).
    /// @param fry vertical beam position [cm].
    /// @param[out] delta reconstructed delta [%].
    /// @param[out] dydz reconstructed yptar.
    /// @param[out] dxdz reconstructed xptar.
    /// @param[out] y reconstructed ytar [cm].
    void reconstructLastTrack(double fry, double& delta, double& dydz, double& dxdz, double& y) override;

    /// @return the kernel (owned by the plug-in).
    McHrs& kernel() { return *k_; }
    /// @return the kernel (owned by the plug-in).
    const McHrs& kernel() const { return *k_; }

private:
    std::unique_ptr<McHrs> k_;
    double ctau_{0};
    std::optional<bool> msOverride_, wcsOverride_, decayOverride_;
    mutable std::map<std::string, long> cntMap_;
};

/// "hrsr": HRS-R (arm 3).
class HrsrSpectrometer : public HrsSpectrometer {
public:
    /// Creates the HRS-R plug-in.
    HrsrSpectrometer();
};
/// "hrsl": HRS-L (arm 4).
class HrslSpectrometer : public HrsSpectrometer {
public:
    /// Creates the HRS-L plug-in.
    HrslSpectrometer();
};

/// @}

}  // namespace hrs
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_HRS_HRSSPECTROMETER_H
