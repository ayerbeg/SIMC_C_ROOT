// include/simc/spectrometers/hrs/HrsSpectrometer.h
//
// ISpectrometer "hrsr" (Hall A HRS-R, arm code 3) and "hrsl" (HRS-L, arm code 4), task P14.2: one class with a side
// parameter around the Fortran-faithful kernel hrs::McHrs (kernel/McHrs.h: mc_hrsr.f / mc_hrsl.f, *_hut.f, *_recon.f on
// the shared COSY engine), called exactly as simc.f:montecarlo calls mc_hrsr / mc_hrsl (l.1473-1482 hadron arm,
// l.1739-1748 electron arm: m2 = Mh2 / me2, decay = doing_decay / .false., ms = wcs = mc_smear, pathlen = 0, recon
// fry = the RASTER fry, SpecInput::raster_fry).
//
// JSON config (all optional): "matrix_dir" (holds hrs/hrs_forward_cosy.dat, hrs/hrs_recon_cosy.dat), "ctau",
// "ms_flag"/"wcs_flag"/"decay" (override the SpecInput flags), "use_sieve" (mc_hrs*.f:47), "verbose".
// Stop counters: trials, slit_hor, slit_vert, Q1_in ... D1_out, hut, dc1, dc2, s1, s2, cal, successes (rSTOP_* / lSTOP_*).
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

class HrsSpectrometer : public ISpectrometer {
public:
    explicit HrsSpectrometer(Side side);
    ~HrsSpectrometer() override;
    HrsSpectrometer(const HrsSpectrometer&) = delete;
    HrsSpectrometer& operator=(const HrsSpectrometer&) = delete;

    std::string name() const override;
    void configure(const Json& config, const SimcParameters& params) override;
    void init(const StageContext& ctx) override;
    SpecResult transport(const SpecInput& in, IRandom& rng) override;
    const std::map<std::string, long>& stopCounters() const override;

    McHrs& kernel() { return *k_; }
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
    HrsrSpectrometer();
};
/// "hrsl": HRS-L (arm 4).
class HrslSpectrometer : public HrsSpectrometer {
public:
    HrslSpectrometer();
};

}  // namespace hrs
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_HRS_HRSSPECTROMETER_H
