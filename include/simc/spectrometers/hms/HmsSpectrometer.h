// include/simc/spectrometers/hms/HmsSpectrometer.h
//
// ISpectrometer "hms": the HMS (Hall C High Momentum Spectrometer) plug-in. It adapts the Fortran-faithful kernel
// hms::McHms (kernel/McHms.h: mc_hms.f, mc_hms_coll.f, mc_hms_hut.f, mc_hms_recon.f on the shared COSY engine) to the
// stage interface: SpecInput -> mc_hms argument list -> SpecResult, exactly as simc.f:montecarlo calls mc_hms
// (l.1459-1500: m2 = Mh2, pathlen = 0, ms = wcs = mc_smear, decay = doing_decay, resmult -> ntup%resfac,
// fry = xtar_init).
//
// JSON config (all optional):
//   "matrix_dir"    : directory holding hms/forward_cosy.dat, hms/recon_cosy.dat (default: the repository copies)
//   "using_HMScoll" : bool (default params.using_HMScoll) — step pions/muons through the collimator
//   "ctau"          : cm (default params.ctau)
//   "ms_flag", "wcs_flag", "decay": bool — if present they OVERRIDE the per-call SpecInput flags
//   "use_sieve"     : bool (default false; mc_hms.f:55 mock sieve slit, a Fortran compile-time option)
//   "verbose"       : bool (default true; false = no Fortran write(6,*) messages)
// Stop counters (stopCounters()): trials, slit_hor, slit_vert, slit_oct, Q1_in ... D1_out, hut, dc1, dc2, scin, cal,
// successes, coll (COMMON /spectrometer1/ hSTOP_*).
#ifndef SIMC_SPECTROMETERS_HMS_HMSSPECTROMETER_H
#define SIMC_SPECTROMETERS_HMS_HMSSPECTROMETER_H

#include <map>
#include <memory>
#include <optional>
#include <string>

#include "simc/stages/ISpectrometer.h"

namespace simc {
namespace hms {

class McHms;   // kernel (kernel/McHms.h)

class HmsSpectrometer : public ISpectrometer {
public:
    HmsSpectrometer();
    ~HmsSpectrometer() override;
    HmsSpectrometer(const HmsSpectrometer&) = delete;
    HmsSpectrometer& operator=(const HmsSpectrometer&) = delete;

    std::string name() const override { return "hms"; }
    void configure(const Json& config, const SimcParameters& params) override;
    /// Loads the COSY maps.
    void init(const StageContext& ctx) override;
    SpecResult transport(const SpecInput& in, IRandom& rng) override;
    const std::map<std::string, long>& stopCounters() const override;
    /// mc_hms_recon on the last hut track (trg_track.f:track_to_tgt, P15.2b).
    void reconstructLastTrack(double fry, double& delta, double& dydz, double& dxdz, double& y) override;

    /// The Fortran-faithful kernel (L1 tests, diagnostics).
    McHms& kernel() { return *k_; }
    const McHms& kernel() const { return *k_; }

private:
    std::unique_ptr<McHms> k_;
    bool usingColl_{false};
    double ctau_{0};
    std::optional<bool> msOverride_, wcsOverride_, decayOverride_;
    mutable std::map<std::string, long> cntMap_;
};

}  // namespace hms
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_HMS_HMSSPECTROMETER_H
