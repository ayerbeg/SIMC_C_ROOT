// include/simc/spectrometers/sos/SosSpectrometer.h
//
// ISpectrometer "sos": the SOS (Hall C Short Orbit Spectrometer, arm code 2) plug-in (task P14.1). It adapts the
// Fortran-faithful kernel sos::McSos (kernel/McSos.h: mc_sos.f, mc_sos_hut.f, mc_sos_recon.f on the shared COSY engine)
// to the stage interface, exactly as simc.f:montecarlo calls mc_sos (l.1468-1472 hadron arm: m2 = Mh2, decay =
// doing_decay, resmult -> ntup%resfac; l.1734-1738 electron arm: m2 = me2, decay = .false., resmult -> tmpfact;
// ms = wcs = mc_smear, pathlen = 0, and the recon argument fry = the RASTER fry of simc.f:1349-1355, SpecInput::raster_fry,
// not xtar_init as for HMS/SHMS).
//
// JSON config (all optional):
//   "matrix_dir"    : directory holding sos/forward_cosy.dat, sos/recon_cosy.dat (default: the repository copies)
//   "ctau"          : cm (default params.ctau)
//   "ms_flag", "wcs_flag", "decay": bool — if present they OVERRIDE the per-call SpecInput flags
//   "use_sieve"     : bool (default false; mc_sos.f:50 mock sieve slit, a Fortran compile-time option)
//   "verbose"       : bool (default true; false = no Fortran write(6,*) messages)
// Stop counters (stopCounters()): trials, slit_hor, slit_vert, slit_oct, quad_in, quad_mid, quad_out, bm01_in, bm01_out,
// bm02_in, bm02_out, exit, hut, dc1, dc2, scin, cal, successes (COMMON /spectrometer2/ sSTOP_*).
#ifndef SIMC_SPECTROMETERS_SOS_SOSSPECTROMETER_H
#define SIMC_SPECTROMETERS_SOS_SOSSPECTROMETER_H

#include <map>
#include <memory>
#include <optional>
#include <string>

#include "simc/stages/ISpectrometer.h"

namespace simc {
namespace sos {

class McSos;   // kernel (kernel/McSos.h)

class SosSpectrometer : public ISpectrometer {
public:
    SosSpectrometer();
    ~SosSpectrometer() override;
    SosSpectrometer(const SosSpectrometer&) = delete;
    SosSpectrometer& operator=(const SosSpectrometer&) = delete;

    std::string name() const override { return "sos"; }
    void configure(const Json& config, const SimcParameters& params) override;
    /// Loads the COSY maps.
    void init(const StageContext& ctx) override;
    SpecResult transport(const SpecInput& in, IRandom& rng) override;
    const std::map<std::string, long>& stopCounters() const override;

    /// The Fortran-faithful kernel (L1 tests, diagnostics).
    McSos& kernel() { return *k_; }
    const McSos& kernel() const { return *k_; }

private:
    std::unique_ptr<McSos> k_;
    double ctau_{0};
    std::optional<bool> msOverride_, wcsOverride_, decayOverride_;
    mutable std::map<std::string, long> cntMap_;
};

}  // namespace sos
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_SOS_SOSSPECTROMETER_H
