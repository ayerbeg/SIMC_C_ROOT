// src/spectrometers/sos/SosSpectrometer.cpp — see include/simc/spectrometers/sos/SosSpectrometer.h.
// The ISpectrometer plug-in "sos" around the kernel sos::McSos, and its registration.
#include "simc/spectrometers/sos/SosSpectrometer.h"

#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/spectrometers/sos/kernel/McSos.h"

namespace simc {
namespace sos {

SosSpectrometer::SosSpectrometer() : k_(std::make_unique<McSos>()) {}
SosSpectrometer::~SosSpectrometer() = default;

// Stage configuration: the run parameters first, then the optional JSON keys of the pipeline (matrix_dir, use_sieve,
// and ms_flag/wcs_flag/decay, which override the per-event flags of the event chain; verbose = false silences it).
void SosSpectrometer::configure(const Json& config, const SimcParameters& params) {
    ctau_ = params.ctau;
    if (!config.is_object()) return;
    if (config.contains("matrix_dir")) k_->setMatrixDir(config.at("matrix_dir").get<std::string>());
    if (config.contains("ctau")) ctau_ = config.at("ctau").get<double>();
    if (config.contains("use_sieve")) k_->setUseSieve(config.at("use_sieve").get<bool>());
    if (config.contains("ms_flag")) msOverride_ = config.at("ms_flag").get<bool>();
    if (config.contains("wcs_flag")) wcsOverride_ = config.at("wcs_flag").get<bool>();
    if (config.contains("decay")) decayOverride_ = config.at("decay").get<bool>();
    if (config.contains("verbose") && !config.at("verbose").get<bool>()) k_->setLog(nullptr);
}

void SosSpectrometer::init(const StageContext& ctx) {
    k_->setFortranCompat(ctx.setup.fortran_compat);   // D-031 run option
    k_->loadMaps();
}

const std::map<std::string, long>& SosSpectrometer::stopCounters() const {
    cntMap_.clear();
    const SosCounters& c = k_->counters();
    for (int i = 1; i <= kNumCounters; ++i) cntMap_[counterName(i)] = c[i];
    return cntMap_;
}

// simc.f:montecarlo l.1468-1472 (hadron arm) / l.1734-1738 (electron arm), l.1491-1500 / 1795-1803 (RECON <- arguments).
SpecResult SosSpectrometer::transport(const SpecInput& in, IRandom& rng) {
    McSosArgs a;
    a.p_spec = in.p_spec;
    a.th_spec = in.th_spec;
    a.dpp = in.delta;
    a.x = in.x;
    a.y = in.y;
    a.z = in.z;
    a.dxdz = in.dxdz;
    a.dydz = in.dydz;
    a.m2 = in.m2;
    a.ms_flag = msOverride_ ? *msOverride_ : in.ms_flag;
    a.wcs_flag = wcsOverride_ ? *wcsOverride_ : in.wcs_flag;
    a.decay_flag = decayOverride_ ? *decayOverride_ : in.decay_flag;
    a.fry = in.raster_fry;   // simc.f:1472/1738 pass fry (raster), not xtar_init
    a.pathlen = 0.0;
    a.ctau = ctau_;
    a.decdist = 0.0;        // simc.f:177
    a.Mh2_final = in.m2;    // simc.f:174 (Mh2 for the hadron arm)
    const McSosDetail det = k_->mcSos(a, rng);

    SpecResult r;
    r.ok = a.ok_spec;
    r.stop_reason = det.stop_code == kSuccess ? std::string() : std::string(counterName(det.stop_code));
    if (det.fit_done) {
        r.x_fp = a.x_fp;
        r.dx_fp = a.dx_fp;
        r.y_fp = a.y_fp;
        r.dy_fp = a.dy_fp;
    }
    if (r.ok) {
        r.delta_recon = a.dpp;
        r.dxdz_recon = a.dxdz;
        r.dydz_recon = a.dydz;
        r.y_recon = a.y;
    }
    r.pathlen = a.pathlen;
    if (det.reached_hut) r.resolution_factor = a.resmult;   // mc_sos_hut writes resmult
    r.m2_final = a.m2;
    r.decay_distance = a.decdist;
    r.Mh2_final = a.Mh2_final;
    return r;
}

// Fortran: trg_track.f:827-837 -> mc_sos_recon(delta,dy,dx,y,xxd) with the COMMON /track/ hut values of the last call
void SosSpectrometer::reconstructLastTrack(double fry, double& delta, double& dydz, double& dxdz, double& y) {
    k_->engine().recon(2, k_->track(), fry, delta, dydz, dxdz, y);
}

}  // namespace sos
}  // namespace simc

SIMC_REGISTER(simc::ISpectrometer, "sos", simc::sos::SosSpectrometer)
