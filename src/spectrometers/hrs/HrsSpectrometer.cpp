// src/spectrometers/hrs/HrsSpectrometer.cpp — see include/simc/spectrometers/hrs/HrsSpectrometer.h.
#include "simc/spectrometers/hrs/HrsSpectrometer.h"

#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/spectrometers/hrs/kernel/McHrs.h"

namespace simc {
namespace hrs {

HrsSpectrometer::HrsSpectrometer(Side side) : k_(std::make_unique<McHrs>(side)) {}
HrsSpectrometer::~HrsSpectrometer() = default;
HrsrSpectrometer::HrsrSpectrometer() : HrsSpectrometer(Side::Right) {}
HrslSpectrometer::HrslSpectrometer() : HrsSpectrometer(Side::Left) {}

std::string HrsSpectrometer::name() const { return k_->params().name; }

void HrsSpectrometer::configure(const Json& config, const SimcParameters& params) {
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

void HrsSpectrometer::init(const StageContext& ctx) {
    k_->setFortranCompat(ctx.setup.fortran_compat);   // D-031 run option (D-061 F1, H2)
    k_->loadMaps();
}

const std::map<std::string, long>& HrsSpectrometer::stopCounters() const {
    cntMap_.clear();
    const HrsCounters& c = k_->counters();
    for (int i = 1; i <= kNumCounters; ++i) cntMap_[counterName(i)] = c[i];
    return cntMap_;
}

// simc.f:montecarlo l.1473-1482 (hadron arm) / l.1739-1748 (electron arm).
SpecResult HrsSpectrometer::transport(const SpecInput& in, IRandom& rng) {
    McHrsArgs a;
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
    a.fry = in.raster_fry;   // simc.f passes fry (raster), not xtar_init
    a.pathlen = 0.0;
    a.ctau = ctau_;
    a.decdist = 0.0;
    a.Mh2_final = in.m2;
    const McHrsDetail det = k_->mcHrs(a, rng);

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
    if (det.reached_hut) r.resolution_factor = a.resmult;   // mc_hrs*_hut writes resmult = 1.0
    r.m2_final = a.m2;
    r.decay_distance = a.decdist;
    r.Mh2_final = a.Mh2_final;
    return r;
}

}  // namespace hrs
}  // namespace simc

SIMC_REGISTER(simc::ISpectrometer, "hrsr", simc::hrs::HrsrSpectrometer)
SIMC_REGISTER(simc::ISpectrometer, "hrsl", simc::hrs::HrslSpectrometer)
