// src/spectrometers/calo/CaloSpectrometer.cpp — see include/simc/spectrometers/calo/CaloSpectrometer.h.
#include "simc/spectrometers/calo/CaloSpectrometer.h"

#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/spectrometers/calo/kernel/McCalo.h"

namespace simc {
namespace calo {

CaloSpectrometer::CaloSpectrometer() : k_(std::make_unique<McCalo>()) {}
CaloSpectrometer::~CaloSpectrometer() = default;

void CaloSpectrometer::configure(const Json& config, const SimcParameters& params) {
    drift_ = params.drift_to_cal;   // dbase.f:1104 (default 0.0)
    ctau_ = params.ctau;
    side_ = params.hadron_arm;      // spect_side: unused by mc_calo (D-062 G7)
    if (!config.is_object()) return;
    if (config.contains("drift_to_cal")) drift_ = config.at("drift_to_cal").get<double>();
    if (config.contains("ctau")) ctau_ = config.at("ctau").get<double>();
}

void CaloSpectrometer::init(const StageContext& ctx) { k_->setFortranCompat(ctx.setup.fortran_compat); }

const std::map<std::string, long>& CaloSpectrometer::stopCounters() const {
    cntMap_.clear();
    const CaloCounters& c = k_->counters();
    for (int i = 1; i <= kNumCounters; ++i) cntMap_[counterName(i)] = c[i];
    return cntMap_;
}

// simc.f:1557-1563 (hadron arm, no pi0), 1519-1522/1536-1539 (each pi0 photon), 1761-1765 (electron arm)
SpecResult CaloSpectrometer::transport(const SpecInput& in, IRandom& rng) {
    McCaloArgs a;
    a.p_spec = in.p_spec;
    a.th_spec = in.th_spec;
    a.dpp = in.delta;
    a.x = in.x;
    a.y = in.y;
    a.z = in.z;
    a.dxdz = in.dxdz;
    a.dydz = in.dydz;
    a.m2 = in.m2;
    a.ms_flag = in.ms_flag;
    a.wcs_flag = in.wcs_flag;
    a.decay_flag = in.decay_flag;
    a.frx = in.raster_frx;
    a.fry = in.raster_fry;
    a.pathlen = in.pathlen_start;
    a.spect_side = side_;
    a.drift_to_cal = drift_;
    a.ctau = ctau_;
    a.decdist = 0.0;
    a.Mh2_final = in.m2;
    const McCaloDetail det = k_->mcCalo(a, rng);

    SpecResult r;
    r.ok = a.ok_spec;
    r.stop_reason = det.stop_code == kSuccess ? std::string() : std::string(counterName(det.stop_code));
    if (r.ok) {
        r.x_fp = a.x_fp;
        r.dx_fp = a.dx_fp;
        r.y_fp = a.y_fp;
        r.dy_fp = a.dy_fp;
    }
    // the arguments are returned unchanged (no reconstruction, mc_calo.f:167)
    r.delta_recon = a.dpp;
    r.dxdz_recon = a.dxdz;
    r.dydz_recon = a.dydz;
    r.y_recon = a.y;
    r.pathlen = a.pathlen;
    r.m2_final = a.m2;
    r.decay_distance = a.decdist;
    r.Mh2_final = a.Mh2_final;
    return r;   // resolution_factor: not written by mc_calo (D-062 G7)
}

}  // namespace calo
}  // namespace simc

SIMC_REGISTER(simc::ISpectrometer, "calo", simc::calo::CaloSpectrometer)
