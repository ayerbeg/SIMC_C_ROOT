// src/spectrometers/hms/HmsSpectrometer.cpp — see include/simc/spectrometers/hms/HmsSpectrometer.h.
// The ISpectrometer plug-in "hms" around the kernel hms::McHms, and its registration.
#include "simc/spectrometers/hms/HmsSpectrometer.h"

#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/spectrometers/hms/kernel/McHms.h"

namespace simc {
namespace hms {

HmsSpectrometer::HmsSpectrometer() : k_(std::make_unique<McHms>()) {}
HmsSpectrometer::~HmsSpectrometer() = default;

// Stage configuration: the run parameters first, then the optional JSON keys of the pipeline (matrix_dir, use_sieve,
// and ms_flag/wcs_flag/decay, which override the per-event flags of the event chain; verbose = false silences it).
void HmsSpectrometer::configure(const Json& config, const SimcParameters& params) {
    usingColl_ = params.using_HMScoll;
    ctau_ = params.ctau;
    if (!config.is_object()) return;
    if (config.contains("matrix_dir")) k_->setMatrixDir(config.at("matrix_dir").get<std::string>());
    if (config.contains("using_HMScoll")) usingColl_ = config.at("using_HMScoll").get<bool>();
    if (config.contains("ctau")) ctau_ = config.at("ctau").get<double>();
    if (config.contains("use_sieve")) k_->setUseSieve(config.at("use_sieve").get<bool>());
    if (config.contains("ms_flag")) msOverride_ = config.at("ms_flag").get<bool>();
    if (config.contains("wcs_flag")) wcsOverride_ = config.at("wcs_flag").get<bool>();
    if (config.contains("decay")) decayOverride_ = config.at("decay").get<bool>();
    if (config.contains("verbose") && !config.at("verbose").get<bool>()) k_->setLog(nullptr);
}

void HmsSpectrometer::init(const StageContext& ctx) {
    k_->setFortranCompat(ctx.setup.fortran_compat);   // D-031 run option
    k_->setCompat294f662(ctx.setup.compat294f662());   // task U.R (D-006)
    k_->loadMaps();
}

const std::map<std::string, long>& HmsSpectrometer::stopCounters() const {
    cntMap_.clear();
    const HmsCounters& c = k_->counters();
    for (int i = 1; i <= kNumCounters; ++i) cntMap_[counterName(i)] = c[i];
    return cntMap_;
}

// The ISpectrometer contract: copy the event chain's input into the Fortran argument list of the kernel, run it, and
// return the focal-plane and reconstructed values (only valid where the Fortran had computed them).
// simc.f:montecarlo l.1459-1467 (m2 = Mh2, pathlen = 0, ms = wcs = mc_smear, decay = doing_decay,
// resmult -> ntup%resfac, fry argument = xtar_init) and l.1491-1500 (RECON <- arguments).
SpecResult HmsSpectrometer::transport(const SpecInput& in, IRandom& rng) {
    McHmsArgs a;
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
    a.fry = in.fry;
    a.pathlen = 0.0;
    a.using_HMScoll = usingColl_;
    a.ctau = ctau_;
    a.decdist = 0.0;        // simc.f:177
    a.Mh2_final = in.m2;    // simc.f:174 (Mh2 for the hadron arm)
    const McHmsDetail det = k_->mcHms(a, rng);

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
    if (det.reached_hut) r.resolution_factor = a.resmult;   // mc_hms_hut writes resmult
    r.m2_final = a.m2;
    r.decay_distance = a.decdist;
    r.Mh2_final = a.Mh2_final;
    return r;
}

// Fortran: trg_track.f:827-837 -> mc_hms_recon(delta,dy,dx,y,xxd) with the COMMON /track/ hut values of the last call
void HmsSpectrometer::reconstructLastTrack(double fry, double& delta, double& dydz, double& dxdz, double& y) {
    k_->engine().recon(1, k_->track(), fry, delta, dydz, dxdz, y);
}

}  // namespace hms
}  // namespace simc

SIMC_REGISTER(simc::ISpectrometer, "hms", simc::hms::HmsSpectrometer)
