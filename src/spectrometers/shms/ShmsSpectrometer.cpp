// src/spectrometers/shms/ShmsSpectrometer.cpp — see include/simc/spectrometers/shms/ShmsSpectrometer.h.
// The ISpectrometer plug-in "shms" around the kernel shms::McShms, and its registration.
#include "simc/spectrometers/shms/ShmsSpectrometer.h"

#include <stdexcept>

#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/SpectrometerId.h"
#include "simc/spectrometers/shms/kernel/McShms.h"

namespace simc {
namespace shms {

ShmsSpectrometer::ShmsSpectrometer() : k_(std::make_unique<McShms>()) {}
ShmsSpectrometer::~ShmsSpectrometer() = default;

// Stage configuration: the run parameters first, then the optional JSON keys of the pipeline (matrix_dir, spectr,
// using_SHMScoll, ctau, and ms_flag/wcs_flag/decay, which override the per-event flags of the event chain; verbose = false silences it).
void ShmsSpectrometer::configure(const Json& config, const SimcParameters& params) {
    usingColl_ = params.using_SHMScoll;
    ctau_ = params.ctau;
    if (!config.is_object()) return;
    if (config.contains("matrix_dir")) k_->setMatrixDir(config.at("matrix_dir").get<std::string>());
    // D-032: the SHMS has one number (5). The former "spectr" key (5, or the retired LSA tune) is gone.
    if (config.contains("spectr"))
        throw std::invalid_argument("shms: the \"spectr\" key was removed (D-032: the SHMS is always spectrometer 5)");
    if (config.contains("using_SHMScoll")) usingColl_ = config.at("using_SHMScoll").get<bool>();
    if (config.contains("ctau")) ctau_ = config.at("ctau").get<double>();
    k_->setOptions(config.value("use_sieve", k_->useSieve()), config.value("use_2sieve", k_->use2Sieve()),
                   config.value("use_coll", k_->useColl()), config.value("skip_hb", k_->skipHb()));
    if (config.contains("ms_flag")) msOverride_ = config.at("ms_flag").get<bool>();
    if (config.contains("wcs_flag")) wcsOverride_ = config.at("wcs_flag").get<bool>();
    if (config.contains("decay")) decayOverride_ = config.at("decay").get<bool>();
    if (config.contains("verbose") && !config.at("verbose").get<bool>()) k_->setLog(nullptr);
}

void ShmsSpectrometer::init(const StageContext& ctx) {
    k_->setFortranCompat(ctx.setup.fortran_compat);   // D-031 run option
    k_->setCompat294f662(ctx.setup.compat294f662());   // task U.R (D-006)
    k_->loadMaps(Code(SpectrometerId::SHMS));
}

const std::map<std::string, long>& ShmsSpectrometer::stopCounters() const {
    cntMap_.clear();
    const ShmsCounters& c = k_->counters();
    for (int i = 1; i <= kNumSlots; ++i)
        if (i != kStopIdSlot) cntMap_[counterName(i)] = c[i];
    cntMap_["q2_mid_dup"] = c[kQ2Mid];   // simc.f:519-520 prints q2_mid in the Q3 'mid' slot
    return cntMap_;
}

// simc.f:montecarlo l.1482-1489 (hadron arm) / l.1747-1755 (electron arm).
SpecResult ShmsSpectrometer::transport(const SpecInput& in, IRandom& rng) {
    McShmsArgs a;
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
    a.spectr = Code(SpectrometerId::SHMS);   // mc_shms spectr argument (D-032: the SHMS is 5)
    a.using_SHMScoll = usingColl_;
    a.ctau = ctau_;
    a.decdist = 0.0;
    a.Mh2_final = in.m2;
    const McShmsDetail det = k_->mcShms(a, rng);

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
    if (det.reached_hut) r.resolution_factor = a.resmult;   // mc_shms_hut writes resmult (1.0)
    r.m2_final = a.m2;
    r.decay_distance = a.decdist;
    r.Mh2_final = a.Mh2_final;
    return r;
}

// Fortran: trg_track.f:827-837 -> mc_shms_recon(delta,dy,dx,y,xxd) with the COMMON /track/ hut values of the last call
void ShmsSpectrometer::reconstructLastTrack(double fry, double& delta, double& dydz, double& dxdz, double& y) {
    k_->engine().recon(Code(SpectrometerId::SHMS), k_->track(), fry, delta, dydz, dxdz, y);
}

}  // namespace shms
}  // namespace simc

SIMC_REGISTER(simc::ISpectrometer, "shms", simc::shms::ShmsSpectrometer)
