// include/simc/spectrometers/hms/HmsSpectrometer.h
//
// ISpectrometer "hms": the HMS single-arm Monte Carlo (task P6.2), on top of the shared COSY engine (P6.1).
//   HmsSpectrometer::mcHms   Fortran: hms/mc_hms.f:mc_hms (l.1-441) + hit_dipole (l.445-492)
//   mcHmsColl (private)      Fortran: hms/mc_hms_coll.f (l.1-147)
//   mcHmsHut  (private)      Fortran: hms/mc_hms_hut.f  (l.1-608), with cern/lfit.f (simc::cern::lfit)
//   recon                    Fortran: hms/mc_hms_recon.f -> optics::CosyEngine::recon(1, ...)
//   pionCollAbsorb           Fortran: hms/pion_coll_absorb.f (l.1-95)
// Constants: hms/apertures_hms.inc, the PARAMETERs of mc_hms.f / mc_hms_coll.f / mc_hms_hut.f (copied verbatim).
//
// Units and coordinates as in the Fortran: TRANSPORT (x down, y left, z downstream), cm, slopes, dp/p in %,
// momenta MeV/c, m2 MeV^2. The ray enters at the plane z = 0 through the target centre (simc.f:montecarlo
// l.1446-1449 drifts it there); `fry` must be what simc.f passes, i.e. xtar_init = the x at z = 0.
//
// State that mirrors Fortran COMMON / SAVE (one HmsSpectrometer object = one Fortran process):
//   * counters()   : COMMON /spectrometer1/ (struct_hms.inc), 22 integer*4 counters hSTOP_*;
//   * track()      : COMMON /track/ (spectrometers.inc) after the last call (where the ray stopped);
//   * the COSY maps (transp_init on the first call, mc_hms.f:185-190; recon file on the first recon);
//   * mc_hms_coll's static thick_temp (-fno-automatic), enerloss_new's statics (inside the owned SimcMaterial).
// The event-level COMMON variables that project/transp use for decay (ctau, decdist, Mh2_final) are passed
// explicitly through McHmsArgs (simc.f:174-177 resets decdist and Mh2_final at the start of every event).
//
// Random numbers (same order and number as the Fortran):
//   * 1 grnd() per ray that reaches the hut (resmult draw, mc_hms_hut.f:292), even with smearing off;
//   * ms_flag: musc = 2 gauss1(99), musc_ext = 4 gauss1(99); wcs_flag: 2 gauss1(99) per wire plane (both drawn,
//     one of them discarded); a success with ms=wcs=on uses 202 gauss1;
//   * decay: grnd() in project/transp (optics::project); collimator: 1 grnd per hit step for pions
//     (pion_coll_absorb) + musc + enerloss_new gauss1(10).
//
// Reproduced FORTRAN-QUIRKs (see the .cpp for file:line):
//   * the calorimeter aperture test is commented out (mc_hms_hut.f:582-588): hSTOP_cal is never incremented;
//   * the second hodoscope planes are tested with the FIRST hodoscope limits/offsets (mc_hms_hut.f:537-551);
//   * multiple scattering in the Cherenkov exit foil is applied twice (mc_hms_hut.f:523-529);
//   * lfit ignores DC hits that are exactly 0.0 as float (see Lfit.h);
//   * the focal-plane values are real*4 (lfit) converted back to real*8, and the recon runs on those;
//   * mc_hms_coll increments hSTOP_slit_* for every 1/20 step that is outside the octagon (several counts per
//     ray possible, even for rays that finally pass), and hSTOP_coll when the ray is absorbed / ranges out.
#ifndef SIMC_SPECTROMETERS_HMS_HMSSPECTROMETER_H
#define SIMC_SPECTROMETERS_HMS_HMSSPECTROMETER_H

#include <array>
#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <string>

#include "simc/material/SimcMaterial.h"
#include "simc/optics/CosyEngine.h"
#include "simc/stages/ISpectrometer.h"

namespace simc {
namespace hms {

/// Index of each counter in COMMON /spectrometer1/ (struct_hms.inc l.66-74), 1-based as in the harness tables.
/// The "stop code" of a ray is the index of the counter incremented at the `goto 500` that ended it (0 = success).
enum StopCode : int {
    kSuccess = 0,
    kTrials = 1,
    kSlitHor = 2,
    kSlitVert = 3,
    kSlitOct = 4,
    kQ1In = 5,
    kQ1Mid = 6,
    kQ1Out = 7,
    kQ2In = 8,
    kQ2Mid = 9,
    kQ2Out = 10,
    kQ3In = 11,
    kQ3Mid = 12,
    kQ3Out = 13,
    kD1In = 14,
    kD1Out = 15,
    kHut = 16,
    kDc1 = 17,
    kDc2 = 18,
    kScin = 19,
    kCal = 20,
    kSuccesses = 21,
    kColl = 22,
};
constexpr int kNumCounters = 22;

/// Key of counter `idx` (1..22) in the summary JSON / .geni (docs/OUTPUT_FORMAT.md "Stop counters"):
/// trials, slit_hor, slit_vert, slit_oct, Q1_in, ..., D1_out, hut, dc1, dc2, scin, cal, successes, coll.
const char* counterName(int idx);

/// Fortran: COMMON /spectrometer1/ hSTOP_* (integer*4 in the Fortran; long here).
struct HmsCounters {
    std::array<long, kNumCounters + 1> c{};   ///< c[idx], idx = StopCode 1..22 (c[0] unused)
    long& operator[](int idx) { return c.at(static_cast<std::size_t>(idx)); }
    long operator[](int idx) const { return c.at(static_cast<std::size_t>(idx)); }
};

/// The argument list of mc_hms (mc_hms.f:1-4), in/out exactly as in the Fortran, plus the COMMON variables
/// of the decay code. Fields marked "out, unchanged unless written" keep the caller's value when the routine
/// does not reach the statement that writes them (as the Fortran arguments do).
struct McHmsArgs {
    double p_spec{0};    ///< in: central momentum (MeV/c)
    double th_spec{0};   ///< in: central angle (rad); not used (mc_hms.f:153 xt = th_spec)
    double dpp{0};       ///< in: dp/p (%);                out if ok: reconstructed delta (%)
    double x{0}, y{0}, z{0};   ///< in: cm at z = 0 plane;  out if ok: y = reconstructed ytar (cm)
    double dxdz{0}, dydz{0};   ///< in: slopes;             out if ok: dxdz = recon xptar, dydz = recon yptar
    double x_fp{0}, dx_fp{0}, y_fp{0}, dy_fp{0};   ///< out after both DCs: fitted focal plane (real*4 values)
    double m2{0};        ///< in/out: mass^2 (MeV^2), changed on decay
    bool ms_flag{false}, wcs_flag{false}, decay_flag{false};
    double resmult{0};   ///< out when the hut is reached: DC resolution factor (1 or 2)
    double fry{0};       ///< in: simc.f passes xtar_init (cm)
    bool ok_spec{false};   ///< out
    double pathlen{0};   ///< in/out: accumulated (caller sets 0)
    bool using_HMScoll{false};
    // COMMON (simulate.inc) used by project/transp for decay:
    double ctau{0};        ///< in: c*tau (cm)
    double decdist{0};     ///< in/out: decay distance accumulator (simc.f:177 sets 0 per event)
    double Mh2_final{0};   ///< in/out: set on decay (simc.f:174 sets Mh2 per event)
};

/// Diagnostics of one mc_hms call (what the reference tables record besides the arguments).
struct McHmsDetail {
    int stop_code{kSuccess};   ///< StopCode of the counter incremented at the stop (0 = success)
    bool reached_hut{false};   ///< hSTOP_hut incremented
    bool fit_done{false};      ///< x_fp..dy_fp written (both chambers passed)
    float xfp4{0}, dxfp4{0}, yfp4{0}, dyfp4{0};   ///< the real*4 lfit results (valid if fit_done)
};

/// Fortran: hms/pion_coll_absorb.f:pion_coll_absorb(ppi,thick,transmission) — returns transmission.
double pionCollAbsorb(double ppi, double thick);

/// Fortran: hms/mc_hms.f:hit_dipole(x,y) (l.445-492): true if (x,y) is OUTSIDE the dipole aperture.
bool hitDipole(double x, double y);

class HmsSpectrometer : public ISpectrometer {
public:
    /// Default directory with hms/forward_cosy.dat and hms/recon_cosy.dat (the repository copies of the Fortran
    /// files, data/matrices/simc_fortran), fixed at build time.
    static std::string defaultMatrixDir();

    HmsSpectrometer();

    std::string name() const override { return "hms"; }
    /// JSON (all optional):
    ///   "matrix_dir"   : directory holding hms/forward_cosy.dat, hms/recon_cosy.dat (default defaultMatrixDir())
    ///   "using_HMScoll": bool (default params.using_HMScoll) — step pions/muons through the collimator
    ///   "ctau"         : cm (default params.ctau)
    ///   "ms_flag", "wcs_flag", "decay": bool — if present they OVERRIDE the per-call SpecInput flags
    ///                    (simc.f passes ms_flag = wcs_flag = mc_smear, decay_flag = doing_decay for the hadron arm)
    ///   "use_sieve"    : bool (default false; mc_hms.f:55 mock sieve slit, a Fortran compile-time option)
    void configure(const Json& config, const SimcParameters& params) override;
    /// Loads the COSY maps (the Fortran does it on the first call; loading here only moves the file reading).
    void init(const StageContext& ctx) override;

    SpecResult transport(const SpecInput& in, IRandom& rng) override;
    const std::map<std::string, long>& stopCounters() const override;

    /// Fortran: mc_hms(...). Updates `a` like the Fortran arguments, the counters and track().
    McHmsDetail mcHms(McHmsArgs& a, IRandom& rng);

    const HmsCounters& counters() const { return cnt_; }
    void resetCounters() { cnt_ = HmsCounters{}; }
    /// COMMON /track/ after the last call.
    const optics::Track& track() const { return track_; }

    /// Share an engine (e.g. with the SHMS); by default the object owns one. Maps are loaded on demand.
    void setEngine(std::shared_ptr<optics::CosyEngine> eng) { eng_ = std::move(eng); }
    optics::CosyEngine& engine() { return *eng_; }
    void setMatrixDir(const std::string& d) { matrixDir_ = d; }
    void setUseSieve(bool b) { useSieve_ = b; }
    void setUsingHMScoll(bool b) { usingColl_ = b; }
    void setCtau(double c) { ctau_ = c; }
    /// Where the Fortran write(6,*) messages go (default std::cout; nullptr = silent). Also used by the
    /// owned material model and the engine.
    void setLog(std::ostream* os);
    material::SimcMaterial& materialModel() { return mat_; }

private:
    int body(McHmsArgs& a, optics::DecayState& d, IRandom& rng, McHmsDetail& det);
    bool mcHmsColl(optics::DecayState& d, double p_spec, double& pathlen, IRandom& rng);
    /// Returns 0 (ok_hut) or the StopCode of the hut counter incremented (dc1, dc2, scin).
    int mcHmsHut(optics::DecayState& d, McHmsArgs& a, double zinit, IRandom& rng, McHmsDetail& det);
    void ensureForward();
    void ensureRecon();
    int stop(int code) {
        ++cnt_[code];
        return code;
    }

    std::shared_ptr<optics::CosyEngine> eng_;
    material::SimcMaterial mat_;
    HmsCounters cnt_;
    optics::Track track_;
    mutable std::map<std::string, long> cntMap_;
    std::string matrixDir_;
    std::ostream* log_;
    bool useSieve_{false};
    bool usingColl_{false};
    double ctau_{0};
    std::optional<bool> msOverride_, wcsOverride_, decayOverride_;
    double thick_temp_{0};   ///< mc_hms_coll.f local thick_temp (static under -fno-automatic)
};

}  // namespace hms
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_HMS_HMSSPECTROMETER_H
