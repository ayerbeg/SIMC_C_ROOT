/**
 * @file
 * @ingroup spectrometers
 * @brief simc::hms::McHms: the Fortran-faithful HMS single-arm Monte Carlo kernel (mc_hms.f and its hut, collimator, recon).
 *
 * KERNEL header (not installed; included only by src/spectrometers/hms
 * and the L1 tests). The public plug-in is HmsSpectrometer.h (ISpectrometer "hms").
 *
 * McHms: the HMS single-arm Monte Carlo as the Fortran routines (task P6.2; split from the plug-in in P7.4), on top of
 * the shared COSY engine (P6.1). Fortran argument lists and COMMON blocks are mirrored here (McHmsArgs, HmsCounters).
 *   McHms::mcHms             Fortran: hms/mc_hms.f:mc_hms (l.1-441) + hit_dipole (l.445-492)
 *   mcHmsColl (private)      Fortran: hms/mc_hms_coll.f (l.1-147)
 *   mcHmsHut  (private)      Fortran: hms/mc_hms_hut.f  (l.1-608), with cern/lfit.f (simc::cern::lfit)
 *   recon                    Fortran: hms/mc_hms_recon.f -> optics::CosyEngine::recon(1, ...)
 *   pionCollAbsorb           Fortran: hms/pion_coll_absorb.f (l.1-95)
 * Constants: hms/apertures_hms.inc, the PARAMETERs of mc_hms.f / mc_hms_coll.f / mc_hms_hut.f (copied verbatim).
 *
 * Units and coordinates as in the Fortran: TRANSPORT (x down, y left, z downstream), cm, slopes, dp/p in %,
 * momenta MeV/c, m2 MeV^2. The ray enters at the plane z = 0 through the target centre (simc.f:montecarlo
 * l.1446-1449 drifts it there); `fry` must be what simc.f passes, i.e. xtar_init = the x at z = 0.
 *
 * State that mirrors Fortran COMMON / SAVE (one McHms object = one Fortran process):
 *   * counters()   : COMMON /spectrometer1/ (struct_hms.inc), 22 integer*4 counters hSTOP_*;
 *   * track()      : COMMON /track/ (spectrometers.inc) after the last call (where the ray stopped);
 *   * the COSY maps (transp_init on the first call, mc_hms.f:185-190; recon file on the first recon);
 *   * mc_hms_coll's static thick_temp (-fno-automatic), enerloss_new's statics (inside the owned SimcMaterial).
 * The event-level COMMON variables that project/transp use for decay (ctau, decdist, Mh2_final) are passed
 * explicitly through McHmsArgs (simc.f:174-177 resets decdist and Mh2_final at the start of every event).
 *
 * Random numbers (same order and number as the Fortran):
 *   * 1 grnd() per ray that reaches the hut (resmult draw, mc_hms_hut.f:292), even with smearing off;
 *   * ms_flag: musc = 2 gauss1(99), musc_ext = 4 gauss1(99); wcs_flag: 2 gauss1(99) per wire plane (both drawn,
 *     one of them discarded); a success with ms=wcs=on uses 202 gauss1;
 *   * decay: grnd() in project/transp (optics::project); collimator: 1 grnd per hit step for pions
 *     (pion_coll_absorb) + musc + enerloss_new gauss1(10).
 *
 * Fortran behaviour (see the .cpp for file:line and the marker; docs/FORTRAN_DIFFERENCES.md):
 *   * the calorimeter aperture test is commented out (mc_hms_hut.f:582-588): hSTOP_cal is never incremented;
 *   * the second hodoscope planes are tested with the FIRST hodoscope limits/offsets (mc_hms_hut.f:537-551);
 *   * multiple scattering in the Cherenkov exit foil is applied twice (mc_hms_hut.f:523-529);
 *   * lfit ignores DC hits that are exactly 0.0 as float (see Lfit.h);
 *   * the focal-plane values are real*4 (lfit) converted back to real*8, and the recon runs on those;
 *   * FORTRAN-BUG(D-021), corrected unless fortran_compat: mc_hms_coll increments hSTOP_slit_* for every 1/20 step
 *     outside the octagon (several counts per ray, even for rays that finally pass); corrected: at most once per ray
 *     and counter. hSTOP_coll counts the rays absorbed / ranged out (both modes).
 *   * FORTRAN-BUG(D-037 F4), corrected unless fortran_compat: a particle whose energy loss in a collimator step is
 *     clamped by enerloss_new (it ranges out) stops in the collimator (hSTOP_coll); the Fortran lets it continue.
 *   * FORTRAN-BUG(D-037 F7, F9) in project/transp (CosyEngine.cpp): drift remainder after a decay; only pi/K decay.
 *   * FORTRAN-DIFF-CANDIDATE(D-037 F3, F5): collimator musc/eloss regardless of ms_flag; pion_coll_absorb parameters.
 *   * Validated bitwise (fortran_compat) against reference/functions/mc_hms_coll, mc_hms_decay, pion_coll_absorb and
 *     the replay of every recorded mc call of runs/p6f2_pi_{hms,shms}_cap (task P6.4, tests/l1/colldecay.cmake).
 */
#ifndef SIMC_SPECTROMETERS_HMS_KERNEL_MCHMS_H
#define SIMC_SPECTROMETERS_HMS_KERNEL_MCHMS_H

#include <array>
#include <iosfwd>
#include <memory>
#include <string>

#include "simc/material/SimcMaterial.h"
#include "simc/optics/CosyEngine.h"

namespace simc {
namespace hms {

/// @addtogroup spectrometers
/// @{

/// Index of each counter in COMMON /spectrometer1/ (struct_hms.inc l.66-74), 1-based as in the harness tables.
/// The "stop code" of a ray is the index of the counter incremented at the `goto 500` that ended it (0 = success).
enum StopCode : int {
    kSuccess = 0,  ///< the ray passed (no stop).
    kTrials = 1,  ///< counter of the trials.
    kSlitHor = 2,  ///< stop at the SlitHor aperture.
    kSlitVert = 3,  ///< stop at the SlitVert aperture.
    kSlitOct = 4,  ///< stop at the SlitOct aperture.
    kQ1In = 5,  ///< stop at the Q1In aperture.
    kQ1Mid = 6,  ///< stop at the Q1Mid aperture.
    kQ1Out = 7,  ///< stop at the Q1Out aperture.
    kQ2In = 8,  ///< stop at the Q2In aperture.
    kQ2Mid = 9,  ///< stop at the Q2Mid aperture.
    kQ2Out = 10,  ///< stop at the Q2Out aperture.
    kQ3In = 11,  ///< stop at the Q3In aperture.
    kQ3Mid = 12,  ///< stop at the Q3Mid aperture.
    kQ3Out = 13,  ///< stop at the Q3Out aperture.
    kD1In = 14,  ///< stop at the D1In aperture.
    kD1Out = 15,  ///< stop at the D1Out aperture.
    kHut = 16,  ///< stop before the detector hut.
    kDc1 = 17,  ///< stop at drift chamber 1.
    kDc2 = 18,  ///< stop at drift chamber 2.
    kScin = 19,  ///< stop at the scintillators.
    kCal = 20,  ///< stop at the calorimeter.
    kSuccesses = 21,  ///< counter of the successes.
    kColl = 22,  ///< stop in the collimator.
};
/// Number of HMS stop counters.
constexpr int kNumCounters = 22;

/// Key of counter `idx` (1..22) in the summary JSON / .geni (docs/OUTPUT_FORMAT.md "Stop counters"):
/// trials, slit_hor, slit_vert, slit_oct, Q1_in, ..., D1_out, hut, dc1, dc2, scin, cal, successes, coll.
/// @param idx counter index (1..22).
/// @return the key (a string literal).
const char* counterName(int idx);

/// Fortran: COMMON /spectrometer1/ hSTOP_* (integer*4 in the Fortran; long here).
struct HmsCounters {
    std::array<long, kNumCounters + 1> c{};   ///< c[idx], idx = StopCode 1..22 (c[0] unused)
    /// @param idx counter index.
    /// @return the counter (modifiable).
    long& operator[](int idx) { return c.at(static_cast<std::size_t>(idx)); }
    /// @param idx counter index.
    /// @return the counter.
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
    bool ms_flag{false}, wcs_flag{false}, decay_flag{false};  ///< in: multiple scattering, wire-chamber smearing, decay switches.
    double resmult{0};   ///< out when the hut is reached: DC resolution factor (1 or 2)
    double fry{0};       ///< in: simc.f passes xtar_init (cm)
    bool ok_spec{false};   ///< out
    double pathlen{0};   ///< in/out: accumulated (caller sets 0)
    bool using_HMScoll{false};  ///< in: use the collimator instead of the fixed slit.
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
/// @param ppi pion momentum [MeV/c].
/// @param thick collimator thickness traversed [cm].
/// @return the transmission probability.
double pionCollAbsorb(double ppi, double thick);

/// CollTrace bookkeeping of one pion_coll_absorb call / one enerloss_new call in the collimator stepping (shared with
/// the SHMS collimator).
/// @param c the trace to update.
/// @param p momentum [MeV/c].
/// @param thick thickness [cm].
/// @param trans transmission.
/// @param u the uniform random number used.
void traceAbsorb(optics::CollTrace& c, double p, double thick, double trans, double u);
/// Record one enerloss_new call of the collimator stepping.
/// @param c the trace to update.
/// @param len step length [cm].
/// @param epart particle energy [MeV].
/// @param mpart particle mass [MeV].
/// @param eloss energy loss [MeV].
void traceEloss(optics::CollTrace& c, double len, double epart, double mpart, double eloss);

/// Fortran: hms/mc_hms.f:hit_dipole(x,y) (l.445-492): true if (x,y) is OUTSIDE the dipole aperture.
/// @param x,y position [cm].
/// @return true if outside the aperture.
bool hitDipole(double x, double y);

/// @brief The HMS single-arm Monte Carlo with the Fortran COMMON counters and SAVEd state.
///
/// Owns its material model; the COSY engine may be shared. Not thread-safe (one object per worker).
class McHms {
public:
    /// Default directory with hms/forward_cosy.dat and hms/recon_cosy.dat (the repository copies of the Fortran
    /// files, data/matrices/simc_fortran), fixed at build time.
    /// @return the directory.
    static std::string defaultMatrixDir();

    /// Creates the kernel with its own COSY engine (maps loaded on demand).
    McHms();

    /// Fortran: mc_hms(...). Updates `a` like the Fortran arguments, the counters and track().
    /// @param[in,out] a the arguments.
    /// @param rng the random generator.
    /// @return the diagnostics of the call.
    McHmsDetail mcHms(McHmsArgs& a, IRandom& rng);
    /// Load the COSY maps now (the Fortran does it on the first call; loading earlier only moves the file reading).
    void loadMaps() {
        ensureForward();
        ensureRecon();
    }

    /// @return the stop counters.
    const HmsCounters& counters() const { return cnt_; }
    /// Reset all stop counters to 0.
    void resetCounters() { cnt_ = HmsCounters{}; }
    /// COMMON /track/ after the last call.
    /// @return the track.
    const optics::Track& track() const { return track_; }

    /// Share an engine (e.g. with the SHMS); by default the object owns one. Maps are loaded on demand.
    /// @param eng the engine (shared ownership).
    void setEngine(std::shared_ptr<optics::CosyEngine> eng) { eng_ = std::move(eng); }
    /// @return the COSY engine.
    optics::CosyEngine& engine() { return *eng_; }
    /// Directory of the COSY map files (default defaultMatrixDir()).
    /// @param d the directory.
    void setMatrixDir(const std::string& d) { matrixDir_ = d; }
    /// Use the sieve slit.
    /// @param b true to use the sieve.
    void setUseSieve(bool b) { useSieve_ = b; }
    /// D-031 run option, kept in the (possibly shared) COSY engine: true = keep the Fortran bugs of the transport
    /// (truncated raddeg, transp.f:158 kaon decay mass, mc_hms_coll slit counters per step). Default false (corrected).
    /// Call it after setEngine().
    /// @param b true for the Fortran behaviour.
    void setFortranCompat(bool b) { eng_->setFortranCompat(b); }
    /// @return the D-031 run option.
    bool fortranCompat() const { return eng_->fortranCompat(); }
    /// Where the Fortran write(6,*) messages go (default std::cout; nullptr = silent). Also used by the
    /// owned material model and the engine.
    /// @param os the stream (not owned), or nullptr.
    void setLog(std::ostream* os);
    /// Optional diagnostics (optics/TransportTrace.h), reset at the start of every call; nullptr = off (default).
    /// @param t the trace (not owned), or nullptr.
    void setTransportTrace(optics::TransportTrace* t) { transportTrace_ = t; }
    /// Optional collimator diagnostics (reset at the start of every call); nullptr = off.
    /// @param t the trace (not owned), or nullptr.
    void setCollTrace(optics::CollTrace* t) { collTrace_ = t; }
    /// COMMON /save_coll/ save_Hcollsuccess: set to 1.5 by a collimator stop (mc_hms.f), never reset or read in SIMC
    /// (D-037 F10: dead state, mirrored only for the L1 tables).
    /// @return save_Hcollsuccess.
    double saveCollSuccess() const { return saveCollSuccess_; }
    /// @param v save_Hcollsuccess.
    void setSaveCollSuccess(double v) { saveCollSuccess_ = v; }
    /// @return the material model (owned by the kernel).
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

    std::shared_ptr<optics::CosyEngine> eng_;   ///< shared on purpose: several spectrometers may use one engine
    material::SimcMaterial mat_;
    HmsCounters cnt_;
    optics::Track track_;
    std::string matrixDir_;
    std::ostream* log_;
    bool useSieve_{false};
    optics::TransportTrace* transportTrace_{nullptr};
    optics::CollTrace* collTrace_{nullptr};
    double saveCollSuccess_{0};
    double thick_temp_{0};   ///< mc_hms_coll.f local thick_temp (static under -fno-automatic)
};

/// @}

}  // namespace hms
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_HMS_KERNEL_MCHMS_H
