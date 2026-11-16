/**
 * @file
 * @ingroup spectrometers
 * @brief simc::shms::McShms: the Fortran-faithful SHMS single-arm Monte Carlo kernel (mc_shms.f and its hut, collimator, recon).
 *
 * KERNEL header (not installed; included only by
 * src/spectrometers/shms and the L1 tests). The public plug-in is ShmsSpectrometer.h (ISpectrometer "shms").
 *
 * McShms: the SHMS single-arm Monte Carlo as the Fortran routines (task P6.3; split from the plug-in in P7.4), on top
 * of the shared COSY engine (P6.1). Fortran argument lists and COMMON blocks are mirrored here (McShmsArgs,
 * ShmsCounters, ShmsTrace).
 *   McShms::mcShms            Fortran: shms/mc_shms.f:mc_shms (l.1-1108)
 *   mcShmsColl (private)      Fortran: shms/mc_shms_coll.f (l.1-149), with hms/pion_coll_absorb.f (hms::pionCollAbsorb)
 *   mcShmsHut  (private)      Fortran: shms/mc_shms_hut.f  (l.1-458) + shms/hut.inc, with cern/lfit.f (simc::cern::lfit)
 *   recon                     Fortran: shms/mc_shms_recon.f -> optics::CosyEngine::recon(5, ...)
 * Constants: shms/apertures_shms.inc, shms/hut.inc and the PARAMETERs / DATA of mc_shms.f (copied verbatim).
 *
 * Units and coordinates as in the Fortran: TRANSPORT (x down, y left, z downstream), cm, slopes, dp/p in %,
 * MeV/c, MeV^2. The ray enters at z = 0 through the target centre; `fry` is what simc.f passes (xtar_init).
 *
 * State that mirrors Fortran COMMON / SAVE (one McShms object = one Fortran process):
 *   * counters()  : COMMON /spectrometer_shms/ (struct_shms.inc), 52 integer*4: the shmsSTOP_* counters, the
 *                   COMMON variable stop_id at index 51 (kept separately, see stopId()) and sSTOP_coll;
 *   * track()     : COMMON /track/ xs..dpps after the last call; trace(): the other /track/ variables that
 *                   mc_shms writes (x_hb_in ... y_d_out aperture positions, xs_fp ... dydzs_fp);
 *   * the COSY maps (transp_init on the first call; recon file on the first recon);
 *   * mc_shms_coll's static thick_temp; enerloss_new's statics (owned SimcMaterial).
 * Decay COMMON variables (ctau, decdist, Mh2_final) are passed explicitly through McShmsArgs.
 *
 * Random numbers (same order and number as the Fortran):
 *   * flags off: none (mc_shms_hut sets resmult = 1.0 and draws nothing);
 *   * ms_flag: musc = 2 gauss1(99), musc_ext = 4 gauss1(99); wcs_flag: 2 gauss1(99) per wire plane;
 *     202 gauss1 for a success (and for a stop on the fitted-track calorimeter cut), 196 s2, 154 dc2, 84 dc1;
 *   * decay: grnd() in project/transp; collimator (using_SHMScoll, pions/muons): as mc_hms_coll.
 *
 * Fortran behaviour (see the .cpp for file:line and the marker; docs/FORTRAN_DIFFERENCES.md):
 *   * the drifts use the hard-coded mc_shms.f lengths (39.47, 60.68, 307.95) that differ from the COSY drift
 *     classes (D-015 F5); the ray is projected to the Cherenkov entrance zd_fp + hcer_1_zentrance = 16.25 cm;
 *   * stop_id: reset to 0 on every call, set to 1..31 by the magnet apertures and 18..25 by the hut, NOT set by
 *     the collimator (slit_* / coll stops leave 0); the hut reuses 18..25, which the magnets also use;
 *   * the hut stop counters do not follow their names: the S1 y cut and x cut both increment s1 (stop_id 20/21),
 *     the second-hodoscope y cut increments s3 (22) and its x cut s2 (23); `cal` counts two different cuts,
 *     the true position at z = 341 (stop_id 24) and the fitted track 5 cm inside (stop_id 25);
 *   * never incremented anywhere: targ_hor/vert/oct, BP_in/BP_out, cer1/cer2;
 *   * the hut's `spec` ntuple array is real*4 in mc_shms_hut but real*8 in mc_shms (type mismatch through the
 *     argument); it is never read or written out, so it has no observable effect and is not ported;
 *   * the spectr argument must be 5 (SHMS, the only SHMS number, D-032). The Fortran STOPs "LSA tune no longer in use!"
 *     for any other value, but only on the FIRST call (mc_shms.f:335-347); the port throws simc::FortranStop on every
 *     call (FORTRAN-BUG(D-025/D-032): a later call with another value would run on unloaded maps). The plug-in always
 *     passes 5 and the input rejects the retired arm code (core/SpectrometerId.h).
 *   * FORTRAN-BUG(D-021/D-037 F1, F4), corrected unless fortran_compat: mc_shms_coll counts slit_* once per ray (the
 *     Fortran per 1/20 step) and stops a particle that ranges out in the collimator (sSTOP_coll).
 *   * Validated bitwise (fortran_compat) against reference/functions/mc_shms_coll, mc_shms_decay and the capture
 *     replays (task P6.4).
 */
#ifndef SIMC_SPECTROMETERS_SHMS_KERNEL_MCSHMS_H
#define SIMC_SPECTROMETERS_SHMS_KERNEL_MCSHMS_H

#include <array>
#include <iosfwd>
#include <memory>
#include <string>

#include "simc/material/SimcMaterial.h"
#include "simc/optics/CosyEngine.h"

namespace simc {
namespace shms {

/// @addtogroup spectrometers
/// @{

/// Index of each variable in COMMON /spectrometer_shms/ (struct_shms.inc l.89-107), 1-based as in the harness
/// tables. The stop code of a ray is the index of the counter incremented at the stop (0 = success).
enum StopCode : int {
    kSuccess = 0,  ///< the ray passed (no stop).
    kTrials = 1,  ///< counter of the trials.
    kSuccesses = 2,  ///< counter of the successes.
    kTargHor = 3,  ///< stop at the TargHor aperture.
    kTargVert = 4,  ///< stop at the TargVert aperture.
    kTargOct = 5,  ///< stop at the TargOct aperture.
    kSlitHor = 6,  ///< stop at the SlitHor aperture.
    kSlitVert = 7,  ///< stop at the SlitVert aperture.
    kSlitOct = 8,  ///< stop at the SlitOct aperture.
    kHBIn = 9,  ///< stop at the HBIn aperture.
    kHBMen = 10,  ///< stop at the HBMen aperture.
    kHBMex = 11,  ///< stop at the HBMex aperture.
    kHBOut = 12,  ///< stop at the HBOut aperture.
    kQ1In = 13,  ///< stop at the Q1In aperture.
    kQ1Men = 14,  ///< stop at the Q1Men aperture.
    kQ1Mid = 15,  ///< stop at the Q1Mid aperture.
    kQ1Mex = 16,  ///< stop at the Q1Mex aperture.
    kQ1Out = 17,  ///< stop at the Q1Out aperture.
    kQ2In = 18,  ///< stop at the Q2In aperture.
    kQ2Men = 19,  ///< stop at the Q2Men aperture.
    kQ2Mid = 20,  ///< stop at the Q2Mid aperture.
    kQ2Mex = 21,  ///< stop at the Q2Mex aperture.
    kQ2Out = 22,  ///< stop at the Q2Out aperture.
    kQ3In = 23,  ///< stop at the Q3In aperture.
    kQ3Men = 24,  ///< stop at the Q3Men aperture.
    kQ3Mid = 25,  ///< stop at the Q3Mid aperture.
    kQ3Mex = 26,  ///< stop at the Q3Mex aperture.
    kQ3Out = 27,  ///< stop at the Q3Out aperture.
    kD1In = 28,  ///< stop at the D1In aperture.
    kD1Flr = 29,  ///< stop at the D1Flr aperture.
    kD1Men = 30,  ///< stop at the D1Men aperture.
    kD1Mid1 = 31,   ///< ... kD1Mid7 = 37
    kD1Mex = 38,  ///< stop at the D1Mex aperture.
    kD1Out = 39,  ///< stop at the D1Out aperture.
    kBPIn = 40,  ///< stop at the BPIn aperture.
    kBPOut = 41,  ///< stop at the BPOut aperture.
    kHut = 42,  ///< stop before the detector hut.
    kDc1 = 43,  ///< stop at drift chamber 1.
    kDc2 = 44,  ///< stop at drift chamber 2.
    kCer1 = 45,  ///< stop at Cherenkov 1.
    kCer2 = 46,  ///< stop at Cherenkov 2.
    kS1 = 47,  ///< stop at scintillator plane S1.
    kS2 = 48,  ///< stop at scintillator plane S2.
    kS3 = 49,  ///< stop at scintillator plane S3.
    kCal = 50,  ///< stop at the calorimeter.
    kStopIdSlot = 51,   ///< COMMON variable stop_id (not a counter)
    kColl = 52,         ///< sSTOP_coll
};
/// Number of counter slots of COMMON /spectrometer_shms/ (including the stop_id slot).
constexpr int kNumSlots = 52;

/// Summary/.geni key of slot `idx` (docs/OUTPUT_FORMAT.md): trials, successes, targ_hor, ..., hb_in, ..., q1_in,
/// ..., D1_in, D1_flr, D1_men, D1_mid1..7, D1_mex, D1_out, BP_in, BP_out, hut, dc1, dc2, cer1, cer2, s1, s2, s3,
/// cal, "stop_id" (51), coll.
/// @param idx slot index (1..52).
/// @return the key (a string literal).
const char* counterName(int idx);

/// Fortran: COMMON /spectrometer_shms/ (slot 51 = stop_id of the last call).
struct ShmsCounters {
    std::array<long, kNumSlots + 1> c{};   ///< c[idx], idx 1..52 (c[0] unused)
    /// @param idx slot index.
    /// @return the counter (modifiable).
    long& operator[](int idx) { return c.at(static_cast<std::size_t>(idx)); }
    /// @param idx slot index.
    /// @return the counter.
    long operator[](int idx) const { return c.at(static_cast<std::size_t>(idx)); }
};

/// Argument list of mc_shms (mc_shms.f:1-4) plus the decay COMMON variables; in/out as in the Fortran.
struct McShmsArgs {
    double p_spec{0};    ///< in (MeV/c)
    double th_spec{0};   ///< in (rad), unused
    double dpp{0};       ///< in (%); out if ok: recon delta
    double x{0}, y{0}, z{0};   ///< in (cm); out if ok: y = recon ytar
    double dxdz{0}, dydz{0};   ///< in; out if ok: dxdz = recon xptar, dydz = recon yptar
    double x_fp{0}, dx_fp{0}, y_fp{0}, dy_fp{0};   ///< out after the fit (real*4 values), else unchanged
    double m2{0};        ///< in/out (decay)
    bool ms_flag{false}, wcs_flag{false}, decay_flag{false};  ///< in: multiple scattering, wire-chamber smearing, decay switches.
    double resmult{0};   ///< out when the hut is entered (always 1.0)
    double fry{0};       ///< in: xtar_init
    bool ok_spec{false};  ///< out: the ray passed.
    double pathlen{0};   ///< in/out
    int spectr{5};       ///< 5 = SHMS (any other value -> FortranStop, mc_shms.f:335-340)
    bool using_SHMScoll{false};  ///< in: use the collimator.
    double ctau{0}, decdist{0}, Mh2_final{0};   ///< COMMON (simulate.inc), see hms/kernel/McHms.h
};

/// Diagnostics of one mc_shms call.
struct McShmsDetail {
    int stop_code{kSuccess};   ///< slot of the counter incremented at the stop (0 = success)
    int stop_id{0};            ///< COMMON stop_id after the call
    bool reached_hut{false};  ///< the hut was entered.
    bool fit_done{false};  ///< x_fp..dy_fp written.
    float xfp4{0}, dxfp4{0}, yfp4{0}, dyfp4{0};  ///< the real*4 lfit results (valid if fit_done).
};

/// The other COMMON /track/ variables written by mc_shms (spectrometers.inc): aperture-plane positions and the
/// track at the hut entrance. Positions of apertures not reached keep the value of an earlier call.
struct ShmsTrace {
    static constexpr int kN = 31;  ///< number of aperture planes recorded.
    static const char* const kNames[kN];   ///< hb_in hb_men hb_mex hb_out q1_in ... q3_out d_in d_flr d_men d_m1..d_m7 d_mex d_out
    std::array<double, kN> x{}, y{};  ///< x and y at each aperture plane [cm].
    double xs_fp{0}, ys_fp{0}, dxdzs_fp{0}, dydzs_fp{0};  ///< track at the hut entrance (positions [cm], slopes).
};

/// @brief The SHMS single-arm Monte Carlo with the Fortran COMMON counters and SAVEd state.
///
/// Owns its material model; the COSY engine may be shared. Not thread-safe (one object per worker).
class McShms {
public:
    /// Default directory with the SHMS COSY maps, fixed at build time.
    /// @return the directory.
    static std::string defaultMatrixDir();

    /// Creates the kernel with its own COSY engine (maps loaded on demand).
    McShms();

    /// Fortran: mc_shms(...).
    /// @param[in,out] a the arguments.
    /// @param rng the random generator.
    /// @return the diagnostics of the call.
    McShmsDetail mcShms(McShmsArgs& a, IRandom& rng);
    /// Load the COSY maps of spectrometer `spectr` now (5 = SHMS; anything else -> FortranStop, mc_shms.f:335-347).
    /// @param spectr spectrometer number (5).
    void loadMaps(int spectr) {
        ensureForward(spectr);
        ensureRecon(spectr);
    }

    /// @return the stop counters.
    const ShmsCounters& counters() const { return cnt_; }
    /// @return COMMON stop_id of the last call.
    int stopId() const { return static_cast<int>(cnt_[kStopIdSlot]); }
    /// Reset all counter slots to 0.
    void resetCounters() { cnt_ = ShmsCounters{}; }
    /// @return COMMON /track/ after the last call.
    const optics::Track& track() const { return track_; }
    /// @return the aperture-plane positions of the last call.
    const ShmsTrace& trace() const { return trace_; }

    /// Share a COSY engine (by default the object owns one). Maps are loaded on demand.
    /// @param eng the engine (shared ownership).
    void setEngine(std::shared_ptr<optics::CosyEngine> eng) { eng_ = std::move(eng); }
    /// @return the COSY engine.
    optics::CosyEngine& engine() { return *eng_; }
    /// Directory of the COSY map files (default defaultMatrixDir()).
    /// @param d the directory.
    void setMatrixDir(const std::string& d) { matrixDir_ = d; }
    /// mc_shms.f compile-time options (mock sieves, collimator, skip the HB): defaults false/false/true/false.
    /// @param use_sieve mock sieve.
    /// @param use_2sieve second mock sieve.
    /// @param use_coll collimator.
    /// @param skip_hb skip the HB magnet.
    void setOptions(bool use_sieve, bool use_2sieve, bool use_coll, bool skip_hb) {
        useSieve_ = use_sieve;
        use2Sieve_ = use_2sieve;
        useColl_ = use_coll;
        skipHb_ = skip_hb;
    }
    /// @return the use_sieve option.
    bool useSieve() const { return useSieve_; }
    /// D-031 run option, kept in the (possibly shared) COSY engine: true = keep the Fortran bugs of the transport
    /// (truncated raddeg, transp.f:158 kaon decay mass, mc_hms_coll slit counters per step). Default false (corrected).
    /// Call it after setEngine().
    /// @param b true for the Fortran behaviour.
    void setFortranCompat(bool b) { eng_->setFortranCompat(b); }
    /// @return the D-031 run option.
    bool fortranCompat() const { return eng_->fortranCompat(); }
    /// @return the use_2sieve option.
    bool use2Sieve() const { return use2Sieve_; }
    /// @return the use_coll option.
    bool useColl() const { return useColl_; }
    /// @return the skip_hb option.
    bool skipHb() const { return skipHb_; }
    /// Where the Fortran write(6,*) messages go (default std::cout; nullptr = silent).
    /// @param os the stream (not owned), or nullptr.
    void setLog(std::ostream* os);
    /// Optional diagnostics (optics/TransportTrace.h), reset at the start of every call; nullptr = off (default).
    /// @param t the trace (not owned), or nullptr.
    void setTransportTrace(optics::TransportTrace* t) { transportTrace_ = t; }
    /// Optional collimator diagnostics (reset at the start of every call); nullptr = off.
    /// @param t the trace (not owned), or nullptr.
    void setCollTrace(optics::CollTrace* t) { collTrace_ = t; }
    /// COMMON /save_coll/ save_Scollsuccess: set to 1.5 by a collimator stop (mc_shms.f), never reset or read in SIMC
    /// (D-037 F10: dead state, mirrored only for the L1 tables).
    /// @return save_Scollsuccess.
    double saveCollSuccess() const { return saveCollSuccess_; }
    /// @param v save_Scollsuccess.
    void setSaveCollSuccess(double v) { saveCollSuccess_ = v; }
    /// @return the material model (owned by the kernel).
    material::SimcMaterial& materialModel() { return mat_; }

private:
    int body(McShmsArgs& a, optics::DecayState& d, IRandom& rng, McShmsDetail& det);
    bool mcShmsColl(optics::DecayState& d, double p_spec, double& pathlen, IRandom& rng);
    /// Returns 0 (ok_hut) or the slot of the counter of the hut stop (stop_id set as in the Fortran).
    int mcShmsHut(optics::DecayState& d, McShmsArgs& a, double zinit, IRandom& rng, McShmsDetail& det);
    void ensureForward(int spectr);
    void ensureRecon(int spectr);
    int stop(int code, int id) {
        ++cnt_[code];
        cnt_[kStopIdSlot] = id;
        return code;
    }
    int stop(int code) {   // stop_id unchanged (collimator stops)
        ++cnt_[code];
        return code;
    }

    std::shared_ptr<optics::CosyEngine> eng_;   ///< shared on purpose: several spectrometers may use one engine
    material::SimcMaterial mat_;
    ShmsCounters cnt_;
    optics::Track track_;
    ShmsTrace trace_;
    std::string matrixDir_;
    std::ostream* log_;
    bool useSieve_{false}, use2Sieve_{false}, useColl_{true}, skipHb_{false};
    // COMMON /hutflag/ cer_flag, vac_flag: mc_shms.f:257-258 sets .true./.false. on every call.
    bool cerFlag_{true}, vacFlag_{false};
    optics::TransportTrace* transportTrace_{nullptr};
    optics::CollTrace* collTrace_{nullptr};
    double saveCollSuccess_{0};
    double thick_temp_{0};   ///< mc_shms_coll.f local thick_temp (static)
};

/// @}

}  // namespace shms
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_SHMS_KERNEL_MCSHMS_H
