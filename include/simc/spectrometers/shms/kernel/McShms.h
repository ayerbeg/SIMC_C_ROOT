// include/simc/spectrometers/shms/kernel/McShms.h — KERNEL header (not installed; included only by
// src/spectrometers/shms and the L1 tests). The public plug-in is ShmsSpectrometer.h (ISpectrometer "shms").
//
// McShms: the SHMS single-arm Monte Carlo as the Fortran routines (task P6.3; split from the plug-in in P7.4), on top
// of the shared COSY engine (P6.1). Fortran argument lists and COMMON blocks are mirrored here (McShmsArgs,
// ShmsCounters, ShmsTrace).
//   McShms::mcShms            Fortran: shms/mc_shms.f:mc_shms (l.1-1108)
//   mcShmsColl (private)      Fortran: shms/mc_shms_coll.f (l.1-149), with hms/pion_coll_absorb.f (hms::pionCollAbsorb)
//   mcShmsHut  (private)      Fortran: shms/mc_shms_hut.f  (l.1-458) + shms/hut.inc, with cern/lfit.f (simc::cern::lfit)
//   recon                     Fortran: shms/mc_shms_recon.f -> optics::CosyEngine::recon(5, ...)
// Constants: shms/apertures_shms.inc, shms/hut.inc and the PARAMETERs / DATA of mc_shms.f (copied verbatim).
//
// Units and coordinates as in the Fortran: TRANSPORT (x down, y left, z downstream), cm, slopes, dp/p in %,
// MeV/c, MeV^2. The ray enters at z = 0 through the target centre; `fry` is what simc.f passes (xtar_init).
//
// State that mirrors Fortran COMMON / SAVE (one McShms object = one Fortran process):
//   * counters()  : COMMON /spectrometer_shms/ (struct_shms.inc), 52 integer*4: the shmsSTOP_* counters, the
//                   COMMON variable stop_id at index 51 (kept separately, see stopId()) and sSTOP_coll;
//   * track()     : COMMON /track/ xs..dpps after the last call; trace(): the other /track/ variables that
//                   mc_shms writes (x_hb_in ... y_d_out aperture positions, xs_fp ... dydzs_fp);
//   * the COSY maps (transp_init on the first call; recon file on the first recon);
//   * mc_shms_coll's static thick_temp; enerloss_new's statics (owned SimcMaterial).
// Decay COMMON variables (ctau, decdist, Mh2_final) are passed explicitly through McShmsArgs.
//
// Random numbers (same order and number as the Fortran):
//   * flags off: none (mc_shms_hut sets resmult = 1.0 and draws nothing);
//   * ms_flag: musc = 2 gauss1(99), musc_ext = 4 gauss1(99); wcs_flag: 2 gauss1(99) per wire plane;
//     202 gauss1 for a success (and for a stop on the fitted-track calorimeter cut), 196 s2, 154 dc2, 84 dc1;
//   * decay: grnd() in project/transp; collimator (using_SHMScoll, pions/muons): as mc_hms_coll.
//
// Reproduced FORTRAN-QUIRKs (see the .cpp for file:line):
//   * the drifts use the hard-coded mc_shms.f lengths (39.47, 60.68, 307.95) that differ from the COSY drift
//     classes (D-015 F5); the ray is projected to the Cherenkov entrance zd_fp + hcer_1_zentrance = 16.25 cm;
//   * stop_id: reset to 0 on every call, set to 1..31 by the magnet apertures and 18..25 by the hut, NOT set by
//     the collimator (slit_* / coll stops leave 0); the hut reuses 18..25, which the magnets also use;
//   * the hut stop counters do not follow their names: the S1 y cut and x cut both increment s1 (stop_id 20/21),
//     the second-hodoscope y cut increments s3 (22) and its x cut s2 (23); `cal` counts two different cuts,
//     the true position at z = 341 (stop_id 24) and the fitted track 5 cm inside (stop_id 25);
//   * never incremented anywhere: targ_hor/vert/oct, BP_in/BP_out, cer1/cer2;
//   * the hut's `spec` ntuple array is real*4 in mc_shms_hut but real*8 in mc_shms (type mismatch through the
//     argument); it is never read or written out, so it has no observable effect and is not ported;
//   * spectr = 6 (SHMS LSA) -> Fortran STOP, but the Fortran checks it only on the FIRST call; the port throws
//     simc::FortranStop for spectr != 5 on every call (a later spectr = 6 call would use unloaded maps).
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

/// Index of each variable in COMMON /spectrometer_shms/ (struct_shms.inc l.89-107), 1-based as in the harness
/// tables. The stop code of a ray is the index of the counter incremented at the stop (0 = success).
enum StopCode : int {
    kSuccess = 0,
    kTrials = 1,
    kSuccesses = 2,
    kTargHor = 3,
    kTargVert = 4,
    kTargOct = 5,
    kSlitHor = 6,
    kSlitVert = 7,
    kSlitOct = 8,
    kHBIn = 9,
    kHBMen = 10,
    kHBMex = 11,
    kHBOut = 12,
    kQ1In = 13,
    kQ1Men = 14,
    kQ1Mid = 15,
    kQ1Mex = 16,
    kQ1Out = 17,
    kQ2In = 18,
    kQ2Men = 19,
    kQ2Mid = 20,
    kQ2Mex = 21,
    kQ2Out = 22,
    kQ3In = 23,
    kQ3Men = 24,
    kQ3Mid = 25,
    kQ3Mex = 26,
    kQ3Out = 27,
    kD1In = 28,
    kD1Flr = 29,
    kD1Men = 30,
    kD1Mid1 = 31,   ///< ... kD1Mid7 = 37
    kD1Mex = 38,
    kD1Out = 39,
    kBPIn = 40,
    kBPOut = 41,
    kHut = 42,
    kDc1 = 43,
    kDc2 = 44,
    kCer1 = 45,
    kCer2 = 46,
    kS1 = 47,
    kS2 = 48,
    kS3 = 49,
    kCal = 50,
    kStopIdSlot = 51,   ///< COMMON variable stop_id (not a counter)
    kColl = 52,         ///< sSTOP_coll
};
constexpr int kNumSlots = 52;

/// Summary/.geni key of slot `idx` (docs/OUTPUT_FORMAT.md): trials, successes, targ_hor, ..., hb_in, ..., q1_in,
/// ..., D1_in, D1_flr, D1_men, D1_mid1..7, D1_mex, D1_out, BP_in, BP_out, hut, dc1, dc2, cer1, cer2, s1, s2, s3,
/// cal, "stop_id" (51), coll.
const char* counterName(int idx);

/// Fortran: COMMON /spectrometer_shms/ (slot 51 = stop_id of the last call).
struct ShmsCounters {
    std::array<long, kNumSlots + 1> c{};   ///< c[idx], idx 1..52 (c[0] unused)
    long& operator[](int idx) { return c.at(static_cast<std::size_t>(idx)); }
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
    bool ms_flag{false}, wcs_flag{false}, decay_flag{false};
    double resmult{0};   ///< out when the hut is entered (always 1.0)
    double fry{0};       ///< in: xtar_init
    bool ok_spec{false};
    double pathlen{0};   ///< in/out
    int spectr{5};       ///< 5 = SHMS (6 = LSA -> STOP)
    bool using_SHMScoll{false};
    double ctau{0}, decdist{0}, Mh2_final{0};   ///< COMMON (simulate.inc), see hms/kernel/McHms.h
};

struct McShmsDetail {
    int stop_code{kSuccess};   ///< slot of the counter incremented at the stop (0 = success)
    int stop_id{0};            ///< COMMON stop_id after the call
    bool reached_hut{false};
    bool fit_done{false};
    float xfp4{0}, dxfp4{0}, yfp4{0}, dyfp4{0};
};

/// The other COMMON /track/ variables written by mc_shms (spectrometers.inc): aperture-plane positions and the
/// track at the hut entrance. Positions of apertures not reached keep the value of an earlier call.
struct ShmsTrace {
    static constexpr int kN = 31;
    static const char* const kNames[kN];   ///< hb_in hb_men hb_mex hb_out q1_in ... q3_out d_in d_flr d_men d_m1..d_m7 d_mex d_out
    std::array<double, kN> x{}, y{};
    double xs_fp{0}, ys_fp{0}, dxdzs_fp{0}, dydzs_fp{0};
};

class McShms {
public:
    static std::string defaultMatrixDir();

    McShms();

    /// Fortran: mc_shms(...).
    McShmsDetail mcShms(McShmsArgs& a, IRandom& rng);
    /// Load the COSY maps of spectrometer `spectr` now (5 = SHMS; 6 = LSA -> FortranStop, mc_shms.f:335-347).
    void loadMaps(int spectr) {
        ensureForward(spectr);
        ensureRecon(spectr);
    }

    const ShmsCounters& counters() const { return cnt_; }
    int stopId() const { return static_cast<int>(cnt_[kStopIdSlot]); }
    void resetCounters() { cnt_ = ShmsCounters{}; }
    const optics::Track& track() const { return track_; }
    const ShmsTrace& trace() const { return trace_; }

    void setEngine(std::shared_ptr<optics::CosyEngine> eng) { eng_ = std::move(eng); }
    optics::CosyEngine& engine() { return *eng_; }
    void setMatrixDir(const std::string& d) { matrixDir_ = d; }
    /// mc_shms.f compile-time options (mock sieves, collimator, skip the HB): defaults false/false/true/false.
    void setOptions(bool use_sieve, bool use_2sieve, bool use_coll, bool skip_hb) {
        useSieve_ = use_sieve;
        use2Sieve_ = use_2sieve;
        useColl_ = use_coll;
        skipHb_ = skip_hb;
    }
    bool useSieve() const { return useSieve_; }
    bool use2Sieve() const { return use2Sieve_; }
    bool useColl() const { return useColl_; }
    bool skipHb() const { return skipHb_; }
    void setLog(std::ostream* os);
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
    double thick_temp_{0};   ///< mc_shms_coll.f local thick_temp (static)
};

}  // namespace shms
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_SHMS_KERNEL_MCSHMS_H
