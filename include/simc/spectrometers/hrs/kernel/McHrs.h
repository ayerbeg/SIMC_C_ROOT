/**
 * @file
 * @ingroup spectrometers
 * @brief simc::hrs::McHrs: the Fortran-faithful HRS-R / HRS-L single-arm Monte Carlo kernel (mc_hrsr.f, mc_hrsl.f).
 *
 * KERNEL header (not installed; included only by src/spectrometers/hrs
 * and the L1 tests). The public plug-ins are in HrsSpectrometer.h (ISpectrometer "hrsr", "hrsl").
 *
 * McHrs: the Hall A High Resolution Spectrometers (HRS-R = arm 3, HRS-L = arm 4) single-arm Monte Carlo as the Fortran
 * routines (task P14.2), ONE implementation with a side parameter. A line-by-line diff of hrsl/ and hrsr/ shows that
 * the two differ only in (mc_hrsl.f vs mc_hrsr.f):
 *   * spectr 4 / 3 (l.31), the counter COMMON /spectrometer4/ lSTOP_* / /spectrometer3/ rSTOP_* and /g_dump_all_l/ /
 *     /g_dump_all_r/ (same layout);
 *   * h_exit 3.335 / 3.340 (l.58), z_entr 110.9d0 / 110.0d0 (l.67);
 *   * the VDC y offset subtracted from y_fp after the recon input is set: 0.78 / 0.48 (mc_hrsl.f:525 / mc_hrsr.f:524);
 * the hut routines (mc_hrsl_hut.f / mc_hrsr_hut.f) and the recon routines differ only in names, and the COSY files
 * (hrs_forward_cosy.dat, hrs_recon_cosy.dat) are byte-identical: one copy in data/matrices/simc_fortran/hrs/.
 *   McHrs::mcHrs            Fortran: hrsl/mc_hrsl.f:mc_hrsl (l.1-551), hrsr/mc_hrsr.f:mc_hrsr
 *   mcHrsHut (private)      Fortran: hrsl/mc_hrsl_hut.f (l.1-498), hrsr/mc_hrsr_hut.f, with cern/lfit.f
 *   recon                   Fortran: mc_hrs*_recon.f -> optics::CosyEngine::recon(3|4, ...)
 *
 * Units and coordinates as in the Fortran: TRANSPORT (x down, y left, z downstream), cm, slopes, dp/p in %, MeV/c, MeV^2.
 * `fry` is the RASTER fry of simc.f (l.1349-1355; simc.f:1477/1482/1743/1748), as for the SOS.
 *
 * State that mirrors Fortran COMMON / SAVE (one McHrs object = one Fortran process for one side):
 *   * counters()   : COMMON /spectrometer3|4/ (struct_hrs*.inc), 21 integer*4 counters;
 *   * dump()       : COMMON /g_dump_all_r|l/ stop_where, x_stop, y_stop (g_dump_all_events.inc): written only at a stop,
 *                    otherwise they keep the value of an earlier call (mirrored, not used by SIMC);
 *   * track()      : COMMON /track/ after the last call;
 *   * the COSY maps (first call).
 *
 * Random numbers: none per call without smearing/decay (the hut sets resmult = 1.0 itself, mc_hrs*_hut.f:243);
 * ms = wcs = on: success 176 gauss1(99), dc1 stop 70, dc2 140, s1 146; decay: grnd in project/transp.
 *
 * REAL(16) quirks (-fdefault-real-8 promotes the d0 literals to binary128; reference/functions/hrsl/(case)/README.txt):
 *   * FORTRAN-BUG(D-061 F1): mc_hrs*.f:248, 259, 307, 318, 388, 455, 466 pass transport lengths like 62.75333333d0 to the
 *     real*8 dummy zd of transp, which reads the first 8 bytes of the 16-byte constant: zd = garbage (-1.3e-41 ...
 *     -3.1e+232). transp uses zd for the pathlength (FP path = -3.134e+232 for every ray past Q2), the decay test and the
 *     length warning (printed at every call). Corrected: the intended lengths. fortran_compat: the garbage values.
 *   * FORTRAN-BUG(D-061 F1e): mc_hrs*_hut.f:293, 356 call rotate_haxis(45.0d0,...): the dummy reads +0.0 (low 8 bytes of
 *     the quad 45.0), so the VDC acceptance test is done WITHOUT the 45 degree rotation. Corrected: 45 degrees.
 *   * z_entr = 110.9d0 / 110.0d0 + z_off and z_exit = z_entr + 8.0d0 are PARAMETERs declared real*8: the quad values are
 *     rounded to real*8 (no garbage); mc_hrs*_recon.f:118 compares with the quad 1.d-30 (CosyEngine::recon).
 * Other Fortran behaviour (FORTRAN_DIFFERENCES.md): the two circular apertures before the slit box count as slit_hor
 * (stop_where 20, 21); the calorimeter cuts are commented out (hSTOP_cal never incremented); y_fp carries the VDC offset,
 * the recon does not.
 */
#ifndef SIMC_SPECTROMETERS_HRS_KERNEL_MCHRS_H
#define SIMC_SPECTROMETERS_HRS_KERNEL_MCHRS_H

#include <array>
#include <iosfwd>
#include <memory>
#include <string>

#include "simc/material/SimcMaterial.h"
#include "simc/optics/CosyEngine.h"

namespace simc {
namespace hrs {

/// @addtogroup spectrometers
/// @{

/// The side of the HRS pair: Right = HRS-R (SIMC arm 3), Left = HRS-L (arm 4).
enum class Side { Right /**< HRS-R (arm 3) */, Left /**< HRS-L (arm 4) */ };

/// Index of each counter in COMMON /spectrometer3|4/ (struct_hrs*.inc), 1-based in declaration order.
enum StopCode : int {
    kSuccess = 0,  ///< the ray passed (no stop).
    kTrials = 1,  ///< counter of the trials.
    kSlitHor = 2,  ///< stop at the SlitHor aperture.
    kSlitVert = 3,  ///< stop at the SlitVert aperture.
    kQ1In = 4,  ///< stop at the Q1In aperture.
    kQ1Mid = 5,  ///< stop at the Q1Mid aperture.
    kQ1Out = 6,  ///< stop at the Q1Out aperture.
    kQ2In = 7,  ///< stop at the Q2In aperture.
    kQ2Mid = 8,  ///< stop at the Q2Mid aperture.
    kQ2Out = 9,  ///< stop at the Q2Out aperture.
    kQ3In = 10,  ///< stop at the Q3In aperture.
    kQ3Mid = 11,  ///< stop at the Q3Mid aperture.
    kQ3Out = 12,  ///< stop at the Q3Out aperture.
    kD1In = 13,  ///< stop at the D1In aperture.
    kD1Out = 14,  ///< stop at the D1Out aperture.
    kHut = 15,  ///< stop before the detector hut.
    kDc1 = 16,  ///< stop at drift chamber 1.
    kDc2 = 17,  ///< stop at drift chamber 2.
    kS1 = 18,  ///< stop at scintillator plane S1.
    kS2 = 19,  ///< stop at scintillator plane S2.
    kCal = 20,  ///< stop at the calorimeter.
    kSuccesses = 21,  ///< counter of the successes.
};
/// Number of HRS stop counters.
constexpr int kNumCounters = 21;

/// Key of counter `idx` (1..21): trials, slit_hor, slit_vert, Q1_in, Q1_mid, Q1_out, Q2_in, Q2_mid, Q2_out, Q3_in,
/// Q3_mid, Q3_out, D1_in, D1_out, hut, dc1, dc2, s1, s2, cal, successes.
/// @param idx counter index (1..21).
/// @return the key (a string literal).
const char* counterName(int idx);

/// Fortran: COMMON /spectrometer3|4/ STOP counters (long here).
struct HrsCounters {
    std::array<long, kNumCounters + 1> c{};  ///< c[idx], idx = StopCode 1..21 (c[0] unused).
    /// @param idx counter index.
    /// @return the counter (modifiable).
    long& operator[](int idx) { return c.at(static_cast<std::size_t>(idx)); }
    /// @param idx counter index.
    /// @return the counter.
    long operator[](int idx) const { return c.at(static_cast<std::size_t>(idx)); }
};

/// COMMON /g_dump_all_r|l/ (g_dump_all_events.inc): where the last stopped ray stopped.
struct HrsDump {
    double stop_where{0} /**< @copydoc y_stop */, x_stop{0} /**< @copydoc y_stop */, y_stop{0};  ///< z of the stop [cm] and x, y of the ray there [cm].
};

/// The side constants (the only differences between mc_hrsl.f and mc_hrsr.f).
struct HrsSideParams {
    int spectr;          ///< 3 (R) / 4 (L)
    double h_exit;       ///< 3.340 / 3.335
    double z_entr;       ///< real*8(110.0d0 + z_off) / real*8(110.9d0 + z_off)
    double z_exit;       ///< real*8(z_entr + 8.0d0)
    double vdc_offset;   ///< 0.48 / 0.78
    const char* name;    ///< "hrsr" / "hrsl"
};
/// @param s the side.
/// @return the side constants (static table).
const HrsSideParams& sideParams(Side s);

/// The transp length arguments of mc_hrs*.f (classes 2, 3, 5, 6, 8, 10, 11): the intended value (the decimal literal) and
/// what transp reads from the REAL(16) constant (its low 8 bytes, D-061 F1).
struct TranspLength {
    int cls;  ///< COSY class.
    double intended;  ///< the decimal literal [cm].
    double compat;  ///< what the Fortran transp reads (D-061 F1) [cm].
};
/// @return the seven transp lengths (static table).
const std::array<TranspLength, 7>& transpLengths();
/// binary128 (gfortran REAL(16) = the d0 literals under -fdefault-real-8).
__extension__ typedef __float128 Quad;
/// What a real*8 dummy reads from the binary128 constant `q` (x86-64: the first 8 bytes = the low mantissa bits).
/// @param q the binary128 constant.
/// @return the double formed by its low 8 bytes.
double low8(Quad q);

/// The argument list of mc_hrs*(...) (mc_hrsl.f:1-3), as McSosArgs.
struct McHrsArgs {
    double p_spec{0} /**< @copydoc th_spec */, th_spec{0};  ///< in: central momentum [MeV/c] and angle [rad].
    double dpp{0};  ///< in: dp/p [%]; out if ok: reconstructed delta [%].
    double x{0} /**< @copydoc z */, y{0} /**< @copydoc z */, z{0};  ///< in: position at z = 0 [cm]; out if ok: y = reconstructed ytar [cm].
    double dxdz{0} /**< @copydoc dydz */, dydz{0};  ///< in: slopes; out if ok: reconstructed xptar, yptar.
    double x_fp{0} /**< @copydoc dy_fp */, dx_fp{0} /**< @copydoc dy_fp */, y_fp{0} /**< @copydoc dy_fp */, dy_fp{0};  ///< out after both DCs: fitted focal plane.
    double m2{0};  ///< in/out: mass squared [MeV^2] (changed on decay).
    bool ms_flag{false} /**< @copydoc decay_flag */, wcs_flag{false} /**< @copydoc decay_flag */, decay_flag{false};  ///< in: multiple scattering, wire-chamber smearing, decay switches.
    double resmult{0};  ///< out when the hut is reached: DC resolution factor.
    double fry{0};  ///< in: vertical raster position [cm].
    bool ok_spec{false};  ///< out: the ray passed.
    double pathlen{0};  ///< in/out: accumulated path length [cm].
    double ctau{0} /**< @copydoc Mh2_final */, decdist{0} /**< @copydoc Mh2_final */, Mh2_final{0};   ///< COMMON (simulate.inc) of the decay code
};

/// Diagnostics of one mc_hrs* call.
struct McHrsDetail {
    int stop_code{kSuccess};  ///< StopCode of the counter incremented (0 = success).
    bool reached_hut{false};  ///< the hut was reached.
    bool fit_done{false};  ///< x_fp..dy_fp written.
};

/// @brief The HRS single-arm Monte Carlo of one side with the Fortran COMMON counters and SAVEd state.
///
/// Owns its material model; the COSY engine may be shared. Not thread-safe (one object per worker).
class McHrs {
public:
    /// Default directory with the HRS COSY maps, fixed at build time.
    /// @return the directory.
    static std::string defaultMatrixDir();

    /// @param side the HRS side.
    explicit McHrs(Side side);

    /// Fortran: mc_hrsr / mc_hrsl(...).
    /// @param[in,out] a the arguments.
    /// @param rng the random generator.
    /// @return the diagnostics of the call.
    McHrsDetail mcHrs(McHrsArgs& a, IRandom& rng);
    /// Load the COSY maps now (the Fortran does it on the first call).
    void loadMaps() {
        ensureForward();
        ensureRecon();
    }

    /// @return the side.
    Side side() const { return side_; }
    /// @return the side constants.
    const HrsSideParams& params() const { return sideParams(side_); }
    /// @return the stop counters.
    const HrsCounters& counters() const { return cnt_; }
    /// @return COMMON /g_dump_all_*/.
    const HrsDump& dump() const { return dump_; }
    /// @param d COMMON /g_dump_all_*/.
    void setDump(const HrsDump& d) { dump_ = d; }
    /// @return COMMON /track/ after the last call.
    const optics::Track& track() const { return track_; }
    /// @param t COMMON /track/ (for mcHrsHutDirect()).
    void setTrack(const optics::Track& t) { track_ = t; }

    /// Share a COSY engine (by default the object owns one).
    /// @param eng the engine (shared ownership).
    void setEngine(std::shared_ptr<optics::CosyEngine> eng) { eng_ = std::move(eng); }
    /// @return the COSY engine.
    optics::CosyEngine& engine() { return *eng_; }
    /// @param d directory of the COSY map files.
    void setMatrixDir(const std::string& d) { matrixDir_ = d; }
    /// mc_hrs*.f:47 'logical use_sieve /.false./'.
    /// @param b true to use the mock sieve.
    void setUseSieve(bool b) { useSieve_ = b; }
    /// D-031 run option, kept in the (possibly shared) COSY engine.
    /// @param b true for the Fortran behaviour.
    void setFortranCompat(bool b) { eng_->setFortranCompat(b); }
    /// Task U.R: compat target 294f662 for transp (D-006 kaon decay daughter mass Mk); see CosyEngine::setCompat294f662.
    /// @param b true for the compat target 294f662.
    void setCompat294f662(bool b) { eng_->setCompat294f662(b); }
    /// @return the D-031 run option.
    bool fortranCompat() const { return eng_->fortranCompat(); }
    /// @param os stream for the Fortran messages (not owned), or nullptr.
    void setLog(std::ostream* os);
    /// @param t optional diagnostics (not owned), or nullptr.
    void setTransportTrace(optics::TransportTrace* t) { transportTrace_ = t; }

    /// L1 hook (reference/functions/hrs*/hut_hrs*): mc_hrs*_hut alone on track().
    /// @param d the decay state.
    /// @param a the arguments.
    /// @param zinit zinit as mc_hrs* passes it to the hut routine [cm].
    /// @param rng the random generator.
    /// @param det the diagnostics.
    /// @return 0 (ok_hut) or the StopCode of the hut counter.
    int mcHrsHutDirect(optics::DecayState& d, McHrsArgs& a, double zinit, IRandom& rng, McHrsDetail& det) {
        return mcHrsHut(d, a, zinit, rng, det);
    }
    /// The zd passed to transp for class `cls` (D-061 F1: garbage with fortran_compat).
    /// @param cls COSY class.
    /// @return zd [cm].
    double transpLength(int cls) const;

private:
    int body(McHrsArgs& a, optics::DecayState& d, IRandom& rng, McHrsDetail& det);
    int mcHrsHut(optics::DecayState& d, McHrsArgs& a, double zinit, IRandom& rng, McHrsDetail& det);
    void ensureForward();
    void ensureRecon();
    int stop(int code, double where, double x, double y) {
        ++cnt_[code];
        dump_.stop_where = where;
        dump_.x_stop = x;
        dump_.y_stop = y;
        return code;
    }

    Side side_;
    std::shared_ptr<optics::CosyEngine> eng_;
    material::SimcMaterial mat_;
    HrsCounters cnt_;
    HrsDump dump_;
    optics::Track track_;
    std::string matrixDir_;
    std::ostream* log_;
    bool useSieve_{false};
    optics::TransportTrace* transportTrace_{nullptr};
};

/// @}

}  // namespace hrs
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_HRS_KERNEL_MCHRS_H
