// include/simc/spectrometers/hrs/kernel/McHrs.h — KERNEL header (not installed; included only by src/spectrometers/hrs
// and the L1 tests). The public plug-ins are in HrsSpectrometer.h (ISpectrometer "hrsr", "hrsl").
//
// McHrs: the Hall A High Resolution Spectrometers (HRS-R = arm 3, HRS-L = arm 4) single-arm Monte Carlo as the Fortran
// routines (task P14.2), ONE implementation with a side parameter. A line-by-line diff of hrsl/ and hrsr/ shows that
// the two differ only in (mc_hrsl.f vs mc_hrsr.f):
//   * spectr 4 / 3 (l.31), the counter COMMON /spectrometer4/ lSTOP_* / /spectrometer3/ rSTOP_* and /g_dump_all_l/ /
//     /g_dump_all_r/ (same layout);
//   * h_exit 3.335 / 3.340 (l.58), z_entr 110.9d0 / 110.0d0 (l.67);
//   * the VDC y offset subtracted from y_fp after the recon input is set: 0.78 / 0.48 (mc_hrsl.f:525 / mc_hrsr.f:524);
// the hut routines (mc_hrsl_hut.f / mc_hrsr_hut.f) and the recon routines differ only in names, and the COSY files
// (hrs_forward_cosy.dat, hrs_recon_cosy.dat) are byte-identical: one copy in data/matrices/simc_fortran/hrs/.
//   McHrs::mcHrs            Fortran: hrsl/mc_hrsl.f:mc_hrsl (l.1-551), hrsr/mc_hrsr.f:mc_hrsr
//   mcHrsHut (private)      Fortran: hrsl/mc_hrsl_hut.f (l.1-498), hrsr/mc_hrsr_hut.f, with cern/lfit.f
//   recon                   Fortran: mc_hrs*_recon.f -> optics::CosyEngine::recon(3|4, ...)
//
// Units and coordinates as in the Fortran: TRANSPORT (x down, y left, z downstream), cm, slopes, dp/p in %, MeV/c, MeV^2.
// `fry` is the RASTER fry of simc.f (l.1349-1355; simc.f:1477/1482/1743/1748), as for the SOS.
//
// State that mirrors Fortran COMMON / SAVE (one McHrs object = one Fortran process for one side):
//   * counters()   : COMMON /spectrometer3|4/ (struct_hrs*.inc), 21 integer*4 counters;
//   * dump()       : COMMON /g_dump_all_r|l/ stop_where, x_stop, y_stop (g_dump_all_events.inc): written only at a stop,
//                    otherwise they keep the value of an earlier call (mirrored, not used by SIMC);
//   * track()      : COMMON /track/ after the last call;
//   * the COSY maps (first call).
//
// Random numbers: none per call without smearing/decay (the hut sets resmult = 1.0 itself, mc_hrs*_hut.f:243);
// ms = wcs = on: success 176 gauss1(99), dc1 stop 70, dc2 140, s1 146; decay: grnd in project/transp.
//
// REAL(16) quirks (-fdefault-real-8 promotes the d0 literals to binary128; reference/functions/hrsl/*/README.txt):
//   * FORTRAN-BUG(D-061 F1): mc_hrs*.f:248, 259, 307, 318, 388, 455, 466 pass transport lengths like 62.75333333d0 to the
//     real*8 dummy zd of transp, which reads the first 8 bytes of the 16-byte constant: zd = garbage (-1.3e-41 ...
//     -3.1e+232). transp uses zd for the pathlength (FP path = -3.134e+232 for every ray past Q2), the decay test and the
//     length warning (printed at every call). Corrected: the intended lengths. fortran_compat: the garbage values.
//   * FORTRAN-BUG(D-061 F1e): mc_hrs*_hut.f:293, 356 call rotate_haxis(45.0d0,...): the dummy reads +0.0 (low 8 bytes of
//     the quad 45.0), so the VDC acceptance test is done WITHOUT the 45 degree rotation. Corrected: 45 degrees.
//   * z_entr = 110.9d0 / 110.0d0 + z_off and z_exit = z_entr + 8.0d0 are PARAMETERs declared real*8: the quad values are
//     rounded to real*8 (no garbage); mc_hrs*_recon.f:118 compares with the quad 1.d-30 (CosyEngine::recon).
// Other Fortran behaviour (FORTRAN_DIFFERENCES.md): the two circular apertures before the slit box count as slit_hor
// (stop_where 20, 21); the calorimeter cuts are commented out (hSTOP_cal never incremented); y_fp carries the VDC offset,
// the recon does not.
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

/// The side of the HRS pair: Right = HRS-R (SIMC arm 3), Left = HRS-L (arm 4).
enum class Side { Right, Left };

/// Index of each counter in COMMON /spectrometer3|4/ (struct_hrs*.inc), 1-based in declaration order.
enum StopCode : int {
    kSuccess = 0,
    kTrials = 1,
    kSlitHor = 2,
    kSlitVert = 3,
    kQ1In = 4,
    kQ1Mid = 5,
    kQ1Out = 6,
    kQ2In = 7,
    kQ2Mid = 8,
    kQ2Out = 9,
    kQ3In = 10,
    kQ3Mid = 11,
    kQ3Out = 12,
    kD1In = 13,
    kD1Out = 14,
    kHut = 15,
    kDc1 = 16,
    kDc2 = 17,
    kS1 = 18,
    kS2 = 19,
    kCal = 20,
    kSuccesses = 21,
};
constexpr int kNumCounters = 21;

/// Key of counter `idx` (1..21): trials, slit_hor, slit_vert, Q1_in, Q1_mid, Q1_out, Q2_in, Q2_mid, Q2_out, Q3_in,
/// Q3_mid, Q3_out, D1_in, D1_out, hut, dc1, dc2, s1, s2, cal, successes.
const char* counterName(int idx);

struct HrsCounters {
    std::array<long, kNumCounters + 1> c{};
    long& operator[](int idx) { return c.at(static_cast<std::size_t>(idx)); }
    long operator[](int idx) const { return c.at(static_cast<std::size_t>(idx)); }
};

/// COMMON /g_dump_all_r|l/ (g_dump_all_events.inc): where the last stopped ray stopped.
struct HrsDump {
    double stop_where{0}, x_stop{0}, y_stop{0};
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
const HrsSideParams& sideParams(Side s);

/// The transp length arguments of mc_hrs*.f (classes 2, 3, 5, 6, 8, 10, 11): the intended value (the decimal literal) and
/// what transp reads from the REAL(16) constant (its low 8 bytes, D-061 F1).
struct TranspLength {
    int cls;
    double intended;
    double compat;
};
const std::array<TranspLength, 7>& transpLengths();
/// binary128 (gfortran REAL(16) = the d0 literals under -fdefault-real-8).
__extension__ typedef __float128 Quad;
/// What a real*8 dummy reads from the binary128 constant `q` (x86-64: the first 8 bytes = the low mantissa bits).
double low8(Quad q);

/// The argument list of mc_hrs*(...) (mc_hrsl.f:1-3), as McSosArgs.
struct McHrsArgs {
    double p_spec{0}, th_spec{0};
    double dpp{0};
    double x{0}, y{0}, z{0};
    double dxdz{0}, dydz{0};
    double x_fp{0}, dx_fp{0}, y_fp{0}, dy_fp{0};
    double m2{0};
    bool ms_flag{false}, wcs_flag{false}, decay_flag{false};
    double resmult{0};
    double fry{0};
    bool ok_spec{false};
    double pathlen{0};
    double ctau{0}, decdist{0}, Mh2_final{0};   ///< COMMON (simulate.inc) of the decay code
};

struct McHrsDetail {
    int stop_code{kSuccess};
    bool reached_hut{false};
    bool fit_done{false};
};

class McHrs {
public:
    static std::string defaultMatrixDir();

    explicit McHrs(Side side);

    /// Fortran: mc_hrsr / mc_hrsl(...).
    McHrsDetail mcHrs(McHrsArgs& a, IRandom& rng);
    void loadMaps() {
        ensureForward();
        ensureRecon();
    }

    Side side() const { return side_; }
    const HrsSideParams& params() const { return sideParams(side_); }
    const HrsCounters& counters() const { return cnt_; }
    const HrsDump& dump() const { return dump_; }
    void setDump(const HrsDump& d) { dump_ = d; }
    const optics::Track& track() const { return track_; }
    void setTrack(const optics::Track& t) { track_ = t; }

    void setEngine(std::shared_ptr<optics::CosyEngine> eng) { eng_ = std::move(eng); }
    optics::CosyEngine& engine() { return *eng_; }
    void setMatrixDir(const std::string& d) { matrixDir_ = d; }
    /// mc_hrs*.f:47 `logical use_sieve /.false./`.
    void setUseSieve(bool b) { useSieve_ = b; }
    void setFortranCompat(bool b) { eng_->setFortranCompat(b); }
    bool fortranCompat() const { return eng_->fortranCompat(); }
    void setLog(std::ostream* os);
    void setTransportTrace(optics::TransportTrace* t) { transportTrace_ = t; }

    /// L1 hook (reference/functions/hrs*/hut_hrs*): mc_hrs*_hut alone on track().
    int mcHrsHutDirect(optics::DecayState& d, McHrsArgs& a, double zinit, IRandom& rng, McHrsDetail& det) {
        return mcHrsHut(d, a, zinit, rng, det);
    }
    /// The zd passed to transp for class `cls` (D-061 F1: garbage with fortran_compat).
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

}  // namespace hrs
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_HRS_KERNEL_MCHRS_H
