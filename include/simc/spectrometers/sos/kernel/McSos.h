// include/simc/spectrometers/sos/kernel/McSos.h — KERNEL header (not installed; included only by src/spectrometers/sos
// and the L1 tests). The public plug-in is SosSpectrometer.h (ISpectrometer "sos").
//
// McSos: the SOS (Hall C Short Orbit Spectrometer) single-arm Monte Carlo as the Fortran routines (task P14.1), on top
// of the shared COSY engine. Fortran argument lists and COMMON blocks are mirrored here (McSosArgs, SosCounters).
//   McSos::mcSos             Fortran: sos/mc_sos.f:mc_sos (l.1-394)
//   mcSosHut (private)       Fortran: sos/mc_sos_hut.f:mc_sos_hut (l.1-548), with cern/lfit.f (simc::cern::lfit)
//   recon                    Fortran: sos/mc_sos_recon.f -> optics::CosyEngine::recon(2, ...)
// Constants: sos/apertures_sos.inc and the PARAMETERs of mc_sos.f / mc_sos_hut.f (copied verbatim; -fdefault-real-8
// makes every real literal real*8, so they are plain doubles here).
// Matrix files: sos/forward_cosy.dat (transp.f:302, 10 classes, mc_sos.f:158 STOPs otherwise) and sos/recon_cosy.dat
// (mc_sos_recon.f:73); byte copies in data/matrices/simc_fortran/sos (MD5SUMS).
//
// Units and coordinates as in the Fortran: TRANSPORT (x down, y left, z downstream), cm, slopes, dp/p in %,
// momenta MeV/c, m2 MeV^2. The ray enters at the plane z = 0 through the target centre (simc.f:montecarlo
// l.1446-1449). `fry` is what simc.f passes to mc_sos: the RASTER fry (simc.f:1349-1355: -main%target%rastery if
// correct_raster, else 0; simc.f:1472/1738), NOT xtar_init as for HMS/SHMS.
//
// State that mirrors Fortran COMMON / SAVE (one McSos object = one Fortran process):
//   * counters()   : COMMON /spectrometer2/ (struct_sos.inc), 18 integer*4 counters sSTOP_*;
//   * track()      : COMMON /track/ (spectrometers.inc) after the last call (where the ray stopped);
//   * the COSY maps (transp_init on the first call, mc_sos.f:155-160; recon file on the first recon).
// The event-level COMMON variables that project/transp use for decay (ctau, decdist, Mh2_final) are passed
// explicitly through McSosArgs (simc.f:174-177 resets decdist and Mh2_final at the start of every event).
//
// Random numbers (same order and number as the Fortran):
//   * 1 grnd() per ray that reaches the hut (resmult draw, mc_sos_hut.f:251), even with smearing off;
//   * ms_flag: musc = 2 gauss1(99), musc_ext = 4 gauss1(99); wcs_flag: 2 gauss1(99) per wire plane (both drawn,
//     one of them discarded). A success with ms = wcs = on draws 1 grnd and 180 gauss1 (dc1 stop 66, dc2 stop 132, scin stop 176);
//   * decay: grnd() in project/transp (optics::project, CosyEngine::transp);
//   * the mock sieve (use_sieve, a compile-time .false. in mc_sos.f:50): 2 grnd per call.
//
// Fortran behaviour (see the .cpp for file:line and the marker; docs/FORTRAN_DIFFERENCES.md, D-059 P14.1):
//   * the calorimeter cuts are commented out (mc_sos_hut.f:519-525, 535-539): sSTOP_cal is never incremented;
//   * the last 7.62 cm to the new exit window are a straight line WITHOUT project (mc_sos.f:339-340): no pathlength,
//     no decay check, zs unchanged;
//   * the hut starts at zinit = -3.206267 cm, the exit foil is at -3.22 cm, so the first drift is negative and is
//     clamped to 0.001 cm (mc_sos_hut.f:283-284; the code comment calls the 3.22 cm "double counted");
//   * lfit ignores DC hits that are exactly 0.0 as float (see Lfit.h); the focal-plane values are real*4 (lfit)
//     converted back to real*8, and the recon runs on those;
//   * mc_sos_recon.f replaces EVERY hut value |v| <= 1e-30 by 1e-30 (HMS/SHMS: only fry), see CosyEngine::recon;
//   * the hut is entered even if the particle decayed; x_fp..dy_fp are written only if both chambers are passed.
#ifndef SIMC_SPECTROMETERS_SOS_KERNEL_MCSOS_H
#define SIMC_SPECTROMETERS_SOS_KERNEL_MCSOS_H

#include <array>
#include <iosfwd>
#include <memory>
#include <string>

#include "simc/material/SimcMaterial.h"
#include "simc/optics/CosyEngine.h"

namespace simc {
namespace sos {

/// Index of each counter in COMMON /spectrometer2/ (struct_sos.inc l.36-52), 1-based in declaration order.
/// The "stop code" of a ray is the index of the counter incremented at the `goto 500` that ended it (0 = success).
enum StopCode : int {
    kSuccess = 0,
    kTrials = 1,
    kSlitHor = 2,
    kSlitVert = 3,
    kSlitOct = 4,
    kQuadIn = 5,
    kQuadMid = 6,
    kQuadOut = 7,
    kBm01In = 8,
    kBm01Out = 9,
    kBm02In = 10,
    kBm02Out = 11,
    kExit = 12,
    kHut = 13,
    kDc1 = 14,
    kDc2 = 15,
    kScin = 16,
    kCal = 17,
    kSuccesses = 18,
};
constexpr int kNumCounters = 18;

/// Key of counter `idx` (1..18) in the summary JSON / .geni: trials, slit_hor, slit_vert, slit_oct, quad_in, quad_mid,
/// quad_out, bm01_in, bm01_out, bm02_in, bm02_out, exit, hut, dc1, dc2, scin, cal, successes (the sSTOP_* names).
const char* counterName(int idx);

/// Fortran: COMMON /spectrometer2/ sSTOP_* (integer*4 in the Fortran; long here).
struct SosCounters {
    std::array<long, kNumCounters + 1> c{};   ///< c[idx], idx = StopCode 1..18 (c[0] unused)
    long& operator[](int idx) { return c.at(static_cast<std::size_t>(idx)); }
    long operator[](int idx) const { return c.at(static_cast<std::size_t>(idx)); }
};

/// The argument list of mc_sos (mc_sos.f:1-3), in/out exactly as in the Fortran, plus the COMMON variables of the
/// decay code. Fields marked "out" keep the caller's value when the routine does not reach the statement that
/// writes them (as the Fortran arguments do).
struct McSosArgs {
    double p_spec{0};    ///< in: central momentum (MeV/c)
    double th_spec{0};   ///< in: central angle (rad); not used (mc_sos.f:118 xt = th_spec)
    double dpp{0};       ///< in: dp/p (%);                out if ok: reconstructed delta (%)
    double x{0}, y{0}, z{0};   ///< in: cm at z = 0 plane;  out if ok: y = reconstructed ytar (cm)
    double dxdz{0}, dydz{0};   ///< in: slopes;             out if ok: dxdz = recon xptar, dydz = recon yptar
                               ///< (with use_sieve also changed before the transport, mc_sos.f:135-136)
    double x_fp{0}, dx_fp{0}, y_fp{0}, dy_fp{0};   ///< out after both DCs: fitted focal plane (real*4 values)
    double m2{0};        ///< in/out: mass^2 (MeV^2), changed on decay
    bool ms_flag{false}, wcs_flag{false}, decay_flag{false};
    double resmult{0};   ///< out when the hut is reached: DC resolution factor (1 or 2)
    double fry{0};       ///< in: simc.f passes the raster fry (cm), see the file header
    bool ok_spec{false};   ///< out
    double pathlen{0};   ///< in/out: accumulated (caller sets 0)
    // COMMON (simulate.inc) used by project/transp for decay:
    double ctau{0};        ///< in: c*tau (cm)
    double decdist{0};     ///< in/out: decay distance accumulator (simc.f:177 sets 0 per event)
    double Mh2_final{0};   ///< in/out: set on decay (simc.f:174 sets Mh2 per event)
};

/// Diagnostics of one mc_sos call (what the reference tables record besides the arguments).
struct McSosDetail {
    int stop_code{kSuccess};   ///< StopCode of the counter incremented at the stop (0 = success)
    bool reached_hut{false};   ///< sSTOP_hut incremented
    bool fit_done{false};      ///< x_fp..dy_fp written (both chambers passed)
    float xfp4{0}, dxfp4{0}, yfp4{0}, dyfp4{0};   ///< the real*4 lfit results (valid if fit_done)
};

class McSos {
public:
    /// Default directory with sos/forward_cosy.dat and sos/recon_cosy.dat (the repository copies of the Fortran
    /// files, data/matrices/simc_fortran), fixed at build time.
    static std::string defaultMatrixDir();

    McSos();

    /// Fortran: mc_sos(...). Updates `a` like the Fortran arguments, the counters and track().
    McSosDetail mcSos(McSosArgs& a, IRandom& rng);
    /// Load the COSY maps now (the Fortran does it on the first call; loading earlier only moves the file reading).
    void loadMaps() {
        ensureForward();
        ensureRecon();
    }

    const SosCounters& counters() const { return cnt_; }
    void resetCounters() { cnt_ = SosCounters{}; }
    /// COMMON /track/ after the last call.
    const optics::Track& track() const { return track_; }

    /// Share an engine (e.g. with the HMS); by default the object owns one. Maps are loaded on demand.
    void setEngine(std::shared_ptr<optics::CosyEngine> eng) { eng_ = std::move(eng); }
    optics::CosyEngine& engine() { return *eng_; }
    void setMatrixDir(const std::string& d) { matrixDir_ = d; }
    /// mc_sos.f:50 `logical use_sieve /.false./` (mock sieve slit; a compile-time option in the Fortran).
    void setUseSieve(bool b) { useSieve_ = b; }
    /// D-031 run option, kept in the (possibly shared) COSY engine: true = keep the Fortran bugs of the transport
    /// (truncated raddeg of rotate_haxis, transp.f:158 kaon decay mass, project.f decay remainder). Default false.
    void setFortranCompat(bool b) { eng_->setFortranCompat(b); }
    bool fortranCompat() const { return eng_->fortranCompat(); }
    /// Where the Fortran write(6,*) messages go (default std::cout; nullptr = silent). Also used by the
    /// owned material model and the engine.
    void setLog(std::ostream* os);
    /// Optional diagnostics (optics/TransportTrace.h), reset at the start of every call; nullptr = off (default).
    void setTransportTrace(optics::TransportTrace* t) { transportTrace_ = t; }
    material::SimcMaterial& materialModel() { return mat_; }

    /// L1 hook (reference/functions/sos/hut_sos): mc_sos_hut alone. The ray is the COMMON /track/ (setTrack); `d`
    /// carries m2, p, dflag, decay_flag and the decay COMMON; `a` the flags, resmult, pathlen and the focal-plane outputs.
    /// Returns 0 (ok_hut) or the StopCode of the hut counter incremented.
    int mcSosHutDirect(optics::DecayState& d, McSosArgs& a, double zinit, IRandom& rng, McSosDetail& det) {
        return mcSosHut(d, a, zinit, rng, det);
    }
    void setTrack(const optics::Track& t) { track_ = t; }

private:
    int body(McSosArgs& a, optics::DecayState& d, IRandom& rng, McSosDetail& det);
    /// Returns 0 (ok_hut) or the StopCode of the hut counter incremented (dc1, dc2, scin).
    int mcSosHut(optics::DecayState& d, McSosArgs& a, double zinit, IRandom& rng, McSosDetail& det);
    void ensureForward();
    void ensureRecon();
    int stop(int code) {
        ++cnt_[code];
        return code;
    }

    std::shared_ptr<optics::CosyEngine> eng_;
    material::SimcMaterial mat_;
    SosCounters cnt_;
    optics::Track track_;
    std::string matrixDir_;
    std::ostream* log_;
    bool useSieve_{false};
    optics::TransportTrace* transportTrace_{nullptr};
};

}  // namespace sos
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_SOS_KERNEL_MCSOS_H
