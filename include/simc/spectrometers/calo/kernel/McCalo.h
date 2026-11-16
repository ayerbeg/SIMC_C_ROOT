/**
 * @file
 * @ingroup spectrometers
 * @brief simc::calo::McCalo: the Fortran-faithful calorimeter arm kernel (calo/mc_calo.f, mc_calo_recon.f).
 *
 * KERNEL header (not installed; included only by
 * src/spectrometers/calo and the L1 tests). The public plug-in is CaloSpectrometer.h (ISpectrometer "calo").
 *
 * McCalo: the calorimeter "spectrometer" of SIMC (arm codes 7 = HMS side, 8 = SHMS side; NPS-like), task P14.3.
 *   McCalo::mcCalo   Fortran: calo/mc_calo.f:mc_calo (l.1-191)
 *   mcCaloRecon      Fortran: calo/mc_calo_recon.f (l.1-60) — NOT called by SIMC as built (mc_calo.f:167 commented out);
 *                    ported as a library function for its L1 table only.
 * What the routine does as built (reference/functions/calo/mc_calo/README.txt): copy the ray to COMMON /track/, ONE
 * project() by drift_to_cal (decay check only with decay_flag), stop if |ys| > 30.75 (slit_hor) else |xs| > 36.9
 * (slit_vert); otherwise x_fp = xs, y_fp = ys, dx_fp = dxdzs, dy_fp = dydzs. dpp/dxdz/dydz/y are returned UNCHANGED (no
 * reconstruction). FORTRAN-QUIRK(D-062 G7): ms_flag, wcs_flag, resmult, frx, fry, th_spec, using_tgt_field, prot_zbeam
 * and spect_side change nothing; no random numbers unless decay_flag.
 * State: COMMON /spectrometer_calo/ caloSTOP_{trials, slit_hor, slit_vert, successes} (struct_calo.inc); COMMON /track/.
 */
#ifndef SIMC_SPECTROMETERS_CALO_KERNEL_MCCALO_H
#define SIMC_SPECTROMETERS_CALO_KERNEL_MCCALO_H

#include <array>
#include <iosfwd>

#include "simc/optics/CosyEngine.h"

namespace simc {
namespace calo {

/// @addtogroup spectrometers
/// @{

/// Index of each counter in COMMON /spectrometer_calo/ (struct_calo.inc), 1-based.
enum StopCode : int {
    kSuccess = 0,     ///< the ray passed (no stop).
    kTrials = 1,      ///< counter of the trials.
    kSlitHor = 2,     ///< stop outside the horizontal extent of the calorimeter.
    kSlitVert = 3,    ///< stop outside the vertical extent of the calorimeter.
    kSuccesses = 4    ///< counter of the successes.
};
/// Number of calorimeter stop counters.
constexpr int kNumCounters = 4;
/// "trials", "slit_hor", "slit_vert", "successes" (the .geni block 'Calo Trials' / 'Extent hor/vert' / 'Successes').
/// @param idx counter index (1..4).
/// @return the key (a string literal).
const char* counterName(int idx);

/// Fortran: COMMON /spectrometer_calo/ counters (long here).
struct CaloCounters {
    std::array<long, kNumCounters + 1> c{};  ///< c[idx], idx = StopCode 1..4 (c[0] unused).
    /// @param idx counter index.
    /// @return the counter (modifiable).
    long& operator[](int idx) { return c.at(static_cast<std::size_t>(idx)); }
    /// @param idx counter index.
    /// @return the counter.
    long operator[](int idx) const { return c.at(static_cast<std::size_t>(idx)); }
};

/// The argument list of mc_calo (mc_calo.f:1-5) plus the decay COMMON variables.
struct McCaloArgs {
    double p_spec{0}, th_spec{0};  ///< in: central momentum [MeV/c] and angle [rad].
    double dpp{0};  ///< in: dp/p [%]; out: reconstructed delta [%].
    double x{0}, y{0}, z{0};  ///< in: position [cm]; out: reconstructed ytar in y [cm].
    double dxdz{0}, dydz{0};  ///< in: slopes; out: reconstructed xptar, yptar.
    double x_fp{0}, dx_fp{0}, y_fp{0}, dy_fp{0};   ///< out on success only
    double m2{0};                                  ///< in/out (decay)
    bool ms_flag{false}, wcs_flag{false}, decay_flag{false};  ///< in: multiple scattering, smearing, decay switches.
    double resmult{0};   ///< not written (D-062 G7)
    double frx{0}, fry{0};  ///< in: beam (raster) position [cm].
    bool ok_spec{false};  ///< out: the ray hit the calorimeter.
    double pathlen{0};   ///< in/out: accumulated
    bool using_tgt_field{false};  ///< in: polarized-target field in use.
    double prot_zbeam{0};  ///< in: prot_zbeam argument of mc_calo (only in commented-out code of the Fortran).
    int spect_side{0};  ///< in: 7 (HMS side) or 8 (SHMS side).
    double drift_to_cal{0};  ///< in: drift distance from the target to the calorimeter face [cm].
    double ctau{0}, decdist{0}, Mh2_final{0};  ///< COMMON (simulate.inc) of the decay code.
};

/// Diagnostics of one mc_calo call.
struct McCaloDetail {
    int stop_code{kSuccess};  ///< StopCode of the counter incremented (0 = success).
};

/// Fortran: calo/mc_calo_recon.f(delta_p,delta_t,delta_phi,y_tgt,fry,delta_y,delta_z,drift_to_cal) on COMMON /track/.
/// @param t the track at the calorimeter.
/// @param fry vertical beam position [cm].
/// @param delta_y,delta_z vertex offsets used by the reconstruction [cm].
/// @param drift_to_cal drift distance to the calorimeter [cm].
/// @param[out] delta_p reconstructed delta [%].
/// @param[out] delta_t reconstructed yptar.
/// @param[out] delta_phi reconstructed xptar.
/// @param[out] y_tgt reconstructed ytar [cm].
void mcCaloRecon(const optics::Track& t, double fry, double delta_y, double delta_z, double drift_to_cal,
                 double& delta_p, double& delta_t, double& delta_phi, double& y_tgt);

/// @brief The calorimeter arm (smearing and extent test) with the Fortran COMMON counters. Not thread-safe.
class McCalo {
public:
    /// Creates the kernel with zero counters.
    McCalo() = default;
    /// Fortran: mc_calo(...).
    /// @param[in,out] a the arguments.
    /// @param rng the random generator.
    /// @return the diagnostics of the call.
    McCaloDetail mcCalo(McCaloArgs& a, IRandom& rng);
    /// @return the stop counters.
    const CaloCounters& counters() const { return cnt_; }
    /// @return COMMON /track/ after the last call.
    const optics::Track& track() const { return track_; }
    /// D-031 run option (project.f decay corrections D-037 F7/F9).
    /// @param b true for the Fortran behaviour.
    void setFortranCompat(bool b) { compat_ = b; }
    /// @return the D-031 run option.
    bool fortranCompat() const { return compat_; }

private:
    CaloCounters cnt_;
    optics::Track track_;
    bool compat_{false};
};

/// @}

}  // namespace calo
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_CALO_KERNEL_MCCALO_H
