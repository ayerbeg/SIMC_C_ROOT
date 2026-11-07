// include/simc/spectrometers/calo/kernel/McCalo.h — KERNEL header (not installed; included only by
// src/spectrometers/calo and the L1 tests). The public plug-in is CaloSpectrometer.h (ISpectrometer "calo").
//
// McCalo: the calorimeter "spectrometer" of SIMC (arm codes 7 = HMS side, 8 = SHMS side; NPS-like), task P14.3.
//   McCalo::mcCalo   Fortran: calo/mc_calo.f:mc_calo (l.1-191)
//   mcCaloRecon      Fortran: calo/mc_calo_recon.f (l.1-60) — NOT called by SIMC as built (mc_calo.f:167 commented out);
//                    ported as a library function for its L1 table only.
// What the routine does as built (reference/functions/calo/mc_calo/README.txt): copy the ray to COMMON /track/, ONE
// project() by drift_to_cal (decay check only with decay_flag), stop if |ys| > 30.75 (slit_hor) else |xs| > 36.9
// (slit_vert); otherwise x_fp = xs, y_fp = ys, dx_fp = dxdzs, dy_fp = dydzs. dpp/dxdz/dydz/y are returned UNCHANGED (no
// reconstruction). FORTRAN-QUIRK(D-062 G7): ms_flag, wcs_flag, resmult, frx, fry, th_spec, using_tgt_field, prot_zbeam
// and spect_side change nothing; no random numbers unless decay_flag.
// State: COMMON /spectrometer_calo/ caloSTOP_{trials, slit_hor, slit_vert, successes} (struct_calo.inc); COMMON /track/.
#ifndef SIMC_SPECTROMETERS_CALO_KERNEL_MCCALO_H
#define SIMC_SPECTROMETERS_CALO_KERNEL_MCCALO_H

#include <array>
#include <iosfwd>

#include "simc/optics/CosyEngine.h"

namespace simc {
namespace calo {

/// Index of each counter in COMMON /spectrometer_calo/ (struct_calo.inc), 1-based.
enum StopCode : int { kSuccess = 0, kTrials = 1, kSlitHor = 2, kSlitVert = 3, kSuccesses = 4 };
constexpr int kNumCounters = 4;
/// "trials", "slit_hor", "slit_vert", "successes" (the .geni block 'Calo Trials' / 'Extent hor/vert' / 'Successes').
const char* counterName(int idx);

struct CaloCounters {
    std::array<long, kNumCounters + 1> c{};
    long& operator[](int idx) { return c.at(static_cast<std::size_t>(idx)); }
    long operator[](int idx) const { return c.at(static_cast<std::size_t>(idx)); }
};

/// The argument list of mc_calo (mc_calo.f:1-5) plus the decay COMMON variables.
struct McCaloArgs {
    double p_spec{0}, th_spec{0};
    double dpp{0};
    double x{0}, y{0}, z{0};
    double dxdz{0}, dydz{0};
    double x_fp{0}, dx_fp{0}, y_fp{0}, dy_fp{0};   ///< out on success only
    double m2{0};                                  ///< in/out (decay)
    bool ms_flag{false}, wcs_flag{false}, decay_flag{false};
    double resmult{0};   ///< not written (D-062 G7)
    double frx{0}, fry{0};
    bool ok_spec{false};
    double pathlen{0};   ///< in/out: accumulated
    bool using_tgt_field{false};
    double prot_zbeam{0};
    int spect_side{0};
    double drift_to_cal{0};
    double ctau{0}, decdist{0}, Mh2_final{0};
};

struct McCaloDetail {
    int stop_code{kSuccess};
};

/// Fortran: calo/mc_calo_recon.f(delta_p,delta_t,delta_phi,y_tgt,fry,delta_y,delta_z,drift_to_cal) on COMMON /track/.
void mcCaloRecon(const optics::Track& t, double fry, double delta_y, double delta_z, double drift_to_cal,
                 double& delta_p, double& delta_t, double& delta_phi, double& y_tgt);

class McCalo {
public:
    McCalo() = default;
    /// Fortran: mc_calo(...).
    McCaloDetail mcCalo(McCaloArgs& a, IRandom& rng);
    const CaloCounters& counters() const { return cnt_; }
    const optics::Track& track() const { return track_; }
    /// D-031 run option (project.f decay corrections D-037 F7/F9).
    void setFortranCompat(bool b) { compat_ = b; }
    bool fortranCompat() const { return compat_; }

private:
    CaloCounters cnt_;
    optics::Track track_;
    bool compat_{false};
};

}  // namespace calo
}  // namespace simc

#endif  // SIMC_SPECTROMETERS_CALO_KERNEL_MCCALO_H
