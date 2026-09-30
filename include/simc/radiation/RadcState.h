// include/simc/radiation/RadcState.h
//
// State of the SIMC radiative-correction chain (task P5.1).
//
// The Fortran keeps everything in COMMON blocks:
//   radc.inc  /radccom/  run flags (rad_flag, extrad_flag, doing_tail, ...), run limits (etatzai, Egamma*_max)
//                        and the per-event quantities (frac, lambda, bt, hardcorfac, c, g, ... Egamma_used/min/max,
//                        ntail, rad_proton_this_ev)
//             /radctest/ lim1..3, phot1..3
//   brem.inc  /bremcom/  produce_output, exponentiate, include_hard, calculate_spence
// plus SAVEd locals (-fno-automatic) whose stale values matter (D-019):
//   radc_init_ev: dsoft_prime (stale when intcor_mode = 0), peaked_rad_weight: dsoft_intmin/dsoft_intmax
//   (stale when rad_flag = 1 and extrad_flag = 3), lambda_dave: warned (DATA .false.).
//
// Here the run-constant part is RadcRun (filled once from RunSetup, or directly by the L1 tests), and the
// per-event part plus the statics is RadcState, owned by the radiator (one instance = one Fortran process).
//
// INDEXING: arrays whose Fortran lower bound is 1 (frac(3), lambda(3), bt(2), Egamma_used/min/max(3),
// lim/phot) are stored with one extra leading element so that the C++ index equals the Fortran index
// (element [0] is unused and stays 0). Arrays with lower bound 0 (c_int(0:3), c_ext(0:3), c(0:4), g(0:4))
// are stored as is.
#ifndef SIMC_RADIATION_RADCSTATE_H
#define SIMC_RADIATION_RADCSTATE_H

#include <array>
#include <stdexcept>
#include <string>

#include "simc/core/RunSetup.h"

namespace simc {
namespace radiation {

/// A Fortran STOP inside the radiation code (radc.f:215, 697; init.f:628).
class FortranStop : public std::runtime_error {
public:
    explicit FortranStop(const std::string& what) : std::runtime_error(what) {}
};

/// Fortran: brem.inc COMMON /bremcom/ (set by init.f:radc_init l.638-647).
struct BremFlags {
    bool produce_output{false};    ///< = debug(1)
    bool exponentiate{false};      ///< = (use_expon == 1)
    bool include_hard{true};
    bool calculate_spence{true};
};

/// Run-constant inputs of the radiation routines (radc.inc flags and limits, simulate.inc values they read).
struct RadcRun {
    // ---- radc.inc, fixed by dbase_read / radc_init / limits_init --------------------------------
    int rad_flag{0};                     ///< dbase.f:228 (0 for doing_phsp)
    int extrad_flag{0};                  ///< after the radc_init default (init.f:621-626)
    int intcor_mode{0};
    int use_expon{0};
    int one_tail{0};
    bool use_offshell_rad{false};
    bool hardwired_rad{false};           ///< dbase.f:521-525
    std::array<bool, 4> doing_tail{};    ///< [1..3] = doing_tail(1:3) (dbase.f:494-513); [0] unused
    double etatzai{0};                   ///< init.f:634
    double Egamma_tot_max{0};            ///< init.f:limits_init
    double Egamma1_max{0}, Egamma2_max{0}, Egamma3_max{0};
    BremFlags brem;                      ///< /bremcom/
    // ---- simulate.inc values read by radc.f / init.f ----------------------------------------------
    double Egamma_gen_max{0};            ///< [MeV]
    double dE_edge_test{0};              ///< [MeV]
    bool doing_eep{false}, doing_hyd_elast{false}, doing_deuterium{false}, doing_heavy{false};
    bool doing_pion{false}, doing_kaon{false}, doing_rho{false}, doing_semi{false};
    double targ_A{0};                    ///< targ%A (schwinger)
    double Mh{0}, Mh2{0};                ///< [MeV], [MeV^2]
    double spec_e_P{0}, spec_p_P{0};     ///< spec%e%P, spec%p%P [MeV/c]
    Cut edge_e_E, edge_p_E, edge_Em;     ///< edge%e%E, edge%p%E, edge%Em [MeV]
    Cut VERTEXedge_Em, VERTEXedge_Pm, VERTEXedge_Trec;
    Cut gen_sumEgen;                     ///< gen%sumEgen [MeV]

    /// Take every field from a built RunSetup (BuildRunSetup already applied radc_init).
    static RadcRun FromSetup(const RunSetup& s);
};

/// Per-event COMMON /radccom/ + /radctest/ quantities and the SAVEd locals of the radiation routines.
struct RadcState {
    // ---- /radccom/ per-event part (radc_init_ev, basicrad_init_ev) --------------------------------
    double etta{0};
    std::array<double, 4> frac{};        ///< frac(1:3)
    std::array<double, 4> lambda{};      ///< lambda(1:3)
    std::array<double, 3> bt{};          ///< bt(1:2)
    double hardcorfac{0};
    std::array<double, 4> c_int{};       ///< c_int(0:3)
    std::array<double, 4> c_ext{};       ///< c_ext(0:3)
    std::array<double, 5> c{};           ///< c(0:4)
    double g_int{0}, g_ext{0};
    std::array<double, 5> g{};           ///< g(0:4)
    bool rad_proton_this_ev{false};
    // ---- /radccom/ generate_rad part --------------------------------------------------------------
    std::array<double, 4> Egamma_used{}; ///< Egamma_used(1:3)
    std::array<double, 4> Egamma_min{};  ///< Egamma_min(1:3): NOT reset by generate_rad (stale entries kept)
    std::array<double, 4> Egamma_max{};  ///< Egamma_max(1:3): idem
    int ntail{0};
    // ---- /radctest/ ----------------------------------------------------------------------------------
    std::array<double, 4> lim{};         ///< lim1..lim3
    std::array<double, 4> phot{};        ///< phot1..phot3
    // ---- SAVEd locals (-fno-automatic) -----------------------------------------------------------------
    double initev_dsoft{0}, initev_dhard{0}, initev_dsoft_prime{0};   ///< init.f:radc_init_ev locals
    double prw_dsoft_intmin{0}, prw_dsoft_intmax{0};                  ///< radc.f:peaked_rad_weight locals
    bool lambda_warned{false};                                        ///< radc.f:lambda_dave DATA warned
};

}  // namespace radiation
}  // namespace simc

#endif  // SIMC_RADIATION_RADCSTATE_H
