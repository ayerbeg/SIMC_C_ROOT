/**
 * @file
 * @ingroup radiation
 * @brief Kernel state of the radiation code: COMMON /bremcom/, the radc.inc run constants and /radccom/ per-event state.
 *
 * KERNEL header (not installed): included only by src/radiation, the radiation plug-in and the L1 tests.
 *
 * State of the SIMC radiative-correction chain (task P5.1).
 *
 * The Fortran keeps everything in COMMON blocks:
 *   radc.inc  /radccom/  run flags (rad_flag, extrad_flag, doing_tail, ...), run limits (etatzai, Egamma*_max)
 *                        and the per-event quantities (frac, lambda, bt, hardcorfac, c, g, ... Egamma_used/min/max,
 *                        ntail, rad_proton_this_ev)
 *             /radctest/ lim1..3, phot1..3
 *   brem.inc  /bremcom/  produce_output, exponentiate, include_hard, calculate_spence
 * plus SAVEd locals (-fno-automatic) whose stale values matter (D-019):
 *   radc_init_ev: dsoft_prime (stale when intcor_mode = 0), peaked_rad_weight: dsoft_intmin/dsoft_intmax
 *   (stale when rad_flag = 1 and extrad_flag = 3), lambda_dave: warned (DATA .false.).
 *
 * Here the run-constant part is RadcRun (filled once from RunSetup, or directly by the L1 tests), and the
 * per-event part plus the statics is RadcState, owned by the radiator (one instance = one Fortran process).
 *
 * INDEXING: arrays whose Fortran lower bound is 1 (frac(3), lambda(3), bt(2), Egamma_used/min/max(3),
 * lim/phot) are stored with one extra leading element so that the C++ index equals the Fortran index
 * (element [0] is unused and stays 0). Arrays with lower bound 0 (c_int(0:3), c_ext(0:3), c(0:4), g(0:4))
 * are stored as is.
 */
#ifndef SIMC_RADIATION_RADCSTATE_H
#define SIMC_RADIATION_RADCSTATE_H

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

#include "simc/core/RunSetup.h"
#include "simc/core/Errors.h"

namespace simc {
namespace radiation {


/// @addtogroup radiation
/// @{

/// Fortran: brem.inc COMMON /bremcom/ (set by init.f:radc_init l.638-647).
struct BremFlags {
    bool produce_output{false};    ///< = debug(1)
    bool exponentiate{false};      ///< = (use_expon == 1)
    bool include_hard{true};  ///< include the hard-photon correction.
    bool calculate_spence{true};  ///< compute the Spence function (else its approximation).
};

/// Run-constant inputs of the radiation routines (radc.inc flags and limits, simulate.inc values they read).
struct RadcRun {
    // ---- radc.inc, fixed by dbase_read / radc_init / limits_init --------------------------------
    int rad_flag{0};                     ///< dbase.f:228 (0 for doing_phsp)
    int extrad_flag{0};                  ///< after the radc_init default (init.f:621-626)
    int intcor_mode{0};  ///< intcor_mode of the input (radc.inc).
    int use_expon{0};  ///< exponentiate the soft-photon factor.
    int one_tail{0};  ///< tail selection of the input: 0 = all tails, 1..3 = only that tail, -1..-3 = all but that tail (dbase.f:494-499).
    bool use_offshell_rad{false};  ///< use_offshell_rad of the input (radc.inc).
    bool hardwired_rad{false};           ///< dbase.f:521-525
    std::array<bool, 4> doing_tail{};    ///< [1..3] = doing_tail(1:3) (dbase.f:494-513); [0] unused
    double etatzai{0};                   ///< init.f:634
    double Egamma_tot_max{0};            ///< init.f:limits_init
    double Egamma1_max{0} /**< @copydoc Egamma3_max */, Egamma2_max{0} /**< @copydoc Egamma3_max */, Egamma3_max{0};  ///< maximum photon energy per tail [MeV].
    BremFlags brem;                      ///< /bremcom/
    // ---- simulate.inc values read by radc.f / init.f ----------------------------------------------
    double Egamma_gen_max{0};            ///< [MeV]
    double dE_edge_test{0};              ///< [MeV]
    bool doing_eep{false} /**< @copydoc doing_heavy */, doing_hyd_elast{false} /**< @copydoc doing_heavy */, doing_deuterium{false} /**< @copydoc doing_heavy */, doing_heavy{false};  ///< reaction flags read by the radiation code.
    bool doing_pion{false} /**< @copydoc doing_semi */, doing_kaon{false} /**< @copydoc doing_semi */, doing_rho{false} /**< @copydoc doing_semi */, doing_semi{false};  ///< reaction flags read by the radiation code.
    bool doing_delta{false};  ///< doing_delta flag.
    /// D-048 F2 (FORTRAN-BUG): tail-1 Egamma limits of doing_delta set by CASE 4 (pion production); false = the Fortran
    /// (no branch for doing_delta, radc.f:288-293). = doing_delta && !fortran_compat.
    bool delta_tail1_case4{false};
    double targ_A{0};                    ///< targ%A (schwinger)
    double Mh{0} /**< @copydoc Mh2 */, Mh2{0};                ///< [MeV], [MeV^2]
    double spec_e_P{0} /**< @copydoc spec_p_P */, spec_p_P{0};     ///< spec%e%P, spec%p%P [MeV/c]
    Cut edge_e_E /**< @copydoc edge_Em */, edge_p_E /**< @copydoc edge_Em */, edge_Em;     ///< edge%e%E, edge%p%E, edge%Em [MeV]
    Cut VERTEXedge_Em /**< @copydoc VERTEXedge_Trec */, VERTEXedge_Pm /**< @copydoc VERTEXedge_Trec */, VERTEXedge_Trec;  ///< vertex cuts on Em, Pm, Trec [MeV].
    Cut gen_sumEgen;                     ///< gen%sumEgen [MeV]

    /// Take every field from a built RunSetup (BuildRunSetup already applied radc_init).
    /// @param s the run set-up.
    /// @return the run constants.
    static RadcRun FromSetup(const RunSetup& s);
};

/// Memo of std::tgamma (task P15.4, performance): radc_init_ev/basicrad_init_ev and peaked_rad_weight call gamma() about
/// 29 times per try on 12 distinct arguments (1+lambda(i), 1+bt(i), 1+g(i), ...), and libm tgamma was 20-56% of the run
/// time. tgamma is a pure function, so returning the value computed earlier for the SAME argument (compared bit for bit,
/// so -0/+0 and NaN payloads are distinct keys) gives bitwise-identical results. 16 slots, round-robin replacement: the
/// values of one try stay in it until the try's last use. One instance per RadcState (per radiator, per worker thread).
struct TgammaMemo {
    static constexpr int kSlots = 16;  ///< number of memo slots.
    std::array<std::uint64_t, kSlots> key{};  ///< bit patterns of the arguments.
    std::array<double, kSlots> val{};  ///< tgamma values.
    std::array<bool, kSlots> used{};  ///< slot in use.
    int next{0};  ///< next slot to replace.
    /// @param x the argument.
    /// @return std::tgamma(x), bit-identical.
    double operator()(double x) {
        std::uint64_t k;
        std::memcpy(&k, &x, sizeof k);
        for (int i = 0; i < kSlots; ++i)
            if (used[i] && key[i] == k) return val[i];
        const double y = std::tgamma(x);
        key[next] = k;
        val[next] = y;
        used[next] = true;
        next = (next + 1) % kSlots;
        return y;
    }
};

/// Per-event COMMON /radccom/ + /radctest/ quantities and the SAVEd locals of the radiation routines.
struct RadcState {
    // ---- /radccom/ per-event part (radc_init_ev, basicrad_init_ev) --------------------------------
    double etta{0};  ///< eta (radc_init_ev).
    std::array<double, 4> frac{};        ///< frac(1:3)
    std::array<double, 4> lambda{};      ///< lambda(1:3)
    std::array<double, 3> bt{};          ///< bt(1:2)
    double hardcorfac{0};  ///< hard-photon correction factor.
    std::array<double, 4> c_int{};       ///< c_int(0:3)
    std::array<double, 4> c_ext{};       ///< c_ext(0:3)
    std::array<double, 5> c{};           ///< c(0:4)
    double g_int{0} /**< @copydoc g_ext */, g_ext{0};  ///< internal and external parts of g.
    std::array<double, 5> g{};           ///< g(0:4)
    bool rad_proton_this_ev{false};  ///< the hadron radiates in this event.
    // ---- /radccom/ generate_rad part --------------------------------------------------------------
    std::array<double, 4> Egamma_used{}; ///< Egamma_used(1:3)
    std::array<double, 4> Egamma_min{};  ///< Egamma_min(1:3): NOT reset by generate_rad (stale entries kept)
    std::array<double, 4> Egamma_max{};  ///< Egamma_max(1:3): idem
    int ntail{0};  ///< tail chosen for this event (0 = all).
    // ---- /radctest/ ----------------------------------------------------------------------------------
    std::array<double, 4> lim{};         ///< lim1..lim3
    std::array<double, 4> phot{};        ///< phot1..phot3
    // ---- SAVEd locals (-fno-automatic) -----------------------------------------------------------------
    double initev_dsoft{0} /**< @copydoc initev_dsoft_prime */, initev_dhard{0} /**< @copydoc initev_dsoft_prime */, initev_dsoft_prime{0};   ///< init.f:radc_init_ev locals
    double prw_dsoft_intmin{0} /**< @copydoc prw_dsoft_intmax */, prw_dsoft_intmax{0};                  ///< radc.f:peaked_rad_weight locals
    bool lambda_warned{false};                                        ///< radc.f:lambda_dave DATA warned
    // ---- not Fortran state: the gamma() memo (P15.4); mutable because extrad_phi reads a const RadcState ------------
    mutable TgammaMemo gamma;  ///< the tgamma memo.
};

/// @}

}  // namespace radiation
}  // namespace simc

#endif  // SIMC_RADIATION_RADCSTATE_H
