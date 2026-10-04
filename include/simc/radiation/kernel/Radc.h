// include/simc/radiation/kernel/Radc.h
// KERNEL header (not installed): included only by src/radiation, the radiation plug-in and the L1 tests.
//
// The SIMC radiative-correction routines as free functions (task P5.1), ported operation by operation:
//
//   basicrad          radc.f:3-88       BASICRAD photon energy generation (1 grnd)
//   gammaRadc         radc.f:92-116     DEAD CODE in SIMC (D-019): the callers bind to the intrinsic = std::tgamma.
//                                       Ported only for the L1 table; nothing in the chain calls it.
//   generateRad       radc.f:120-519    radiate the event (calls basicrad, complete_ev, peaked_rad_weight/extrad_phi)
//   peakedRadWeight   radc.f:523-646
//   extradFriedrich   radc.f:650-664
//   extradPhi         radc.f:668-707
//   schwinger         radc.f:711-739
//   spen              radc.f:743-764
//   lambdaDave        radc.f:768-821
//   radcInit          init.f:576-651    (the same logic is applied by BuildRunSetup; this copy serves the L1 table)
//   radcInitEv        init.f:655-728
//   basicradInitEv    init.f:732-813
//
// Fortran behaviour kept (D-019, P5.F1 report Q3-Q7; FORTRAN-QUIRK = documented, FORTRAN-DIFF-CANDIDATE = open
// question, docs/FORTRAN_DIFFERENCES.md), marked at the code:
//   * gamma() is libm tgamma everywhere (std::tgamma).
//   * radc_init_ev with intcor_mode = 0 uses the stale SAVEd dsoft_prime (RadcState::initev_dsoft_prime).
//   * peaked_rad_weight: rad_flag = 0 result independent of Egamma; rad_flag = 1 & extrad_flag = 3 uses the stale
//     SAVEd dsoft_intmin/max (RadcState::prw_dsoft_intmin/max).
//   * generate_rad overwrites rad_weight with peaked_rad_weight (rad_flag <= 1); the complete_ev inside tail 1
//     changes hardcorfac/frac/g/c before the final division (radc.f:512).
//   * one_tail != 0: a disabled tail can still be drawn -> no photon, gen_weight = 1/hardcorfac.
//   * doing_tail(3) is passed as doing_proton to lambda_dave (init.f:693).
//   * Egamma_min/max(1:3) are COMMON and are not reset: stale values of tails that are not processed are kept
//     (and, for D(e,e'p) with rad_flag = 3, Egamma_min(1) is decremented by dE_edge_test from its stale value).
//
// RNG (P5.F1 §6): generate_rad with rad_flag <= 1 draws 1 uniform for the tail + 1 in basicrad when the chosen tail
// is enabled and has room; rad_flag = 2: 1 + 1; rad_flag = 3: up to 3. The gaussians come only from complete_ev.
#ifndef SIMC_RADIATION_RADC_H
#define SIMC_RADIATION_RADC_H

#include <iosfwd>

#include "simc/core/Event.h"
#include "simc/radiation/kernel/RadcState.h"
#include "simc/stages/IRadiator.h"

namespace simc {
class IRandom;

namespace radiation {

/// Fortran: radc.f:basicrad (l.3-88). itail is in/out (0 -> 4 on entry, back to 0 on the normal exit only).
/// Reads st.g, st.c. Consumes 1 uniform only when it generates.
void basicrad(const RadcState& st, int& itail, double Egamma_lo, double Egamma_hi, double& Egamma, double& weight,
              double& val_reciprocal, IRandom& rng);

/// Fortran: radc.f:gamma (l.92-116). Dead code in SIMC (see header).
double gammaRadc(double x);

/// Fortran: radc.f:peaked_rad_weight (l.523-646). Reads the per-event COMMON in `st`; updates its SAVEd locals.
double peakedRadWeight(const RadcRun& run, RadcState& st, const Event& vertex, double Egamma, double emin, double emax,
                       double basicrad_val_reciprocal, double basicrad_weight, std::ostream* log = nullptr);

/// Fortran: radc.f:extrad_friedrich (l.650-664). etatzai is the COMMON value.
void extradFriedrich(double etatzai, double Ei, double Ecutoff, double trad, double& dbrem, double& dbrem_prime);

/// Fortran: radc.f:extrad_phi (l.668-707). Reads extrad_flag, etatzai (run) and bt, g (state).
/// Throws FortranStop for itail = 0 with extrad_flag = 3 (radc.f:697).
double extradPhi(const RadcRun& run, const RadcState& st, int itail, double E1, double E2, double Egamma);

/// Fortran: radc.f:schwinger (l.711-739). Reads vertex Q2, nu, e%theta, Ein, e%E; COMMON targ%A, etta, use_expon.
double schwinger(double Ecutoff, const Event& vertex, bool include_hard, double targ_A, double etta, int use_expon,
                 double& dsoft, double& dhard);

/// Fortran: radc.f:spen (l.743-764).
double spen(double x);

/// Fortran: radc.f:lambda_dave (l.768-821). `warned` is the DATA/SAVE flag (warning printed once per process).
double lambdaDave(int itail, int plus_flag, bool doing_proton, double e1, double e2, double e3, double p3, double th,
                  bool& warned, std::ostream* log = nullptr);

/// Outputs of init.f:radc_init that are not simple copies.
struct RadcInitResult {
    int extrad_flag{0};
    double etatzai{0};
    BremFlags brem;
};
/// Fortran: init.f:radc_init (l.576-651). Throws FortranStop for extrad_flag < 0.
RadcInitResult radcInit(double targ_Z, double targ_L1, double targ_L2, int rad_flag, int extrad_flag, int use_expon,
                        bool debug1);

/// Fortran: init.f:radc_init_ev (l.655-728). Reads main.target.teff(1:2) and the vertex; writes the state.
void radcInitEv(const RadcRun& run, RadcState& st, const MainEvent& main, const Event& vertex,
                std::ostream* log = nullptr);

/// Fortran: init.f:basicrad_init_ev (l.732-813). Reads st.lambda, st.bt, st.g(4); writes g(0:3), c_*, g_*.
void basicradInitEv(RadcState& st, double e1, double e2, double e3);

/// Fortran: radc.f:generate_rad (l.120-519). vertex%Ein is reduced by the tail-1 photon and the vertex is
/// recompleted by `complete_ev` (event.f:complete_ev, which must call radcInitEv on this same state).
/// Returns success; main.gen_weight is multiplied on success; main.ntup.radphot/radarm are set.
bool generateRad(const RadcRun& run, RadcState& st, MainEvent& main, Event& vertex, Event& orig, IRandom& rng,
                 const CompleteEvFn& complete_ev, std::ostream* log = nullptr);

}  // namespace radiation
}  // namespace simc

#endif  // SIMC_RADIATION_RADC_H
