/**
 * @file
 * @ingroup physics
 * @brief Semi-inclusive (e,e'h)X cross section, a bit-for-bit port of semi_physics.f and semi_dilution.f.
 *
 * Semi-inclusive (e,e'h)X cross section of the Fortran SIMC, a bit-for-bit port of semi_physics.f and semi_dilution.f
 * (task P12.1):
 *   peepiX        semi_physics.f:1-631    d5sigma/dE_e' dOmega_e' dE_h dOmega_h [ub/MeV^2/sr^2] (D. Gaskell 2003, P. Bosted
 *                                         2020-2021): CTEQ5M PDFs x fragmentation functions (pi+-/pi0: P. Bosted's 9/20/2021
 *                                         fit in the Accardi z' variable; K+-: fDSS NLO) x a Gaussian pt2 dependence
 *                                         x the inclusive F1/F2 of F1F2IN21 (Christy 2021) x the EMC (slacemcfit) and
 *                                         hadron attenuation (Rhad_global) factors; and, when the decay is not
 *                                         simulated, the survival probability of the hadron (COMMON decdist).
 *   slacemcfit    semi_physics.f:634-650  SLAC EMC ratio fit (A > 2), 1 otherwise
 *   rhadGlobal    semi_physics.f:652-680  hadron multiplicity ratio fit R_h(A, z), 1 for A = 1
 *   semiDilution  semi_dilution.f:1-270   dilution factor (Binnewies FFs). DEAD CODE in SIMC: its only call,
 *                                         event.f:1524, is commented out (ntup%dilu = 1.0, not in the ntuple). Ported
 *                                         for completeness; no plug-in calls it.
 *
 * event.f:complete_main (l.1520-1522) calls peepiX(vertex, vertex0, main, survivalprob, .FALSE.) and sets
 * main%sigcc_recon = 1; the "central" variant (doing_cent = .TRUE., l.1523) is commented out, so the doing_cent branch
 * (sigc_flag, sigc_nbin, sigc_kin_min/max/ind) is dead code too; it is ported (peepiX's doing_cent argument).
 *
 * Numerics: semi_physics.f is compiled with -fdefault-real-8 and has no D-exponent literal: every operation is REAL(8)
 * (semi_physics.o converts to binary128 only for the -ff2c function-return convention, which is exact). x**2 -> x*x,
 * x**3 -> x*(x*x), the x**2..x**8 chain of slacemcfit as gcc's powi (x5 = x2*x3, x6 = x3*x3, x7 = x3*x4, x8 = x4*x4),
 * real exponents -> libm pow; operation order of the source (D-024). Units: MeV, MeV^2, rad as the Fortran event; the
 * cross section in ub/MeV^2/sr^2 (semi_physics.f:565, the jacobian 2 p_h^2 cos(theta_pq)/nu converts dz dpt2 to
 * dE_h dOmega_h).
 */
#ifndef SIMC_PHYSICS_SIDIS_SIDISCROSSSECTION_H
#define SIMC_PHYSICS_SIDIS_SIDISCROSSSECTION_H

#include <ostream>

#include "simc/core/Event.h"
#include "simc/physics/sidis/Cteq5Pdf.h"
#include "simc/physics/sidis/Fdss.h"

namespace simc {
namespace sidis {

/// @addtogroup physics
/// @{

/// Fortran: semi_physics.f:slacemcfit(A, x) (l.634-650)
/// @param A target A.
/// @param x Bjorken x.
/// @return the EMC ratio of the SLAC fit.
double slacemcfit(double A, double x);
/// Fortran: semi_physics.f:Rhad_global(A, z) (l.652-680)
/// @param A target A.
/// @param z hadron energy fraction z.
/// @return the hadron attenuation ratio.
double rhadGlobal(double A, double z);

/// The run constants that peepiX reads from COMMON blocks, and the D-031 switches (false = the Fortran, fortran_compat).
struct PeepiXConfig {
    double targA{0}, targZ{0};        ///< targ%A, targ%Z
    double Mtar_struck{0};            ///< targ%Mtar_struck [MeV] (Mp for doing_semi, dbase.f:301-302)
    bool doing_semipi{false}, doing_semika{false}, doing_pizero{false}, doing_hplus{false};  ///< SIDIS flags of the input.
    bool do_fermi{false};             ///< after dbase.f:202-206 (forced off for hydrogen)
    bool doing_decay{false};  ///< doing_decay of the input.
    int hadron_arm{1};                ///< spectrometer code (1 HMS, 2 SOS, 3 HRSR, 4 HRSL, 5 SHMS, ...)
    double ctau{0};                   ///< hadron c*tau [cm] (input)
    double spec_p_P{0};               ///< spec%p%P [MeV/c]
    double Mh2{0};                    ///< Mh2 [MeV^2]
    int sigc_flag{0}, sigc_nbin{0};   ///< central cross section binning (doing_cent only, dead code)
    double sigc_kin_min{0}, sigc_kin_max{0}, sigc_kin_ind{0};  ///< doing_cent binning (dead code).
    /// FORTRAN-BUG(D-051 F1): dsigdz = 0 when sum_sq <= 0 (all PDFs 0: x_Bj clamped to 1, or Q < Alambda); the
    /// Fortran divides 0/0 (semi_physics.f:508-513, the old guard l.378-400 is commented out) and returns NaN.
    bool sum_sq_guard{false};
    /// FORTRAN-BUG(D-051 F2): a non-positive Fermi-motion x_Bj (2 P.q <= 0, semi_physics.f:247-248) returns 0 (early
    /// return as below the threshold) instead of the Ctq5Pdf STOP (Ctq5Pdf.f:80-83) that ends the run.
    bool nonpositive_xbj_zero{false};
    /// FORTRAN-BUG(D-051 F4): pi0 fragmentation function of dbar = the pi+/pi- average (the Fortran averages d1 twice
    /// and leaves db = the pi+ or pi- value, semi_physics.f:490-495) and mhad = Mpi0 for doing_pizero (l.227 uses the
    /// charged Mpi while Mpi_gev, l.149, and dbase.f Mh use Mpi0).
    bool pizero_fixed{false};
    /// FORTRAN-BUG(D-051 F6, as D-043 F5): survival probability with beta*gamma of p_hadron = spec%p%P(1+delta/100).
    bool survival_p_hadron{false};
    /// FORTRAN-BUG(D-051 F6, as D-043 F6): zaero = 0 for the arms without an aerogel distance (SHMS) and the slope
    /// correction zaero*sqrt(1+dx^2+dy^2).
    bool survival_path_fixed{false};
};

/// SAVEd locals of peepiX (-fno-automatic): the first-call flags, zaero (only assigned for hadron_arm 1-4) and the
/// doing_cent kinematics pt2/zhad (kept from the previous call when sigc_flag is neither 0 nor 1).
struct PeepiXState {
    bool first{true}, firstqe{true}, first_cent{true};  ///< first-call flags.
    double zaero{0};  ///< aerogel distance [cm].
    double pt2{0}, zhad{0};  ///< doing_cent kinematics.
};

/// Intermediate values of one peepiX call (optional diagnostics for the L1 tables; never changes the numerics).
struct PeepiXTrace {
    bool threshold_return{false};     ///< mmpi2 below threshold: returned 0 at l.280-284
    double Mpi_gev{0}, mhad{0}, mtar{0};  ///< peepiX locals (semi_physics.f names; GeV).
    double xbj_raw{0};                ///< x_Bj before the clamp to 1
    double xbj{0}, Q2gev{0}, Qgev{0}, pt2gev{0}, wsq{0}, w{0}, mmpi2{0}, zhad{0}, pt2{0}, Ehad{0}, phad{0}, cthpq{0};  ///< peepiX locals: kinematics.
    double u{0}, ubar{0}, d{0}, dbar{0}, s{0}, sbar{0}, sum_sq{0};  ///< parton distributions of the proton and their sum of squares.
    double uA{0}, ubarA{0}, dA{0}, dbarA{0}, sA{0}, sbarA{0};  ///< parton distributions of the target.
    double Qgev_after{0};             ///< Qgev after the six Ctq5Pdf calls (Alambda if Q < Alambda)
    double sin2th2{0};  ///< sin^2(theta/2).
    double u1{0}, ub{0}, d1{0}, db{0}, s1{0}, sb{0};   ///< fragmentation functions (z D) used
    double xp{0}, zp{0}, sv{0}, yf{0}, yu{0};  ///< peepiX locals (semi_physics.f names).
    double dsigdz{0}, b{0}, rhad{0}, sighad{0};  ///< peepiX locals: cross-section factors.
    double F1p{0}, F2p{0}, F1n{0}, F2n{0}, emc{0}, F1{0}, F2{0}, W1{0}, W2{0}, sige{0}, sigsemi{0};  ///< peepiX locals: structure functions and inclusive cross section.
    double jacobian{0}, fac{0}, sigma{0};  ///< peepiX locals: Jacobian, factor and result.
    double pathlen{0}, betak{0}, gammak{0}, survival{0};  ///< peepiX locals: decay path and survival probability.
};

/// Fortran: semi_physics.f:peepiX(vertex, vertex0, main, survivalprob, doing_cent) (l.1-631).
/// Reads vertex (nu, q, uq, Q2, Ein, e%E, ue, pt2, zhad, theta_pq), vertex0 and main%theta_pq (doing_cent only),
/// main.pfermi (COMMON /pfermi/, do_fermi), main%FP%p (path, dx, dy) and main%SP%p%delta (!doing_decay).
/// Writes main%davejac (or main%johnjac for doing_cent), ntup%sigcm = sighad, ntup%xfermi (do_fermi, !doing_cent),
/// and when !doing_decay the survival probability into `survivalprob` and the COMMON decdist (main.ntup.decdist).
/// The early return below the mmpi2 threshold (l.280-284) leaves all of these unchanged.
/// `pdf` must be loaded (SetCtq5(1)); `fdss` (fDSS kaon NLO grid) is needed for doing_semika. Messages (xbj > 1, PDF
/// and FF range warnings, the doing_cent banner) go to `log`. Neither doing_semipi nor doing_semika: FortranStop
/// (l.231-232). Compat: a negative x_Bj STOPs in Ctq5Pdf (FortranStop, D-051 F2).
/// @param vertex the vertex.
/// @param vertex0 the central vertex (doing_cent only).
/// @param main the main event record (written).
/// @param[in,out] survivalprob the survival probability.
/// @param doing_cent central cross section mode.
/// @param cfg the run constants.
/// @param pdf the CTEQ5 set.
/// @param fdss the fDSS kaon grid (not owned), or nullptr for pions.
/// @param st the SAVEd state.
/// @param log stream for the messages (not owned), or nullptr.
/// @param[out] trace intermediate values (optional, may be nullptr).
/// @return the cross section.
double peepiX(const Event& vertex, const Event& vertex0, MainEvent& main, double& survivalprob, bool doing_cent,
              const PeepiXConfig& cfg, Cteq5Pdf& pdf, const FdssGrid* fdss, PeepiXState& st,
              std::ostream* log = nullptr, PeepiXTrace* trace = nullptr);

/// SAVEd locals of semi_dilution.
struct SemiDilutionState {
    bool first{true}, kfirst{true};  ///< first-call flags.
    double N{0}, a1{0}, a2{0};   ///< kept when neither doing_semipi nor doing_semika
};
/// Fortran: semi_dilution.f:semi_dilution(vertex, main) (l.1-270). DEAD CODE in SIMC (event.f:1524 commented out).
/// Reads vertex (Ein, Q2, nu, pt2, zhad, e%E), doing_semipi/semika/hplus/deutsemi. Returns the dilution factor
/// (999 if a proton or neutron cross section is not positive).
/// @param vertex the vertex.
/// @param doing_semipi,doing_semika,doing_hplus,doing_deutsemi SIDIS flags.
/// @param pdf the CTEQ5 set.
/// @param st the SAVEd state.
/// @param log stream for the messages (not owned), or nullptr.
/// @return the dilution factor.
double semiDilution(const Event& vertex, bool doing_semipi, bool doing_semika, bool doing_hplus, bool doing_deutsemi,
                    Cteq5Pdf& pdf, SemiDilutionState& st, std::ostream* log = nullptr);

/// @}

}  // namespace sidis
}  // namespace simc

#endif  // SIMC_PHYSICS_SIDIS_SIDISCROSSSECTION_H
