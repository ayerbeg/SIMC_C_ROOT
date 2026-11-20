/**
 * @file
 * @ingroup core
 * @brief simc::RunSetup: every constant of a run derived once, before the event loop, from the input parameters.
 *
 * The set-up converts the raw input (simc::SimcParameters, in `.inp` units) into the internal units and derives the
 * reaction flags, the spectrometer settings, the target, the generation limits and acceptance edges, the radiative
 * set-up, the spectral function and the normalisation constants. It is filled by simc::BuildRunSetup
 * (setup/SetupBuilder.h); RunSetup::FromParameters fills only the direct conversions (parameters, flags, run name,
 * spectrometer angles, histogram axes). Every stage reads it (read-only) through simc::StageContext::setup.
 *
 * Units are the internal units after conversion: MeV, MeV/c, cm, rad, g/cm^2, delta in %. `params` keeps the raw
 * input untouched. The structs mirror the derived types and COMMON blocks of the original program; each field names
 * its origin.
 * @note Origin: Fortran SIMC dbase.f:dbase_read (post-read logic), init.f:min_max_init / target_init / limits_init /
 *       radc_init / theory_init, simc.f:95-101; types from structures.inc, target.inc, simulate.inc, radc.inc,
 *       brem.inc.
 */
#ifndef SIMC_CORE_RUNSETUP_H
#define SIMC_CORE_RUNSETUP_H

#include <array>
#include <string>
#include <vector>

#include "simc/core/SimcParameters.h"

namespace simc {

/** @addtogroup core
 *  @{ */

/**
 * @brief Which reaction family (and target type) the run simulates.
 *
 * The family flags (doing_eep, doing_pion, doing_kaon, doing_delta, doing_semi, doing_rho, doing_phsp) come from the
 * input; the target-specific flags (hydrogen `hyd*`, deuterium `deut*`, heavier `he*`) are derived from the target A.
 * At most one family is set for a valid input.
 * @note Origin: Fortran SIMC dbase.f:dbase_read l.124-234 (the doing_* logicals).
 */
struct ReactionFlags {
    /// (e,e'p): any target, hydrogen elastic, deuterium, heavier nucleus.
    bool doing_eep{false}, doing_hyd_elast{false} /**< @copydoc doing_eep */, doing_deuterium{false} /**< @copydoc doing_eep */, doing_heavy{false} /**< @copydoc doing_eep */;
    /// Exclusive pion production: any target, on hydrogen, on deuterium, on a heavier nucleus.
    bool doing_pion{false}, doing_hydpi{false} /**< @copydoc doing_pion */, doing_deutpi{false} /**< @copydoc doing_pion */, doing_hepi{false} /**< @copydoc doing_pion */;
    /// Exclusive kaon production: any target, on hydrogen, on deuterium, on a heavier nucleus.
    bool doing_kaon{false}, doing_hydkaon{false} /**< @copydoc doing_kaon */, doing_deutkaon{false} /**< @copydoc doing_kaon */, doing_hekaon{false} /**< @copydoc doing_kaon */;
    /// Delta production: any target, on hydrogen, on deuterium, on a heavier nucleus.
    bool doing_delta{false}, doing_hyddelta{false} /**< @copydoc doing_delta */, doing_deutdelta{false} /**< @copydoc doing_delta */, doing_hedelta{false} /**< @copydoc doing_delta */;
    /// Semi-inclusive (SIDIS) production: any hadron, pion, kaon.
    bool doing_semi{false}, doing_semipi{false} /**< @copydoc doing_semi */, doing_semika{false} /**< @copydoc doing_semi */;
    /// Semi-inclusive production on hydrogen, deuterium, a heavier nucleus.
    bool doing_hydsemi{false}, doing_deutsemi{false} /**< @copydoc doing_hydsemi */, doing_hesemi{false} /**< @copydoc doing_hydsemi */;
    /// Exclusive rho production: any target, on hydrogen, on deuterium, on a heavier nucleus.
    bool doing_rho{false}, doing_hydrho{false} /**< @copydoc doing_rho */, doing_deutrho{false} /**< @copydoc doing_rho */, doing_herho{false} /**< @copydoc doing_rho */;
    /// Phase-space (acceptance) mode; neutral-pion detection; positive hadron; Fermi motion of the struck nucleon.
    bool doing_phsp{false}, doing_pizero{false} /**< @copydoc doing_phsp */, doing_hplus{false} /**< @copydoc doing_phsp */, do_fermi{false} /**< @copydoc doing_phsp */;
    /// Polarised-target magnetic field tracking is on.
    bool using_tgt_field{false};
    /// A(e,e')A nuclear elastic scattering (upstream SIMC commit 4df7351, dbase.f:216-218): the raw input flag. With it
    /// doing_eep is true and none of doing_hyd_elast / doing_deuterium / doing_heavy is set.
    bool doing_nuc_elast{false};
};

/**
 * @brief Derives the reaction flags from the input parameters.
 * @param p The input parameters.
 * @param problems If not null, receives one message per condition where the original program stops (inconsistent
 *        family / target combinations); the flags are still returned.
 * @return The derived flags.
 * @note Origin: Fortran SIMC dbase.f l.124-234 (the STOP checks become messages in @p problems).
 */
ReactionFlags DeriveReactionFlags(const SimcParameters& p, std::vector<std::string>* problems = nullptr);

/**
 * @brief The reaction title line printed in the `.hist` report (e.g. "H(e,e'p)").
 * @param f The reaction flags.
 * @param p The input parameters (hadron species, target).
 * @return The title text.
 * @note Origin: Fortran SIMC simc.f:report l.664-790.
 */
std::string ReactionTitle(const ReactionFlags& f, const SimcParameters& p);

/**
 * @brief One axis of the diagnostic 1-D histograms (50 bins starting at min, of width bin).
 * @note Origin: Fortran SIMC histograms.inc (min, bin; nHbins = 50).
 */
struct HistAxis {
    double min{0};   ///< Lower edge of the first bin (unit of the histogrammed quantity).
    double bin{1};   ///< Bin width (same unit).
};

/**
 * @brief Axes of the generated / generated-inside / reconstructed diagnostic histograms (all three share them).
 * @note Origin: Fortran SIMC H%gen, H%geni, H%RECON; init.f:limits_init l.515-570.
 */
struct DiagHistAxes {
    /// Axes of one arm.
    struct Arm {
        HistAxis delta /**< @copydoc xptar */, yptar /**< @copydoc xptar */, xptar;   ///< delta [%], yptar and xptar [rad].
    };
    /// Electron arm (e) and hadron arm (p).
    Arm e, p /**< @copydoc e */;
    /// Missing energy [MeV] and missing momentum [MeV/c].
    HistAxis Em, Pm /**< @copydoc Em */;
};

// ---------------------------------------------------------------------------------------------
// Derived types (structures.inc / target.inc). Defaults = init.f:min_max_init (l.921-1193):
// every min/max pair that min_max_init touches starts at -/+1e10; everything else is zero
// (zero-initialised COMMON).
// ---------------------------------------------------------------------------------------------

constexpr double kMinMaxInit = 1.0e10;   ///< Fortran: min_max_init xmax (init.f:936-937)

/**
 * @brief A [min, max] interval (a cut, a limit or an acceptance edge), unbounded (-/+1e10) by default.
 * @note Origin: Fortran SIMC structures.inc cutstype, initialised by init.f:min_max_init.
 */
struct Cut {
    double min{-kMinMaxInit};   ///< Lower bound (unit of the quantity).
    double max{kMinMaxInit};    ///< Upper bound (unit of the quantity).
};

/**
 * @brief A "slop": the extra range added to a generation limit to cover resolution and radiation effects.
 * @note Origin: Fortran SIMC structures.inc slop_item (lo, hi, used).
 */
struct SlopItem {
    /// lo, hi: run accumulators of the smallest / largest value seen (0 at start); used: the slop applied.
    double lo{0}, hi{0} /**< @copydoc lo */, used{0} /**< @copydoc lo */;
};

/**
 * @brief Average and extreme values of a target quantity (energy loss, thickness, Coulomb shift).
 * @note Origin: Fortran SIMC target.inc target_function.
 */
struct TargetFunction {
    /// Average, minimum and maximum over the target (unit of the quantity).
    double ave{0}, min{0} /**< @copydoc ave */, max{0} /**< @copydoc ave */;
};

/// Fortran: target.inc type target_info (COMMON /target_info/ targ), after dbase.f conversions.
struct TargetInfo {
    double A{0} /**< @copydoc N */, Z{0} /**< @copydoc N */, N{0};          ///< dbase.f:250 N = A - Z
    double mass_amu{0};               ///< [amu] input
    double M{0};                      ///< [MeV] dbase.f:251 (+ forced Mp/Md l.256-266)
    double mrec_amu{0};               ///< [amu] input
    double Mrec{0};                   ///< [MeV] dbase.f:252 (+ l.256-272, 410-412, 437-439, 445-450)
    double rho{0};                    ///< [g/cm^3]
    double thick{0};                  ///< [g/cm^2] dbase.f:452 thick/1000
    double angle{0};                  ///< [rad] dbase.f:454 (forced 0 for liquid, l.458-461)
    double abundancy{0};              ///< [%]
    double length{0};                 ///< [cm] dbase.f:453 thick/rho
    double zoffset{0} /**< @copydoc yoffset */, xoffset{0} /**< @copydoc yoffset */, yoffset{0};   ///< [cm]
    double X0{0};                     ///< [g/cm^2] init.f:target_init l.33-42
    double X0_cm{0};                  ///< [cm] init.f:43
    double L1{0} /**< @copydoc L2 */, L2{0};              ///< init.f:19-24
    double fr1{0} /**< @copydoc fr2 */, fr2{0};            ///< [cm] raster
    /// Raster pattern (1 bedpost, 2 circle, 3 flat square) and target can type (1 beer can, 2 pudding can, 3 cryo17).
    int fr_pattern{0}, can{0} /**< @copydoc fr_pattern */;
    std::array<TargetFunction, 3> Eloss{};   ///< [MeV] 1 beam, 2 e', 3 hadron. ave: target_init (hook),
                                             ///< min/max: extreme_trip_thru_target (hook); zeroed if !using_Eloss
    std::array<TargetFunction, 3> teff{};    ///< [rad. len.] same origin (never zeroed)
    TargetFunction Coulomb;                  ///< [MeV] init.f:target_init l.62-82
    double Coulomb_constant{0};              ///< [MeV] init.f:75/79
    std::array<double, 3> musc_max{};        ///< [rad] extreme_target_musc (hook); zeroed if !mc_smear (init.f:254-258)
    double musc_nsig_max{0};                 ///< extreme_target_musc (hook)
    double Mtar_struck{0} /**< @copydoc Mrec_struck */, Mrec_struck{0};   ///< [MeV] dbase.f:277-450
};

/// Fortran: structures.inc spectrometer offsets (after dbase.f:472-475 xptar/yptar mrad -> rad).
struct SpecOffset {
    double x{0} /**< @copydoc z */, y{0} /**< @copydoc z */, z{0};          ///< [cm]
    double xptar{0} /**< @copydoc yptar */, yptar{0};        ///< [rad] (NB: the .hist prints these under 'mr' without *1000)
};

/// Fortran: structures.inc spec%e / spec%p (after dbase.f:238-262).
struct SpecArm {
    double P{0};                      ///< [MeV/c]
    double theta{0};                  ///< [rad] dbase.f:239/242 abs(theta)/degrad
    /// cos(theta) and sin(theta) of the central angle.
    double cos_th{0}, sin_th{0} /**< @copydoc cos_th */;
    double phi{0};                    ///< [rad] dbase.f:247-262: 3pi/2 (beam-right arms, SpectrometerId::IsBeamRight) or pi/2
    SpecOffset offset;                ///< Mispointing offsets of the arm.
};

/// Fortran: arm cuts in double_arm_cuts (SPedge%e / SPedge%p): delta [%], yptar/xptar [rad], z.
struct ArmCuts {
    /// delta [%], yptar, xptar [rad], z [cm] acceptance edges.
    Cut delta, yptar /**< @copydoc delta */, xptar /**< @copydoc delta */, z /**< @copydoc delta */;
};

/// Fortran: edge_arm (E, yptar, xptar).
struct EdgeArm {
    /// Energy [MeV] and target slopes yptar, xptar [rad].
    Cut E, yptar /**< @copydoc E */, xptar /**< @copydoc E */;
};

/// Fortran: edge_true (edge, VERTEXedge).
struct EdgeTrue {
    /// Electron arm (e) and hadron arm (p).
    EdgeArm e, p /**< @copydoc e */;
    /// Missing energy [MeV], missing momentum [MeV/c], recoil mass [MeV], recoil kinetic energy [MeV] (of the system
    /// and of the struck nucleon).
    Cut Em, Pm /**< @copydoc Em */, Mrec /**< @copydoc Em */, Trec /**< @copydoc Em */, Trec_struck /**< @copydoc Em */;
};

/// Fortran: arm_limits (gen%e / gen%p).
struct GenArm {
    /// delta [%], yptar, xptar [rad], energy E [MeV].
    Cut delta, yptar /**< @copydoc delta */, xptar /**< @copydoc delta */, E /**< @copydoc delta */;
};

/// Fortran: gen_limits (gen).
struct GenLimits {
    /// Electron arm (e) and hadron arm (p).
    GenArm e, p /**< @copydoc e */;
    /// Sum of the generated energies [MeV] and recoil kinetic energy [MeV].
    Cut sumEgen, Trec /**< @copydoc sumEgen */;
    double xwid{0} /**< @copydoc ywid */, ywid{0};          ///< [cm] input gen%xwid/ywid
};

/// Fortran: sloptype (slop%total, slop%MC).
struct Slops {
    /// Slops of the missing energy [MeV] and missing momentum [MeV/c].
    struct Total {
        SlopItem Em /**< @copydoc Pm */, Pm;   ///< Missing-energy slop [MeV] and missing-momentum slop [MeV/c].
    } total;               ///< Slops on the missing quantities.
    /// Slops of one arm.
    struct Arm {
        SlopItem delta /**< @copydoc xptar */, yptar /**< @copydoc xptar */, xptar;   ///< delta [%], yptar and xptar [rad].
    };
    /// Monte Carlo (resolution) slops of the two arms.
    struct MC {
        Arm e /**< @copydoc p */, p;   ///< Electron arm (e) and hadron arm (p).
    } MC;           ///< Monte Carlo (resolution) slops per arm.
};

/// Run-level logicals and codes after dbase.f/init.f adjustments (simulate.inc, radc.inc, brem.inc).
/// Where the Fortran overwrites an input variable, the effective value is here; params keeps the raw one.
struct RunSwitches {
    bool using_E_arm_montecarlo{false};  ///< dbase.f:533 spect_mode /= 1 and /= -1
    bool using_P_arm_montecarlo{false};  ///< dbase.f:534 spect_mode /= 1 and /= -2
    bool using_Eloss{false};             ///< Energy loss in the target is simulated (input using_Eloss).
    bool correct_Eloss{false};           ///< dbase.f:542-545 forced .false. if !using_Eloss
    bool using_Coulomb{false};           ///< dbase.f:547 forced .false. for nint(Z)=1
    bool mc_smear{false};                ///< Multiple scattering and detector smearing are simulated (input mc_smear).
    bool using_rad{false};               ///< Radiative corrections are on (input using_rad).
    int rad_flag{0};                     ///< dbase.f:228 forced 0 by doing_phsp
    int extrad_flag{0};                  ///< init.f:radc_init l.619-627 default 3 (rad_flag 0) or 1 (1..3)
    std::array<bool, 3> doing_tail{};    ///< dbase.f:500-519
    bool hardwired_rad{false};           ///< dbase.f:521-525 Egamma_gen_max > 0.01
    // brem.inc flags set by radc_init (init.f:634-645)
    bool produce_output{false};          ///< = debug(1)
    bool exponentiate{false};            ///< use_expon == 1
    bool include_hard{true};             ///< Include the hard-photon correction in the radiative weight (always true).
    bool calculate_spence{true};         ///< Evaluate the Spence function in the radiative weight (always true).
};

/// Fortran: theory_init (init.f:828-920) — spectral-function table for D(e,e'p) / A(e,e'p).
struct TheoryData {
    /// Missing-momentum axis of one distribution (structures.inc axis).
    struct PmAxis {
        /// Lower edge, upper edge and bin width [MeV/c].
        double min{0}, max{0} /**< @copydoc min */, bin{0} /**< @copydoc min */;
        int n{0};   ///< Number of bins.
    };
    std::string file;                   ///< theory_file (h2.theory, c12.theory, ...)
    bool loaded{false};                 ///< True when the file was read.
    int nrhoPm{0};                      ///< Number of momentum distributions (shells) in the file.
    double E_Fermi{0};                  ///< Fermi energy of the nucleus [MeV] (from the file).
    std::vector<double> nprot /**< @copydoc Em_int */, Em /**< @copydoc Em_int */, Emsig /**< @copydoc Em_int */, bs_norm /**< @copydoc Em_int */, Em_int;   ///< per distribution m (0-based)
    std::vector<PmAxis> Pm;                                  ///< Pm_theory(m)
    std::vector<std::vector<double>> theory;                 ///< theory(m, n) (0-based both)
};

/// The spectral function of D(e,e'p) / A(e,e'p) (task P8.1): which kind and which file, as chosen at start-up
/// (dbase.f:555-558 theory_init, dbase.f:592-621 sf_lookup_init), possibly named explicitly in the JSON
/// (pipeline.reaction.spectral_function, SetupOptions::sf_file).
struct SpectralFunctionChoice {
    std::string type;        ///< "theory" (init.f:theory_init), "benhar" (sf_lookup_init), "" = none (H, phsp, ...)
    std::string file;        ///< file name, looked up in RunSetup::data_dir
    bool fallback{false};    ///< D-033 F1: the Fortran default file for an A without its own file (C-12)
    bool explicit_file{false};   ///< the file was named in the input (no fallback warning)
    bool proton{true};       ///< Benhar: proton column SFp (Mtar_struck = Mp) or neutron SFn (dbase.f:610-620)
};

/// Fortran: dbase.f:571-596 momentum distribution (pval, mprob, nump) for deut/he pion/kaon/semi.
struct MomentumDistribution {
    std::string file;                   ///< File name of the distribution (in RunSetup::data_dir).
    std::vector<double> pval /**< @copydoc mprob */, mprob;    ///< mprob normalised to mprob(nump)
};

/// Fortran: radc.inc quantities fixed by radc_init / limits_init (per-event ones live elsewhere).
struct RadSetup {
    double etatzai{0};                  ///< init.f:630 (12+(Z+1)/(Z*L1+L2))/9
    double Egamma_tot_max{0};           ///< init.f:limits_init l.408-418
    double Egamma1_max{0} /**< @copydoc Egamma3_max */, Egamma2_max{0} /**< @copydoc Egamma3_max */, Egamma3_max{0};   ///< init.f:419-421 (0 if tail off)
};

/// Layout of the diagnostic ntuple (tree "h10"), task U.2 (D-082).
/// Upstream: the layout of upstream SIMC 59ccdf6 (NtupleInit.f / results_write.f after e4dd175), the DEFAULT;
/// Legacy: the layout of SIMC 294f662 (and 20d46fc for doing_nuc_elast). CLI --tree-layout, JSON "output": {"layout"}.
enum class TreeLayout { Upstream /**< upstream 59ccdf6 (default) */, Legacy /**< 294f662 */ };

/**
 * @brief Every constant of a run, derived from the input before the event loop (see the file description).
 *
 * Built once by simc::BuildRunSetup and then read-only: shared by all stages and, in the parallel event loop, by all
 * workers (concurrent reads are safe). Copyable value type.
 */
struct RunSetup {
    SimcParameters params;          ///< raw input (.inp units)
    ReactionFlags flags;            ///< Reaction family and target type.
    std::string run_name{"simc"};   ///< Fortran `base` (file stem)
    std::string output_dir{"."};    ///< Fortran 'outfiles/' and 'worksim/'
    /// Run option `fortran_compat` (D-031, task P7.5). false (default) = the clear Fortran bugs listed in
    /// docs/translation/markdown/FORTRAN_DIFFERENCES.md are CORRECTED (code marked `FORTRAN-BUG(D-xxx)`); true = the exact Fortran
    /// behaviour, used by the bitwise L1/replay/regression tests. Set once from SetupOptions::fortran_compat by
    /// BuildRunSetup and read by every stage through StageContext::setup (no global).
    bool fortran_compat{false};
    /// Ntuple layout of the tree h10 (task U.2): upstream 59ccdf6 (default) or legacy 294f662. Set from
    /// SetupOptions::tree_layout by BuildRunSetup.
    TreeLayout tree_layout{TreeLayout::Upstream};

    // dbase.f:239-244 (kept for the P1.4 sinks; identical to spec.e/p.theta, sin_th, cos_th)
    double e_theta{0} /**< @copydoc e_cos_th */, e_sin_th{0} /**< @copydoc e_cos_th */, e_cos_th{0};   ///< rad
    double p_theta{0} /**< @copydoc p_cos_th */, p_sin_th{0} /**< @copydoc p_cos_th */, p_cos_th{0};   ///< rad (hadron arm)

    DiagHistAxes hist;              ///< init.f:limits_init l.515-570

    // ---------------- filled by BuildRunSetup (P2.1) ----------------
    double Mh{0} /**< @copydoc Mh2_final */, Mh2{0} /**< @copydoc Mh2_final */, Mh2_final{0};   ///< [MeV] dbase.f:136-219, 221-222
    double sign_hadron{0};                ///< dbase.f:277-446 (0 if no branch sets it, e.g. doing_phsp)
    double Ebeam{0};                      ///< [MeV]
    double dEbeam{0};                     ///< [MeV] dbase.f:238 Ebeam*dEbeam/100 (full width)
    double Ebeam_vertex_ave{0};           ///< [MeV] init.f:283 Ebeam + Coulomb.ave - Eloss(1).ave
    struct { SpecArm e, p; } spec;        ///< dbase.f:238-262, 472-475
    TargetInfo targ;                      ///< The target after the dbase.f/init.f conversions.
    double targ_Bangle{0} /**< @copydoc targ_Bphi */, targ_Bphi{0};  ///< [rad] dbase.f:455-456
    struct { ArmCuts e, p; } SPedge;      ///< dbase.f:479-490 then + slop%MC (init.f:209-220)
    struct { Cut Em, Pm; } cuts;          ///< [MeV] dbase.f:536-540 (±1e6 for pion/kaon/delta or Em.min=Em.max)
    EdgeTrue edge;                        ///< init.f:limits_init
    EdgeTrue VERTEXedge;                  ///< init.f:limits_init
    GenLimits gen;                        ///< init.f:limits_init (+ dbase.f:524-531 /1000 of the ±1e10 defaults)
    Slops slop;                           ///< init.f:limits_init
    RunSwitches sw;                       ///< Effective run switches (after the set-up adjustments).
    RadSetup rad;                         ///< Run constants of the radiative corrections.
    TheoryData theory;                    ///< only for doing_deuterium / doing_heavy (no Benhar SF)
    SpectralFunctionChoice sf;            ///< D(e,e'p) / A(e,e'p) spectral function (type, file)
    std::string data_dir;                 ///< SetupOptions::data_dir: the Fortran working-directory data files
    MomentumDistribution pdist;           ///< only for deut/he pion, kaon, semi
    double targetfac{0};                  ///< simc.f:94-95
    double luminosity{0};                 ///< [ub^-1] simc.f:101 EXPER%charge/targetfac

    /// True when filled by BuildRunSetup (false for RunSetup::FromParameters only).
    bool built{false};
    /// Provenance of the energy-loss/multiple-scattering numbers (see setup/TargetLossProvider.h).
    bool target_losses_from_provider{false};
    /// Fortran STOP conditions reached (the Fortran aborts: the fields after it are not meaningful).
    std::vector<std::string> stops;
    /// Informational messages the Fortran prints (WARNING/NOTE lines of dbase.f/init.f).
    std::vector<std::string> messages;
    /// Pieces that are not ported yet and were skipped (e.g. sf_lookup_init).
    std::vector<std::string> not_ported;
    /**
     * @brief Tells whether the set-up is usable.
     * @return True if no stop condition was reached (#stops is empty).
     */
    bool ok() const { return stops.empty(); }

    /**
     * @brief A partial set-up: the parameters, the reaction flags and the spectrometer angles (direct conversions).
     *
     * For tools and sinks that need only those; a run uses simc::BuildRunSetup.
     * @param p The input parameters (copied into #params).
     * @param run_name The run name (file stem of the outputs).
     * @return The set-up with #built false.
     */
    static RunSetup FromParameters(const SimcParameters& p, const std::string& run_name = "simc");
};

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_RUNSETUP_H
