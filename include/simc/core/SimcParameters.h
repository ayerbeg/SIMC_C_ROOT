/**
 * @file
 * @ingroup core
 * @brief simc::SimcParameters: every input variable of a SIMC `.inp` file, and the parameter table that drives the
 *        `.inp` reader, the JSON (de)serialisation and the reference tables.
 *
 * The struct has exactly one field for each of the 106 input variables of the original program (107 since upstream 4df7351 added doing_nuc_elast), grouped like its
 * COMMON blocks / derived types.
 *
 * Conventions:
 *  - Values are stored EXACTLY as written in the `.inp`, in the `.inp` units (MeV, MeV/c, deg, mrad, mg/cm^2, %, ...).
 *    No conversion is applied here: that is the job of simc::RunSetup. simc::EffectiveValues (input/EffectiveValues.h)
 *    shows the converted values.
 *  - Defaults are the TRUE defaults of the original program when a variable is absent from the `.inp`: zero / false /
 *    blank, except SPedge.*.*.min/max = -/+1e10 and cuts.Em.min/max = -/+1e10 (set before the input is read), and
 *    use_benhar_sf = false, random_state_file = blank, random_seed = 0, extra_dbase_file = blank.
 *  - The `*_int` integer flags of the original that are converted to logicals (int > 0 means true) are stored as
 *    `bool`. Integer codes that are tested numerically (use_expon, Nntu, rad_flag, ...) stay `int`.
 *  - The JSON path of each field is its original name split on '%' (e.g. "targ%thick" -> {"targ":{"thick":...}});
 *    see simc::ParameterTable().
 *
 * @note Origin: Fortran SIMC dbase.f:regallvars (l.967-1128, 106 CTP variables; + doing_nuc_elast of upstream 4df7351) and dbase.f:convert_to_logical
 *       (l.1129-1166). The third argument of the CTP regparm*() calls is a title pointer, NOT a default (e.g.
 *       doing_hplus, sigc_nbin and sigc_kin_max default to 0); source: the reference input-variable table (D-005).
 */
#ifndef SIMC_CORE_SIMCPARAMETERS_H
#define SIMC_CORE_SIMCPARAMETERS_H

#include <array>
#include <string>
#include <vector>

#include "simc/core/Json.h"

namespace simc {

/** @addtogroup core
 *  @{ */

/// Fortran: cutstype / rangetype (structures.inc) — a min/max pair.
struct ParamRange {
    double min{0.0};   ///< Lower bound (unit of the parameter).
    double max{0.0};   ///< Upper bound (unit of the parameter).
};

/**
 * @brief Every input variable of a SIMC run, in `.inp` units (see the file description for the conventions).
 *
 * Filled by simc::InpReader (from `.inp` files) or simc::FromJson (from JSON); converted into the run constants by
 * simc::BuildRunSetup. Field documentation: original name, [unit or kind], meaning. Plain copyable value type.
 */
struct SimcParameters {
    // ------------------------------------------------------------------
    // dbase.f:984-992 ("RESTSW" section of regallvars)
    // ------------------------------------------------------------------
    bool mc_smear{false};            ///< Fortran: mc_smear (mc_smear_int) [flag] target/hut mult. scattering + DC smearing
    int electron_arm{0};             ///< Fortran: electron_arm [code] simc::SpectrometerId (core/SpectrometerId.h); 6 rejected (D-032)
    int hadron_arm{0};               ///< Fortran: hadron_arm [code] same codes as electron_arm
    int use_first_cer{0};            ///< Fortran: use_first_cer [flag] registered, only printed in .hist
    double transparency{0.0};        ///< Fortran: transparency [dimensionless] A(e,e'p) Benhar SF weight factor
    bool use_benhar_sf{false};       ///< Fortran: use_benhar_sf (use_benhar_sf_int) [flag]; dbase.f:69 sets .false.
    std::string extra_dbase_file;    ///< Fortran: extra_dbase_file [character*60 file stem] 'infiles/'//name(.inp), dbase.f:105-117
    std::string random_state_file;   ///< Fortran: random_state_file [character*80 path] restore RNG state (simc.f:106)
    int random_seed{0};              ///< Fortran: random_seed [int] -1 => time() (dbase.f:119)

    // ------------------------------------------------------------------
    // dbase.f:996-1010 (EXPERIMENT)
    // ------------------------------------------------------------------
    double Ebeam{0.0};               ///< Fortran: Ebeam [MeV]
    double dEbeam{0.0};              ///< Fortran: dEbeam [% of Ebeam, full width, flat]; dbase.f:238 -> MeV
    /// Experiment-level quantities (Fortran EXPER).
    struct Exper {
        double charge{0.0};          ///< Fortran: EXPER%charge [mC]
    } EXPER;                         ///< Fortran: EXPER
    bool doing_kaon{false};          ///< Fortran: doing_kaon (doing_kaon_int) [flag]
    int which_kaon{0};               ///< Fortran: which_kaon [code] 0/1/2 Lambda/Sigma0/Sigma-, 10/11/12 coherent
    bool doing_pion{false};          ///< Fortran: doing_pion (doing_pion_int) [flag]
    int which_pion{0};               ///< Fortran: which_pion [code] 0/1 pi+/pi-, 2/3 Delta final states, 10/11 coherent
    bool doing_delta{false};         ///< Fortran: doing_delta (doing_delta_int) [flag]
    bool doing_semi{false};          ///< Fortran: doing_semi (doing_semi_int) [flag]
    bool doing_hplus{false};         ///< Fortran: doing_hplus (doing_hplus_int) [flag] (title arg '1' is NOT a default)
    bool doing_pizero{false};        ///< Fortran: doing_pizero (doing_pizero_int) [flag]
    bool doing_nuc_elast{false};     ///< Fortran (upstream 4df7351): doing_nuc_elast (doing_nuc_elast_int) [flag] A(e,e')A nuclear elastic
    int pizero_ngamma{0};            ///< Fortran: pizero_ngamma [count] 1 or 2 photons required
    bool doing_rho{false};           ///< Fortran: doing_rho (doing_rho_int) [flag]
    bool doing_decay{false};         ///< Fortran: doing_decay (doing_decay_int) [flag]
    double ctau{0.0};                ///< Fortran: ctau [cm] hadron c*tau

    // ------------------------------------------------------------------
    // dbase.f:1014 (DEBUG)
    // ------------------------------------------------------------------
    std::array<bool, 6> debug{};     ///< Fortran: debug(1:6) (debug_int(6)) [flags]; CTP debug(1) = element 0

    // ------------------------------------------------------------------
    // dbase.f:1018-1035 (TARGET), target.inc
    // ------------------------------------------------------------------
    /// The target as written in the input (Fortran targ, target.inc).
    struct Target {
        double A{0.0};               ///< Fortran: targ%A [nucleons]
        double Z{0.0};               ///< Fortran: targ%Z [charge]
        double mass_amu{0.0};        ///< Fortran: targ%mass_amu [amu]
        double mrec_amu{0.0};        ///< Fortran: targ%mrec_amu [amu]
        double rho{0.0};             ///< Fortran: targ%rho [g/cm^3]
        double thick{0.0};           ///< Fortran: targ%thick [mg/cm^2]; dbase.f:451 -> g/cm^2
        double xoffset{0.0};         ///< Fortran: targ%xoffset [cm] +x = beam right
        double yoffset{0.0};         ///< Fortran: targ%yoffset [cm] +y = up
        int fr_pattern{0};           ///< Fortran: targ%fr_pattern [code] 1=bedpost 2=circle 3=flat square
        double fr1{0.0};             ///< Fortran: targ%fr1 [cm]
        double fr2{0.0};             ///< Fortran: targ%fr2 [cm]
        double zoffset{0.0};         ///< Fortran: targ%zoffset [cm]
        double angle{0.0};           ///< Fortran: targ%angle [deg]; dbase.f:453 -> rad
        double abundancy{0.0};       ///< Fortran: targ%abundancy [%]
        int can{0};                  ///< Fortran: targ%can [code] 1=beer can 2=pudding can 3=cryo17
    } targ;                          ///< Fortran: targ
    double targ_Bangle{0.0};         ///< Fortran: targ_Bangle [deg] polarized target field angle
    double targ_Bphi{0.0};           ///< Fortran: targ_Bphi [deg]
    double targ_pol{0.0};            ///< Fortran: targ_pol [dimensionless]

    // ------------------------------------------------------------------
    // dbase.f:1039-1063 (E_ARM_MAIN, P_ARM_MAIN, *_OFFSET)
    // ------------------------------------------------------------------
    /// Mispointing offsets of a spectrometer arm.
    struct Offset {
        double x{0.0};               ///< Fortran: spec%?%offset%x [cm]
        double y{0.0};               ///< Fortran: spec%?%offset%y [cm]
        double z{0.0};               ///< Fortran: spec%?%offset%z [cm]
        double xptar{0.0};           ///< Fortran: spec%?%offset%xptar [mrad]; dbase.f:472-475 -> rad
        double yptar{0.0};           ///< Fortran: spec%?%offset%yptar [mrad]
    };
    /// Central setting of a spectrometer arm.
    struct Arm {
        double P{0.0};               ///< Fortran: spec%?%P [MeV/c] central momentum
        double theta{0.0};           ///< Fortran: spec%?%theta [deg]; dbase.f:239-244 abs()/degrad
        Offset offset;               ///< Fortran: spec%?%offset
    };
    /// The two spectrometer arms (Fortran spec).
    struct Spec {
        Arm e;                       ///< Fortran: spec%e
        Arm p;                       ///< Fortran: spec%p
    } spec;                          ///< Fortran: spec
    /// Beam spot size (Fortran gen%xwid, gen%ywid).
    struct Gen {
        double xwid{0.0};            ///< Fortran: gen%xwid [cm, 1 sigma]
        double ywid{0.0};            ///< Fortran: gen%ywid [cm, 1 sigma]
    } gen;                           ///< Fortran: gen (beam spot widths)

    // ------------------------------------------------------------------
    // dbase.f:1067-1069 (MISC2INT, radc.inc)
    // ------------------------------------------------------------------
    int use_expon{0};                ///< Fortran: use_expon [flag, tested ==0/1]
    int one_tail{0};                 ///< Fortran: one_tail [code] 0=all 1=e 2=e' 3=p -3=all but p ...
    int intcor_mode{0};              ///< Fortran: intcor_mode [code]

    // ------------------------------------------------------------------
    // dbase.f:1073-1104 (SIMULATE)
    // ------------------------------------------------------------------
    int ngen{0};                     ///< Fortran: ngen [events] >0 successes, <0 |ngen| tries
    bool hard_cuts{false};           ///< Fortran: hard_cuts (hard_cuts_int) [flag]
    bool using_rad{false};           ///< Fortran: using_rad (using_rad_int) [flag]
    int spect_mode{0};               ///< Fortran: spect_mode [code] 0=e+p -1=p only -2=e only 1=none
    bool doing_phsp{false};          ///< Fortran: doing_phsp (doing_phsp_int) [flag] (acceptance mode, D-008)
    /// Cuts on the reconstructed event (Fortran cuts).
    struct Cuts {
        ParamRange Em{-1.0e10, 1.0e10};  ///< Fortran: cuts%Em%min/max [MeV]; default init.f:min_max_init
    } cuts;                          ///< Fortran: cuts
    bool using_Eloss{false};         ///< Fortran: using_Eloss (using_Eloss_int) [flag]
    bool correct_Eloss{false};       ///< Fortran: correct_Eloss (correct_eloss_int) [flag]
    bool correct_raster{false};      ///< Fortran: correct_raster (correct_raster_int) [flag]
    bool using_HMScoll{false};       ///< Fortran: using_HMScoll (using_HMScoll_int) [flag]
    bool using_SHMScoll{false};      ///< Fortran: using_SHMScoll (using_SHMScoll_int) [flag]
    int deForest_flag{0};            ///< Fortran: deForest_flag [code] 0=sigcc1 1=sigcc2 -1=sigcc1 on-shell
    int rad_flag{0};                 ///< Fortran: rad_flag [code] (radc.inc)
    int extrad_flag{0};              ///< Fortran: extrad_flag [code] 0 => 3 or 1 (init.f:621-626)
    std::array<double, 3> lambda{};  ///< Fortran: lambda(1:3) [dimensionless] (overwritten at init.f:692-695, D-005)
    bool using_cit_generation{false};///< Fortran: using_cit_generation (using_cit_generation_int) [flag]
    int Nntu{0};                     ///< Fortran: Nntu [flag, tested >0] write ntuple
    bool using_Coulomb{false};       ///< Fortran: using_Coulomb (using_Coulomb_int) [flag]; forced off for Z=1 (dbase.f:547)
    double dE_edge_test{0.0};        ///< Fortran: dE_edge_test [MeV]
    bool use_offshell_rad{false};    ///< Fortran: use_offshell_rad (use_offshell_rad_int) [flag]
    double Egamma_gen_max{0.0};      ///< Fortran: Egamma_gen_max [MeV] >0.01 => hardwired_rad
    bool do_fermi{false};            ///< Fortran: do_fermi (do_fermi_int) [flag]
    double pt_b_param{0.0};          ///< Fortran: pt_b_param [GeV^-2]
    int sigc_flag{0};                ///< Fortran: sigc_flag [code]
    int sigc_nbin{0};                ///< Fortran: sigc_nbin [count] (title arg 10 is NOT a default; D-005)
    double sigc_kin_min{0.0};        ///< Fortran: sigc_kin_min [z or pt2 GeV^2]
    double sigc_kin_max{0.0};        ///< Fortran: sigc_kin_max [z or pt2 GeV^2] (title arg 1.0 is NOT a default)
    double sigc_kin_ind{0.0};        ///< Fortran: sigc_kin_ind [pt2 GeV^2 or z]
    bool using_tgt_field{false};     ///< Fortran: using_tgt_field (using_tgt_field_int) [flag]
    std::string tgt_field_file;      ///< Fortran: tgt_field_file [character*60 path] (Fortran default: NUL bytes)
    double drift_to_cal{0.0};        ///< Fortran: drift_to_cal [cm]

    // ------------------------------------------------------------------
    // dbase.f:1108-1122 (E_ARM_ACCEPT, P_ARM_ACCEPT)
    // ------------------------------------------------------------------
    /// Acceptance edges of one arm applied to the generated and reconstructed target quantities.
    struct ArmEdge {
        ParamRange delta{-1.0e10, 1.0e10};  ///< Fortran: SPedge%?%delta%min/max [%]
        ParamRange yptar{-1.0e10, 1.0e10};  ///< Fortran: SPedge%?%yptar%min/max [mrad]
        ParamRange xptar{-1.0e10, 1.0e10};  ///< Fortran: SPedge%?%xptar%min/max [mrad]
    };
    /// Acceptance edges of both arms (Fortran SPedge).
    struct SPEdge {
        ArmEdge e;                   ///< Fortran: SPedge%e
        ArmEdge p;                   ///< Fortran: SPedge%p
    } SPedge;                        ///< Fortran: SPedge
};

// ----------------------------------------------------------------------
// Parameter table: one descriptor per regallvars variable (in regallvars order).
// It drives the .inp reader, JSON (de)serialisation, the completeness test and
// the effective-value table. Nothing else needs to know the field list.
// ----------------------------------------------------------------------

/// Storage kind as registered with CTP (regparmint/double/string/intarray/doublearray).
enum class ParamKind {
    Int,          ///< regparmint, stored as int
    Flag,         ///< regparmint on a *_int variable converted by convert_to_logical: stored as bool
    Double,       ///< regparmdouble
    String,       ///< regparmstring (Fortran fixed-length character*N)
    FlagArray,    ///< regparmintarray converted to logical (debug)
    DoubleArray   ///< regparmdoublearray (lambda)
};

/**
 * @brief Description of one input variable: name, storage kind, units, origin and typed accessors into a
 *        simc::SimcParameters.
 *
 * The table of descriptors (simc::ParameterTable) is the single list of the input variables: the `.inp` reader, the
 * JSON (de)serialisation, the completeness test and the effective-value table iterate over it. The strings are static
 * literals; descriptors are immutable and safe to share between threads.
 */
struct ParamDescriptor {
    const char* name;        ///< Fortran/CTP name as registered (case preserved), e.g. "targ%thick"
    ParamKind kind;          ///< Storage kind (selects the accessor to use).
    int size;                ///< array length (debug 6, lambda 3), string length (60/80), else 1
    const char* units;       ///< units as written in the .inp
    const char* fortran;     ///< Fortran target variable and dbase.f registration line
    // Accessors (never null for the matching kind)
    int* (*asInt)(SimcParameters&);          ///< Int: pointer to the field in the given parameters (else null)
    bool* (*asBool)(SimcParameters&);        ///< Flag: the bool; FlagArray: first element
    double* (*asDouble)(SimcParameters&);    ///< Double: the value; DoubleArray: first element
    std::string* (*asString)(SimcParameters&);   ///< String: pointer to the field (else null)

    /**
     * @brief The JSON path of the parameter.
     * @return The components of #name split on '%' (e.g. {"targ", "thick"}).
     */
    std::vector<std::string> jsonPath() const;
};

/**
 * @brief The table of all input variables.
 * @return All 107 descriptors (106 + upstream doing_nuc_elast), in the registration order of the original program (dbase.f:regallvars); a static
 *         table, built on first use (thread-safe).
 */
const std::vector<ParamDescriptor>& ParameterTable();

/**
 * @brief Looks up an input variable by its original name.
 * @param fortran_name The name as written in a `.inp` file, e.g. "targ%thick"; matched case-insensitively (as the
 *        original CTP parser does).
 * @return The descriptor, or nullptr if no variable has this name.
 */
const ParamDescriptor* FindParameter(const std::string& fortran_name);

/**
 * @brief Serialises the parameters to JSON.
 * @param p The parameters.
 * @return A JSON object keyed by the original names split on '%' (nested objects), one entry per variable.
 */
Json ToJson(const SimcParameters& p);

/**
 * @brief Reads parameters from JSON (the inverse of simc::ToJson).
 *
 * Keys are matched case-insensitively against the original names. Missing keys keep the value already in @p p (a
 * default-constructed simc::SimcParameters holds the original defaults).
 * @param j A JSON object in the layout of simc::ToJson.
 * @param p The parameters to update.
 * @param warnings If not null, receives one message per unknown key or type problem (never fatal).
 * @throws std::invalid_argument for a rejected input (electron_arm / hadron_arm = 6, the retired SHMS LSA tune; D-032).
 */
void FromJson(const Json& j, SimcParameters& p, std::vector<std::string>* warnings = nullptr);

/**
 * @brief The value of one parameter as text (for tables and tests).
 * @param d The descriptor of the parameter.
 * @param p The parameters to read.
 * @return The value formatted for its kind (arrays as a comma-separated list).
 */
std::string FormatValue(const ParamDescriptor& d, const SimcParameters& p);

/** @} */

}  // namespace simc

#endif  // SIMC_CORE_SIMCPARAMETERS_H
