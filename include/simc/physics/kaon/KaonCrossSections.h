/**
 * @file
 * @ingroup physics
 * @brief Kaon electroproduction cross sections, a bit-for-bit port of physics_kaon.f and cern/fint.f.
 *
 * Kaon electroproduction cross sections, a bit-for-bit port of physics_kaon.f and cern/fint.f (task P10.1):
 *   peeK            physics_kaon.f:1-168     d5sigma/dEe'dOmega_e'dOmega_K [ub/MeV/sr^2] + the survival probability
 *                                            of the kaon when its decay is not simulated (doing_decay = 0)
 *   sigFactorized   physics_kaon.f:173-237   factorized model of D. Koltenuk (the model SIMC uses for the weight)
 *   eekeek          physics_kaon.f:241-355   Saghai model, K+ Lambda (ntup%sigcm1, the "saghai" ntuple column)
 *   eekeeks         physics_kaon.f:357-477   Saghai model, K+ Sigma tables (every Mrec_struck >= 1150 MeV)
 *   phspwght        physics_kaon.f:479-544   HMS/SOS phase-space weight (dead code in SIMC; ported for the L1 table)
 *   fint            cern/fint.f              CERNLIB multilinear interpolation (KERNLIB E104)
 *   SaghaiTables    dbase.f:624-679          COMMON /sigkaon/ (simulate.inc:188-206): weightc, weightd and the
 *                                            Saghai amplitude tables, read from weightc.dat, weightd.dat,
 *                                            saghai_proton.dat, saghai_sigma0.dat (copies in data/kaon)
 *
 * Numerics. The Fortran is compiled with -fdefault-real-8 -ff2c (D-010, D-024): every literal is double, x**2 is x*x,
 * the operation order is kept. Two build-level facts of the Fortran objects are reproduced (D-043 F1, F4):
 *  * physics_kaon.f declares `real*4 fint`, but cern/fint.f's implicitly typed FUNCTION FINT is REAL(8) and is
 *    returned f2c-style as binary128; the REAL*4 caller reads the low 32 bits of that binary128, which are always 0.
 *    AS BUILT every fint() value seen by eekeek/eekeeks/phspwght is +0.0 (FORTRAN-BUG D-043 F1, see SaghaiOptions).
 *  * DOUBLE COMPLEX z1a..z8a, ur, ui, dble() and the literal 1.d+06 become kind 16 under -fdefault-real-8: the
 *    amplitude sums of eekeek/eekeeks are evaluated in binary128 (__float128, glibc hypotf128 for abs()), each
 *    real*8 subexpression in double, exactly as gfortran types every operation (D-043 F4).
 *
 * Units: the event in MeV / MeV^2 / rad as the Fortran; sigFactorized takes q2 [MeV^2], w [MeV], t [MeV^2], pk [MeV/c],
 * mrec [MeV]; eekeek/eekeeks take ss = s/1e6 [GeV^2], q22 = Q^2 [MeV^2], angles in rad.
 */
#ifndef SIMC_PHYSICS_KAON_KAONCROSSSECTIONS_H
#define SIMC_PHYSICS_KAON_KAONCROSSSECTIONS_H

#include <array>
#include <ostream>
#include <string>
#include <vector>

#include "simc/core/Event.h"

namespace simc {
namespace kaon {

// --------------------------------------------------------------------------------------------- CERNLIB fint
/// @addtogroup physics
/// @{

/// Fortran: cern/fint.f FUNCTION FINT(NARG, ARG, NENT, ENT, TABLE) (KERNLIB E104, C. Letertre / B. Schorr).
/// Multilinear interpolation (linear extrapolation outside the grid) in NARG <= 5 dimensions. ARG, ENT, TABLE are REAL*4
/// as in CERNLIB; X, H, ETA, WEIGHT and the result are REAL(8) as compiled by SIMC (-fdefault-real-8, implicit typing):
/// this is the value cern/fint.o computes. ENT(LMIN+1)-ENT(LMIN) is a REAL*4 difference (both operands REAL*4).
/// NARG outside 1..5: the CERNLIB error branch (KERMTR E104.1); returns 0 (never reached by SIMC).
/// @param narg number of dimensions (1-5).
/// @param arg the point (narg values).
/// @param nent number of grid nodes per dimension.
/// @param ent the grid nodes, dimension after dimension.
/// @param table the table (column-major).
/// @return the interpolated value.
double fint(int narg, const float* arg, const int* nent, const float* ent, const float* table);

// --------------------------------------------------------------------------------------------- /sigkaon/ tables
/// COMMON /sigkaon/ (simulate.inc:188-206), REAL*4 as in the Fortran, column-major (first index fastest) so that the
/// arrays can be handed to fint exactly as the Fortran does (TABLE(I) with the linear index of fint).
///
/// Read-only after load(); thread-safe for reading.
class SaghaiTables {
public:
    static constexpr int kNsP = 10 /**< @copydoc kNa */, kNqP = 11 /**< @copydoc kNa */, kNa = 19;   ///< zrff*(10,11,19): s, Q2, angle
    static constexpr int kNsS = 20 /**< @copydoc kNqS */, kNqS = 10;   ///< zsrff*(20,10,19): s, Q2 (angle as kNa)
    static const std::array<const char*, 4> kFiles;        ///< weightc.dat weightd.dat saghai_proton.dat saghai_sigma0.dat

    /// dbase.f:624-679. `sigma0Corrected` = false reads saghai_sigma0.dat exactly as dbase.f:662-677 (one line off,
    /// D-043 F2); true reads it record by record after the truncated first line (FORTRAN-BUG D-043 F2, see the .cpp).
    /// Warnings of the corrected reader go to `msgs`. Returns false and sets `err` if a file is missing or short.
    /// @param dir the data directory.
    /// @param sigma0Corrected see above.
    /// @param[out] msgs warnings (optional, may be nullptr).
    /// @param[out] err error text (optional, may be nullptr).
    /// @return false on a missing or short file.
    bool load(const std::string& dir, bool sigma0Corrected, std::vector<std::string>* msgs = nullptr,
              std::string* err = nullptr);
    /// @return true after a successful load().
    bool loaded() const { return !weightc.empty(); }

    std::vector<float> weightc;                 ///< weightc(20,50)
    std::vector<float> weightd;                 ///< weightd(8,40,30)
    std::array<std::vector<float>, 6> zrff /**< @copydoc ziff */, ziff;     ///< zrff1..6, ziff1..6 (10,11,19)
    std::array<std::vector<float>, 6> zsrff /**< @copydoc zsiff */, zsiff;   ///< zsrff1..6, zsiff1..6 (20,10,19)

    /// 1-based accessors (Fortran indices)
    /// Element weightc(j, i).
    /// @param j,i 1-based indices.
    /// @return the value.
    float wc(int j, int i) const { return weightc[static_cast<std::size_t>((i - 1) * 20 + (j - 1))]; }
    /// Element weightd(k, j, i).
    /// @param k,j,i 1-based indices.
    /// @return the value.
    float wd(int k, int j, int i) const { return weightd[static_cast<std::size_t>(((i - 1) * 40 + (j - 1)) * 8 + (k - 1))]; }
    /// Linear index into zrff/ziff.
    /// @param is,iq,ia 1-based s, Q2, angle indices.
    /// @return the 0-based linear index.
    static std::size_t idxP(int is, int iq, int ia) {
        return static_cast<std::size_t>(((ia - 1) * kNqP + (iq - 1)) * kNsP + (is - 1));
    }
    /// Linear index into zsrff/zsiff.
    /// @param is,iq,ia 1-based s, Q2, angle indices.
    /// @return the 0-based linear index.
    static std::size_t idxS(int is, int iq, int ia) {
        return static_cast<std::size_t>(((ia - 1) * kNqS + (iq - 1)) * kNsS + (is - 1));
    }
};

/// The D-031 switches of the Saghai model (all false = the Fortran as built, fortran_compat).
struct SaghaiOptions {
    /// FORTRAN-BUG(D-043 F1): use the value fint computes (rounded to REAL*4, as the `real*4 fint` declaration of the
    /// callers intends) instead of the +0.0 the as-built ABI mismatch delivers.
    bool fint_value{false};
    /// FORTRAN-BUG(D-043 F3): eekeek s-grid nodes 9 and 10 = 5.1 and 5.4 GeV^2 (the s of the 9th and 10th blocks of
    /// saghai_proton.dat) instead of 2.6 + 0.3*(i-1) = 5.0 and 5.3.
    bool proton_grid_fixed{false};
};

/// Optional diagnostics of eekeek/eekeeks (L1 tables): the REAL*4 fint arguments px(1:3) and the 12 fint values the
/// routine used, zf(1,k) (real parts) and zf(2,k) (imaginary parts).
struct SaghaiTrace {
    float px[3]{};  ///< the fint arguments px(1:3).
    double zr[6]{} /**< @copydoc zi */, zi[6]{};  ///< the fint values: real and imaginary parts.
};

/// Fortran: physics_kaon.f:eekeek (l.241-355). Mrec_struck = targ%Mrec_struck (COMMON) [MeV]. Returns sigma_eep.
/// @param ss s = W^2 [GeV^2].
/// @param q22 Q2 as peeK passes it (vertex%q2) [MeV^2].
/// @param angl CM angle theta_cm as peeK passes it [rad].
/// @param theta theta_pq [rad].
/// @param phi azimuthal angle (peeK passes its never-assigned phi = 0) [rad].
/// @param epsi virtual photon polarization.
/// @param Mrec_struck recoil (hyperon) mass [MeV].
/// @param tab the Saghai tables.
/// @param opt the D-031 switches.
/// @param[out] trace diagnostics (optional, may be nullptr).
/// @return sigma_eep.
double eekeek(double ss, double q22, double angl, double theta, double phi, double epsi, double Mrec_struck,
              const SaghaiTables& tab, const SaghaiOptions& opt, SaghaiTrace* trace = nullptr);
/// Fortran: physics_kaon.f:eekeeks (l.357-477). Returns sigma_eep.
/// @param ss s = W^2 [GeV^2].
/// @param q22 Q2 as peeK passes it (vertex%q2) [MeV^2].
/// @param angl CM angle theta_cm as peeK passes it [rad].
/// @param theta theta_pq [rad].
/// @param phi azimuthal angle (peeK passes its never-assigned phi = 0) [rad].
/// @param epsi virtual photon polarization.
/// @param Mrec_struck recoil (hyperon) mass [MeV].
/// @param tab the Saghai tables.
/// @param opt the D-031 switches.
/// @param[out] trace diagnostics (optional, may be nullptr).
/// @return sigma_eep.
double eekeeks(double ss, double q22, double angl, double theta, double phi, double epsi, double Mrec_struck,
               const SaghaiTables& tab, const SaghaiOptions& opt, SaghaiTrace* trace = nullptr);

/// Result of physics_kaon.f:phspwght: the weight, and the value the fint call delivered (before the clamps).
struct PhspwghtResult {
    double weight{0};      ///< the phase-space weight
    double fint_seen{0};   ///< the fint value used by the routine (0 as built)
};
/// Fortran: physics_kaon.f:phspwght (l.479-544), NOT called by SIMC (dead code). electron_arm 1 + hadron_arm 2 (fint
/// over weightc) or 2 + 1 (weightd); any other pair prints three lines to `log` and executes a bare STOP
/// (simc::FortranStop).
/// @param delta,theta,phi phase-space coordinates of the routine.
/// @param electron_arm electron arm code.
/// @param hadron_arm hadron arm code.
/// @param tab the tables.
/// @param opt the D-031 switches.
/// @param log stream for the messages (not owned), or nullptr.
/// @return the weight.
/// @throws FortranStop for an unsupported arm pair.
PhspwghtResult phspwght(double delta, double theta, double phi, int electron_arm, int hadron_arm,
                        const SaghaiTables& tab, const SaghaiOptions& opt, std::ostream* log = nullptr);

// --------------------------------------------------------------------------------------------- factorized model
/// SAVEd local of sig_factorized (-fno-automatic): fact_w is assigned only when w**2 != 0.
struct SigFactorizedState {
    double fact_w{0};  ///< the W factor of the previous call.
};
/// Fortran: physics_kaon.f:sig_factorized (l.173-237). The 3-line warning (|pktest-pk|/pk > 0.001) goes to `log`.
/// FORTRAN-DIFF-CANDIDATE(D-043 F10): the proton mass in nu/qcm for every target (Sigma- on the neutron, coherent).
/// @param q2 Q2 [MeV^2].
/// @param w W [MeV].
/// @param t t as peeK passes it (main%t) [MeV^2].
/// @param pk kaon CM momentum [MeV/c].
/// @param mrec recoil mass [MeV].
/// @param st the SAVEd state.
/// @param log stream for the warning (not owned), or nullptr.
/// @return the factorized cross section.
double sigFactorized(double q2, double w, double t, double pk, double mrec, SigFactorizedState& st,
                     std::ostream* log = nullptr);

// --------------------------------------------------------------------------------------------- peeK
/// The run constants that peeK reads from COMMON blocks, and the D-031 switches.
struct PeekConfig {
    double Mtar_struck{0};      ///< targ%Mtar_struck [MeV]
    double Mrec_struck{0};      ///< targ%Mrec_struck [MeV] (hyperon or hypernucleus)
    bool doing_decay{false};    ///< dbase.f doing_decay
    int hadron_arm{1};          ///< dbase.f hadron_arm (spectrometer code)
    double ctau{0};             ///< kaon c*tau [cm]
    double spec_p_P{0};         ///< spec%p%P, central hadron momentum [MeV/c]
    double Mh2{0};              ///< Mh2 [MeV^2]
    /// FORTRAN-BUG(D-043 F5): beta*gamma from the kaon momentum p_kaon (not the central momentum).
    bool survival_p_kaon{false};
    /// FORTRAN-BUG(D-043 F6): zaero assigned for every arm (0 for the SHMS and the other arms without an aerogel
    /// distance) and the slope correction sqrt(1+dx^2+dy^2).
    bool survival_path_fixed{false};
    SaghaiOptions saghai;  ///< Saghai model switches.
};
/// SAVEd locals of peeK (-fno-automatic): zaero keeps its value between calls (hadron_arm 5/6 never assign it), phi is
/// never assigned (0), and the sig_factorized state.
struct PeekState {
    double zaero{0};  ///< aerogel distance of the previous call [cm].
    double phi{0};  ///< never assigned (0).
    SigFactorizedState sf;  ///< sig_factorized state.
};

/// Fortran: physics_kaon.f:peeK (l.1-168). Reads vertex (Ein, nu, q, Q2, e%E, p%P, p%E, up, uq), main%t, epsilon,
/// theta_pq, phi_pq, FP%p (path, dx, dy), SP%p%delta and main.pfermi (COMMON /pfermi/). Writes main%thetacm, phicm,
/// pcm, davejac, johnjac, wcm, ntup%sigcm1 (Saghai), sigcm2 = sigcm (factorized) and, when !doing_decay, the survival
/// probability into `survivalprob` and the COMMON decdist (main.ntup.decdist), physics_kaon.f:164.
/// The |beta| > 1 warning of transform_to_cm and the sig_factorized warning go to `log` (if non-null).
/// @param vertex the vertex.
/// @param main the main event record (written).
/// @param[in,out] survivalprob the survival probability.
/// @param cfg the run constants.
/// @param tab the Saghai tables.
/// @param st the SAVEd state.
/// @param log stream for the warnings (not owned), or nullptr.
/// @return the cross section.
double peeK(const Event& vertex, MainEvent& main, double& survivalprob, const PeekConfig& cfg,
            const SaghaiTables& tab, PeekState& st, std::ostream* log = nullptr);

/// @}

}  // namespace kaon
}  // namespace simc

#endif  // SIMC_PHYSICS_KAON_KAONCROSSSECTIONS_H
