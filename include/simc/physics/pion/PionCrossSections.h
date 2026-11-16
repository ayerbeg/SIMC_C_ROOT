/**
 * @file
 * @ingroup physics
 * @brief Exclusive pion electroproduction cross sections, a bit-for-bit port of physics_pion.f.
 *
 * Exclusive pion electroproduction cross sections, a bit-for-bit port of physics_pion.f (task P9.1):
 *   peepi           physics_pion.f:1-194     d5sigma/dEe'dOmega_e'dOmega_pi [ub/MeV/sr^2]; the model SIMC uses
 *   sigParam2021    physics_pion.f:738-805   P. Bosted 2021 fit (two exclfit calls for pi0), COMMON /excmn/
 *   exclfit         physics_pion.f:807-854
 *   sigmaid         physics_pion.f:577-736   MAID2007 tables (maid*.dat), nearest bin, used by peepi for W < 2 GeV
 *   sigBlok         physics_pion.f:404-465   } not called by SIMC (their calls in peepi are commented out,
 *   sigParam04      physics_pion.f:467-516   } physics_pion.f:114-123); ported for completeness and validated
 *   sigParam3000    physics_pion.f:519-573   } against the P9.F1 L1 tables
 *   SigMultipole    physics_pion.f:197-398   multipole model (dead code too, physics_pion.f:157-163), COMPLEX*8
 * transform_to_cm (jacobians.f) is NOT duplicated here: peepi calls reactions::transformToCm (task P4.1).
 *
 * Numerics (D-010, D-024): the Fortran is compiled with -fdefault-real-8, so every literal is double; x**2 is x*x,
 * real exponents use pow, and the Fortran operation order is kept everywhere. sig_multipole keeps the single
 * precision of its COMPLEX*8 amplitudes (D-038 F12): see SigMultipole.
 *
 * Units: the cross-section functions take GeV units (t, Q2, s in GeV^2) exactly as the Fortran routines; peepi takes
 * the event in MeV (vertex, main) like the Fortran.
 */
#ifndef SIMC_PHYSICS_PION_PIONCROSSSECTIONS_H
#define SIMC_PHYSICS_PION_PIONCROSSSECTIONS_H

#include <array>
#include <complex>
#include <ostream>
#include <string>
#include <vector>

#include "simc/core/Event.h"

namespace simc {
namespace pion {

// --------------------------------------------------------------------------------------------- simple fits
/// @addtogroup physics
/// @{

/// Fortran: physics_pion.f:sig_blok (l.404-465). Arguments as the Fortran (t = -t_Mandelstam >= 0 [GeV^2],
/// q2_gev, s_gev [GeV^2], mtar_gev [GeV]); returns dsigma/dt/dphi_cm [ub/MeV^2/rad].
/// @param thetacm,phicm CM angles of the pion [rad].
/// @param t -t (Mandelstam) [GeV^2].
/// @param q2_gev Q2 [GeV^2].
/// @param s_gev s = W^2 [GeV^2].
/// @param eps virtual photon polarization epsilon.
/// @param mtar_gev target nucleon mass [GeV].
/// @param which_pion which_pion of the input.
/// @return dsigma/dt/dphi_cm [ub/MeV^2/rad].
double sigBlok(double thetacm, double phicm, double t, double q2_gev, double s_gev, double eps, double mtar_gev,
               int which_pion);
/// Fortran: physics_pion.f:sig_param04 (l.467-516). FORTRAN-QUIRK(D-038 F13): pole at |t| = 0.02 GeV^2.
/// @param thetacm,phicm CM angles of the pion [rad].
/// @param t -t (Mandelstam) [GeV^2].
/// @param q2_gev Q2 [GeV^2].
/// @param s_gev s = W^2 [GeV^2].
/// @param eps virtual photon polarization epsilon.
/// @param mtar_gev target nucleon mass [GeV].
/// @param which_pion which_pion of the input.
/// @return dsigma/dt/dphi_cm [ub/MeV^2/rad].
double sigParam04(double thetacm, double phicm, double t, double q2_gev, double s_gev, double eps, double mtar_gev,
                  int which_pion);
/// Fortran: physics_pion.f:sig_param_3000 (l.519-573).
/// @param thetacm,phicm CM angles of the pion [rad].
/// @param t -t (Mandelstam) [GeV^2].
/// @param q2_gev Q2 [GeV^2].
/// @param s_gev s = W^2 [GeV^2].
/// @param eps virtual photon polarization epsilon.
/// @param mtar_gev target nucleon mass [GeV].
/// @param which_pion which_pion of the input.
/// @return dsigma/dt/dphi_cm [ub/MeV^2/rad].
double sigParam3000(double thetacm, double phicm, double t, double q2_gev, double s_gev, double eps, double mtar_gev,
                    int which_pion);

// --------------------------------------------------------------------------------------------- 2021 fit
/// The outputs of physics_pion.f:exclfit.
struct ExclfitResult {
    double sigl{0}, sigt{0}, sigtt{0}, siglt{0}, sig{0};  ///< L, T, TT, LT components and the total cross section.
};
/// Fortran: physics_pion.f:exclfit (l.807-854). p = p(1:17) at [0..16]. mtar_gev is the local 0.938 and pi the
/// literal 3.1415928 of the Fortran (FORTRAN-DIFF-CANDIDATE D-038 F4: the nucleon mass also for coherent production).
/// @param t t as peepi passes it (main%t/1e6) [GeV^2].
/// @param thetacm,phicm CM angles [rad].
/// @param q2_gev Q2 [GeV^2].
/// @param s_gev s = W^2 [GeV^2].
/// @param eps virtual photon polarization.
/// @param p fit parameters.
/// @param fpifact pion form factor scale.
/// @return the cross-section components.
ExclfitResult exclfit(double t, double thetacm, double phicm, double q2_gev, double s_gev, double eps,
                      const std::array<double, 17>& p, double fpifact);

/// Fortran: physics_pion.f:738-786, the parameter vectors pp (pi+) and pm (pi-) of sig_param_2021.
extern const std::array<double, 17> kParam2021PiPlus;
/// The pi- parameter vector of sig_param_2021.
extern const std::array<double, 17> kParam2021PiMinus;

/// Fortran: COMMON /excmn/ sigl, sigt, sigttv, sigltv (physics_pion.f:745), the components of the last exclfit call.
struct Excmn {
    double sigl{0}, sigt{0}, sigttv{0}, sigltv{0};  ///< L, T, TT, LT components of the last exclfit call.
};
/// Fortran: physics_pion.f:sig_param_2021 (l.738-805). pi0 (doing_pizero): the average of the pi+ and pi- fits with
/// fpifact = 0; pi- for which_pion 1, 11, 3; pi+ otherwise. `excmn` (optional) receives COMMON /excmn/ after the call
/// (for pi0: the values of the second, pi-, exclfit call).
/// @param thcm,phicm CM angles [rad].
/// @param t t as peepi passes it (main%t/1e6) [GeV^2].
/// @param q2 Q2 [GeV^2].
/// @param wsq W^2 [GeV^2].
/// @param eps virtual photon polarization.
/// @param which_pion which_pion of the input.
/// @param doing_pizero pi0 production.
/// @param[out] excmn COMMON /excmn/ after the call (optional, may be nullptr).
/// @return the cross section.
double sigParam2021(double thcm, double phicm, double t, double q2, double wsq, double eps, int which_pion,
                    bool doing_pizero, Excmn* excmn = nullptr);

// --------------------------------------------------------------------------------------------- MAID
/// The four MAID2007 tables of physics_pion.f:sigmaid, maidtbl(4,25,46,23,18) (l.597). Read once (the Fortran reads
/// them on the first sigmaid call, l.616-637; SIMC++ reads them when the model is initialised).
///
/// Read-only after load(); thread-safe for reading.
class MaidTables {
public:
    static constexpr int kNipi = 4, kNq = 25, kNw = 46, kNth = 23, kNcol = 18;  ///< table dimensions: channels, Q2, W, theta bins and columns.
    /// File names in the order of ipi = 1..4 (physics_pion.f:620-623).
    static const std::array<const char*, 4> kFiles;

    /// Read <dir>/maidpi0p.dat, maidpi0n.dat, maidpipn.dat, maidpimp.dat with the Fortran format (f11.6,18f8.4).
    /// Returns false and sets `err` if a file is missing or short.
    /// @param dir the data directory.
    /// @param[out] err error text (optional, may be nullptr).
    /// @return false on a missing or short file.
    bool load(const std::string& dir, std::string* err = nullptr);
    /// @return true after a successful load().
    bool loaded() const { return !tbl_.empty(); }
    /// maidtbl(ipi, iq, iw, ith, j), 1-based as the Fortran.
    /// @param ipi channel (1-4).
    /// @param iq,iw,ith Q2, W and theta bins (1-based).
    /// @param j column (1-based).
    /// @return the table value.
    double at(int ipi, int iq, int iw, int ith, int j) const {
        return tbl_[static_cast<std::size_t>(((((ipi - 1) * kNq + (iq - 1)) * kNw + (iw - 1)) * kNth + (ith - 1)) * kNcol +
                                             (j - 1))];
    }

private:
    std::vector<double> tbl_;
};

/// Outputs of physics_pion.f:sigmaid. Only sig0 is returned by the Fortran routine; sigz, sigez and sige are its
/// local polarisation terms (computed for completeness).
struct SigmaidResult {
    double sig0{0}, sigz{0}, sigez{0}, sige{0};  ///< dsigma/dOmega* [ub/sr] and the local polarisation terms.
};
/// The nearest-bin indices of sigmaid (physics_pion.f:659-669), 1-based.
struct MaidBin {
    int iq{1}, iw{1}, ith{1};  ///< Q2, W and theta bins.
};
/// Fortran: physics_pion.f:659-669. FORTRAN-QUIRK(D-038 F14): nearest (lower) bin, no interpolation; only the theta
/// rows 1..6 of the 23 stored per (Q2, W) are reachable.
/// @param q2 Q2 [GeV^2].
/// @param w W [GeV].
/// @param costh cos(theta*).
/// @return the bins.
MaidBin sigmaidBin(double q2, double w, double costh);
/// Fortran: physics_pion.f:sigmaid (l.577-736). ipi 1 pi0 p, 2 pi0 n, 3 pi+ n, 4 pi- p; q2 [GeV^2], w [GeV], e0 [GeV],
/// costh = cos(theta*), phi [rad]. Returns sig0 = dsigma/dOmega* [ub/sr] (0 for the early returns l.645-650).
/// @param ipi channel (see above).
/// @param q2 Q2 [GeV^2].
/// @param w W [GeV].
/// @param e0 beam energy [GeV].
/// @param costh cos(theta*).
/// @param phi azimuthal angle [rad].
/// @param maid the MAID tables.
/// @return the cross section and polarisation terms.
SigmaidResult sigmaid(int ipi, double q2, double w, double e0, double costh, double phi, const MaidTables& maid);

// --------------------------------------------------------------------------------------------- multipoles
/// Fortran: physics_pion.f:sig_multipole (l.197-398), with its SAVEd state (one object = one Fortran process).
///
/// FORTRAN-DIFF-CANDIDATE(D-038 F12), reproduced in both modes (Question: are the single precision and the reuse of
/// the previous call's amplitudes intended? The routine is not called by SIMC):
///  * the amplitudes are COMPLEX*8 (single precision) although the build uses -fdefault-real-8. They are kept as
///    std::complex<float>; every mixed expression follows the Fortran kind rules (a real*8 operand promotes the
///    expression to COMPLEX(8), the assignment rounds back to COMPLEX(4)), with the arithmetic written out
///    component by component as gfortran compiles it (-fcx-fortran-rules, real operands not promoted to complex);
///  * the interpolation assigns amplitudes only strictly inside a (Q2, nu) cell: at exact nodes or outside the table
///    the amplitudes of the PREVIOUS call are used (rows must be evaluated in order);
///  * `nu` is an in/out argument: nu >= 500 prints a warning and sets nu = 499.9;
///  * the conversion factor 197.327/mh is fixed at the first call (SAVEd), as are the Clebsch-Gordan coefficients.
class SigMultipole {
public:
    /// Read MULTIPOL.file (the Fortran reads it at the first call, l.255-286). Returns false and sets `err` on error.
    /// @param path the MULTIPOL file.
    /// @param[out] err error text (optional, may be nullptr).
    /// @return false on error.
    bool load(const std::string& path, std::string* err = nullptr);
    /// @return true after a successful load().
    bool loaded() const { return loaded_; }

    /// One sig_multipole call [ub/MeV^2/rad]. Messages (nu >= 500) go to `log` if non-null.
    /// @param[in,out] nu photon energy in the nucleon rest frame [MeV] (clamped to 499.9).
    /// @param efer energy of the struck nucleon [MeV].
    /// @param qstar virtual photon momentum in the CM [MeV].
    /// @param tfcos cosine of the CM angle.
    /// @param ppicm pion CM momentum [MeV].
    /// @param s_gev s [GeV^2].
    /// @param thetacm,phicm CM angles [rad].
    /// @param eps virtual photon polarization.
    /// @param final_state final state (1-based).
    /// @param q2_gev Q2 [GeV^2].
    /// @param mtar_gev target mass [GeV].
    /// @param pfer Fermi momentum [MeV/c].
    /// @param mh pion mass [MeV].
    /// @param log stream for the messages (not owned), or nullptr.
    /// @return the cross section.
    double operator()(double& nu, double efer, double qstar, double tfcos, double ppicm, double s_gev, double thetacm,
                      double phicm, double eps, int final_state, double q2_gev, double mtar_gev, double pfer, double mh,
                      std::ostream* log = nullptr);

private:
    bool loaded_{false};
    bool first_{true};
    double root_two_{0}, mpi_to_microbarn_{0};
    double cleb_[2][3]{};                            // CLEB(2,3)
    double q2_table_[8]{}, nu_table_[10]{};
    std::complex<float> at_[7][3][8][10]{};          // AT(7,3,8,10)
    std::complex<float> amplitude_[7]{};             // SAVEd (zero-initialised static storage in the Fortran)
};

// --------------------------------------------------------------------------------------------- peepi
/// The run constants that peepi reads from COMMON blocks.
struct PeepiConfig {
    int which_pion{0};           ///< dbase.f which_pion
    bool doing_pizero{false};    ///< dbase.f doing_pizero
    double Mtar_struck{0};       ///< targ%Mtar_struck [MeV]
    bool fortran_compat{false};  ///< D-031 (only D-038 F10 depends on it)
};

/// Fortran: physics_pion.f:peepi (l.1-194). d5sigma/dEe'dOmega_e'dOmega_pi [ub/MeV/sr^2].
/// Reads vertex (Ein, nu, q, Q2, e%E, p%P, p%E, up, uq), main%t, main%epsilon, main%phi_pq and main.pfermi (COMMON
/// /pfermi/). Writes main%thetacm, phicm, pcm, davejac, johnjac, wcm and main.ntup.sigcm1, sigcm2, sigcm, as the
/// Fortran. The |beta| > 1 warning of transform_to_cm goes to `log` (if non-null).
///
/// FORTRAN-BUG(D-038 F10): physics_pion.f:133-154 assigns ntup%sigcm2 only for W < 2000 MeV, so for W >= 2 GeV it
/// keeps the value of an earlier event. Corrected (default): sigcm2 = 0 for W >= 2 GeV (not evaluated);
/// cfg.fortran_compat keeps the stale value (main.ntup.sigcm2 is left unchanged).
/// @param vertex the vertex.
/// @param main the main event record (written).
/// @param cfg the run constants.
/// @param maid the MAID tables.
/// @param log stream for the warning (not owned), or nullptr.
/// @return the cross section.
double peepi(const Event& vertex, MainEvent& main, const PeepiConfig& cfg, const MaidTables& maid,
             std::ostream* log = nullptr);

/// @}

}  // namespace pion
}  // namespace simc

#endif  // SIMC_PHYSICS_PION_PIONCROSSSECTIONS_H
