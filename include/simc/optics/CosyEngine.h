/**
 * @file
 * @ingroup optics
 * @brief simc::optics::CosyEngine: the spectrometer-independent COSY transport engine (transp, project, rotate, recon).
 *
 * Spectrometer-independent COSY transport engine (task P6.1). Ports:
 *   transp()        Fortran: shared/transp.f:transp (l.1-281)          one forward class (map or drift)
 *   loadForward()   Fortran: shared/transp.f:transp_init / transp_init_file (entries, l.285-478)
 *   recon()         Fortran: hms/mc_hms_recon.f (l.1-146), shms/mc_shms_recon.f (l.1-166), sos/mc_sos_recon.f (l.1-150)
 *   project()       Fortran: shared/project.f (l.1-127)                 field-free drift
 *   rotateHaxis()   Fortran: shared/rotate_haxis.f                      (raddeg = 0.017453292, D-015)
 *   rotateVaxis()   Fortran: shared/rotate_vaxis.f
 * shared/locforunt.f only finds a free Fortran I/O unit; it has no C++ counterpart (std::ifstream).
 *
 * Units and coordinates as in the Fortran (TRANSPORT: x down/dispersive, y left, z downstream):
 * positions cm, slopes dx/dz, dy/dz, dpps in % (the transp.f comment "Fractional" is wrong), path cm.
 * Spectrometer ids as in SIMC: 1 HMS, 2 SOS, 3 HRSR, 4 HRSL, 5 SHMS (6 = SHMS LSA -> STOP).
 *
 * Numerics: Fortran real**integer*2 is __builtin_powi -> libgcc __powidf2 (repeated squaring), NOT pow();
 * powi() below reproduces it bit for bit. Everything is real*8.
 */
#ifndef SIMC_OPTICS_COSYENGINE_H
#define SIMC_OPTICS_COSYENGINE_H

#include <array>
#include <iosfwd>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include "simc/core/Errors.h"
#include "simc/optics/TransportTrace.h"

namespace simc {

class IRandom;

namespace optics {


/// @addtogroup optics
/// @{

/// Fortran: COMMON /track/ (spectrometers.inc l.23-25): the ray carried through transp/project/rotate/recon.
struct Track {
    double xs{0} /**< @copydoc zs */, ys{0} /**< @copydoc zs */, zs{0};   ///< cm
    double dxdzs{0} /**< @copydoc dydzs */, dydzs{0};    ///< slopes
    double dpps{0};               ///< delta p/p in %
};

/// Decay bookkeeping of transp/project (their arguments decay_flag, dflag, m2, ph plus the COMMON variables
/// ctau, decdist, Mh2_final). Pass nullptr to transp()/project() for decay_flag = .false.
struct DecayState {
    bool decay_flag{false};   ///< check for decay
    bool dflag{false};        ///< particle has already decayed
    double m2{0};             ///< mass^2 of the particle (MeV^2), updated on decay
    double ph{0};             ///< momentum (MeV/c), updated on decay
    double ctau{0};           ///< Fortran COMMON ctau (cm)
    double decdist{0};        ///< Fortran COMMON decdist (cm), accumulated
    double Mh2_final{0};      ///< Fortran COMMON Mh2_final (MeV^2), set on decay
    IRandom* rng{nullptr};    ///< grnd() source (required if decay_flag)
    /// D-031 run option for project() (transp() uses CosyEngine::fortranCompat()): true = the Fortran behaviour,
    /// false = the D-037 corrections (F7: project.f:113 drift remainder with the pre-decay slopes; F9: only pi/K decay).
    bool fortran_compat{false};
    TransportTrace* trace{nullptr};   ///< optional diagnostics (TransportTrace.h); no effect on the numerics
};

/// D-037 F9: the particles project/transp can decay (|sqrt(m2) - Mpi| < 2 or |sqrt(m2) - Mk| < 2, the tests of
/// project.f:78-79 / transp.f:152-153 / 245-246).
/// @param m2 particle mass squared [MeV^2].
/// @return true for a pion or a kaon.
bool decaysInFlight(double m2);

/// Fortran real**integer (libgcc __powidf2, used by gfortran for x**n with a variable integer n).
/// @param x the base.
/// @param m the integer exponent.
/// @return x**m with the libgcc rounding.
double powi(double x, int m);

/// Fortran: shared/project.f:project(x_new,y_new,z_drift,decay_flag,dflag,m2,ph,pathlen).
/// x_new/y_new are usually t.xs/t.ys (the Fortran passes the COMMON variables).
/// @param[in,out] t the track (slopes; dpps on decay).
/// @param[in,out] x_new x position to update [cm] (usually t.xs).
/// @param[in,out] y_new y position to update [cm] (usually t.ys).
/// @param z_drift drift length [cm].
/// @param[in,out] pathlen accumulated path length [cm].
/// @param decay decay state (not owned), or nullptr for no decay.
void project(Track& t, double& x_new, double& y_new, double z_drift, double& pathlen, DecayState* decay = nullptr);
/// project(xs, ys, ...) — the usual call.
/// @param[in,out] t the track.
/// @param z_drift drift length [cm].
/// @param[in,out] pathlen accumulated path length [cm].
/// @param decay decay state (not owned), or nullptr for no decay.
inline void project(Track& t, double z_drift, double& pathlen, DecayState* decay = nullptr) {
    project(t, t.xs, t.ys, z_drift, pathlen, decay);
}

/// Fortran: shared/rotate_haxis.f:rotate_haxis(rotang,xp0,yp0). rotang in degrees; uses t.dxdzs, t.dydzs.
/// `fortran_compat`: true = the truncated raddeg of the Fortran (kRotateRaddeg), false = pi/180 (FORTRAN-BUG D-015).
/// @param t the track.
/// @param rotang rotation angle [deg].
/// @param[out] xp0 rotated x [cm].
/// @param[out] yp0 rotated y [cm].
/// @param fortran_compat see above.
void rotateHaxis(const Track& t, double rotang, double& xp0, double& yp0, bool fortran_compat);
/// Fortran: shared/rotate_vaxis.f:rotate_vaxis(rotang,xp0,yp0).
/// @param t the track.
/// @param rotang rotation angle [deg].
/// @param[out] xp0 rotated x [cm].
/// @param[out] yp0 rotated y [cm].
/// @param fortran_compat true: the truncated raddeg of the Fortran.
void rotateVaxis(const Track& t, double rotang, double& xp0, double& yp0, bool fortran_compat);

/// Fortran rotate_*axis.f: parameter (raddeg=0.017453292) — truncated degree->radian factor (D-015 F9), used with
/// fortran_compat only.
constexpr double kRotateRaddeg = 0.017453292;

/// @brief COSY forward maps (transp) and reconstruction maps (recon) of the spectrometers, keyed by spectrometer number.
///
/// Loading is not thread-safe; after loading, transp() and recon() only read the maps (transp() writes the
/// caller's Track/DecayState). Owns the maps.
class CosyEngine {
public:
    /// transp.f: max_class = 41 (spectrometers.inc), max_elements = 1000, coeff_min = 1.0e-14.
    static constexpr int kMaxClass = 41;  ///< maximum number of COSY classes.
    static constexpr int kMaxElements = 1000;  ///< maximum number of matrix elements per class.
    static constexpr double kCoeffMin = 1.0e-14;  ///< coefficients below this are dropped.

    /// One forward matrix element: coefficients for X, A, Y, B, dL; exponents of x, a, y, b, delta
    /// (the TOF digit is dropped at read time, transp.f:377-391).
    struct ForwardTerm {
        std::array<double, 5> coeff{};  ///< coefficients.
        std::array<int, 5> expon{};  ///< exponents.
    };
    /// One COSY class (element) of a forward map.
    struct ForwardClass {
        std::vector<ForwardTerm> terms;  ///< matrix elements.
        double length{0};      ///< !LENGTH: comment in cm (0 if absent)   (transp.f local `length`)
        bool adrift{true};     ///< COMMON /driftdistances/ adrift
        double driftdist{0};   ///< COMMON /driftdistances/ driftdist (cm)
        std::array<int, 5> max_expon{};   ///< largest exponent of each variable in `terms` (P15.4 power tables)
    };
    /// One recon element: coefficients for delta_phi (xptar), y_tgt, delta_t (yptar), delta_p;
    /// exponents of xfp, xpfp, yfp, ypfp, fry.
    struct ReconTerm {
        std::array<double, 4> coeff{};  ///< coefficients.
        std::array<int, 5> expon{};  ///< exponents.
    };

    /// Creates an engine without maps (log on std::cout).
    CosyEngine();

    /// Default forward file of transp_init (relative to the SIMC directory), transp.f:300-314.
    /// @param spectr spectrometer number.
    /// @return the relative path.
    static std::string defaultForwardFile(int spectr);
    /// Default recon file (mc_hms_recon.f / mc_sos_recon.f / mc_shms_recon.f). 1 (HMS), 2 (SOS) and 5 (SHMS) are ported.
    /// @param spectr spectrometer number.
    /// @return the relative path.
    static std::string defaultReconFile(int spectr);

    /// Fortran: transp_init_file(file, n_classes) for spectrometer `spectr`. Returns n_classes.
    /// @param spectr spectrometer number.
    /// @param path the forward map file.
    /// @return the number of classes.
    /// @throws std::runtime_error (or FortranStop) if the file cannot be read.
    int loadForward(int spectr, const std::string& path);
    /// Fortran: transp_init(spectr, n_classes) with the default file below `simcDir`.
    /// @param spectr spectrometer number.
    /// @param simcDir the SIMC data directory.
    /// @return the number of classes.
    int initForward(int spectr, const std::string& simcDir);
    /// Fortran: the firsttime block of mc_hms_recon / mc_shms_recon. Returns n_terms.
    /// @param spectr spectrometer number.
    /// @param path the recon map file.
    /// @return the number of terms.
    int loadRecon(int spectr, const std::string& path);
    /// Load the default recon file below @p simcDir.
    /// @param spectr spectrometer number.
    /// @param simcDir the SIMC data directory.
    /// @return the number of terms.
    int initRecon(int spectr, const std::string& simcDir);

    /// @param spectr spectrometer number.
    /// @return true if a forward map is loaded.
    bool hasForward(int spectr) const { return fwd_.count(spectr) != 0; }
    /// @param spectr spectrometer number.
    /// @return true if a recon map is loaded.
    bool hasRecon(int spectr) const { return rec_.count(spectr) != 0; }
    /// @param spectr spectrometer number.
    /// @return the number of forward classes.
    int nClasses(int spectr) const;
    /// @param spectr spectrometer number.
    /// @param cls class number (1-based).
    /// @return the class.
    const ForwardClass& forwardClass(int spectr, int cls) const;   ///< cls 1-based
    /// @param spectr spectrometer number.
    /// @param cls class number (1-based).
    /// @return true if the class is a drift.
    bool adrift(int spectr, int cls) const { return forwardClass(spectr, cls).adrift; }
    /// @param spectr spectrometer number.
    /// @param cls class number (1-based).
    /// @return the drift length [cm].
    double driftdist(int spectr, int cls) const { return forwardClass(spectr, cls).driftdist; }
    /// @param spectr spectrometer number.
    /// @return the recon matrix elements.
    const std::vector<ReconTerm>& reconTerms(int spectr) const;
    /// @param spectr spectrometer number.
    /// @return the largest total order of the recon map.
    int reconMaxOrder(int spectr) const;

    /// Fortran: transp(spectr,class,decay_flag,dflag,m2,ph,zd,pathlen). Updates t (xs, dxdzs, ys, dydzs; dpps
    /// only on decay) and adds zd + delta_z to pathlen. Returns delta_z (= -sum(5)).
    /// @param spectr spectrometer number.
    /// @param cls class number (1-based).
    /// @param[in,out] t the track.
    /// @param zd drift length of the class [cm].
    /// @param[in,out] pathlen accumulated path length [cm].
    /// @param decay decay state (not owned), or nullptr for no decay.
    /// @return delta_z [cm].
    double transp(int spectr, int cls, Track& t, double zd, double& pathlen, DecayState* decay = nullptr);

    /// Fortran: mc_hms_recon / mc_shms_recon(delta_p,delta_t,delta_phi,y_tgt,fry[,spectr]) on the focal-plane
    /// track t (xs, dxdzs, ys, dydzs). Outputs: delta_p [%], delta_t = yptar, delta_phi = xptar [slopes], y_tgt [cm].
    /// FORTRAN-QUIRK (D-015): |fry/100| <= 1e-30 is replaced by 1e-30 (avoid 0**0); for the SOS (spectr 2,
    /// mc_sos_recon.f:120-122) every one of the five hut values.
    /// @param spectr spectrometer number.
    /// @param t the focal-plane track.
    /// @param fry vertical beam position (raster) [cm].
    /// @param[out] delta_p delta [%].
    /// @param[out] delta_t yptar.
    /// @param[out] delta_phi xptar.
    /// @param[out] y_tgt ytar [cm].
    void recon(int spectr, const Track& t, double fry, double& delta_p, double& delta_t, double& delta_phi,
               double& y_tgt) const;

    /// Where the Fortran write(6,*) warnings go (default std::cout; nullptr = silent).
    /// @param os the stream (not owned), or nullptr.
    void setLog(std::ostream* os) { log_ = os; }
    /// D-031 run option: true = keep the Fortran bugs of the transport code (transp.f:158 kaon decay mass, D-006;
    /// truncated raddeg of rotate_*axis.f, D-015 F9, used by the chains). Default false (corrected).
    /// @param b true for the Fortran behaviour.
    void setFortranCompat(bool b) { fortranCompat_ = b; }
    /// @return the D-031 run option.
    bool fortranCompat() const { return fortranCompat_; }

private:
    std::map<int, std::vector<ForwardClass>> fwd_;
    std::map<int, std::vector<ReconTerm>> rec_;
    std::map<int, int> recMaxOrder_;
    std::map<int, std::array<int, 5>> recMaxExpon_;   ///< largest exponent of each hut variable (P15.4 power tables)
    std::ostream* log_;
    bool fortranCompat_{false};
};

namespace detail {
/// Fortran formatted input of a real with Gw.d editing (blanks ignored, implied decimal point if none).
/// @param field the input field.
/// @param d the d of Gw.d.
/// @return the value.
double readFortranReal(const std::string& field, int d);
/// Fortran I1 input (blank -> 0).
/// @param field the input field.
/// @return the value.
int readFortranInt(const std::string& field);
}  // namespace detail

/// @}

}  // namespace optics
}  // namespace simc

#endif  // SIMC_OPTICS_COSYENGINE_H
