// include/simc/optics/CosyEngine.h
//
// Spectrometer-independent COSY transport engine (task P6.1). Ports:
//   transp()        Fortran: shared/transp.f:transp (l.1-281)          one forward class (map or drift)
//   loadForward()   Fortran: shared/transp.f:transp_init / transp_init_file (entries, l.285-478)
//   recon()         Fortran: hms/mc_hms_recon.f (l.1-146), shms/mc_shms_recon.f (l.1-166), sos/mc_sos_recon.f (l.1-150)
//   project()       Fortran: shared/project.f (l.1-127)                 field-free drift
//   rotateHaxis()   Fortran: shared/rotate_haxis.f                      (raddeg = 0.017453292, D-015)
//   rotateVaxis()   Fortran: shared/rotate_vaxis.f
// shared/locforunt.f only finds a free Fortran I/O unit; it has no C++ counterpart (std::ifstream).
//
// Units and coordinates as in the Fortran (TRANSPORT: x down/dispersive, y left, z downstream):
// positions cm, slopes dx/dz, dy/dz, dpps in % (the transp.f comment "Fractional" is wrong), path cm.
// Spectrometer ids as in SIMC: 1 HMS, 2 SOS, 3 HRSR, 4 HRSL, 5 SHMS (6 = SHMS LSA -> STOP).
//
// Numerics: Fortran real**integer*2 is __builtin_powi -> libgcc __powidf2 (repeated squaring), NOT pow();
// powi() below reproduces it bit for bit. Everything is real*8.
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


/// Fortran: COMMON /track/ (spectrometers.inc l.23-25): the ray carried through transp/project/rotate/recon.
struct Track {
    double xs{0}, ys{0}, zs{0};   ///< cm
    double dxdzs{0}, dydzs{0};    ///< slopes
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
bool decaysInFlight(double m2);

/// Fortran real**integer (libgcc __powidf2, used by gfortran for x**n with a variable integer n).
double powi(double x, int m);

/// Fortran: shared/project.f:project(x_new,y_new,z_drift,decay_flag,dflag,m2,ph,pathlen).
/// x_new/y_new are usually t.xs/t.ys (the Fortran passes the COMMON variables).
void project(Track& t, double& x_new, double& y_new, double z_drift, double& pathlen, DecayState* decay = nullptr);
/// project(xs, ys, ...) — the usual call.
inline void project(Track& t, double z_drift, double& pathlen, DecayState* decay = nullptr) {
    project(t, t.xs, t.ys, z_drift, pathlen, decay);
}

/// Fortran: shared/rotate_haxis.f:rotate_haxis(rotang,xp0,yp0). rotang in degrees; uses t.dxdzs, t.dydzs.
/// `fortran_compat`: true = the truncated raddeg of the Fortran (kRotateRaddeg), false = pi/180 (FORTRAN-BUG D-015).
void rotateHaxis(const Track& t, double rotang, double& xp0, double& yp0, bool fortran_compat);
/// Fortran: shared/rotate_vaxis.f:rotate_vaxis(rotang,xp0,yp0).
void rotateVaxis(const Track& t, double rotang, double& xp0, double& yp0, bool fortran_compat);

/// Fortran rotate_*axis.f: parameter (raddeg=0.017453292) — truncated degree->radian factor (D-015 F9), used with
/// fortran_compat only.
constexpr double kRotateRaddeg = 0.017453292;

class CosyEngine {
public:
    /// transp.f: max_class = 41 (spectrometers.inc), max_elements = 1000, coeff_min = 1.0e-14.
    static constexpr int kMaxClass = 41;
    static constexpr int kMaxElements = 1000;
    static constexpr double kCoeffMin = 1.0e-14;

    /// One forward matrix element: coefficients for X, A, Y, B, dL; exponents of x, a, y, b, delta
    /// (the TOF digit is dropped at read time, transp.f:377-391).
    struct ForwardTerm {
        std::array<double, 5> coeff{};
        std::array<int, 5> expon{};
    };
    struct ForwardClass {
        std::vector<ForwardTerm> terms;
        double length{0};      ///< !LENGTH: comment in cm (0 if absent)   (transp.f local `length`)
        bool adrift{true};     ///< COMMON /driftdistances/ adrift
        double driftdist{0};   ///< COMMON /driftdistances/ driftdist (cm)
    };
    /// One recon element: coefficients for delta_phi (xptar), y_tgt, delta_t (yptar), delta_p;
    /// exponents of xfp, xpfp, yfp, ypfp, fry.
    struct ReconTerm {
        std::array<double, 4> coeff{};
        std::array<int, 5> expon{};
    };

    CosyEngine();

    /// Default forward file of transp_init (relative to the SIMC directory), transp.f:300-314.
    static std::string defaultForwardFile(int spectr);
    /// Default recon file (mc_hms_recon.f / mc_sos_recon.f / mc_shms_recon.f). 1 (HMS), 2 (SOS) and 5 (SHMS) are ported.
    static std::string defaultReconFile(int spectr);

    /// Fortran: transp_init_file(file, n_classes) for spectrometer `spectr`. Returns n_classes.
    int loadForward(int spectr, const std::string& path);
    /// Fortran: transp_init(spectr, n_classes) with the default file below `simcDir`.
    int initForward(int spectr, const std::string& simcDir);
    /// Fortran: the firsttime block of mc_hms_recon / mc_shms_recon. Returns n_terms.
    int loadRecon(int spectr, const std::string& path);
    int initRecon(int spectr, const std::string& simcDir);

    bool hasForward(int spectr) const { return fwd_.count(spectr) != 0; }
    bool hasRecon(int spectr) const { return rec_.count(spectr) != 0; }
    int nClasses(int spectr) const;
    const ForwardClass& forwardClass(int spectr, int cls) const;   ///< cls 1-based
    bool adrift(int spectr, int cls) const { return forwardClass(spectr, cls).adrift; }
    double driftdist(int spectr, int cls) const { return forwardClass(spectr, cls).driftdist; }
    const std::vector<ReconTerm>& reconTerms(int spectr) const;
    int reconMaxOrder(int spectr) const;

    /// Fortran: transp(spectr,class,decay_flag,dflag,m2,ph,zd,pathlen). Updates t (xs, dxdzs, ys, dydzs; dpps
    /// only on decay) and adds zd + delta_z to pathlen. Returns delta_z (= -sum(5)).
    double transp(int spectr, int cls, Track& t, double zd, double& pathlen, DecayState* decay = nullptr);

    /// Fortran: mc_hms_recon / mc_shms_recon(delta_p,delta_t,delta_phi,y_tgt,fry[,spectr]) on the focal-plane
    /// track t (xs, dxdzs, ys, dydzs). Outputs: delta_p [%], delta_t = yptar, delta_phi = xptar [slopes], y_tgt [cm].
    /// FORTRAN-QUIRK (D-015): |fry/100| <= 1e-30 is replaced by 1e-30 (avoid 0**0); for the SOS (spectr 2,
    /// mc_sos_recon.f:120-122) every one of the five hut values.
    void recon(int spectr, const Track& t, double fry, double& delta_p, double& delta_t, double& delta_phi,
               double& y_tgt) const;

    /// Where the Fortran write(6,*) warnings go (default std::cout; nullptr = silent).
    void setLog(std::ostream* os) { log_ = os; }
    /// D-031 run option: true = keep the Fortran bugs of the transport code (transp.f:158 kaon decay mass, D-006;
    /// truncated raddeg of rotate_*axis.f, D-015 F9, used by the chains). Default false (corrected).
    void setFortranCompat(bool b) { fortranCompat_ = b; }
    bool fortranCompat() const { return fortranCompat_; }

private:
    std::map<int, std::vector<ForwardClass>> fwd_;
    std::map<int, std::vector<ReconTerm>> rec_;
    std::map<int, int> recMaxOrder_;
    std::ostream* log_;
    bool fortranCompat_{false};
};

namespace detail {
/// Fortran formatted input of a real with Gw.d editing (blanks ignored, implied decimal point if none).
double readFortranReal(const std::string& field, int d);
/// Fortran I1 input (blank -> 0).
int readFortranInt(const std::string& field);
}  // namespace detail

}  // namespace optics
}  // namespace simc

#endif  // SIMC_OPTICS_COSYENGINE_H
