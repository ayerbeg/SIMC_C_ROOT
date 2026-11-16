/**
 * @file
 * @ingroup optics
 * @brief simc::optics::TargetField: tracking through the polarized-target magnetic field (trg_track.f).
 *
 * Polarized-target field tracking (task P15.2b). Ports trg_track.f (simc_gfortran, GEN/E93-026 target tracking by
 * G. Warren and M. Muehlbauer):
 *   init()          Fortran: trg_track.f:trgInit (l.243-339)            field map + rotation angle of each arm
 *   field()         Fortran: trg_track.f:trgField (l.342-444)           bilinear map interpolation, frame rotation
 *   deriv()         Fortran: trg_track.f:trgDeriv (l.449-488)           dv/dt = factor * (v x B)
 *   rk4()           Fortran: trg_track.f:trgRK4 (l.490-532)             classical RK4 step
 *   trackToPlane()  Fortran: trg_track.f:trgTrackToPlane (l.155-238)    track to the plane a x + b y + c z + d = 0
 *   trackFromTgt()  Fortran: trg_track.f:track_from_tgt (l.591-664)     vertex -> field-free image at z = 100 cm
 *   trackToTgt()    Fortran: trg_track.f:track_to_tgt (l.738-874)       iterative target reconstruction in the field
 * The arm angles of simc.f:120-155 are set up by driver::TargetFieldArmAngles (driver/TargetFieldAngles.h).
 *
 * The routines are deterministic (no random numbers, P15.F2 report §1). Coordinates are the TRANSPORT frame of the arm
 * (x down, y left, z along the spectrometer), cm, cm/ns, T, MeV. The state of the Fortran COMMON blocks is held here:
 *   /trgFieldStrength/ (Bz, Br, zz, rr), /trgFieldAngles_e/, /trgFieldAngles_p/ (TargetField members),
 *   /trgConversionFactor/ factor (factor(): set at the entry of every trackToPlane call with ok = true; T9).
 *
 * FORTRAN-QUIRKs kept (D-067 Categories B/D, docs/translation/markdown/FORTRAN_DIFFERENCES.md §17): pi180 = 3.141592653/180 (T7); field zero
 * for z' >= 100 or r >= 100 cm and for non-finite coordinates (T8: INT() of NaN/Inf/huge = -2147483648); signed dist0 in
 * the step limit (T4); convergence test on the start-of-iteration distance (T5); factor = 90/E (T6); the restart from
 * the wrong-direction step after the direction test (T10); phi treated as 0 (T17).
 */
#ifndef SIMC_OPTICS_TARGETFIELD_H
#define SIMC_OPTICS_TARGETFIELD_H

#include <array>
#include <functional>
#include <string>
#include <vector>

namespace simc {
namespace optics {

/// @addtogroup optics
/// @{

/// @brief The polarized-target field map and the Runge-Kutta tracking between the vertex and z = 0 (trg_track.f).
///
/// Holds the Fortran COMMON blocks and SAVEd locals of trg_track.f. Not thread-safe; one object per worker.
class TargetField {
public:
    static constexpr int nz = 51;  ///< number of z points of the map.
    static constexpr int nr = 51;  ///< number of r points of the map.
    using Vec = std::array<double, 9>;   ///< u(1..9) of the Fortran: x,y,z [cm], dx/dt,dy/dt,dz/dt [cm/ns], (7-9 unused)

    /// Fortran: trgInit(map,theta_e,phi_e,theta_p,phi_p). map "0" = zero field, " " (or empty) = the uniform 5 T test
    /// field (r <= 16 cm, z <= 26 cm), anything else = the file name of the map (7 numbers per line, r outer, z inner;
    /// columns 3 = Bz, 4 = Br). Angles in degrees. Throws std::runtime_error if the file cannot be read.
    /// @param map the map selection (see above).
    /// @param theta_e,phi_e angles of the electron arm [deg].
    /// @param theta_p,phi_p angles of the hadron arm [deg].
    void init(const std::string& map, double theta_e, double phi_e, double theta_p, double phi_p);

    /// Fortran: trgField(x_,B_,spect); spect -1 = electron arm, +1 = hadron arm (any other value: the angles of the
    /// previous call, as the SAVEd locals of the Fortran).
    /// @param x_ position [cm].
    /// @param[out] B_ field [T].
    /// @param spect arm (-1 electron, +1 hadron).
    void field(const double x_[3], double B_[3], int spect);
    /// Fortran: trgDeriv(u,dudt,spect) with the current factor().
    /// @param u the state vector.
    /// @param[out] dudt its time derivative.
    /// @param spect arm (-1 electron, +1 hadron).
    void deriv(const Vec& u, Vec& dudt, int spect);
    /// Fortran: trgRK4(u0,u1,h,spect).
    /// @param u0 the state at the start.
    /// @param[out] u1 the state after the step.
    /// @param h time step [ns].
    /// @param spect arm (-1 electron, +1 hadron).
    void rk4(const Vec& u0, Vec& u1, double h, int spect);
    /// Fortran: trgTrackToPlane(u,E,dl,a,b,c,d,ok,spect). E = energy [MeV] times the charge sign.
    /// @param[in,out] u the state vector.
    /// @param E energy times the charge sign [MeV].
    /// @param dl step length [cm].
    /// @param a,b,c,d the plane a*x + b*y + c*z + d = 0.
    /// @param[out] ok the plane was reached.
    /// @param spect arm (-1 electron, +1 hadron).
    void trackToPlane(Vec& u, double E, double dl, double a, double b, double c, double d, bool& ok, int spect);
    /// Fortran: track_from_tgt(x,y,z,dx,dy,mom,mass,spect,ok). mom [MeV/c] with the charge sign.
    /// @param[in,out] x,y,z position [cm] (vertex in, position at z = 0 out).
    /// @param[in,out] dx,dy slopes.
    /// @param mom momentum with the charge sign [MeV/c].
    /// @param mass particle mass [MeV].
    /// @param spect arm (-1 electron, +1 hadron).
    /// @param[out] ok the track succeeded.
    void trackFromTgt(double& x, double& y, double& z, double& dx, double& dy, double mom, double mass, int spect,
                      bool& ok);

    /// mc_<arm>_recon(delta,dy,dx,y,fry) on the focal-plane track of the last spectrometer call (COMMON /track/).
    using ReconFn = std::function<void(double fry, double& delta, double& dy, double& dx, double& y)>;
    /// Fortran: track_to_tgt(delta,y,dx,dy,frx,fry,mom,mass,ctheta,stheta,spect,ok,arm). `recon` is called for arm
    /// 1..6 (mc_hms/sos/hrsr/hrsl/shms_recon dispatch, trg_track.f:827-837); for any other arm no reconstruction is
    /// made. mom is in/out (INTENT(INOUT) in the Fortran); `n_recon` (optional) returns the number of iterations.
    /// @param[in,out] delta reconstructed delta [%].
    /// @param[in,out] y reconstructed ytar [cm].
    /// @param[in,out] dx,dy reconstructed slopes.
    /// @param frx,fry beam position [cm].
    /// @param[in,out] mom momentum [MeV/c].
    /// @param mass particle mass [MeV].
    /// @param ctheta,stheta cosine and sine of the arm angle.
    /// @param spect arm (-1 electron, +1 hadron).
    /// @param[out] ok the tracking succeeded.
    /// @param arm the arm code.
    /// @param recon the reconstruction of the arm.
    /// @param[out] n_recon number of iterations (optional, may be nullptr).
    void trackToTgt(double& delta, double& y, double& dx, double& dy, double frx, double fry, double& mom, double mass,
                    double ctheta, double stheta, int spect, bool& ok, int arm, const ReconFn& recon,
                    int* n_recon = nullptr);

    /// COMMON /trgConversionFactor/ factor.
    /// @return the conversion factor.
    double factor() const { return factor_; }
    /// @param f the conversion factor.
    void setFactor(double f) { factor_ = f; }

    /// The trgInit state (L1 tests): angles of arm e (spect -1) / p (+1) and the map arrays (1-based iz, ir).
    struct Angles {
        double theta{0} /**< @copydoc cphi */, stheta{0} /**< @copydoc cphi */, ctheta{0} /**< @copydoc cphi */, phi{0} /**< @copydoc cphi */, sphi{0} /**< @copydoc cphi */, cphi{0};  ///< angle [rad] and its sine and cosine, for theta and phi.
    };
    /// @return the electron-arm angles.
    const Angles& anglesE() const { return e_; }
    /// @return the hadron-arm angles.
    const Angles& anglesP() const { return p_; }
    /// @param iz,ir 1-based map indices.
    /// @return Bz [T].
    double Bz(int iz, int ir) const { return Bz_[idx(iz, ir)]; }
    /// @param iz,ir 1-based map indices.
    /// @return Br [T].
    double Br(int iz, int ir) const { return Br_[idx(iz, ir)]; }
    /// @param i 1-based index.
    /// @return the z of the map point [cm].
    double zz(int i) const { return zz_[static_cast<std::size_t>(i - 1)]; }
    /// @param i 1-based index.
    /// @return the r of the map point [cm].
    double rr(int i) const { return rr_[static_cast<std::size_t>(i - 1)]; }

    /// Fortran: trg_track.f:273 PARAMETER pi180 = 3.141592653/180. (truncated pi, T7)
    static constexpr double pi180 = 3.141592653 / 180.;

private:
    static std::size_t idx(int iz, int ir) {   // column-major B_field_z(iz,ir), 1-based
        return static_cast<std::size_t>((ir - 1) * nz + (iz - 1));
    }
    std::vector<double> Bz_ = std::vector<double>(nz * nr, 0.0), Br_ = std::vector<double>(nz * nr, 0.0);
    std::array<double, nz> zz_{};
    std::array<double, nr> rr_{};
    Angles e_, p_;
    Angles last_;          ///< trgField: SAVEd B_tht.. of the previous call (spect not +-1)
    double factor_{0.0};   ///< COMMON /trgConversionFactor/
};

/// @}

}  // namespace optics
}  // namespace simc

#endif  // SIMC_OPTICS_TARGETFIELD_H
