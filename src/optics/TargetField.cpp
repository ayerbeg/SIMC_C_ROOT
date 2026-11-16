// src/optics/TargetField.cpp — see include/simc/optics/TargetField.h.
// The polarized-target magnetic field (using_tgt_field): the field map (or the uniform test field), Runge-Kutta
// tracking of a charged particle through it, and the iterative reconstruction of the target quantities with the
// field (trg_track.f). Called by MonteCarlo before and after the spectrometer transport.
// Fortran: simc_gfortran/trg_track.f (compiled -fdefault-real-8 -fno-automatic): every literal is REAL*8 (T18), the
// locals are SAVEd. Operation order kept for bitwise results.
#include "simc/optics/TargetField.h"

#include <climits>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>


namespace simc {
namespace optics {

namespace {
constexpr double cc = 29.9792458;   // trg_track.f:619, 759 speed of light [cm/ns]

/// Fortran INT(x) for a REAL*8 -> INTEGER*4 on x86-64 (cvttsd2si): truncation; NaN, +-Inf and |x| >= 2^31 give
/// -2147483648 (T8: such coordinates produce a zero field).
int fint(double x) {
    if (!(x > -2147483649.0 && x < 2147483648.0)) return INT_MIN;
    return static_cast<int>(x);
}
/// Fortran SIGN(1.,x) (gfortran: copysign, so SIGN(1.,-0.) = -1, T10).
double sgn(double x) { return std::copysign(1.0, x); }
}  // namespace

// Fortran: trg_track.f:trgInit (l.243-339)
void TargetField::init(const std::string& map, double theta_e, double phi_e, double theta_p, double phi_p) {
    e_.theta = theta_e;
    e_.stheta = std::sin(theta_e * pi180);
    e_.ctheta = std::cos(theta_e * pi180);
    p_.theta = theta_p;
    p_.stheta = std::sin(theta_p * pi180);
    p_.ctheta = std::cos(theta_p * pi180);
    // Note: for performance reasons B_phi is always treated 0 in trgField (T17)
    e_.phi = phi_e;
    e_.sphi = std::sin(phi_e * pi180);
    e_.cphi = std::cos(phi_e * pi180);
    p_.phi = phi_p;
    p_.sphi = std::sin(phi_p * pi180);
    p_.cphi = std::cos(phi_p * pi180);

    const bool blank = map.find_first_not_of(' ') == std::string::npos;   // Fortran map .EQ. ' '
    std::string trimmed = map;
    while (!trimmed.empty() && trimmed.back() == ' ') trimmed.pop_back();
    if (trimmed == "0") {   // zero field option
        for (int ir = 1; ir <= nr; ++ir) {
            rr_[static_cast<std::size_t>(ir - 1)] = 2. * static_cast<double>(ir - 1);
            zz_[static_cast<std::size_t>(ir - 1)] = 2. * static_cast<double>(ir - 1);
            for (int iz = 1; iz <= nz; ++iz) {
                Br_[idx(iz, ir)] = 0.;
                Bz_[idx(iz, ir)] = 0.;
            }
        }
    } else if (!blank) {   // read in numerical field map: READ (1,*) xx,xx,B_field_z(iz,ir),B_field_r(iz,ir),xx,xx,xx
        std::ifstream in(trimmed);
        if (!in) throw std::runtime_error("trgInit: cannot open the target field map '" + trimmed + "'");
        for (int ir = 1; ir <= nr; ++ir) {
            rr_[static_cast<std::size_t>(ir - 1)] = 2. * static_cast<double>(ir - 1);
            zz_[static_cast<std::size_t>(ir - 1)] = 2. * static_cast<double>(ir - 1);
            for (int iz = 1; iz <= nz; ++iz) {
                std::string line;
                if (!std::getline(in, line))
                    throw std::runtime_error("trgInit: target field map '" + trimmed + "' ends early");
                std::istringstream ls(line);
                std::string tok[7];
                for (auto& t : tok) ls >> t;
                if (tok[3].empty()) throw std::runtime_error("trgInit: bad line in the target field map: " + line);
                Bz_[idx(iz, ir)] = std::strtod(tok[2].c_str(), nullptr);
                Br_[idx(iz, ir)] = std::strtod(tok[3].c_str(), nullptr);
            }
        }
    } else {   // uniform 5T field over 26 cm in z and 16 cm in r
        for (int ir = 1; ir <= nr; ++ir) {
            rr_[static_cast<std::size_t>(ir - 1)] = 2. * static_cast<double>(ir - 1);
            zz_[static_cast<std::size_t>(ir - 1)] = 2. * static_cast<double>(ir - 1);
        }
        for (int ir = 1; ir <= nr; ++ir) {
            for (int iz = 1; iz <= nz; ++iz) {
                Br_[idx(iz, ir)] = 0.;
                if (rr_[static_cast<std::size_t>(ir - 1)] <= 16. && zz_[static_cast<std::size_t>(iz - 1)] <= 26.)
                    Bz_[idx(iz, ir)] = 5.0;
                else
                    Bz_[idx(iz, ir)] = 0.0;
            }
        }
    }
}

// Fortran: trg_track.f:trgField (l.342-444)
void TargetField::field(const double x_[3], double B_[3], int spect) {
    // set up angles for electron or proton spectrometer
    if (spect == -1)
        last_ = e_;
    else if (spect == 1)
        last_ = p_;
    const double B_stht = last_.stheta, B_ctht = last_.ctheta;

    // rotate to coordinates with z' along field direction
    double x[3];
    x[0] = x_[0];
    x[1] = B_stht * x_[2] + B_ctht * x_[1];
    x[2] = B_ctht * x_[2] - B_stht * x_[1];
    // compute zylinder coordinates
    const double z = std::fabs(x[2]);
    const double r = std::sqrt(x[0] * x[0] + x[1] * x[1]);

    // interpolate the field map
    const double dzz = zz_[1] - zz_[0], drr = rr_[1] - rr_[0];
    const int i = fint((z - zz_[0]) / dzz) + 1;
    const int j = fint((r - rr_[0]) / drr) + 1;
    // (i+1 computed in 64 bits: for i = INT_MAX the Fortran INTEGER*4 i+1 is never reached, INT() caps at INT_MIN)
    if ((static_cast<long long>(i) + 1 > nz) || (i < 1) || (static_cast<long long>(j) + 1 > nr) || (j < 1)) {
        B_[0] = 0.;
        B_[1] = 0.;
        B_[2] = 0.;
        return;
    }
    double B[3];
    // calculate the Bz component
    const double az = ((z - zz_[static_cast<std::size_t>(i - 1)]) / dzz);
    const double ar = ((r - rr_[static_cast<std::size_t>(j - 1)]) / drr);
    double a0 = az * (Bz_[idx(i + 1, j)] - Bz_[idx(i, j)]) + Bz_[idx(i, j)];
    double a1 = az * (Bz_[idx(i + 1, j + 1)] - Bz_[idx(i, j + 1)]) + Bz_[idx(i, j + 1)];
    B[2] = (ar * (a1 - a0) + a0);
    if (r > 0.) {
        // calculate the Bx,By components
        a0 = az * (Br_[idx(i + 1, j)] - Br_[idx(i, j)]) + Br_[idx(i, j)];
        a1 = az * (Br_[idx(i + 1, j + 1)] - Br_[idx(i, j + 1)]) + Br_[idx(i, j + 1)];
        B[1] = (ar * (a1 - a0) + a0) / r;
        if (x[2] < 0.) B[1] = -B[1];
        B[0] = B[1] * x[0];
        B[1] = B[1] * x[1];
        // transform B field to lab. system
        B_[0] = B[0];
        B_[1] = -B_stht * B[2] + B_ctht * B[1];
        B_[2] = B_ctht * B[2] + B_stht * B[1];
    } else {
        B_[0] = 0.;
        B_[1] = -B_stht * B[2];
        B_[2] = B_ctht * B[2];
    }
}

// Fortran: trg_track.f:trgDeriv (l.449-488)
void TargetField::deriv(const Vec& u, Vec& dudt, int spect) {
    double B[3];
    field(u.data(), B, spect);
    // These are just the velocities
    dudt[0] = u[3];
    dudt[1] = u[4];
    dudt[2] = u[5];
    // This is just (v_vec X B_vec)
    dudt[6] = u[4] * B[2] - u[5] * B[1];
    dudt[7] = u[5] * B[0] - u[3] * B[2];
    dudt[8] = u[3] * B[1] - u[4] * B[0];
    // This is just (v_vec X B_vec) * factor
    dudt[3] = dudt[6] * factor_;
    dudt[4] = dudt[7] * factor_;
    dudt[5] = dudt[8] * factor_;
}

// Fortran: trg_track.f:trgRK4 (l.490-532). ut(7..9) are never used by trgDeriv.
void TargetField::rk4(const Vec& u0in, Vec& u1, double h, int spect) {
    const Vec u0 = u0in;   // (u1 may alias u0 as in trgTrack; the Fortran reads u0(i) before writing u1(i))
    Vec ut{}, dudt{}, dut{}, dum{};
    const double hh = h * 0.5;
    const double h6 = h / 6.;
    deriv(u0, dudt, spect);
    for (int i = 0; i < 6; ++i) ut[i] = u0[i] + hh * dudt[i];
    deriv(ut, dut, spect);
    for (int i = 0; i < 6; ++i) ut[i] = u0[i] + hh * dut[i];
    deriv(ut, dum, spect);
    for (int i = 0; i < 6; ++i) {
        ut[i] = u0[i] + h * dum[i];
        dum[i] = dut[i] + dum[i];
    }
    deriv(ut, dut, spect);
    for (int i = 0; i < 6; ++i) u1[i] = u0[i] + h6 * (dudt[i] + dut[i] + 2. * dum[i]);
}

// Fortran: trg_track.f:trgTrackToPlane (l.155-238)
void TargetField::trackToPlane(Vec& u, double E, double dl, double a, double b, double c, double d, bool& ok,
                               int spect) {
    if (!ok) return;   // (T9: factor keeps its previous value)
    const double n = 1 / std::sqrt(a * a + b * b + c * c);
    const double an = a * n, bn = b * n, cn = c * n, dn = d * n;
    factor_ = 90. / E;   // FORTRAN-QUIRK(D-067 T6): 90 instead of c^2 * 1e-16 / 1e6 = 89.8755178737
    double ts = -dl / std::sqrt(u[3] * u[3] + u[4] * u[4] + u[5] * u[5]);

    double dist0 = u[0] * an + u[1] * bn + u[2] * cn + dn;
    const double maxdist = std::fmax(std::fabs(dist0) * 4., 1.0);

    // check for the tracking direction
    Vec u0{}, u1{};
    rk4(u, u1, ts, spect);
    double dist1 = u1[0] * an + u1[1] * bn + u1[2] * cn + dn;
    // FORTRAN-QUIRK(D-067 T10): ts is flipped but the loop below starts from the wrong-direction step u1
    if ((sgn(dist0) == sgn(dist1)) && (std::fabs(dist0) < std::fabs(dist1))) ts = -ts;

    // track through the intersection plane
    int steps = 0;
    // FORTRAN-QUIRK(D-067 T4, Category B): the signed dist0 (negative for the -sin(theta) arms: 100 iterations)
    const int max_steps_raw = fint(std::fmax(dist0, 10. * dl) / dl);
    const long long max_steps = static_cast<long long>(max_steps_raw) * 10;   // (INTEGER*4 product; wraps only for INT_MIN)
    const int max_steps32 = static_cast<int>(static_cast<unsigned int>(max_steps & 0xffffffffLL));
    if (sgn(dist0) == sgn(dist1)) {
        dist1 = dist0;
        while ((sgn(dist0) == sgn(dist1)) && ok) {
            rk4(u1, u0, ts, spect);
            dist0 = u0[0] * an + u0[1] * bn + u0[2] * cn + dn;
            if (sgn(dist0) == sgn(dist1)) {
                rk4(u0, u1, ts, spect);
                dist1 = u1[0] * an + u1[1] * bn + u1[2] * cn + dn;
            }
            ok = (std::fabs(dist1) < maxdist) && steps < max_steps32;
            steps = steps + 1;
        }
    } else {
        for (int i = 0; i < 6; ++i) u0[i] = u[i];
    }
    if (ok) {
        // calculate the intersection point: ((u1-u0)*dist0)/(dist0-dist1), left to right as the Fortran
        for (int i = 0; i < 6; ++i) u[i] = u0[i] + (u1[i] - u0[i]) * dist0 / (dist0 - dist1);
    }
}

// Fortran: trg_track.f:track_from_tgt (l.591-664)
void TargetField::trackFromTgt(double& x, double& y, double& z, double& dx, double& dy, double mom, double mass,
                               int spect, bool& ok) {
    const double vel = std::fabs(mom) / std::sqrt(mom * mom + mass * mass) * cc;
    const double eng = sgn(mom) * std::sqrt(mom * mom + mass * mass);
    Vec vT{};
    vT[0] = x;
    vT[1] = y;
    vT[2] = z;
    vT[5] = vel / std::sqrt(1 + dx * dx + dy * dy);
    vT[3] = dx * vT[5];
    vT[4] = dy * vT[5];
    // for debugging, run track first to z=0 (T11: the result is only the start of the next tracking)
    ok = true;
    trackToPlane(vT, eng, 1., 0., 0., 1., 0., ok, spect);
    // track through magnetic field to z=100 cm plane
    ok = true;
    trackToPlane(vT, eng, 1., 0., 0., 1., -100., ok, spect);
    // translate back into SIMC variables
    x = vT[0];
    y = vT[1];
    z = vT[2];
    dx = vT[3] / vT[5];
    dy = vT[4] / vT[5];
}

// Fortran: trg_track.f:track_to_tgt (l.738-874)
void TargetField::trackToTgt(double& delta, double& y, double& dx, double& dy, double frx, double fry, double& mom,
                             double mass, double ctheta, double stheta, int spect, bool& ok, int arm,
                             const ReconFn& recon, int* n_recon) {
    // copy vertical offset into another variable
    double xx = -fry;
    double vel = std::fabs(mom) / std::sqrt(mom * mom + mass * mass) * cc;
    double eng = sgn(mom) * std::sqrt(mom * mom + mass * mass);
    const double mom_0 = mom / (1.e0 + delta / 100.e0);

    Vec vT{}, vTx{};
    vT[0] = -fry + 100. * dx;
    vT[1] = y + 100. * dy;
    vT[2] = 100.;
    vT[5] = vel / std::sqrt(1 + dy * dy + dx * dx);
    vT[3] = dx * vT[5];
    vT[4] = dy * vT[5];
    // and track into the magnetic field to the beam plane (perp. to y)
    trackToPlane(vT, eng, 1., 0., -ctheta, stheta, frx, ok, spect);

    int n = 0;
    double delx = 1.;
    // FORTRAN-QUIRK(D-067 T5, Category B): delx is the distance at the START of the iteration
    while ((delx > .0001) && (n < 10) && ok) {
        delx = std::fabs(-fry - vT[0]);
        // track to the z=0 plane to find a correction for the x-offset
        vTx[0] = -fry;
        for (int i = 1; i < 6; ++i) vTx[i] = vT[i];
        trackToPlane(vT, eng, 1., 0., 0., 1., 0., ok, spect);
        trackToPlane(vTx, eng, 1., 0., 0., 1., 0., ok, spect);
        xx = xx + std::fmin(1., std::fmax(-1., (vTx[0] - vT[0])));
        const double xxd = xx;
        // now find a better approximation: mc_<arm>_recon(delta,dy,dx,y,xxd)
        if (arm >= 1 && arm <= 6 && recon) recon(xxd, delta, dy, dx, y);
        mom = mom_0 * (1.e0 + delta / 100.e0);
        vel = std::fabs(mom) / std::sqrt(mom * mom + mass * mass) * cc;
        eng = sgn(mom) * std::sqrt(mom * mom + mass * mass);
        // drift to a field free region and calculate the velocities
        vT[0] = xx + 100. * dx;
        vT[1] = y + 100. * dy;
        vT[2] = 100.;
        vT[5] = vel / std::sqrt(1 + dy * dy + dx * dx);
        vT[3] = dx * vT[5];
        vT[4] = dy * vT[5];
        // and track into the magnetic field to the beam plane (perp. to y)
        trackToPlane(vT, eng, 1., 0., -ctheta, stheta, frx, ok, spect);
        n = n + 1;
    }
    if (delx > .2) ok = false;
    // calculate the result in HMS coordinates
    dy = vT[4] / vT[5];
    dx = vT[3] / vT[5];
    y = vT[1];
    if (n_recon) *n_recon = n;
}

}  // namespace optics
}  // namespace simc
