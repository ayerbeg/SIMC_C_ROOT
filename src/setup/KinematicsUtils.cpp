// src/setup/KinematicsUtils.cpp
//
// Bit-for-bit ports (the Fortran is compiled with -fdefault-real-8: every literal is a double;
// x**2 is x*x; expressions are evaluated left to right as written). No FMA contraction is
// assumed on either side (x86-64 baseline, no -march).
#include "simc/setup/KinematicsUtils.h"

#include <cmath>
#include <iomanip>

namespace simc {

namespace {
constexpr double kPi = 3.141592653589793;   // constants.inc pi
}

// Fortran: event.f:physics_angles (l.1572-1616)
void PhysicsAngles(double theta0, double phi0, double dx, double dy, double& theta, double& phi,
                   std::ostream* log) {
    const double costh = std::cos(theta0);
    const double sinth = std::sin(theta0);
    const double sinph = std::sin(phi0);
    const double cosph = std::cos(phi0);
    const double r = std::sqrt(1.0 + dx * dx + dy * dy);

    if (std::abs(cosph) > 0.0001) {   // phi not at +/- pi/2 (event.f:1596-1599)
        if (log) {
            *log << " theta,phi will be incorrect if phi0 <> pi/2 or 3*pi/2\n";
            *log << std::setprecision(17) << " phi0=   " << phi0 << "      =   " << phi0 * 180 / kPi
                 << "      degrees\n";
        }
    }

    const double tmp = (costh - dy * sinth * sinph) / r;
    if (std::abs(tmp) > 1 && log) *log << " tmp=" << std::setprecision(17) << tmp << "\n";
    theta = std::acos((costh - dy * sinth * sinph) / r);
    if (dx != 0.0) {
        phi = std::atan((dy * costh + sinth * sinph) / dx);   // gives -90 to 90 deg.
        if (phi <= 0) phi = phi + kPi;                        // make 0 to 180 deg.
        if (sinph < 0.) phi = phi + kPi;                      // add pi to phi for HMS
    } else {
        phi = phi0;
    }
}

// Fortran: event.f:spectrometer_angles (l.1618-1647)
// FORTRAN-QUIRK(D-012 F1): dy = sqrt(1/c^2 - 1 - dx^2) is NaN or ~3e-9 when the true dy is 0.
void SpectrometerAngles(double theta0, double phi0, double& dx, double& dy, double theta, double phi) {
    const double x = std::sin(theta) * std::cos(phi);
    const double y = std::sin(theta) * std::sin(phi);
    const double z = std::cos(theta);
    const double x0 = std::sin(theta0) * std::cos(phi0);
    const double y0 = std::sin(theta0) * std::sin(phi0);
    const double z0 = std::cos(theta0);

    const double cos_dtheta = x * x0 + y * y0 + z * z0;
    dx = x / cos_dtheta;
    dy = std::sqrt(1 / (cos_dtheta * cos_dtheta) - 1. - dx * dx);

    const double y_event = y / cos_dtheta;   // projected to plane perp. to spectrometer.
    if (y_event < y0) dy = -dy;
}

// Fortran: loren.f:loren (l.1-26)
LorenResult Loren(double gam, double bx, double by, double bz, double e, double x, double y, double z) {
    LorenResult r;
    const double gam1 = gam * gam / (1. + gam);
    r.ef1 = gam * (e - bx * x - by * y - bz * z);
    r.pxf = (1 + gam1 * (bx * bx)) * x + gam1 * bx * (by * y + bz * z) - gam * bx * e;
    r.pyf = (1 + gam1 * (by * by)) * y + gam1 * by * (bx * x + bz * z) - gam * by * e;
    r.pzf = (1 + gam1 * (bz * bz)) * z + gam1 * bz * (by * y + bx * x) - gam * bz * e;
    r.pf1 = std::sqrt(r.pxf * r.pxf + r.pyf * r.pyf + r.pzf * r.pzf);
    return r;
}

}  // namespace simc
