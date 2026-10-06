// src/material/SimcMaterial.cpp — see include/simc/material/SimcMaterial.h.
#include "simc/material/SimcMaterial.h"

#include <cmath>
#include <iostream>
#include <string>

#include "simc/core/Registry.h"
#include "simc/material/FortranFormat.h"
#include "simc/material/MaterialConstants.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace material {

SimcMaterial::SimcMaterial() : log_(&std::cout) {}

// Fortran: enerloss_new.f:enerloss_new (l.1-85). Operation order kept for bitwise results.
double SimcMaterial::energyLoss(double len, double dens, double zeff, double aeff, double epart, double mpart,
                                int typeflag, IRandom& rng) {
    const double me = me_eloss;
    double Eloss;
    const double thick = len * dens;
    const double gamma = epart / mpart;
    const double beta = std::sqrt(1. - 1. / (gamma * gamma));

    double I;   // ionization potential in MeV
    if (zeff == 1) {
        I = 21.8e-06;
    } else {
        I = (16. * std::pow(zeff, 0.9)) * 1.0e-06;
    }

    const double hnup = 28.816e-06 * std::sqrt(dens * zeff / aeff);   // plasma frequency
    const double log10bg = std::log(beta * gamma) / std::log(10.);
    const double CO = std::log(hnup) - std::log(I) + 0.5;

    // DJG density effect correction
    double denscorr;
    if (log10bg < 0.) {
        denscorr = 0.;
    } else if (log10bg < 3.) {
        const double t = 3. - log10bg;
        denscorr = CO + std::log(10.) * log10bg + std::fabs(CO / 27.) * (t * t * t);
    } else if (log10bg < 4.7) {
        denscorr = CO + std::log(10.) * log10bg;
    } else {
        denscorr = CO + std::log(10.) * 4.7;
    }

    if (thick <= 0.) {
        Eloss = 0.;
    } else {
        const double b2 = beta * beta;
        // (Eloss_mp, the old formula, is computed and overwritten in the Fortran: no side effects.)
        const double Eloss_mp_new = 0.1536e-03 * zeff / aeff * thick / b2 *
                                    (std::log(me / (I * I)) + 1.063 + 2. * std::log(gamma * beta) +
                                     std::log(0.1536 * zeff / aeff * thick / b2) - b2 - denscorr);
        const double Eloss_mp = Eloss_mp_new * 1000.;   // convert to MeV
        const double chsi = 0.307075 / 2. * zeff / aeff * thick / b2;
        if (typeflag == 1) {
            x_ = std::fabs(rng.gauss1(kEnerlossNsig));
        } else if (typeflag == 2) {
            x_ = 3;
        } else if (typeflag == 3) {
            x_ = 0.0067;
        } else if (typeflag == 4) {
            x_ = 1;
        }
        // FORTRAN-QUIRK(D-012 F3): any other typeflag keeps x from the previous call.
        double lambda;
        if (x_ > 0.0) {
            lambda = -2.0 * std::log(x_);
        } else {
            lambda = 100000.;   // FORTRAN-QUIRK(D-012 F2): replaces -2 ln 0
        }
        Eloss = lambda * chsi + Eloss_mp;
    }

    lastClamped_ = Eloss > (epart - mpart);
    if (Eloss > (epart - mpart)) {
        Eloss = (epart - mpart) - 0.0000001;
        ++numerr_;
        if (numerr_ <= 10 && log_) {
            *log_ << " Eloss>Total KE; forcing Eloss=KE\n";
            if (numerr_ == 10) *log_ << "      FURTHER ELOSS ERRORS SUPPRESSED\n";
        }
    }
    return Eloss;
}

// Fortran: shared/musc.f:musc (l.1-58)
void SimcMaterial::musc(double m2, double p, double rad_len, double& dth, double& dph, IRandom& rng) {
    constexpr double Es = 13.6;         // MeV
    constexpr double epsilon = 0.088;
    if (rad_len == 0) return;
    if (p < 25. && log_)
        *log_ << " Momentum passed to musc.f should be in MeV, but p= " << fortran::listReal(p) << "\n";

    const double beta = p / std::sqrt(m2 + p * p);
    const double theta_sigma = Es / p / beta * std::sqrt(rad_len) * (1 + epsilon * std::log10(rad_len / (beta * beta)));

    // Compute new trajectory angles (units are rad)
    dth = dth + theta_sigma * rng.gauss1(kMuscNsig);
    dph = dph + theta_sigma * rng.gauss1(kMuscNsig);
}

// Fortran: shared/musc_ext.f:musc_ext (l.1-54). Note the argument order (dph before dth, y before x).
void SimcMaterial::muscExt(double m2, double p, double rad_len, double x_len, double& dph, double& dth, double& y,
                           double& x, IRandom& rng) {
    constexpr double Es = 13.6;
    constexpr double epsilon = 0.088;
    if (rad_len == 0) return;
    if (x_len <= 0 || rad_len < 0) {
        if (log_) {
            *log_ << " x_len or rad_len < 0 in musc_ext.f\n";
            *log_ << " This is bad.  Really bad.  Dont even ask how bad it is.\n";
            *log_ << " Just fix it now.\n";
        }
        throw FortranStop("musc_ext.f: x_len <= 0 or rad_len < 0 (Fortran STOP); x_len=" + std::to_string(x_len) +
                          " rad_len=" + std::to_string(rad_len));
    }
    if (p < 10. && log_)
        *log_ << " Momentum passed to musc_ext.f should be in MeV, but p= " << fortran::listReal(p) << "\n";

    const double beta = p / std::sqrt(m2 + p * p);
    const double theta_sigma = Es / p / beta * std::sqrt(rad_len) * (1 + epsilon * std::log10(rad_len / (beta * beta)));

    // Compute new trajectory angles and displacements (units are rad and cm)
    double g1 = rng.gauss1(kMuscNsig);   // gaussian, truncated at 99 sigma
    double g2 = rng.gauss1(kMuscNsig);
    dth = dth + theta_sigma * g1;
    x = x + theta_sigma * x_len * g2 / std::sqrt(12.) + theta_sigma * x_len * g1 / 2.;

    g1 = rng.gauss1(kMuscNsig);
    g2 = rng.gauss1(kMuscNsig);
    dph = dph + theta_sigma * g1;
    y = y + theta_sigma * x_len * g2 / std::sqrt(12.) + theta_sigma * x_len * g1 / 2.;
}

}  // namespace material
}  // namespace simc

SIMC_REGISTER(simc::IMaterialModel, "simc_material", simc::material::SimcMaterial)
