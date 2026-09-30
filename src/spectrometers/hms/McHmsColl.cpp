// src/spectrometers/hms/McHmsColl.cpp
// Fortran: hms/mc_hms_coll.f:mc_hms_coll (l.1-147) — step pions/muons through the 6.3 cm HMS collimator in
// 20 steps, with hadronic absorption (pion_coll_absorb.f), multiple scattering and energy loss wherever the
// ray is outside the octagon opening. Only called by mc_hms when using_HMScoll and 100 < sqrt(m2) < 200 MeV.
#include <cmath>

#include "simc/spectrometers/hms/HmsSpectrometer.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace hms {

// Fortran: hms/pion_coll_absorb.f (l.1-95). Reaction cross sections (DG + TH extension > 2 GeV).
double pionCollAbsorb(double ppi, double thick) {
    constexpr double Navagadro = 6.0221367e+23;
    constexpr double mpi = 139.56995;
    static const double T[14] = {85.0,   125.0,   165.0,   205.0,   245.0,   315.0,  584.02,
                                 711.95, 870.12,  1227.57, 1446.58, 1865.29, 2858.0, 4159.0};
    static const double sigreac[14] = {26.03, 84.47, 117.3, 117.4, 101.9, 69.58, 42.5,
                                       44.7,  47.9,  46.5,  45.2,  39.6,  35.34, 33.15};
    static const double qreac[14] = {0.948, 0.659, 0.56,  0.5342, 0.5452, 0.60796, 0.699,
                                     0.689, 0.679, 0.683, 0.688,  0.704,  0.7483,  0.7705};

    const double mate_dens = 17.0;   // Collimator density (g/cm^3)
    // mate_t = 6.35 (total aluminum seen by pion) is set but unused in the Fortran
    const double mate_A = 171.57;

    const double Epi = std::sqrt(ppi * ppi + mpi * mpi);
    const double Tpi = Epi - mpi;

    double sigA = 0.0;
    // GH: loop to 13 (T(i+1) must exist); no exit, at most one interval matches
    for (int i = 0; i < 13; ++i) {
        if ((Tpi > T[i]) && (Tpi <= T[i + 1])) {
            const double Thi = T[i + 1];
            const double Tlo = T[i];
            const double sighi = sigreac[i + 1];
            const double siglo = sigreac[i];
            const double qhi = qreac[i + 1];
            const double qlo = qreac[i];
            const double sigAhi = sighi * std::pow(mate_A, qhi);
            const double sigAlo = siglo * std::pow(mate_A, qlo);
            sigA = (sigAlo * (Thi - Tpi) + sigAhi * (Tpi - Tlo)) / (Thi - Tlo);
            sigA = sigA * 1.e-27;   // convert to cm^2
        }
    }
    // GH: if pion kinetic energy is greater than T(14)=4159 MeV, then use T(14) parameters
    if (Tpi > T[13]) {
        sigA = sigreac[13] * std::pow(mate_A, qreac[13]);
        sigA = sigA * 1.e-27;
    }

    double sum = 0.0;
    const double lambdai = mate_dens * Navagadro * sigA / mate_A;
    sum = sum + thick * lambdai;
    return std::exp(-sum);
}

// Fortran: hms/mc_hms_coll.f:mc_hms_coll(m2,p,p_spec,decay_flag,dflag,success,coll_flag,pathlen)
// m2 = d.m2, p = d.ph (mc_hms's local p, modified in place). Returns `success`.
bool HmsSpectrometer::mcHmsColl(optics::DecayState& d, double p_spec, double& pathlen, IRandom& rng) {
    constexpr double h_entr = 4.575;
    constexpr double v_entr = 11.646;
    constexpr double h_exit = 4.759;
    constexpr double v_exit = 12.114;
    constexpr double x_off = 0.000;
    constexpr double y_off = +0.028;
    constexpr int nstep = 20;
    constexpr double coll_dens = 17.0;
    constexpr double zcoll = 69.45;
    constexpr double acoll = 171.56797;
    constexpr double coll_radl = 0.41753;

    optics::Track& t = track_;
    double& m2 = d.m2;
    double& p = d.ph;

    double eloss_tot = 0.;
    double eloss = 0.;
    const double step_size = 6.30 / nstep;
    bool coll_flag = false;   // (argument; not used by mc_hms)
    double h_step = h_entr;
    double v_step = v_entr;
    double epart = std::sqrt(p * p + m2);
    (void)coll_flag;
    (void)eloss_tot;

    // Now start stepping
    for (int n = 1; n <= nstep; ++n) {
        bool step_flag = false;
        if (std::fabs(t.ys - y_off) > h_step) {
            coll_flag = true;
            step_flag = true;
            cnt_[kSlitHor] = cnt_[kSlitHor] + 1;
            thick_temp_ = step_size;
        }
        if (std::fabs(t.xs - x_off) > v_step) {
            coll_flag = true;
            step_flag = true;
            cnt_[kSlitVert] = cnt_[kSlitVert] + 1;
            thick_temp_ = step_size;
        }
        if (std::fabs(t.xs - x_off) > (-v_step / h_step * std::fabs(t.ys - y_off) + 3 * v_step / 2)) {
            coll_flag = true;
            step_flag = true;
            cnt_[kSlitOct] = cnt_[kSlitOct] + 1;
            thick_temp_ = step_size;
        }
        const double thick = thick_temp_;

        // If we hit the collimator, check for hadronic reaction, do eloss, mult. scat. for first step
        if (step_flag) {
            if (m2 > 12000. && m2 < 20000.) {   // pion - if muon, don't check for hadronic interaction
                const double trans = pionCollAbsorb(p, thick);
                const double rantemp = rng.uniform();
                if (rantemp > trans) return false;
            }
            const double coll_radw = thick / coll_radl;
            mat_.musc(m2, p, coll_radw, t.dydzs, t.dxdzs, rng);
            epart = std::sqrt(p * p + m2);
            eloss = mat_.energyLoss(thick, coll_dens, zcoll, acoll, epart, std::sqrt(m2), 1, rng);
            epart = epart - eloss;
            eloss_tot = eloss_tot + eloss;
            if (epart < std::sqrt(m2)) return false;
            p = std::sqrt(epart * epart - m2);
            t.dpps = 100. * (p / p_spec - 1.);
        }

        optics::project(t, step_size, pathlen, &d);   // drift and decay

        h_step = h_step + (h_exit - h_entr) / nstep;
        v_step = v_step + (v_exit - v_entr) / nstep;
    }
    return true;
}

}  // namespace hms
}  // namespace simc
