// src/spectrometers/hms/McHmsColl.cpp
// Fortran: hms/mc_hms_coll.f:mc_hms_coll (l.1-147) — step pions/muons through the 6.3 cm HMS collimator in
// 20 steps, with hadronic absorption (pion_coll_absorb.f), multiple scattering and energy loss wherever the
// ray is outside the octagon opening. Only called by mc_hms when using_HMScoll and 100 < sqrt(m2) < 200 MeV.
#include <cmath>

#include "simc/spectrometers/hms/kernel/McHms.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace hms {

namespace {
// FORTRAN-QUIRK: `sigA = sigA*1.d-27` (pion_coll_absorb.f:76, 82): 1.d-27 is a DOUBLE PRECISION literal, which
// -fdefault-real-8 promotes to real*16, so the product is formed in quad precision with the quad value of 1e-27 and
// rounded to real*8 (objdump of pion_coll_absorb.o: __extenddftf2, __multf3, __trunctfdf2). The double product differs
// by 1 ulp for 7 of the 788 rows of reference/functions/pion_coll_absorb. 1e-27 = 1/(5**27 * 2**27): 5**27 fits in 64
// bits, so the quad denominator is exact and the quotient is the correctly rounded quad 1e-27, as the compiler's.
double timesQuad1em27(double x) {
    __extension__ typedef __float128 quad;
    const quad e27 = static_cast<quad>(1) / (static_cast<quad>(7450580596923828125ULL) * static_cast<quad>(134217728));
    return static_cast<double>(static_cast<quad>(x) * e27);
}
}  // namespace

// Fortran: hms/pion_coll_absorb.f (l.1-92). Reaction cross sections (DG + TH extension > 2 GeV).
// FORTRAN-DIFF-CANDIDATE(D-037 F5): the routine uses its own mpi = 139.56995 (constants.inc Mpi = 139.57018), a
// transmission that jumps from 1 (Tpi <= 85 MeV, sigA = 0) to 0.938 per 0.315 cm at the first node, the unused mate_t,
// and the HMS collimator material (17 g/cm3, A = 171.57) also for the SHMS collimator. Question: intended parameters
// (and should Tpi <= 85 MeV extrapolate)? Kept in both modes.
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
            sigA = timesQuad1em27(sigA);   // convert to cm^2 (FORTRAN-QUIRK above)
        }
    }
    // GH: if pion kinetic energy is greater than T(14)=4159 MeV, then use T(14) parameters
    if (Tpi > T[13]) {
        sigA = sigreac[13] * std::pow(mate_A, qreac[13]);
        sigA = timesQuad1em27(sigA);
    }

    double sum = 0.0;
    const double lambdai = mate_dens * Navagadro * sigA / mate_A;
    sum = sum + thick * lambdai;
    return std::exp(-sum);
}

void traceAbsorb(optics::CollTrace& c, double p, double thick, double trans, double u) {
    ++c.n_absorb_calls;
    c.ab_p_last = p;
    c.ab_thick_last = thick;
    c.ab_trans_last = trans;
    c.ab_u_last = u;
    if (u > trans) c.absorbed = true;
    if (c.calls) c.calls->push_back({1, {p, thick, trans, u, u > trans ? 1. : 0.}});
}

void traceEloss(optics::CollTrace& c, double len, double epart, double mpart, double eloss) {
    ++c.n_eloss_calls;
    c.eloss_sum = c.eloss_sum + eloss;
    if (c.calls) c.calls->push_back({2, {len, epart, mpart, eloss, 0.}});
}

// Fortran: hms/mc_hms_coll.f:mc_hms_coll(m2,p,p_spec,decay_flag,dflag,success,coll_flag,pathlen)
// m2 = d.m2, p = d.ph (mc_hms's local p, modified in place). Returns `success`.
bool McHms::mcHmsColl(optics::DecayState& d, double p_spec, double& pathlen, IRandom& rng) {
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

    // FORTRAN-BUG(D-021): original Fortran hms/mc_hms_coll.f:83-104: `hSTOP_slit_hor = hSTOP_slit_hor + 1` (and _vert,
    // _oct) inside the 20-step loop, for EVERY step outside the opening; corrected because these are per-ray stop
    // counters (.geni "Slit hor/vert/corners"): a ray crossing the collimator material was counted up to 20 times per
    // counter, even when it finally passed. The corrected code counts each kind at most once per ray, using the
    // slit_h_flag/slit_v_flag/slit_o_flag flags that the Fortran sets but never reads. fortran_compat counts per step.
    const bool perStep = fortranCompat();
    bool slit_h_flag = false, slit_v_flag = false, slit_o_flag = false;
    auto count = [&](int code, bool& flag) {
        if (perStep || !flag) cnt_[code] = cnt_[code] + 1;
        flag = true;
    };

    // Now start stepping
    for (int n = 1; n <= nstep; ++n) {
        bool step_flag = false;
        if (std::fabs(t.ys - y_off) > h_step) {
            coll_flag = true;
            step_flag = true;
            count(kSlitHor, slit_h_flag);
            thick_temp_ = step_size;
        }
        if (std::fabs(t.xs - x_off) > v_step) {
            coll_flag = true;
            step_flag = true;
            count(kSlitVert, slit_v_flag);
            thick_temp_ = step_size;
        }
        if (std::fabs(t.xs - x_off) > (-v_step / h_step * std::fabs(t.ys - y_off) + 3 * v_step / 2)) {
            coll_flag = true;
            step_flag = true;
            count(kSlitOct, slit_o_flag);
            thick_temp_ = step_size;
        }
        const double thick = thick_temp_;

        // If we hit the collimator, check for hadronic reaction, do eloss, mult. scat. for first step
        if (step_flag) {
            if (m2 > 12000. && m2 < 20000.) {   // pion - if muon, don't check for hadronic interaction
                const double trans = pionCollAbsorb(p, thick);
                const double rantemp = rng.uniform();
                if (collTrace_) traceAbsorb(*collTrace_, p, thick, trans, rantemp);
                if (rantemp > trans) return false;
            }
            // FORTRAN-DIFF-CANDIDATE(D-037 F3): hms/mc_hms_coll.f:115-118: musc (2 gauss1) and the Landau-sampled enerloss_new
            // (typeflag 1, 1 gauss1) are called whatever ms_flag / mc_smear say, so the collimator scatters and
            // straggles even with mc_smear = 0. Question: is the collimator scattering intended physics (it is the
            // purpose of the stepping) or should it follow ms_flag? Kept in both modes.
            const double coll_radw = thick / coll_radl;
            mat_.musc(m2, p, coll_radw, t.dydzs, t.dxdzs, rng);
            epart = std::sqrt(p * p + m2);
            eloss = mat_.energyLoss(thick, coll_dens, zcoll, acoll, epart, std::sqrt(m2), 1, rng);
            if (collTrace_) traceEloss(*collTrace_, thick, epart, std::sqrt(m2), eloss);
            epart = epart - eloss;
            eloss_tot = eloss_tot + eloss;
            if (epart < std::sqrt(m2)) return false;
            // FORTRAN-BUG(D-037 F4): original Fortran hms/mc_hms_coll.f:121-124: `if(epart.lt.sqrt(m2)) goto 500` is
            // never true, because enerloss_new (enerloss_new.f:75-77) clamps the loss to (E-m)-1e-7: a particle that
            // ranges out in the collimator continues with p = sqrt((m+1e-7)**2-m**2) ~ 5e-3 MeV/c (up to a path of
            // 4.65e9 cm in the Fortran reference tables). Corrected: a clamped loss stops the particle (counted in hSTOP_coll), as
            // the branch intends. fortran_compat keeps the Fortran.
            if (!fortranCompat() && mat_.lastEnergyLossClamped()) {
                if (collTrace_) collTrace_->ranged_out = true;
                return false;
            }
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
