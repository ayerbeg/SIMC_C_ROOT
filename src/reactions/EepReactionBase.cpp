// src/reactions/EepReactionBase.cpp — see include/simc/reactions/EepReactionBase.h.
// Fortran: simc_gfortran/event.f. Operation order follows the Fortran source; x**2 is x*x (D-024), x**3 is (x*x)*x
// (__powidf2).
#include "simc/reactions/EepReactionBase.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "simc/core/Errors.h"
#include "simc/core/Event.h"
#include "simc/core/SimcConstants.h"
#include "simc/material/FortranFormat.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace reactions {

namespace {
inline double sq(double x) { return x * x; }
// gfortran MAX/MIN of two reals: mvar = a; if (b > mvar || isnan(mvar)) mvar = b  (resp. b < mvar)
inline double fmaxF(double a, double b) { return (b > a || std::isnan(a)) ? b : a; }
inline double fminF(double a, double b) { return (b < a || std::isnan(a)) ? b : a; }
}  // namespace

EepReactionBase::EepReactionBase(std::string name) : name_(std::move(name)), log_(&std::cout) {}

void EepReactionBase::configure(const Json& config, const SimcParameters& params) {
    targ_pol_ = params.targ_pol;
    if (config.is_object() && config.contains("verbose") && !config.at("verbose").get<bool>()) log_ = nullptr;
}

void EepReactionBase::init(const StageContext& ctx) {
    s_ = &ctx.setup;
    target_ = &ctx.requireTarget("reaction '" + name_ + "'");
    radiator_ = &ctx.requireRadiator("reaction '" + name_ + "'");
    physics_ = ctx.physics;   // required by completeMain only
    mainCtx_ = MainContext::FromSetup(ctx.setup);
}

const RunSetup& EepReactionBase::setup() const {
    if (!s_) throw std::logic_error("reaction '" + name_ + "': init() was not called");
    return *s_;
}

CompleteEvFn EepReactionBase::completeEvFn() {
    return [this](MainEvent& m, Event& v, IRandom& r) { return completeEv(m, v, r); };
}

// Fortran: event.f:complete_ev (l.432-1054), branches doing_hyd_elast / doing_deuterium / doing_heavy / doing_phsp and the
// common part. (The pion/kaon/delta/rho/semi branches l.611-903 are refused by the plug-ins' init.)
bool EepReactionBase::completeEv(MainEvent& main, Event& vertex, IRandom& rng) {
    using constants::Me;
    using constants::Mp;
    using constants::pi;
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;
    const double Mh2 = s.Mh2;

    main.jacobian = 1.0;   // l.498

    // l.504-511
    vertex.ue.x = std::sin(vertex.e.theta) * std::cos(vertex.e.phi);
    vertex.ue.y = std::sin(vertex.e.theta) * std::sin(vertex.e.phi);
    vertex.ue.z = std::cos(vertex.e.theta);
    if (!f.doing_hyd_elast && !f.doing_rho) {
        vertex.up.x = std::sin(vertex.p.theta) * std::cos(vertex.p.phi);
        vertex.up.y = std::sin(vertex.p.theta) * std::sin(vertex.p.phi);
        vertex.up.z = std::cos(vertex.p.theta);
    }

    // l.517-524
    if (f.doing_hyd_elast) {
        vertex.e.E = vertex.Ein * s.Mh / (s.Mh + vertex.Ein * (1. - vertex.ue.z));
        // FORTRAN-QUIRK(D-024): the only failure; at theta = 0 Ein*Mh/Mh can round 1 ulp above Ein.
        if (vertex.e.E > vertex.Ein) return false;
        vertex.e.P = vertex.e.E;
        vertex.e.delta = (vertex.e.P - s.spec.e.P) * 100. / s.spec.e.P;
    }

    // l.530-539
    vertex.nu = vertex.Ein - vertex.e.E;
    vertex.Q2 = 2 * vertex.Ein * vertex.e.E * (1. - vertex.ue.z);
    vertex.q = std::sqrt(vertex.Q2 + sq(vertex.nu));
    vertex.xbj = vertex.Q2 / 2. / Mp / vertex.nu;
    vertex.uq.x = -vertex.e.P * vertex.ue.x / vertex.q;
    vertex.uq.y = -vertex.e.P * vertex.ue.y / vertex.q;
    vertex.uq.z = (vertex.Ein - vertex.e.P * vertex.ue.z) / vertex.q;
    // FORTRAN-QUIRK(D-024): NaN (theta = 0, q = 0) does not trigger the STOP.
    if (std::fabs(sq(vertex.uq.x) + sq(vertex.uq.y) + sq(vertex.uq.z) - 1) > 0.01)
        throw FortranStop("Error in q vector normalization (event.f:538-539)");

    if (f.doing_hyd_elast) {   // p = q (l.546-563)
        vertex.Em = 0.0;
        vertex.Pm = 0.0;
        vertex.Mrec = 0.0;
        vertex.up.x = vertex.uq.x;
        vertex.up.y = vertex.uq.y;
        vertex.up.z = vertex.uq.z;
        vertex.p.P = vertex.q;
        vertex.p.theta = std::acos(vertex.up.z);
        if (std::fabs(vertex.up.x / std::sin(vertex.p.theta)) > 1 && log_)
            *log_ << " cos(phi)=" << fortran::listReal(vertex.up.x / std::sin(vertex.p.theta)) << "\n";
        vertex.p.phi = std::atan2(vertex.up.y, vertex.up.x);
        if (vertex.p.phi < 0.) vertex.p.phi = vertex.p.phi + 2. * pi;
        SpectrometerAngles(s.spec.p.theta, s.spec.p.phi, vertex.p.xptar, vertex.p.yptar, vertex.p.theta,
                           vertex.p.phi, s.fortran_compat);   // FORTRAN-BUG(D-028 F2) corrected inside unless compat
        vertex.p.E = std::sqrt(sq(vertex.p.P) + Mh2);
        vertex.p.delta = (vertex.p.P - s.spec.p.P) * 100. / s.spec.p.P;
    } else if (f.doing_deuterium) {   // need Ep, and a jacobian (l.565-609)
        vertex.Em = s.targ.Mtar_struck + s.targ.Mrec - s.targ.M;   // =2.2249 MeV
        vertex.Mrec = s.targ.M - s.targ.Mtar_struck + vertex.Em;   // =targ.Mrec
        const double a = -1. * vertex.q * (vertex.uq.x * vertex.up.x + vertex.uq.y * vertex.up.y + vertex.uq.z * vertex.up.z);
        const double b = sq(vertex.q);
        const double c = vertex.nu + s.targ.M;
        const double t = sq(c) - b + Mh2 - sq(vertex.Mrec);
        const double QA = 4. * (sq(a) - sq(c));
        const double QB = 4. * c * t;
        const double QC = -4. * sq(a) * Mh2 - sq(t);
        const double radical = sq(QB) - 4. * QA * QC;
        if (radical < 0) return false;
        vertex.p.E = (-QB - std::sqrt(radical)) / 2. / QA;
        // Check for two solutions, but only print warning if both are within event generation limits (l.588-593).
        const double Ehad2 = (-QB + std::sqrt(radical)) / 2. / QA;
        if (Ehad2 > s.edge.p.E.min && Ehad2 < s.edge.p.E.max && ntried_ <= 5000 && log_) {
            *log_ << " The low-momentum solution to E_hadron is within the spectrometer generation\n"
                  << " limits.  If it is in the acceptance, Its the end of the world as we know it!!!\n"
                  << " E_hadron solns=" << fortran::listReal(vertex.p.E) << fortran::listReal(Ehad2) << "\n";
        }
        if (vertex.p.E <= s.Mh) return false;
        vertex.p.P = std::sqrt(sq(vertex.p.E) - Mh2);
        vertex.p.delta = (vertex.p.P - s.spec.p.P) * 100. / s.spec.p.P;
        // ........ the Jacobian here is |dEp'/dEm|
        main.jacobian = (t * (c - vertex.p.E) + 2 * c * vertex.p.E * (vertex.p.E - c)) /
                        (2 * (sq(a) - sq(c)) * vertex.p.E + c * t);
        main.jacobian = std::fabs(main.jacobian);
    } else if (f.doing_phsp) {   // l.709-716
        vertex.p.P = s.spec.p.P;
        vertex.p.E = std::sqrt(Mh2 + sq(vertex.p.P));
        vertex.p.delta = (vertex.p.P - s.spec.p.P) * 100. / s.spec.p.P;
    }
    // doing_heavy: no p-side block (E_p, P_p, delta were generated); pion/kaon/delta/rho/semi (l.611-903): not ported

    // l.914-926
    vertex.Pmx = vertex.p.P * vertex.up.x - vertex.q * vertex.uq.x;
    vertex.Pmy = vertex.p.P * vertex.up.y - vertex.q * vertex.uq.y;
    vertex.Pmz = vertex.p.P * vertex.up.z - vertex.q * vertex.uq.z;
    vertex.Pmiss = std::sqrt(sq(vertex.Pmx) + sq(vertex.Pmy) + sq(vertex.Pmz));
    vertex.Emiss = vertex.nu + s.targ.M - vertex.p.E;
    const double oop_x = -vertex.uq.y;
    const double oop_y = vertex.uq.x;
    vertex.PmPar = (vertex.Pmx * vertex.uq.x + vertex.Pmy * vertex.uq.y + vertex.Pmz * vertex.uq.z);
    vertex.PmOop = (vertex.Pmx * oop_x + vertex.Pmy * oop_y) / std::sqrt(sq(oop_x) + sq(oop_y));
    // FORTRAN-QUIRK(D-024): for H elastic vertex%Pm is 0 (l.547), so PmPer = 0. For D/A see FORTRAN-BUG(D-033 F5) below.
    vertex.PmPer = std::sqrt(fmaxF(0.e0, sq(vertex.Pm) - sq(vertex.PmPar) - sq(vertex.PmOop)));

    // l.937-955
    if (f.doing_hyd_elast) {
        vertex.Trec = 0.0;
    } else if (f.doing_deuterium) {
        vertex.Pm = vertex.Pmiss;
        vertex.Trec = std::sqrt(sq(vertex.Mrec) + sq(vertex.Pm)) - vertex.Mrec;
    } else if (f.doing_heavy) {
        vertex.Pm = vertex.Pmiss;
        // FORTRAN-BUG(D-033 F7): original Fortran event.f:944: `vertex%Mrec = sqrt(vertex%Emiss**2-vertex%Pmiss**2)`
        // without a guard: Emiss**2 < Pmiss**2 gives NaN Mrec/Em/Trec and the try still succeeds (a NaN Em even passes
        // the VERTEXedge tests of generate_rad, NaN comparisons being false); corrected because such a kinematics has
        // no physical recoil system: the try fails.
        const double mm2 = sq(vertex.Emiss) - sq(vertex.Pmiss);
        if (!s.fortran_compat && !(mm2 >= 0)) return false;
        vertex.Mrec = std::sqrt(mm2);
        vertex.Em = s.targ.Mtar_struck + vertex.Mrec - s.targ.M;
        vertex.Trec = std::sqrt(sq(vertex.Mrec) + sq(vertex.Pm)) - vertex.Mrec;
    }
    if ((f.doing_deuterium || f.doing_heavy) && !s.fortran_compat) {
        // FORTRAN-BUG(D-033 F5): original Fortran event.f:926: `vertex%PmPer = sqrt( max(0.e0, vertex%Pm**2 - ...) )` is
        // computed BEFORE vertex%Pm is assigned (l.940/943), i.e. from the previous try's Pm; corrected because PmPer
        // is the component of this event's missing momentum: recomputed with the current Pm.
        vertex.PmPer = std::sqrt(fmaxF(0.e0, sq(vertex.Pm) - sq(vertex.PmPar) - sq(vertex.PmOop)));
    }

    // l.964-977: krel only for deut/he pion/kaon. FORTRAN-QUIRK(D-033 F8): ntup%krel = the SAVEd local krel, never
    // assigned on the (e,e'p) branches (0).
    main.ntup.krel = krel_;

    // l.1013-1024
    double r = std::sqrt(1. + sq(vertex.e.yptar) + sq(vertex.e.xptar));   // r=cos(theta-theta0)
    main.jacobian = main.jacobian / (r * r * r);   // r**3 (powi: (r*r)*r)
    // FORTRAN-DIFF-CANDIDATE(D-033 F4): event.f:1020-1024 applies the hadron 1/cos**3(theta-theta0) jacobian only for
    // doing_heavy/pion/kaon/delta/semi, NOT for D(e,e'p), although the proton angles are generated in xptar/yptar and
    // genvol contains dOmega_p (simc.f:388). Question for an expert: is the missing jacobian for D intended?
    // Kept as the Fortran in both modes.
    if (f.doing_heavy || f.doing_pion || f.doing_kaon || f.doing_delta || f.doing_semi) {
        r = std::sqrt(1. + sq(vertex.p.yptar) + sq(vertex.p.xptar));
        main.jacobian = main.jacobian / (r * r * r);
    }

    // l.1029-1040
    {
        const ITarget::Traversal te = target_->tripThruTarget(ITarget::kElectron, main.target.z - s.targ.zoffset,
                                                              vertex.e.E, vertex.e.theta, Me, ITarget::kSample, rng);
        main.target.Eloss[1] = te.eloss;
        main.target.teff[1] = te.radlen;
        const ITarget::Traversal tp = target_->tripThruTarget(ITarget::kHadron, main.target.z - s.targ.zoffset,
                                                              vertex.p.E, vertex.p.theta, s.Mh, ITarget::kSample, rng);
        main.target.Eloss[2] = tp.eloss;
        main.target.teff[2] = tp.radlen;
    }
    if (!s.sw.using_Eloss) {
        main.target.Eloss[1] = 0.0;
        main.target.Eloss[2] = 0.0;
    }

    radiator_->initEvent(main, vertex);   // l.1045 radc_init_ev
    return true;
}

// Fortran: event.f:complete_recon_ev (l.1056-1359)
bool EepReactionBase::completeRecon(Event& recon) {
    using constants::Mp;
    using constants::pi;
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;
    double qx, qy, qz, px, py, pz, dummy, new_x_x, new_x_y, new_x_z, new_y_x, new_y_y, new_y_z, p_new_x, p_new_y;

    recon.Ein = s.Ebeam_vertex_ave;                          // l.1087
    recon.Ein = s.Ebeam_vertex_ave - s.targ.Coulomb.ave;     // l.1096 (overwrites)

    // l.1099-1104
    recon.ue.x = std::sin(recon.e.theta) * std::cos(recon.e.phi);
    recon.ue.y = std::sin(recon.e.theta) * std::sin(recon.e.phi);
    recon.ue.z = std::cos(recon.e.theta);
    recon.up.x = std::sin(recon.p.theta) * std::cos(recon.p.phi);
    recon.up.y = std::sin(recon.p.theta) * std::sin(recon.p.phi);
    recon.up.z = std::cos(recon.p.theta);

    // l.1110-1115
    recon.nu = recon.Ein - recon.e.E;
    recon.Q2 = 2 * recon.Ein * recon.e.E * (1 - recon.ue.z);
    recon.q = std::sqrt(recon.Q2 + sq(recon.nu));
    recon.uq.x = -recon.e.P * recon.ue.x / recon.q;
    recon.uq.y = -recon.e.P * recon.ue.y / recon.q;
    recon.uq.z = (recon.Ein - recon.e.P * recon.ue.z) / recon.q;

    // l.1124-1127 (W with the proton mass)
    const double W2 = sq(Mp) + 2. * Mp * recon.nu - recon.Q2;
    recon.W = std::sqrt(std::fabs(W2)) * W2 / std::fabs(W2);
    recon.xbj = recon.Q2 / 2. / Mp / recon.nu;

    // l.1132-1136. FORTRAN-QUIRK(D-024): phsp overwrites the reconstructed hadron momentum.
    if (f.doing_phsp) {
        recon.p.P = s.spec.p.P;
        recon.p.E = std::sqrt(s.Mh2 + sq(recon.p.P));
    }

    // l.1140-1141
    recon.epsilon = 1. / (1. + 2. * (1 + sq(recon.nu) / recon.Q2) * sq(std::tan(recon.e.theta / 2.)));
    recon.theta_pq =
        std::acos(fminF(1.0, recon.up.x * recon.uq.x + recon.up.y * recon.uq.y + recon.up.z * recon.uq.z));

    // l.1160-1186
    qx = -recon.uq.y;
    qy = recon.uq.x;
    qz = recon.uq.z;
    px = -recon.up.y;
    py = recon.up.x;
    pz = recon.up.z;
    dummy = std::sqrt((sq(qx) + sq(qy)) * (sq(qx) + sq(qy) + sq(qz)));
    new_x_x = -qx * qz / dummy;
    new_x_y = -qy * qz / dummy;
    new_x_z = (sq(qx) + sq(qy)) / dummy;
    dummy = std::sqrt(sq(qx) + sq(qy));
    new_y_x = qy / dummy;
    new_y_y = -qx / dummy;
    new_y_z = 0.0;
    p_new_x = px * new_x_x + py * new_x_y + pz * new_x_z;
    p_new_y = px * new_y_x + py * new_y_y + pz * new_y_z;
    if ((sq(p_new_x) + sq(p_new_y)) == 0.)
        recon.phi_pq = 0.0;
    else
        recon.phi_pq = std::acos(p_new_x / std::sqrt(sq(p_new_x) + sq(p_new_y)));
    if (p_new_y < 0.) recon.phi_pq = 2 * pi - recon.phi_pq;

    if (f.using_tgt_field) {   // l.1189-1271 polarized-target angles
        const double targ_Bangle = s.targ_Bangle;
        qx = -recon.uq.y;
        qy = recon.uq.x;
        qz = recon.uq.z;
        double targx = -targ_pol_ * std::sin(std::fabs(targ_Bangle));
        double targy = 0.0;
        double targz = targ_pol_ * std::cos(std::fabs(targ_Bangle));
        dummy = std::sqrt((sq(qx) + sq(qy)) * (sq(qx) + sq(qy) + sq(qz)));
        new_x_x = -qx * qz / dummy;
        new_x_y = -qy * qz / dummy;
        new_x_z = (sq(qx) + sq(qy)) / dummy;
        dummy = std::sqrt(sq(qx) + sq(qy));
        new_y_x = qy / dummy;
        new_y_y = -qx / dummy;
        new_y_z = 0.0;
        p_new_x = targx * new_x_x + targy * new_x_y + targz * new_x_z;
        p_new_y = targx * new_y_x + targy * new_y_y + targz * new_y_z;
        recon.phi_targ = std::atan2(p_new_y, p_new_x);
        if (recon.phi_targ < 0.) recon.phi_targ = 2. * pi + recon.phi_targ;

        qx = -recon.uq.y;
        qy = recon.uq.x;
        qz = recon.uq.z;
        targx = -targ_pol_ * std::sin(std::fabs(targ_Bangle));
        targy = 0.0;
        targz = targ_pol_ * std::cos(std::fabs(targ_Bangle));
        px = -recon.up.y;
        py = recon.up.x;
        pz = recon.up.z;
        dummy = std::sqrt(sq(qy * pz - qz * py) + sq(qz * px - qx * pz) + sq(qx * py - qy * px));
        new_y_x = (qy * pz - qz * py) / dummy;
        new_y_y = (qz * px - qx * pz) / dummy;
        new_y_z = (qx * py - qy * px) / dummy;
        dummy = std::sqrt(sq(new_y_y * qz - new_y_z * qy) + sq(new_y_z * qx - new_y_x * qz) +
                          sq(new_y_x * qy - new_y_y * qx));
        new_x_x = (new_y_y * qz - new_y_z * qy) / dummy;
        new_x_y = (new_y_z * qx - new_y_x * qz) / dummy;
        new_x_z = (new_y_x * qy - new_y_y * qx) / dummy;
        const double targ_new_x = targx * new_x_x + targy * new_x_y + targz * new_x_z;
        const double targ_new_y = targx * new_y_x + targy * new_y_y + targz * new_y_z;
        recon.beta = std::atan2(targ_new_y, targ_new_x);
        if (recon.beta < 0.) recon.beta = 2. * pi + recon.beta;
        recon.phi_s = recon.phi_pq - recon.phi_targ;
        if (recon.phi_s < 0.) recon.phi_s = 2 * pi + recon.phi_s;
        recon.phi_c = recon.phi_pq + recon.phi_targ;
        if (recon.phi_c > 2. * pi) recon.phi_c = recon.phi_c - 2 * pi;
        if (recon.phi_c < 0.) recon.phi_c = 2 * pi + recon.phi_c;
        dummy = std::sqrt((sq(qx) + sq(qy) + sq(qz))) * std::sqrt((sq(targx) + sq(targy) + sq(targz)));
        recon.theta_tarq = std::acos((qx * targx + qy * targy + qz * targz) / dummy);
    }

    // l.1291-1301
    recon.Pmx = recon.p.P * recon.up.x - recon.q * recon.uq.x;
    recon.Pmy = recon.p.P * recon.up.y - recon.q * recon.uq.y;
    recon.Pmz = recon.p.P * recon.up.z - recon.q * recon.uq.z;
    recon.Pm = std::sqrt(sq(recon.Pmx) + sq(recon.Pmy) + sq(recon.Pmz));
    const double oop_x = -recon.uq.y;
    const double oop_y = recon.uq.x;
    recon.PmPar = (recon.Pmx * recon.uq.x + recon.Pmy * recon.uq.y + recon.Pmz * recon.uq.z);
    recon.PmOop = (recon.Pmx * oop_x + recon.Pmy * oop_y) / std::sqrt(sq(oop_x) + sq(oop_y));
    recon.PmPer = std::sqrt(fmaxF(0.e0, sq(recon.Pm) - sq(recon.PmPar) - sq(recon.PmOop)));

    // l.1303-1320
    if (f.doing_pion || f.doing_kaon || f.doing_delta || f.doing_rho || f.doing_semi) {
        recon.Em = recon.nu + s.targ.Mtar_struck - recon.p.E;
        const double mm2 = sq(recon.Em) - sq(recon.Pm);
        const double mm = std::sqrt(std::fabs(mm2)) * std::fabs(mm2) / mm2;
        const double mmA2 = sq(recon.nu + s.targ.M - recon.p.E) - sq(recon.Pm);
        const double mmA = std::sqrt(std::fabs(mmA2)) * std::fabs(mmA2) / mmA2;
        const double t = recon.Q2 - s.Mh2 + 2 * (recon.nu * recon.p.E - recon.p.P * recon.q * std::cos(recon.theta_pq));
        reconNtup_ = {mm, mmA, t};
    }
    if (f.doing_semi || f.doing_rho) {
        recon.zhad = recon.p.E / recon.nu;
        recon.pt2 = sq(recon.p.P) * (1.0 - sq(std::cos(recon.theta_pq)));
    }

    // l.1325-1335
    if (f.doing_hyd_elast) {
        recon.Trec = 0.0;
        recon.Em = recon.nu + s.targ.M - recon.p.E - recon.Trec;
    } else if (f.doing_deuterium || f.doing_heavy) {
        recon.Trec = std::sqrt(sq(recon.Pm) + sq(s.targ.Mrec)) - s.targ.Mrec;
        recon.Em = recon.nu + s.targ.Mtar_struck - recon.p.E - recon.Trec;
    } else if (f.doing_pion || f.doing_kaon || f.doing_delta || f.doing_rho) {
        recon.Em = recon.nu + s.targ.Mtar_struck - recon.p.E;
    }
    return true;
}

// Fortran: event.f:complete_main (l.1363-1569). The reaction's part (main%SF_weight, l.1397-1430) is spectralWeight();
// the rest is reactions::assembleWeight (tgtweight = 1 for (e,e'p)).
bool EepReactionBase::completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0,
                                   const Event& recon) {
    if (!physics_) throw std::logic_error("reaction '" + name_ + "': no physics model (StageContext::physics)");
    main.SF_weight = spectralWeight(vertex);
    return assembleWeight(force_sigcc, main, vertex, vertex0, recon, *physics_, mainCtx_);
}

// Fortran: simc.f:383-396
double EepReactionBase::genVolume() const {
    using constants::pi;
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;
    const double domega_e = (s.gen.e.yptar.max - s.gen.e.yptar.min) * (s.gen.e.xptar.max - s.gen.e.xptar.min);
    double domega_p;
    if (f.doing_rho)
        domega_p = 4. * pi;
    else
        domega_p = (s.gen.p.yptar.max - s.gen.p.yptar.min) * (s.gen.p.xptar.max - s.gen.p.xptar.min);
    double genvol = domega_e;
    if (f.doing_deuterium || f.doing_heavy || f.doing_pion || f.doing_kaon || f.doing_delta || f.doing_rho ||
        f.doing_semi)
        genvol = genvol * domega_p * (s.gen.e.E.max - s.gen.e.E.min);
    if (f.doing_heavy || f.doing_semi) genvol = genvol * (s.gen.p.E.max - s.gen.p.E.min);
    return genvol;
}

}  // namespace reactions
}  // namespace simc
