// src/reactions/EventChain.cpp — see include/simc/reactions/EventChain.h.
#include "simc/reactions/EventChain.h"

#include <cmath>
#include <cstdio>
#include <ostream>
#include <stdexcept>
#include <utility>

#include "simc/core/Errors.h"
#include "simc/core/SimcConstants.h"
#include "simc/setup/KinematicsUtils.h"

namespace simc {
namespace reactions {

namespace {
inline double sq(double x) { return x * x; }
inline double fmaxF(double a, double b) { return (b > a || std::isnan(a)) ? b : a; }
inline double fminF(double a, double b) { return (b < a || std::isnan(a)) ? b : a; }
}  // namespace

bool generateEvent(IBeam& beam, IReaction& reaction, MainEvent& main, Event& vertex, Event& orig, IRandom& rng) {
    beam.generate(main, vertex, rng);                    // event.f:150-212
    return reaction.generate(main, vertex, orig, rng);   // event.f:214-428
}

MainContext MainContext::FromSetup(const RunSetup& s) {
    MainContext c;
    c.flags = s.flags;
    c.using_Coulomb = s.sw.using_Coulomb;
    c.Coulomb_ave = s.targ.Coulomb.ave;
    c.Ebeam = s.Ebeam;
    return c;
}

// Fortran: event.f:complete_main (l.1432-1563), after the SF weight (l.1397-1430, IReaction::completeMain)
bool assembleWeight(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& vertex0, const Event& recon,
                    IPhysicsModel& model, const MainContext& ctx, double tgtweight) {
    (void)vertex0;   // used only by the semi-inclusive model (peepiX)
    const ReactionFlags& f = ctx.flags;
    // ... if we have come up with weight<=0, quit now unless FORCE_SIGCC is on (l.1436). NaN passes, as in the Fortran.
    if (main.SF_weight <= 0 && !force_sigcc) return false;

    // Cross section (l.1444-1529)
    if (f.doing_phsp) {
        main.sigcc = 1.0;
        main.sigcc_recon = 1.0;
    } else if (f.doing_hyd_elast || f.doing_deuterium || f.doing_heavy) {
        main.sigcc = model.sigcc(vertex, main);
        main.sigcc_recon = model.sigccRecon(recon, main);
    } else {
        throw NotPorted("complete_main: pion/kaon/delta/rho/semi cross sections not ported");
    }

    // Coulomb focusing factor (l.1533-1535)
    if (ctx.using_Coulomb) main.sigcc = main.sigcc * sq(1.0 + ctx.Coulomb_ave / ctx.Ebeam);

    // l.1557-1560
    main.weight = main.SF_weight * main.jacobian * main.gen_weight * main.sigcc;
    main.weight = main.weight * tgtweight;
    return true;
}

// Fortran: event.f:update_range (l.95-108)
void updateRange(double val, Range& range) {
    range.lo = fminF(range.lo, val);
    range.hi = fmaxF(range.hi, val);
}

// Fortran: event.f:update_slop_range (l.112-125)
void updateSlopRange(double val, SlopItem& slop) {
    slop.lo = fminF(slop.lo, val);
    slop.hi = fmaxF(slop.hi, val);
}

// P15.4: parallel event loop (EventChain.h).
void mergeRange(Range& into, const Range& w) {
    into.lo = fminF(into.lo, w.lo);
    into.hi = fmaxF(into.hi, w.hi);
}

void mergeSlopRange(SlopItem& into, const SlopItem& w) {
    into.lo = fminF(into.lo, w.lo);
    into.hi = fmaxF(into.hi, w.hi);
}

void mergeContrib(Contrib& a, const Contrib& b) {
    auto arm = [](ArmRange& x, const ArmRange& y) {
        mergeRange(x.delta, y.delta);
        mergeRange(x.yptar, y.yptar);
        mergeRange(x.xptar, y.xptar);
        mergeRange(x.z, y.z);
    };
    arm(a.gen.e, b.gen.e);
    arm(a.gen.p, b.gen.p);
    mergeRange(a.gen.Trec, b.gen.Trec);
    mergeRange(a.gen.sumEgen, b.gen.sumEgen);
    mergeRange(a.tru.e.E, b.tru.e.E);
    mergeRange(a.tru.e.yptar, b.tru.e.yptar);
    mergeRange(a.tru.e.xptar, b.tru.e.xptar);
    mergeRange(a.tru.p.E, b.tru.p.E);
    mergeRange(a.tru.p.yptar, b.tru.p.yptar);
    mergeRange(a.tru.p.xptar, b.tru.p.xptar);
    mergeRange(a.tru.Em, b.tru.Em);
    mergeRange(a.tru.Pm, b.tru.Pm);
    mergeRange(a.tru.Trec, b.tru.Trec);
    arm(a.SP.e, b.SP.e);
    arm(a.SP.p, b.SP.p);
    mergeRange(a.vertex.Trec, b.vertex.Trec);
    mergeRange(a.vertex.Em, b.vertex.Em);
    mergeRange(a.vertex.Pm, b.vertex.Pm);
    for (std::size_t i = 0; i < a.rad.Egamma.size(); ++i) mergeRange(a.rad.Egamma[i], b.rad.Egamma[i]);
    mergeRange(a.rad.Egamma_total, b.rad.Egamma_total);
}

void mergeSlops(Slops& a, const Slops& b) {
    mergeSlopRange(a.total.Em, b.total.Em);
    mergeSlopRange(a.total.Pm, b.total.Pm);
    for (auto [x, y] : {std::pair<Slops::Arm*, const Slops::Arm*>{&a.MC.e, &b.MC.e}, {&a.MC.p, &b.MC.p}}) {
        mergeSlopRange(x->delta, y->delta);
        mergeSlopRange(x->yptar, y->yptar);
        mergeSlopRange(x->xptar, y->xptar);
    }
}

// Fortran: event.f:limits_update (l.1-91)
void limitsUpdate(const MainEvent& main, const Event& vertex, const Event& orig, const Event& recon,
                  const ReactionFlags& f, const std::array<double, 3>& Egamma_used, Contrib& contrib, Slops& slop,
                  bool rho_gen_from_orig) {
    // FORTRAN-BUG(D-057 F18): original Fortran event.f:22-24 update_range(vertex%p%delta/yptar/xptar, contrib%gen%p ...):
    // for doing_rho vertex%p is the rho with its lab theta/phi in the slope fields ('Limiting VERTEX P arm' lines of
    // the .hist); corrected: the decay pion (orig%p).
    const ArmFull& hp = GeneratedHadronArm(rho_gen_from_orig, vertex, orig);
    // GENERATION values
    updateRange(vertex.e.delta, contrib.gen.e.delta);
    updateRange(vertex.e.yptar, contrib.gen.e.yptar);
    updateRange(vertex.e.xptar, contrib.gen.e.xptar);
    updateRange(hp.delta, contrib.gen.p.delta);
    updateRange(hp.yptar, contrib.gen.p.yptar);
    updateRange(hp.xptar, contrib.gen.p.xptar);
    updateRange(main.Trec, contrib.gen.Trec);
    if (f.doing_deuterium || f.doing_pion || f.doing_kaon || f.doing_delta || f.doing_rho)
        updateRange(vertex.e.E - main.Ein_shift, contrib.gen.sumEgen);
    else
        updateRange(vertex.e.E + vertex.p.E - main.Ein_shift, contrib.gen.sumEgen);

    // TRUE values (FORTRAN-QUIRK: orig%e%xptar is updated twice, l.36/38)
    updateRange(orig.e.E - main.Ee_shift, contrib.tru.e.E);
    updateRange(orig.e.xptar, contrib.tru.e.xptar);
    updateRange(orig.e.yptar, contrib.tru.e.yptar);
    updateRange(orig.e.xptar, contrib.tru.e.xptar);
    updateRange(orig.p.E, contrib.tru.p.E);
    updateRange(orig.p.yptar, contrib.tru.p.yptar);
    updateRange(orig.p.xptar, contrib.tru.p.xptar);
    updateRange(orig.Em - main.Ein_shift + main.Ee_shift, contrib.tru.Em);
    updateRange(orig.Pm, contrib.tru.Pm);
    updateRange(orig.Trec, contrib.tru.Trec);

    // SPECTROMETER values
    updateRange(main.SP.e.delta, contrib.SP.e.delta);
    updateRange(main.SP.e.yptar, contrib.SP.e.yptar);
    updateRange(main.SP.e.xptar, contrib.SP.e.xptar);
    updateRange(main.SP.p.delta, contrib.SP.p.delta);
    updateRange(main.SP.p.yptar, contrib.SP.p.yptar);
    updateRange(main.SP.p.xptar, contrib.SP.p.xptar);

    // VERTEX values
    updateRange(vertex.Trec, contrib.vertex.Trec);
    updateRange(vertex.Em, contrib.vertex.Em);
    updateRange(vertex.Pm, contrib.vertex.Pm);

    // RADIATION
    for (int i = 0; i < 3; ++i) updateRange(Egamma_used[i], contrib.rad.Egamma[i]);   // Fortran Egamma_used(i+1)
    updateRange(Egamma_used[0] + Egamma_used[1] + Egamma_used[2], contrib.rad.Egamma_total);

    // MC slops
    updateSlopRange(main.RECON.e.delta - main.SP.e.delta, slop.MC.e.delta);
    updateSlopRange(main.RECON.e.yptar - main.SP.e.yptar, slop.MC.e.yptar);
    updateSlopRange(main.RECON.e.xptar - main.SP.e.xptar, slop.MC.e.xptar);
    updateSlopRange(main.RECON.p.delta - main.SP.p.delta, slop.MC.p.delta);
    updateSlopRange(main.RECON.p.yptar - main.SP.p.yptar, slop.MC.p.yptar);
    updateSlopRange(main.RECON.p.xptar - main.SP.p.xptar, slop.MC.p.xptar);
    // total slops
    updateSlopRange(recon.Em - (orig.Em - main.Ein_shift + main.Ee_shift), slop.total.Em);
    updateSlopRange(std::fabs(recon.Pm) - std::fabs(orig.Pm), slop.total.Pm);
}

// Fortran: jacobians.f:transform_to_cm (l.1-282)
CmTransform transformToCm(const Event& vertex, const MainEvent& main, double pfer, double pferx, double pfery,
                          double pferz, double efer, std::ostream* log) {
    using constants::pi;
    CmTransform o;
    const Event& v = vertex;

    double tcos = v.up.x * v.uq.x + v.up.y * v.uq.y + v.up.z * v.uq.z;
    if (tcos - 1. > 0. && tcos - 1. < 1.e-8) tcos = 1.0;
    const double tsin = std::sqrt(1. - sq(tcos));
    double tfcos = pferx * v.uq.x + pfery * v.uq.y + pferz * v.uq.z;
    if (tfcos - 1. > 0. && tfcos - 1. < 1.e-8) tfcos = 1.0;
    const double tfsin = std::sqrt(1. - sq(tfcos));
    const double cospq = std::cos(main.phi_pq);
    const double sinpq = std::sin(main.phi_pq);

    const double qx = -v.uq.y, qy = v.uq.x, qz = v.uq.z;
    const double px = -pfery, py = pferx, pz = pferz;
    double dummy = std::sqrt((sq(qx) + sq(qy)) * (sq(qx) + sq(qy) + sq(qz)));
    const double tmp_x_x = -qx * qz / dummy;
    const double tmp_x_y = -qy * qz / dummy;
    const double tmp_x_z = (sq(qx) + sq(qy)) / dummy;
    dummy = std::sqrt(sq(qx) + sq(qy));
    const double tmp_y_x = qy / dummy;
    const double tmp_y_y = -qx / dummy;
    const double tmp_y_z = 0.0;
    const double p_tmp_x = pfer * (px * tmp_x_x + py * tmp_x_y + pz * tmp_x_z);
    const double p_tmp_y = pfer * (px * tmp_y_x + py * tmp_y_y + pz * tmp_y_z);
    if (p_tmp_x == 0.)
        o.phiqn = 0.;
    else
        o.phiqn = std::atan2(p_tmp_y, p_tmp_x);
    if (o.phiqn < 0.) o.phiqn = o.phiqn + 2. * pi;

    const double pbeam = v.Ein;
    const double beam_tmpx = pbeam * tmp_x_z;
    const double beam_tmpy = pbeam * tmp_y_z;
    const double beam_tmpz = pbeam * v.uq.z;

    o.bstar = std::sqrt(sq(v.q + pfer * tfcos) + sq(pfer * tfsin)) / (efer + v.nu);
    if (std::fabs(o.bstar) > 1. && log) {
        char buf[128];
        std::snprintf(buf, sizeof buf, "q,nu,efer,beta=%15.5f%15.5f%15.5f%15.5f\n", v.q, v.nu, efer, o.bstar);
        *log << buf;
    }
    o.gstar = 1. / std::sqrt(1. - sq(o.bstar));
    o.bstarz = (v.q + pfer * tfcos) / (efer + v.nu);
    o.bstarx = p_tmp_x / (efer + v.nu);
    o.bstary = p_tmp_y / (efer + v.nu);

    LorenResult L = Loren(o.gstar, o.bstarx, o.bstary, o.bstarz, v.Ein, beam_tmpx, beam_tmpy, beam_tmpz);
    o.ebeamcm = L.ef1; o.pbeamcmx = L.pxf; o.pbeamcmy = L.pyf; o.pbeamcmz = L.pzf; o.pbeamcm = L.pf1;
    const double zero = 0.e0;
    L = Loren(o.gstar, o.bstarx, o.bstary, o.bstarz, v.nu, zero, zero, v.q);
    o.nustar = L.ef1; o.qstarx = L.pxf; o.qstary = L.pyf; o.qstarz = L.pzf; o.qstar = L.pf1;

    const double phadz = v.p.P * tcos;
    const double phadx = v.p.P * tsin * cospq;
    const double phady = v.p.P * tsin * sinpq;
    L = Loren(o.gstar, o.bstarx, o.bstary, o.bstarz, v.p.E, phadx, phady, phadz);
    o.ehadcm = L.ef1; o.phadcmx = L.pxf; o.phadcmy = L.pyf; o.phadcmz = L.pzf; o.phadcm = L.pf1;
    o.thetacm = std::acos((o.phadcmx * o.qstarx + o.phadcmy * o.qstary + o.phadcmz * o.qstarz) / o.phadcm / o.qstar);

    const double ptarz = pfer * tfcos, ptarx = p_tmp_x, ptary = p_tmp_y;
    L = Loren(o.gstar, o.bstarx, o.bstary, o.bstarz, efer, ptarx, ptary, ptarz);
    o.etarcm = L.ef1; o.ptarcmx = L.pxf; o.ptarcmy = L.pyf; o.ptarcmz = L.pzf; o.ptarcm = L.pf1;

    dummy = std::sqrt(sq(o.qstary * o.pbeamcmz - o.qstarz * o.pbeamcmy) + sq(o.qstarz * o.pbeamcmx - o.qstarx * o.pbeamcmz) +
                      sq(o.qstarx * o.pbeamcmy - o.qstary * o.pbeamcmx));
    const double tmp2_y_x = (o.qstary * o.pbeamcmz - o.qstarz * o.pbeamcmy) / dummy;
    const double tmp2_y_y = (o.qstarz * o.pbeamcmx - o.qstarx * o.pbeamcmz) / dummy;
    const double tmp2_y_z = (o.qstarx * o.pbeamcmy - o.qstary * o.pbeamcmx) / dummy;
    dummy = std::sqrt(sq(tmp2_y_y * o.qstarz - tmp2_y_z * o.qstary) + sq(tmp2_y_z * o.qstarx - tmp2_y_x * o.qstarz) +
                      sq(tmp2_y_x * o.qstary - tmp2_y_y * o.qstarx));
    const double tmp2_x_x = (tmp2_y_y * o.qstarz - tmp2_y_z * o.qstary) / dummy;
    const double tmp2_x_y = (tmp2_y_z * o.qstarx - tmp2_y_x * o.qstarz) / dummy;
    const double tmp2_x_z = (tmp2_y_x * o.qstary - tmp2_y_y * o.qstarx) / dummy;
    const double tmp2_z_x = o.qstarx / o.qstar;
    const double tmp2_z_y = o.qstary / o.qstar;
    const double tmp2_z_z = o.qstarz / o.qstar;
    const double phadcm_tmp2x = o.phadcmx * tmp2_x_x + o.phadcmy * tmp2_x_y + o.phadcmz * tmp2_x_z;
    const double phadcm_tmp2y = o.phadcmx * tmp2_y_x + o.phadcmy * tmp2_y_y + o.phadcmz * tmp2_y_z;
    const double phadcm_tmp2z = o.phadcmx * tmp2_z_x + o.phadcmy * tmp2_z_y + o.phadcmz * tmp2_z_z;
    (void)phadcm_tmp2z;
    o.phicm = std::atan2(phadcm_tmp2y, phadcm_tmp2x);
    if (o.phicm < 0.) o.phicm = 2. * pi + o.phicm;

    const double P = v.p.P, E = v.p.E, nu = v.nu, q = v.q;
    const double g = o.gstar, bx = o.bstarx, by = o.bstary, bz = o.bstarz, b2 = sq(o.bstar);
    const double psign = std::cos(o.phiqn) * cospq + std::sin(o.phiqn) * sinpq;
    const double square_root = q + pfer * tfcos - P * tcos;
    const double dp_dcos_num = P + (sq(P) * tcos - psign * pfer * P * tfsin * tcos / tsin) / square_root;
    const double dp_dcos_den = ((nu + efer - E) * P / E + P * sq(tsin) - psign * pfer * tfsin * tsin) / square_root - tcos;
    const double dp_dcos = dp_dcos_num / dp_dcos_den;
    const double dp_dphi_num =
        pfer * P * tsin * tfsin * (std::cos(o.phiqn) * sinpq - std::sin(o.phiqn) * cospq) / square_root;
    const double dp_dphi_den = tcos + (pfer * tsin * tfsin * psign - P * sq(tsin) - (nu + efer - E) * P / E) / square_root;
    const double dp_dphi = dp_dphi_num / dp_dphi_den;
    const double dt_dcos_lab = 2. * (q * P + (q * tcos - nu * P / E) * dp_dcos);
    const double dt_dphi_lab = 2. * (q * tcos - nu * P / E) * dp_dphi;

    const double dpxdphi = P * tsin * (-sinpq + (g - 1.) * bx / b2 * (by * cospq - bx * sinpq)) +
                           ((o.phadcmx + g * bx * E) / P - g * bx * P / E) * dp_dphi;
    const double dpydphi = P * tsin * (cospq + (g - 1.) * by / b2 * (by * cospq - bx * sinpq)) +
                           ((o.phadcmy + g * by * E) / P - g * by * P / E) * dp_dphi;
    const double dpzdphi = P * (g - 1.) / b2 * bz * tsin * (by * cospq - bx * sinpq) +
                           ((o.phadcmz + g * bz * E) / P - g * bz * P / E) * dp_dphi;
    const double dpxdcos = -P * tcos / tsin * (cospq + (g - 1.) * bx / b2 * (bx * cospq + by * sinpq - bz * tsin / tcos)) +
                           ((o.phadcmx + g * bx * E) / P - g * bx * P / E) * dp_dcos;
    const double dpydcos = -P * tcos / tsin * (sinpq + (g - 1.) * by / b2 * (bx * cospq + by * sinpq - bz * tsin / tcos)) +
                           ((o.phadcmy + g * by * E) / P - g * by * P / E) * dp_dcos;
    const double dpzdcos = P * (1. - (g - 1.) / b2 * bz * tcos / tsin * (bx * cospq + by * sinpq - tsin / tcos * bz)) +
                           ((o.phadcmz + g * bz * E) / P - g * bz * P / E) * dp_dcos;

    const double dpxnewdphi = dpxdphi * tmp2_x_x + dpydphi * tmp2_x_y + dpzdphi * tmp2_x_z;
    const double dpynewdphi = dpxdphi * tmp2_y_x + dpydphi * tmp2_y_y + dpzdphi * tmp2_y_z;
    const double dphicmdphi =
        (dpynewdphi * phadcm_tmp2x - phadcm_tmp2y * dpxnewdphi) / (sq(phadcm_tmp2x) + sq(phadcm_tmp2y));
    const double dpxnewdcos = dpxdcos * tmp2_x_x + dpydcos * tmp2_x_y + dpzdcos * tmp2_x_z;
    const double dpynewdcos = dpxdcos * tmp2_y_x + dpydcos * tmp2_y_y + dpzdcos * tmp2_y_z;
    const double dphicmdcos =
        (dpynewdcos * phadcm_tmp2x - phadcm_tmp2y * dpxnewdcos) / (sq(phadcm_tmp2x) + sq(phadcm_tmp2y));

    o.jacobian = std::fabs(dt_dcos_lab * dphicmdphi - dt_dphi_lab * dphicmdcos);
    o.jac_old = 2 * (efer - 2 * pferz * pfer * E / P * tcos) * (q + pferz * pfer) * P /
                    (efer + nu - (q + pferz * pfer) * E / P * tcos) -
                2 * P * pfer;
    return o;
}

}  // namespace reactions
}  // namespace simc
