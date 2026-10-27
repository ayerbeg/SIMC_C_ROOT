// src/physics/rho/RhoPhysics.cpp — rho_physics.f:peerho, generate_rho.f:generate_rho, rho_decay.f:rho_decay (task P13.1).
// See include/simc/physics/rho/RhoPhysics.h for the numerics and the Fortran findings R1..R7.
#include "simc/physics/rho/RhoPhysics.h"

#include <cmath>
#include <cstdio>
#include <string>

#include "simc/core/Errors.h"
#include "simc/core/SimcConstants.h"
#include "simc/setup/KinematicsUtils.h"

namespace simc {
namespace rho {

namespace {
__extension__ typedef __float128 quad;
inline double sq(double x) { return x * x; }
using constants::pi;

// x/1.d+06 and x/1.d3 with a REAL(16) literal: binary128 division (__divtf3), rounded to REAL*8 on assignment.
inline double divQuad(double x, double d) { return static_cast<double>(static_cast<quad>(x) / static_cast<quad>(d)); }

// Fortran e12.3 / e10.2 / e12.4 edit descriptors for the fort.222 diagnostics (R7): 0.ddddE+xx right-aligned.
std::string fe(double v, int w, int d) {
    char buf[64];
    if (v == 0.0 || !std::isfinite(v)) {
        std::snprintf(buf, sizeof buf, "%*.*E", w, d, v);
        return buf;
    }
    int e = static_cast<int>(std::floor(std::log10(std::fabs(v)))) + 1;
    double m = v / std::pow(10.0, e);
    // rounding may carry to 1.0
    char mb[32];
    std::snprintf(mb, sizeof mb, "%.*f", d, std::fabs(m));
    if (mb[0] == '1') {
        ++e;
        m = v / std::pow(10.0, e);
        std::snprintf(mb, sizeof mb, "%.*f", d, std::fabs(m));
    }
    std::snprintf(buf, sizeof buf, "%s%s%c%c%02d", v < 0 ? "-" : "", mb, 'E', e < 0 ? '-' : '+', std::abs(e));
    std::string s(buf);
    if (static_cast<int>(s.size()) < w) s = std::string(static_cast<std::size_t>(w) - s.size(), ' ') + s;
    return s;
}
}  // namespace

double rhoR(double Q2) { return 0.33 * std::pow(Q2 / constants::Mrho2, 0.61); }

// ============================================================================================ peerho
// Fortran: rho_physics.f:peerho (l.1-402)
double peerho(const Event& vertex, MainEvent& main, const PeerhoConfig& cfg, std::ostream* diag, PeerhoTrace* trace) {
    using constants::Mrho2;
    FermiMomentum& f = main.pfermi;   // COMMON /pfermi/
    const double pfer = f.pfer, pferx = f.pferx, pfery = f.pfery, pferz = f.pferz;
    const double Mts = cfg.Mtar_struck;

    // l.77-81
    const double Q2_g = vertex.Q2 / 1000000.;
    const double phipq = main.phi_pq;
    const double cospq = std::cos(phipq);
    const double sinpq = std::sin(phipq);

    // l.88-93: energy of the struck nucleon
    double efer;
    if (cfg.efer_from_common) {
        // FORTRAN-BUG(P13.1 R1): original Fortran rho_physics.f:88-92: `efer = sqrt(pfer**2+targ%Mtar_struck**2)`,
        // `if(doing_deutpi.or.doing_hepi) then efer = targ%M-sqrt(Mn2+pfer**2) ...` (efer is COMMON /pfermi/, so the
        // on-shell value is also written back); corrected because the comment (l.83-86) asks for "the same assumptions
        // that go into calculating the pion angle/momentum (in event.f)", but the test uses the PION flags, which are
        // never set in a rho run: D/3He rho events were generated with the off-shell efer of event.f:366-371. The
        // COMMON efer of the generation is used (and kept). Identical for H (pfer = 0).
        efer = f.efer;
    } else {
        efer = std::sqrt(sq(pfer) + sq(Mts));
        if (cfg.doing_deutpi || cfg.doing_hepi) {
            efer = cfg.targ_M - std::sqrt(constants::Mn2 + sq(pfer));
            if (cfg.doing_hepi) efer = efer - constants::Mp;
        }
        f.efer = efer;   // COMMON /pfermi/ efer is overwritten (Fortran)
    }

    // l.99-105
    double tcos = vertex.up.x * vertex.uq.x + vertex.up.y * vertex.uq.y + vertex.up.z * vertex.uq.z;
    if (tcos - 1. > 0. && tcos - 1. < 1.e-8) tcos = 1.0;
    const double tsin = std::sqrt(1. - sq(tcos));

    double tfcos = pferx * vertex.uq.x + pfery * vertex.uq.y + pferz * vertex.uq.z;
    if (tfcos - 1. > 0. && tfcos - 1. < 1.e-8) tfcos = 1.0;
    const double tfsin = std::sqrt(1. - sq(tfcos));

    // l.107: `2*(...)`: integer 2
    const double epsi = 1. / (1. + 2. * (1. + sq(vertex.nu) / vertex.Q2) * sq(std::tan(vertex.e.theta / 2.)));

    // l.109-117
    double s = sq(vertex.nu + efer) - sq(vertex.q + pfer * tfcos) - sq(pfer * tfsin);
    const double nu_eq = (s - (sq(efer) - sq(pfer))) / 2. / (efer - pfer * tfcos);   // equiv. photon energy (unused)

    s = s / 1.e6;   // (GeV)**2
    // FORTRAN-QUIRK(P13.1 R5): main%w and main%t are overwritten in GeV units.
    main.w = std::sqrt(s);

    // FORTRAN-DIFF-CANDIDATE(P13.1 R4): t (and tmin, l.233) use the nominal Mrho2, not the Breit-Wigner mass of the
    // event. Question: intended (model defined at the pole mass) or should Mh2 be used as in complete_ev's main%t?
    double t = vertex.Q2 - Mrho2 + 2. * vertex.nu * vertex.p.E - 2. * vertex.p.P * vertex.q * tcos;
    t = t / 1.e6;   // (GeV/c)**2
    main.t = t;

    // l.129-154: pfer in the 'q' system (replay coordinates x = right, y = down)
    const double qx = -vertex.uq.y;
    const double qy = vertex.uq.x;
    const double qz = vertex.uq.z;
    const double px = -pfery;
    const double py = pferx;
    const double pz = pferz;

    double dummy = std::sqrt((sq(qx) + sq(qy)) * (sq(qx) + sq(qy) + sq(qz)));
    double new_x_x = -qx * qz / dummy;
    double new_x_y = -qy * qz / dummy;
    double new_x_z = (sq(qx) + sq(qy)) / dummy;

    dummy = std::sqrt(sq(qx) + sq(qy));
    double new_y_x = qy / dummy;
    double new_y_y = -qx / dummy;
    double new_y_z = 0.0;

    const double p_new_x = pfer * (px * new_x_x + py * new_x_y + pz * new_x_z);
    const double p_new_y = pfer * (px * new_y_x + py * new_y_y + pz * new_y_z);

    double phiqn;
    if (p_new_x == 0.)
        phiqn = 0.;
    else
        phiqn = std::atan2(p_new_y, p_new_x);
    if (phiqn < 0.) phiqn = phiqn + 2. * pi;

    // l.158-161: beam in the 'q' system
    const double pbeam = std::sqrt(sq(vertex.Ein) - constants::Me2);
    const double beam_newx = pbeam * new_x_z;
    const double beam_newy = pbeam * new_y_z;
    const double beam_newz = pbeam * vertex.uq.z;

    // l.163-168
    const double bstar = std::sqrt(sq(vertex.q + pfer * tfcos) + sq(pfer * tfsin)) / (efer + vertex.nu);
    const double gstar = 1. / std::sqrt(1. - sq(bstar));

    const double bstarz = (vertex.q + pfer * tfcos) / (efer + vertex.nu);
    const double bstarx = p_new_x / (efer + vertex.nu);
    const double bstary = p_new_y / (efer + vertex.nu);

    // l.172-174: virtual photon to the CM
    const double zero = 0.e0;
    const LorenResult Lq = Loren(gstar, bstarx, bstary, bstarz, vertex.nu, zero, zero, vertex.q);
    const double qstarx = Lq.pxf, qstary = Lq.pyf, qstarz = Lq.pzf, qstar = Lq.pf1;

    // l.178-184: rho to the CM
    const double ppiz = vertex.p.P * tcos;
    const double ppix = vertex.p.P * tsin * cospq;
    const double ppiy = vertex.p.P * tsin * sinpq;
    const LorenResult Lp = Loren(gstar, bstarx, bstary, bstarz, vertex.p.E, ppix, ppiy, ppiz);
    const double ppicmx = Lp.pxf, ppicmy = Lp.pyf, ppicmz = Lp.pzf, ppicm = Lp.pf1;
    const double thetacm = std::acos((ppicmx * qstarx + ppicmy * qstary + ppicmz * qstarz) / ppicm / qstar);
    main.pcm = ppicm;

    // l.188-189: beam to the CM
    const LorenResult Lb = Loren(gstar, bstarx, bstary, bstarz, vertex.Ein, beam_newx, beam_newy, beam_newz);
    const double pbeamcmx = Lb.pxf, pbeamcmy = Lb.pyf, pbeamcmz = Lb.pzf;

    // l.198-218: y = qcm x pbeamcm, x = y x qcm, z = qcm
    dummy = std::sqrt(sq(qstary * pbeamcmz - qstarz * pbeamcmy) + sq(qstarz * pbeamcmx - qstarx * pbeamcmz) +
                      sq(qstarx * pbeamcmy - qstary * pbeamcmx));
    new_y_x = (qstary * pbeamcmz - qstarz * pbeamcmy) / dummy;
    new_y_y = (qstarz * pbeamcmx - qstarx * pbeamcmz) / dummy;
    new_y_z = (qstarx * pbeamcmy - qstary * pbeamcmx) / dummy;

    dummy = std::sqrt(sq(new_y_y * qstarz - new_y_z * qstary) + sq(new_y_z * qstarx - new_y_x * qstarz) +
                      sq(new_y_x * qstary - new_y_y * qstarx));
    new_x_x = (new_y_y * qstarz - new_y_z * qstary) / dummy;
    new_x_y = (new_y_z * qstarx - new_y_x * qstarz) / dummy;
    new_x_z = (new_y_x * qstary - new_y_y * qstarx) / dummy;

    const double new_z_x = qstarx / qstar;
    const double new_z_y = qstary / qstar;
    const double new_z_z = qstarz / qstar;

    const double ppicm_newx = ppicmx * new_x_x + ppicmy * new_x_y + ppicmz * new_x_z;
    const double ppicm_newy = ppicmx * new_y_x + ppicmy * new_y_y + ppicmz * new_y_z;
    const double ppicm_newz = ppicmx * new_z_x + ppicmy * new_z_y + ppicmz * new_z_z;
    (void)ppicm_newz;   // only in the commented-out dcoscm lines (l.370-371)

    // l.220-228 (sintcm, costcm, sinto2cm, costo2cm are unused)
    double phicm = std::atan2(ppicm_newy, ppicm_newx);
    if (phicm < 0.) phicm = 2. * 3.141592654 + phicm;   // FORTRAN-QUIRK: the literal 3.141592654
    main.thetacm = thetacm;
    main.phicm = phicm;

    // l.233-238: tmin (overall minus sign: t above is -t)
    const double mg = Mts / 1000.;
    const double tmin =
        -(sq((-Q2_g - Mrho2 / 1.e6 - sq(mg) + sq(mg)) / (2. * std::sqrt(s))) - sq((qstar - ppicm) / 1000.));
    const double tprime = std::fabs(t - tmin);

    // l.242: W dependence from photoproduction data
    const double sig0 = 41.263 / std::pow(vertex.nu / 1000.0, 0.4765);   // microbarns

    // l.250-254
    double R = 0.33 * std::pow(vertex.Q2 / Mrho2, 0.61);
    if (R < 0.) {
        if (diag) *diag << "r" << fe(R, 12, 3) << fe(vertex.Q2, 12, 3) << fe(Mrho2, 12, 3) << "\n";
        R = 0.;
    }

    // l.259-264
    double sigt = sig0 * (1.0 + epsi * R) * std::pow(Mrho2 / (vertex.Q2 + Mrho2), 2.575);
    if (sigt < 0.) {
        if (diag)
            *diag << "sigt" << fe(sigt, 12, 3) << fe(sig0, 12, 3) << fe(epsi, 12, 3) << fe(R, 12, 3) << fe(Mrho2, 12, 3)
                  << fe(vertex.Q2, 12, 3) << "\n";
        sigt = 0.;
    }

    // l.269-276: b parameter from the coherence length c*delta tau [fm]
    const double cdeltatau = constants::hbarc / (std::sqrt(sq(vertex.nu) + vertex.Q2 + Mrho2) - vertex.nu);
    double brho;
    if (cdeltatau < 2.0) {
        brho = 4.4679 + 8.6106 * std::log10(cdeltatau);
        if (brho < 1.0) brho = 1.0;
    } else {
        brho = 7.0;
    }

    // l.278-285
    const double sig219 = sigt * brho * std::exp(-brho * tprime) / 2.0 / pi;   // ub/GeV**2/rad
    double sig = divQuad(sig219, 1.e6);   // sig219/1.d+06: REAL(16) literal
    sig = sig * 2. * qstar * ppicm;

    // l.289-297: virtual photon flux
    double gtpr = constants::alpha / 2. / (pi * pi) * vertex.e.E / vertex.Ein * (s - sq(mg)) / 2. /
                  ((efer - pfer * tfcos) / 1000.) / Q2_g / (1. - epsi);
    if (gtpr <= 0.) {
        if (diag)
            *diag << "gtpr" << fe(gtpr, 10, 2) << fe(efer, 10, 2) << fe(pfer, 10, 2) << fe(vertex.Ein, 10, 2)
                  << fe(s, 10, 2) << fe(Mts, 10, 2) << fe(tfcos, 10, 2) << fe(Q2_g, 10, 2) << "\n";
        gtpr = 0.;
    }

    // l.299-300 (s5lab is not used)
    const double factor = Mts / (Mts + vertex.nu - vertex.q * vertex.p.E / vertex.p.P * tcos);

    // l.306-324
    const double psign = std::cos(phiqn) * cospq + std::sin(phiqn) * sinpq;
    const double P = vertex.p.P, E = vertex.p.E;

    const double square_root = vertex.q + pfer * tfcos - P * tcos;
    const double dp_dcos_num = P + (sq(P) * tcos - psign * pfer * P * tfsin * tcos / tsin) / square_root;
    const double dp_dcos_den =
        ((vertex.nu + efer - E) * P / E + P * sq(tsin) - psign * pfer * tfsin * tsin) / square_root - tcos;
    const double dp_dcos = dp_dcos_num / dp_dcos_den;

    const double dp_dphi_num =
        pfer * P * tsin * tfsin * (std::cos(phiqn) * sinpq - std::sin(phiqn) * cospq) / square_root;
    const double dp_dphi_den = tcos + (pfer * tsin * tfsin * psign - P * sq(tsin) - (vertex.nu + efer - E) * P / E) /
                                          square_root;
    const double dp_dphi = dp_dphi_num / dp_dphi_den;

    const double dt_dcos_lab = 2. * (vertex.q * P + (vertex.q * tcos - vertex.nu * P / E) * dp_dcos);
    const double dt_dphi_lab = 2. * (vertex.q * tcos - vertex.nu * P / E) * dp_dphi;

    // l.329-353: dphicm/dphi and dphicm/dcos(theta)
    const double bs2 = sq(bstar);
    const double dpxdphi = P * tsin * (-sinpq + (gstar - 1.) * bstarx / bs2 * (bstary * cospq - bstarx * sinpq)) +
                           ((ppicmx + gstar * bstarx * E) / P - gstar * bstarx * P / E) * dp_dphi;
    const double dpydphi = P * tsin * (cospq + (gstar - 1.) * bstary / bs2 * (bstary * cospq - bstarx * sinpq)) +
                           ((ppicmy + gstar * bstary * E) / P - gstar * bstary * P / E) * dp_dphi;
    const double dpzdphi = P * (gstar - 1.) / bs2 * bstarz * tsin * (bstary * cospq - bstarx * sinpq) +
                           ((ppicmz + gstar * bstarz * E) / P - gstar * bstarz * P / E) * dp_dphi;

    const double dpxdcos =
        -P * tcos / tsin * (cospq + (gstar - 1.) * bstarx / bs2 * (bstarx * cospq + bstary * sinpq - bstarz * tsin / tcos)) +
        ((ppicmx + gstar * bstarx * E) / P - gstar * bstarx * P / E) * dp_dcos;
    const double dpydcos =
        -P * tcos / tsin * (sinpq + (gstar - 1.) * bstary / bs2 * (bstarx * cospq + bstary * sinpq - bstarz * tsin / tcos)) +
        ((ppicmy + gstar * bstary * E) / P - gstar * bstary * P / E) * dp_dcos;
    const double dpzdcos =
        P * (1. - (gstar - 1.) / bs2 * bstarz * tcos / tsin * (bstarx * cospq + bstary * sinpq - tsin / tcos * bstarz)) +
        ((ppicmz + gstar * bstarz * E) / P - gstar * bstarz * P / E) * dp_dcos;

    // l.355-367
    const double dpxnewdphi = dpxdphi * new_x_x + dpydphi * new_x_y + dpzdphi * new_x_z;
    const double dpynewdphi = dpxdphi * new_y_x + dpydphi * new_y_y + dpzdphi * new_y_z;
    const double dphicmdphi = (dpynewdphi * ppicm_newx - ppicm_newy * dpxnewdphi) / (sq(ppicm_newx) + sq(ppicm_newy));

    const double dpxnewdcos = dpxdcos * new_x_x + dpydcos * new_x_y + dpzdcos * new_x_z;
    const double dpynewdcos = dpxdcos * new_y_x + dpydcos * new_y_y + dpzdcos * new_y_z;
    const double dphicmdcos = (dpynewdcos * ppicm_newx - ppicm_newy * dpxnewdcos) / (sq(ppicm_newx) + sq(ppicm_newy));
    // (dpznewdphi, dpznewdcos: only in the commented-out lines l.370-371)

    // l.373-379
    const double jacobian = std::fabs(dt_dcos_lab * dphicmdphi - dt_dphi_lab * dphicmdcos);
    main.davejac = jacobian;
    main.johnjac = 2. * (efer - 2. * pferz * pfer * E / P * tcos) * (vertex.q + pferz * pfer) * P /
                       (efer + vertex.nu - (vertex.q + pferz * pfer) * E / P * tcos) -
                   2. * P * pfer;

    // l.383-395: davesig = gtpr*sig (the jacobian is NOT applied, l.381 is commented out)
    const double davesig = gtpr * sig;
    double sigma_eerho = divQuad(davesig, 1.e3);   // davesig/1.d3: REAL(16) literal; ub/GeV-sr --> ub/MeV-sr
    bool crazy = false;
    if (sigma_eerho > 1.E10 || sigma_eerho < 0.) {
        if (diag) *diag << "crazy sigma_rho" << fe(gtpr, 12, 4) << fe(sig, 12, 4) << fe(sigma_eerho, 12, 4) << "\n";
        sigma_eerho = 0.;
        crazy = true;
    }
    main.ntup.sigcm = sig;   // sig_cm

    if (trace) {
        trace->efer = efer;
        trace->tcos = tcos;
        trace->tsin = tsin;
        trace->tfcos = tfcos;
        trace->tfsin = tfsin;
        trace->epsi = epsi;
        trace->s = s;
        trace->nu_eq = nu_eq;
        trace->t = t;
        trace->tmin = tmin;
        trace->tprime = tprime;
        trace->qstar = qstar;
        trace->ppicm = ppicm;
        trace->thetacm = thetacm;
        trace->phicm = phicm;
        trace->phiqn = phiqn;
        trace->sig0 = sig0;
        trace->R = R;
        trace->sigt = sigt;
        trace->cdeltatau = cdeltatau;
        trace->brho = brho;
        trace->sig219 = sig219;
        trace->sig = sig;
        trace->gtpr = gtpr;
        trace->factor = factor;
        trace->jacobian = jacobian;
        trace->johnjac = main.johnjac;
        trace->davesig = davesig;
        trace->crazy = crazy;
    }
    return sigma_eerho;
}

// ============================================================================================ generate_rho
// Fortran: generate_rho.f:generate_rho (l.1-128)
GenerateRhoResult generateRho(Event& vertex, const FermiMomentum& pf, RhoMassState& gnrl, double& ntup_rhomass,
                              const GenerateRhoConfig& cfg, IRandom& rng) {
    GenerateRhoResult res;
    const double efer = pf.efer, pfer = pf.pfer;

    // l.49-56
    const double nu = vertex.nu;
    const double Q2 = vertex.Q2;
    const double q = vertex.q;
    const double qux = vertex.uq.x;
    const double quy = vertex.uq.y;
    const double quz = vertex.uq.z;

    // l.59-69: velocity of the photon-nucleon CM (minus sign: boost from the CM back to the lab)
    const double bx = -(q * qux + pf.pferx * pfer) / (nu + efer);
    const double by = -(q * quy + pf.pfery * pfer) / (nu + efer);
    const double bz = -(q * quz + pf.pferz * pfer) / (nu + efer);

    const double betacm = std::sqrt(sq(bx) + sq(by) + sq(bz));
    if (betacm > 1.0) {
        res.success = false;
        return res;
    }
    const double gammacm = 1. / std::sqrt(1.0 - sq(betacm));

    // l.72. FORTRAN-DIFF-CANDIDATE(P13.1 R3): "W2 in cm - including possible Fermi motion effects", but the
    // -2 q.pfer - pfer**2 terms of (nu+efer)**2 - |q+pfer|**2 are missing (D/3He only; the boost above includes them).
    // Question: intended approximation? Kept.
    const double s = -Q2 + sq(efer) + 2. * nu * efer;

    // l.75-80: non-relativistic Breit-Wigner rho mass (Gamma = 150.2 MeV, cut at +-500 MeV)
    gnrl.Mh = constants::Mrho;
    gnrl.Mh = gnrl.Mh + 0.5 * 150.2 * std::tan((2. * rng.uniform() - 1.) * std::atan(2. * 500. / 150.2));
    res.n_uniform = 1;
    gnrl.Mh2 = gnrl.Mh * gnrl.Mh;
    ntup_rhomass = gnrl.Mh;
    if (cfg.mh2_final_fixed) {
        // FORTRAN-BUG(P13.1 R2): original Fortran generate_rho.f:80: `Mh2_final = Mh`; corrected because Mh2_final is
        // a mass squared (dbase.f:222, rho_decay.f:282). Not observable in a run (simc.f:174 / rho_decay reset it).
        gnrl.Mh2_final = gnrl.Mh2;
    } else {
        gnrl.Mh2_final = gnrl.Mh;
    }

    // l.82-90
    const double Erhocm = (s + gnrl.Mh2 - sq(cfg.Mrec_struck)) / 2. / std::sqrt(s);
    if (Erhocm < gnrl.Mh) {
        res.success = false;
        return res;
    }
    const double Prhocm = std::sqrt(sq(Erhocm) - gnrl.Mh2);

    // l.94-101: flat in phi and cos(theta) in the CM
    const double rph = rng.uniform() * 2. * pi;
    const double rth1 = rng.uniform() * 2. - 1.;
    res.n_uniform = 3;
    const double rth = std::acos(rth1);

    const double pxr = Prhocm * std::sin(rth) * std::cos(rph);
    const double pyr = Prhocm * std::sin(rth) * std::sin(rph);
    const double pzr = Prhocm * std::cos(rth);
    const double er = Erhocm;

    // l.109: boost to the lab
    const LorenResult L = Loren(gammacm, bx, by, bz, er, pxr, pyr, pzr);
    const double ef = L.ef1, pxf = L.pxf, pyf = L.pyf, pzf = L.pzf, pfm = L.pf1;

    // l.112-124
    vertex.p.delta = 100. * (pfm / cfg.spec_p_P - 1.);
    vertex.p.P = pfm;
    vertex.p.E = ef;
    vertex.up.x = pxf / pfm;
    vertex.up.y = pyf / pfm;
    vertex.up.z = pzf / pfm;
    // "These aren't really used for anything": lab theta and phi
    vertex.p.xptar = std::acos(pzf / pfm);
    vertex.p.yptar = std::atan2(pyf, pxf);

    res.success = true;
    return res;
}

// ============================================================================================ rho_decay
// Fortran: rho_decay.f:rho_decay (l.132-291)
RhoDecayResult rhoDecay(Event& orig, double p_spec, double epsilon, bool success_in, RhoMassState& gnrl,
                        NtupVars& ntup, const RhoDecayConfig& cfg, IRandom& rng, std::ostream* log) {
    RhoDecayResult res;
    res.success = success_in;   // FORTRAN-QUIRK(P13.1 R6): never set to .true. by rho_decay

    // l.179: R = sigma_L/sigma_T (the same parameterisation as rho_physics.f:250)
    const double R_rho = 0.33 * std::pow(orig.Q2 / constants::Mrho2, 0.61);

    // l.182-185 (dlen = ctaurho*beta*gamma is not used)
    const double ph = orig.p.P;
    const double beta = ph / std::sqrt(sq(ph) + gnrl.Mh2);
    const double gamma = 1. / std::sqrt(1. - beta * beta);

    // l.188-195: phi* flat; cos(theta*) by accept-reject from sin**2 + 2 eps R cos**2.
    // FORTRAN-DIFF-CANDIDATE(P13.1 R9): theta*, phi* are taken w.r.t. the LAB z and x axes in the rho rest frame (the
    // pion is boosted with -beta*up below), not w.r.t. the helicity axis (rho direction in the photon-nucleon CM) in
    // which sin**2 + 2 eps R cos**2 is the s-channel-helicity-conserving decay distribution. Question: intended
    // approximation? Kept.
    const double rph = rng.uniform() * 2. * pi;
    res.n_uniform = 1;
    double rth;
    for (;;) {
        const double rth1 = rng.uniform() * 2. - 1.;
        rth = std::acos(rth1);
        const double norm = (1.0 + 2.0 * epsilon * R_rho) * rng.uniform();
        res.n_uniform += 2;
        ++res.n_trials;
        const double dist = sq(std::sin(rth)) + 2.0 * epsilon * R_rho * sq(std::cos(rth));
        if (!(dist < norm)) break;   // `if(dist.lt.norm) goto 100`
    }
    ntup.rhotheta = rth;

    // l.200-208: the pion in the rho rest frame
    const double er = ntup.rhomass / 2.0;
    if (er < constants::Mpi) {
        res.success = false;
        return res;
    }
    const double pr = std::sqrt(sq(er) - constants::Mpi2);
    const double pxr = pr * std::sin(rth) * std::cos(rph);
    const double pyr = pr * std::sin(rth) * std::sin(rph);
    const double pzr = pr * std::cos(rth);

    // l.212-215: boost to the lab
    const double bx = -beta * orig.up.x;
    const double by = -beta * orig.up.y;
    const double bz = -beta * orig.up.z;
    const LorenResult L = Loren(gamma, bx, by, bz, er, pxr, pyr, pzr);
    const double pxf = L.pxf, pyf = L.pyf, pzf = L.pzf, pf = L.pf1;

    // l.223-237: is the pion going toward the spectrometer? (rotate by theta_spec about x)
    const int arm = cfg.hadron_arm;
    const bool right = arm == 1 || arm == 3;
    const bool left = arm == 2 || arm == 4 || arm == 5;
    const double th = cfg.spec_p_theta;
    double pzprime;
    if (right)
        pzprime = pzf * std::cos(th) - pyf * std::sin(th);
    else if (left)
        pzprime = pzf * std::cos(th) + pyf * std::sin(th);
    else
        throw FortranStop("rho_decay.f:228-229: Unknown spectrometer setup dude! (hadron_arm = " + std::to_string(arm) +
                          ")");
    if (pzprime < 0.0) {
        res.success = false;
        return res;
    }

    // l.241-259: xptar and yptar.
    // FORTRAN-DIFF-CANDIDATE(P13.1 R8): xptar = tan(th_oop) = px/sqrt(py**2+pz**2) is dx/dz' * cos(in-plane offset)
    // (the exact spectrometer slope is px/pz'); -0.1 % at 45 mr. Question: approximation accepted? Kept.
    const double th_oop = std::asin(pxf / pf);
    const double cos_th_inp = pzf / pf / std::cos(th_oop);
    double th_inp = std::acos(cos_th_inp);
    if (right) {
        if (pyf < 0.0)
            th_inp = th - th_inp;
        else
            th_inp = th + th_inp;
    } else {
        if (pyf > 0.0)
            th_inp = th_inp - th;
        else
            th_inp = th_inp + th;
    }

    // l.263-267 ("This should never happen!!!")
    if ((th_inp > pi / 2.) || (th_oop > pi / 2.)) {
        if (log) *log << " Decay pion going backwards - fahgettaboutit!\n";
        res.success = false;
        return res;
    }

    // l.269-288
    orig.up.x = pxf / pf;
    orig.up.y = pyf / pf;
    orig.up.z = pzf / pf;
    orig.p.xptar = std::tan(th_oop);
    orig.p.yptar = std::tan(th_inp);
    orig.p.delta = 100. * (pf / p_spec - 1.);
    orig.p.P = pf;

    gnrl.Mh = constants::Mpi;
    gnrl.Mh2 = constants::Mpi2;
    gnrl.Mh2_final = gnrl.Mh2;

    orig.p.E = std::sqrt(sq(orig.p.P) + gnrl.Mh2);

    PhysicsAngles(cfg.spec_p_theta, cfg.spec_p_phi, orig.p.xptar, orig.p.yptar, orig.p.theta, orig.p.phi, log);
    return res;
}

}  // namespace rho
}  // namespace simc
