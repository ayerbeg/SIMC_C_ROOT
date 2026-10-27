// tests/l1/test_l1_rho.cpp — task P13.1: the exclusive rho0 kernels (rho_physics.f:peerho, generate_rho.f, rho_decay.f)
// and the IPhysicsModel "rho".
//
// Modes:
//   sanity                     kinematic identities, RNG order and counts, failure paths, the R1/R2 switches, the
//                              decay frame, the plug-in (no Fortran tables needed)
//   fortran_diff               P13.1 R1/R2 in both modes
//   setup <pion .inp>          the rho set-up from a pion H input switched to doing_rho: flags, Mh, doing_decay accepted,
//                              rho + pion refused, the default pipeline (physics "rho"), checkSupported / model init
// Bitwise L1 against the P13.F1 tables (reference/functions/rho/) is added when the tables exist.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include "simc/core/Errors.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SimcParameters.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/driver/Pipeline.h"
#include "simc/input/InpReader.h"
#include "simc/physics/rho/RhoModel.h"
#include "simc/physics/rho/RhoPhysics.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/setup/SetupBuilder.h"
#include "simc/stages/IPhysicsModel.h"

using namespace simc;
namespace K = simc::constants;

namespace {
int nfail = 0;
void check(bool ok, const std::string& what) {
    std::printf("  [%s] %s\n", ok ? " ok " : "FAIL", what.c_str());
    if (!ok) ++nfail;
}
bool near(double a, double b, double rel) { return std::fabs(a - b) <= rel * std::max(std::fabs(a), std::fabs(b)); }
double sq(double x) { return x * x; }

// An electron-scattering vertex: beam along z, scattered electron in the HMS (beam right: -y) at theta_e.
Event electronVertex(double Ein, double Ee, double theta_e) {
    Event v;
    v.Ein = Ein;
    v.e.E = Ee;
    v.e.P = Ee;
    v.e.theta = theta_e;
    v.ue = {0., -std::sin(theta_e), std::cos(theta_e)};
    const double qx = -Ee * v.ue.x, qy = -Ee * v.ue.y, qz = Ein - Ee * v.ue.z;
    v.nu = Ein - Ee;
    v.q = std::sqrt(qx * qx + qy * qy + qz * qz);
    v.Q2 = v.q * v.q - v.nu * v.nu;
    v.uq = {qx / v.q, qy / v.q, qz / v.q};
    return v;
}

// complete_ev's phi_pq (event.f:743-790), enough for the test: angle of up around uq in the replay frame
double phiPq(const Event& v) {
    const double qx = -v.uq.y, qy = v.uq.x, qz = v.uq.z;
    const double px = -v.up.y, py = v.up.x, pz = v.up.z;
    const double dummy = std::sqrt((qx * qx + qy * qy) * (qx * qx + qy * qy + qz * qz));
    const double nxx = -qx * qz / dummy, nxy = -qy * qz / dummy, nxz = (qx * qx + qy * qy) / dummy;
    const double d2 = std::sqrt(qx * qx + qy * qy);
    const double nyx = qy / d2, nyy = -qx / d2;
    const double pnx = px * nxx + py * nxy + pz * nxz, pny = px * nyx + py * nyy;
    double phi = std::atan2(pny, pnx);
    if (phi < 0) phi += 2 * K::pi;
    return phi;
}

int sanity() {
    std::printf("rho sanity (P13.1)\n");
    const Event v0 = electronVertex(10600., 6000., 12.0 / K::degrad);
    FermiMomentum pfH;
    pfH.efer = K::Mp;   // event.f:335 efer = targ%Mtar_struck, pfer = 0
    rho::GenerateRhoConfig gc;
    gc.Mrec_struck = K::Mp;
    gc.spec_p_P = 4000.;

    // ---- generate_rho at the pole mass (u1 = 0.5 -> tan(0) = 0)
    {
        Event v = v0;
        rho::RhoMassState g{K::Mrho, K::Mrho2, K::Mrho2};
        double rm = 0;
        ScriptedRandom rng;
        rng.load({}, {0.5, 0.3, 0.7, 0.123});
        const auto r = rho::generateRho(v, pfH, g, rm, gc, rng);
        check(r.success && r.n_uniform == 3 && rng.nUniform() == 3, "generate_rho: success, 3 grnd (BW, phi, cos)");
        check(g.Mh == K::Mrho && rm == K::Mrho && g.Mh2 == K::Mrho * K::Mrho && g.Mh2_final == g.Mh2,
              "generate_rho: u1 = 0.5 gives exactly Mrho; Mh2, Mh2_final (corrected R2), ntup%rhomass");
        check(near(sq(v.p.E) - sq(v.p.P), g.Mh2, 1e-9), "generate_rho: rho on shell (E^2 - P^2 = Mh^2)");
        const double rx = -v.p.P * v.up.x + v.q * v.uq.x, ry = -v.p.P * v.up.y + v.q * v.uq.y,
                     rz = -v.p.P * v.up.z + v.q * v.uq.z;
        const double Erec = v.nu + K::Mp - v.p.E;
        check(near(sq(Erec) - (rx * rx + ry * ry + rz * rz), K::Mp2, 1e-9),
              "generate_rho: recoil proton on shell (energy-momentum conservation)");
        check(v.p.delta == 100. * (v.p.P / gc.spec_p_P - 1.), "generate_rho: delta");
        check(v.p.xptar == std::acos(v.up.z * v.p.P / v.p.P) || near(std::cos(v.p.xptar), v.up.z, 1e-12),
              "generate_rho: xptar = lab theta");
        // the CM angles: generated in the CM with the LAB axes, cos(theta*) = 2*0.7-1 w.r.t. the lab z axis
        const double gb = (v.nu + K::Mp), bcx = v.q * v.uq.x / gb, bcy = v.q * v.uq.y / gb, bcz = v.q * v.uq.z / gb;
        const double gcm = 1. / std::sqrt(1. - (bcx * bcx + bcy * bcy + bcz * bcz));
        const LorenResult cm = Loren(gcm, bcx, bcy, bcz, v.p.E, v.p.P * v.up.x, v.p.P * v.up.y, v.p.P * v.up.z);
        check(near(cm.pzf / cm.pf1, 0.4, 1e-9) && near(std::atan2(cm.pyf, cm.pxf) + 2 * K::pi, 2 * K::pi * 0.3 + 2 * K::pi, 1e-9),
              "generate_rho: in the CM (lab axes) cos(theta*) = 2 u3 - 1, phi* = 2 pi u2");
    }

    // ---- Breit-Wigner limits and Mh2_final (R2) in compat
    {
        Event v = v0;
        rho::RhoMassState g{};
        double rm = 0;
        ScriptedRandom rng;
        rng.load({}, {0.0, 0.3, 0.7});
        rho::GenerateRhoConfig gcc = gc;
        gcc.mh2_final_fixed = false;
        rho::generateRho(v, pfH, g, rm, gcc, rng);
        check(near(g.Mh, K::Mrho - 500., 1e-12), "generate_rho: u1 = 0 gives Mrho - 500 MeV (BW cut)");
        check(g.Mh2_final == g.Mh, "generate_rho: fortran_compat keeps Mh2_final = Mh (R2)");
    }

    // ---- below threshold: 1 grnd, the mass is written anyway
    {
        Event v = electronVertex(2000., 1500., 20.0 / K::degrad);   // W ~ 1.25 GeV < Mp + Mrho
        rho::RhoMassState g{K::Mpi, K::Mpi2, K::Mpi2};
        double rm = -1;
        ScriptedRandom rng;
        rng.load({}, {0.5, 0.3, 0.7});
        const auto r = rho::generateRho(v, pfH, g, rm, gc, rng);
        check(!r.success && r.n_uniform == 1 && rng.nUniform() == 1 && rm == K::Mrho && g.Mh == K::Mrho,
              "generate_rho: below threshold -> fail after 1 grnd, Mh/ntup%rhomass already written");
    }

    // ---- |beta_cm| > 1: no grnd
    {
        Event v = v0;
        v.q = 1.0e5;   // unphysical: q >> nu + efer
        rho::RhoMassState g{K::Mrho, K::Mrho2, K::Mrho2};
        double rm = -1;
        ScriptedRandom rng;
        rng.load({}, {0.5});
        const auto r = rho::generateRho(v, pfH, g, rm, gc, rng);
        check(!r.success && r.n_uniform == 0 && rng.nUniform() == 0 && rm == -1 && g.Mh == K::Mrho,
              "generate_rho: beta_cm > 1 -> fail before any grnd, nothing written");
    }

    // ---- peerho on the generated H event
    Event vr = v0;
    rho::RhoMassState gr{K::Mrho, K::Mrho2, K::Mrho2};
    double rhomass = 0;
    {
        ScriptedRandom rng;
        rng.load({}, {0.55, 0.3, 0.8});
        rho::generateRho(vr, pfH, gr, rhomass, gc, rng);
    }
    MainEvent m0;
    m0.pfermi = pfH;
    m0.phi_pq = phiPq(vr);
    rho::PeerhoConfig pc;
    pc.Mtar_struck = K::Mp;
    pc.targ_M = K::Mp;
    {
        MainEvent m = m0;
        rho::PeerhoTrace tr;
        std::ostringstream diag;
        const double sig = rho::peerho(vr, m, pc, &diag, &tr);
        check(std::isfinite(sig) && sig > 0 && !tr.crazy && diag.str().empty(), "peerho: finite, positive, no fort.222 line");
        const double W2 = K::Mp2 + 2 * K::Mp * vr.nu - vr.Q2;
        check(near(m.w * m.w * 1e6, W2, 1e-12), "peerho: main%w = sqrt(s) in GeV (R5)");
        check(m.t == (vr.Q2 - K::Mrho2 + 2. * vr.nu * vr.p.E - 2. * vr.p.P * vr.q * tr.tcos) / 1.e6,
              "peerho: main%t in GeV^2 with the nominal Mrho2 (R4)");
        check(m.ntup.sigcm == tr.sig && m.davejac == tr.jacobian && m.pcm == tr.ppicm && m.thetacm == tr.thetacm &&
                  m.phicm == tr.phicm,
              "peerho: main%davejac, pcm, thetacm, phicm, ntup%sigcm written");
        check(m.phicm >= 0 && m.phicm < 2 * 3.141592654 + 1e-12, "peerho: phicm in [0, 2*3.141592654)");
        // pcm: the rho CM momentum must match the generation (pfer = 0: exact two-body)
        const double W = std::sqrt(W2);
        const double Ecm = (W2 + gr.Mh2 - K::Mp2) / 2. / W;
        check(near(tr.ppicm, std::sqrt(Ecm * Ecm - gr.Mh2), 1e-9), "peerho: p_rho in the CM = the generated two-body momentum");
        // theta_cm = angle between p_rho and q in the CM
        {
            const double gb = (vr.nu + K::Mp), bcx = vr.q * vr.uq.x / gb, bcy = vr.q * vr.uq.y / gb,
                         bcz = vr.q * vr.uq.z / gb;
            const double gcm = 1. / std::sqrt(1. - (bcx * bcx + bcy * bcy + bcz * bcz));
            const LorenResult cm = Loren(gcm, bcx, bcy, bcz, vr.p.E, vr.p.P * vr.up.x, vr.p.P * vr.up.y, vr.p.P * vr.up.z);
            const double c = (cm.pxf * vr.uq.x + cm.pyf * vr.uq.y + cm.pzf * vr.uq.z) / cm.pf1;
            check(near(std::cos(tr.thetacm), c, 1e-9), "peerho: theta_cm = angle(p_rho*, q*) of the generated event");
        }
        // the REAL(16) division: sigma = davesig/1.d3 within the double rounding of a binary128 quotient
        check(sig == tr.davesig / 1000. || near(sig, tr.davesig / 1000., 2.3e-16),
              "peerho: davesig/1.d3 (binary128) = davesig/1000 to 1 ulp");
        check(tr.davesig == tr.gtpr * tr.sig, "peerho: davesig = gtpr*sig (no jacobian, l.381-383)");
        // H: R1 has no effect; the Fortran mode writes the on-shell efer (= Mp) back
        MainEvent mc = m0;
        rho::PeerhoConfig pcc = pc;
        pcc.efer_from_common = false;
        const double sigc = rho::peerho(vr, mc, pcc);
        check(sigc == sig && mc.pfermi.efer == K::Mp && mc.davejac == m.davejac, "peerho: H, fortran_compat == corrected (R1)");
        // phi_pq: with pfer = 0 the cross section does not depend on it
        MainEvent mp = m0;
        mp.phi_pq = m0.phi_pq + 1.0;
        const double sigp = rho::peerho(vr, mp, pc);
        check(near(sigp, sig, 1e-12), "peerho: H, independent of phi_pq");
        // a crazy value is clamped: Q2 tiny -> gtpr huge
        Event vz = vr;
        vz.Q2 = 1e-30;
        MainEvent mz = m0;
        std::ostringstream dz;
        rho::PeerhoTrace tz;
        const double sz = rho::peerho(vz, mz, pc, &dz, &tz);
        check(sz == 0. && tz.crazy && dz.str().rfind("crazy sigma_rho", 0) == 0,
              "peerho: sigma > 1e10 -> 0 with the fort.222 'crazy sigma_rho' line");
    }

    // ---- peerho with Fermi motion (D-like): R1 changes efer
    {
        MainEvent m = m0;
        m.pfermi.pfer = 150.;
        m.pfermi.pferx = 0.48;
        m.pfermi.pfery = 0.6;
        m.pfermi.pferz = 0.64;
        const double Md = K::Md;
        m.pfermi.efer = Md - std::sqrt(K::Mn2 + sq(150.));   // event.f:366 (m_spec = Mn)
        MainEvent mc = m;
        rho::PeerhoConfig pd = pc, pdc = pc;
        pdc.efer_from_common = false;
        rho::PeerhoTrace t1, t2;
        const double s1 = rho::peerho(vr, m, pd, nullptr, &t1);
        const double s2 = rho::peerho(vr, mc, pdc, nullptr, &t2);
        check(t1.efer == Md - std::sqrt(K::Mn2 + sq(150.)) && m.pfermi.efer == t1.efer,
              "peerho R1 corrected: the generated off-shell efer is used and kept");
        check(t2.efer == std::sqrt(sq(150.) + K::Mp2) && mc.pfermi.efer == t2.efer && s1 != s2,
              "peerho R1 fortran_compat: on-shell efer, written back to COMMON /pfermi/");
        check(std::isfinite(s1) && std::isfinite(s2) && std::isfinite(m.davejac) && std::isfinite(m.johnjac),
              "peerho: finite with Fermi motion (jacobians included)");
    }

    // ---- rho_decay
    rho::RhoDecayConfig dc;
    // spectrometer along the rho: beam-left arm (SHMS) if the rho goes to +y
    const double th_rho = std::acos(vr.up.z / std::sqrt(sq(vr.up.y) + sq(vr.up.z)));
    dc.hadron_arm = vr.up.y > 0 ? 5 : 1;
    dc.spec_p_theta = th_rho;
    dc.spec_p_phi = vr.up.y > 0 ? K::pi / 2 : 3 * K::pi / 2;
    const double eps = 0.8;
    const double R = rho::rhoR(vr.Q2);
    {
        Event o = vr;
        rho::RhoMassState g = gr;
        NtupVars nt;
        nt.rhomass = rhomass;
        ScriptedRandom rng;
        // phi* = 2 pi 0.25; trial 1: cos = 2*1-1 = 1 (dist = 2 eps R), norm = (1+2 eps R)*0.999 -> reject when
        // 2 eps R < 999; trial 2: cos = 0 (dist = 1), norm = (1+2 eps R)*0.01 -> accept
        rng.load({}, {0.25, 1.0, 0.999, 0.5, 0.01, 0.77});
        std::ostringstream log;
        const auto r = rho::rhoDecay(o, 4000., eps, true, g, nt, dc, rng, &log);
        check(2 * eps * R < 999., "rho_decay: test set-up (2 eps R small)");
        check(r.success && r.n_trials == 2 && r.n_uniform == 5 && rng.nUniform() == 5,
              "rho_decay: success; RNG 1 (phi) + 2 per trial, 2 trials");
        check(nt.rhotheta == std::acos(2. * 0.5 - 1.), "rho_decay: ntup%rhotheta = accepted theta*");
        check(g.Mh == K::Mpi && g.Mh2 == K::Mpi2 && g.Mh2_final == K::Mpi2, "rho_decay: /gnrl/ Mh, Mh2, Mh2_final = pion");
        check(o.p.E == std::sqrt(o.p.P * o.p.P + K::Mpi2), "rho_decay: E = sqrt(P^2 + Mpi2)");
        check(o.p.delta == 100. * (o.p.P / 4000. - 1.), "rho_decay: delta w.r.t. p_spec");
        // boost the lab pion back into the rho rest frame: energy rhomass/2
        const double beta = vr.p.P / std::sqrt(sq(vr.p.P) + gr.Mh2);
        const double gam = 1. / std::sqrt(1. - beta * beta);
        const LorenResult b = Loren(gam, beta * vr.up.x, beta * vr.up.y, beta * vr.up.z, o.p.E, o.p.P * o.up.x,
                                    o.p.P * o.up.y, o.p.P * o.up.z);
        check(near(b.ef1, rhomass / 2., 1e-9), "rho_decay: pion energy in the rho frame = rhomass/2");
        // R9: theta* is measured from the LAB z axis in the rho rest frame (not the helicity axis), phi* from lab x
        check(std::fabs(b.pzf / b.pf1) < 1e-9 && near(std::atan2(b.pyf, b.pxf), 2 * K::pi * 0.25, 1e-9),
              "rho_decay: in the rho frame theta* = pi/2 and phi* = pi/2 w.r.t. the lab axes (R9)");
        // slopes: yptar exact in-plane, xptar = tan(asin(px/p)) (R8: = dx/dz * cos(in-plane offset))
        const double alpha = std::acos(o.up.z / std::sqrt(sq(o.up.y) + sq(o.up.z)));   // in-plane angle from the beam
        const bool sameSide = (dc.hadron_arm == 5) == (o.up.y > 0);
        check(sameSide && near(std::fabs(o.p.yptar), std::fabs(std::tan(alpha - th_rho)), 1e-9),
              "rho_decay: |yptar| = tan of the in-plane angle to the spectrometer axis");
        check(near(o.p.xptar, o.up.x / std::sqrt(sq(o.up.y) + sq(o.up.z)), 1e-12),
              "rho_decay: xptar = tan(asin(px/p)) = px/sqrt(py^2+pz^2) (R8, Fortran)");
        check(std::fabs(std::cos(o.p.theta) - o.up.z) < 2e-3, "rho_decay: physics_angles theta ~ the pion direction");
        check(log.str().empty(), "rho_decay: no message");
        // success_in = false is returned unchanged (R6)
        Event o2 = vr;
        rho::RhoMassState g2 = gr;
        NtupVars n2;
        n2.rhomass = rhomass;
        rng.load({}, {0.25, 1.0, 0.999, 0.5, 0.01});
        const auto r2 = rho::rhoDecay(o2, 4000., eps, false, g2, n2, dc, rng);
        check(!r2.success && o2.p.P == o.p.P && g2.Mh == K::Mpi, "rho_decay: success_in = .false. is kept (R6), event decayed");
    }
    {
        // rhomass < 2 Mpi: fail after the angle loop; rhotheta written, gnrl unchanged
        Event o = vr;
        rho::RhoMassState g = gr;
        NtupVars nt;
        nt.rhomass = 270.;
        ScriptedRandom rng;
        rng.load({}, {0.25, 0.5, 0.0});
        const auto r = rho::rhoDecay(o, 4000., eps, true, g, nt, dc, rng);
        check(!r.success && r.n_uniform == 3 && nt.rhotheta == std::acos(0.) && g.Mh == gr.Mh && o.p.P == vr.p.P,
              "rho_decay: rhomass/2 < Mpi -> fail after the angle draw, event unchanged");
    }
    {
        // spectrometer on the other side at 90 degrees: pzprime < 0
        Event o = vr;
        rho::RhoMassState g = gr;
        NtupVars nt;
        nt.rhomass = rhomass;
        rho::RhoDecayConfig d2 = dc;
        d2.hadron_arm = vr.up.y > 0 ? 1 : 5;
        d2.spec_p_theta = K::pi / 2;
        ScriptedRandom rng;
        rng.load({}, {0.25, 0.5, 0.0});
        const auto r = rho::rhoDecay(o, 4000., eps, true, g, nt, d2, rng);
        check(!r.success && o.p.P == vr.p.P, "rho_decay: pion away from the spectrometer (pzprime < 0) -> fail");
        // arm 7 -> STOP
        d2.hadron_arm = 7;
        rng.load({}, {0.25, 0.5, 0.0});
        bool stopped = false;
        try {
            rho::rhoDecay(o, 4000., eps, true, g, nt, d2, rng);
        } catch (const FortranStop& e) {
            stopped = std::string(e.what()).find("Unknown spectrometer setup dude!") != std::string::npos;
        }
        check(stopped, "rho_decay: hadron_arm 7 -> FortranStop 'Unknown spectrometer setup dude!'");
    }
    check(rho::rhoR(K::Mrho2) == 0.33, "R(Q2 = Mrho2) = 0.33");

    std::printf("%s (%d failure%s)\n", nfail ? "FAIL" : "PASS", nfail, nfail == 1 ? "" : "s");
    return nfail ? 1 : 0;
}

// P13.1 R1 and R2 (FORTRAN-BUG, corrected by default): both modes on a D-like and an H event.
int fortranDiff() {
    std::printf("rho fortran_diff P13.1 R1/R2\n");
    const Event v0 = electronVertex(10600., 6000., 12.0 / K::degrad);
    rho::GenerateRhoConfig gc;
    gc.Mrec_struck = K::Mp;
    gc.spec_p_P = 4000.;
    // R2
    for (const bool fixed : {true, false}) {
        Event v = v0;
        FermiMomentum pf;
        pf.efer = K::Mp;
        rho::RhoMassState g{};
        double rm = 0;
        ScriptedRandom rng;
        rng.load({}, {0.6, 0.3, 0.7});
        gc.mh2_final_fixed = fixed;
        rho::generateRho(v, pf, g, rm, gc, rng);
        check(fixed ? g.Mh2_final == g.Mh * g.Mh : g.Mh2_final == g.Mh,
              std::string("R2 generate_rho.f:80 Mh2_final ") + (fixed ? "= Mh2 (corrected)" : "= Mh (compat)"));
    }
    // R1: H identical, D differs
    Event v = v0;
    FermiMomentum pfH;
    pfH.efer = K::Mp;
    rho::RhoMassState g{};
    double rm = 0;
    ScriptedRandom rng;
    rng.load({}, {0.6, 0.3, 0.7});
    gc.mh2_final_fixed = true;
    rho::generateRho(v, pfH, g, rm, gc, rng);
    rho::PeerhoConfig corr, compat;
    corr.Mtar_struck = compat.Mtar_struck = K::Mp;
    corr.targ_M = compat.targ_M = K::Mp;
    compat.efer_from_common = false;
    MainEvent mH;
    mH.pfermi = pfH;
    mH.phi_pq = phiPq(v);
    MainEvent a = mH, b = mH;
    check(rho::peerho(v, a, corr) == rho::peerho(v, b, compat) && a.davejac == b.davejac && a.pfermi.efer == b.pfermi.efer,
          "R1 H: corrected == compat (bitwise)");
    MainEvent mD = mH;
    mD.pfermi.pfer = 200.;
    mD.pfermi.pferx = 0.;
    mD.pfermi.pfery = 0.6;
    mD.pfermi.pferz = 0.8;
    mD.pfermi.efer = K::Md - std::sqrt(K::Mn2 + 200. * 200.);
    MainEvent c = mD, d = mD;
    const double sc = rho::peerho(v, c, corr), sd = rho::peerho(v, d, compat);
    std::printf("  D-like event: corrected %.10g, compat %.10g (ratio %.6f)\n", sc, sd, sc / sd);
    check(sc != sd && c.pfermi.efer == mD.pfermi.efer && d.pfermi.efer == std::sqrt(200. * 200. + K::Mp2),
          "R1 D: corrected keeps the generated efer, compat writes the on-shell efer back");
    std::printf("%s (%d failure%s)\n", nfail ? "FAIL" : "PASS", nfail, nfail == 1 ? "" : "s");
    return nfail ? 1 : 0;
}

int setup(const std::string& inp) {
    std::printf("rho set-up (P13.1) from %s\n", inp.c_str());
    SimcParameters p;
    InpReader rd(p);
    if (!rd.loadRun(inp)) {
        std::printf("SKIP: cannot read %s\n", inp.c_str());
        return 77;
    }
    p.doing_pion = false;
    p.doing_rho = true;
    p.doing_hplus = true;
    for (const bool compat : {false, true}) {
        SetupOptions so;
        so.fortran_compat = compat;
        RunSetup s;
        bool threw = false;
        try {
            s = BuildRunSetup(p, so);
        } catch (const std::exception& e) {
            threw = true;
            std::printf("  %s\n", e.what());
        }
        const std::string m = compat ? " (compat)" : " (corrected)";
        check(!threw && s.ok(), "BuildRunSetup accepts doing_rho with doing_decay = 1" + m);
        if (threw || !s.ok()) continue;
        check(s.flags.doing_rho && s.flags.doing_hydrho && !s.flags.doing_deutrho && !s.flags.doing_eep,
              "flags doing_rho, doing_hydrho" + m);
        check(s.Mh == K::Mrho && s.Mh2 == K::Mrho * K::Mrho && s.targ.Mtar_struck == K::Mp && s.targ.Mrec_struck == K::Mp &&
                  s.sign_hadron == 1.0,
              "Mh = Mrho, Mtar_struck = Mrec_struck = Mp, sign_hadron (dbase.f:210, 315-321)" + m);
        const driver::PipelineConfig pc = driver::PipelineConfig::Defaults(s);
        check(pc.physics.type == "rho" && pc.reaction.type.empty(),
              "default pipeline: physics \"rho\", reaction not ported yet (P13.2)" + m);
        auto model = Registry<IPhysicsModel>::create("rho");
        model->configure(Json::object(), s.params);
        bool ok = true;
        try {
            model->init(s);
        } catch (const std::exception& e) {
            ok = false;
            std::printf("  %s\n", e.what());
        }
        check(ok, "RhoModel::init accepts the set-up" + m);
        auto* rm = dynamic_cast<rho::RhoModel*>(model.get());
        check(rm && rm->config().efer_from_common == !compat && rho::makeGenerateRhoConfig(s).mh2_final_fixed == !compat,
              "fortran_compat selects R1/R2" + m);
        const auto dcfg = rho::makeRhoDecayConfig(s);
        check(dcfg.hadron_arm == p.hadron_arm && dcfg.spec_p_theta == s.spec.p.theta && dcfg.spec_p_phi == s.spec.p.phi,
              "makeRhoDecayConfig" + m);
        for (const int arm : {2, 3, 4, 7, 8}) {
            RunSetup s2 = s;
            s2.params.hadron_arm = arm;
            auto mdl = Registry<IPhysicsModel>::create("rho");
            mdl->configure(Json::object(), s2.params);
            std::string kind;
            try {
                mdl->init(s2);
                kind = "none";
            } catch (const NotPorted&) {
                kind = "NotPorted";
            } catch (const std::invalid_argument&) {
                kind = "invalid_argument";
            }
            check(kind == (arm <= 4 ? "NotPorted" : "invalid_argument"),
                  "RhoModel::init hadron_arm " + std::to_string(arm) + " -> " + kind + m);
        }
        RunSetup s3 = s;
        s3.flags.doing_rho = false;
        check(!rho::checkSupported(s3).empty(), "checkSupported: doing_rho not set" + m);
        // rho + pion: refused
        SimcParameters q = p;
        q.doing_pion = true;
        bool refused = false;
        try {
            BuildRunSetup(q, so);
        } catch (const std::invalid_argument&) {
            refused = true;
        }
        check(refused, "BuildRunSetup refuses doing_rho together with doing_pion" + m);
    }
    std::printf("%s (%d failure%s)\n", nfail ? "FAIL" : "PASS", nfail, nfail == 1 ? "" : "s");
    return nfail ? 1 : 0;
}
}  // namespace

int main(int argc, char** argv) {
    const std::string what = argc > 1 ? argv[1] : "";
    try {
        if (what == "sanity") return sanity();
        if (what == "fortran_diff") return fortranDiff();
        if (what == "setup" && argc > 2) return setup(argv[2]);
    } catch (const std::exception& e) {
        std::printf("ERROR: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "usage: test_l1_rho sanity | fortran_diff | setup <pion .inp>\n");
    return 2;
}
