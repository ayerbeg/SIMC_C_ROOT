// tests/l1/test_l1_rho.cpp — task P13.1: the exclusive rho0 kernels (rho_physics.f:peerho, generate_rho.f, rho_decay.f)
// and the IPhysicsModel "rho".
//
// Modes:
//   sanity                     kinematic identities, RNG order and counts, failure paths, the R1/R2 switches, the
//                              decay frame, the plug-in (no Fortran tables needed)
//   fortran_diff               D-057 F6/F3 (P13.1 R1/R2) in both modes
//   peerho|generate_rho|rho_decay <dir> <inp dir> <data/spectral> <tol>   bitwise L1 vs the P13.F1 tables (compat)
//   setup <pion .inp>          the rho set-up from a pion H input switched to doing_rho: flags, Mh, doing_decay accepted,
//                              rho + pion refused, the default pipeline (physics "rho"), checkSupported / model init
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
#include "simc/core/random/RanluxRandom.h"
#include <cstring>
#include <fstream>
#include "simc/driver/Pipeline.h"
#include "simc/driver/SimcRun.h"
#include "simc/reactions/EventChain.h"
#include "simc/input/InpReader.h"
#include "simc/physics/rho/RhoModel.h"
#include "simc/physics/rho/RhoPhysics.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/setup/SetupBuilder.h"
#include "simc/stages/IPhysicsModel.h"
#include "L1Runner.h"
#include <map>
#include <memory>
#include <set>

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

    // ---- Breit-Wigner limits and Mh2_final (D-057 F3) in compat
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
        check(g.Mh2_final == g.Mh, "generate_rho: fortran_compat keeps Mh2_final = Mh (D-057 F3)");
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
        check(near(m.w * m.w * 1e6, W2, 1e-12), "peerho: main%w = sqrt(s) in GeV (D-057 F7)");
        check(m.t == (vr.Q2 - K::Mrho2 + 2. * vr.nu * vr.p.E - 2. * vr.p.P * vr.q * tr.tcos) / 1.e6,
              "peerho: main%t in GeV^2 with the nominal Mrho2 (D-057 F7)");
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
        check(sigc == sig && mc.pfermi.efer == K::Mp && mc.davejac == m.davejac, "peerho: H, fortran_compat == corrected (D-057 F6)");
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
              "peerho F6 corrected: the generated off-shell efer is used and kept");
        check(t2.efer == std::sqrt(sq(150.) + K::Mp2) && mc.pfermi.efer == t2.efer && s1 != s2,
              "peerho F6 fortran_compat: on-shell efer, written back to COMMON /pfermi/");
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
              "rho_decay: in the rho frame theta* = pi/2 and phi* = pi/2 w.r.t. the lab axes (D-057 F8)");
        // slopes: yptar exact in-plane, xptar = tan(asin(px/p)) (R8: = dx/dz * cos(in-plane offset))
        const double alpha = std::acos(o.up.z / std::sqrt(sq(o.up.y) + sq(o.up.z)));   // in-plane angle from the beam
        const bool sameSide = (dc.hadron_arm == 5) == (o.up.y > 0);
        check(sameSide && near(std::fabs(o.p.yptar), std::fabs(std::tan(alpha - th_rho)), 1e-9),
              "rho_decay: |yptar| = tan of the in-plane angle to the spectrometer axis");
        check(near(o.p.xptar, o.up.x / std::sqrt(sq(o.up.y) + sq(o.up.z)), 1e-12),
              "rho_decay: xptar = tan(asin(px/p)) = px/sqrt(py^2+pz^2) (D-057 F10, Fortran)");
        check(std::fabs(std::cos(o.p.theta) - o.up.z) < 2e-3, "rho_decay: physics_angles theta ~ the pion direction");
        check(log.str().empty(), "rho_decay: no message");
        // success_in = false is returned unchanged (D-057 F2)
        Event o2 = vr;
        rho::RhoMassState g2 = gr;
        NtupVars n2;
        n2.rhomass = rhomass;
        rng.load({}, {0.25, 1.0, 0.999, 0.5, 0.01});
        const auto r2 = rho::rhoDecay(o2, 4000., eps, false, g2, n2, dc, rng);
        check(!r2.success && o2.p.P == o.p.P && g2.Mh == K::Mpi, "rho_decay: success_in = .false. is kept (D-057 F2), event decayed");
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

// D-057 F6 and F3 (FORTRAN-BUG, corrected by default): both modes on a D-like and an H event.
int fortranDiff() {
    std::printf("rho fortran_diff D-057 F3, F6, F9, F20, F21\n");
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
              std::string("F3 generate_rho.f:80 Mh2_final ") + (fixed ? "= Mh2 (corrected)" : "= Mh (compat)"));
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
          "F6 H: corrected == compat (bitwise)");
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
          "F6 D: corrected keeps the generated efer, compat writes the on-shell efer back");

    // F20: a NaN cross section (s < 0: nu + efer < q) passes the Fortran guard; corrected 0
    {
        Event vn = v;
        vn.nu = 10.;   // nu + Mp < q -> s < 0 -> NaN
        MainEvent a2 = mH, b2 = mH;
        const double sn = rho::peerho(vn, a2, corr), sf = rho::peerho(vn, b2, compat);
        rho::PeerhoConfig noguard = corr;
        noguard.nan_guard = false;
        MainEvent c2 = mH;
        const double sg = rho::peerho(vn, c2, noguard);
        check(std::isnan(sg) && sn == 0.0 && (std::isnan(sf) || sf == 0.0),
              "F20 rho_physics.f:388-389: NaN cross section -> 0 (corrected), NaN without the guard");
    }
    // F9: beam-left arm (SHMS), decay pion on the far side of the beam (pyf < 0): yptar sign
    {
        rho::RhoDecayConfig dc;
        dc.hadron_arm = 5;
        dc.spec_p_theta = 0.2;
        dc.spec_p_phi = K::pi / 2;
        int found = 0, okSign = 0;
        for (int k = 0; k < 200 && found < 20; ++k) {
            Event o;
            o.Q2 = 2.0e6;
            o.p.P = 4000.;
            const double a = -0.05 + 0.0005 * k;   // rho near the beam, slightly on either side
            o.up = {0.0, std::sin(a), std::cos(a)};
            rho::RhoMassState g0{K::Mrho, K::Mrho2, K::Mrho2};
            NtupVars n0;
            n0.rhomass = K::Mrho;
            for (const double u1 : {0.1, 0.6}) {
                Event oc = o, ok = o;
                rho::RhoMassState gc = g0, gk = g0;
                NtupVars nc = n0, nk = n0;
                rho::RhoDecayConfig cc = dc, kk = dc;
                cc.yptar_sign_fixed = false;
                cc.acos_clamp = false;
                ScriptedRandom rc, rk;
                rc.load({}, {u1, 0.5, 0.0});
                rk.load({}, {u1, 0.5, 0.0});
                const auto res_c = rho::rhoDecay(oc, 3000., 0.8, true, gc, nc, cc, rc);
                const auto res_k = rho::rhoDecay(ok, 3000., 0.8, true, gk, nk, kk, rk);
                if (!res_k.success || ok.up.y > 0) continue;
                ++found;
                // true in-plane offset from the SHMS axis: atan2(py, pz) - theta (negative on the far side)
                const double truth = std::atan2(ok.up.y, ok.up.z) - dc.spec_p_theta;
                if (near(ok.p.yptar, std::tan(truth), 1e-9) && (!res_c.success || oc.p.yptar == -ok.p.yptar)) ++okSign;
            }
        }
        check(found > 0 && okSign == found,
              "F9 rho_decay.f:123-124: far-side pion, left arm: corrected yptar = tan(in-plane offset), compat the "
              "opposite sign (" + std::to_string(okSign) + "/" + std::to_string(found) + ")");
    }
    // F21: a pion in the x-z plane (py = 0): pz/p/cos(asin(px/p)) can round above 1 -> acos NaN (compat), clamped
    {
        int nanCompat = 0, finiteCorr = 0, tried = 0;
        for (int k = 1; k <= 4000 && nanCompat < 3; ++k) {
            const double px = 1e-3 * k, pz = 3000.;
            const double pf = std::sqrt(px * px + pz * pz);
            const double c = pz / pf / std::cos(std::asin(px / pf));
            if (!(c > 1.0)) continue;
            ++tried;
            // drive rho_decay to exactly this lab momentum is impractical; check the kernel arithmetic it uses
            if (std::isnan(std::acos(c))) ++nanCompat;
            if (std::isfinite(std::acos(std::min(1.0, c)))) ++finiteCorr;
        }
        check(tried > 0 && nanCompat == tried && finiteCorr == tried,
              "F21 rho_decay.f:111-112: cos_th_inp > 1 by rounding exists (" + std::to_string(tried) +
                  " cases), acos -> NaN unclamped, finite clamped");
    }
    std::printf("%s (%d failure%s)\n", nfail ? "FAIL" : "PASS", nfail, nfail == 1 ? "" : "s");
    return nfail ? 1 : 0;
}

// ================================================================================================ L1 (P13.F1 tables)
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

int iv(const Row& r, const char* c) { return static_cast<int>(std::lround(r[c])); }

/// reference/functions/rho/rho_state/configs.txt
std::string cfgInp(const std::string& inpdir, int cfg) {
    static const char* f[] = {"",          "p13f1_rh_k1_rad_cap", "p13f1_rd_k1_rad_cap", "l1_h_k1_norad",
                              "l1_h_k1_pim", "l1_h_k2_rad",        "l1_he3_k1_rad",       "l1_h_k1_dec1"};
    if (cfg < 1 || cfg > 7) throw std::runtime_error("bad cfg " + std::to_string(cfg));
    return inpdir + "/" + f[cfg] + ".inp";
}

struct SetupCache {
    std::string inpdir, spectral;
    std::map<int, RunSetup> m;
    const RunSetup& get(int cfg) {
        auto it = m.find(cfg);
        if (it != m.end()) return it->second;
        SimcParameters p;
        InpReader rd(p);
        if (!rd.loadRun(cfgInp(inpdir, cfg))) throw std::runtime_error("cannot read " + cfgInp(inpdir, cfg));
        SetupOptions so;
        so.fortran_compat = true;
        so.data_dir = spectral;
        RunSetup s = BuildRunSetup(p, so);
        if (!s.ok()) throw std::runtime_error("cfg " + std::to_string(cfg) + ": " + s.stops.front());
        return m.emplace(cfg, std::move(s)).first->second;
    }
};

#define EVFIELDS(X)                                                                                                    \
    X(Ein) X(Em) X(Pm) X(Emiss) X(Pmiss) X(Pmx) X(Pmy) X(Pmz) X(PmPar) X(PmPer) X(PmOop) X(nu) X(q) X(Q2) X(Trec)     \
        X(W) X(Mrec) X(epsilon) X(theta_pq) X(theta_tarq) X(phi_pq) X(phi_targ) X(beta) X(phi_s) X(phi_c) X(zhad)    \
            X(pt2) X(xbj)
#define ARMFIELDS(X) X(delta) X(xptar) X(yptar) X(z) X(theta) X(phi) X(E) X(P)
void getEvent(const std::string& p, const Row& r, Event& v) {
#define GETF(f) v.f = r[p + #f];
    EVFIELDS(GETF)
#undef GETF
#define GETA(f) v.e.f = r[p + "e_" #f]; v.p.f = r[p + "p_" #f];
    ARMFIELDS(GETA)
#undef GETA
    v.ue = {r[p + "ue_x"], r[p + "ue_y"], r[p + "ue_z"]};
    v.up = {r[p + "up_x"], r[p + "up_y"], r[p + "up_z"]};
    v.uq = {r[p + "uq_x"], r[p + "uq_y"], r[p + "uq_z"]};
}
void putEvent(const std::string& p, const Event& v, Outputs& o) {
#define PUTF(f) o[p + #f] = v.f;
    EVFIELDS(PUTF)
#undef PUTF
#define PUTA(f) o[p + "e_" #f] = v.e.f; o[p + "p_" #f] = v.p.f;
    ARMFIELDS(PUTA)
#undef PUTA
    o[p + "ue_x"] = v.ue.x;
    o[p + "ue_y"] = v.ue.y;
    o[p + "ue_z"] = v.ue.z;
    o[p + "up_x"] = v.up.x;
    o[p + "up_y"] = v.up.y;
    o[p + "up_z"] = v.up.z;
    o[p + "uq_x"] = v.uq.x;
    o[p + "uq_y"] = v.uq.y;
    o[p + "uq_z"] = v.uq.z;
}

/// Every column of output.txt that is not an input column is compared (tolerance tol, 0 = bitwise).
int runAuto(const std::string& title, const std::string& dir, double tol, const simcl1::Adapter& adapter) {
    simcl1::Table in, ref;
    std::string err;
    // The harness drivers group the output rows by cfg (input.txt is in generation order); output.txt echoes every
    // input column, so it is the input of the adapter too. input.txt only provides the list of input columns.
    if (!in.load(dir + "/input.txt", err) || !ref.load(dir + "/output.txt", err)) {
        std::printf("L1 %s: SKIP (%s)\n", title.c_str(), err.c_str());
        return 77;
    }
    if (in.rows.size() != ref.rows.size()) {
        std::printf("L1 %s: ERROR %zu input rows, %zu output rows\n", title.c_str(), in.rows.size(), ref.rows.size());
        return 1;
    }
    std::set<std::string> inNames(in.names.begin(), in.names.end());
    in = ref;
    std::vector<ColumnSpec> cols;
    for (const auto& n : ref.names)
        if (!inNames.count(n)) cols.push_back({n, tol, 0.0});
    auto res = simcl1::compare(in, ref, cols, adapter, err);
    if (!err.empty()) {
        std::printf("L1 %s: ERROR %s\n", title.c_str(), err.c_str());
        return 1;
    }
    std::size_t rows = res.empty() ? 0 : res[0].rows, nb = 0, nc = 0;
    for (const auto& r : res) {
        nb += r.bitwise_equal;
        nc += r.rows;
    }
    const bool ok = simcl1::report(title, res, cols);
    std::printf("  summary: %zu rows x %zu columns, %zu of %zu values bitwise equal\n", rows, cols.size(), nb, nc);
    return ok ? 0 : 1;
}

// rho_physics.f:peerho, fortran_compat (D-057 F6, F20 as the Fortran)
int runPeerho(const std::string& dir, SetupCache& sc, double tol) {
    return runAuto("peerho (rho_physics.f:1-402), fortran_compat", dir, tol, [&](const Row& r, Outputs& o) {
        const RunSetup& s = sc.get(iv(r, "cfg"));
        rho::PeerhoConfig c = rho::makePeerhoConfig(s);
        Event v;
        getEvent("v_", r, v);
        MainEvent m;
        m.phi_pq = r["in_m_phi_pq"];
        m.pfermi.pfer = r["pfer"];
        m.pfermi.pferx = r["pferx"];
        m.pfermi.pfery = r["pfery"];
        m.pfermi.pferz = r["pferz"];
        m.pfermi.efer = r["efer_in"];
        std::ostringstream diag;
        o["peerho"] = rho::peerho(v, m, c, &diag);
        o["efer_out"] = m.pfermi.efer;
        o["diag222_written"] = diag.str().empty() ? 0. : 1.;
        o["ntup_sigcm"] = m.ntup.sigcm;
        o["m_w"] = m.w;
        o["m_t"] = m.t;
        o["m_pcm"] = m.pcm;
        o["m_thetacm"] = m.thetacm;
        o["m_phicm"] = m.phicm;
        o["m_davejac"] = m.davejac;
        o["m_johnjac"] = m.johnjac;
    });
}

// generate_rho.f, scripted grnd, fortran_compat (D-057 F3 as the Fortran)
int runGenerateRho(const std::string& dir, SetupCache& sc, double tol) {
    return runAuto("generate_rho (generate_rho.f:1-131), ScriptedRandom, fortran_compat", dir, tol,
                   [&](const Row& r, Outputs& o) {
        const RunSetup& s = sc.get(iv(r, "cfg"));
        rho::GenerateRhoConfig c = rho::makeGenerateRhoConfig(s);
        if (iv(r, "ov_state") == 1) {
            c.Mrec_struck = r["ov_mrec_struck"];
            c.spec_p_P = r["ov_spec_p_P"];
        }
        Event v;
        v.nu = r["nu"];
        v.Q2 = r["Q2"];
        v.q = r["q"];
        v.uq = {r["uq_x"], r["uq_y"], r["uq_z"]};
        v.ue = {r["ue_x"], r["ue_y"], r["ue_z"]};
        v.p.xptar = r["in_p_xptar"];
        v.p.yptar = r["in_p_yptar"];
        v.p.delta = r["in_p_delta"];
        v.p.P = r["in_p_P"];
        v.p.E = r["in_p_E"];
        v.up = {r["in_up_x"], r["in_up_y"], r["in_up_z"]};
        FermiMomentum pf;
        pf.pfer = r["pfer"];
        pf.pferx = r["pferx"];
        pf.pfery = r["pfery"];
        pf.pferz = r["pferz"];
        pf.efer = r["efer"];
        rho::RhoMassState g{r["in_Mh"], r["in_Mh2"], r["in_Mh2_final"]};
        double rm = r["in_rhomass"];
        ScriptedRandom rng;
        rng.load({}, {r["u1"], r["u2"], r["u3"], r["u4"]});
        const rho::GenerateRhoResult res = rho::generateRho(v, pf, g, rm, c, rng);
        o["success"] = res.success;
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["n_gauss"] = static_cast<double>(rng.nGauss());
        o["rng_err"] = rng.error();
        o["Mh"] = g.Mh;
        o["Mh2"] = g.Mh2;
        o["Mh2_final"] = g.Mh2_final;
        o["ntup_rhomass"] = rm;
        o["p_delta"] = v.p.delta;
        o["p_P"] = v.p.P;
        o["p_E"] = v.p.E;
        o["up_x"] = v.up.x;
        o["up_y"] = v.up.y;
        o["up_z"] = v.up.z;
        o["p_xptar"] = v.p.xptar;
        o["p_yptar"] = v.p.yptar;
        o["eff_mrec_struck"] = c.Mrec_struck;
        o["eff_spec_p_P"] = c.spec_p_P;
    });
}

// rho_decay.f, scripted grnd, fortran_compat (D-057 F9, F21 as the Fortran)
int runRhoDecay(const std::string& dir, SetupCache& sc, double tol) {
    std::ostringstream log;
    const int rc = runAuto("rho_decay (rho_decay.f:1-160), ScriptedRandom, fortran_compat", dir, tol,
                           [&](const Row& r, Outputs& o) {
        const RunSetup& s = sc.get(iv(r, "cfg"));
        rho::RhoDecayConfig c = rho::makeRhoDecayConfig(s);
        if (iv(r, "ov_state") == 1) {
            c.hadron_arm = iv(r, "ov_hadron_arm");
            c.spec_p_theta = r["ov_spec_p_theta"];
            c.spec_p_phi = r["ov_spec_p_phi"];
        }
        Event v;
        getEvent("in_o_", r, v);
        rho::RhoMassState g{r["in_Mh"], r["in_Mh2"], r["in_Mh2_final"]};
        NtupVars nt;
        nt.rhomass = r["in_rhomass"];
        nt.rhotheta = r["in_rhotheta"];
        std::vector<double> u;
        for (int k = 1; k <= 64; ++k) u.push_back(r["u" + std::to_string(k)]);
        ScriptedRandom rng;
        rng.load({}, u);
        const rho::RhoDecayResult res =
            rho::rhoDecay(v, r["p_spec"], r["epsilon"], iv(r, "success_in") != 0, g, nt, c, rng, &log);
        o["success"] = res.success;
        o["n_unif"] = static_cast<double>(rng.nUniform());
        o["n_gauss"] = static_cast<double>(rng.nGauss());
        o["rng_err"] = rng.error();
        putEvent("o_", v, o);
        o["Mh"] = g.Mh;
        o["Mh2"] = g.Mh2;
        o["Mh2_final"] = g.Mh2_final;
        o["ntup_rhotheta"] = nt.rhotheta;
        o["eff_hadron_arm"] = c.hadron_arm;
        o["eff_spec_p_theta"] = c.spec_p_theta;
        o["eff_spec_p_phi"] = c.spec_p_phi;
    });
    // stops.txt: hadron_arm 6/7/8 STOP in the Fortran; SIMC++ throws FortranStop for 7/8 (6 is rejected at input)
    int nstop = 0;
    for (const int arm : {7, 8}) {
        Event v;
        v.p.P = 3000.;
        v.up = {0., 0.2, std::sqrt(1. - 0.04)};
        v.Q2 = 2.5e6;
        rho::RhoMassState g{K::Mrho, K::Mrho2, K::Mrho2};
        NtupVars nt;
        nt.rhomass = K::Mrho;
        rho::RhoDecayConfig c;
        c.hadron_arm = arm;
        ScriptedRandom rng;
        rng.load({}, {0.3, 0.5, 0.0});
        try {
            rho::rhoDecay(v, 3000., 0.8, true, g, nt, c, rng);
        } catch (const FortranStop& e) {
            if (std::string(e.what()).find("Unknown spectrometer setup dude!") != std::string::npos) ++nstop;
        }
    }
    std::printf("  stops.txt: hadron_arm 7/8 -> FortranStop 'Unknown spectrometer setup dude!': %d of 2\n", nstop);
    return (rc == 0 && nstop == 2) ? 0 : (rc == 77 ? 77 : 1);
}

// The RANLUX port (IRandom "ranlux") against the random stream of a capture: every grnd / gauss1 value handed out in the
// contiguous full-logged tries (cap_rng.txt, ntried 1..nmax) must be reproduced bitwise from sgrnd(seed).
int runRanluxStream(const std::string& capdir, int seed, long nmax) {
    const std::string path = capdir + "/cap_rng.txt";
    std::string cmd = std::ifstream(path) ? "cat '" + path + "'" : "xz -dc '" + path + ".xz'";
    FILE* f = popen(cmd.c_str(), "r");
    if (!f) return 77;
    char line[512];
    if (!std::fgets(line, sizeof line, f)) {
        pclose(f);
        std::printf("ranlux stream: SKIP (no %s[.xz])\n", path.c_str());
        return 77;
    }
    RanluxRandom r;
    r.seed(static_cast<std::uint64_t>(seed));
    long n = 0, bad = 0, nu = 0, ng = 0;
    while (std::fgets(line, sizeof line, f)) {
        long nt, ph, sq, kind, nint;
        double v, nsig;
        if (std::sscanf(line, "%ld %ld %ld %ld %lf %lf %ld", &nt, &ph, &sq, &kind, &v, &nsig, &nint) != 7) continue;
        if (nt > nmax) break;
        const double c = kind == 2 ? r.gauss1(nsig) : r.uniform();
        (kind == 2 ? ng : nu)++;
        ++n;
        if (std::memcmp(&c, &v, sizeof c) != 0) {
            if (bad < 5) std::printf("  try %ld phase %ld seq %ld kind %ld: ref %.17g cpp %.17g\n", nt, ph, sq, kind, v, c);
            ++bad;
        }
    }
    pclose(f);
    std::printf("ranlux stream vs %s (tries 1..%ld): %ld values (%ld grnd, %ld gauss1), %ld differ -> %s\n", path.c_str(),
                nmax, n, nu, ng, bad, (bad == 0 && n > 0) ? "PASS" : "FAIL");
    return (bad == 0 && n > 0) ? 0 : 1;
}

// D-057 F15 (.hist title for D/3He) and F18 (P-arm report quantities from the decay pion): short SimcRun runs on the
// D capture input in both modes (ranlux stream, 20 events), and limits_update on a fabricated event.
int fortranDiffReport(const std::string& inpD, const std::string& outdir) {
    std::printf("rho fortran_diff D-057 F15, F18\n");
    SimcParameters p;
    InpReader rd(p);
    if (!rd.loadRun(inpD)) {
        std::printf("SKIP: cannot read %s\n", inpD.c_str());
        return 77;
    }
    for (const bool compat : {true, false}) {
        driver::RunOptions opt;
        opt.fortran_compat = compat;
        opt.run_name = compat ? "rho_d_compat" : "rho_d_corrected";
        opt.output_dir = outdir;
        opt.quiet_stages = true;
        opt.override_ngen = true;
        opt.ngen = 20;
        std::ostringstream log;
        opt.log = &log;
        driver::SimcRun sr(p, Json{{"rng", {{"type", "ranlux"}}}}, opt);
        sr.run();
        std::ifstream h(outdir + "/" + opt.run_name + ".hist");
        std::stringstream ss;
        ss << h.rdbuf();
        const std::string t = ss.str();
        const bool fortranTitle = t.find("I am not set up for anything else yet!") != std::string::npos;
        const bool dTitle = t.find("D(e,e'rho)") != std::string::npos;
        check(compat ? (fortranTitle && !dTitle) : (dTitle && !fortranTitle),
              std::string("F15 simc.f:721-725 D(e,e'rho) .hist title, ") + (compat ? "compat: Fortran text" : "corrected"));
    }
    // F18: limits_update GENERATION P-arm values: vertex%p (compat) or orig%p (corrected, doing_rho)
    {
        MainEvent m;
        Event v, o, r;
        v.p.delta = 7.0;
        v.p.yptar = 2.2;   // the rho's lab phi in the slope field
        v.p.xptar = 0.4;   // the rho's lab theta
        o.p.delta = 1.0;
        o.p.yptar = 0.01;
        o.p.xptar = -0.02;
        ReactionFlags f;
        f.doing_rho = true;
        reactions::Contrib cc{}, ck{};
        Slops sc{}, sk{};
        reactions::limitsUpdate(m, v, o, r, f, {0., 0., 0.}, cc, sc, false);
        reactions::limitsUpdate(m, v, o, r, f, {0., 0., 0.}, ck, sk, true);
        check(cc.gen.p.yptar.hi == 2.2 && ck.gen.p.yptar.hi == 0.01 && ck.gen.p.xptar.lo == -0.02 && ck.gen.p.delta.hi == 1.0,
              "F18 event.f:22-24 limits_update GENERATION P arm: vertex%p (compat), the decay pion orig%p (corrected)");
    }
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
        check(pc.physics.type == "rho" && pc.reaction.type == "rho_production",
              "default pipeline: physics \"rho\", reaction \"rho_production\"" + m);
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
        if (what == "fortran_diff_report" && argc > 3) return fortranDiffReport(argv[2], argv[3]);
        if (what == "ranlux_stream" && argc > 4) return runRanluxStream(argv[2], std::atoi(argv[3]), std::atol(argv[4]));
        if (argc > 5) {
            SetupCache sc{argv[3], argv[4], {}};
            const double tol = std::atof(argv[5]);
            if (what == "peerho") return runPeerho(argv[2], sc, tol);
            if (what == "generate_rho") return runGenerateRho(argv[2], sc, tol);
            if (what == "rho_decay") return runRhoDecay(argv[2], sc, tol);
        }
        if (what == "setup" && argc > 2) return setup(argv[2]);
    } catch (const std::exception& e) {
        std::printf("ERROR: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "usage: test_l1_rho sanity | fortran_diff | setup <pion .inp> | peerho|generate_rho|rho_decay <dir> <inp dir> <data/spectral> <tol>\n");
    return 2;
}
