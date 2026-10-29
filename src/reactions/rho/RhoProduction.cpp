// src/reactions/rho/RhoProduction.cpp — see include/simc/reactions/RhoProduction.h (task P13.2).
// Fortran: simc_gfortran/event.f (generate, complete_ev, complete_main: doing_rho branches), generate_rho.f, rho_decay.f.
// Operation order follows the Fortran source; x**2 is x*x, r**3 is (r*r)*r (D-024).
#include "simc/reactions/RhoProduction.h"

#include <cmath>
#include <stdexcept>

#include "simc/core/Errors.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/material/FortranFormat.h"
#include "simc/physics/rho/RhoModel.h"
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

RhoProduction::RhoProduction(std::string name) : PionProduction(std::move(name), false) { doingRho_ = true; }

void RhoProduction::init(const StageContext& ctx) {
    PionProduction::init(ctx);
    const RunSetup& s = ctx.setup;
    const auto problems = rho::checkSupported(s);
    if (!problems.empty()) {
        for (const auto& p : problems)
            if (rho::isNotPortedProblem(p)) throw NotPorted("reaction '" + name_ + "': " + p);
        throw std::invalid_argument("reaction '" + name_ + "': " + problems.front());
    }
    gcfg_ = rho::makeGenerateRhoConfig(s);
    dcfg_ = rho::makeRhoDecayConfig(s);
    reconMh2_ = s.Mh2;   // calculate_central's complete_recon_ev: Mh2 = Mrho2 (dbase.f:221)
    fixStaleAngle_ = !s.fortran_compat;
    decayOnlyOnSuccess_ = !s.fortran_compat;
}

void RhoProduction::initGnrl(MainEvent& main) const {
    if (main.gnrl_per_event) return;
    const RunSetup& s = setup();
    main.gnrl = HadronMass{s.Mh, s.Mh2, s.Mh2_final};
    main.gnrl_per_event = true;
}

// Fortran: event.f:generate (l.214-428) after the beam part (StandardBeam), doing_rho branches.
bool RhoProduction::generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) {
    using constants::Mn;
    using constants::Mp;
    using constants::pi;
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;
    const GenLimits& gen = s.gen;
    initGnrl(main);
    main.gnrl.Mh2_final = main.gnrl.Mh2;   // simc.f:174 (per try, before generate)

    bool success = false;
    main.gen_weight = 1.0;   // l.216-217
    lastCompleteEvOk_ = lastRadSuccess_ = lastDecayCalled_ = lastDecaySuccess_ = false;
    lastDecayNUnif_ = 0;

    // Generate Electron Angles (l.276-277). No hadron angles for doing_rho (l.280-286).
    vertex.e.yptar = gen.e.yptar.min + rng.uniform() * (gen.e.yptar.max - gen.e.yptar.min);
    vertex.e.xptar = gen.e.xptar.min + rng.uniform() * (gen.e.xptar.max - gen.e.xptar.min);

    // Generate Electron Energy (l.300-316): within the sumEgen window as for pion production
    double Emin = gen.e.E.min;
    double Emax = gen.e.E.max;
    Emin = fmaxF(Emin, gen.sumEgen.min);
    Emax = fminF(Emax, gen.sumEgen.max);
    if (Emin > Emax) return false;   // goto 100 (no rho_decay)
    main.gen_weight = main.gen_weight * (Emax - Emin) / (gen.e.E.max - gen.e.E.min);
    vertex.e.E = Emin + rng.uniform() * (Emax - Emin);
    vertex.e.P = vertex.e.E;
    vertex.e.delta = 100. * (vertex.e.P - s.spec.e.P) / s.spec.e.P;

    // physics angles (l.322-325). FORTRAN-QUIRK(D-057 F1): the hadron slopes are the previous rho's lab theta/phi
    // (generate_rho.f:123-124), see completeEv for the correction.
    PhysicsAngles(s.spec.e.theta, s.spec.e.phi, vertex.e.xptar, vertex.e.yptar, vertex.e.theta, vertex.e.phi, log_);
    PhysicsAngles(s.spec.p.theta, s.spec.p.phi, vertex.p.xptar, vertex.p.yptar, vertex.p.theta, vertex.p.phi, log_);

    // Fermi block (l.329-372). COMMON /pfermi/ = main.pfermi.
    simc::FermiMomentum& F = main.pfermi;
    F.pfer = 0.0;
    F.pferx = 0.0;
    F.pfery = 0.0;
    F.pferz = 0.0;
    vertex.Em = 0.0;
    F.efer = s.targ.Mtar_struck;
    if (f.doing_deutrho || f.doing_herho) {
        // FORTRAN-DIFF-CANDIDATE(D-057 F4): dbase.f:563/594 load no momentum distribution (nump = 0) and no 3He spectral
        // function for doing_deutrho/herho: the loop below stops at ii = 1 and pferlo = pferhi = 0 (the unloaded
        // pval), so pfer = 0 for every event although 4 uniforms are drawn (a nucleon at rest; kept in both modes).
        const double ranprob = rng.uniform();
        (void)ranprob;                         // do while (ranprob.gt.mprob(ii) .and. ii.lt.nump): ii = 1
        const double pferlo = 0;               // ii.eq.1
        const double pferhi = (0. + 0.) / 2;   // (pval(1)+pval(2))/2 of the unloaded table (ii.ne.nump = 0)
        F.pfer = pferlo + (pferhi - pferlo) * rng.uniform();
        const double ranth1 = rng.uniform() * 2. - 1.0;
        const double ranth = std::acos(ranth1);
        const double ranph = rng.uniform() * 2. * pi;
        F.pferx = std::sin(ranth) * std::cos(ranph);
        F.pfery = std::sin(ranth) * std::sin(ranph);
        F.pferz = std::cos(ranth);
        if (f.doing_deutrho) {   // Em = binding energy (l.362-366)
            vertex.Em = Mp + Mn - s.targ.M;
            const double m_spec = s.targ.M - s.targ.Mtar_struck + vertex.Em;
            F.efer = s.targ.M - std::sqrt(sq(m_spec) + sq(F.pfer));
        }
        if (f.doing_herho) {   // l.368-371: generate_em(pfer, vertex%Em) on the EMPTY 3He table (sf_lookup.f:195-245)
            // FORTRAN-DIFF-CANDIDATE(D-057 F4): generate_em draws 2 uniforms (ranprob and the Em smearing) and returns
            // Em = Emval(ind) + dEm(ind)*(u - 0.5) = 0 from the unloaded table (the Fortran reads its arrays out of bounds,
            // sf_lookup Pmval(0)); SIMC++ evaluates the same result without reading the table (both modes).
            (void)rng.uniform();
            vertex.Em = 0.0 + 0.0 * (rng.uniform() - 0.5);
            const double m_spec = s.targ.M - s.targ.Mtar_struck + vertex.Em;   // = M^*_{A-1}
            F.efer = s.targ.M - std::sqrt(sq(m_spec) + sq(F.pfer));
        }
    }

    // l.380-393
    success = completeEv(main, vertex, rng);
    main.sigcc = 1.0;
    if (!success) return false;   // goto 100 (no rho_decay)
    lastCompleteEvOk_ = true;
    main.Trec = vertex.Trec;

    // l.400-415
    if (s.sw.using_rad) {
        // FORTRAN-DIFF-CANDIDATE(D-024): orig keeps the previous try's values when generate_rad fails.
        success = radiator_->radiate(main, vertex, orig, rng, completeEvFn());
        lastRadSuccess_ = success;
    } else {
        success = true;
        orig = vertex;
    }

    // l.417-422: decay the rho before the spectrometers
    success = decay(main, orig, success, rng);
    reconMh2_ = main.gnrl.Mh2;
    return success;
}

// Fortran: event.f:420-422
bool RhoProduction::decay(MainEvent& main, Event& orig, bool success, IRandom& rng) {
    if (decayOnlyOnSuccess_ && !success) {
        // FORTRAN-BUG(D-057 F2): original Fortran event.f:420-422 `if(doing_rho) then call rho_decay(orig,spec%p%P,
        // main%epsilon,success)` for every try that reaches it, also after a failed generate_rad: rho_decay then decays
        // the PREVIOUS try's orig (3.4 % of the tries), draws 1+2k uniforms and may set Mh = Mpi; corrected because a
        // failed try has nothing to decay. No weight changes; the random stream differs from the Fortran's after the
        // first such try (those uniforms are no longer drawn), and Mh/Mh2/ntup%rhotheta/orig keep their values.
        return false;
    }
    const rho::RhoDecayResult r =
        rho::rhoDecay(orig, setup().spec.p.P, main.epsilon, success, main.gnrl, main.ntup, dcfg_, rng, log_);
    lastDecayCalled_ = true;
    lastDecaySuccess_ = r.success;
    lastDecayNUnif_ = r.n_uniform;
    return r.success;
}

// Fortran: event.f:complete_ev (l.432-1052), doing_rho branches.
bool RhoProduction::completeEv(MainEvent& main, Event& vertex, IRandom& rng) {
    using constants::Me;
    using constants::Mp;
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;
    const simc::FermiMomentum& F = main.pfermi;
    initGnrl(main);

    main.jacobian = 1.0;   // l.498

    // l.504-511 (vertex%up is not set for doing_rho)
    vertex.ue.x = std::sin(vertex.e.theta) * std::cos(vertex.e.phi);
    vertex.ue.y = std::sin(vertex.e.theta) * std::sin(vertex.e.phi);
    vertex.ue.z = std::cos(vertex.e.theta);

    // The q vector (l.528-539)
    vertex.nu = vertex.Ein - vertex.e.E;
    vertex.Q2 = 2 * vertex.Ein * vertex.e.E * (1. - vertex.ue.z);
    vertex.q = std::sqrt(vertex.Q2 + sq(vertex.nu));
    vertex.xbj = vertex.Q2 / 2. / Mp / vertex.nu;
    vertex.uq.x = -vertex.e.P * vertex.ue.x / vertex.q;
    vertex.uq.y = -vertex.e.P * vertex.ue.y / vertex.q;
    vertex.uq.z = (vertex.Ein - vertex.e.P * vertex.ue.z) / vertex.q;
    if (std::fabs(sq(vertex.uq.x) + sq(vertex.uq.y) + sq(vertex.uq.z) - 1) > 0.01)
        throw FortranStop("Error in q vector normalization (event.f:538-539)");

    // l.701-708: generate the rho in 4pi in the CM
    const rho::GenerateRhoResult g = rho::generateRho(vertex, F, main.gnrl, main.ntup.rhomass, gcfg_, rng);
    if (!g.success) return false;

    // l.723-786 (t and tmin with the COMMON Mh2 = the BW mass of this try)
    pionMainKinematics(main, vertex, main.gnrl.Mh2);
    // (using_tgt_field angles l.789-896: refused in init)

    // l.914-926. FORTRAN-QUIRK: vertex%Pm is never assigned on the rho branch (PmPer uses the stale value, 0 in a run).
    vertex.Pmx = vertex.p.P * vertex.up.x - vertex.q * vertex.uq.x;
    vertex.Pmy = vertex.p.P * vertex.up.y - vertex.q * vertex.uq.y;
    vertex.Pmz = vertex.p.P * vertex.up.z - vertex.q * vertex.uq.z;
    vertex.Pmiss = std::sqrt(sq(vertex.Pmx) + sq(vertex.Pmy) + sq(vertex.Pmz));
    vertex.Emiss = vertex.nu + s.targ.M - vertex.p.E;
    const double oop_x = -vertex.uq.y;
    const double oop_y = vertex.uq.x;
    vertex.PmPar = (vertex.Pmx * vertex.uq.x + vertex.Pmy * vertex.uq.y + vertex.Pmz * vertex.uq.z);
    vertex.PmOop = (vertex.Pmx * oop_x + vertex.Pmy * oop_y) / std::sqrt(sq(oop_x) + sq(oop_y));
    vertex.PmPer = std::sqrt(fmaxF(0.e0, sq(vertex.Pm) - sq(vertex.PmPar) - sq(vertex.PmOop)));

    // l.947-951 (vertex%Mrec, vertex%Pm never assigned on the rho branch)
    if (f.doing_hydrho) {
        vertex.Trec = 0.0;
    } else if (f.doing_deutrho || f.doing_herho) {
        vertex.Trec = std::sqrt(sq(vertex.Mrec) + sq(vertex.Pm)) - vertex.Mrec;
    }

    // krel (l.960-977): SAVEd local, never set for H
    if (f.doing_deutrho) {
        const double mm2 = sq(vertex.Emiss) - sq(vertex.Pmiss);
        if (mm2 < 0 && log_)
            *log_ << " BAD MM!!!!! Emiss,Pmiss=" << fortran::listReal(vertex.Emiss) << fortran::listReal(vertex.Pmiss)
                  << "\n";
        const double MM = std::sqrt(fmaxF(0.e0, mm2));
        krel_ = std::sqrt(fmaxF(0.e0, sq(MM) - 4. * sq(s.targ.Mrec_struck)));
    } else if (f.doing_herho) {
        const double mm2 = sq(vertex.Emiss) - sq(vertex.Pmiss);
        if (mm2 < 0 && log_)
            *log_ << " BAD MM!!!!! Emiss,Pmiss=" << fortran::listReal(vertex.Emiss) << fortran::listReal(vertex.Pmiss)
                  << "\n";
        const double krelx = vertex.Pmx + 1.5 * F.pferx * F.pfer;
        const double krely = vertex.Pmy + 1.5 * F.pfery * F.pfer;
        const double krelz = vertex.Pmz + 1.5 * F.pferz * F.pfer;
        krel_ = std::sqrt(sq(krelx) + sq(krely) + sq(krelz));
        if (vertex.Em < 6.0) krel_ = -krel_;   // bound state test???
    }
    main.ntup.krel = krel_;

    // l.1013-1022: the electron 1/cos**3 only ("we generate rho's in 4pi ... we don't need no stinkin' Jacobian")
    const double r = std::sqrt(1. + sq(vertex.e.yptar) + sq(vertex.e.xptar));
    main.jacobian = main.jacobian / (r * r * r);

    if (fixStaleAngle_) {
        // FORTRAN-BUG(D-057 F1): original Fortran event.f:325 + 1034: `call physics_angles(spec%p%theta,spec%p%phi,
        // vertex%p%xptar,vertex%p%yptar,vertex%p%theta,vertex%p%phi)` runs in generate on the slopes left by the
        // PREVIOUS try's generate_rho (which stores the rho's lab theta and phi there, generate_rho.f:123-124), and
        // trip_thru_target(3) uses that vertex%p%theta; corrected because the energy loss and the thickness (teff(3),
        // used for the decay pion's multiple scattering) belong to this rho: its lab angles are used.
        vertex.p.theta = vertex.p.xptar;   // acos(pz/p)
        vertex.p.phi = vertex.p.yptar;     // atan2(py, px)
    }

    // l.1031-1040: trip_thru_target 2 and 3 (tail 3 with the rho energy and the BW mass)
    {
        const ITarget::Traversal te = target_->tripThruTarget(ITarget::kElectron, main.target.z - s.targ.zoffset,
                                                              vertex.e.E, vertex.e.theta, Me, ITarget::kSample, rng);
        main.target.Eloss[1] = te.eloss;
        main.target.teff[1] = te.radlen;
        const ITarget::Traversal tp = target_->tripThruTarget(ITarget::kHadron, main.target.z - s.targ.zoffset, vertex.p.E,
                                                              vertex.p.theta, main.gnrl.Mh, ITarget::kSample, rng);
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

// Fortran: event.f:complete_main (l.1363-1569), doing_rho branch.
bool RhoProduction::completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& /*vertex0*/,
                                 const Event& recon) {
    if (!physics_) throw std::logic_error("reaction '" + name_ + "': no physics model (StageContext::physics)");
    const RunSetup& s = setup();
    // ntup%mm, mmA, t of complete_recon_ev (COMMON /ntuple_variables/)
    main.ntup.mm = reconNtup_.mm;
    main.ntup.mmA = reconNtup_.mmA;
    main.ntup.t = reconNtup_.t;

    main.SF_weight = 1.0;   // l.1397-1398
    if (main.SF_weight <= 0 && !force_sigcc) return false;

    main.sigcc = physics_->sigccMain(vertex, main);         // l.1516 peerho(vertex, main)
    main.sigcc_recon = physics_->sigccRecon(recon, main);   // l.1517: 1.0
    const double tgtweight = s.targ.Z + s.targ.N;           // l.1518

    // Coulomb focusing factor (l.1533-1535)
    if (mainCtx_.using_Coulomb) main.sigcc = main.sigcc * sq(1.0 + mainCtx_.Coulomb_ave / mainCtx_.Ebeam);

    // l.1557-1559
    main.weight = main.SF_weight * main.jacobian * main.gen_weight * main.sigcc;
    main.weight = main.weight * tgtweight;
    return true;
}

SIMC_REGISTER(IReaction, "rho_production", RhoProduction)

}  // namespace reactions
}  // namespace simc
