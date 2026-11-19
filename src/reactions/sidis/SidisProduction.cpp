// src/reactions/sidis/SidisProduction.cpp — see include/simc/reactions/SidisProduction.h.
// The "sidis_production" reaction, semi-inclusive (e,e'pi)X / (e,e'K)X: the hadron momentum is generated (no two-body
// solution) and the weight is the peepiX cross section.
// Fortran: simc_gfortran/event.f (generate, complete_ev, complete_main: doing_semi branches).
// Operation order follows the Fortran source; x**2 is x*x, r**3 is (r*r)*r (D-024).
#include "simc/reactions/SidisProduction.h"

#include <cmath>
#include <stdexcept>

#include "simc/core/Errors.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/physics/pion/PizeroDecay.h"
#include "simc/physics/sidis/SidisModel.h"
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

SidisProduction::SidisProduction(std::string name) : PionProduction(std::move(name), false) { doingSemi_ = true; }

// Fortran: event.f:generate (l.214-428) after the beam part (StandardBeam), doing_semi branches.
bool SidisProduction::generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) {
    using constants::Mn;
    using constants::Mp;
    using constants::pi;
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;
    const GenLimits& gen = s.gen;

    bool success = false;
    main.gen_weight = 1.0;

    // Generate Electron Angles (all cases) (l.276-277)
    vertex.e.yptar = gen.e.yptar.min + rng.uniform() * (gen.e.yptar.max - gen.e.yptar.min);
    vertex.e.xptar = gen.e.xptar.min + rng.uniform() * (gen.e.xptar.max - gen.e.xptar.min);
    // Generate Hadron Angles (l.280-286)
    vertex.p.yptar = gen.p.yptar.min + rng.uniform() * (gen.p.yptar.max - gen.p.yptar.min);
    vertex.p.xptar = gen.p.xptar.min + rng.uniform() * (gen.p.xptar.max - gen.p.xptar.min);

    // Generate Hadron Momentum (semi-inclusive production) (l.289-298)
    {
        const double Emin = fmaxF(gen.p.E.min, gen.sumEgen.min - gen.e.E.max);
        const double Emax = fminF(gen.p.E.max, gen.sumEgen.max - gen.e.E.min);
        if (Emin > Emax) return false;   // goto 100
        main.gen_weight = main.gen_weight * (Emax - Emin) / (gen.p.E.max - gen.p.E.min);
        vertex.p.E = Emin + rng.uniform() * (Emax - Emin);
        vertex.p.P = std::sqrt(sq(vertex.p.E) - s.Mh2);
        vertex.p.delta = 100. * (vertex.p.P - s.spec.p.P) / s.spec.p.P;
    }

    // Generate Electron Energy (l.301-317): no sumEgen adjustment on the semi branch
    {
        const double Emin = gen.e.E.min;
        const double Emax = gen.e.E.max;
        if (Emin > Emax) return false;
        main.gen_weight = main.gen_weight * (Emax - Emin) / (gen.e.E.max - gen.e.E.min);
        vertex.e.E = Emin + rng.uniform() * (Emax - Emin);
        vertex.e.P = vertex.e.E;
        vertex.e.delta = 100. * (vertex.e.P - s.spec.e.P) / s.spec.e.P;
    }

    // physics angles (l.322-325)
    PhysicsAngles(s.spec.e.theta, s.spec.e.phi, vertex.e.xptar, vertex.e.yptar, vertex.e.theta, vertex.e.phi, log_);
    PhysicsAngles(s.spec.p.theta, s.spec.p.phi, vertex.p.xptar, vertex.p.yptar, vertex.p.theta, vertex.p.phi, log_);

    // Fermi momentum (l.329-372): for semi only doing_deutsemi (drawn whatever do_fermi; peepiX uses it only with
    // do_fermi). FORTRAN-DIFF-CANDIDATE(D-051 F10): no Fermi momentum for A > 2 (doing_hesemi is not in the list).
    simc::FermiMomentum& F = main.pfermi;
    F.pfer = 0.0;
    F.pferx = 0.0;
    F.pfery = 0.0;
    F.pferz = 0.0;
    vertex.Em = 0.0;
    F.efer = s.targ.Mtar_struck;
    if (f.doing_deutsemi) {
        const std::vector<double>& pval = s.pdist.pval;
        const std::vector<double>& mprob = s.pdist.mprob;
        const int nump = static_cast<int>(pval.size());
        if (nump < 1) throw std::logic_error("reaction '" + name_ + "': no momentum distribution in the set-up");
        auto PV = [&](int i) { return pval[static_cast<std::size_t>(i - 1)]; };
        auto MP = [&](int i) { return mprob[static_cast<std::size_t>(i - 1)]; };
        const double ranprob = rng.uniform();
        int ii = 1;
        while (ranprob > MP(ii) && ii < nump) ii = ii + 1;
        double pferlo, pferhi;
        if (ii == 1)
            pferlo = 0;
        else
            pferlo = (PV(ii - 1) + PV(ii)) / 2;
        if (ii == nump)
            pferhi = PV(nump);
        else
            pferhi = (PV(ii) + PV(ii + 1)) / 2;
        F.pfer = pferlo + (pferhi - pferlo) * rng.uniform();
        const double ranth1 = rng.uniform() * 2. - 1.0;
        const double ranth = std::acos(ranth1);
        const double ranph = rng.uniform() * 2. * pi;
        F.pferx = std::sin(ranth) * std::cos(ranph);
        F.pfery = std::sin(ranth) * std::sin(ranph);
        F.pferz = std::cos(ranth);
        vertex.Em = Mp + Mn - s.targ.M;   // Em = binding energy
        const double m_spec = s.targ.M - s.targ.Mtar_struck + vertex.Em;
        F.efer = s.targ.M - std::sqrt(sq(m_spec) + sq(F.pfer));
    }

    // l.380-393
    success = completeEv(main, vertex, rng);
    main.sigcc = 1.0;
    if (!success) return false;
    main.Trec = vertex.Trec;

    // l.400-415
    if (s.sw.using_rad) {
        success = radiator_->radiate(main, vertex, orig, rng, completeEvFn());
    } else {
        success = true;
        orig = vertex;
    }
    return success;
}

// Fortran: event.f:complete_ev (l.432-1054), doing_semi branch.
bool SidisProduction::completeEv(MainEvent& main, Event& vertex, IRandom& rng) {
    using constants::Me;
    using constants::Mp;
    using constants::Mpi0;
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;

    main.jacobian = 1.0;   // l.498

    // l.504-511
    vertex.ue.x = std::sin(vertex.e.theta) * std::cos(vertex.e.phi);
    vertex.ue.y = std::sin(vertex.e.theta) * std::sin(vertex.e.phi);
    vertex.ue.z = std::cos(vertex.e.theta);
    vertex.up.x = std::sin(vertex.p.theta) * std::cos(vertex.p.phi);
    vertex.up.y = std::sin(vertex.p.theta) * std::sin(vertex.p.phi);
    vertex.up.z = std::cos(vertex.p.theta);

    // The q vector (l.530-539)
    vertex.nu = vertex.Ein - vertex.e.E;
    vertex.Q2 = 2 * vertex.Ein * vertex.e.E * (1. - vertex.ue.z);
    vertex.q = std::sqrt(vertex.Q2 + sq(vertex.nu));
    vertex.xbj = vertex.Q2 / 2. / Mp / vertex.nu;
    vertex.uq.x = -vertex.e.P * vertex.ue.x / vertex.q;
    vertex.uq.y = -vertex.e.P * vertex.ue.y / vertex.q;
    vertex.uq.z = (vertex.Ein - vertex.e.P * vertex.ue.z) / vertex.q;
    if (std::fabs(sq(vertex.uq.x) + sq(vertex.uq.y) + sq(vertex.uq.z) - 1) > 0.01)
        throw FortranStop("Error in q vector normalization (event.f:538-539)");

    // semi: no p-side block (the hadron was generated). Pion/kaon stuff (l.723-786), pizero_decay (l.899-901)
    pionMainKinematics(main, vertex);
    if (f.using_tgt_field) polarizedTargetAngles(main, vertex);   // l.796-896
    if (f.doing_pizero) {
        const pion::PizeroDecayResult g = pion::pizeroDecay(vertex, s.Mh, rng, log_);
        if (g.success) {
            main.ntup.gamma1 = g.gamma1;
            main.ntup.gamma2 = g.gamma2;
        }
    }

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
    vertex.PmPer = std::sqrt(fmaxF(0.e0, sq(vertex.Pm) - sq(vertex.PmPar) - sq(vertex.PmOop)));

    // l.952-954 (Trec and Mrec are not assigned on the semi branch)
    vertex.Pm = vertex.Pmiss;
    vertex.Em = vertex.Emiss;
    if (!s.fortran_compat) {
        // FORTRAN-BUG(D-051 F9, as D-033 F5): original Fortran event.f:922: `vertex%PmPer = sqrt( max(0.e0, vertex%Pm**2 -
        // ...) )` is computed BEFORE vertex%Pm = vertex%Pmiss (l.952), i.e. from the previous try's Pm; corrected because
        // PmPer is the component of this event's missing momentum.
        vertex.PmPer = std::sqrt(fmaxF(0.e0, sq(vertex.Pm) - sq(vertex.PmPar) - sq(vertex.PmOop)));
    }

    // l.961-977: krel is not assigned on the semi branch; FORTRAN-QUIRK: ntup%krel = the SAVEd local krel
    main.ntup.krel = krel_;

    // l.979-988: the missing mass must be above the two-pion threshold.
    // FORTRAN-DIFF-CANDIDATE(D-051 F8): (Mp+Mpi0)^2 for every hadron (also kaons), with Mtar_struck, no Fermi motion.
    if (sq(s.targ.Mtar_struck + vertex.nu - vertex.p.E) - sq(vertex.Pmiss) < sq(Mp + Mpi0)) return false;

    // l.990-996
    vertex.zhad = vertex.p.E / vertex.nu;
    const double c = std::cos(main.theta_pq);
    vertex.pt2 = sq(vertex.p.P) * (1.0 - c * c);
    if (vertex.zhad > 1.0) {
        // upstream event.f:996 (e4dd175) prints this line for every such try
        if (log_) *log_ << " WARNING! z=E/nu>1.0. This should never happen. Setting success to .false.\n";
        return false;
    }

    // l.1013-1024: 1/cos**3(theta-theta0) for both arms
    double r = std::sqrt(1. + sq(vertex.e.yptar) + sq(vertex.e.xptar));
    main.jacobian = main.jacobian / (r * r * r);
    r = std::sqrt(1. + sq(vertex.p.yptar) + sq(vertex.p.xptar));
    main.jacobian = main.jacobian / (r * r * r);

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

// Fortran: event.f:complete_main (l.1363-1569), doing_semi branch.
bool SidisProduction::completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& /*vertex0*/,
                                   const Event& recon) {
    if (!physics_) throw std::logic_error("reaction '" + name_ + "': no physics model (StageContext::physics)");
    const RunSetup& s = setup();
    // ntup%mm, mmA, t of complete_recon_ev (event.f:1323-1335, COMMON /ntuple_variables/)
    main.ntup.mm = reconNtup_.mm;
    main.ntup.mmA = reconNtup_.mmA;
    main.ntup.t = reconNtup_.t;

    main.SF_weight = 1.0;   // l.1397-1398
    if (main.SF_weight <= 0 && !force_sigcc) return false;   // l.1436

    const double tgtweight = 1.0;   // l.1444 (not set on the semi branch: the nucleon sum is inside peepiX)
    main.sigcc = physics_->sigccMain(vertex, main);   // l.1521 peepiX(vertex, vertex0, main, survivalprob, .FALSE.)
    // complete_main's SAVEd survivalprob: peepiX sets it when the decay is not simulated (and not below its threshold)
    if (auto* sm = dynamic_cast<sidis::SidisModel*>(physics_))
        survivalprob_ = sm->survivalProbability();
    else if (!s.params.doing_decay)
        survivalprob_ = main.ntup.decdist;
    main.sigcc_recon = physics_->sigccRecon(recon, main);   // l.1522: 1.0
    main.ntup.dilu = 1.0;                                     // l.1525

    // Coulomb focusing factor (l.1533-1535)
    if (mainCtx_.using_Coulomb) {
        const double c = 1.0 + mainCtx_.Coulomb_ave / mainCtx_.Ebeam;
        main.sigcc = main.sigcc * (c * c);
    }

    // l.1557-1561
    main.weight = main.SF_weight * main.jacobian * main.gen_weight * main.sigcc;
    main.weight = main.weight * tgtweight;
    if (s.flags.doing_semika && !s.params.doing_decay) main.weight = main.weight * survivalprob_;
    return true;
}

SIMC_REGISTER(IReaction, "sidis_production", SidisProduction)

}  // namespace reactions
}  // namespace simc
