// src/reactions/pion/PionProduction.cpp — see include/simc/reactions/PionProduction.h (task P9.2).
// Fortran: simc_gfortran/event.f (generate, complete_ev, complete_main: doing_pion branches), simc.f:calculate_central.
// Operation order follows the Fortran source; x**2 is x*x, x**3 is (x*x)*x (__powidf2), D-024.
#include "simc/reactions/PionProduction.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "simc/core/Errors.h"
#include "simc/core/Event.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/material/FortranFormat.h"
#include "simc/physics/nuclear/SpectralFunctionTables.h"
#include "simc/physics/pion/PizeroDecay.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/stages/IRandom.h"

#ifndef SIMC_PION_DATA_DIR
#define SIMC_PION_DATA_DIR "data/pion"
#endif

namespace simc {
namespace reactions {

namespace {
inline double sq(double x) { return x * x; }
// gfortran MAX/MIN of two reals: mvar = a; if (b > mvar || isnan(mvar)) mvar = b  (resp. b < mvar)
inline double fmaxF(double a, double b) { return (b > a || std::isnan(a)) ? b : a; }
inline double fminF(double a, double b) { return (b < a || std::isnan(a)) ? b : a; }
std::string joinPath(const std::string& d, const std::string& f) {
    if (d.empty() || f.empty() || f[0] == '/') return f;
    return d.back() == '/' ? d + f : d + "/" + f;
}
}  // namespace

// ---------------------------------------------------------------------------------------------------------- Delta
// Fortran: event.f:complete_ev l.478-496 (first call): read(55,'(f8.3,f12.5)') w, prob for i = 1..1000.
bool DeltaMassTable::load(const std::string& path, std::vector<std::string>* msgs, std::string* err) {
    mass.clear();
    probsum.clear();
    std::ifstream f(path);
    if (!f) {
        if (err) *err = "cannot open " + path;
        return false;
    }
    auto field = [](const std::string& line, std::size_t pos, std::size_t len) {
        if (pos >= line.size()) return 0.0;   // PAD mode: blank = 0
        const std::string s = line.substr(pos, len);
        return s.find_first_not_of(' ') == std::string::npos ? 0.0 : std::strtod(s.c_str(), nullptr);
    };
    double probtot = 0.;
    std::string line;
    for (int i = 1; i <= 1000; ++i) {
        if (!std::getline(f, line)) {
            if (err) *err = path + ": fewer than 1000 records";
            return false;
        }
        const double w = field(line, 0, 8);
        const double prob = field(line, 8, 12);
        mass.push_back(w);
        probtot = probtot + prob;
        probsum.push_back(probtot);
    }
    for (int i = 1; i <= 1000; ++i) {
        probsum[static_cast<std::size_t>(i - 1)] = probsum[static_cast<std::size_t>(i - 1)] / probtot;
        if ((i / 50) * 50 == i && msgs) {   // write(6,'(''Delta'',i5,2f8.3)')
            char b[64];
            std::snprintf(b, sizeof b, "Delta%5d%8.3f%8.3f", i, mass[static_cast<std::size_t>(i - 1)],
                          probsum[static_cast<std::size_t>(i - 1)]);
            msgs->push_back(b);
        }
    }
    return true;
}

// Fortran: event.f:621-631
double DeltaMassTable::sample(double rn) const {
    double diffmin = 10000.;
    double m = 0.;
    for (std::size_t i = 0; i < probsum.size(); ++i) {
        if (std::fabs(rn - probsum[i]) < diffmin) {
            diffmin = std::fabs(rn - probsum[i]);
            m = mass[i] * 1000.;   // in MeV
        }
    }
    return m;
}

// ---------------------------------------------------------------------------------------------------------- plug-in
PionProduction::PionProduction(std::string name, bool kaon) : EepReactionBase(std::move(name)), kaon_(kaon) {}
PionProduction::~PionProduction() = default;

void PionProduction::configure(const Json& config, const SimcParameters& params) {
    EepReactionBase::configure(config, params);
    dataDir_ = SIMC_PION_DATA_DIR;
    if (config.is_object()) {
        if (config.contains("data_dir")) dataDir_ = config.at("data_dir").get<std::string>();
        if (config.contains("verbose")) verbose_ = config.at("verbose").get<bool>();
    }
}

void PionProduction::init(const StageContext& ctx) {
    const RunSetup& s = ctx.setup;
    const ReactionFlags& f = s.flags;
    if (doingDelta_ ? !f.doing_delta : (kaon_ ? !f.doing_kaon : !f.doing_pion))
        throw std::invalid_argument("reaction '" + name_ + "': the input is not " +
                                    (doingDelta_ ? std::string("(e,e'p)pi in the Delta region (doing_delta)")
                                            : (kaon_ ? std::string("kaon production (doing_kaon)")
                                                     : std::string("exclusive pion production (doing_pion)"))));
    if (f.using_tgt_field || s.params.using_tgt_field)
        throw NotPorted("reaction '" + name_ + "': using_tgt_field (polarized-target angles) not ported");
    EepReactionBase::init(ctx);
    if (!verbose_) log_ = nullptr;
    mrecStruck_ = s.targ.Mrec_struck;
    const int w = s.params.which_pion;
    if (w == 2 || w == 3) {
        std::vector<std::string> msgs;
        std::string err;
        if (!delta_.load(joinPath(dataDir_, "delta_relativistic_bw.inp"), &msgs, &err))
            throw std::invalid_argument("reaction '" + name_ + "': " + err + " (JSON key \"data_dir\" of the reaction)");
        if (log_)
            for (const auto& m : msgs) *log_ << m << "\n";
    }
    // dbase.f:594-621 sf_lookup_init for generate_em. doing_hedelta: see SetupBuilder (D-048 F1, the Fortran never loads
    // the table for the Delta; the A(e,e'p)pi runs only exist in corrected mode).
    if (f.doing_hepi || f.doing_hekaon || f.doing_hedelta) {
        sf_ = std::make_unique<nuclear::SfTable>();
        sf_->setFortranCompat(s.fortran_compat);
        std::vector<std::string> msgs;
        std::string err;
        if (!sf_->load(joinPath(s.data_dir, s.sf.file), s.sf.proton, &msgs, &err))
            throw std::invalid_argument("reaction '" + name_ + "': spectral function " + s.sf.file + ": " + err);
        if (log_)
            for (const auto& m : msgs) *log_ << m << "\n";
    }
}

// Fortran: event.f:generate (l.214-428), after the beam part (StandardBeam). doing_pion branches.
bool PionProduction::generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) {
    using constants::Mn;
    using constants::Mp;
    using constants::pi;
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;
    const GenLimits& gen = s.gen;

    // l.216-217
    bool success = false;
    main.gen_weight = 1.0;

    // Generate Electron Angles (all cases) (l.276-277)
    vertex.e.yptar = gen.e.yptar.min + rng.uniform() * (gen.e.yptar.max - gen.e.yptar.min);
    vertex.e.xptar = gen.e.xptar.min + rng.uniform() * (gen.e.xptar.max - gen.e.xptar.min);
    // Generate Hadron Angles (l.280-286)
    vertex.p.yptar = gen.p.yptar.min + rng.uniform() * (gen.p.yptar.max - gen.p.yptar.min);
    vertex.p.xptar = gen.p.xptar.min + rng.uniform() * (gen.p.xptar.max - gen.p.xptar.min);

    // Generate Electron Energy (l.300-315): sumEgen = generated electron energy for pion production
    double Emin = gen.e.E.min;
    double Emax = gen.e.E.max;
    Emin = fmaxF(Emin, gen.sumEgen.min);
    Emax = fminF(Emax, gen.sumEgen.max);
    if (Emin > Emax) return false;   // goto 100
    main.gen_weight = main.gen_weight * (Emax - Emin) / (gen.e.E.max - gen.e.E.min);
    vertex.e.E = Emin + rng.uniform() * (Emax - Emin);
    vertex.e.P = vertex.e.E;
    vertex.e.delta = 100. * (vertex.e.P - s.spec.e.P) / s.spec.e.P;

    // physics angles (l.322-325)
    PhysicsAngles(s.spec.e.theta, s.spec.e.phi, vertex.e.xptar, vertex.e.yptar, vertex.e.theta, vertex.e.phi, log_);
    PhysicsAngles(s.spec.p.theta, s.spec.p.phi, vertex.p.xptar, vertex.p.yptar, vertex.p.theta, vertex.p.phi, log_);

    // Generate Fermi Momentum and Em for A(e,e'pi) (l.329-372). COMMON /pfermi/ = main.pfermi.
    simc::FermiMomentum& F = main.pfermi;
    F.pfer = 0.0;
    F.pferx = 0.0;
    F.pfery = 0.0;
    F.pferz = 0.0;
    vertex.Em = 0.0;
    F.efer = s.targ.Mtar_struck;   // used for pion/kaon xsec calcs.
    // event.f:336-338 doing_deutpi.or.doing_hepi.or.doing_deutkaon.or.doing_hekaon.or.doing_deutdelta.or.doing_hedelta
    // (kaon: P10.2; delta: P11.1); l.362 and l.368 select the D / A branch with the same flags.
    const bool deut = f.doing_deutpi || f.doing_deutkaon || f.doing_deutdelta,
               he = f.doing_hepi || f.doing_hekaon || f.doing_hedelta;
    if (deut || he) {
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
        if (deut) {   // Em = binding energy
            vertex.Em = Mp + Mn - s.targ.M;
            const double m_spec = s.targ.M - s.targ.Mtar_struck + vertex.Em;   // = Mn(Mp) for pi+(-)
            F.efer = s.targ.M - std::sqrt(sq(m_spec) + sq(F.pfer));
        }
        if (he) {
            vertex.Em = sf_->generateEm(F.pfer, rng);   // call generate_em(pfer, vertex%Em)
            const double m_spec = s.targ.M - s.targ.Mtar_struck + vertex.Em;   // = M^*_{A-1}
            F.efer = s.targ.M - std::sqrt(sq(m_spec) + sq(F.pfer));
        }
    }

    // l.380-393
    success = completeEv(main, vertex, rng);
    main.sigcc = 1.0;
    if (!success) return false;
    main.Trec = vertex.Trec;

    // l.400-415
    if (s.sw.using_rad) {
        // FORTRAN-DIFF-CANDIDATE(D-024): orig keeps the previous try's values when generate_rad fails.
        success = radiator_->radiate(main, vertex, orig, rng, completeEvFn());
    } else {
        success = true;   // (the doing_heavy VERTEXedge test does not apply)
        orig = vertex;
    }
    return success;
}

// Fortran: event.f:720-786 (pion/kaon block of complete_ev, without the polarized-target angles)
void PionProduction::pionMainKinematics(MainEvent& main, const Event& vertex) const {
    using constants::pi;
    const RunSetup& s = setup();
    const double Mts = s.targ.Mtar_struck;
    const double W2 = sq(Mts) + 2. * Mts * vertex.nu - vertex.Q2;
    main.w = std::sqrt(std::fabs(W2)) * W2 / std::fabs(W2);
    main.epsilon = 1. / (1. + 2. * (1 + sq(vertex.nu) / vertex.Q2) * sq(std::tan(vertex.e.theta / 2.)));
    main.theta_pq = std::acos(vertex.up.x * vertex.uq.x + vertex.up.y * vertex.uq.y + vertex.up.z * vertex.uq.z);
    main.t = vertex.Q2 - s.Mh2 + 2 * vertex.nu * vertex.p.E - 2 * vertex.p.P * vertex.q * std::cos(main.theta_pq);
    main.tmin = vertex.Q2 - s.Mh2 + 2 * vertex.p.E * vertex.nu - 2 * vertex.p.P * vertex.q;
    main.q2 = vertex.Q2;

    // phi between the scattering and the reaction plane, in the 'replay' coordinate system (l.761-786)
    const double qx = -vertex.uq.y, qy = vertex.uq.x, qz = vertex.uq.z;
    const double px = -vertex.up.y, py = vertex.up.x, pz = vertex.up.z;
    double dummy = std::sqrt((sq(qx) + sq(qy)) * (sq(qx) + sq(qy) + sq(qz)));
    const double new_x_x = -qx * qz / dummy;
    const double new_x_y = -qy * qz / dummy;
    const double new_x_z = (sq(qx) + sq(qy)) / dummy;
    dummy = std::sqrt(sq(qx) + sq(qy));
    const double new_y_x = qy / dummy;
    const double new_y_y = -qx / dummy;
    const double new_y_z = 0.0;
    const double p_new_x = px * new_x_x + py * new_x_y + pz * new_x_z;
    const double p_new_y = px * new_y_x + py * new_y_y + pz * new_y_z;
    main.phi_pq = std::atan2(p_new_y, p_new_x);
    if (main.phi_pq < 0.e0) main.phi_pq = main.phi_pq + 2. * pi;
}

// Fortran: event.f:complete_ev (l.432-1054), doing_pion branches.
bool PionProduction::completeEv(MainEvent& main, Event& vertex, IRandom& rng) {
    using constants::Me;
    using constants::Mp;
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;
    const double Mh2 = s.Mh2;
    const simc::FermiMomentum& F = main.pfermi;

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

    // pion branch (l.611-698)
    // Delta final states: generate the Delta mass (relativistic Breit-Wigner, l.620-631).
    // FORTRAN-DIFF-CANDIDATE(D-038 F9): the COMMON targ%Mrec_struck is overwritten on every call, also by the second
    // complete_ev call of radiation tail 1 (radc.f:324). Question: is the re-generation intended? Kept in both modes.
    const int which_pion = s.params.which_pion;
    if (which_pion == 2 || which_pion == 3) mrecStruck_ = delta_.sample(rng.uniform());

    vertex.Pm = F.pfer;   // vertex%Em generated at beginning.
    vertex.Mrec = s.targ.M - s.targ.Mtar_struck + vertex.Em;
    double a = -1. * vertex.q * (vertex.uq.x * vertex.up.x + vertex.uq.y * vertex.up.y + vertex.uq.z * vertex.up.z);
    double b = sq(vertex.q);
    double c = vertex.nu + s.targ.M;
    // For nuclei, correct for fermi motion and missing energy (l.648-654)
    // l.646-647, l.949-950, l.964-968 (the same flags incl. doing_deutdelta / doing_hedelta)
    const bool deut = f.doing_deutpi || f.doing_deutkaon || f.doing_deutdelta,
               he = f.doing_hepi || f.doing_hekaon || f.doing_hedelta;
    if (deut || he) {
        a = a - std::fabs(F.pfer) * (F.pferx * vertex.up.x + F.pfery * vertex.up.y + F.pferz * vertex.up.z);
        b = b + sq(F.pfer) +
            2 * vertex.q * std::fabs(F.pfer) * (F.pferx * vertex.uq.x + F.pfery * vertex.uq.y + F.pferz * vertex.uq.z);
        c = vertex.nu + F.efer;
    }
    const double t = sq(c) - b + Mh2 - sq(mrecStruck_);
    const double QA = 4. * (sq(a) - sq(c));
    const double QB = 4. * c * t;
    const double QC = -4. * sq(a) * Mh2 - sq(t);
    const double radical = sq(QB) - 4. * QA * QC;
    if (radical < 0) return false;
    vertex.p.E = (-QB - std::sqrt(radical)) / 2. / QA;
    if (vertex.p.E < 0.0) return false;
    const double Ehad2 = (-QB + std::sqrt(radical)) / 2. / QA;
    if (f.doing_delta) {   // choose one of the two solutions (l.682-684)
        // FORTRAN-DIFF-CANDIDATE(D-048 F5): the second root is taken with probability 0.5 and no factor 1/0.5 enters
        // gen_weight or the jacobian, so when both roots are physical each is generated at half rate with full weight.
        // Question: is a factor 2 missing? Kept in both modes (one grnd, the Fortran random-number order).
        if (rng.uniform() > 0.5) vertex.p.E = Ehad2;
    } else if (Ehad2 > s.edge.p.E.min && Ehad2 < s.edge.p.E.max && ntried_ <= 5000 && log_) {
        // verify that 'backwards' soln. is no good (l.685-689)
        *log_ << " The low-momentum solution to E_hadron is within the spectrometer generation\n"
              << " limits.  If it is in the acceptance, Its the end of the world as we know it!!!\n"
              << " E_hadron solns=" << fortran::listReal(vertex.p.E) << fortran::listReal(Ehad2) << "\n";
    }
    const double E_rec = c - vertex.p.E;      // energy of recoil system
    if (E_rec <= mrecStruck_) return false;   // non-physical solution
    if (vertex.p.E <= s.Mh) return false;
    vertex.p.P = std::sqrt(sq(vertex.p.E) - Mh2);
    vertex.p.delta = (vertex.p.P - s.spec.p.P) * 100. / s.spec.p.P;

    // Compute some pion and kaon stuff (l.720-786)
    pionMainKinematics(main, vertex);
    // (using_tgt_field angles l.789-896: refused in init)
    if (f.doing_pizero) {   // l.899-901: decay the pi0 into 2 photons (2 grnd)
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

    // l.937-955
    if (f.doing_hydpi || f.doing_hydkaon || f.doing_hyddelta) {
        vertex.Trec = 0.0;
    } else if (deut || he) {
        vertex.Trec = std::sqrt(sq(vertex.Mrec) + sq(vertex.Pm)) - vertex.Mrec;
    }

    // krel (l.961-977). FORTRAN-QUIRK(D-038 F11): SAVEd local, never set by the H/coherent branch.
    if (deut) {
        const double mm2 = sq(vertex.Emiss) - sq(vertex.Pmiss);
        if (mm2 < 0 && log_)
            *log_ << " BAD MM!!!!! Emiss,Pmiss=" << fortran::listReal(vertex.Emiss) << fortran::listReal(vertex.Pmiss)
                  << "\n";
        const double MM = std::sqrt(fmaxF(0.e0, mm2));
        krel_ = std::sqrt(fmaxF(0.e0, sq(MM) - 4. * sq(mrecStruck_)));
    } else if (he) {
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

// Fortran: event.f:complete_main (l.1363-1569), doing_pion branch.
bool PionProduction::completeMain(bool force_sigcc, MainEvent& main, const Event& vertex, const Event& /*vertex0*/,
                                  const Event& recon) {
    if (!physics_) throw std::logic_error("reaction '" + name_ + "': no physics model (StageContext::physics)");
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;
    // ntup%mm, mmA, t of complete_recon_ev (event.f:1323-1335, COMMON /ntuple_variables/)
    main.ntup.mm = reconNtup_.mm;
    main.ntup.mmA = reconNtup_.mmA;
    main.ntup.t = reconNtup_.t;

    main.SF_weight = 1.0;   // l.1397-1398
    if (main.SF_weight <= 0 && !force_sigcc) return false;   // l.1436

    double tgtweight = 1.0;
    main.sigcc = physics_->sigccMain(vertex, main);   // l.1459 peepi(vertex, main)
    // Clebsch-Gordan factors for the Delta final states (l.1464-1494)
    const int w = s.params.which_pion;
    if (w == 2) {
        if (f.doing_hydpi) {
            if (f.doing_pizero)
                main.sigcc = 0.55 * main.sigcc;   // g* p -> pi0 Delta+
            else
                main.sigcc = 0.4 * main.sigcc;    // (pi+ Delta0)/(pi+ n)
        } else if (f.doing_deutpi) {
            if (f.doing_pizero)
                main.sigcc = 0.55 * main.sigcc;
            else
                main.sigcc = 0.4 * main.sigcc + 0.8 * main.sigcc;   // + (pi+ Delta-)/(pi+ n)
        }
    } else if (w == 3) {
        if (f.doing_hydpi) {
            if (f.doing_pizero)
                main.sigcc = 0;                   // can't do gamma* n -> pi0 Delta0 for hydpi
            else
                main.sigcc = 0.55 * main.sigcc;   // (pi- Delta++)/(pi- p)
        } else if (f.doing_deutpi) {
            if (f.doing_pizero)
                main.sigcc = 0.99 * main.sigcc;   // g* n -> pi0 Delta0
            else
                main.sigcc = 0.55 * main.sigcc + 0.99 * main.sigcc;   // + (pi- Delta+)/(pi- p)
        }
    }
    main.sigcc_recon = physics_->sigccRecon(recon, main);   // l.1505: 1.0
    // FORTRAN-DIFF-CANDIDATE(D-038 F4): "OK for coherent???" (event.f:1496): Z or N also for coherent production.
    if (w == 1 || w == 11)
        tgtweight = s.targ.N;
    else
        tgtweight = s.targ.Z;

    // Coulomb focusing factor (l.1533-1535)
    if (mainCtx_.using_Coulomb) main.sigcc = main.sigcc * sq(1.0 + mainCtx_.Coulomb_ave / mainCtx_.Ebeam);

    // l.1557-1560
    main.weight = main.SF_weight * main.jacobian * main.gen_weight * main.sigcc;
    main.weight = main.weight * tgtweight;   // correct for #/nucleons involved
    return true;
}

// Fortran: simc.f:calculate_central l.1236-1290
void PionProduction::prepareCentral(MainEvent& main0, const Event& vertex0) {
    const RunSetup& s = setup();
    if (s.fortran_compat) return;
    // FORTRAN-BUG(D-038 F5): original Fortran simc.f:1236-1290 calls only complete_recon_ev(vertex0) before
    // complete_main(.true., main0, ...), so the central peepi call sees main0%epsilon = main0%t = main0%phi_pq = 0
    // (main0 is a SAVEd local): sigma_L ~ |t| = 0 and epsilon = 0, and the .hist CENTRAL.sigcc is not the cross section
    // at the central kinematics; corrected because CENTRAL.sigcc is reported as the central cross section: main0 gets
    // the central W, epsilon, theta_pq, t, tmin, q2, phi_pq of vertex0 (the complete_ev formulas, event.f:720-786).
    pionMainKinematics(main0, vertex0);
}

// simc.f:844 (report): the COMMON targ%Mrec_struck as complete_ev left it (FORTRAN-DIFF-CANDIDATE D-038 F9).
void PionProduction::finalTargetState(TargetInfo& targ) const {
    const int w = setup().params.which_pion;
    if (w == 2 || w == 3) targ.Mrec_struck = mrecStruck_;
}

SIMC_REGISTER(IReaction, "pion_production", PionProduction)

}  // namespace reactions
}  // namespace simc
