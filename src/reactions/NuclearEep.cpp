// src/reactions/NuclearEep.cpp — see include/simc/reactions/NuclearEep.h.
// The "deuterium_eep" and "nucleus_eep" reactions: D(e,e'p) and A(e,e'p) quasi-elastic knock-out, where the missing
// energy and momentum are generated (flat, within the VERTEXedge limits) and weighted by the spectral function
// (ISpectralFunction: the theory momentum distribution or the Benhar table).
// Fortran: simc_gfortran/event.f:generate (D and heavy branches). Operation order follows the Fortran; x**2 is x*x.
#include "simc/reactions/NuclearEep.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "simc/core/Errors.h"
#include "simc/core/Event.h"
#include "simc/core/Registry.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/stages/IRandom.h"

namespace simc {
namespace reactions {

namespace {
// gfortran MAX/MIN of two reals: mvar = a; if (b > mvar || isnan(mvar)) mvar = b  (resp. b < mvar)
inline double fmaxF(double a, double b) { return (b > a || std::isnan(a)) ? b : a; }
inline double fminF(double a, double b) { return (b < a || std::isnan(a)) ? b : a; }
}  // namespace

NuclearEep::NuclearEep(std::string name) : EepReactionBase(std::move(name)) {}

void NuclearEep::configure(const Json& config, const SimcParameters& params) {
    EepReactionBase::configure(config, params);
    params_ = params;
    sfConfig_ = Json::object();
    if (config.is_object() && config.contains("spectral_function")) {
        const Json& c = config.at("spectral_function");
        if (c.is_string())
            sfConfig_["type"] = c.get<std::string>();
        else if (c.is_object())
            sfConfig_ = c;
        else
            throw std::invalid_argument("reaction '" + name_ + "': spectral_function must be a name or an object");
    }
    if (config.is_object() && config.contains("verbose") && !sfConfig_.contains("verbose"))
        sfConfig_["verbose"] = config.at("verbose");
}

void NuclearEep::init(const StageContext& ctx) {
    const RunSetup& s = ctx.setup;
    const ReactionFlags& f = s.flags;
    if (!f.doing_deuterium && !f.doing_heavy)
        throw std::invalid_argument("reaction '" + name_ + "': the input is not D(e,e'p) or A(e,e'p) (doing_deuterium / "
                                    "doing_heavy)");
    EepReactionBase::init(ctx);
    // The spectral function (dbase.f:555-621 file choice in the set-up; RunSetup::sf). Type from the input unless the
    // JSON names one; it must agree with the set-up, which already read the theory file / chose the Benhar file.
    std::string type = s.sf.type;
    if (sfConfig_.contains("type")) {
        const std::string t = sfConfig_.at("type").get<std::string>();
        if ((t == "theory" || t == "benhar") && t != type)
            throw std::invalid_argument("reaction '" + name_ + "': spectral_function type '" + t +
                                        "' differs from the set-up's '" + type +
                                        "' (use_benhar_sf; SimcRun passes the JSON choice to the set-up)");
        type = t;   // a user plug-in (Registry<ISpectralFunction>) gets the set-up's choice through its StageContext
    }
    Json c = sfConfig_;
    c.erase("type");
    if (!c.contains("file")) c["file"] = s.sf.file;
    sf_ = Registry<ISpectralFunction>::create(type);
    sf_->configure(c, params_);
    sf_->init(ctx);
}

// Fortran: event.f:complete_main l.1399-1430 (Benhar: l.1399-1401; theory file: l.1402-1430)
double NuclearEep::spectralWeight(const Event& vertex) {
    if (!sf_) throw std::logic_error("reaction '" + name_ + "': init() was not called");
    return sf_->weight(vertex);
}

// Fortran: event.f:generate (l.214-428), after the beam part (StandardBeam). D(e,e'p) and A(e,e'p) branches.
bool NuclearEep::generate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng) {
    const RunSetup& s = setup();
    const ReactionFlags& f = s.flags;
    const GenLimits& gen = s.gen;

    // l.216-217
    bool success = false;
    main.gen_weight = 1.0;

    // Generate Electron Angles (all cases) (l.276-277)
    vertex.e.yptar = gen.e.yptar.min + rng.uniform() * (gen.e.yptar.max - gen.e.yptar.min);
    vertex.e.xptar = gen.e.xptar.min + rng.uniform() * (gen.e.xptar.max - gen.e.xptar.min);

    // Generate Hadron Angles (all but H(e,e'p)) (l.280-286)
    vertex.p.yptar = gen.p.yptar.min + rng.uniform() * (gen.p.yptar.max - gen.p.yptar.min);
    vertex.p.xptar = gen.p.xptar.min + rng.uniform() * (gen.p.xptar.max - gen.p.xptar.min);

    // Generate Hadron Momentum (A(e,e'p)) (l.289-297)
    double Emin, Emax;
    if (f.doing_heavy) {
        Emin = fmaxF(gen.p.E.min, gen.sumEgen.min - gen.e.E.max);
        Emax = fminF(gen.p.E.max, gen.sumEgen.max - gen.e.E.min);
        if (Emin > Emax) return false;   // goto 100
        main.gen_weight = main.gen_weight * (Emax - Emin) / (gen.p.E.max - gen.p.E.min);
        vertex.p.E = Emin + rng.uniform() * (Emax - Emin);
        vertex.p.P = std::sqrt(vertex.p.E * vertex.p.E - s.Mh2);
        vertex.p.delta = 100. * (vertex.p.P - s.spec.p.P) / s.spec.p.P;
    }

    // Generate Electron Energy (all but hydrogen elastic) (l.300-315)
    Emin = gen.e.E.min;
    Emax = gen.e.E.max;
    if (f.doing_deuterium) {
        Emin = fmaxF(Emin, gen.sumEgen.min);
        Emax = fminF(Emax, gen.sumEgen.max);
    } else if (f.doing_heavy) {   // A(e,e'p)
        Emin = fmaxF(Emin, gen.sumEgen.min - vertex.p.E);
        Emax = fminF(Emax, gen.sumEgen.max - vertex.p.E);
    }
    if (Emin > Emax) return false;   // goto 100
    main.gen_weight = main.gen_weight * (Emax - Emin) / (gen.e.E.max - gen.e.E.min);
    vertex.e.E = Emin + rng.uniform() * (Emax - Emin);
    vertex.e.P = vertex.e.E;
    vertex.e.delta = 100. * (vertex.e.P - s.spec.e.P) / s.spec.e.P;

    // physics angles (l.322-325)
    PhysicsAngles(s.spec.e.theta, s.spec.e.phi, vertex.e.xptar, vertex.e.yptar, vertex.e.theta, vertex.e.phi, log_);
    PhysicsAngles(s.spec.p.theta, s.spec.p.phi, vertex.p.xptar, vertex.p.yptar, vertex.p.theta, vertex.p.phi, log_);

    // l.329-335 (/pfermi/, Em; the Fermi-motion generation l.336-372 is for A(e,e'pi/K) only: none for (e,e'p))
    fermi_.pfer = 0.0;
    fermi_.pferx = 0.0;
    fermi_.pfery = 0.0;
    fermi_.pferz = 0.0;
    vertex.Em = 0.0;
    fermi_.efer = s.targ.Mtar_struck;

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
        success = true;
        if (f.doing_heavy)
            success = (vertex.Em >= s.VERTEXedge.Em.min && vertex.Em <= s.VERTEXedge.Em.max &&
                       vertex.Pm >= s.VERTEXedge.Pm.min && vertex.Pm <= s.VERTEXedge.Pm.max);
        if (success) orig = vertex;
    }
    return success;
}

namespace {
struct DeuteriumEep : NuclearEep {
    DeuteriumEep() : NuclearEep("deuterium_eep") {}
};
}  // namespace

SIMC_REGISTER(IReaction, "nucleus_eep", NuclearEep)
SIMC_REGISTER(IReaction, "deuterium_eep", DeuteriumEep)

}  // namespace reactions
}  // namespace simc
