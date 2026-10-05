// src/physics/nuclear/SpectralFunctions.cpp — see include/simc/physics/nuclear/SpectralFunctions.h.
#include "simc/physics/nuclear/SpectralFunctions.h"

#include <iostream>
#include <stdexcept>

#include "simc/core/Event.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcParameters.h"

namespace simc {
namespace nuclear {

namespace {
std::string joinPath(const std::string& dir, const std::string& file) {
    if (dir.empty() || (!file.empty() && file[0] == '/')) return file;
    return dir.back() == '/' ? dir + file : dir + "/" + file;
}
std::string configFile(const Json& c) {
    return (c.is_object() && c.contains("file")) ? c.at("file").get<std::string>() : std::string();
}
}  // namespace

// ---------------------------------------------------------------------------------------------------- "theory"
void TheorySpectralFunction::configure(const Json& config, const SimcParameters&) { file_ = configFile(config); }

void TheorySpectralFunction::init(const StageContext& ctx) {
    const RunSetup& s = ctx.setup;
    if (!s.theory.loaded)
        throw std::invalid_argument("spectral function 'theory': the set-up read no theory file (init.f:theory_init is "
                                    "called for D(e,e'p) and for A(e,e'p) without use_benhar_sf)");
    if (!file_.empty() && file_ != s.theory.file)
        throw std::invalid_argument("spectral function 'theory': file '" + file_ + "' differs from the file read by the "
                                    "set-up ('" + s.theory.file + "'); pass it through SetupOptions::sf_file");
    t_ = &s.theory;
    heavy_ = s.flags.doing_heavy;
}

// Fortran: event.f:complete_main l.1402-1430
double TheorySpectralFunction::weight(const Event& vertex) {
    if (!t_) throw std::logic_error("spectral function 'theory': init() was not called");
    return TheorySfWeight(*t_, heavy_, vertex.Em, vertex.Pm);
}

// ---------------------------------------------------------------------------------------------------- "benhar"
BenharSpectralFunction::BenharSpectralFunction() : log_(&std::cout) {}

void BenharSpectralFunction::configure(const Json& config, const SimcParameters& params) {
    file_ = configFile(config);
    transparency_ = params.transparency;
    if (config.is_object() && config.contains("transparency")) {
        transparency_ = config.at("transparency").get<double>();
    }
    if (config.is_object() && config.contains("verbose") && !config.at("verbose").get<bool>()) log_ = nullptr;
}

// Fortran: dbase.f:592-621 -> sf_lookup.f:sf_lookup_init
void BenharSpectralFunction::init(const StageContext& ctx) {
    const RunSetup& s = ctx.setup;
    if (!s.flags.doing_heavy)
        throw std::invalid_argument("spectral function 'benhar': only for A(e,e'p) with A > 2 (event.f:1399; D(e,e'p) "
                                    "always uses h2.theory)");
    if (file_.empty()) file_ = s.sf.file;
    if (file_.empty())
        throw std::invalid_argument("spectral function 'benhar': no file (set use_benhar_sf or name the file)");
    std::vector<std::string> msgs;
    std::string err;
    table_.setFortranCompat(s.fortran_compat);
    if (!table_.load(joinPath(s.data_dir, file_), s.sf.proton, &msgs, &err)) throw std::invalid_argument(err);
    if (log_)
        for (const auto& m : msgs) *log_ << m << "\n";
    Z_ = s.targ.Z;
}

// Fortran: event.f:1399-1401   call sf_lookup_diff(vertex%Em, vertex%Pm, weight); SF_weight = targ%Z*transparency*weight
double BenharSpectralFunction::weight(const Event& vertex) {
    const double w = table_.lookupDiff(vertex.Em, vertex.Pm);
    return Z_ * transparency_ * w;
}

}  // namespace nuclear
}  // namespace simc

SIMC_REGISTER(simc::ISpectralFunction, "theory", simc::nuclear::TheorySpectralFunction)
SIMC_REGISTER(simc::ISpectralFunction, "benhar", simc::nuclear::BenharSpectralFunction)
