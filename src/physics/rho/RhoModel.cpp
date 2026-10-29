// src/physics/rho/RhoModel.cpp — IPhysicsModel "rho" (task P13.1). See include/simc/physics/rho/RhoModel.h.
#include "simc/physics/rho/RhoModel.h"

#include <iostream>
#include <stdexcept>

#include "simc/core/Errors.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"

namespace simc {
namespace rho {

namespace {
const char* const kNotPortedTag = "[not ported]";
}

std::vector<std::string> checkSupported(const RunSetup& s) {
    std::vector<std::string> v;
    const ReactionFlags& f = s.flags;
    if (!f.doing_rho)
        v.push_back("rho: doing_rho is not set (the \"rho\" cross section is rho_physics.f:peerho; with doing_phsp the "
                    "Fortran switches doing_rho off, dbase.f:227-234)");
    const int arm = s.params.hadron_arm;
    if (arm == 2 || arm == 3 || arm == 4)
        v.push_back(std::string(kNotPortedTag) + " rho: hadron_arm = " + std::to_string(arm) +
                    " (SOS/HRS): these spectrometers are not ported yet (task P14); use 1 = HMS or 5 = SHMS");
    else if (arm != 1 && arm != 5)
        v.push_back("rho: hadron_arm = " + std::to_string(arm) +
                    ": rho_decay.f:92-99 handles only the arms 1-5 (the Fortran STOPs with 'Unknown spectrometer "
                    "setup dude!' at the first decay); use 1 = HMS or 5 = SHMS");
    return v;
}

bool isNotPortedProblem(const std::string& problem) { return problem.rfind(kNotPortedTag, 0) == 0; }

PeerhoConfig makePeerhoConfig(const RunSetup& s) {
    PeerhoConfig c;
    c.Mtar_struck = s.targ.Mtar_struck;
    c.targ_M = s.targ.M;
    c.doing_deutpi = s.flags.doing_deutpi;   // rho_physics.f:89 (false in a rho run)
    c.doing_hepi = s.flags.doing_hepi;
    c.efer_from_common = !s.fortran_compat;   // D-057 F6
    c.nan_guard = !s.fortran_compat;          // D-057 F20
    return c;
}

GenerateRhoConfig makeGenerateRhoConfig(const RunSetup& s) {
    GenerateRhoConfig c;
    c.Mrec_struck = s.targ.Mrec_struck;
    c.spec_p_P = s.spec.p.P;
    c.mh2_final_fixed = !s.fortran_compat;   // D-057 F3
    return c;
}

RhoDecayConfig makeRhoDecayConfig(const RunSetup& s) {
    RhoDecayConfig c;
    c.hadron_arm = s.params.hadron_arm;
    c.spec_p_theta = s.spec.p.theta;
    c.spec_p_phi = s.spec.p.phi;
    c.yptar_sign_fixed = !s.fortran_compat;   // D-057 F9
    c.acos_clamp = !s.fortran_compat;         // D-057 F21
    return c;
}

void RhoModel::configure(const Json& config, const SimcParameters&) {
    if (config.is_object() && config.contains("verbose")) verbose_ = config.at("verbose").get<bool>();
}

void RhoModel::init(const StageContext& ctx) {
    const RunSetup& s = ctx.setup;
    const auto problems = checkSupported(s);
    if (!problems.empty()) {
        std::string m;
        bool notPorted = false;
        for (const auto& p : problems) {
            if (isNotPortedProblem(p)) {
                notPorted = true;
                m += (m.empty() ? "" : "; ") + p.substr(std::string(kNotPortedTag).size() + 1);
            } else {
                m += (m.empty() ? "" : "; ") + p;
            }
        }
        if (notPorted) throw NotPorted(m);
        throw std::invalid_argument(m);
    }
    cfg_ = makePeerhoConfig(s);
    if (!(cfg_.Mtar_struck > 0))
        throw std::invalid_argument("rho: targ%Mtar_struck is not set (the \"rho\" cross section needs a doing_rho set-up)");
    diag_.str("");
}

double RhoModel::sigcc(const Event& vertex, const MainEvent& main) {
    MainEvent m = main;
    return peerho(vertex, m, cfg_, nullptr);
}

double RhoModel::sigccMain(const Event& vertex, MainEvent& main) {
    std::ostringstream d;
    const double r = peerho(vertex, main, cfg_, &d);
    if (!d.str().empty()) {
        // keep the first ~100 diagnostic lines (the Fortran writes every one to fort.222)
        if (diag_.tellp() < 8000) diag_ << d.str();
        if (verbose_) std::cout << " rho (fort.222): " << d.str();
    }
    return r;
}

}  // namespace rho
}  // namespace simc

SIMC_REGISTER(simc::IPhysicsModel, "rho", simc::rho::RhoModel)
