// src/physics/proton/ProtonModels.cpp
//
// IPhysicsModel plug-ins for (e,e'p) on top of ProtonCrossSections (physics_proton.f).
// Fortran dispatch: event.f:complete_main l.1446-1456
//     doing_hyd_elast            -> main%sigcc = sigep(vertex),    sigcc_recon = sigep(recon)
//     doing_deuterium/heavy      -> main%sigcc = deForest(vertex), sigcc_recon = deForest(recon)
// Choosing the plug-in from the reaction flags is the driver's job; the names are:
//     "sigep", "deforest" (flag = SimcParameters::deForest_flag, or JSON "deForest_flag"),
//     "deforest_cc1" (0), "deforest_cc2" (+1), "deforest_cc1_onshell" (-1).
// JSON config of the deforest models (all optional):
//     { "deForest_flag": int (only for "deforest"), "Mh": MeV (default Mp, dbase.f:217) }
#include <stdexcept>
#include <string>

#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SimcParameters.h"
#include "simc/physics/proton/ProtonCrossSections.h"
#include "simc/stages/IPhysicsModel.h"

namespace simc {
namespace proton {

/// Fortran: physics_proton.f:sigep, called at event.f:1451-1452.
class SigEPModel : public IPhysicsModel {
public:
    std::string name() const override { return "sigep"; }
    void configure(const Json&, const SimcParameters&) override {}
    void init(RunSetup&) override {}
    double sigcc(const Event& vertex, const MainEvent&) override { return sigep(vertex); }
    double sigccRecon(const Event& recon, const MainEvent&) override { return sigep(recon); }
};

/// Fortran: physics_proton.f:deForest, called at event.f:1455-1456.
/// fixedFlag < -1 means "take the flag from the parameters / JSON" (plug-in "deforest").
class DeForestModel : public IPhysicsModel {
public:
    explicit DeForestModel(const char* name = "deforest", int fixedFlag = -99) : name_(name), fixed_(fixedFlag) {}
    std::string name() const override { return name_; }
    void configure(const Json& cfg, const SimcParameters& params) override {
        st_.deForest_flag = params.deForest_flag;
        if (fixed_ >= -1) {
            st_.deForest_flag = fixed_;
        } else if (cfg.is_object() && cfg.contains("deForest_flag")) {
            st_.deForest_flag = cfg.at("deForest_flag").get<int>();
        }
        if (st_.deForest_flag < -1 || st_.deForest_flag > 1)
            throw std::invalid_argument("deforest: deForest_flag must be -1, 0 or 1 (got " +
                                        std::to_string(st_.deForest_flag) + ")");
        // Fortran: dbase.f:217,224  Mh = Mp; Mh2 = Mh*Mh (the hadron is a proton for (e,e'p)).
        double Mh = constants::Mp;
        if (cfg.is_object() && cfg.contains("Mh")) Mh = cfg.at("Mh").get<double>();
        st_.Mh2 = Mh * Mh;
    }
    void init(RunSetup&) override {}
    void beginEvent(long nevent) override { st_.nevent = nevent; }
    double sigcc(const Event& vertex, const MainEvent&) override { return deForest(vertex, st_); }
    double sigccRecon(const Event& recon, const MainEvent&) override { return deForest(recon, st_); }

    const DeForestState& state() const { return st_; }

private:
    std::string name_;
    int fixed_;
    DeForestState st_;
};

struct DeForestCC1 : DeForestModel { DeForestCC1() : DeForestModel("deforest_cc1", 0) {} };
struct DeForestCC2 : DeForestModel { DeForestCC2() : DeForestModel("deforest_cc2", 1) {} };
struct DeForestCC1OnShell : DeForestModel { DeForestCC1OnShell() : DeForestModel("deforest_cc1_onshell", -1) {} };

}  // namespace proton
}  // namespace simc

SIMC_REGISTER(simc::IPhysicsModel, "sigep", simc::proton::SigEPModel)
SIMC_REGISTER(simc::IPhysicsModel, "deforest", simc::proton::DeForestModel)
SIMC_REGISTER(simc::IPhysicsModel, "deforest_cc1", simc::proton::DeForestCC1)
SIMC_REGISTER(simc::IPhysicsModel, "deforest_cc2", simc::proton::DeForestCC2)
SIMC_REGISTER(simc::IPhysicsModel, "deforest_cc1_onshell", simc::proton::DeForestCC1OnShell)
