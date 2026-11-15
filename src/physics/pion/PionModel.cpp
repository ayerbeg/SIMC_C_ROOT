// src/physics/pion/PionModel.cpp — IPhysicsModel "pion". See include/simc/physics/pion/PionModel.h.
// The pion-electroproduction cross section of SIMC as a physics-model plug-in: physics_pion.f:peepi
// (PionCrossSections.cpp) with the model chosen by which_pion and the MAID tables read from the data tree. The
// reaction (reactions/pion/PionProduction.cpp) calls it in complete_main.
#include "simc/physics/pion/PionModel.h"

#include "simc/core/DataPath.h"

#include <iostream>
#include <stdexcept>

#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"


namespace simc {
namespace pion {

void PionModel::configure(const Json& config, const SimcParameters&) {
    data_dir_ = simc::DataDir("pion");
    if (config.is_object()) {
        if (config.contains("data_dir")) data_dir_ = config.at("data_dir").get<std::string>();
        if (config.contains("verbose")) verbose_ = config.at("verbose").get<bool>();
    }
}

// The model reads its run constants from the set-up and loads the MAID tables once.
void PionModel::init(const StageContext& ctx) {
    const RunSetup& s = ctx.setup;
    cfg_.which_pion = s.params.which_pion;
    cfg_.doing_pizero = s.flags.doing_pizero;
    cfg_.Mtar_struck = s.targ.Mtar_struck;
    cfg_.fortran_compat = s.fortran_compat;
    if (!(cfg_.Mtar_struck > 0))
        throw std::invalid_argument("pion: targ%Mtar_struck is not set (the \"pion\" cross section needs a pion set-up, "
                                    "doing_pion)");
    std::string err;
    if (!maid_.load(data_dir_, &err))
        throw std::invalid_argument("pion: cannot read the MAID tables (physics_pion.f:sigmaid): " + err +
                                    " (JSON key \"data_dir\" of the physics model)");
}

// sigcc: the side-effect-free evaluation (IPhysicsModel contract): the Fortran routine writes event fields and SAVEd
// state, so it runs on copies here. sigccMain is the complete_main call that keeps those writes.
double PionModel::sigcc(const Event& vertex, const MainEvent& main) {
    MainEvent m = main;
    return peepi(vertex, m, cfg_, maid_, verbose_ ? &std::cout : nullptr);
}

double PionModel::sigccMain(const Event& vertex, MainEvent& main) {
    return peepi(vertex, main, cfg_, maid_, verbose_ ? &std::cout : nullptr);
}

}  // namespace pion
}  // namespace simc

SIMC_REGISTER(simc::IPhysicsModel, "pion", simc::pion::PionModel)
