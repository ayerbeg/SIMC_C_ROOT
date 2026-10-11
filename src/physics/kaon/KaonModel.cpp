// src/physics/kaon/KaonModel.cpp — IPhysicsModel "kaon" (task P10.1). See include/simc/physics/kaon/KaonModel.h.
#include "simc/physics/kaon/KaonModel.h"

#include <iostream>
#include <stdexcept>

#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"

#ifndef SIMC_KAON_DATA_DIR
#define SIMC_KAON_DATA_DIR "data/kaon"
#endif

namespace simc {
namespace kaon {

void KaonModel::configure(const Json& config, const SimcParameters&) {
    data_dir_ = SIMC_KAON_DATA_DIR;
    if (config.is_object()) {
        if (config.contains("data_dir")) data_dir_ = config.at("data_dir").get<std::string>();
        if (config.contains("verbose")) verbose_ = config.at("verbose").get<bool>();
    }
}

void KaonModel::init(const StageContext& ctx) {
    const RunSetup& s = ctx.setup;
    cfg_.Mtar_struck = s.targ.Mtar_struck;
    cfg_.Mrec_struck = s.targ.Mrec_struck;
    cfg_.doing_decay = s.params.doing_decay;
    cfg_.hadron_arm = s.params.hadron_arm;
    cfg_.ctau = s.params.ctau;
    cfg_.spec_p_P = s.spec.p.P;
    cfg_.Mh2 = s.Mh2;
    const bool corr = !s.fortran_compat;
    cfg_.survival_p_kaon = corr;        // D-043 F5
    cfg_.survival_path_fixed = corr;    // D-043 F6
    cfg_.saghai.fint_value = corr;      // D-043 F1
    cfg_.saghai.proton_grid_fixed = corr;   // D-043 F3
    if (!(cfg_.Mtar_struck > 0))
        throw std::invalid_argument("kaon: targ%Mtar_struck is not set (the \"kaon\" cross section needs a kaon set-up, "
                                    "doing_kaon)");
    std::string err;
    std::vector<std::string> msgs;
    if (!tab_.load(data_dir_, corr, &msgs, &err))   // D-043 F2 (saghai_sigma0.dat records) in corrected mode
        throw std::invalid_argument("kaon: cannot read the kaon model tables (dbase.f:624-679): " + err +
                                    " (JSON key \"data_dir\" of the physics model)");
    if (verbose_)
        for (const auto& m : msgs) std::cout << " kaon: " << m << "\n";
}

double KaonModel::sigcc(const Event& vertex, const MainEvent& main) {
    MainEvent m = main;
    PeekState st = st_;
    double sp = survivalprob_;
    return peeK(vertex, m, sp, cfg_, tab_, st, verbose_ ? &std::cout : nullptr);
}

double KaonModel::sigccMain(const Event& vertex, MainEvent& main) {
    return peeK(vertex, main, survivalprob_, cfg_, tab_, st_, verbose_ ? &std::cout : nullptr);
}

}  // namespace kaon
}  // namespace simc

SIMC_REGISTER(simc::IPhysicsModel, "kaon", simc::kaon::KaonModel)
