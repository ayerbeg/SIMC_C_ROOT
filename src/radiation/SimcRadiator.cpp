// src/radiation/SimcRadiator.cpp — see include/simc/radiation/SimcRadiator.h.
#include "simc/radiation/SimcRadiator.h"

#include <iostream>

#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"

namespace simc {
namespace radiation {

SimcRadiator::SimcRadiator() : log_(&std::cout) {}

void SimcRadiator::configure(const Json& config, const SimcParameters&) {
    if (config.is_object() && config.contains("verbose") && !config.at("verbose").get<bool>()) log_ = nullptr;
}

void SimcRadiator::init(RunSetup& setup) {
    run_ = RadcRun::FromSetup(setup);
    st_ = RadcState{};   // zero-initialised COMMON and SAVEd locals, as at program start
}

void SimcRadiator::initEvent(MainEvent& main, const Event& vertex) { radcInitEv(run_, st_, main, vertex, log_); }

bool SimcRadiator::radiate(MainEvent& main, Event& vertex, Event& orig, IRandom& rng, const CompleteEvFn& complete_ev) {
    return generateRad(run_, st_, main, vertex, orig, rng, complete_ev, log_);
}

double SimcRadiator::peakedRadWeight(const Event& vertex, double Egamma, double emin, double emax,
                                     double basicrad_val_reciprocal, double basicrad_weight) {
    return radiation::peakedRadWeight(run_, st_, vertex, Egamma, emin, emax, basicrad_val_reciprocal, basicrad_weight,
                                      log_);
}

SIMC_REGISTER(IRadiator, "simc_rad", SimcRadiator)

}  // namespace radiation
}  // namespace simc
