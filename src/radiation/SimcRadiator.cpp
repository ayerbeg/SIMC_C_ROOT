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

void SimcRadiator::init(const StageContext& ctx) {
    run_ = RadcRun::FromSetup(ctx.setup);
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

// Fortran: simc.f:calculate_central l.1268-1281 (central%rad = the /radccom/ values after radc_init_ev)
RadiationValues SimcRadiator::currentValues() const {
    RadiationValues v;
    v.hardcorfac = st_.hardcorfac;
    v.etatzai = run_.etatzai;
    v.g_int = st_.g_int;
    v.g_ext = st_.g_ext;
    for (int i = 1; i <= 3; ++i) {   // Fortran index i -> [i-1]
        v.frac[i - 1] = st_.frac[i];
        v.lambda[i - 1] = st_.lambda[i];
        if (i < 3) v.bt[i - 1] = st_.bt[i];
    }
    for (int i = 0; i <= 3; ++i) {
        v.c_int[i] = st_.c_int[i];
        v.c_ext[i] = st_.c_ext[i];
        v.c[i] = st_.c[i];
        v.g[i] = st_.g[i];
    }
    return v;
}

std::array<double, 3> SimcRadiator::photonEnergies() const {
    return {st_.Egamma_used[1], st_.Egamma_used[2], st_.Egamma_used[3]};
}

bool SimcRadiator::hadronRadiates() const { return st_.rad_proton_this_ev; }

SIMC_REGISTER(IRadiator, "simc_rad", SimcRadiator)

}  // namespace radiation
}  // namespace simc
