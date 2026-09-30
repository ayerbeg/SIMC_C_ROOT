// src/setup/HallCTargetLosses.cpp — see include/simc/setup/HallCTargetLosses.h.
#include "simc/setup/HallCTargetLosses.h"

#include <sstream>

#include "simc/material/HallCCryoTarget.h"

namespace simc {

namespace {
void copyIn(const TargetInfo& t, material::TargetInfo& m) {
    m.A = t.A; m.Z = t.Z; m.N = t.N;
    m.mass_amu = t.mass_amu; m.M = t.M; m.mrec_amu = t.mrec_amu; m.Mrec = t.Mrec;
    m.rho = t.rho; m.thick = t.thick; m.angle = t.angle; m.abundancy = t.abundancy; m.length = t.length;
    m.zoffset = t.zoffset; m.X0 = t.X0; m.X0_cm = t.X0_cm; m.L1 = t.L1; m.L2 = t.L2;
    m.fr1 = t.fr1; m.fr2 = t.fr2; m.xoffset = t.xoffset; m.yoffset = t.yoffset;
    m.fr_pattern = t.fr_pattern; m.can = t.can;
    for (int i = 0; i < 3; ++i) {
        m.Eloss[i] = {t.Eloss[i].ave, t.Eloss[i].min, t.Eloss[i].max};
        m.teff[i] = {t.teff[i].ave, t.teff[i].min, t.teff[i].max};
        m.musc_max[i] = t.musc_max[i];
    }
    m.Coulomb = {t.Coulomb.ave, t.Coulomb.min, t.Coulomb.max};
    m.Coulomb_constant = t.Coulomb_constant;
    m.musc_nsig_max = t.musc_nsig_max;
    m.Mtar_struck = t.Mtar_struck; m.Mrec_struck = t.Mrec_struck;
}
void splitLines(const std::string& text, std::vector<std::string>& out) {
    std::istringstream is(text);
    std::string l;
    while (std::getline(is, l)) {
        const auto a = l.find_first_not_of(' ');
        if (a != std::string::npos) out.push_back(l.substr(a));
    }
}
}  // namespace

HallCTargetLosses::HallCTargetLosses(int electron_arm, int hadron_arm)
    : target_(std::make_unique<material::HallCCryoTarget>()) {
    target_->setArms(electron_arm, hadron_arm);
}
HallCTargetLosses::~HallCTargetLosses() = default;

// Fortran: init.f:target_init (l.1-88) via HallCCryoTarget::targetInit
void HallCTargetLosses::targetInit(const TargetInitArgs& a, TargetInfo& t, std::vector<std::string>& messages) {
    std::ostringstream log;
    target_->setLog(&log);
    copyIn(t, target_->targ());
    target_->targetInit(a.using_Eloss, a.using_Coulomb, a.the_cent, a.thp_cent, a.Pp_cent, a.Mh2, a.Ebeam,
                        a.Pe_cent);
    const material::TargetInfo& m = target_->targ();
    t.L1 = m.L1;
    t.L2 = m.L2;
    t.X0 = m.X0;
    t.X0_cm = m.X0_cm;
    for (int i = 0; i < 3; ++i) {
        t.Eloss[i].ave = m.Eloss[i].ave;
        t.teff[i].ave = m.teff[i].ave;
    }
    t.Coulomb = {m.Coulomb.ave, m.Coulomb.min, m.Coulomb.max};
    t.Coulomb_constant = m.Coulomb_constant;
    target_->setLog(nullptr);
    splitLines(log.str(), messages);
}

// Fortran: target.f:extreme_trip_thru_target (+ extreme_target_musc), called at init.f:245
void HallCTargetLosses::extremeLosses(const ExtremeTripArgs& x, TargetInfo& t) {
    using L = material::Limits;
    target_->setLog(nullptr);
    copyIn(t, target_->targ());
    target_->extremeTripThruTarget(x.Ebeam, L{x.the.min, x.the.max}, L{x.thp.min, x.thp.max},
                                   L{x.pe.min, x.pe.max}, L{x.pp.min, x.pp.max}, L{x.z.min, x.z.max}, x.m);
    const material::TargetInfo& m = target_->targ();
    for (int i = 0; i < 3; ++i) {
        t.Eloss[i].min = m.Eloss[i].min;
        t.Eloss[i].max = m.Eloss[i].max;
        t.teff[i].min = m.teff[i].min;
        t.teff[i].max = m.teff[i].max;
        t.musc_max[i] = m.musc_max[i];
    }
    t.musc_nsig_max = m.musc_nsig_max;
}

}  // namespace simc
