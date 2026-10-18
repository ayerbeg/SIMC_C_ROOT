// src/physics/sidis/SidisModel.cpp — IPhysicsModel "sidis" (task P12.1). See include/simc/physics/sidis/SidisModel.h.
#include "simc/physics/sidis/SidisModel.h"

#include <iostream>
#include <stdexcept>

#include "simc/core/Errors.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcParameters.h"

#ifndef SIMC_SIDIS_DATA_DIR
#define SIMC_SIDIS_DATA_DIR "data/sidis"
#endif

namespace simc {
namespace sidis {

std::vector<std::string> checkSupported(const RunSetup& s) {
    std::vector<std::string> v;
    const ReactionFlags& f = s.flags;
    if (!f.doing_semi) v.push_back("sidis: doing_semi is not set (the \"sidis\" cross section is semi_physics.f:peepiX)");
    if (f.doing_semi && f.doing_semipi == f.doing_semika)
        v.push_back("sidis: doing_semi needs doing_pion = 1 (pi+-, (e,e'pi)X) or doing_kaon = 1 (K+-, (e,e'K)X) "
                    "(dbase.f:125-132; the Fortran STOPs with 'I don't have ANY idea what (e,e'm)X we're doing!!!')");
    if (f.doing_semi && f.doing_pizero)
        v.push_back("sidis: doing_pizero with doing_semi (pi0 SIDIS) is not supported yet: the pi0 photons are detected "
                    "with the calorimeter arm (task P14)");
    if (f.doing_semi && f.doing_phsp)
        v.push_back("sidis: doing_phsp with doing_semi: the Fortran weight of a phase-space run uses sigcc = 1 "
                    "(event.f:1446-1448), not peepiX; the \"sidis\" model does not apply");
    return v;
}

std::vector<std::string> ignoredInputNotes(const RunSetup& s) {
    std::vector<std::string> v;
    const SimcParameters& p = s.params;
    if (p.pt_b_param != 0.0)
        v.push_back("NOTE (SIMC++): pt_b_param is read but not used by semi_physics.f (l.142 is commented out; the pt2 "
                    "slope is b = 1/(0.120 z^2 + 0.200) GeV^-2, l.525)");
    if (p.sigc_flag != 0 || p.sigc_nbin != 0 || p.sigc_kin_min != 0.0 || p.sigc_kin_max != 0.0 || p.sigc_kin_ind != 0.0)
        v.push_back("NOTE (SIMC++): sigc_flag/sigc_nbin/sigc_kin_min/sigc_kin_max/sigc_kin_ind are read but not used: "
                    "the central cross section call peepiX(..., .TRUE.) is commented out (event.f:1523)");
    return v;
}

PeepiXConfig makeConfig(const RunSetup& s) {
    PeepiXConfig c;
    c.targA = s.targ.A;
    c.targZ = s.targ.Z;
    c.Mtar_struck = s.targ.Mtar_struck;
    c.doing_semipi = s.flags.doing_semipi;
    c.doing_semika = s.flags.doing_semika;
    c.doing_pizero = s.flags.doing_pizero;
    c.doing_hplus = s.flags.doing_hplus;
    c.do_fermi = s.flags.do_fermi;
    c.doing_decay = s.params.doing_decay;
    c.hadron_arm = s.params.hadron_arm;
    c.ctau = s.params.ctau;
    c.spec_p_P = s.spec.p.P;
    c.Mh2 = s.Mh2;
    c.sigc_flag = s.params.sigc_flag;
    c.sigc_nbin = s.params.sigc_nbin;
    c.sigc_kin_min = s.params.sigc_kin_min;
    c.sigc_kin_max = s.params.sigc_kin_max;
    c.sigc_kin_ind = s.params.sigc_kin_ind;
    const bool corr = !s.fortran_compat;
    c.pizero_db_fixed = corr;       // P12.1 F1
    c.survival_p_hadron = corr;     // P12.1 F3
    c.survival_path_fixed = corr;   // P12.1 F4
    return c;
}

void SidisModel::configure(const Json& config, const SimcParameters&) {
    data_dir_ = SIMC_SIDIS_DATA_DIR;
    if (config.is_object()) {
        if (config.contains("data_dir")) data_dir_ = config.at("data_dir").get<std::string>();
        if (config.contains("verbose")) verbose_ = config.at("verbose").get<bool>();
    }
}

void SidisModel::init(const StageContext& ctx) {
    const RunSetup& s = ctx.setup;
    const auto problems = checkSupported(s);
    if (!problems.empty()) {
        std::string m;
        for (const auto& p : problems) m += (m.empty() ? "" : "; ") + p;
        if (s.flags.doing_semi && s.flags.doing_pizero) throw NotPorted(m);
        throw std::invalid_argument(m);
    }
    cfg_ = makeConfig(s);
    if (!(cfg_.Mtar_struck > 0))
        throw std::invalid_argument("sidis: targ%Mtar_struck is not set (the \"sidis\" cross section needs a doing_semi "
                                    "set-up)");
    std::string err;
    if (!pdf_.load(data_dir_, 1, &err))   // SetCtq5(iset = 1), semi_physics.f:121, 300
        throw std::invalid_argument("sidis: cannot read the CTEQ5M table: " + err +
                                    " (JSON key \"data_dir\" of the physics model)");
    if (cfg_.doing_semika && !fdss_.load(data_dir_, 2, 1, &err))   // fDSS IH = 2, IO = 1 (semi_physics.f:497-505)
        throw std::invalid_argument("sidis: cannot read the fDSS kaon grid: " + err +
                                    " (JSON key \"data_dir\" of the physics model)");
    st_ = PeepiXState{};
    survivalprob_ = 0;
    if (verbose_)
        for (const auto& n : ignoredInputNotes(s)) std::cout << " sidis: " << n << "\n";
}

double SidisModel::sigcc(const Event& vertex, const MainEvent& main) {
    MainEvent m = main;
    PeepiXState st = st_;
    double sp = survivalprob_;
    const bool warn = pdf_.warn, first = pdf_.first;   // the SAVEd PDF warning flags are restored
    const double r = peepiX(vertex, vertex, m, sp, false, cfg_, pdf_, cfg_.doing_semika ? &fdss_ : nullptr, st,
                            verbose_ ? &std::cout : nullptr);
    pdf_.warn = warn;
    pdf_.first = first;
    return r;
}

double SidisModel::sigccMain(const Event& vertex, MainEvent& main) {
    return peepiX(vertex, vertex, main, survivalprob_, false, cfg_, pdf_, cfg_.doing_semika ? &fdss_ : nullptr, st_,
                  verbose_ ? &std::cout : nullptr);
}

}  // namespace sidis
}  // namespace simc

SIMC_REGISTER(simc::IPhysicsModel, "sidis", simc::sidis::SidisModel)
