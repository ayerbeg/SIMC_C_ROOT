// tests/l1/test_l1_physics_proton.cpp
//
// L1 tests of the proton cross sections (src/physics/proton) against the Fortran harness
// tables reference/functions/{sigmott,fofa_best_fit,sigep,deforest}/ (D-010: |rel| <= 1e-12).
//
// Usage: test_l1_physics_proton <routine> <table_dir> [rel_tol]
//   routine: sigmott | fofa_best_fit | sigep | deforest | plugins
// Exit: 0 pass, 1 fail, 77 table missing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "simc/core/Registry.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SimcParameters.h"
#include "simc/physics/proton/ProtonCrossSections.h"
#include "simc/stages/IPhysicsModel.h"

using namespace simc;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

// Fill the type(event) fields read by deForest from the table columns (README.txt of deforest).
Event deforestEvent(const Row& r) {
    Event ev;
    ev.Q2 = r["Q2"];
    ev.q = r["q"];
    ev.nu = r["nu"];
    ev.Pm = r["Pm"];
    ev.Pmx = r["Pmx"];
    ev.Pmy = r["Pmy"];
    ev.Pmz = r["Pmz"];
    ev.e.E = r["e_E"];
    ev.e.theta = r["e_theta"];
    ev.p.E = r["p_E"];
    ev.p.P = r["p_P"];
    ev.uq = {r["uq_x"], r["uq_y"], r["uq_z"]};
    ev.up = {r["up_x"], r["up_y"], r["up_z"]};
    return ev;
}

std::vector<std::string> warningLines(std::istream& in) {
    std::vector<std::string> v;
    std::string l;
    while (std::getline(in, l))
        if (l.find("WARNING") != std::string::npos) v.push_back(l);
    return v;
}

// Harness state hstate_eep_h: Mh = Mp, Mh2 = Mh*Mh.
constexpr double kMh2 = constants::Mp * constants::Mp;

int runDeforest(const std::string& dir, double tol) {
    std::ostringstream log;
    int rc = simcl1::run("deforest (" + dir + ")", dir, {{"deforest", tol}}, [&](const Row& r, Outputs& o) {
        proton::DeForestState st;
        st.deForest_flag = static_cast<int>(std::lround(r["deForest_flag"]));   // driver uses nint
        st.Mh2 = kMh2;
        st.nevent = static_cast<long>(r.number());   // driver: nevent = row number
        st.log = &log;
        o["deforest"] = proton::deForest(deforestEvent(r), st);
    });
    if (rc != 0) return rc;
    // The sin_gamma<0 clamp (physics_proton.f:74-78): same WARNING lines as the Fortran driver.log.
    std::ifstream flog(dir + "/driver.log");
    if (!flog) {
        std::printf("  driver.log missing: warning comparison skipped\n");
        return rc;
    }
    std::istringstream cl(log.str());
    auto fw = warningLines(flog), cw = warningLines(cl);
    std::size_t same = 0;
    for (std::size_t i = 0; i < fw.size() && i < cw.size(); ++i) same += (fw[i] == cw[i]);
    std::printf("  WARNING lines: Fortran %zu, C++ %zu, identical text %zu\n", fw.size(), cw.size(), same);
    if (fw.size() != cw.size() || same != fw.size()) {
        for (std::size_t i = 0; i < fw.size() || i < cw.size(); ++i) {
            std::printf("   F: %s\n   C: %s\n", i < fw.size() ? fw[i].c_str() : "-", i < cw.size() ? cw[i].c_str() : "-");
            if (i > 5) break;
        }
        std::printf("  => FAIL (warnings)\n");
        return 1;
    }
    return 0;
}

// The registered plug-ins give the free-function values (sigcc and sigccRecon) and pick the flag.
int runPlugins() {
    int fail = 0;
    auto expect = [&](bool ok, const char* what) {
        if (!ok) {
            std::printf("FAIL: %s\n", what);
            ++fail;
        }
    };
    Event ev;
    ev.Ein = 10600.;
    ev.e.E = 8360.989;
    ev.e.theta = 12.5 / constants::degrad;
    ev.Q2 = 2. * ev.Ein * ev.e.E * (1. - std::cos(ev.e.theta));
    ev.q = 2500.;
    ev.nu = ev.Ein - ev.e.E;
    ev.Pm = 120.;
    ev.Pmx = 50.;
    ev.Pmy = 60.;
    ev.Pmz = 90.;
    ev.p.E = 3000.;
    ev.p.P = 2850.;
    ev.uq = {0.1, -0.2, std::sqrt(1. - 0.05)};
    ev.up = {0.12, -0.18, std::sqrt(1. - 0.0144 - 0.0324)};
    MainEvent main;
    SimcParameters params;

    auto sig = Registry<IPhysicsModel>::create("sigep");
    sig->configure(Json::object(), params);
    expect(sig->sigcc(ev, main) == proton::sigep(ev), "sigep plug-in sigcc == sigep()");
    expect(sig->sigccRecon(ev, main) == proton::sigep(ev), "sigep plug-in sigccRecon == sigep()");

    struct Case { const char* name; int paramFlag; Json cfg; int expectFlag; };
    std::vector<Case> cases = {{"deforest", 0, Json::object(), 0},
                               {"deforest", 1, Json::object(), 1},
                               {"deforest", -1, Json::object(), -1},
                               {"deforest", 0, Json{{"deForest_flag", 1}}, 1},
                               {"deforest_cc1", 1, Json::object(), 0},
                               {"deforest_cc2", 0, Json::object(), 1},
                               {"deforest_cc1_onshell", 0, Json::object(), -1}};
    for (const auto& c : cases) {
        params.deForest_flag = c.paramFlag;
        auto m = Registry<IPhysicsModel>::create(c.name);
        m->configure(c.cfg, params);
        m->beginEvent(7);
        proton::DeForestState st;
        st.deForest_flag = c.expectFlag;
        st.Mh2 = kMh2;
        double ref = proton::deForest(ev, st);
        std::string what = std::string(c.name) + " flag " + std::to_string(c.expectFlag);
        expect(m->sigcc(ev, main) == ref, (what + " sigcc").c_str());
        expect(m->sigccRecon(ev, main) == ref, (what + " sigccRecon").c_str());
        expect(m->name() == c.name, (what + " name").c_str());
    }
    bool threw = false;
    try {
        auto m = Registry<IPhysicsModel>::create("deforest");
        m->configure(Json{{"deForest_flag", 2}}, params);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    expect(threw, "deforest rejects deForest_flag=2");
    std::printf("plugins: %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s sigmott|fofa_best_fit|sigep|deforest|plugins <table_dir> [rel_tol]\n", argv[0]);
        return 2;
    }
    const std::string routine = argv[1];
    if (routine == "plugins") return runPlugins();
    if (argc < 3) return 2;
    const std::string dir = argv[2];
    const double tol = argc > 3 ? std::atof(argv[3]) : 1e-12;

    if (routine == "sigmott")
        return simcl1::run("sigmott (" + dir + ")", dir, {{"sigmott", tol}}, [](const Row& r, Outputs& o) {
            o["sigmott"] = proton::sigMott(r["e0"], r["theta"], r["Q2"]);
        });
    if (routine == "fofa_best_fit")
        return simcl1::run("fofa_best_fit (" + dir + ")", dir, {{"GE", tol}, {"GM", tol}}, [](const Row& r, Outputs& o) {
            double GE = 0, GM = 0;
            proton::fofa_best_fit(r["qsquar"], GE, GM);
            o["GE"] = GE;
            o["GM"] = GM;
        });
    if (routine == "sigep")
        return simcl1::run("sigep (" + dir + ")", dir, {{"sigep", tol}}, [](const Row& r, Outputs& o) {
            Event v;
            v.Ein = r["Ein"];
            v.e.E = r["e_E"];
            v.e.theta = r["e_theta"];
            v.Q2 = r["Q2"];
            o["sigep"] = proton::sigep(v);
        });
    if (routine == "deforest") return runDeforest(dir, tol);
    std::fprintf(stderr, "unknown routine %s\n", routine.c_str());
    return 2;
}
