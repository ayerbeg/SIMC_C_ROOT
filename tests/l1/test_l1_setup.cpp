// tests/l1/test_l1_setup.cpp
//
// L1 tests of the P2.1 setup utilities against the Fortran harness tables:
//   physics_angles      reference/functions/physics_angles/       (event.f:1572-1616)
//   spectrometer_angles reference/functions/spectrometer_angles/  (event.f:1618-1647, NaN rows included)
//   loren               reference/functions/loren/                (loren.f)
//   target_init_nonloss reference/functions/target_init/  columns that need no energy loss
//                       (targ_thick, targ_length, L1, L2, X0, X0_cm, Coulomb_ave, Coulomb_const)
//                       computed through BuildRunSetup (dbase.f conversions + init.f:target_init via the default HallCCryoTarget provider).
//
// Usage: test_l1_setup <routine> <table_dir> [rel_tol] [data_dir (target_init_nonloss: *.theory files)]
// Exit: 0 pass, 1 fail, 77 table missing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

#include "L1Runner.h"
#include "simc/core/SimcConstants.h"
#include "simc/setup/KinematicsUtils.h"
#include "simc/setup/SetupBuilder.h"

using namespace simc;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

int countLines(std::istream& in, const std::string& needle) {
    int n = 0;
    std::string l;
    while (std::getline(in, l))
        if (l.find(needle) != std::string::npos) ++n;
    return n;
}

int runPhysicsAngles(const std::string& dir, double tol) {
    std::ostringstream log;
    int rc = simcl1::run("physics_angles (" + dir + ")", dir, {{"theta", tol}, {"phi", tol}},
                         [&](const Row& r, Outputs& o) {
                             double th = 0, ph = 0;
                             PhysicsAngles(r["theta0"], r["phi0"], r["dx"], r["dy"], th, ph, &log);
                             o["theta"] = th;
                             o["phi"] = ph;
                         });
    if (rc != 0) return rc;
    // The phi0 /= pi/2, 3pi/2 warning (event.f:1596-1599): same number of occurrences as driver.log.
    std::ifstream flog(dir + "/driver.log");
    if (!flog) return rc;
    std::istringstream cl(log.str());
    const std::string w = "theta,phi will be incorrect";
    const int nf = countLines(flog, w), nc = countLines(cl, w);
    std::printf("  warning '%s': Fortran %d, C++ %d -> %s\n", w.c_str(), nf, nc, nf == nc ? "OK" : "MISMATCH");
    return nf == nc ? 0 : 1;
}

int runSpectrometerAngles(const std::string& dir, double tol) {
    // dy is compared too: the Fortran NaN / ~3e-9 rows (dy = 0 ill-conditioning, D-012 F1) must be
    // reproduced bit for bit (NaN == NaN counts as equal).
    return simcl1::run("spectrometer_angles (" + dir + ")", dir, {{"dx", tol}, {"dy", tol}},
                       [&](const Row& r, Outputs& o) {
                           double dx = 0, dy = 0;
                           SpectrometerAngles(r["theta0"], r["phi0"], dx, dy, r["theta"], r["phi"], /*fortran_compat=*/true);
                           o["dx"] = dx;
                           o["dy"] = dy;
                       });
}

int runLoren(const std::string& dir, double tol) {
    return simcl1::run("loren (" + dir + ")", dir,
                       {{"ef1", tol}, {"pxf", tol}, {"pyf", tol}, {"pzf", tol}, {"pf1", tol}},
                       [&](const Row& r, Outputs& o) {
                           const LorenResult x =
                               Loren(r["gam"], r["bx"], r["by"], r["bz"], r["e"], r["x"], r["y"], r["z"]);
                           o["ef1"] = x.ef1;
                           o["pxf"] = x.pxf;
                           o["pyf"] = x.pyf;
                           o["pzf"] = x.pzf;
                           o["pf1"] = x.pf1;
                       });
}

// target_init non energy-loss part, through the full BuildRunSetup (harness tstate.f applies exactly
// the dbase.f conversions, then calls target_init).
int runTargetInitNonLoss(const std::string& dir, double tol, const std::string& data_dir) {
    const std::vector<ColumnSpec> cols = {{"targ_thick", tol}, {"targ_length", tol}, {"L1", tol},
                                          {"L2", tol},         {"X0", tol},          {"X0_cm", tol},
                                          {"Coulomb_ave", tol}, {"Coulomb_const", tol}};
    return simcl1::run("target_init non-eloss (" + dir + ")", dir, cols, [&](const Row& r, Outputs& o) {
        SimcParameters p;
        p.targ.A = r["targ_A"];
        p.targ.Z = r["targ_Z"];
        p.targ.rho = r["targ_rho"];
        p.targ.thick = r["targ_thick_mg"];
        p.targ.angle = r["targ_angle_deg"];
        p.targ.can = static_cast<int>(std::lround(r["targ_can"]));
        p.targ.mass_amu = r["targ_A"];   // not used by these columns
        p.targ.abundancy = 100.;
        p.electron_arm = static_cast<int>(std::lround(r["electron_arm"]));
        p.hadron_arm = static_cast<int>(std::lround(r["hadron_arm"]));
        p.Ebeam = r["Ebeam"];
        p.spec.e.theta = r["the_cent_deg"];
        p.spec.p.theta = r["thp_cent_deg"];
        p.spec.e.P = r["Pe_cent"];
        p.spec.p.P = r["Pp_cent"];
        p.using_Eloss = r["using_Eloss"] != 0;
        p.using_Coulomb = r["using_Coulomb"] != 0;
        SetupOptions opt;
        opt.fortran_compat = true;   // D-031: Fortran references -> exact Fortran behaviour
        opt.data_dir = data_dir;   // theory files for the C12 rows (dbase.f:555-558)
        const RunSetup s = BuildRunSetup(p, opt);
        if (!s.ok()) std::printf("  row %zu: setup stopped: %s\n", r.number(), s.stops.front().c_str());
        o["targ_thick"] = s.targ.thick;
        o["targ_length"] = s.targ.length;
        o["L1"] = s.targ.L1;
        o["L2"] = s.targ.L2;
        o["X0"] = s.targ.X0;
        o["X0_cm"] = s.targ.X0_cm;
        o["Coulomb_ave"] = s.targ.Coulomb.ave;
        o["Coulomb_const"] = s.targ.Coulomb_constant;
    });
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s physics_angles|spectrometer_angles|loren|target_init_nonloss <dir> [tol]\n",
                     argv[0]);
        return 1;
    }
    const std::string what = argv[1], dir = argv[2];
    const double tol = argc > 3 ? std::atof(argv[3]) : 1e-12;
    if (what == "physics_angles") return runPhysicsAngles(dir, tol);
    if (what == "spectrometer_angles") return runSpectrometerAngles(dir, tol);
    if (what == "loren") return runLoren(dir, tol);
    if (what == "target_init_nonloss") return runTargetInitNonLoss(dir, tol, argc > 4 ? argv[4] : "");
    std::fprintf(stderr, "unknown routine %s\n", what.c_str());
    return 1;
}
