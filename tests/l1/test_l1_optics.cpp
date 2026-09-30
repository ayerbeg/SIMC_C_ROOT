// tests/l1/test_l1_optics.cpp
//
// L1 tests of the COSY optics engine (task P6.1) against the Fortran harness tables
// reference/functions/{transp_hms,transp_shms,transp_hms_chain,transp_shms_chain,recon_hms,recon_shms}/.
//
// Usage: test_l1_optics <mode> <matrix_dir> <table_dir> [rel_tol]
//   transp  : one class applied alone (both spectrometers are loaded, as the Fortran driver does)
//   chain   : aperture-free forward chain (every step) + pure COSY chain
//   recon   : mc_hms_recon / mc_shms_recon
//   files <matrix_dir> <simc_dir> : the repository copies of the COSY files equal the Fortran ones
//   selftest <matrix_dir>         : parser details, powi, class counts, decay hook
// <matrix_dir> holds hms/forward_cosy.dat, hms/recon_cosy.dat, shms/shms_forward.dat, shms/shms_recon.dat
// (the default relative paths of transp_init / mc_*_recon).
// Exit: 0 pass, 1 fail, 77 table missing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "simc/core/random/ScriptedRandom.h"
#include "simc/optics/CosyEngine.h"
#include "simc/optics/ForwardChains.h"

using namespace simc;
using namespace simc::optics;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

double g_tol = 1e-12;

bool loadEngine(CosyEngine& eng, const std::string& mdir, bool recon) {
    try {
        int nh = eng.initForward(1, mdir);
        int ns = eng.initForward(5, mdir);
        std::printf("  n_classes HMS = %d, SHMS = %d\n", nh, ns);
        if (nh != 12 || ns != 32) {   // mc_hms.f:188, mc_shms.f:344
            std::printf("  wrong number of transport classes\n");
            return false;
        }
        if (recon) {
            int rh = eng.initRecon(1, mdir), rs = eng.initRecon(5, mdir);
            std::printf("  recon terms HMS = %d (max order %d), SHMS = %d (max order %d)\n", rh, eng.reconMaxOrder(1), rs,
                        eng.reconMaxOrder(5));
        }
    } catch (const std::exception& e) {
        std::printf("  cannot load COSY files from %s: %s\n", mdir.c_str(), e.what());
        return false;
    }
    return true;
}

int runTransp(const std::string& mdir, const std::string& dir) {
    CosyEngine eng;
    std::ostringstream log;
    eng.setLog(&log);
    if (!loadEngine(eng, mdir, false)) return 1;
    std::vector<ColumnSpec> c;
    for (const char* n : {"xs_out", "dxdzs_out", "ys_out", "dydzs_out", "dpps_out", "pathlen", "delta_z", "driftdist"})
        c.push_back({n, g_tol, 0.0});
    c.push_back({"adrift", 0.0, 0.0});
    int rc = simcl1::run("transp (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        const int spectr = static_cast<int>(std::lround(r["spectr"]));
        const int cls = static_cast<int>(std::lround(r["class"]));
        const double zd = r["zd"];
        Track t;
        t.xs = r["xs"];
        t.dxdzs = r["dxdzs"];
        t.ys = r["ys"];
        t.dydzs = r["dydzs"];
        t.dpps = r["dpps"];
        double pathlen = 0.;
        eng.transp(spectr, cls, t, zd, pathlen, nullptr);
        o["xs_out"] = t.xs;
        o["dxdzs_out"] = t.dxdzs;
        o["ys_out"] = t.ys;
        o["dydzs_out"] = t.dydzs;
        o["dpps_out"] = t.dpps;
        o["pathlen"] = pathlen;
        o["delta_z"] = pathlen - zd;   // drv_transp.f:66
        o["adrift"] = eng.adrift(spectr, cls) ? 1. : 0.;
        o["driftdist"] = eng.driftdist(spectr, cls);
    });
    if (!log.str().empty()) {
        std::printf("  unexpected warnings (the Fortran driver.log has none):\n%s", log.str().c_str());
        return 1;
    }
    return rc;
}

int runChain(const std::string& mdir, const std::string& dir) {
    CosyEngine eng;
    std::ostringstream log;
    eng.setLog(&log);
    if (!loadEngine(eng, mdir, false)) return 1;
    simcl1::Table ref;
    std::string err;
    if (!ref.load(dir + "/output.txt", err)) {
        std::printf("chain: SKIP (%s)\n", err.c_str());
        return 77;
    }
    // Output columns: every <step>_* column plus cosy_*.
    std::vector<ColumnSpec> c;
    for (std::size_t i = 6; i < ref.names.size(); ++i) c.push_back({ref.names[i], g_tol, 0.0});
    int rc = simcl1::run("chain (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        const int spectr = static_cast<int>(std::lround(r["spectr"]));
        Track t0;
        t0.xs = r["x"];
        t0.ys = r["y"];
        t0.dxdzs = r["dxdz"];
        t0.dydzs = r["dydz"];
        t0.dpps = r["dpp"];
        auto steps = spectr == 1 ? collectHmsChain(eng, t0) : collectShmsChain(eng, t0);
        for (const auto& s : steps) {
            o[s.name + "_x"] = s.track.xs;
            o[s.name + "_y"] = s.track.ys;
            o[s.name + "_dxdz"] = s.track.dxdzs;
            o[s.name + "_dydz"] = s.track.dydzs;
            o[s.name + "_path"] = s.path;
            o[s.name + "_xt"] = s.xt;
            o[s.name + "_yt"] = s.yt;
        }
        Track t = t0;
        double path = 0.;
        cosyChain(eng, spectr, t, path);
        o["cosy_x"] = t.xs;
        o["cosy_dxdz"] = t.dxdzs;
        o["cosy_y"] = t.ys;
        o["cosy_dydz"] = t.dydzs;
        o["cosy_path"] = path;
    });
    if (!log.str().empty()) {
        std::printf("  unexpected warnings:\n%s", log.str().c_str());
        return 1;
    }
    return rc;
}

int runRecon(const std::string& mdir, const std::string& dir) {
    CosyEngine eng;
    if (!loadEngine(eng, mdir, true)) return 1;
    std::vector<ColumnSpec> c;
    for (const char* n : {"dpp_recon", "dph_recon", "dth_recon", "y_recon"}) c.push_back({n, g_tol, 0.0});
    return simcl1::run("recon (" + dir + ")", dir, c, [&](const Row& r, Outputs& o) {
        const int spectr = static_cast<int>(std::lround(r["spectr"]));
        Track t;
        t.xs = r["xfp"];
        t.dxdzs = r["xpfp"];
        t.ys = r["yfp"];
        t.dydzs = r["ypfp"];
        double dp, dt, dph, yt;
        eng.recon(spectr, t, r["fry"], dp, dt, dph, yt);
        o["dpp_recon"] = dp;
        o["dph_recon"] = dph;
        o["dth_recon"] = dt;
        o["y_recon"] = yt;
    });
}

std::string slurp(const std::string& p, bool& ok) {
    std::ifstream f(p, std::ios::binary);
    ok = static_cast<bool>(f);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

int runFiles(const std::string& mdir, const std::string& simc) {
    int fail = 0;
    for (const char* f : {"hms/forward_cosy.dat", "hms/recon_cosy.dat", "shms/shms_forward.dat", "shms/shms_recon.dat"}) {
        bool a, b;
        std::string x = slurp(mdir + "/" + f, a), y = slurp(simc + "/" + f, b);
        if (!b) {
            std::printf("  %s: Fortran copy missing, SKIP\n", f);
            return 77;
        }
        bool same = a && x == y;
        std::printf("  %-24s %zu bytes  %s\n", f, x.size(), same ? "identical" : "DIFFERENT");
        if (!same) ++fail;
    }
    return fail ? 1 : 0;
}

int runSelftest(const std::string& mdir) {
    int fail = 0;
    auto expect = [&](bool ok, const std::string& what) {
        std::printf("  %-72s %s\n", what.c_str(), ok ? "ok" : "FAIL");
        if (!ok) ++fail;
    };
    // Fortran Gw.d input editing
    expect(detail::readFortranReal("  0.1234567E-01", 7) == 0.1234567e-01, "G14.7 with exponent");
    expect(detail::readFortranReal("   1234567", 7) == 0.1234567, "G14.7 without decimal point: implied d");
    expect(detail::readFortranReal(" -0.1234567-100", 7) == -0.1234567e-100, "exponent without letter");
    expect(detail::readFortranReal("", 7) == 0.0 && detail::readFortranInt(" ") == 0, "blank fields read as 0");
    expect(detail::readFortranReal(" 1 . 5 ", 7) == 1.5, "embedded blanks ignored (BLANK='NULL')");
    // powi = libgcc __powidf2
    expect(powi(3.0, 0) == 1.0 && powi(0.0, 0) == 1.0, "powi(x,0) = 1");
    {
        double x = 1.1, y = x;           // __powidf2(1.1, 5): y=x; x=x*x (n=2); x=x*x (n=1) -> y*=x
        double x2 = x * x, x4 = x2 * x2;
        y = y * x4;
        expect(powi(1.1, 5) == y, "powi(1.1,5) = x*(x^2)^2 (repeated squaring)");
        expect(powi(2.0, -2) == 0.25, "powi negative exponent");
    }
    CosyEngine eng;
    std::ostringstream log;
    eng.setLog(&log);
    bool loaded = loadEngine(eng, mdir, true);
    expect(loaded, "COSY files load (HMS 12 classes, SHMS 32)");
    if (loaded) {
        expect(eng.reconTerms(1).size() == 461 && eng.reconTerms(5).size() == 461, "recon: 461 terms each (P6.F1)");
        int nd = 0;
        for (int k = 1; k <= 12; ++k) nd += eng.adrift(1, k);
        expect(nd == 5, "HMS: classes 1,4,7,10,12 are drifts");
        expect(eng.driftdist(1, 12) == 625.0, "HMS class 12 drift = 625.0 cm");
        expect(log.str().empty(), "no transp_init warnings");
        // decay hook: a pion with ctau tiny decays in the first project(); uses the RNG
        ScriptedRandom rng;
        rng.load({}, {0.5, 0.25, 0.5});
        DecayState d;
        d.decay_flag = true;
        d.m2 = 139.57018 * 139.57018;
        d.ph = 1000.;
        d.ctau = 1e-3;
        d.rng = &rng;
        Track t;
        t.dpps = 0.;
        double path = 0.;
        project(t, 100., path, &d);
        expect(d.dflag && rng.nUniform() == 3 && std::fabs(d.Mh2_final - 105.6583755 * 105.6583755) < 1e-9,
               "project: decay branch (pi -> mu), 3 grnd, Mh2_final = Mmu^2");
        bool thr = false;
        try {
            CosyEngine::defaultForwardFile(6);
        } catch (const FortranStop&) {
            thr = true;
        }
        expect(thr, "spectr 6 (SHMS LSA) -> FortranStop");
    }
    std::printf("selftest: %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <mode> <matrix_dir> [table_dir|simc_dir] [rel_tol]\n", argv[0]);
        return 2;
    }
    const std::string mode = argv[1], mdir = argv[2];
    if (mode == "selftest") return runSelftest(mdir);
    if (argc < 4) return 2;
    const std::string dir = argv[3];
    if (argc > 4) g_tol = std::atof(argv[4]);
    if (mode == "transp") return runTransp(mdir, dir);
    if (mode == "chain") return runChain(mdir, dir);
    if (mode == "recon") return runRecon(mdir, dir);
    if (mode == "files") return runFiles(mdir, dir);
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 2;
}
