// tests/l1/test_l1_trg.cpp
//
// L1 tests of the polarized-target field tracking optics::TargetField (task P15.2b, trg_track.f) against the Fortran
// harness tables reference/functions/trg_track/<table>/ (P15.F2). Tolerance 0 (bitwise) for every column.
//
// Usage: test_l1_trg <mode> <table_dir> <field_dir> <matrix_dir>
//   modes: trginit, trginit_map, trgfield, trgderiv, trgrk4, trgtracktoplane, track_from_tgt, track_to_tgt
//          files <field_dir> <simc_dir> : data/field/trg_field_map.dat equals the Fortran file
//          replay <capture_dir> <field_dir> <matrix_dir> : every trgTrackToPlane / track_from_tgt / track_to_tgt call
//                 of the capture reference/runs/p15f2_poltar_cap (records cap_plane/cap_from/cap_to.txt.xz, via xz -dc)
// <field_dir> holds trg_field_map.dat (map_cfg 1); <matrix_dir> the COSY recon files (hms/, sos/, hrs/, shms/) for
// the mc_*_recon calls of track_to_tgt (whose COMMON /track/ hut values are the track_* input columns).
// Exit: 0 pass, 1 fail, 77 table missing.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "simc/optics/CosyEngine.h"
#include "simc/optics/TargetField.h"

using namespace simc;
using namespace simc::optics;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

namespace {

std::string g_fieldDir;

/// The trgInit call of the Fortran driver for the configuration of a row (only when it changes).
struct FieldCfg {
    TargetField f;
    double key[5]{-1e300, 0, 0, 0, 0};
    void set(const Row& r) {
        const double k[5] = {r["map_cfg"], r["theta_e_deg"], r["phi_e_deg"], r["theta_p_deg"], r["phi_p_deg"]};
        bool same = true;
        for (int i = 0; i < 5; ++i) same = same && (std::memcmp(&k[i], &key[i], sizeof(double)) == 0);
        if (same) return;
        for (int i = 0; i < 5; ++i) key[i] = k[i];
        const int m = static_cast<int>(std::lround(k[0]));
        const std::string map = m == 1 ? g_fieldDir + "/trg_field_map.dat" : (m == 2 ? std::string("0") : std::string(" "));
        f.init(map, k[1], k[2], k[3], k[4]);
    }
};

std::vector<ColumnSpec> exact(std::initializer_list<const char*> names) {
    std::vector<ColumnSpec> c;
    for (const char* n : names) c.push_back({n, 0.0, 0.0});
    return c;
}

TargetField::Vec readU(const Row& r, const char* prefix = "u") {
    TargetField::Vec u{};
    for (int i = 0; i < 6; ++i) u[static_cast<std::size_t>(i)] = r[std::string(prefix) + std::to_string(i + 1)];
    return u;
}

int spectOf(const Row& r) { return static_cast<int>(std::lround(r["spect"])); }

int runTrginit(const std::string& dir) {
    FieldCfg F;
    return simcl1::run("trgInit", dir,
                       exact({"B_theta_e", "B_stheta_e", "B_ctheta_e", "B_phi_e", "B_sphi_e", "B_cphi_e", "B_theta_p",
                              "B_stheta_p", "B_ctheta_p", "B_phi_p", "B_sphi_p", "B_cphi_p", "sum_Bz", "sum_Br",
                              "sum_zz", "sum_rr"}),
                       [&](const Row& r, Outputs& o) {
                           F.key[0] = -1e300;
                           F.set(r);
                           const auto& e = F.f.anglesE();
                           const auto& p = F.f.anglesP();
                           o["B_theta_e"] = e.theta, o["B_stheta_e"] = e.stheta, o["B_ctheta_e"] = e.ctheta;
                           o["B_phi_e"] = e.phi, o["B_sphi_e"] = e.sphi, o["B_cphi_e"] = e.cphi;
                           o["B_theta_p"] = p.theta, o["B_stheta_p"] = p.stheta, o["B_ctheta_p"] = p.ctheta;
                           o["B_phi_p"] = p.phi, o["B_sphi_p"] = p.sphi, o["B_cphi_p"] = p.cphi;
                           double sz = 0, sr = 0, szz = 0, srr = 0;
                           for (int ir = 1; ir <= TargetField::nr; ++ir)   // column-major array order
                               for (int iz = 1; iz <= TargetField::nz; ++iz) {
                                   sz += F.f.Bz(iz, ir);
                                   sr += F.f.Br(iz, ir);
                               }
                           for (int i = 1; i <= TargetField::nz; ++i) szz += F.f.zz(i), srr += F.f.rr(i);
                           o["sum_Bz"] = sz, o["sum_Br"] = sr, o["sum_zz"] = szz, o["sum_rr"] = srr;
                       });
}

// output.txt has one row per map node (map_cfg iz ir zz rr B_field_z B_field_r): compared directly.
int runTrginitMap(const std::string& dir) {
    simcl1::Table out;
    std::string err;
    if (!out.load(dir + "/output.txt", err)) {
        std::printf("SKIP: %s\n", err.c_str());
        return 77;
    }
    TargetField f[3];
    f[0].init(g_fieldDir + "/trg_field_map.dat", 0, 0, 0, 0);
    f[1].init("0", 0, 0, 0, 0);
    f[2].init(" ", 0, 0, 0, 0);
    std::size_t bad = 0;
    for (std::size_t i = 0; i < out.rows.size(); ++i) {
        const Row r(out, i);
        const int m = static_cast<int>(std::lround(r["map_cfg"]));
        const int iz = static_cast<int>(std::lround(r["iz"])), ir = static_cast<int>(std::lround(r["ir"]));
        const TargetField& t = f[m - 1];
        const double v[4] = {t.zz(iz), t.rr(ir), t.Bz(iz, ir), t.Br(iz, ir)};
        const double ref[4] = {r["zz"], r["rr"], r["B_field_z"], r["B_field_r"]};
        for (int k = 0; k < 4; ++k)
            if (std::memcmp(&v[k], &ref[k], sizeof(double)) != 0) {
                if (bad < 10) std::printf("  row %zu col %d: C++ %.17g Fortran %.17g\n", i + 1, k, v[k], ref[k]);
                ++bad;
            }
    }
    std::printf("trgInit field arrays (%s): %zu rows x 4 values, %zu differing (bitwise)\n", dir.c_str(), out.rows.size(),
                bad);
    std::printf("%s\n", bad == 0 && !out.rows.empty() ? "PASS" : "FAIL");
    return bad == 0 && !out.rows.empty() ? 0 : 1;
}

int runTrgfield(const std::string& dir) {
    FieldCfg F;
    return simcl1::run("trgField", dir, exact({"B1", "B2", "B3"}), [&](const Row& r, Outputs& o) {
        F.set(r);
        const double x[3] = {r["x1"], r["x2"], r["x3"]};
        double B[3];
        F.f.field(x, B, spectOf(r));
        o["B1"] = B[0], o["B2"] = B[1], o["B3"] = B[2];
    });
}

int runTrgderiv(const std::string& dir) {
    FieldCfg F;
    return simcl1::run("trgDeriv", dir,
                       exact({"dudt1", "dudt2", "dudt3", "dudt4", "dudt5", "dudt6", "dudt7", "dudt8", "dudt9"}),
                       [&](const Row& r, Outputs& o) {
                           F.set(r);
                           F.f.setFactor(r["factor"]);
                           TargetField::Vec d{};
                           F.f.deriv(readU(r), d, spectOf(r));
                           for (int i = 0; i < 9; ++i) o["dudt" + std::to_string(i + 1)] = d[static_cast<std::size_t>(i)];
                       });
}

int runTrgrk4(const std::string& dir) {
    FieldCfg F;
    return simcl1::run("trgRK4", dir, exact({"u1_out", "u2_out", "u3_out", "u4_out", "u5_out", "u6_out"}),
                       [&](const Row& r, Outputs& o) {
                           F.set(r);
                           F.f.setFactor(r["factor"]);
                           TargetField::Vec u1{};
                           F.f.rk4(readU(r), u1, r["h"], spectOf(r));
                           for (int i = 0; i < 6; ++i)
                               o["u" + std::to_string(i + 1) + "_out"] = u1[static_cast<std::size_t>(i)];
                       });
}

int runTrackToPlane(const std::string& dir) {
    FieldCfg F;
    return simcl1::run("trgTrackToPlane", dir,
                       exact({"u1_out", "u2_out", "u3_out", "u4_out", "u5_out", "u6_out", "ok_out", "factor_out"}),
                       [&](const Row& r, Outputs& o) {
                           F.set(r);
                           F.f.setFactor(0.0);   // the driver's value before the call (T9)
                           TargetField::Vec u = readU(r);
                           bool ok = r["ok_in"] != 0.0;
                           F.f.trackToPlane(u, r["E"], r["dl"], r["a"], r["b"], r["c"], r["d"], ok, spectOf(r));
                           for (int i = 0; i < 6; ++i)
                               o["u" + std::to_string(i + 1) + "_out"] = u[static_cast<std::size_t>(i)];
                           o["ok_out"] = ok ? 1.0 : 0.0;
                           o["factor_out"] = F.f.factor();
                       });
}

int runTrackFrom(const std::string& dir) {
    FieldCfg F;
    return simcl1::run("track_from_tgt", dir,
                       exact({"x_out", "y_out", "z_out", "dx_out", "dy_out", "ok_out", "factor_out"}),
                       [&](const Row& r, Outputs& o) {
                           F.set(r);
                           F.f.setFactor(0.0);
                           double x = r["x"], y = r["y"], z = r["z"], dx = r["dx"], dy = r["dy"];
                           bool ok = false;
                           F.f.trackFromTgt(x, y, z, dx, dy, r["mom"], r["mass"], spectOf(r), ok);
                           o["x_out"] = x, o["y_out"] = y, o["z_out"] = z, o["dx_out"] = dx, o["dy_out"] = dy;
                           o["ok_out"] = ok ? 1.0 : 0.0;
                           o["factor_out"] = F.f.factor();
                       });
}

int runTrackTo(const std::string& dir, const std::string& mdir) {
    CosyEngine eng;
    eng.setFortranCompat(true);
    eng.setLog(nullptr);
    try {
        eng.initRecon(1, mdir);
        eng.initRecon(2, mdir);
        eng.loadRecon(3, mdir + "/hrs/hrs_recon_cosy.dat");
        eng.loadRecon(4, mdir + "/hrs/hrs_recon_cosy.dat");
        eng.initRecon(5, mdir);
    } catch (const std::exception& e) {
        std::printf("cannot load the recon matrices from %s: %s\n", mdir.c_str(), e.what());
        return 1;
    }
    FieldCfg F;
    return simcl1::run("track_to_tgt", dir,
                       exact({"delta_out", "y_out", "dx_out", "dy_out", "mom_out", "ok_out", "n_recon", "factor_out"}),
                       [&](const Row& r, Outputs& o) {
                           F.set(r);
                           F.f.setFactor(0.0);
                           Track t;
                           t.xs = r["track_xs"];
                           t.ys = r["track_ys"];
                           t.dxdzs = r["track_dxdzs"];
                           t.dydzs = r["track_dydzs"];
                           const int arm = static_cast<int>(std::lround(r["arm"]));
                           double delta = r["delta"], y = r["y"], dx = r["dx"], dy = r["dy"], mom = r["mom"];
                           bool ok = r["ok_in"] != 0.0;
                           int n = 0;
                           F.f.trackToTgt(delta, y, dx, dy, r["frx"], r["fry"], mom, r["mass"], r["ctheta"],
                                          r["stheta"], spectOf(r), ok, arm,
                                          [&](double fry, double& d, double& yp, double& xp, double& yt) {
                                              eng.recon(arm, t, fry, d, yp, xp, yt);
                                          },
                                          &n);
                           o["delta_out"] = delta, o["y_out"] = y, o["dx_out"] = dx, o["dy_out"] = dy;
                           o["mom_out"] = mom;
                           o["ok_out"] = ok ? 1.0 : 0.0;
                           o["n_recon"] = n;
                           o["factor_out"] = F.f.factor();
                       });
}

// ---- capture replay (reference/runs/p15f2_poltar_cap/capture: every call of the 1000-event capture, config A)
struct XzTable {
    std::vector<std::string> names;
    FILE* pipe{nullptr};
    std::vector<double> row;
    bool open(const std::string& path) {
        pipe = popen(("xz -dc '" + path + "'").c_str(), "r");
        if (!pipe) return false;
        char buf[8192];
        if (!fgets(buf, sizeof buf, pipe)) return false;
        std::istringstream h(buf);
        std::string n;
        while (h >> n) names.push_back(n);
        return !names.empty();
    }
    bool next() {
        static char buf[1 << 16];
        if (!fgets(buf, sizeof buf, pipe)) return false;
        row.clear();
        char* p = buf;
        char* e = nullptr;
        for (std::size_t i = 0; i < names.size(); ++i) {
            row.push_back(std::strtod(p, &e));
            p = e;
        }
        return true;
    }
    double operator[](const std::string& n) const {
        for (std::size_t i = 0; i < names.size(); ++i)
            if (names[i] == n) return row[i];
        throw std::out_of_range("column " + n);
    }
    ~XzTable() {
        if (pipe) pclose(pipe);
    }
};

bool same(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

int runReplay(const std::string& cap, const std::string& mdir) {
    TargetField f;
    f.init(g_fieldDir + "/trg_field_map.dat", 6.4749999999999986E+001, 0, 9.2699999809265137E+001, 0);   // config A
    CosyEngine eng;
    eng.setFortranCompat(true);
    eng.setLog(nullptr);
    eng.initRecon(1, mdir);
    eng.initRecon(5, mdir);
    int fails = 0;
    {   // every trgTrackToPlane call (factor_out compared when ok_in = 1: stale otherwise, T9)
        XzTable t;
        if (!t.open(cap + "/cap_plane.txt.xz")) {
            std::printf("SKIP: %s/cap_plane.txt.xz missing\n", cap.c_str());
            return 77;
        }
        long n = 0, bad = 0;
        while (t.next()) {
            TargetField::Vec u{};
            for (int i = 0; i < 6; ++i) u[static_cast<std::size_t>(i)] = t["u" + std::to_string(i + 1) + "_in"];
            bool ok = t["ok_in"] != 0;
            const bool ok_in = ok;
            f.trackToPlane(u, t["E"], t["dl"], t["a"], t["b"], t["c"], t["d"], ok, static_cast<int>(t["spect"]));
            bool d = (ok ? 1.0 : 0.0) != t["ok_out"] || (ok_in && !same(f.factor(), t["factor_out"]));
            for (int i = 0; i < 6; ++i) d = d || !same(u[static_cast<std::size_t>(i)], t["u" + std::to_string(i + 1) + "_out"]);
            if (d && bad++ < 5) std::printf("  trgTrackToPlane call %ld differs\n", n + 1);
            ++n;
        }
        std::printf("trgTrackToPlane (all captured calls): %ld compared, %ld differing\n", n, bad);
        fails += bad != 0 || n == 0;
    }
    {
        XzTable t;
        if (!t.open(cap + "/cap_from.txt.xz")) return 1;
        long n = 0, bad = 0;
        while (t.next()) {
            double x = t["x_in"], y = t["y_in"], z = t["z_in"], dx = t["dx_in"], dy = t["dy_in"];
            bool ok = false;
            f.trackFromTgt(x, y, z, dx, dy, t["mom_in"], t["mass_in"], static_cast<int>(t["spect_in"]), ok);
            const bool d = !same(x, t["x_out"]) || !same(y, t["y_out"]) || !same(z, t["z_out"]) || !same(dx, t["dx_out"]) ||
                           !same(dy, t["dy_out"]) || (ok ? 1.0 : 0.0) != t["ok_out"];
            if (d && bad++ < 5) std::printf("  track_from_tgt call %ld differs\n", n + 1);
            ++n;
        }
        std::printf("track_from_tgt (all captured calls): %ld compared, %ld differing\n", n, bad);
        fails += bad != 0 || n == 0;
    }
    {
        XzTable t;
        if (!t.open(cap + "/cap_to.txt.xz")) return 1;
        long n = 0, bad = 0, skipped = 0;
        while (t.next()) {
            if (t["ok_in"] == 0) {   // immediate return in SIMC (the spectrometer failed)
                ++skipped;
                continue;
            }
            Track tr;
            tr.xs = t["track_xs"];
            tr.ys = t["track_ys"];
            tr.dxdzs = t["track_dxdzs"];
            tr.dydzs = t["track_dydzs"];
            const int arm = static_cast<int>(t["arm_in"]);
            double delta = t["delta_in"], y = t["y_in"], dx = t["dx_in"], dy = t["dy_in"], mom = t["mom_in"];
            bool ok = true;
            int nr = 0;
            f.trackToTgt(delta, y, dx, dy, t["frx_in"], t["fry_in"], mom, t["mass_in"], t["ctheta_in"], t["stheta_in"],
                         static_cast<int>(t["spect_in"]), ok, arm,
                         [&](double fry, double& d, double& yp, double& xp, double& yt) {
                             eng.recon(arm, tr, fry, d, yp, xp, yt);
                         },
                         &nr);
            const bool d = !same(delta, t["delta_out"]) || !same(y, t["y_out"]) || !same(dx, t["dx_out"]) ||
                           !same(dy, t["dy_out"]) || !same(mom, t["mom_out"]) || (ok ? 1.0 : 0.0) != t["ok_out"] ||
                           nr != static_cast<int>(t["n_recon"]);
            if (d && bad++ < 5) std::printf("  track_to_tgt call %ld differs\n", n + skipped + 1);
            ++n;
        }
        std::printf("track_to_tgt (captured calls with ok_in = 1): %ld compared, %ld differing (%ld ok_in = 0 calls)\n", n,
                    bad, skipped);
        fails += bad != 0 || n == 0;
    }
    std::printf("%s\n", fails == 0 ? "PASS" : "FAIL");
    return fails == 0 ? 0 : 1;
}

bool slurp(const std::string& p, std::string& s) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return false;
    s.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::printf("usage: test_l1_trg <mode> <table_dir> <field_dir> [matrix_dir] | files <field_dir> <simc_dir>\n");
        return 1;
    }
    const std::string mode = argv[1];
    if (mode == "files") {
        std::string a, b;
        if (!slurp(std::string(argv[3]) + "/trg_field_map.dat", b)) {
            std::printf("SKIP: Fortran map missing\n");
            return 77;
        }
        if (!slurp(std::string(argv[2]) + "/trg_field_map.dat", a)) return 1;
        std::printf("data/field/trg_field_map.dat %s the Fortran file (%zu bytes)\n", a == b ? "equals" : "DIFFERS FROM",
                    a.size());
        return a == b ? 0 : 1;
    }
    const std::string dir = argv[2];
    g_fieldDir = argv[3];
    const std::string mdir = argc > 4 ? argv[4] : "";
    try {
        if (mode == "trginit") return runTrginit(dir);
        if (mode == "trginit_map") return runTrginitMap(dir);
        if (mode == "trgfield") return runTrgfield(dir);
        if (mode == "trgderiv") return runTrgderiv(dir);
        if (mode == "trgrk4") return runTrgrk4(dir);
        if (mode == "trgtracktoplane") return runTrackToPlane(dir);
        if (mode == "track_from_tgt") return runTrackFrom(dir);
        if (mode == "track_to_tgt") return runTrackTo(dir, mdir);
        if (mode == "replay") return runReplay(dir, mdir);
    } catch (const std::exception& e) {
        std::printf("error: %s\n", e.what());
        return 1;
    }
    std::printf("unknown mode %s\n", mode.c_str());
    return 1;
}
