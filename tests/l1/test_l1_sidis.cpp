// tests/l1/test_l1_sidis.cpp
//
// Tests of the SIDIS physics (task P12.1): semi_physics.f, semi_dilution.f, cteq5/Ctq5Pdf.f, fdss/fdss.f and the free
// nucleon part of F1F2IN21_v1.0.f.
//
//   files <data/sidis> <fortran dir>   md5 of the copies (MD5SUMS) + byte equality with simc_gfortran/cteq5/cteq5m.tbl
//                                      and simc_gfortran/fdss/KANLO.GRID
//   tables <data/sidis>                table reading: CTEQ5M header/grid/point count, fDSS grid records and the
//                                      interpolation arrays (spot values of the files, reconstruction of PARTON)
//   sanity <data/sidis>                physical sanity of every routine (valence sum rules, FF charge conjugation,
//                                      F2 ranges, EMC and attenuation limits, peepiX pi+/pi-/K+/D/threshold, the D-031
//                                      switches, semi_dilution in [0,1]) and the plug-in "sidis" (registry, refused
//                                      input combinations)
// The bitwise L1 comparisons against the P12.F1 Fortran tables (reference/functions/sidis/) are added when the tables
// exist. Exit: 0 pass, 1 fail, 77 data missing.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "L1Runner.h"
#include "Md5.h"
#include "simc/core/Errors.h"
#include "simc/core/Registry.h"
#include "simc/core/RunSetup.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SimcParameters.h"
#include "simc/physics/sidis/Cteq5Pdf.h"
#include "simc/physics/sidis/F1F2In21.h"
#include "simc/physics/sidis/Fdss.h"
#include "simc/physics/sidis/SidisCrossSection.h"
#include "simc/physics/sidis/SidisModel.h"
#include "simc/stages/IPhysicsModel.h"

using namespace simc;

namespace {

int g_fail = 0, g_checks = 0;
bool g_verbose = std::getenv("SIMC_TEST_VERBOSE") != nullptr;   // print every check
void check(bool ok, const std::string& what) {
    ++g_checks;
    if (g_verbose && ok) std::printf("  ok: %s\n", what.c_str());
    if (!ok) {
        ++g_fail;
        std::printf("  FAIL: %s\n", what.c_str());
    }
}
std::string fmt(const char* f, double a, double b = 0, double c = 0) {
    char buf[256];
    std::snprintf(buf, sizeof buf, f, a, b, c);
    return buf;
}

// ------------------------------------------------------------------------------------------------- files
int runFiles(const std::string& dataDir, const std::string& fortranDir) {
    std::ifstream sums(dataDir + "/MD5SUMS");
    if (!sums) return 77;
    const std::vector<std::pair<std::string, std::string>> where = {{"cteq5m.tbl", "cteq5/cteq5m.tbl"},
                                                                    {"KANLO.GRID", "fdss/KANLO.GRID"}};
    int bad = 0, n = 0;
    std::string md5, name;
    while (sums >> md5 >> name) {
        ++n;
        std::string a, b, rel;
        for (const auto& w : where)
            if (w.first == name) rel = w.second;
        if (rel.empty()) {
            std::printf("  %s: not a SIDIS data file\n", name.c_str());
            ++bad;
            continue;
        }
        if (!simcl1::slurp(dataDir + "/" + name, a)) {
            std::printf("  missing copy %s\n", name.c_str());
            ++bad;
            continue;
        }
        const std::string got = simcl1::md5hex(a);
        const bool okmd5 = got == md5;
        std::string cmp = "Fortran file missing (not compared)";
        bool okF = true;
        if (simcl1::slurp(fortranDir + "/" + rel, b)) {
            okF = a == b;
            cmp = okF ? "identical to the Fortran file" : "DIFFERS from the Fortran file";
        }
        std::printf("  %-12s md5 %s %s, %s %s\n", name.c_str(), got.c_str(), okmd5 ? "ok" : "MISMATCH", rel.c_str(),
                    cmp.c_str());
        if (!okmd5 || !okF) ++bad;
    }
    std::printf("sidis data files: %d listed, %d bad -> %s\n", n, bad, bad || n != 2 ? "FAIL" : "PASS");
    return (bad || n != 2) ? 1 : 0;
}

// ------------------------------------------------------------------------------------------------- tables
int runTables(const std::string& dataDir) {
    sidis::Cteq5Pdf p;
    std::string err;
    if (!p.load(dataDir, 1, &err)) {
        std::printf("cteq5m: %s\n", err.c_str());
        return 77;
    }
    check(p.Nx == 90 && p.Nt == 13 && p.NfMx == 5, "CTEQ5M NX NT NfMx = 90 13 5");
    check(p.Iorder == 2 && p.Nfl == 5 && p.Al == 0.2260 && p.Alambda == p.Al, "CTEQ5M Ordr Nfl lambda = 2 5 0.2260");
    check(p.Amass[3] == 1.3 && p.Amass[4] == 4.5 && p.Amass[5] == 180.0, "CTEQ5M quark masses");
    check(p.Qini == 1.0 && p.Qmax == 1.0e4 && p.Xmin == 1.0e-5, "CTEQ5M QINI QMAX XMIN");
    check(p.XV[0] == 0.0 && p.XV[1] == 9.99983E-06 && p.XV[90] == 1.0, "CTEQ5M XV(0) XV(1) XV(90)");
    check(p.QL[0] == std::log(1.0 / 0.2260) && p.QL[13] == std::log(1.0e4 / 0.2260), "CTEQ5M QL = ln(Q/Al)");
    check(p.points_expected() == 91 * 14 * 8 && p.points_read() == p.points_expected(),
          "CTEQ5M parton table: " + std::to_string(p.points_read()) + " of " + std::to_string(p.points_expected()) +
              " values read");
    check(p.UPD[0] == 0.0 && p.UPD[static_cast<std::size_t>(p.points_expected() - 1)] == 0.0,
          "CTEQ5M UPD(1) and UPD(Npts) (x = 0 / x = 1 nodes) are 0");
    // a few values in the middle: monotonic grids
    bool mono = true;
    for (int i = 1; i <= p.Nx; ++i) mono = mono && p.XV[static_cast<std::size_t>(i)] > p.XV[static_cast<std::size_t>(i - 1)];
    for (int i = 1; i <= p.Nt; ++i) mono = mono && p.QL[static_cast<std::size_t>(i)] > p.QL[static_cast<std::size_t>(i - 1)];
    check(mono, "CTEQ5M x and Q grids increasing");
    std::string e2;
    sidis::Cteq5Pdf bad;
    check(!bad.load(dataDir, 10, &e2) && e2.find("Invalid Iset") != std::string::npos, "SetCtq5(10) refused");
    check(!bad.load(dataDir + "/nonexistent", 1, &e2), "missing table refused");

    sidis::FdssGrid g;
    if (!g.load(dataDir, 2, 1, &err)) {
        std::printf("KANLO.GRID: %s\n", err.c_str());
        return 1;
    }
    auto P = [&](int i, int n, int m) {
        return g.parton[static_cast<std::size_t>(((m - 1) * sidis::FdssGrid::NQ + (n - 1)) * 9 + (i - 1))];
    };
    // first record (M = 1, N = 1) and last record (M = 34, N = 24) of KANLO.GRID
    check(P(1, 1, 1) == 2.372E-03 && P(2, 1, 1) == 1.600E-03 && P(3, 1, 1) == 8.859E-02 && P(9, 1, 1) == -8.699E-02,
          "KANLO.GRID record 1");
    check(P(1, 24, 34) == 4.894E-04 && P(8, 24, 34) == 0.0 && P(9, 24, 34) == -1.568E-03, "KANLO.GRID record 816");
    // XUTOTF(IX,IQ) = PARTON(1,IQ,IX)/(XB1**4*XB0**0.5); the x = 1 column is 0
    const double xb0 = sidis::FdssGrid::kXB[0], xb1 = 1 - xb0;
    check(g.F[0][0] == P(1, 1, 1) / (xb1 * xb1 * (xb1 * xb1) * std::pow(xb0, 0.5)), "XUTOTF(1,1)");
    check(g.F[3][34] == 0.0 && g.F[5][sidis::FdssGrid::NX * 24 - 1] == 0.0, "XCTOTF/XGF(NX, IQ) = 0");
    check(g.ARRF[0] == std::log(0.01) && g.ARRF[34] == 0.0 && g.ARRF[35] == 0.0 && g.ARRF[58] == std::log(1.0e5),
          "ARRF = log(XB), log(QS)");
    check(sidis::FdssGrid::kQS[5] == 6.4, "QS(6) = 6.4D0 rounded to REAL(8)");
    check(!g.load(dataDir, 5, 1, &err) && err.find("WRONG SET") != std::string::npos, "fDSS IH = 5 refused");
    std::printf("sidis tables: %d checks, %d failed -> %s\n", g_checks, g_fail, g_fail ? "FAIL" : "PASS");
    return g_fail ? 1 : 0;
}

// ------------------------------------------------------------------------------------------------- sanity helpers
// A SIDIS vertex: beam Ein, x, Q2 [MeV^2], z, pt2 [MeV^2], hadron mass mh. Directions: q along +z, e' in the x-z plane.
Event sidisVertex(double Ein, double x, double Q2, double z, double pt2, double mh) {
    using constants::Mp;
    Event v;
    v.Ein = Ein;
    v.Q2 = Q2;
    v.nu = Q2 / (2. * Mp * x);
    v.q = std::sqrt(v.nu * v.nu + Q2);
    v.e.E = Ein - v.nu;
    v.e.P = v.e.E;
    v.uq = {0, 0, 1};
    v.ue = {0.3, 0, std::sqrt(1 - 0.09)};
    v.zhad = z;
    v.pt2 = pt2;
    const double ph = std::sqrt((z * v.nu) * (z * v.nu) - mh * mh);
    v.theta_pq = std::asin(std::sqrt(pt2) / ph);
    v.p.E = z * v.nu;
    v.p.P = ph;
    v.xbj = x;
    return v;
}

sidis::PeepiXConfig cfgFor(bool pion, bool plus, double A, double Z) {
    sidis::PeepiXConfig c;
    c.targA = A;
    c.targZ = Z;
    c.Mtar_struck = constants::Mp;
    c.doing_semipi = pion;
    c.doing_semika = !pion;
    c.doing_hplus = plus;
    c.doing_decay = true;
    c.hadron_arm = 5;
    return c;
}

double integrate(const std::function<double(double)>& f, double lo, double hi, int n) {   // in ln x, Simpson
    const double a = std::log(lo), b = std::log(hi), h = (b - a) / n;
    double s = 0;
    for (int i = 0; i <= n; ++i) {
        const double x = std::exp(a + i * h);
        const double w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
        s += w * f(x) * x;
    }
    return s * h / 3;
}

// ------------------------------------------------------------------------------------------------- sanity
int runSanity(const std::string& dataDir) {
    sidis::Cteq5Pdf pdf;
    std::string err;
    if (!pdf.load(dataDir, 1, &err)) {
        std::printf("cteq5m: %s\n", err.c_str());
        return 77;
    }
    sidis::FdssGrid kgrid;
    if (!kgrid.load(dataDir, 2, 1, &err)) {
        std::printf("KANLO.GRID: %s\n", err.c_str());
        return 77;
    }
    std::ostringstream log;

    // polint: exact for a quadratic
    {
        const double xa[3] = {0.1, 0.4, 0.9}, ya[3] = {0.1 * 0.1 + 1, 0.4 * 0.4 + 1, 0.9 * 0.9 + 1};
        double dy;
        const double y = sidis::polint(xa, ya, 3, 0.55, dy);
        check(std::abs(y - (0.55 * 0.55 + 1)) < 1e-14, fmt("polint quadratic %.17g", y));
    }
    // CTEQ5M at Q = 2 GeV: valence sum rules int (u - ubar) = 2, int (d - dbar) = 1; u > d > ubar at x = 0.3
    {
        double Q = 2.0;
        auto val = [&](int ip) {
            return [&, ip](double x) {
                double q1 = Q, q2 = Q;
                return pdf.pdf(ip, x, q1) - pdf.pdf(-ip, x, q2);
            };
        };
        const double uv = integrate(val(1), 1e-5, 0.999, 4000), dv = integrate(val(2), 1e-5, 0.999, 4000);
        check(std::abs(uv - 2.0) < 0.03, fmt("CTEQ5M int u_v dx = %.4f (2)", uv));
        check(std::abs(dv - 1.0) < 0.03, fmt("CTEQ5M int d_v dx = %.4f (1)", dv));
        double q = 2.0;
        const double u = pdf.pdf(1, 0.3, q), d = pdf.pdf(2, 0.3, q), ub = pdf.pdf(-1, 0.3, q), s = pdf.pdf(3, 0.3, q);
        check(u > d && d > ub && ub > 0 && s > 0 && s == pdf.pdf(-3, 0.3, q), fmt("CTEQ5M x=0.3: u %.4g d %.4g ub %.4g", u, d, ub));
        // momentum sum rule incl. gluon (Iparton 0) and heavy flavours: ~1
        auto mom = [&](double x) {
            double sum = 0;
            for (int ip = -5; ip <= 5; ++ip) {
                double qq = 10.0;
                sum += pdf.pdf(ip, x, qq);
            }
            return x * sum;
        };
        const double msum = integrate(mom, 1e-5, 0.999, 4000);
        check(std::abs(msum - 1.0) < 0.03, fmt("CTEQ5M momentum sum at Q = 10 GeV = %.4f (1)", msum));
        // out of range: Iparton 6 -> 0 (warning once), Q < Alambda -> Q = Alambda, x > 1 -> STOP
        double qq = 2.0;
        check(pdf.pdf(6, 0.3, qq, &log) == 0.0 && !pdf.warn, "Ctq5Pdf Iparton 6 = 0 (Warn)");
        double ql = 0.1;
        pdf.pdf(1, 0.3, ql, &log);
        check(ql == pdf.Alambda, "Ctq5Pdf sets Q = Alambda");
        bool stopped = false;
        try {
            double q3 = 2.0;
            pdf.pdf(1, 1.5, q3);
        } catch (const FortranStop&) {
            stopped = true;
        }
        check(stopped, "Ctq5Pdf x > 1 STOP");
    }
    // ffint: table node values are reproduced, linear between nodes
    {
        const double ent[5] = {0, 1, 2, 0, 1};
        const int nent[2] = {3, 2};
        double table[6];
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 3; ++i) table[j * 3 + i] = 1 + 2 * i + 10 * j;   // f = 1 + 2x + 10y
        const double a1[2] = {1, 0}, a2[2] = {1.5, 0.25}, a3[2] = {2.5, 1.5};
        check(sidis::ffint(2, a1, nent, ent, table) == 3, "ffint node");
        check(std::abs(sidis::ffint(2, a2, nent, ent, table) - (1 + 3 + 2.5)) < 1e-14, "ffint bilinear");
        check(std::abs(sidis::ffint(2, a3, nent, ent, table) - (1 + 5 + 15)) < 1e-13, "ffint linear extrapolation");
    }
    // fDSS kaons: K+ favoured u > ubar, s-bar favoured; K- = charge conjugate; average for IC = 0
    {
        const auto kp = kgrid.eval(1, 0.5, 4.0), km = kgrid.eval(-1, 0.5, 4.0), k0 = kgrid.eval(0, 0.5, 4.0);
        check(kp.U > kp.UB && kp.SB > kp.S && kp.U > 0 && kp.SB > 0, fmt("fDSS K+ z=0.5 Q2=4: u %.4g ubar %.4g sbar %.4g", kp.U, kp.UB, kp.SB));
        check(km.U == kp.UB && km.UB == kp.U && km.S == kp.SB && km.D == kp.DB, "fDSS K- = conjugate of K+");
        check(k0.U == (kp.UB + kp.U) / 2. && k0.U == k0.UB, "fDSS IC = 0 average");
        const auto lowq = kgrid.eval(1, 0.5, 2.0);
        check(lowq.U > kp.U, "fDSS z D_u^K+ at z = 0.5 decreases with Q2");
        std::ostringstream w;
        kgrid.eval(1, 0.03, 0.5, &w);
        check(w.str().find("X OUT OF RANGE 0.0") != std::string::npos && w.str().find("Q2 OUT OF RANGE") != std::string::npos,
              "fDSS range warnings");
        bool stopped = false;
        try {
            kgrid.eval(2, 0.5, 4.0);
        } catch (const FortranStop&) {
            stopped = true;
        }
        check(stopped, "fDSS IC = 2 WRONG CHARGE STOP");
    }
    // F1F2IN21: proton/neutron F2 in the DIS region, Delta peak, Callan-Gross-like 2xF1 <= F2 (1+4M^2x^2/Q2)
    {
        double F1p, F2p, F1n, F2n;
        sidis::f1f2in21(1., 1., 4.0, 8.0, F1p, F2p);
        sidis::f1f2in21(0., 1., 4.0, 8.0, F1n, F2n);
        check(F2p > 0.2 && F2p < 0.45 && F2n > 0.1 && F2n < F2p, fmt("F1F2IN21 Q2=4 W2=8: F2p %.4f F2n %.4f", F2p, F2n));
        check(F1p > 0 && F1n > 0, "F1F2IN21 F1 > 0");
        double F1d, F2d, F1o, F2o;
        sidis::f1f2in21(1., 1., 0.5, 1.232 * 1.232, F1d, F2d);
        sidis::f1f2in21(1., 1., 0.5, 1.45 * 1.45, F1o, F2o);
        check(F2d > F2o, fmt("F1F2IN21 Delta peak F2(1.232) %.4f > F2(1.45) %.4f", F2d, F2o));
        double f1t, f2t;
        sidis::f1f2in21(1., 1., 1.0, 1.0, f1t, f2t);   // below the pion threshold
        check(f1t == 0 && f2t == 0, "F1F2IN21 below threshold = 0");
        bool np = false;
        try {
            sidis::f1f2in21(1., 2., 4., 8., f1t, f2t);
        } catch (const NotPorted&) {
            np = true;
        }
        check(np, "F1F2IN21 A = 2 NotPorted");
    }
    // slacemcfit, Rhad_global
    {
        check(sidis::slacemcfit(1., 0.3) == 1.0 && sidis::slacemcfit(2., 0.3) == 1.0, "slacemcfit(A <= 2) = 1");
        const double e12 = sidis::slacemcfit(12., 0.3), e12b = sidis::slacemcfit(12., 0.6);
        check(e12 > 0.95 && e12 < 1.02 && e12b < e12, fmt("slacemcfit C12 x=0.3 %.4f, x=0.6 %.4f (EMC slope)", e12, e12b));
        check(sidis::rhadGlobal(1., 0.5) == 1.0, "Rhad_global(1) = 1");
        const double r12 = sidis::rhadGlobal(12., 0.5), r63 = sidis::rhadGlobal(63., 0.5);
        check(r12 < 1 && r12 > 0.8 && r63 < r12, fmt("Rhad_global C %.4f Cu %.4f", r12, r63));
    }
    // peepiX: H(e,e'pi+-)X at E = 10.6 GeV, x = 0.35, Q2 = 4 GeV2, z = 0.5, pt2 = 0.1 GeV2
    {
        using constants::Mk;
        using constants::Mpi;
        const Event v = sidisVertex(10600., 0.35, 4.0e6, 0.5, 0.1e6, Mpi);
        sidis::PeepiXState st;
        MainEvent m;
        double sp = 0;
        sidis::PeepiXTrace tr;
        auto cp = cfgFor(true, true, 1, 1), cm = cfgFor(true, false, 1, 1);
        const double sigp = sidis::peepiX(v, v, m, sp, false, cp, pdf, nullptr, st, &log, &tr);
        check(sigp > 0 && std::isfinite(sigp), fmt("peepiX H pi+ = %.5g ub/MeV^2/sr^2", sigp));
        check(m.davejac == tr.jacobian * 1000.0 && m.ntup.sigcm == tr.sighad, "peepiX writes davejac, ntup%sigcm");
        check(std::abs(tr.xbj - 0.35) < 1e-12 && std::abs(tr.cthpq - std::cos(v.theta_pq)) == 0, "peepiX xbj, cthpq");
        check(std::abs(tr.wsq - (0.93827231 * 0.93827231 + 4.0 * (1 / 0.35 - 1))) < 1e-12, "peepiX W^2");
        const double sigm = sidis::peepiX(v, v, m, sp, false, cm, pdf, nullptr, st, &log, &tr);
        check(sigm > 0 && sigp / sigm > 1.5 && sigp / sigm < 5, fmt("peepiX H pi+/pi- = %.3f", sigp / sigm));
        // deuteron: pi+ and pi- closer
        auto dp = cfgFor(true, true, 2, 1), dm = cfgFor(true, false, 2, 1);
        const double sdp = sidis::peepiX(v, v, m, sp, false, dp, pdf, nullptr, st, &log, &tr);
        const double sdm = sidis::peepiX(v, v, m, sp, false, dm, pdf, nullptr, st, &log, &tr);
        check(sdp / sdm < sigp / sigm && sdp > sigp, fmt("peepiX D pi+/pi- = %.3f", sdp / sdm));
        // pt2 slope b = 1/(0.12 z^2+0.2): sigma(pt2)/sigma(pt2') = exp(-b (pt2-pt2')) at fixed jacobian
        const double b = 1. / (0.120 * 0.25 + 0.200);
        check(std::abs(tr.b - b) < 1e-15, "peepiX b(z)");
        // kaons
        const Event vk = sidisVertex(10600., 0.35, 4.0e6, 0.5, 0.1e6, Mk);
        auto kp = cfgFor(false, true, 1, 1), km = cfgFor(false, false, 1, 1);
        const double skp = sidis::peepiX(vk, vk, m, sp, false, kp, pdf, &kgrid, st, &log, &tr);
        const double skm = sidis::peepiX(vk, vk, m, sp, false, km, pdf, &kgrid, st, &log, &tr);
        check(skp > 0 && skm > 0 && skp > 2 * skm, fmt("peepiX H K+ %.4g, K- %.4g", skp, skm));
        check(skp < sigp, "peepiX K+ < pi+");
        // threshold: z -> 1 gives mmpi2 below (M + mh)^2: 0, davejac unchanged
        Event vt = sidisVertex(10600., 0.35, 4.0e6, 0.98, 0.01e6, Mpi);
        m.davejac = -7;
        const double s0 = sidis::peepiX(vt, vt, m, sp, false, cp, pdf, nullptr, st, &log, &tr);
        check(s0 == 0 && tr.threshold_return && m.davejac == -7, "peepiX threshold return");
        // unknown hadron: STOP
        auto none = cp;
        none.doing_semipi = false;
        bool stopped = false;
        try {
            sidis::peepiX(v, v, m, sp, false, none, pdf, nullptr, st, &log, &tr);
        } catch (const FortranStop&) {
            stopped = true;
        }
        check(stopped, "peepiX unknown hadron STOP");
        // pi0 (P12.1 F1): compat uses db = u1 (pi+ value), corrected the average; d1 is the average in both
        auto c0 = cp;
        c0.doing_pizero = true;
        sidis::PeepiXTrace t0, t1;
        sidis::peepiX(v, v, m, sp, false, c0, pdf, nullptr, st, &log, &t0);
        c0.pizero_db_fixed = true;
        sidis::peepiX(v, v, m, sp, false, c0, pdf, nullptr, st, &log, &t1);
        check(t0.db == t0.yf && t1.db == (t1.yf + t1.yu) / 2. && t0.d1 == t1.d1 && t0.u1 == t1.u1 && t0.dsigdz != t1.dsigdz,
              "pi0 dbar FF: compat yf, corrected average (P12.1 F1)");
        // Fermi motion: xfermi written, fac
        auto cf = dp;
        cf.do_fermi = true;
        MainEvent mf;
        mf.pfermi = {100., 0., 0., 1., std::sqrt(constants::Mp * constants::Mp + 1e4)};
        sidis::peepiX(v, v, mf, sp, false, cf, pdf, nullptr, st, &log, &tr);
        check(mf.ntup.xfermi == tr.xbj && tr.xbj != 0.35 && tr.fac != 1.0, "peepiX do_fermi: xfermi, flux factor");
        // survival probability (doing_decay = 0), HMS (arm 1, zaero = -331.491) and SHMS (arm 5)
        auto cs = kp;
        cs.doing_decay = false;
        cs.ctau = 371.2;
        cs.spec_p_P = 3000.;
        cs.Mh2 = Mk * Mk;
        cs.hadron_arm = 1;
        MainEvent ms;
        ms.FP.p.path = 2000.;
        ms.FP.p.dx = 0.05;
        ms.FP.p.dy = 0.02;
        ms.SP.p.delta = 5.;
        sidis::PeepiXState ss;
        double spc = 0;
        sidis::peepiX(vk, vk, ms, spc, false, cs, pdf, &kgrid, ss, &log, &tr);
        const double bk = 3000. / std::sqrt(9e6 + Mk * Mk), gk = 1 / std::sqrt(1 - bk * bk);
        check(spc > 0 && spc < 1 && ms.ntup.decdist == spc && ss.zaero == -331.491 &&
                  spc == 1. / std::exp((2000. + -331.491 * (1 + 0.0025 + 0.0004)) / (371.2 * bk * gk)),
              fmt("survival HMS compat %.6f", spc));
        cs.survival_p_hadron = cs.survival_path_fixed = true;
        double spr = 0;
        sidis::peepiX(vk, vk, ms, spr, false, cs, pdf, &kgrid, ss, &log, &tr);
        check(spr > spc && spr < 1, fmt("survival HMS corrected %.6f (higher momentum, shorter path)", spr));
        cs.hadron_arm = 5;
        cs.survival_p_hadron = cs.survival_path_fixed = false;
        ss.zaero = -183.;   // stale value of an earlier arm-3 call
        sidis::peepiX(vk, vk, ms, spc, false, cs, pdf, &kgrid, ss, &log, &tr);
        check(ss.zaero == -183., "SHMS compat keeps the SAVEd zaero (P12.1 F4)");
        cs.survival_path_fixed = true;
        sidis::peepiX(vk, vk, ms, spc, false, cs, pdf, &kgrid, ss, &log, &tr);
        check(ss.zaero == 0. && tr.pathlen == 2000., "SHMS corrected zaero = 0 (P12.1 F4)");
        // doing_cent (dead code): z binning
        auto cc = cp;
        cc.sigc_flag = 0;
        cc.sigc_nbin = 10;
        cc.sigc_kin_min = 0.;
        cc.sigc_kin_max = 1.;
        cc.sigc_kin_ind = 0.1;
        MainEvent mc;
        sidis::PeepiXState sc;
        const double scent = sidis::peepiX(v, v, mc, sp, true, cc, pdf, nullptr, sc, &log, &tr);
        check(scent > 0 && mc.johnjac == tr.jacobian * 1000.0 && mc.davejac == 0 && std::abs(tr.zhad - 0.45) < 1e-15 &&
                  tr.pt2 == 0.1e6 && !sc.first_cent,
              fmt("peepiX doing_cent z bin centre %.3f", tr.zhad));
    }
    // semi_dilution (dead code): in (0, 1) for hydrogen pi+
    {
        const Event v = sidisVertex(10600., 0.35, 4.0e6, 0.5, 0.1e6, constants::Mpi);
        sidis::SemiDilutionState sd;
        const double dil = sidis::semiDilution(v, true, false, true, false, pdf, sd, &log);
        check(dil > 0 && dil < 1, fmt("semi_dilution H pi+ = %.4f", dil));
        const double dilk = sidis::semiDilution(v, false, true, true, false, pdf, sd, &log);
        check(dilk == dil && !sd.kfirst, "semi_dilution kaons use the pion FFs");
    }
    // the plug-in "sidis"
    {
        auto model = Registry<IPhysicsModel>::create("sidis");
        check(model && model->name() == "sidis", "Registry creates \"sidis\"");
        SimcParameters p;
        RunSetup s = RunSetup::FromParameters(p);
        s.flags.doing_semi = true;
        s.flags.doing_semipi = true;
        s.flags.doing_hplus = true;
        s.targ.A = 1;
        s.targ.Z = 1;
        s.targ.Mtar_struck = constants::Mp;
        s.params.doing_decay = true;
        s.params.hadron_arm = 5;
        model->configure(Json{{"data_dir", dataDir}, {"verbose", false}}, p);
        model->init(StageContext(s));
        const Event v = sidisVertex(10600., 0.35, 4.0e6, 0.5, 0.1e6, constants::Mpi);
        MainEvent m;
        sidis::PeepiXConfig c = sidis::makeConfig(s);
        sidis::PeepiXState st;
        double sp = 0;
        const double ref = sidis::peepiX(v, v, m, sp, false, c, pdf, nullptr, st);
        MainEvent m2;
        check(model->sigcc(v, m2) == ref && m2.davejac == 0, "plug-in sigcc == peepiX, main unchanged");
        check(model->sigccMain(v, m2) == ref && m2.davejac == m.davejac, "plug-in sigccMain writes main");
        check(model->sigccRecon(v, m2) == 1.0, "plug-in sigcc_recon = 1");
        check(c.pizero_db_fixed && c.survival_p_hadron && c.survival_path_fixed, "default mode = corrected");
        s.fortran_compat = true;
        const auto cc = sidis::makeConfig(s);
        check(!cc.pizero_db_fixed && !cc.survival_p_hadron && !cc.survival_path_fixed, "fortran_compat = Fortran");
        // refused combinations
        RunSetup s0 = s;
        s0.flags.doing_pizero = true;
        bool np = false;
        try {
            auto m0 = Registry<IPhysicsModel>::create("sidis");
            m0->configure(Json{{"data_dir", dataDir}, {"verbose", false}}, p);
            m0->init(StageContext(s0));
        } catch (const NotPorted&) {
            np = true;
        }
        check(np, "pi0 SIDIS refused (NotPorted)");
        RunSetup s1 = s;
        s1.flags.doing_semipi = false;
        bool ia = false;
        try {
            auto m1 = Registry<IPhysicsModel>::create("sidis");
            m1->configure(Json{{"data_dir", dataDir}, {"verbose", false}}, p);
            m1->init(StageContext(s1));
        } catch (const std::invalid_argument&) {
            ia = true;
        }
        check(ia, "doing_semi without pion/kaon refused");
        RunSetup s2 = s;
        s2.params.pt_b_param = 4.68;
        check(sidis::ignoredInputNotes(s2).size() == 1, "pt_b_param ignored-input note");
    }
    std::printf("sidis sanity: %d checks, %d failed -> %s\n", g_checks, g_fail, g_fail ? "FAIL" : "PASS");
    return g_fail ? 1 : 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <mode> ... (see the file header)\n", argv[0]);
        return 2;
    }
    const std::string what = argv[1];
    try {
        if (what == "files" && argc > 3) return runFiles(argv[2], argv[3]);
        if (what == "tables" && argc > 2) return runTables(argv[2]);
        if (what == "sanity" && argc > 2) return runSanity(argv[2]);
    } catch (const std::exception& e) {
        std::printf("ERROR: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "bad arguments for mode %s\n", what.c_str());
    return 2;
}
