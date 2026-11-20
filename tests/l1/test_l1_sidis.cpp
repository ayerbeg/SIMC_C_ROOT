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
// Bitwise L1 comparisons against the P12.F1 Fortran tables (reference/functions/sidis/<routine>/, fortran_compat; the
// tolerance argument is |rel|, 0 = bitwise; NaN/Inf must match):
//   ctq5pdf <dir> <data/sidis> <tol>           Ctq5Pdf(ipart, x, Q) in table order (one process): pdf, Q_after
//   ctq5_state <dir> <data/sidis>              COMMON /CtqPar1/ ... after SetCtq5(1) (ctq5_state.txt), bitwise
//   fdss <dir> <data/sidis> <fortran fdss dir> <tol>   fDSS(IH, IC, IO, X, Q2) for every set (KANLO.GRID from
//                                              data/sidis, the unreachable grids from the Fortran fdss/ directory)
//   f1f2in21 <dir> <tol>                       F1F2IN21(Z, A, Q2, W2): F1, F2
//   slacemcfit <dir> <tol>                     slacemcfit(A, x), Rhad_global(A, x)
//   fortran_diff <peepix dir> <sidis_state dir> <inp dirs> <data/sidis> <data/spectral>
//                                              D-051 F1/F4/F6 in both modes on the peepix rows (F2: sanity)
//   peepix <dir> <sidis_state dir> <inp dirs (':'-separated)> <data/sidis> <data/spectral> <tol>
//                                              peepiX with the per-configuration set-up (BuildRunSetup of the cfg's .inp,
//                                              fortran_compat), the row overrides (ov_state) and the SAVEd state; the
//                                              outputs and the harness intermediates hx_* (rows in table order per cfg)
// Exit: 0 pass, 1 fail, 77 data missing.
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
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
#include "simc/input/InpReader.h"
#include "simc/setup/SetupBuilder.h"
#include "simc/physics/sidis/Cteq5Pdf.h"
#include "simc/physics/sidis/F1F2In21.h"
#include "simc/physics/sidis/Fdss.h"
#include "simc/physics/sidis/SidisCrossSection.h"
#include "simc/physics/sidis/SidisModel.h"
#include "simc/stages/IPhysicsModel.h"

using namespace simc;
using simcl1::ColumnSpec;
using simcl1::Outputs;
using simcl1::Row;

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
    check(p.Iorder == 2 && p.Nfl == 5 && p.Al() == 0.2260 && p.Alambda == p.Al(), "CTEQ5M Ordr Nfl lambda = 2 5 0.2260");
    check(p.Amass[3] == 1.3 && p.Amass[4] == 4.5 && p.Amass[5] == 180.0, "CTEQ5M quark masses");
    check(p.Qini == 1.0 && p.Qmax == 1.0e4 && p.Xmin == 1.0e-5, "CTEQ5M QINI QMAX XMIN");
    check(p.XV(0) == 0.0 && p.XV(1) == 9.99983E-06 && p.XV(90) == 1.0, "CTEQ5M XV(0) XV(1) XV(90)");
    check(p.QL(0) == std::log(1.0 / 0.2260) && p.QL(13) == std::log(1.0e4 / 0.2260), "CTEQ5M QL = ln(Q/Al)");
    check(p.points_expected() == 91 * 14 * 8 && p.points_read() == p.points_expected(),
          "CTEQ5M parton table: " + std::to_string(p.points_read()) + " of " + std::to_string(p.points_expected()) +
              " values read");
    check(p.UPD(1) == 0.0 && p.UPD(p.points_expected()) == 0.0,
          "CTEQ5M UPD(1) and UPD(Npts) (x = 0 / x = 1 nodes) are 0");
    // a few values in the middle: monotonic grids
    bool mono = true;
    for (int i = 1; i <= p.Nx; ++i) mono = mono && p.XV(i) > p.XV(i - 1);
    for (int i = 1; i <= p.Nt; ++i) mono = mono && p.QL(i) > p.QL(i - 1);
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
        double qn = 2.0;
        check(std::isnan(pdf.pdf(1, std::nan(""), qn)), "Ctq5Pdf NaN x -> NaN (as the Fortran)");
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
        // pi0 (D-051 F4): compat uses db = u1 (pi+ value) and mhad = Mpi, corrected the average and Mpi0
        auto c0 = cp;
        c0.doing_pizero = true;
        sidis::PeepiXTrace t0, t1;
        sidis::peepiX(v, v, m, sp, false, c0, pdf, nullptr, st, &log, &t0);
        c0.pizero_fixed = true;
        sidis::peepiX(v, v, m, sp, false, c0, pdf, nullptr, st, &log, &t1);
        check(t0.db == t0.yf && t1.db == (t1.yf + t1.yu) / 2. && t0.d1 == t1.d1 && t0.u1 == t1.u1 &&
                  t0.dsigdz != t1.dsigdz && t0.mhad == Mpi && t1.mhad == constants::Mpi0,
              "pi0 dbar FF and mass: compat yf/Mpi, corrected average/Mpi0 (D-051 F4)");
        // D-051 F1: x_Bj > 1 -> clamped to 1 -> every PDF 0 -> compat 0/0 = NaN, corrected 0
        {
            const Event vx = sidisVertex(10600., 1.2, 4.0e6, 0.3, 0.05e6, Mpi);
            auto c1 = cp;
            MainEvent m1;
            sidis::PeepiXTrace tx;
            const double nanv = sidis::peepiX(vx, vx, m1, sp, false, c1, pdf, nullptr, st, &log, &tx);
            c1.sum_sq_guard = true;
            const double zero = sidis::peepiX(vx, vx, m1, sp, false, c1, pdf, nullptr, st, &log, &tx);
            check(tx.xbj == 1.0 && tx.sum_sq == 0 && std::isnan(nanv) && zero == 0.0,
                  fmt("sum_sq = 0: compat %g, corrected %g (D-051 F1)", nanv, zero));
        }
        // D-051 F2: Fermi momentum along q with P.q < 0 -> negative x_Bj: compat Ctq5Pdf STOP, corrected 0
        {
            auto c2 = dp;
            c2.do_fermi = true;
            MainEvent m2;
            m2.pfermi = {6000., 0., 0., 1., std::sqrt(constants::Mp * constants::Mp + 3.6e7)};   // pfer along q
            bool stop2 = false;
            try {
                sidis::peepiX(v, v, m2, sp, false, c2, pdf, nullptr, st, &log, &tr);
            } catch (const FortranStop&) {
                stop2 = true;
            }
            c2.nonpositive_xbj_zero = true;
            m2.davejac = -5;
            const double z2 = sidis::peepiX(v, v, m2, sp, false, c2, pdf, nullptr, st, &log, &tr);
            check(stop2 && z2 == 0.0 && tr.xbj < 0 && m2.ntup.xfermi == tr.xbj && m2.davejac == -5,
                  fmt("negative Fermi x_Bj %g: compat STOP, corrected 0 (D-051 F2)", tr.xbj));
        }
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
        check(ss.zaero == -183., "SHMS compat keeps the SAVEd zaero (D-051 F6)");
        cs.survival_path_fixed = true;
        sidis::peepiX(vk, vk, ms, spc, false, cs, pdf, &kgrid, ss, &log, &tr);
        check(ss.zaero == 0. && tr.pathlen == 2000., "SHMS corrected zaero = 0 (D-051 F6)");
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
        check(c.sum_sq_guard && c.nonpositive_xbj_zero && c.pizero_fixed && c.survival_p_hadron && c.survival_path_fixed,
              "default mode = corrected");
        s.fortran_compat = true;
        const auto cc = sidis::makeConfig(s);
        check(!cc.sum_sq_guard && !cc.nonpositive_xbj_zero && !cc.pizero_fixed && !cc.survival_p_hadron &&
                  !cc.survival_path_fixed,
              "fortran_compat = Fortran");
        // refused combinations
        RunSetup s0 = s;
        s0.flags.doing_pizero = true;
        bool np = true;
        try {
            auto m0 = Registry<IPhysicsModel>::create("sidis");
            m0->configure(Json{{"data_dir", dataDir}, {"verbose", false}}, p);
            m0->init(StageContext(s0));
        } catch (const std::exception&) {
            np = false;
        }
        check(np, "pi0 SIDIS accepted (calorimeter arm, P14.3)");
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


// ================================================================================================= L1 tables
int iv(const Row& r, const char* c) { return static_cast<int>(std::lround(r[c])); }
std::vector<ColumnSpec> specs(const std::vector<std::string>& names, double tol) {
    std::vector<ColumnSpec> v;
    for (const auto& n : names) v.push_back(ColumnSpec{n, tol, 0.0});
    return v;
}
bool sameBits(double a, double b) { return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof a) == 0; }

int runCtq5pdf(const std::string& dir, const std::string& dataDir, double tol) {
    sidis::Cteq5Pdf pdf;
    std::string err;
    if (!pdf.load(dataDir, 1, &err)) {
        std::printf("cteq5m: %s\n", err.c_str());
        return 77;
    }
    std::ostringstream log;
    return simcl1::run("ctq5pdf (Ctq5Pdf.f:Ctq5Pdf, SetCtq5(1))", dir, specs({"pdf", "Q_after"}, tol),
                       [&](const Row& r, Outputs& o) {
                           double q = r["Q"];
                           o["pdf"] = pdf.pdf(iv(r, "ipart"), r["x"], q, &log);
                           o["Q_after"] = q;
                       });
}

int runCtq5State(const std::string& dir, const std::string& dataDir) {
    simcl1::Table t;
    std::string err;
    if (!t.load(dir + "/ctq5_state.txt", err)) {
        std::printf("%s\n", err.c_str());
        return 77;
    }
    sidis::Cteq5Pdf p;
    if (!p.load(dataDir, 1, &err)) {
        std::printf("cteq5m: %s\n", err.c_str());
        return 1;
    }
    Row r(t, 0);
    int n = 0, bad = 0;
    auto c = [&](const std::string& name, double v) {
        if (!r.has(name)) return;
        ++n;
        if (!sameBits(v, r[name])) {
            ++bad;
            std::printf("  %s: ref %.17g cpp %.17g\n", name.c_str(), r[name], v);
        }
    };
    c("Al", p.Al());
    c("Alambda", p.Alambda);
    c("Nfl", p.Nfl);
    c("Iorder", p.Iorder);
    c("Nx", p.Nx);
    c("Nt", p.Nt);
    c("NfMx", p.NfMx);
    c("Qini", p.Qini);
    c("Qmax", p.Qmax);
    c("Xmin", p.Xmin);
    for (int i = 1; i <= 6; ++i) c("Amass_" + std::to_string(i), p.Amass[static_cast<std::size_t>(i - 1)]);
    for (int i = 0; i <= sidis::Cteq5Pdf::MXQ; ++i) c("QL_" + std::to_string(i), p.QL(i));
    for (int i = 0; i <= sidis::Cteq5Pdf::MXX; ++i) c("XV_" + std::to_string(i), p.XV(i));
    std::printf("ctq5_state: %d of %zu columns compared, %d differ -> %s\n", n, t.names.size(), bad,
                bad || n != static_cast<int>(t.names.size()) ? "FAIL" : "PASS");
    return bad || n != static_cast<int>(t.names.size()) ? 1 : 0;
}

int runFdss(const std::string& dir, const std::string& dataDir, const std::string& fortranFdss, double tol) {
    std::map<std::pair<int, int>, sidis::FdssGrid> grids;
    std::string err;
    auto grid = [&](int ih, int io) -> const sidis::FdssGrid& {
        auto k = std::make_pair(ih, io);
        auto it = grids.find(k);
        if (it == grids.end()) {
            sidis::FdssGrid g;
            const bool reach = ih == 2 && io == 1;
            if (!g.load(reach ? dataDir : fortranFdss, ih, io, &err)) throw std::runtime_error(err);
            it = grids.emplace(k, std::move(g)).first;
        }
        return it->second;
    };
    return simcl1::run("fdss (fdss.f:fDSS)", dir,
                       specs({"U", "UB", "D", "DB", "S", "SB", "C", "B", "GL", "n_warn_x", "n_warn_q2"}, tol),
                       [&](const Row& r, Outputs& o) {
                           std::ostringstream log;
                           const auto v = grid(iv(r, "IH"), iv(r, "IO")).eval(iv(r, "IC"), r["X"], r["Q2"], &log);
                           o["U"] = v.U;
                           o["UB"] = v.UB;
                           o["D"] = v.D;
                           o["DB"] = v.DB;
                           o["S"] = v.S;
                           o["SB"] = v.SB;
                           o["C"] = v.C;
                           o["B"] = v.B;
                           o["GL"] = v.GL;
                           const std::string l = log.str();
                           o["n_warn_x"] = l.find("X OUT OF RANGE") != std::string::npos ? 1 : 0;
                           o["n_warn_q2"] = l.find("Q2 OUT OF RANGE") != std::string::npos ? 1 : 0;
                       });
}

int runF1f2in21(const std::string& dir, double tol) {
    return simcl1::run("f1f2in21 (F1F2IN21_v1.0.f, A = 1)", dir, specs({"F1", "F2"}, tol), [&](const Row& r, Outputs& o) {
        double F1 = 0, F2 = 0;
        sidis::f1f2in21(r["Z"], r["A"], r["Q2"], r["W2"], F1, F2);
        o["F1"] = F1;
        o["F2"] = F2;
    });
}

int runSlac(const std::string& dir, double tol) {
    return simcl1::run("slacemcfit / Rhad_global (semi_physics.f)", dir, specs({"slacemcfit", "rhad_global"}, tol),
                       [&](const Row& r, Outputs& o) {
                           o["slacemcfit"] = sidis::slacemcfit(r["A"], r["x"]);
                           o["rhad_global"] = sidis::rhadGlobal(r["A"], r["x"]);
                       });
}

// ------------------------------------------------------------------------------------------------- peepix
Event readEvent(const Row& r, const std::string& p) {
    Event e;
    auto g = [&](const char* n) { return r[p + n]; };
    e.Ein = g("Ein"); e.Em = g("Em"); e.Pm = g("Pm"); e.Emiss = g("Emiss"); e.Pmiss = g("Pmiss");
    e.Pmx = g("Pmx"); e.Pmy = g("Pmy"); e.Pmz = g("Pmz"); e.PmPar = g("PmPar"); e.PmPer = g("PmPer"); e.PmOop = g("PmOop");
    e.nu = g("nu"); e.q = g("q"); e.Q2 = g("Q2"); e.Trec = g("Trec"); e.W = g("W"); e.Mrec = g("Mrec");
    e.epsilon = g("epsilon"); e.theta_pq = g("theta_pq"); e.theta_tarq = g("theta_tarq"); e.phi_pq = g("phi_pq");
    e.phi_targ = g("phi_targ"); e.beta = g("beta"); e.phi_s = g("phi_s"); e.phi_c = g("phi_c");
    e.zhad = g("zhad"); e.pt2 = g("pt2"); e.xbj = g("xbj");
    auto arm = [&](ArmFull& a, const std::string& s) {
        a.delta = r[p + s + "delta"]; a.xptar = r[p + s + "xptar"]; a.yptar = r[p + s + "yptar"]; a.z = r[p + s + "z"];
        a.theta = r[p + s + "theta"]; a.phi = r[p + s + "phi"]; a.E = r[p + s + "E"]; a.P = r[p + s + "P"];
    };
    arm(e.e, "e_");
    arm(e.p, "p_");
    e.ue = {g("ue_x"), g("ue_y"), g("ue_z")};
    e.up = {g("up_x"), g("up_y"), g("up_z")};
    e.uq = {g("uq_x"), g("uq_y"), g("uq_z")};
    return e;
}

std::string findInp(const std::string& dirs, const std::string& name) {
    std::string base = name;
    if (base.find('/') != std::string::npos) base = base.substr(base.find_last_of('/') + 1);
    if (base.size() < 4 || base.substr(base.size() - 4) != ".inp") base += ".inp";
    std::stringstream ss(dirs);
    for (std::string d; std::getline(ss, d, ':');) {
        const std::string p = d + "/" + base;
        if (std::ifstream(p).good()) return p;
    }
    return "";
}

std::map<int, std::string> readConfigs(const std::string& stateDir) {
    std::map<int, std::string> m;
    std::ifstream in(stateDir + "/configs.txt");
    std::string line;
    std::getline(in, line);   // header
    while (std::getline(in, line)) {
        std::istringstream ls(line);
        int cfg;
        std::string inp;
        if (ls >> cfg >> inp) m[cfg] = inp;
    }
    return m;
}

int runPeepix(const std::string& dir, const std::string& stateDir, const std::string& inpDirs, const std::string& dataDir,
              const std::string& spectralDir, double tol) {
    const auto cfgs = readConfigs(stateDir);
    if (cfgs.empty()) {
        std::printf("no %s/configs.txt\n", stateDir.c_str());
        return 77;
    }
    struct Proc {
        RunSetup setup;
        sidis::PeepiXConfig cfg;
        sidis::PeepiXState st;
        sidis::Cteq5Pdf pdf;
        sidis::FdssGrid fdss;
    };
    std::map<int, std::unique_ptr<Proc>> procs;
    sidis::Cteq5Pdf pdf0;
    sidis::FdssGrid fdss0;
    std::string err;
    if (!pdf0.load(dataDir, 1, &err) || !fdss0.load(dataDir, 2, 1, &err)) {
        std::printf("%s\n", err.c_str());
        return 77;
    }
    auto proc = [&](int cfg) -> Proc& {
        auto it = procs.find(cfg);
        if (it != procs.end()) return *it->second;
        const auto c = cfgs.find(cfg);
        if (c == cfgs.end()) throw std::runtime_error("cfg " + std::to_string(cfg) + " not in configs.txt");
        const std::string inp = findInp(inpDirs, c->second);
        if (inp.empty()) throw std::runtime_error("no .inp for cfg " + std::to_string(cfg) + " (" + c->second + ")");
        auto P = std::make_unique<Proc>();
        SimcParameters p;
        InpReader rd(p);
        if (!rd.loadRun(inp)) throw std::runtime_error("cannot read " + inp);
        SetupOptions so;
        so.fortran_compat = true;
        so.compat_target = CompatTarget::U294f662;   // U.R: 294f662 references (frozen tree)
        so.data_dir = spectralDir;
        P->setup = BuildRunSetup(p, so);
        if (!P->setup.ok()) throw std::runtime_error(inp + ": setup STOP " + P->setup.stops.front());
        P->cfg = sidis::makeConfig(P->setup);
        P->pdf = pdf0;   // one Fortran process per cfg: fresh SAVEd Warn/First
        P->fdss = fdss0;
        return *procs.emplace(cfg, std::move(P)).first->second;
    };
    const std::vector<std::string> outs = {"peepix", "survivalprob", "m_davejac", "m_johnjac", "ntup_sigcm",
                                           "ntup_xfermi", "decdist_out"};
    const std::vector<std::string> hx = {
        "hx_thr_return", "hx_targA", "hx_targZ", "hx_targN", "hx_Mpi_gev", "hx_mhad", "hx_mtar", "hx_Ehad", "hx_phad",
        "hx_cthpq", "hx_xbj_raw", "hx_xbj", "hx_Q2gev", "hx_Qgev", "hx_pt2gev", "hx_wsq", "hx_w", "hx_mmpi2", "hx_u",
        "hx_ubar", "hx_d", "hx_dbar", "hx_s", "hx_sbar", "hx_uA", "hx_ubarA", "hx_dA", "hx_dbarA", "hx_sA", "hx_sbarA",
        "hx_sum_sq", "hx_xp", "hx_zp", "hx_sv", "hx_yf", "hx_yu", "hx_u1", "hx_ub", "hx_d1", "hx_db", "hx_s1", "hx_sb",
        "hx_dsigdz", "hx_b", "hx_Rhad", "hx_sighad", "hx_F1p", "hx_F2p", "hx_F1n", "hx_F2n", "hx_slac", "hx_F1",
        "hx_F2", "hx_W1", "hx_W2", "hx_sin2th2", "hx_sige", "hx_sigsemi", "hx_jacobian", "hx_fac", "hx_peepix",
        "hx_Qgev_after"};
    // hx_* are the harness's re-evaluation (hxsemi.f) of the peepiX intermediates: all 0 on the doing_cent rows (not
    // evaluated). hx_pathlen/hx_survival are not compared: the harness uses zaero = 0 for the arms without a zaero
    // branch (5, 6) where peepiX keeps the SAVEd zaero; the SIMC values are the outputs survivalprob/decdist_out.
    std::vector<std::string> all = outs;
    // the hx_* columns present in the table (an absent one is skipped)
    simcl1::Table probe;
    if (probe.load(dir + "/output.txt", err))
        for (const auto& h : hx)
            if (probe.col(h) >= 0) all.push_back(h);
    std::ostringstream log;
    return simcl1::run("peepix (semi_physics.f:peepiX)", dir, specs(all, tol), [&](const Row& r, Outputs& o) {
        Proc& P = proc(iv(r, "cfg"));
        sidis::PeepiXConfig c = P.cfg;
        if (iv(r, "ov_state") == 1) {
            c.doing_hplus = iv(r, "ov_hplus") != 0;
            c.doing_pizero = iv(r, "ov_pizero") != 0;
            c.do_fermi = iv(r, "ov_do_fermi") != 0;
            c.sigc_flag = iv(r, "ov_sigc_flag");
            c.sigc_nbin = iv(r, "ov_sigc_nbin");
            c.sigc_kin_min = r["ov_sigc_kin_min"];
            c.sigc_kin_max = r["ov_sigc_kin_max"];
            c.sigc_kin_ind = r["ov_sigc_kin_ind"];
            c.hadron_arm = iv(r, "ov_hadron_arm");
        }
        const Event v = readEvent(r, "v_");
        const Event v0 = readEvent(r, "v0_");
        MainEvent m;
        m.theta_pq = r["in_m_theta_pq"];
        m.FP.p.path = r["in_m_FP_p_path"];
        m.FP.p.dx = r["in_m_FP_p_dx"];
        m.FP.p.dy = r["in_m_FP_p_dy"];
        m.SP.p.delta = r["in_m_SP_p_delta"];
        m.davejac = r["in_m_davejac"];
        m.johnjac = r["in_m_johnjac"];
        m.ntup.sigcm = r["in_ntup_sigcm"];
        m.ntup.xfermi = r["in_ntup_xfermi"];
        m.ntup.decdist = r["decdist_in"];
        m.pfermi = {r["pfer"], r["pferx"], r["pfery"], r["pferz"], r["efer"]};
        double sp = r["in_survivalprob"];
        sidis::PeepiXTrace t;
        const double res = sidis::peepiX(v, v0, m, sp, iv(r, "doing_cent") != 0, c, P.pdf,
                                         c.doing_semika ? &P.fdss : nullptr, P.st, &log, &t);
        o["peepix"] = res;
        o["survivalprob"] = sp;
        o["m_davejac"] = m.davejac;
        o["m_johnjac"] = m.johnjac;
        o["ntup_sigcm"] = m.ntup.sigcm;
        o["ntup_xfermi"] = m.ntup.xfermi;
        o["decdist_out"] = m.ntup.decdist;
        const std::map<std::string, double> h = {
            {"hx_thr_return", t.threshold_return ? 1 : 0}, {"hx_targA", c.targA}, {"hx_targZ", c.targZ},
            {"hx_targN", c.targA - c.targZ}, {"hx_Mpi_gev", t.Mpi_gev}, {"hx_mhad", t.mhad}, {"hx_mtar", t.mtar},
            {"hx_Ehad", t.Ehad}, {"hx_phad", t.phad}, {"hx_cthpq", t.cthpq}, {"hx_xbj_raw", t.xbj_raw},
            {"hx_xbj", t.xbj}, {"hx_Q2gev", t.Q2gev}, {"hx_Qgev", t.Qgev}, {"hx_pt2gev", t.pt2gev}, {"hx_wsq", t.wsq},
            {"hx_w", t.w}, {"hx_mmpi2", t.mmpi2}, {"hx_u", t.u}, {"hx_ubar", t.ubar}, {"hx_d", t.d},
            {"hx_dbar", t.dbar}, {"hx_s", t.s}, {"hx_sbar", t.sbar}, {"hx_uA", t.uA}, {"hx_ubarA", t.ubarA},
            {"hx_dA", t.dA}, {"hx_dbarA", t.dbarA}, {"hx_sA", t.sA}, {"hx_sbarA", t.sbarA}, {"hx_sum_sq", t.sum_sq},
            {"hx_xp", t.xp}, {"hx_zp", t.zp}, {"hx_sv", t.sv}, {"hx_yf", t.yf}, {"hx_yu", t.yu}, {"hx_u1", t.u1},
            {"hx_ub", t.ub}, {"hx_d1", t.d1}, {"hx_db", t.db}, {"hx_s1", t.s1}, {"hx_sb", t.sb},
            {"hx_dsigdz", t.dsigdz}, {"hx_b", t.b}, {"hx_Rhad", t.rhad}, {"hx_sighad", t.sighad}, {"hx_F1p", t.F1p},
            {"hx_F2p", t.F2p}, {"hx_F1n", t.F1n}, {"hx_F2n", t.F2n}, {"hx_slac", t.emc}, {"hx_F1", t.F1},
            {"hx_F2", t.F2}, {"hx_W1", t.W1}, {"hx_W2", t.W2}, {"hx_sin2th2", t.sin2th2}, {"hx_sige", t.sige},
            {"hx_sigsemi", t.sigsemi}, {"hx_jacobian", t.jacobian}, {"hx_fac", t.fac}, {"hx_peepix", res},
            {"hx_Qgev_after", t.Qgev_after}};
        const bool cent = iv(r, "doing_cent") != 0;
        for (const auto& kv : h) o[kv.first] = cent ? 0.0 : kv.second;
    });
}

// ------------------------------------------------------------------------------------------------- fortran_diff
// D-051 F1/F4/F6 (FORTRAN-BUG, corrected by default) on the peepix table rows: the corrected mode differs from the
// Fortran exactly where a correction applies. Compat = the Fortran (bitwise, as l1_peepix); corrected: peepiX changes
// only on the sum_sq = 0 rows (F1: NaN -> 0) and the pi0 rows (F4: dbar FF, Mpi0), survivalprob/decdist only on
// decay-OFF rows (F6: p_hadron, sqrt, zaero); every other output is unchanged. (F2 needs a negative x_Bj, which STOPs
// the Fortran: not in the table, covered by sidis_sanity.)
int runFortranDiff(const std::string& dir, const std::string& stateDir, const std::string& inpDirs,
                   const std::string& dataDir, const std::string& spectralDir) {
    const auto cfgs = readConfigs(stateDir);
    simcl1::Table in;
    std::string err;
    if (cfgs.empty() || !in.load(dir + "/input.txt", err)) return 77;
    sidis::Cteq5Pdf pdf0;
    sidis::FdssGrid fdss0;
    if (!pdf0.load(dataDir, 1, &err) || !fdss0.load(dataDir, 2, 1, &err)) return 77;
    struct Mode {
        sidis::PeepiXConfig cfg;
        sidis::PeepiXState st;
        sidis::Cteq5Pdf pdf;
    };
    std::map<int, std::array<Mode, 2>> procs;   // [0] compat, [1] corrected
    long nrow = 0, nan0 = 0, nan0chg = 0, pi0 = 0, pi0chg = 0, other_chg = 0, surv = 0, survchg = 0, surv_other = 0, zaero_rows = 0;
    std::ostringstream log;
    for (std::size_t i = 0; i < in.rows.size(); ++i) {
        Row r(in, i);
        const int cfg = iv(r, "cfg");
        auto it = procs.find(cfg);
        if (it == procs.end()) {
            const std::string inp = findInp(inpDirs, cfgs.at(cfg));
            if (inp.empty()) throw std::runtime_error("no .inp for cfg " + std::to_string(cfg));
            SimcParameters p;
            InpReader rd(p);
            if (!rd.loadRun(inp)) throw std::runtime_error("cannot read " + inp);
            std::array<Mode, 2> m;
            for (int k = 0; k < 2; ++k) {
                SetupOptions so;
                so.fortran_compat = k == 0;
                so.data_dir = spectralDir;
                const RunSetup st = BuildRunSetup(p, so);
                m[static_cast<std::size_t>(k)].cfg = sidis::makeConfig(st);
                m[static_cast<std::size_t>(k)].pdf = pdf0;
            }
            it = procs.emplace(cfg, std::move(m)).first;
        }
        double res[2], sp[2], dd[2];
        bool pizero = false, decayOff = false, thr = false, sumsq0 = false;
        int arm = 0;
        for (int k = 0; k < 2; ++k) {
            Mode& M = it->second[static_cast<std::size_t>(k)];
            sidis::PeepiXConfig c = M.cfg;
            if (iv(r, "ov_state") == 1) {
                c.doing_hplus = iv(r, "ov_hplus") != 0;
                c.doing_pizero = iv(r, "ov_pizero") != 0;
                c.do_fermi = iv(r, "ov_do_fermi") != 0;
                c.sigc_flag = iv(r, "ov_sigc_flag");
                c.sigc_nbin = iv(r, "ov_sigc_nbin");
                c.sigc_kin_min = r["ov_sigc_kin_min"];
                c.sigc_kin_max = r["ov_sigc_kin_max"];
                c.sigc_kin_ind = r["ov_sigc_kin_ind"];
                c.hadron_arm = iv(r, "ov_hadron_arm");
            }
            const Event v = readEvent(r, "v_");
            const Event v0 = readEvent(r, "v0_");
            MainEvent m;
            m.theta_pq = r["in_m_theta_pq"];
            m.FP.p.path = r["in_m_FP_p_path"];
            m.FP.p.dx = r["in_m_FP_p_dx"];
            m.FP.p.dy = r["in_m_FP_p_dy"];
            m.SP.p.delta = r["in_m_SP_p_delta"];
            m.ntup.decdist = r["decdist_in"];
            m.pfermi = {r["pfer"], r["pferx"], r["pfery"], r["pferz"], r["efer"]};
            sp[k] = r["in_survivalprob"];
            sidis::PeepiXTrace t;
            res[k] = sidis::peepiX(v, v0, m, sp[k], iv(r, "doing_cent") != 0, c, M.pdf,
                                   c.doing_semika ? &fdss0 : nullptr, M.st, &log, &t);
            dd[k] = m.ntup.decdist;
            pizero = c.doing_pizero && c.doing_semipi;
            decayOff = !c.doing_decay;
            if (k == 0) {
                thr = t.threshold_return;
                sumsq0 = !thr && t.sum_sq == 0.0;
            }
            arm = c.hadron_arm;
        }
        ++nrow;
        const bool peq = sameBits(res[0], res[1]);
        if (sumsq0) {
            ++nan0;
            if (std::isnan(res[0]) && res[1] == 0.0) ++nan0chg;
        } else if (pizero && !thr) {
            ++pi0;
            if (!peq) ++pi0chg;
        } else if (!peq) {
            ++other_chg;
        }
        if (decayOff && !thr) {
            ++surv;
            if (!sameBits(sp[0], sp[1]) && !sameBits(dd[0], dd[1])) ++survchg;
            if (arm == 5 || arm == 6) ++zaero_rows;
        } else if (!sameBits(sp[0], sp[1]) || !sameBits(dd[0], dd[1])) {
            ++surv_other;
        }
    }
    std::printf("fortran_diff D-051 (peepix rows %ld):\n", nrow);
    std::printf("  F1 sum_sq = 0: %ld rows, compat NaN -> corrected 0 on %ld\n", nan0, nan0chg);
    std::printf("  F4 pi0 FF/mass: %ld pi0 rows, peepiX changed on %ld; other rows changed: %ld\n", pi0, pi0chg,
                other_chg);
    std::printf("  F6 survival: %ld decay-OFF rows (%ld with hadron_arm 5/6), survivalprob and decdist changed on %ld; "
                "other rows changed: %ld\n", surv, zaero_rows, survchg, surv_other);
    const bool ok = nan0 > 0 && nan0chg == nan0 && pi0 > 0 && pi0chg == pi0 && other_chg == 0 && surv > 0 && survchg > 0 && surv_other == 0;
    std::printf("fortran_diff D-051 -> %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
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
        if (what == "ctq5pdf" && argc > 4) return runCtq5pdf(argv[2], argv[3], std::atof(argv[4]));
        if (what == "ctq5_state" && argc > 3) return runCtq5State(argv[2], argv[3]);
        if (what == "fdss" && argc > 5) return runFdss(argv[2], argv[3], argv[4], std::atof(argv[5]));
        if (what == "f1f2in21" && argc > 3) return runF1f2in21(argv[2], std::atof(argv[3]));
        if (what == "slacemcfit" && argc > 3) return runSlac(argv[2], std::atof(argv[3]));
        if (what == "fortran_diff" && argc > 6) return runFortranDiff(argv[2], argv[3], argv[4], argv[5], argv[6]);
        if (what == "peepix" && argc > 7)
            return runPeepix(argv[2], argv[3], argv[4], argv[5], argv[6], std::atof(argv[7]));
    } catch (const std::exception& e) {
        std::printf("ERROR: %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "bad arguments for mode %s\n", what.c_str());
    return 2;
}
