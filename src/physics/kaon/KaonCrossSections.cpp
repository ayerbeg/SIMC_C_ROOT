// src/physics/kaon/KaonCrossSections.cpp — see include/simc/physics/kaon/KaonCrossSections.h.
// The kaon-electroproduction cross sections of SIMC: peeK (lab cross section, kaon survival), the factorized model
// sig_factorized, the Saghai tables (eekeek/eekeeks, interpolated with CERNLIB fint), and the phase-space weight.
// Fortran: simc_gfortran/physics_kaon.f (peeK, sig_factorized, eekeek, eekeeks, phspwght), cern/fint.f, dbase.f:624-679.
// Operation order follows the Fortran source; x**2 is x*x (D-024); every literal is double (-fdefault-real-8, D-010).
#include "simc/physics/kaon/KaonCrossSections.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "simc/core/Errors.h"
#include "simc/core/SimcConstants.h"
#include "simc/material/FortranFormat.h"
#include "simc/reactions/EventChain.h"

namespace simc {
namespace kaon {

namespace {

inline double sq(double x) { return x * x; }
// gfortran MAX of two reals: mvar = a; if (b > mvar || isnan(mvar)) mvar = b (as physics_pion / PionProduction)
inline double fmaxF(double a, double b) { return (b > a || std::isnan(a)) ? b : a; }

// ------------------------------------------------------------------------------------------- binary128 (D-043 F4)
__extension__ typedef __float128 quad;
// glibc's hypotf128 (libm): gfortran's ABS of a COMPLEX(16) calls cabsf128 = hypotf128(re, im). Declared under its
// assembler name to stay independent of the C++ declarations of _Float128 in <math.h>.
extern "C" quad simc_hypotf128(quad, quad) __asm__("hypotf128");

/// COMPLEX(16) with the operations gfortran emits (-fcx-fortran-rules: plain products, no NaN recovery).
struct cq {
    quad re, im;
};
inline cq conj(cq a) { return {a.re, -a.im}; }
inline cq add(cq a, cq b) { return {a.re + b.re, a.im + b.im}; }
inline cq sub(cq a, cq b) { return {a.re - b.re, a.im - b.im}; }
// (ar+i ai)(br+i bi): rr = ar*br - ai*bi, ri = ar*bi + ai*br (tree-complex.cc, expand_complex_multiplication)
inline cq mul(cq a, cq b) { return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re}; }
// COMPLEX * REAL (the real operand is "only real": component-wise)
inline cq scal(cq a, quad x) { return {a.re * x, a.im * x}; }
// abs(z)**2
inline quad abs2(cq a) {
    const quad r = simc_hypotf128(a.re, a.im);
    return r * r;
}

/// The Saghai cross section from the 12 amplitudes zf(1:2, 1:6) (physics_kaon.f:317-348 = 449-474), in binary128.
double saghaiSigma(const double zf[2][6], double aflx, double aflxl, double x, double sx, double theta, double phi,
                   double epsi, double qr) {
    // z1a=zf(1,1)*ur+zf(2,1)*ui ...: ur = (1,0), ui = (0,1) are constants, so z = (zf(1,k), zf(2,k)) exactly.
    const cq z1{zf[0][0], zf[1][0]}, z2{zf[0][1], zf[1][1]}, z3{zf[0][2], zf[1][2]}, z4{zf[0][3], zf[1][3]};
    const cq z7{zf[0][4], zf[1][4]}, z8{zf[0][5], zf[1][5]};
    const quad X = x;
    const quad two = 2., half = 0.5;
    // conjg(z1a)*z4a-conjg(z2a)*z3a+conjg(z3a)*z4a*x
    const cq w = add(sub(mul(conj(z1), z4), mul(conj(z2), z3)), scal(mul(conj(z3), z4), X));
    // dsigt00 = aflx*(abs(z1a)**2+abs(z2a)**2+2.*real(conjg(z1a)*z2a)*x+0.5*sx**2*(abs(z3a)**2+abs(z4a)**2+2.*real(w)))
    const quad e = (abs2(z3) + abs2(z4)) + two * w.re;
    const quad inner = ((abs2(z1) + abs2(z2)) + (two * mul(conj(z1), z2).re) * X) + static_cast<quad>(0.5 * (sx * sx)) * e;
    const double dsigt00 = static_cast<double>(static_cast<quad>(aflx) * inner);
    // dsigl00 = aflxl*epsi*(abs(z7a)**2+abs(z8a)**2+2.*real(conjg(z7a)*z8a)*x)
    const double dsigl00 = static_cast<double>(static_cast<quad>(aflxl * epsi) *
                                               ((abs2(z7) + abs2(z8)) + (two * mul(conj(z7), z8).re) * X));
    // dsigp00 = aflx*epsi*(sin(theta))**2*cos(2.*phi)*(0.5*abs(z3a)**2+0.5*abs(z4a)**2+real(w))
    const double st = std::sin(theta);
    const double dsigp00 = static_cast<double>(static_cast<quad>(((aflx * epsi) * (st * st)) * std::cos(2. * phi)) *
                                               ((half * abs2(z3) + half * abs2(z4)) + w.re));
    // dsigi00 = aflx*sqrt(2.*qr**2*epsi*(1.+epsi))*sin(theta)*cos(phi)*
    //           real((z7a)*(conjg(z3a)-conjg(z2a)+conjg(z4a)*x)+z8a*(conjg(z1a)+conjg(z3a)*x+conjg(z4a)))
    const cq u = add(sub(conj(z3), conj(z2)), scal(conj(z4), X));
    const cq v = add(add(conj(z1), scal(conj(z3), X)), conj(z4));
    const quad r = add(mul(z7, u), mul(z8, v)).re;
    const double d = ((aflx * std::sqrt(((2. * (qr * qr)) * epsi) * (1. + epsi))) * st) * std::cos(phi);
    const double dsigi00 = static_cast<double>(static_cast<quad>(d) * r);
    return ((dsigt00 + dsigl00) + dsigp00) + dsigi00;
}

/// The kinematic prologue shared by eekeek/eekeeks (physics_kaon.f:264-280 = 377-393).
struct SaghaiKin {
    double w, skc, q0, q0c, qr, aflx, aflxl, an, x, sx;
};
SaghaiKin saghaiKin(double ss, double q22, double angl, double Mrec) {
    using constants::Mk2;
    using constants::Mp;
    using constants::Mp2;
    using constants::hbarc;
    using constants::pi;
    SaghaiKin k{};
    k.w = std::sqrt(ss) * 1000.;
    double skc2 = sq(sq(k.w) - Mk2 - sq(Mrec)) - 4. * Mk2 * sq(Mrec);
    skc2 = fmaxF(skc2, 0.e0);
    k.skc = std::sqrt(skc2) / 2. / k.w;
    k.q0 = -(-q22 - sq(k.w) + Mp2) / 2. / Mp;
    k.q0c = (-q22 + k.q0 * Mp) / k.w;
    k.qr = std::sqrt(q22) / k.q0c;
    // qv2, qv, qvc (l.271-273) are computed but unused
    k.aflx = k.skc / 2. / k.w / (sq(k.w) - Mp2) * (hbarc * hbarc) * 10000.;
    k.aflxl = k.aflx * (k.qr * k.qr);
    k.an = angl * 180. / pi;
    k.x = std::cos(angl);
    k.sx = std::sin(angl);
    return k;
}

/// The value the REAL*4 caller sees for `fint(...)`: +0.0 as built (D-043 F1), the REAL*4 rounding of the computed
/// value when corrected.
inline double fintSeen(int pn, const float* px, const int* pna, const float* pa, const float* table, bool fixed) {
    if (!fixed) return 0.0;   // FORTRAN-BUG(D-043 F1), see SaghaiOptions::fint_value
    return static_cast<double>(static_cast<float>(fint(pn, px, pna, pa, table)));
}

// ------------------------------------------------------------------------------------------- formatted input
/// A Fortran Ew.d / Fw.d field read into REAL*4 (PAD='YES', BLANK='NULL': blanks ignored, an empty field is 0, without a
/// decimal point d digits are implied). libgfortran converts the digits directly to the kind of the item (strtof).
float readEField(const std::string& rec, std::size_t pos, std::size_t w, int d) {
    std::string s;
    for (std::size_t i = pos; i < pos + w && i < rec.size(); ++i)
        if (rec[i] != ' ') s += rec[i];
    if (s.empty()) return 0.0f;
    for (auto& c : s)
        if (c == 'D' || c == 'd') c = 'E';
    if (s.find('.') != std::string::npos) return std::strtof(s.c_str(), nullptr);
    // implied decimal point: mantissa digits * 10**(exponent - d)
    std::size_t epos = s.find_first_of("eE");
    std::string mant = s.substr(0, epos);
    long ex = 0;
    if (epos != std::string::npos) ex = std::strtol(s.c_str() + epos + 1, nullptr, 10);
    const std::string t = mant + "e" + std::to_string(ex - d);
    return std::strtof(t.c_str(), nullptr);
}

/// A sequential formatted input file, one record per line (Fortran unit semantics for READ statements).
class Records {
public:
    explicit Records(const std::string& path) : f_(path) {}
    bool ok() const { return static_cast<bool>(f_); }
    long line() const { return line_; }
    bool next(std::string& rec) {
        if (!std::getline(f_, rec)) return false;
        if (!rec.empty() && rec.back() == '\r') rec.pop_back();
        ++line_;
        return true;
    }
    /// list-directed READ of n values (the statement starts a new record and continues on the next records until n
    /// values are read; the rest of the last record is skipped)
    bool listRead(int n, std::vector<std::string>& vals) {
        vals.clear();
        std::string rec;
        while (static_cast<int>(vals.size()) < n) {
            if (!next(rec)) return false;
            for (auto& c : rec)
                if (c == ',' || c == '\t') c = ' ';
            std::istringstream ss(rec);
            std::string tok;
            while (static_cast<int>(vals.size()) < n && ss >> tok) vals.push_back(tok);
        }
        return true;
    }
    /// READ(u, '(6e12.4)') of k REAL*4 items from one record
    bool fmtRead6(float* out, int k) {
        std::string rec;
        if (!next(rec)) return false;
        for (int i = 0; i < k; ++i) out[i] = readEField(rec, static_cast<std::size_t>(12 * i), 12, 4);
        return true;
    }

private:
    std::ifstream f_;
    long line_{0};
};

float listFloat(const std::string& tok) {
    std::string s = tok;
    for (auto& c : s)
        if (c == 'D' || c == 'd') c = 'E';
    return std::strtof(s.c_str(), nullptr);
}

std::string joinPath(const std::string& d, const std::string& f) { return d.empty() ? f : (d.back() == '/' ? d + f : d + "/" + f); }

}  // namespace

// ================================================================================================ CERNLIB fint
// Fortran: cern/fint.f (l.10-84). Implicit typing: X, H, ETA, WEIGHT, FINT are REAL(8) (-fdefault-real-8), ARG, ENT,
// TABLE are REAL*4. Arithmetic IF (l.44) as gfortran compiles it: x <= 0 ? (x < 0 ? 12 : 20) : 13 (NaN -> 13).
double fint(int narg, const float* arg, const int* nent, const float* ent, const float* table) {
    double result = 0.;
    if (narg < 1 || narg > 5) return result;   // l.21 GOTO 300: KERMTR('E104.1') (not reached by SIMC)
    int index[32];
    double weight[32];
    int lmax = 0, istep = 1, knots = 1;
    index[0] = 1;
    weight[0] = 1.;
    for (int n = 1; n <= narg; ++n) {   // DO 100
        const double x = arg[n - 1];
        const int ndim = nent[n - 1];
        int loca = lmax;
        const int lmin = lmax + 1;
        lmax = lmax + ndim;
        int ishift = 0;
        double eta = 0.;
        auto E = [&](int i) { return ent[i - 1]; };   // ENT(i), 1-based
        if (ndim > 2) {
            // 10-14: binary search
            int locb = lmax + 1;
            int locc = 0;
            bool exact = false;
            for (;;) {
                locc = (loca + locb) / 2;   // 11
                const double d = x - E(locc);
                if (d < 0) {
                    locb = locc;   // 12
                } else if (d == 0) {
                    exact = true;   // 20
                    break;
                } else {
                    loca = locc;   // 13
                }
                if (!(locb - loca > 1)) break;   // 14
            }
            if (exact) {
                ishift = (locc - lmin) * istep;   // 20
                for (int k = 1; k <= knots; ++k) index[k - 1] = index[k - 1] + ishift;   // 21-22
                istep = istep * ndim;   // 90
                continue;
            }
            loca = std::min(std::max(loca, lmin), lmax - 1);
            ishift = (loca - lmin) * istep;
            eta = (x - E(loca)) / static_cast<double>(static_cast<float>(E(loca + 1) - E(loca)));
        } else {
            if (ndim == 1) continue;   // GOTO 100
            const double h = x - E(lmin);
            if (h == 0.) {
                istep = istep * ndim;   // GOTO 90
                continue;
            }
            ishift = istep;
            if (x - E(lmin + 1) == 0.) {   // GOTO 21
                for (int k = 1; k <= knots; ++k) index[k - 1] = index[k - 1] + ishift;
                istep = istep * ndim;
                continue;
            }
            ishift = 0;
            eta = h / static_cast<double>(static_cast<float>(E(lmin + 1) - E(lmin)));
        }
        // 30-31
        for (int k = 1; k <= knots; ++k) {
            index[k - 1] = index[k - 1] + ishift;
            index[k + knots - 1] = index[k - 1] + istep;
            weight[k + knots - 1] = weight[k - 1] * eta;
            weight[k - 1] = weight[k - 1] - weight[k + knots - 1];
        }
        knots = 2 * knots;
        istep = istep * ndim;   // 90
    }
    for (int k = 1; k <= knots; ++k) {   // DO 200
        const int i = index[k - 1];
        result = result + weight[k - 1] * static_cast<double>(table[i - 1]);
    }
    return result;
}

// ================================================================================================ /sigkaon/
const std::array<const char*, 4> SaghaiTables::kFiles = {"weightc.dat", "weightd.dat", "saghai_proton.dat",
                                                          "saghai_sigma0.dat"};

// Fortran: dbase.f:624-679
bool SaghaiTables::load(const std::string& dir, bool sigma0Corrected, std::vector<std::string>* msgs, std::string* err) {
    auto fail = [&](const std::string& m) {
        if (err) *err = m;
        return false;
    };
    std::vector<std::string> v;
    // weightc.dat: do i=1,50; do j=1,20; read(1,*) weightc(j,i)
    {
        const std::string p = joinPath(dir, kFiles[0]);
        Records r(p);
        if (!r.ok()) return fail("cannot open " + p);
        std::vector<float> a(20 * 50);
        for (int i = 1; i <= 50; ++i)
            for (int j = 1; j <= 20; ++j) {
                if (!r.listRead(1, v)) return fail(p + ": end of file at weightc(" + std::to_string(j) + "," + std::to_string(i) + ")");
                a[static_cast<std::size_t>((i - 1) * 20 + (j - 1))] = listFloat(v[0]);
            }
        weightc.swap(a);
    }
    // weightd.dat: do i=1,30; do j=1,40; do k=1,8; read(1,*) weightd(k,j,i)
    {
        const std::string p = joinPath(dir, kFiles[1]);
        Records r(p);
        if (!r.ok()) return fail("cannot open " + p);
        std::vector<float> a(8 * 40 * 30);
        for (int i = 1; i <= 30; ++i)
            for (int j = 1; j <= 40; ++j)
                for (int k = 1; k <= 8; ++k) {
                    if (!r.listRead(1, v)) return fail(p + ": end of file");
                    a[static_cast<std::size_t>(((i - 1) * 40 + (j - 1)) * 8 + (k - 1))] = listFloat(v[0]);
                }
        weightd.swap(a);
    }
    // saghai_proton.dat (l.646-661): per (iread, iq2) a list-directed header, then per angle 3 records (6e12.4): the
    // 5-value kinematics line (discarded) and the two amplitude lines.
    {
        const std::string p = joinPath(dir, kFiles[2]);
        Records r(p);
        if (!r.ok()) return fail("cannot open " + p);
        for (int k = 0; k < 6; ++k) {
            zrff[static_cast<std::size_t>(k)].assign(kNsP * kNqP * kNa, 0.f);
            ziff[static_cast<std::size_t>(k)].assign(kNsP * kNqP * kNa, 0.f);
        }
        float d5[6], a[6], b[6];
        for (int ir = 1; ir <= kNsP; ++ir)
            for (int iq = 1; iq <= kNqP; ++iq) {
                if (!r.listRead(2, v)) return fail(p + ": end of file (header)");
                for (int ia = 1; ia <= kNa; ++ia) {
                    if (!r.fmtRead6(d5, 5) || !r.fmtRead6(a, 6) || !r.fmtRead6(b, 6)) return fail(p + ": end of file");
                    const std::size_t i = idxP(ir, iq, ia);
                    for (int k = 0; k < 3; ++k) {
                        zrff[static_cast<std::size_t>(k)][i] = a[2 * k];
                        ziff[static_cast<std::size_t>(k)][i] = a[2 * k + 1];
                        zrff[static_cast<std::size_t>(k + 3)][i] = b[2 * k];
                        ziff[static_cast<std::size_t>(k + 3)][i] = b[2 * k + 1];
                    }
                }
            }
    }
    // saghai_sigma0.dat (l.662-677): per (iread, iq2, iang) a list-directed header and 3 records (6e12.4).
    {
        const std::string p = joinPath(dir, kFiles[3]);
        Records r(p);
        if (!r.ok()) return fail("cannot open " + p);
        for (int k = 0; k < 6; ++k) {
            zsrff[static_cast<std::size_t>(k)].assign(kNsS * kNqS * kNa, 0.f);
            zsiff[static_cast<std::size_t>(k)].assign(kNsS * kNqS * kNa, 0.f);
        }
        // the eekeeks grid (physics_kaon.f:402-431) for the header check of the corrected reader
        static const double sNode[kNsS] = {2.851, 2.898, 2.945, 2.991, 3.038, 3.085, 3.132, 3.320, 3.507, 3.695,
                                           3.883, 4.070, 4.258, 4.446, 4.633, 4.821, 5.009, 5.196, 5.384, 5.572};
        static const double qNode[kNqS] = {0.0, 0.250, 0.376, 0.520, 0.750, 1.000, 1.250, 1.500, 1.750, 2.000};
        std::string first;
        if (sigma0Corrected) {
            // FORTRAN-BUG(D-043 F2): original Fortran dbase.f:662-677 reads 4-line records from line 1, but line 1 of
            // saghai_sigma0.dat is a truncated fragment ('-01  0.2000E-03 -0.4000E-02  0.1900E-02'), so every record
            // is read one line off: zsrff1..zsiff3 get the 5-value kinematics line, zsrff4..zsiff6 the true ff1..ff3,
            // and the true ff4..ff6 are lost; corrected because the amplitudes must come from their own lines: the
            // fragment is skipped and each record is read as header ("s Q2"), kinematics line, ff1..ff3, ff4..ff6.
            if (!r.next(first)) return fail(p + ": empty file");
        }
        long badHeader = 0;
        float d5[6], a[6], b[6];
        for (int ir = 1; ir <= kNsS; ++ir)
            for (int iq = 1; iq <= kNqS; ++iq)
                for (int ia = 1; ia <= kNa; ++ia) {
                    if (!r.listRead(2, v)) return fail(p + ": end of file (header)");
                    if (!r.fmtRead6(d5, 5) || !r.fmtRead6(a, 6) || !r.fmtRead6(b, 6)) return fail(p + ": end of file");
                    const std::size_t i = idxS(ir, iq, ia);
                    for (int k = 0; k < 3; ++k) {
                        zsrff[static_cast<std::size_t>(k)][i] = a[2 * k];
                        zsiff[static_cast<std::size_t>(k)][i] = a[2 * k + 1];
                        zsrff[static_cast<std::size_t>(k + 3)][i] = b[2 * k];
                        zsiff[static_cast<std::size_t>(k + 3)][i] = b[2 * k + 1];
                    }
                    if (sigma0Corrected) {
                        const double s = std::strtod(v[0].c_str(), nullptr), q = std::strtod(v[1].c_str(), nullptr);
                        if (std::fabs(s - sNode[ir - 1]) > 1e-4 || std::fabs(q - qNode[iq - 1]) > 1e-4) ++badHeader;
                    }
                }
        // FORTRAN-DIFF-CANDIDATE(D-043 F2): the file content itself starts with 18 stray records (s = 2.898 GeV^2,
        // Q2 = 0.52/0.75) and has 1 record instead of 19 for s = 2.851, Q2 = 0: the 18 slots (iread 1, iq2 1,
        // angles 1-18) keep the stray records as the file orders them. Question: which amplitudes belong to s = 2.851,
        // Q2 = 0? (only the saghai column near the Sigma threshold at Q2 < 0.25 GeV^2 is affected)
        if (sigma0Corrected && badHeader > 0 && msgs)
            msgs->push_back("saghai_sigma0.dat: " + std::to_string(badHeader) +
                            " records do not carry the (s, Q2) of their table node (the stray records at the start of "
                            "the file, D-043 F2); kept as the file orders them");
    }
    return true;
}

// ================================================================================================ Saghai model
// Fortran: physics_kaon.f:eekeek (l.241-355)
double eekeek(double ss, double q22, double angl, double theta, double phi, double epsi, double Mrec_struck,
              const SaghaiTables& tab, const SaghaiOptions& opt, SaghaiTrace* trace) {
    const SaghaiKin k = saghaiKin(ss, q22, angl, Mrec_struck);
    const int pn = 3;
    float px[3];
    px[0] = static_cast<float>(ss);
    px[1] = static_cast<float>(static_cast<quad>(q22) / static_cast<quad>(1.e6));   // q22/1.d+06 in REAL(16)
    px[2] = static_cast<float>(k.an);
    const int pna[3] = {10, 11, 19};
    float pa[50] = {};
    double ps = 2.6, qs = 0.0, as = 0.0;
    for (int i = 1; i <= 10; ++i) {
        pa[i - 1] = static_cast<float>(ps);
        ps = ps + 0.3;
    }
    if (opt.proton_grid_fixed) {
        // FORTRAN-BUG(D-043 F3): original Fortran physics_kaon.f:288-294: pa(i) = 2.6 + 0.3*(i-1) gives the nodes
        // 5.0 and 5.3 GeV^2 for i = 9, 10; corrected because the 9th and 10th blocks of saghai_proton.dat (read as
        // iread 9, 10 by dbase.f:646-661) are s = 5.1 and 5.4 GeV^2.
        pa[8] = static_cast<float>(5.1);
        pa[9] = static_cast<float>(5.4);
    }
    for (int i = 11; i <= 21; ++i) {
        pa[i - 1] = static_cast<float>(qs);
        qs = qs + 0.2;
    }
    for (int i = 22; i <= 40; ++i) {
        pa[i - 1] = static_cast<float>(as);
        as = as + 10.;
    }
    double zf[2][6];
    for (int j = 0; j < 6; ++j) {
        zf[0][j] = fintSeen(pn, px, pna, pa, tab.zrff[static_cast<std::size_t>(j)].data(), opt.fint_value);
    }
    for (int j = 0; j < 6; ++j) {
        zf[1][j] = fintSeen(pn, px, pna, pa, tab.ziff[static_cast<std::size_t>(j)].data(), opt.fint_value);
    }
    if (trace) {
        for (int i = 0; i < 3; ++i) trace->px[i] = px[i];
        for (int j = 0; j < 6; ++j) {
            trace->zr[j] = zf[0][j];
            trace->zi[j] = zf[1][j];
        }
    }
    return saghaiSigma(zf, k.aflx, k.aflxl, k.x, k.sx, theta, phi, epsi, k.qr);
}

// Fortran: physics_kaon.f:eekeeks (l.357-477)
double eekeeks(double ss, double q22, double angl, double theta, double phi, double epsi, double Mrec_struck,
               const SaghaiTables& tab, const SaghaiOptions& opt, SaghaiTrace* trace) {
    const SaghaiKin k = saghaiKin(ss, q22, angl, Mrec_struck);
    const int pn = 3;
    float px[3];
    px[0] = static_cast<float>(ss);
    px[1] = static_cast<float>(static_cast<quad>(q22) / static_cast<quad>(1.e6));
    px[2] = static_cast<float>(k.an);
    const int pna[3] = {20, 10, 19};
    float pa[50] = {};
    static const double sq0[30] = {2.851, 2.898, 2.945, 2.991, 3.038, 3.085, 3.132, 3.320, 3.507, 3.695,
                                   3.883, 4.070, 4.258, 4.446, 4.633, 4.821, 5.009, 5.196, 5.384, 5.572,
                                   0.0,   0.250, 0.376, 0.520, 0.750, 1.000, 1.250, 1.500, 1.750, 2.000};
    for (int i = 1; i <= 30; ++i) pa[i - 1] = static_cast<float>(sq0[i - 1]);
    double as = 0.0;
    for (int i = 31; i <= 49; ++i) {
        pa[i - 1] = static_cast<float>(as);
        as = as + 10.;
    }
    double zf[2][6];
    for (int j = 0; j < 6; ++j) {
        zf[0][j] = fintSeen(pn, px, pna, pa, tab.zsrff[static_cast<std::size_t>(j)].data(), opt.fint_value);
    }
    for (int j = 0; j < 6; ++j) {
        zf[1][j] = fintSeen(pn, px, pna, pa, tab.zsiff[static_cast<std::size_t>(j)].data(), opt.fint_value);
    }
    if (trace) {
        for (int i = 0; i < 3; ++i) trace->px[i] = px[i];
        for (int j = 0; j < 6; ++j) {
            trace->zr[j] = zf[0][j];
            trace->zi[j] = zf[1][j];
        }
    }
    return saghaiSigma(zf, k.aflx, k.aflxl, k.x, k.sx, theta, phi, epsi, k.qr);
}

// Fortran: physics_kaon.f:phspwght (l.479-544)
PhspwghtResult phspwght(double delta, double theta, double phi, int electron_arm, int hadron_arm,
                        const SaghaiTables& tab, const SaghaiOptions& opt, std::ostream* log) {
    PhspwghtResult res;
    float pa[78] = {};
    double weight;
    if (electron_arm == 1 && hadron_arm == 2) {   // e- in HMS, K in SOS
        const int pn = 2;
        const float px[2] = {static_cast<float>(theta), static_cast<float>(phi)};
        const int pna[2] = {20, 50};
        double ts = -0.03325, ps = -0.0735;
        for (int i = 1; i <= 20; ++i) {
            pa[i - 1] = static_cast<float>(ts);
            ts = ts + 0.0035;
        }
        for (int i = 21; i <= 70; ++i) {
            pa[i - 1] = static_cast<float>(ps);
            ps = ps + 0.003;
        }
        weight = fintSeen(pn, px, pna, pa, tab.weightc.data(), opt.fint_value);
    } else if (electron_arm == 2 && hadron_arm == 1) {   // e- in SOS, K in HMS
        const int pn = 3;
        const float px[3] = {static_cast<float>(delta), static_cast<float>(theta), static_cast<float>(phi)};
        const int pna[3] = {8, 40, 30};
        double ts = -0.063375, ps = -0.0435, ds = -17.5;
        for (int i = 1; i <= 8; ++i) {
            pa[i - 1] = static_cast<float>(ds);
            ds = ds + 5.;
        }
        for (int i = 9; i <= 48; ++i) {
            pa[i - 1] = static_cast<float>(ts);
            ts = ts + 0.00325;
        }
        for (int i = 49; i <= 78; ++i) {
            pa[i - 1] = static_cast<float>(ps);
            ps = ps + 0.003;
        }
        weight = fintSeen(pn, px, pna, pa, tab.weightd.data(), opt.fint_value);
    } else {
        if (log) {
            char b[96];
            std::snprintf(b, sizeof b, " electron_arm=%12d  and hadron_arm=%12d", electron_arm, hadron_arm);
            *log << b << "\n"
                 << " eekeek.f has a phase space factor that is only defined for\n"
                 << " hms&sos case.  Need to update for other spectrometers.\n";
        }
        throw FortranStop("physics_kaon.f:539 phspwght: phase space factor only defined for the HMS/SOS case (stop)");
    }
    res.fint_seen = weight;
    weight = fmaxF(weight, 0.01e00);
    weight = fmaxF(100.e00 / weight, 1.0e00);
    res.weight = weight;
    return res;
}

// ================================================================================================ factorized model
// Fortran: physics_kaon.f:sig_factorized (l.173-237)
double sigFactorized(double q2, double w, double t, double pk, double mrec, SigFactorizedState& st, std::ostream* log) {
    using constants::Mk2;
    using constants::Mp;
    // FORTRAN-DIFF-CANDIDATE(D-043 F10): Mp for every target (Sigma- from the neutron, coherent production).
    const double nu = (sq(w) + q2 - sq(Mp)) / 2. / Mp;
    const double q = std::sqrt(q2 + sq(nu));
    const double qcm = q * (Mp / w);
    const double nucm = std::sqrt(sq(qcm) - q2);
    const double tmin = -1. * (Mk2 - q2 - 2 * nucm * std::sqrt(sq(pk) + Mk2) + 2 * qcm * pk);

    // Check center of mass pk, since we can get it from w2 (l.201-206)
    // FORTRAN-DIFF-CANDIDATE(D-043 F12): spurious for the central event of Sigma0/Sigma-/A>1 quasifree runs.
    const double pktest = std::sqrt(sq(sq(w) + Mk2 - sq(mrec)) / 4. / sq(w) - Mk2);
    if (std::fabs((pktest - pk) / pk) > 0.001 && log) {
        *log << " Kaon C.M. momentum passed to sig_factorized does not agree with\n"
             << " the value calculated from W\n"
             << " Passed,calculated=" << fortran::listReal(pk) << fortran::listReal(pktest) << "\n";
    }

    const double q2val = q2 / 1.e6;
    const double w2val = sq(w) / 1.e6;
    const double pkval = pk / 1000.;
    const double tval = t / 1.e6;
    const double tminval = tmin / 1.e6;
    double fact_q, fact_t;
    double& fact_w = st.fact_w;   // SAVEd: assigned only when w**2 != 0
    if (mrec < 1150.) {   // lambda production
        fact_q = 1. / sq(q2val + 2.67);
        fact_t = std::exp(-2.1 * (tval - tminval));
        if (w2val != 0) {   // =0 for central event!!!
            fact_w = 0.959 * 4.1959 * pkval / (std::sqrt(w2val) * (w2val - 0.93827 * 0.93827));
            fact_w = fact_w + (0.18 * (1.72 * 1.72) * (0.10 * 0.10)) /
                                  (sq(w2val - 1.72 * 1.72) + (1.72 * 1.72) * (0.10 * 0.10));
        }
    } else {
        fact_q = 1. / sq(q2val + 0.79);
        fact_t = std::exp(-1.0 * (tval - tminval));
        if (w2val != 0) fact_w = 0.959 * 4.1959 * pkval / (std::sqrt(w2val) * (w2val - 0.93827 * 0.93827));
    }
    return fact_q * fact_t * fact_w;
}

// ================================================================================================ peeK
// Fortran: physics_kaon.f:peeK (l.1-168)
double peeK(const Event& vertex, MainEvent& main, double& survivalprob, const PeekConfig& cfg, const SaghaiTables& tab,
            PeekState& st, std::ostream* log) {
    using constants::pi;
    const FermiMomentum& f = main.pfermi;   // COMMON /pfermi/

    // l.41-59: photon-NUCLEON C.M. frame (jacobians.f)
    const reactions::CmTransform cm =
        reactions::transformToCm(vertex, main, f.pfer, f.pferx, f.pfery, f.pferz, f.efer, log);
    const double pkcm = cm.phadcm;
    const double qstar = cm.qstar;
    // dOmega_cm --> dt dphi_cm jacobian 1/(2 pk_cm q_cm)
    const double jacobian = cm.jacobian / (2. * pkcm * qstar);
    const double jac_old = cm.jac_old / (2. * pkcm * qstar);
    main.thetacm = cm.thetacm;
    main.phicm = cm.phicm;
    main.pcm = pkcm;
    main.davejac = jacobian;
    main.johnjac = jac_old;   // approx. assuming collinear boost.

    // l.66-71
    double tfcos = f.pferx * vertex.uq.x + f.pfery * vertex.uq.y + f.pferz * vertex.uq.z;
    if (tfcos - 1. > 0. && tfcos - 1. < 1.e-8) tfcos = 1.0;
    const double tfsin = std::sqrt(1. - sq(tfcos));
    const double s = sq(vertex.nu + f.efer) - sq(vertex.q + f.pfer * tfcos) - sq(f.pfer * tfsin);
    main.wcm = std::sqrt(s);

    // l.82-87: Saghai model. FORTRAN-DIFF-CANDIDATE(D-043 F8): phi is a never-assigned SAVEd local (0), the theta
    // argument is main%theta_pq. FORTRAN-DIFF-CANDIDATE(D-043 F9): coherent Lambda (Mrec_struck > 1150) uses eekeeks.
    double sigma_saghai;
    if (cfg.Mrec_struck < 1150.)   // lambda production.
        sigma_saghai = eekeek(s / 1.e6, vertex.Q2, main.thetacm, main.theta_pq, st.phi, main.epsilon, cfg.Mrec_struck, tab,
                              cfg.saghai);
    else
        sigma_saghai = eekeeks(s / 1.e6, vertex.Q2, main.thetacm, main.theta_pq, st.phi, main.epsilon, cfg.Mrec_struck,
                               tab, cfg.saghai);
    main.ntup.sigcm1 = sigma_saghai;

    // l.92-98: factorization model (the model used)
    main.ntup.sigcm2 = sigFactorized(vertex.Q2, main.wcm, main.t, pkcm, cfg.Mrec_struck, st.sf, log);
    const double sigma_eek = main.ntup.sigcm2;
    main.ntup.sigcm = sigma_eek;

    // l.108: k_eq from wcm**2 (JRA 2001)
    const double k_eq = (sq(main.wcm) - sq(cfg.Mtar_struck)) / 2. / cfg.Mtar_struck;
    // l.132-135
    const double fac = 1. / (1. - f.pferz * f.pfer / f.efer) * cfg.Mtar_struck / f.efer;
    const double gtpr = constants::alpha / 2. / (pi * pi) * vertex.e.E / vertex.Ein * k_eq / vertex.Q2 / (1. - main.epsilon);
    const double peek = sigma_eek * jacobian * (gtpr * fac);   // ub/MeV^2/rad-->ub/sr-->ub/MeV/sr

    // l.148-165: survival probability when the decay is not simulated
    if (!cfg.doing_decay) {
        if (cfg.hadron_arm == 1) {
            st.zaero = 0.;   // no aerogel yet, use full length.
        } else if (cfg.hadron_arm == 2) {
            st.zaero = -82.8;   // From Rick: aero at 263.2cm,last project=346(cal).
        } else if (cfg.hadron_arm == 3) {
            st.zaero = -183.;   // aero at 130cm,last project=313(S2)
        } else if (cfg.hadron_arm == 4) {
            st.zaero = -183.;
        } else if (cfg.survival_path_fixed) {
            // FORTRAN-BUG(D-043 F6): original Fortran physics_kaon.f:149-158 has no branch for hadron_arm 5 (SHMS): zaero
            // keeps the SAVEd value of an earlier call (0 in a SIMC run, -82.8/-183 after an arm 2/3/4 call); corrected
            // because the path must not depend on earlier calls: 0 (full length to the last plane, as for the HMS).
            st.zaero = 0.;
        }
        double pathlen;
        if (cfg.survival_path_fixed)
            // FORTRAN-BUG(D-043 F6): original Fortran physics_kaon.f:159: zaero*(1+dx**2+dy**2); corrected because the
            // path length along the track is zaero*sqrt(1+dx**2+dy**2) (no effect for the HMS/SHMS, zaero = 0).
            pathlen = main.FP.p.path + st.zaero * std::sqrt(1 + sq(main.FP.p.dx) + sq(main.FP.p.dy));
        else
            pathlen = main.FP.p.path + st.zaero * (1 + sq(main.FP.p.dx) + sq(main.FP.p.dy));
        const double p_kaon = cfg.spec_p_P * (1. + main.SP.p.delta / 100.);
        double betak;
        if (cfg.survival_p_kaon)
            // FORTRAN-BUG(D-043 F5): original Fortran physics_kaon.f:161: betak = spec%p%P/sqrt(spec%p%P**2+Mh2) (the
            // central momentum; p_kaon is computed and unused, against the comment of l.145-146); corrected because
            // the decay probability depends on the momentum of this kaon after radiation/energy loss.
            betak = p_kaon / std::sqrt(sq(p_kaon) + cfg.Mh2);
        else
            betak = cfg.spec_p_P / std::sqrt(sq(cfg.spec_p_P) + cfg.Mh2);
        const double gammak = 1. / std::sqrt(1. - sq(betak));
        survivalprob = 1. / std::exp(pathlen / (cfg.ctau * betak * gammak));
        main.ntup.decdist = survivalprob;   // decdist in ntuple
    }
    return peek;
}

}  // namespace kaon
}  // namespace simc
