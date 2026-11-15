// src/physics/sidis/Fdss.cpp — fdss/fdss.f (the de Florian-Sassot-Stratmann fragmentation functions, grid +
// interpolation in z and Q2). See include/simc/physics/sidis/Fdss.h.
// Operation order follows the Fortran source (D-024); the binary128 parts of fdss.o (REAL(16) literals under
// -fdefault-real-8, _gfortran_pow_r16_i4, __multf3, __subtf3) are reproduced with __float128.
#include "simc/physics/sidis/Fdss.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

#include "simc/core/Errors.h"

namespace simc {
namespace sidis {

namespace {
__extension__ typedef __float128 quad;

// libgfortran pow_r16_i4 (generated/pow_r16_i4.c): REAL(16)**INTEGER(4) by repeated squaring. fdss.o calls it for
// (1.D0-X)**4 and (1.D0-X)**7 (gcc does not expand powi for binary128).
quad powR16I4(quad a, int b) {
    quad pow = 1, x = a;
    int n = b;
    unsigned u;
    if (n != 0) {
        if (n < 0) {
            u = static_cast<unsigned>(-n);
            x = pow / x;
        } else {
            u = static_cast<unsigned>(n);
        }
        for (;;) {
            if (u & 1u) pow *= x;
            u >>= 1;
            if (u)
                x *= x;
            else
                break;
        }
    }
    return pow;
}

// One Ew.d input field (fdss.f:90 FORMAT (9(1PE10.3)): w = 10, d = 3, scale factor k = 1), gfortran rules with the
// default BLANK='NULL' (blanks ignored): an explicit exponent overrides the scale factor; without a decimal point the
// last d digits are the fraction; without an exponent the value is divided by 10**k. Returns false for a bad field.
bool readEField(const std::string& field, double& v) {
    std::string s;
    for (char c : field)
        if (c != ' ') s += c;
    if (s.empty()) {
        v = 0;
        return true;
    }
    std::size_t i = 0;
    std::string sign;
    if (s[i] == '+' || s[i] == '-') sign = s[i++] == '-' ? "-" : "";
    std::string digits, frac;
    bool point = false;
    while (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.')) {
        if (s[i] == '.') {
            if (point) return false;
            point = true;
        } else if (point) {
            frac += s[i];
        } else {
            digits += s[i];
        }
        ++i;
    }
    if (digits.empty() && frac.empty()) return false;
    long expo = 0;
    bool hasExp = false;
    if (i < s.size()) {
        if (s[i] == 'E' || s[i] == 'e' || s[i] == 'D' || s[i] == 'd' || s[i] == 'Q' || s[i] == 'q') ++i;
        if (i >= s.size()) return false;
        std::size_t j = i;
        if (s[j] == '+' || s[j] == '-') ++j;
        if (j >= s.size()) return false;
        for (std::size_t k = j; k < s.size(); ++k)
            if (!std::isdigit(static_cast<unsigned char>(s[k]))) return false;
        expo = std::strtol(s.c_str() + i, nullptr, 10);
        hasExp = true;
    }
    if (!point) {   // implied decimal point: the last d = 3 digits are the fraction
        std::string all = digits;
        while (all.size() < 3) all = "0" + all;
        frac = all.substr(all.size() - 3);
        digits = all.substr(0, all.size() - 3);
    }
    if (!hasExp) expo -= 1;   // scale factor 1P on input without an exponent
    const std::string txt = sign + (digits.empty() ? "0" : digits) + "." + frac + "e" + std::to_string(expo);
    v = std::strtod(txt.c_str(), nullptr);
    return true;
}

bool fail(std::string* err, const std::string& m) {
    if (err) *err = m;
    return false;
}

double q16(long n, long d) { return static_cast<double>(static_cast<quad>(n) / static_cast<quad>(d)); }
}  // namespace

// fdss.f:57-63
// The REAL(16) literals are rounded to binary128 by the compiler, then to REAL(8) by the DATA statement: q(n, d) = n/d
// in binary128 (a correctly rounded division of exact integers = the correctly rounded literal), then to double.
const std::array<double, FdssGrid::NQ> FdssGrid::kQS = {
    q16(1, 1),     q16(125, 100), q16(15, 10),  q16(25, 10),  q16(4, 1),     q16(64, 10),   q16(10, 1),    q16(15, 1),
    q16(25, 1),    q16(40, 1),    q16(64, 1),   q16(100, 1),  q16(180, 1),   q16(320, 1),   q16(580, 1),   q16(1000, 1),
    q16(1800, 1),  q16(3200, 1),  q16(5800, 1), q16(10000, 1), q16(18000, 1), q16(32000, 1), q16(58000, 1), q16(100000, 1)};
const std::array<double, FdssGrid::NX> FdssGrid::kXB = {
    0.01, 0.02,  0.03,  0.04, 0.05,  0.06, 0.07,  0.08,  0.09, 0.095, 0.1,  0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
    0.3,  0.325, 0.35,  0.375, 0.4,  0.45, 0.5,   0.55,  0.6,  0.65,  0.7,  0.75,  0.8,  0.85,  0.9, 0.93,  1.0};

const char* FdssGrid::gridFile(int ih, int io) {
    static const char* const names[4][2] = {{"PILO.GRID", "PINLO.GRID"}, {"KALO.GRID", "KANLO.GRID"},
                                            {"PROLO.GRID", "PRONLO.GRID"}, {"HLO.GRID", "HNLO.GRID"}};
    if (ih < 1 || ih > 4 || io < 0 || io > 1) return nullptr;
    return names[ih - 1][io];
}

// Fortran: fdss.f:fDSS l.70-139 (FINI = 0)
bool FdssGrid::load(const std::string& dir, int ih, int io, std::string* err) {
    const char* g = gridFile(ih, io);
    if (!g) return fail(err, "fdss.f:96-98: WRONG SET (IH = " + std::to_string(ih) + ", IO = " + std::to_string(io) + ")");
    const std::string path = dir + "/" + g;
    std::ifstream in(path);
    if (!in) return fail(err, "cannot open " + path + " (fdss.f:72-95)");
    parton.assign(static_cast<std::size_t>(NPART * NQ * (NX - 1)), 0.0);
    std::string line;
    // DO 15 M = 1, NX-1; DO 15 N = 1, NQ: READ(IIREAD,90) PARTON(1..9,N,M)
    for (int m = 1; m <= NX - 1; ++m) {
        for (int n = 1; n <= NQ; ++n) {
            if (!std::getline(in, line))
                return fail(err, path + ": end of file at record " + std::to_string((m - 1) * NQ + n) +
                                     " (the Fortran READ fails)");
            if (!line.empty() && line.back() == '\r') line.pop_back();
            for (int i = 1; i <= NPART; ++i) {
                const std::size_t b = static_cast<std::size_t>((i - 1) * 10);
                const std::string f = b < line.size() ? line.substr(b, 10) : std::string();   // PAD='YES'
                double v = 0;
                if (!readEField(f, v))
                    return fail(err, path + ": bad E10.3 field '" + f + "' in record " +
                                         std::to_string((m - 1) * NQ + n));
                parton[static_cast<std::size_t>(((m - 1) * NQ + (n - 1)) * NPART + (i - 1))] = v;
            }
        }
    }
    auto P = [&](int i, int n, int m) { return parton[static_cast<std::size_t>(((m - 1) * NQ + (n - 1)) * NPART + (i - 1))]; };
    for (auto& a : F) a.assign(static_cast<std::size_t>(NX * NQ), 0.0);
    auto at = [](int ix, int iq) { return static_cast<std::size_t>((iq - 1) * NX + (ix - 1)); };
    // ARRAYS FOR THE INTERPOLATION SUBROUTINE (l.112-131)
    for (int iq = 1; iq <= NQ; ++iq) {
        for (int ix = 1; ix <= NX - 1; ++ix) {
            const double xb0 = kXB[static_cast<std::size_t>(ix - 1)];
            const double xb1 = static_cast<double>(static_cast<quad>(1) - static_cast<quad>(xb0));   // 1.D0-XB(IX)
            const double x2 = xb1 * xb1, x3 = xb1 * x2, x4 = x2 * x2, x7 = x3 * x4;   // XB1**4, XB1**7 (powi)
            const double p05 = std::pow(xb0, 0.5), p03 = std::pow(xb0, 0.3);
            F[0][at(ix, iq)] = P(1, iq, ix) / (x4 * p05);   // XUTOTF
            F[1][at(ix, iq)] = P(2, iq, ix) / (x4 * p05);   // XDTOTF
            F[2][at(ix, iq)] = P(3, iq, ix) / (x4 * p05);   // XSTOTF
            F[3][at(ix, iq)] = P(4, iq, ix) / (x7 * p03);   // XCTOTF
            F[4][at(ix, iq)] = P(5, iq, ix) / (x7 * p03);   // XBTOTF
            F[5][at(ix, iq)] = P(6, iq, ix) / (x4 * p03);   // XGF
            F[6][at(ix, iq)] = P(7, iq, ix) / (x4 * p05);   // XUVALF
            F[7][at(ix, iq)] = P(8, iq, ix) / (x4 * p05);   // XDVALF
            F[8][at(ix, iq)] = P(9, iq, ix) / (x4 * p05);   // XSVALF
        }
        for (auto& a : F) a[at(NX, iq)] = 0.;
    }
    NA[0] = NX;
    NA[1] = NQ;
    for (int ix = 1; ix <= NX; ++ix) ARRF[static_cast<std::size_t>(ix - 1)] = std::log(kXB[static_cast<std::size_t>(ix - 1)]);
    for (int iq = 1; iq <= NQ; ++iq)
        ARRF[static_cast<std::size_t>(NX + iq - 1)] = std::log(kQS[static_cast<std::size_t>(iq - 1)]);
    loaded_ = true;
    ih_ = ih;
    io_ = io;
    return true;
}

// Fortran: fdss.f:fDSS l.58-69 and 140-200
FdssValues FdssGrid::eval(int ic, double x, double q2, std::ostream* log) const {
    if (!loaded_) throw std::logic_error("FdssGrid::eval before load");
    // CHECK OF X AND Q2 VALUES (X.LT.0.05D0 is a binary128 comparison, __lttf2)
    if ((static_cast<quad>(x) < static_cast<quad>(5) / static_cast<quad>(100)) || (x > 1.0)) {
        if (log) {
            char b[64];
            std::snprintf(b, sizeof b, "%3.1f", x);   // FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE ', f3.1)
            std::string s(b);
            if (s.size() > 3) s = "***";
            *log << "  PARTON INTERPOLATION: X OUT OF RANGE " << s << "\n";
        }
    }
    if ((q2 < 1.) || (q2 > 1.e5)) {
        if (log) *log << "  PARTON INTERPOLATION: Q2 OUT OF RANGE\n";
    }
    // INTERPOLATION (l.141-158)
    double xt[2];
    xt[0] = std::log(x);
    xt[1] = std::log(q2);
    const quad omx = static_cast<quad>(1) - static_cast<quad>(x);   // (1.D0-X), binary128
    const quad p4 = powR16I4(omx, 4), p7 = powR16I4(omx, 7);
    const double x05 = std::pow(x, 0.5), x03 = std::pow(x, 0.3);
    auto term = [&](int k, quad pw, double xp) {
        const double fi = ffint(NARG, xt, NA.data(), ARRF.data(), F[static_cast<std::size_t>(k)].data());
        return static_cast<double>((static_cast<quad>(fi) * pw) * static_cast<quad>(xp));
    };
    const double UTOT = term(0, p4, x05);
    const double DTOT = term(1, p4, x05);
    const double STOT = term(2, p4, x05);
    const double CTOT = term(3, p7, x03);
    const double BTOT = term(4, p7, x03);
    const double GL = term(5, p4, x03);
    const double UVAL = term(6, p4, x05);
    const double DVAL = term(7, p4, x05);
    const double SVAL = term(8, p4, x05);

    const double Up = (UTOT + UVAL) / 2.;
    const double UBp = (UTOT - UVAL) / 2.;
    const double Dp = (DTOT + DVAL) / 2.;
    const double DBp = (DTOT - DVAL) / 2.;
    const double Sp = (STOT + SVAL) / 2.;
    const double SBp = (STOT - SVAL) / 2.;
    const double Cp = CTOT / 2.;
    const double Bp = BTOT / 2.;

    FdssValues r;
    r.GL = GL;
    if (ic == 1) {
        r.U = Up; r.UB = UBp; r.D = Dp; r.DB = DBp; r.S = Sp; r.SB = SBp; r.C = Cp; r.B = Bp;
    } else if (ic == -1) {
        r.U = UBp; r.UB = Up; r.D = DBp; r.DB = Dp; r.S = SBp; r.SB = Sp; r.C = Cp; r.B = Bp;
    } else if (ic == 0) {
        r.U = (UBp + Up) / 2.;
        r.UB = r.U;
        r.D = (DBp + Dp) / 2.;
        r.DB = r.D;
        r.S = (SBp + Sp) / 2.;
        r.SB = r.S;
        r.C = Cp;
        r.B = Bp;
    } else {
        throw FortranStop("fdss.f:194-196 STOP: WRONG CHARGE (IC = " + std::to_string(ic) + ")");
    }
    return r;
}

// Fortran: fdss.f:fFINT (l.206-247), CERNLIB E104. FAC*(1.D0-D(I)) is a binary128 product (REAL(16) literal).
double ffint(int narg, const double* arg, const int* nent, const double* ent, const double* table) {
    if (narg < 1 || narg > 2) throw std::invalid_argument("ffint: NARG must be 1 or 2 (fdss.f:209 D(2), NCOMB(2))");
    double d[3] = {0, 0, 0};   // D(1:2)
    int ncomb[3] = {0, 0, 0}, ient[3] = {0, 0, 0};
    auto ARG = [&](int i) { return arg[i - 1]; };
    auto NENT = [&](int i) { return nent[i - 1]; };
    auto ENT = [&](int j) { return ent[j - 1]; };
    int kd = 1, m = 1, ja = 1;
    for (int i = 1; i <= narg; ++i) {
        ncomb[i] = 1;
        const int jb = ja - 1 + NENT(i);
        int j = ja;
        bool found = false;
        for (j = ja; j <= jb; ++j)
            if (ARG(i) <= ENT(j)) {
                found = true;
                break;
            }
        if (!found) j = jb;
        if (j == ja) j = j + 1;   // 3
        const int jr = j - 1;     // 4
        d[i] = (ENT(j) - ARG(i)) / (ENT(j) - ENT(jr));
        ient[i] = j - ja;
        kd = kd + ient[i] * m;
        m = m * NENT(i);
        ja = jb + 1;   // 5
    }
    double r = 0.;
    for (;;) {
        // 10
        double fac = 1.;
        int iadr = kd, ifadr = 1;
        for (int i = 1; i <= narg; ++i) {
            if (ncomb[i] == 0) {
                fac = fac * d[i];   // 12
                iadr = iadr - ifadr;
            } else {
                fac = static_cast<double>(static_cast<quad>(fac) * (static_cast<quad>(1) - static_cast<quad>(d[i])));
            }
            ifadr = ifadr * NENT(i);   // 15
        }
        r = r + fac * table[iadr - 1];
        int il = narg;
        // 40 / 80
        bool again = false;
        for (;;) {
            if (ncomb[il] != 0) {
                ncomb[il] = 0;
                if (il != narg) {
                    il = il + 1;
                    for (int k = il; k <= narg; ++k) ncomb[k] = 1;
                }
                again = true;   // GO TO 10
                break;
            }
            il = il - 1;   // 80
            if (il == 0) break;
        }
        if (!again) break;
    }
    return r;
}

}  // namespace sidis
}  // namespace simc
