// src/radiation/Brem.cpp — see include/simc/radiation/kernel/Brem.h.
// The internal radiative correction of SIMC (brem.f): brem (electron-only, Mo-Tsai type, with the soft-photon
// exponentiation), bremos (including the radiation of the outgoing hadron), and their helpers inter/inter_prime,
// srad and spence. Called per event by radc_init_ev (RadcInit.cpp) to get the soft and hard correction factors.
// Fortran: simc_gfortran/brem.f. Operation order follows the Fortran source (left-to-right evaluation,
// parentheses kept); comments give the Fortran line numbers.
#include "simc/radiation/kernel/Brem.h"

#include <cmath>
#include <ostream>

#include "simc/material/FortranFormat.h"

namespace simc {
namespace radiation {

namespace {

// brem.f local PARAMETERs (brem l.11-15, bremos l.359-364, inter/inter_prime/spence l.224-225, srad l.257-259)
constexpr double pi = 3.141592653589793;
constexpr double twopi = 2. * pi;
constexpr double am = .93827231;     // brem: proton mass [GeV]
constexpr double mp = .93827231;     // bremos
constexpr double ame = .00051099906;
constexpr double e2 = 1. / 137.0359895;

// PARAMETER**integer is folded by gfortran with MPFR (correctly rounded). ame**2 and am**2 are single
// products (x*x is correctly rounded as well); ame**4 and am**4 need the correctly rounded 4th power, which the
// 64-bit-mantissa long double product gives here (checked bit for bit by the brem/bremos L1 tables).
constexpr double kAme2 = ame * ame;
constexpr double kAm2 = am * am;
const double kAme4 = static_cast<double>(static_cast<long double>(ame) * ame * ame * ame);
const double kAm4 = static_cast<double>(static_cast<long double>(am) * am * am * am);

inline double sq(double x) { return x * x; }

// Fortran list-directed write(6,*) of real*8 values (used only when produce_output).
void wr(std::ostream* log, const char* label, double v) {
    if (log) *log << label << fortran::listReal(v) << "\n";
}
void wr2(std::ostream* log, double a, double b) {
    if (log) *log << fortran::listReal(a) << fortran::listReal(b) << "\n";
}

}  // namespace

// Fortran: brem.f:brem (l.6-214)
double brem(const BremFlags& f, double ein, double eout, double egamma, bool radiate_proton, double& bsoft,
            double& bhard, double& dbsoft, std::ostream* log) {
    double eang, ak, akp, ap, pang, eta, ape, q2, de;
    double aprod, adot, ar1, ar2, alpha;
    double bpi = 0, bpf = 0, bpp = 0, bei, bef, bee, bepii = 0, bepif = 0, bepfi = 0, bepff = 0;
    double dbpi = 0, dbpf = 0, dbpp = 0, dbei, dbef, dbee, dbepii = 0, dbepif = 0, dbepfi = 0, dbepff = 0;
    double b, bz, bzz, db, dbz, dbzz;

    // Convert to GeV (l.32-34)
    ak = ein / 1000.;
    akp = eout / 1000.;
    de = egamma / 1000.;

    // electron scattering angle (l.38-40)
    eang = 2. * std::asin(std::pow(am / (2. * ak) * (ak / akp - 1.), 0.5));
    eta = 1. + 2. * ak * sq(std::sin(eang / 2.)) / am;
    q2 = 4. * ak * akp * sq(std::sin(eang / 2.));

    // final proton energy/momentum/scattering angle (l.44-46)
    ape = am + ak - akp;
    ap = std::sqrt(sq(ape) - kAm2);
    pang = std::acos((ak - akp * std::cos(eang)) / ap);

    if (f.produce_output) {   // l.48-56
        wr(log, " q2   ", q2);
        wr(log, " k\t   ", ak);
        wr(log, " kp   ", akp);
        wr(log, " eang ", eang * 180. / pi);
        wr(log, " ap   ", ap);
        wr(log, " pang ", pang * 180. / pi);
        wr(log, " eta  ", eta);
    }

    // direct initial electron (l.64-66)
    aprod = 1.e0;
    bei = aprod * (-1. / (2. * pi)) * std::log(ak / de);
    dbei = aprod * (1. / (2. * pi * de));

    // direct final electron (l.70-72)
    aprod = 1.e0;
    bef = aprod * (-1. / (2. * pi)) * std::log(akp / de);
    dbef = aprod * (1. / (2. * pi * de));

    // e-e interference (l.76-84)
    aprod = -1.e0;
    adot = ak * akp * (1. - std::cos(eang));
    alpha = 2. * kAme2 - 2. * adot;
    ar1 = 0.5 + std::sqrt(sq(adot) - kAme4) / alpha;
    ar2 = 0.5 - std::sqrt(sq(adot) - kAme4) / alpha;
    bee = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, ak, akp, de);
    dbee = aprod * adot / (pi * alpha * (ar1 - ar2) * de) * (std::log((ar1 - 1) / ar1) - std::log((ar2 - 1) / ar2));
    if (f.produce_output) wr2(log, ar1, ar2);

    if (radiate_proton) {   // l.88-162
        // initial p direct
        aprod = 1.e0;
        bpi = aprod * (-1. / (2. * pi)) * std::log(am / de);
        dbpi = aprod * (1. / (2. * pi * de));

        // final p direct
        aprod = 1.e0;
        bpf = aprod * (-1. / (2. * pi)) * std::log(ape / de);
        dbpf = aprod * (1 / (2. * pi * de));

        // p-p interference (l.104-112)
        aprod = -1.e0;
        adot = am * ape;
        alpha = 2. * kAm2 - 2. * adot;
        ar1 = 0.5 + std::sqrt(sq(adot) - kAm4) / alpha;
        ar2 = 0.5 - std::sqrt(sq(adot) - kAm4) / alpha;
        bpp = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, am, ape, de);
        dbpp = aprod * adot / (pi * alpha * (ar1 - ar2) * de) * (std::log((ar1 - 1) / ar1) - std::log((ar2 - 1) / ar2));
        if (f.produce_output) wr2(log, ar1, ar2);

        constexpr double kAmPlusAme = kAm2 + kAme2;          // am**2+ame**2 (folded)
        constexpr double kAmeAm2 = (ame * am) * (ame * am);  // (ame*am)**2 (folded)

        // ei-pi interference (l.116-124)
        aprod = -1.e0;
        adot = ak * am;
        alpha = kAmPlusAme - 2. * adot;
        ar1 = (kAm2 - adot + std::sqrt(sq(adot) - kAmeAm2)) / alpha;
        ar2 = (kAm2 - adot - std::sqrt(sq(adot) - kAmeAm2)) / alpha;
        bepii = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, ak, am, de);
        dbepii = aprod * adot / (pi * alpha * (ar1 - ar2) * de) * (std::log((ar1 - 1) / ar1) - std::log((ar2 - 1) / ar2));
        if (f.produce_output) wr2(log, ar1, ar2);

        // ef-pf interference (l.128-136)
        aprod = -1.e0;
        adot = akp * ape - akp * ap * std::cos(eang + pang);
        alpha = kAmPlusAme - 2. * adot;
        ar1 = (kAm2 - adot + std::sqrt(sq(adot) - kAmeAm2)) / alpha;
        ar2 = (kAm2 - adot - std::sqrt(sq(adot) - kAmeAm2)) / alpha;
        bepff = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, akp, ape, de);
        dbepff = aprod * adot / (pi * alpha * (ar1 - ar2) * de) * (std::log((ar1 - 1) / ar1) - std::log((ar2 - 1) / ar2));
        if (f.produce_output) wr2(log, ar1, ar2);

        // ei-pf interference (l.140-148)
        aprod = 1.e0;
        adot = ak * ape - ak * ap * std::cos(pang);
        alpha = kAmPlusAme - 2. * adot;
        ar1 = (kAm2 - adot + std::sqrt(sq(adot) - kAmeAm2)) / alpha;
        ar2 = (kAm2 - adot - std::sqrt(sq(adot) - kAmeAm2)) / alpha;
        bepif = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, ak, ape, de);
        dbepif = aprod * adot / (pi * alpha * (ar1 - ar2) * de) * (std::log((ar1 - 1) / ar1) - std::log((ar2 - 1) / ar2));
        if (f.produce_output) wr2(log, ar1, ar2);

        // ef-pi interference (l.152-160)
        aprod = 1.e0;
        adot = akp * am;
        alpha = kAmPlusAme - 2. * adot;
        ar1 = (kAm2 - adot + std::sqrt(sq(adot) - kAmeAm2)) / alpha;
        ar2 = (kAm2 - adot - std::sqrt(sq(adot) - kAmeAm2)) / alpha;
        bepfi = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, akp, am, de);
        dbepfi = aprod * adot / (pi * alpha * (ar1 - ar2) * de) * (std::log((ar1 - 1) / ar1) - std::log((ar2 - 1) / ar2));
        if (f.produce_output) wr2(log, ar1, ar2);
    }

    // All together now! (l.165-183)
    b = 2. * e2 * (bei + bef + bee);
    if (radiate_proton) {
        bzz = 2. * e2 * (bpi + bpf + bpp);
        bz = 2. * e2 * (bepii + bepff + bepif + bepfi);
    } else {
        bzz = 0.0;
        bz = 0.0;
    }
    bsoft = b + bz + bzz;
    bhard = -1. * (e2 / pi) * (-28 / 9. + 13. / 6. * std::log(q2 / kAme2));
    db = 2. * e2 * (dbei + dbef + dbee);
    if (radiate_proton) {
        dbzz = 2. * e2 * (dbpi + dbpf + dbpp);
        dbz = 2. * e2 * (dbepii + dbepff + dbepif + dbepfi);
    } else {
        dbzz = 0.0;
        dbz = 0.0;
    }
    dbsoft = db + dbz + dbzz;

    if (f.produce_output) {   // l.185-200
        if (log) *log << "   ----- results ----- \n";
        wr(log, " b=     ", b);
        wr(log, " bz=    ", bz);
        wr(log, " bzz=   ", bzz);
        wr(log, " bhard= ", bhard);
        wr(log, " total= ", bsoft + bhard);
        wr(log, " exp=   ", 1. - std::exp(-1. * (bsoft)) * (1. - bhard));
        if (log) *log << " \n ultra-relativistic limit\n";
        srad(ak, akp, eang, q2, ame, am, ap, de, f.produce_output, log);
        if (log) *log << "  Schwinger Result\n";
        const double bsch = 2. * e2 / pi *
                            ((std::log(ak / de) - 13. / 12.) * (std::log(q2 / kAme2) - 1.) + 17. / 36. +
                             0.5 * (sq(pi) / 6. - spence(sq(std::cos(eang / 2.)))));
        wr(log, "  b=     ", bsch);
    }

    // the result (l.205-211)
    dbsoft = dbsoft / 1000.;
    double r;
    if (f.exponentiate)
        r = -dbsoft / std::exp(bsoft);
    else
        r = 1. - dbsoft;
    if (f.include_hard) r = r * (1. - bhard);
    return r;
}

// Fortran: brem.f:inter (l.218-247)
double inter(bool calculate_spence, double alpha, double ar1, double ar2, double e1, double e2_, double de) {
    const double de2 = e1 - e2_;
    const double amult = -1. / (alpha * (ar1 - ar2));
    double r = std::log(std::fabs((e2_ / de) + ar1 * (de2 / de))) * std::log(std::fabs((ar1 - 1.) / ar1)) -
               std::log(std::fabs((e2_ / de) + ar2 * (de2 / de))) * std::log(std::fabs((ar2 - 1.) / ar2));
    if (calculate_spence) {
        const double arg1 = (de2 / (e2_ + ar1 * de2)) * (ar1 - 1.);
        const double arg2 = (de2 / (e2_ + ar1 * de2)) * (ar1);
        const double arg3 = (de2 / (e2_ + ar2 * de2)) * (ar2 - 1.);
        const double arg4 = (de2 / (e2_ + ar2 * de2)) * (ar2);
        r = r - spence(arg1) + spence(arg2) + spence(arg3) - spence(arg4);
    }
    r = r * amult / (pi);
    return r;
}

// Fortran: brem.f:srad (l.251-283)
SradValues srad(double ak, double akp, double eang, double q2, double ame_, double am_, double ap, double de,
                bool produce_output, std::ostream* log) {
    constexpr double alpha = 1. / 137.0359895;
    SradValues v;
    const double ep = std::sqrt(sq(ap) + sq(am_));
    const double dhard = -13. / 12. * (std::log(q2 / sq(ame_)) - 1.) + 17. / 36.;
    const double eta = 1. + 2. * ak * (sq(std::sin(eang / 2.))) / am_;
    (void)dhard;   // computed and never used (as in the Fortran)
    v.dsoft = alpha / pi * std::log(ak * akp / sq(de)) * (std::log(q2 / sq(ame_)) - 1.);
    v.dsoftz = alpha / pi * std::log(eta) * (std::log(ak * akp / sq(de)) + std::log(am_ * ep / sq(de)));
    v.dsoftzz = alpha / pi * std::log(am_ * ep / sq(de)) * (std::log(q2 / sq(am_)) - 1.);
    if (v.dsoftzz <= 0) v.dsoftzz = 0.0;
    if (produce_output && log) {
        *log << "  \n";
        wr(log, " b=    ", v.dsoft);
        wr(log, " bz=   ", v.dsoftz);
        wr(log, " bzz=  ", v.dsoftzz);
    }
    return v;
}

// Fortran: brem.f:spence (l.287-320) — the approximation actually compiled (the series is commented out)
double spence(double ax) {
    const double bx = std::fabs(ax);
    if (bx <= 1) return 0.0;
    return -0.5 * sq(std::log(bx));
}

// Fortran: brem.f:bremos (l.344-574)
double bremos(const BremFlags& f, double egamma, double k_ix, double k_iy, double k_iz, double k_fx, double k_fy,
              double k_fz, double p_ix, double p_iy, double p_iz, double p_fx, double p_fy, double p_fz, double p_fe,
              bool radiate_proton, double& bsoft, double& bhard, double& dbsoft, std::ostream* log) {
    struct FourVector {
        double e, x, y, z;
    } k_i{}, k_f{}, p_i{}, p_f{};
    double ami, amf, q2, aprod, adot, ar1, ar2, alpha;
    double bpi = 0, bpf = 0, bpp = 0, bei, bef, bee, bepii = 0, bepif = 0, bepfi = 0, bepff = 0;
    double dbpi = 0, dbpf = 0, dbpp = 0, dbei, dbef, dbee, dbepii = 0, dbepif = 0, dbepfi = 0, dbepff = 0;
    double b, bz, bzz, db, dbz, dbzz;

    // l.389-402: to GeV
    const double de = egamma / 1000.;
    k_i.x = k_ix / 1000.;
    k_i.y = k_iy / 1000.;
    k_i.z = k_iz / 1000.;
    k_f.x = k_fx / 1000.;
    k_f.y = k_fy / 1000.;
    k_f.z = k_fz / 1000.;
    p_i.x = p_ix / 1000.;
    p_i.y = p_iy / 1000.;
    p_i.z = p_iz / 1000.;
    p_f.e = p_fe / 1000.;
    p_f.x = p_fx / 1000.;
    p_f.y = p_fy / 1000.;
    p_f.z = p_fz / 1000.;

    // l.405-409: electron energies (x**0.5 = libm pow, Q9); p_i%e forced to mp
    k_i.e = std::pow(sq(k_i.x) + sq(k_i.y) + sq(k_i.z) + kAme2, 0.5);
    k_f.e = std::pow(sq(k_f.x) + sq(k_f.y) + sq(k_f.z) + kAme2, 0.5);
    p_i.e = mp;

    // l.424-429
    q2 = -1. * (sq(k_f.e - k_i.e) - sq(k_f.x - k_i.x) - sq(k_f.y - k_i.y) - sq(k_f.z - k_i.z));
    if (f.produce_output) wr(log, " q2= ", q2);
    ami = mp;
    amf = std::pow(sq(p_f.e) - sq(p_f.x) - sq(p_f.y) - sq(p_f.z), 0.5);

    // direct initial electron (l.436-438)
    aprod = 1.e0;
    bei = aprod * (-1. / twopi) * std::log(k_i.e / de);
    dbei = aprod * (-1. / twopi) * (-1. / de);

    // direct final electron (l.441-443)
    aprod = 1.e0;
    bef = aprod * (-1. / twopi) * std::log(k_f.e / de);
    dbef = aprod * (-1. / twopi) * (-1. / de);

    // e-e interference (l.446-453)
    aprod = -1.e0;
    adot = k_i.e * k_f.e - k_i.x * k_f.x - k_i.y * k_f.y - k_i.z * k_f.z;
    alpha = 2. * kAme2 - 2. * adot;
    ar1 = 0.5 + std::sqrt(4. * sq(adot) - 4. * kAme4) / (2. * alpha);
    ar2 = 0.5 - std::sqrt(4. * sq(adot) - 4. * kAme4) / (2. * alpha);
    bee = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, k_i.e, k_f.e, de);
    dbee = aprod * adot * interPrime(alpha, ar1, ar2, de);
    if (f.produce_output) wr2(log, ar1, ar2);

    if (radiate_proton) {   // l.457-519
        // initial p direct
        aprod = 1.e0;
        bpi = aprod * (-1. / twopi) * std::log(p_i.e / de);
        dbpi = aprod * (-1. / twopi) * (-1. / de);

        // final p direct
        aprod = 1.e0;
        bpf = aprod * (-1. / twopi) * std::log(p_f.e / de);
        dbpf = aprod * (-1. / twopi) * (-1. / de);

        // p-p interference (l.470-477)
        aprod = -1.e0;
        adot = p_i.e * p_f.e - p_i.x * p_f.x - p_i.y * p_f.y - p_i.z * p_f.z;
        alpha = sq(ami) + sq(amf) - 2. * adot;
        ar1 = (2. * sq(amf) - 2. * adot + std::sqrt(4. * sq(adot) - 4. * sq(ami * amf))) / (2. * alpha);
        ar2 = (2. * sq(amf) - 2. * adot - std::sqrt(4. * sq(adot) - 4. * sq(ami * amf))) / (2. * alpha);
        bpp = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, p_i.e, p_f.e, de);
        dbpp = aprod * adot * interPrime(alpha, ar1, ar2, de);
        if (f.produce_output) wr2(log, ar1, ar2);

        // ei-pi interference (l.480-487)
        aprod = -1.e0;
        adot = k_i.e * p_i.e - k_i.x * p_i.x - k_i.y * p_i.y - k_i.z * p_i.z;
        alpha = sq(ami) + kAme2 - 2. * adot;
        ar1 = (2. * sq(ami) - 2. * adot + std::sqrt(4. * sq(adot) - 4. * sq(ame * ami))) / (2. * alpha);
        ar2 = (2. * sq(ami) - 2. * adot - std::sqrt(4. * sq(adot) - 4. * sq(ame * ami))) / (2. * alpha);
        bepii = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, k_i.e, p_i.e, de);
        dbepii = aprod * adot * interPrime(alpha, ar1, ar2, de);
        if (f.produce_output) wr2(log, ar1, ar2);

        // ef-pf interference (l.490-497)
        aprod = -1.e0;
        adot = k_f.e * p_f.e - k_f.x * p_f.x - k_f.y * p_f.y - k_f.z * p_f.z;
        alpha = sq(amf) + kAme2 - 2. * adot;
        ar1 = (2. * sq(amf) - 2. * adot + std::sqrt(4. * sq(adot) - 4. * sq(ame * amf))) / (2. * alpha);
        ar2 = (2. * sq(amf) - 2. * adot - std::sqrt(4. * sq(adot) - 4. * sq(ame * amf))) / (2. * alpha);
        bepff = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, k_f.e, p_f.e, de);
        dbepff = aprod * adot * interPrime(alpha, ar1, ar2, de);
        if (f.produce_output) wr2(log, ar1, ar2);

        // ei-pf interference (l.500-507)
        aprod = 1.e0;
        adot = k_i.e * p_f.e - k_i.x * p_f.x - k_i.y * p_f.y - k_i.z * p_f.z;
        alpha = sq(amf) + kAme2 - 2. * adot;
        ar1 = (2. * sq(amf) - 2. * adot + std::sqrt(4. * sq(adot) - 4. * sq(ame * amf))) / (2. * alpha);
        ar2 = (2. * sq(amf) - 2. * adot - std::sqrt(4. * sq(adot) - 4. * sq(ame * amf))) / (2. * alpha);
        bepif = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, k_i.e, p_f.e, de);
        dbepif = aprod * adot * interPrime(alpha, ar1, ar2, de);
        if (f.produce_output) wr2(log, ar1, ar2);

        // ef-pi interference (l.510-517)
        aprod = 1.e0;
        adot = k_f.e * p_i.e - k_f.x * p_i.x - k_f.y * p_i.y - k_f.z * p_i.z;
        alpha = sq(ami) + kAme2 - 2. * adot;
        ar1 = (2. * sq(ami) - 2. * adot + std::sqrt(4. * sq(adot) - 4. * sq(ame * ami))) / (2. * alpha);
        ar2 = (2. * sq(ami) - 2. * adot - std::sqrt(4. * sq(adot) - 4. * sq(ame * ami))) / (2. * alpha);
        bepfi = aprod * adot * inter(f.calculate_spence, alpha, ar1, ar2, k_f.e, p_i.e, de);
        dbepfi = aprod * adot * interPrime(alpha, ar1, ar2, de);
        if (f.produce_output) wr2(log, ar1, ar2);
    }

    // All together now! (l.522-542)
    b = 2. * e2 * (bei + bef + bee);
    if (radiate_proton) {
        bzz = 2. * e2 * (bpi + bpf + bpp);
        bz = 2. * e2 * (bepii + bepff + bepif + bepfi);
    } else {
        bzz = 0.0;
        bz = 0.0;
    }
    bsoft = b + bz + bzz;
    bhard = -1. * (e2 / pi) * (-28 / 9. + 13. / 6. * std::log(q2 / kAme2));
    db = 2. * e2 * (dbei + dbef + dbee);
    if (radiate_proton) {
        dbzz = 2. * e2 * (dbpi + dbpf + dbpp);
        dbz = 2. * e2 * (dbepii + dbepff + dbepif + dbepfi);
    } else {
        dbzz = 0.0;
        dbz = 0.0;
    }
    dbsoft = db + dbz + dbzz;

    if (f.produce_output) {   // l.544-560
        if (log) *log << "   ----- results ----- \n";
        wr(log, " b+bhard= ", b + bhard);
        wr(log, " bz=      ", bz);
        wr(log, " bzz=     ", bzz);
        wr(log, " bhard=   ", bhard);
        wr(log, " total=   ", 1 - bsoft - bhard);
        wr(log, " exp=     ", std::exp(-1. * (bsoft)) * (1. - bhard));
        wr(log, " 1-bhard= ", 1 - bhard);
        wr(log, " exps=    ", std::exp(-1. * (bsoft)));
        wr(log, " expse=   ", std::exp(-1. * b));
        wr(log, " expsp=   ", std::exp(-1. * (bz + bzz)));
        if (log) *log << " \n  Schwinger Result\n";
        const double bsch =
            2. * e2 / pi * ((std::log(k_i.e / de) - 13. / 12.) * (std::log(q2 / kAme2) - 1.) + 17. / 36.);
        wr(log, "  b=     ", bsch);
    }

    // l.565-571
    dbsoft = dbsoft / 1000.;
    double r;
    if (f.exponentiate)
        r = -dbsoft / std::exp(bsoft);
    else
        r = 1. - dbsoft;
    if (f.include_hard) r = r * (1. - bhard);
    return r;
}

// Fortran: brem.f:inter_prime (l.581-596)
double interPrime(double alpha, double ar1, double ar2, double de) {
    const double amult = -1. / (alpha * (ar1 - ar2));
    return (-1. / de) * amult / pi * (std::log(std::fabs((ar1 - 1.) / ar1)) - std::log(std::fabs((ar2 - 1.) / ar2)));
}

}  // namespace radiation
}  // namespace simc
