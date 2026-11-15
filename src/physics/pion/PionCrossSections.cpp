// src/physics/pion/PionCrossSections.cpp
//
// physics_pion.f: sig_blok, sig_param04, sig_param_3000, sig_param_2021, exclfit, sigmaid (+ the MAID tables), peepi:
// the pion-electroproduction models of SIMC (which_pion selects the parametrisation) and the transformation of the
// centre-of-mass cross section to the lab with its Jacobian. See include/simc/physics/pion/PionCrossSections.h. sig_multipole is in SigMultipole.cpp.
//
// Numerics (D-010, D-024): all literals are double (-fdefault-real-8); Fortran x**2 -> x*x; real exponents -> pow;
// products and quotients are evaluated left to right as written in the Fortran (a*b/c*d = ((a*b)/c)*d); the
// exponentiation binds tighter than the product (q2*fpi**2 = q2*(fpi*fpi)).
#include "simc/physics/pion/PionCrossSections.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "simc/core/SimcConstants.h"
#include "simc/reactions/EventChain.h"

namespace simc {
namespace pion {

namespace {
inline double sq(double x) { return x * x; }   // Fortran x**2
using constants::pi;

/// The pi- reduction of sigT and sigTT shared by sig_blok, sig_param04, sig_param_3000 (which_pion 1, 11, 3).
inline bool isPiMinus(int which_pion) { return which_pion == 1 || which_pion == 11 || which_pion == 3; }
}  // namespace

// Fortran: physics_pion.f:sig_blok (l.404-465)
double sigBlok(double thetacm, double phicm, double t, double q2_gev, double s_gev, double eps, double mtar_gev,
               int which_pion) {
    // exp(-11.5*abs(t)) is compiled as exp(-(abs(t)*11.5)) (mulsd, then xorpd): same value, and for t = NaN the
    // argument is -NaN (the C++ -11.5*fabs(t) would be +NaN), which glibc exp returns unchanged.
    double sigl = 27.8 * std::exp(-(std::fabs(t) * 11.5));
    double sigt = 10.0 * (5. * std::fabs(t)) * std::exp(-(std::fabs(t) * 5.));
    const double siglt = 0.0 * std::sin(thetacm);
    // -(a)*b as the object code computes it: the product first, then the sign flip (xorpd after mulsd, physics_pion.o);
    // identical for finite values, and gives the Fortran sign of a NaN (D-048 F3 garbage arguments).
    double sigtt = -((4.0 * sigl + 0.5 * sigt) * sq(std::sin(thetacm)));

    if (isPiMinus(which_pion)) {
        sigt = sigt * 0.25 * (1. + 3. * std::exp(-(std::fabs(t) * 10.)));
        sigtt = sigtt * 0.25 * (1. + 3. * std::exp(-(std::fabs(t) * 10.)));
    }
    // mrho_emp = 0.712 is set but not used (l.437)
    const double fpi = 1. / (1. + 1.65 * q2_gev + 0.5 * sq(q2_gev));
    const double fpi2 = sq(fpi);

    // Operand order of the object code (physics_pion.o, sig_blok__): identical values for finite numbers (commuted
    // products/sums), and the Fortran sign of a NaN (x86 returns the first operand's NaN when both are NaN; D-048 F3
    // garbage arguments give NaN for ~0.05 % of the SOS Delta events).
    sigl = ((fpi2 * q2_gev) * sigl) / 0.1215;
    sigt = sigt / (0.3 + q2_gev);
    sigtt = sigtt / (0.3 + q2_gev);

    const double sig219 =
        (((sigl * eps + sigt) + (eps * std::cos(2. * phicm)) * sigtt) + (std::sqrt(2.0 * eps * (1. + eps)) * std::cos(phicm)) * siglt) /
        1.e0;
    double sig = sig219 * 15.333 / sq(s_gev - sq(mtar_gev));
    sig = sig / 2. / pi / 1.e+06;
    return sig;
}

// Fortran: physics_pion.f:sig_param04 (l.467-516)
double sigParam04(double thetacm, double phicm, double t, double q2_gev, double s_gev, double eps, double mtar_gev,
                  int which_pion) {
    const double fpi = 1. / (1. + 1.77 * q2_gev + 0.05 * sq(q2_gev));
    // q2fpi = q2_gev*fpi is set but not used (l.484)
    const double q2fpi2 = q2_gev * sq(fpi);
    // FORTRAN-QUIRK(D-038 F13): pole at |t| = 0.02 GeV^2 (Inf/NaN rows of the L1 table)
    const double polefactor = std::fabs(t) / sq(std::fabs(t) - 0.02);

    // alog(q2_gev): the specific name of the default-real LOG; with -fdefault-real-8 it is the double log.
    const double sigl = 350.0 * q2fpi2 * std::exp((16.0 - 7.5 * std::log(q2_gev)) * (-t));
    double sigt = 4.5 / q2_gev + 2.0 / sq(q2_gev);
    const double siglt = (std::exp(0.79 - 3.4 / std::sqrt(q2_gev) * (t)) + 1.1 - 3.6 / sq(q2_gev)) * std::sin(thetacm);
    double sigtt = -5.0 / sq(q2_gev) * polefactor * sq(std::sin(thetacm));

    if (isPiMinus(which_pion)) {
        sigt = sigt * 0.25 * (1. + 3. * std::exp(-10. * std::fabs(t)));
        sigtt = sigtt * 0.25 * (1. + 3. * std::exp(-10. * std::fabs(t)));
    }

    const double sig219 =
        (sigt + eps * sigl + eps * std::cos(2. * phicm) * sigtt + std::sqrt(2.0 * eps * (1. + eps)) * std::cos(phicm) * siglt) /
        1.e0;
    double sig = sig219 * 8.539 / sq(s_gev - sq(mtar_gev));
    sig = sig / 2. / pi / 1.e+06;
    return sig;
}

// Fortran: physics_pion.f:sig_param_3000 (l.519-573)
double sigParam3000(double thetacm, double phicm, double t, double q2_gev, double s_gev, double eps, double mtar_gev,
                    int which_pion) {
    const double fpi = 1.0 / (1.0 + 1.77 * q2_gev + 0.12 * sq(q2_gev));
    const double q2fpi2 = q2_gev * sq(fpi);

    const double sigl = 19.8 * std::fabs(t) / sq(std::fabs(t) + 0.02) * q2fpi2 * std::exp(-3.66 * std::fabs(t));
    double sigt = 5.96 / q2_gev * std::exp(-0.013 * sq(q2_gev));
    const double siglt = (16.533 / (1.0 + q2_gev)) * std::exp(-5.1437 * std::fabs(t)) * std::sin(thetacm);
    double sigtt = (-178.06 / (1.0 + q2_gev)) * std::exp(-7.1381 * std::fabs(t)) * sq(std::sin(thetacm));

    if (isPiMinus(which_pion)) {
        sigt = sigt * 0.25 * (1.0 + 3.0 * std::exp(-10.0 * std::fabs(t)));
        sigtt = sigtt * 0.25 * (1.0 + 3.0 * std::exp(-10.0 * std::fabs(t)));
    }

    const double sig219 = (sigt + eps * sigl + eps * std::cos(2.0 * phicm) * sigtt +
                           std::sqrt(2.0 * eps * (1.0 + eps)) * std::cos(phicm) * siglt) /
                          1.0;
    double sig = sig219 * 8.539 / sq(s_gev - sq(mtar_gev));
    sig = sig / 2.0 / pi / 1.0e+06;
    return sig;
}

// Fortran: physics_pion.f:exclfit (l.807-854)
ExclfitResult exclfit(double t, double thetacm, double phicm, double q2_gev, double s_gev, double eps,
                      const std::array<double, 17>& pv, double fpifact) {
    // p(k) of the Fortran = p(k)
    auto p = [&pv](int k) { return pv[static_cast<std::size_t>(k - 1)]; };
    // FORTRAN-DIFF-CANDIDATE(D-038 F4): the DATA value 0.938 (l.812), also for coherent production (Mtar_struck = M_A).
    // Question: is the nucleon fit meant for coherent production? Kept in both modes.
    const double mtar_gev = 0.938;
    ExclfitResult r;

    const double fpi = fpifact / (1.0 + p(1) * q2_gev + p(2) * sq(q2_gev));
    const double q2fpi2 = q2_gev * sq(fpi);

    r.sigl = (p(3) + p(15) / q2_gev) * std::fabs(t) / sq(std::fabs(t) + 0.02) * q2fpi2 * std::exp(p(4) * std::fabs(t));
    r.sigl = r.sigl / (std::pow(s_gev, p(11)) + std::pow(std::sqrt(s_gev), p(17)));

    r.sigt = p(5) / q2_gev * std::exp(p(6) * sq(q2_gev));
    r.sigt = r.sigt / (std::pow(s_gev, p(12)) + std::pow(std::sqrt(s_gev), p(16)));
    r.sigt = r.sigt * std::exp(p(14) * std::fabs(t));

    r.siglt = (p(7) / (1.0 + p(10) * q2_gev)) * std::exp(p(8) * std::fabs(t)) * std::sin(thetacm);
    r.siglt = r.siglt / std::pow(s_gev, p(13));

    r.sigtt = (p(9) / (1. + 1.0 * q2_gev)) * std::exp(-7.0 * std::fabs(t)) * sq(std::sin(thetacm));

    const double sig219 = (r.sigt + eps * r.sigl + eps * std::cos(2.0 * phicm) * r.sigtt +
                           std::sqrt(2.0 * eps * (1.0 + eps)) * std::cos(phicm) * r.siglt) /
                          1.0;

    r.sig = sig219 * 8.539 / sq(s_gev - sq(mtar_gev));
    r.sig = r.sig / 2.0 / 3.1415928 / 1.0e+06;   // the Fortran literal 3.1415928 (not constants.inc pi)
    return r;
}

// Fortran: physics_pion.f:746-786 (DATA pp, pm)
const std::array<double, 17> kParam2021PiPlus = {1.60077,  -0.01523, 37.08142, -4.11060, 23.26192, 0.00983,
                                                 0.87073,  -5.77115, -271.08678, 0.13766, -0.00855, 0.27885,
                                                 -1.13212, -1.50415, -6.34766, 0.55769,  -0.01709};
const std::array<double, 17> kParam2021PiMinus = {1.75169,  0.11144,  47.35877, -4.69434, 1.60552,   0.00800,
                                                  0.44194,  -2.29188, -41.67194, 0.69475, 0.02527,   -0.50178,
                                                  -1.22825, -1.16878, 5.75825,   -1.00355, 0.05055};

// Fortran: physics_pion.f:sig_param_2021 (l.738-805)
double sigParam2021(double thcm, double phicm, double t, double q2, double wsq, double eps, int which_pion,
                    bool doing_pizero, Excmn* excmn) {
    double sigexcl3;
    ExclfitResult r;
    if (doing_pizero) {
        const ExclfitResult rp = exclfit(t, thcm, phicm, q2, wsq, eps, kParam2021PiPlus, 0.0);
        r = exclfit(t, thcm, phicm, q2, wsq, eps, kParam2021PiMinus, 0.0);
        sigexcl3 = (r.sig + rp.sig) / 2.0;   // (sigexcl3m + sigexcl3p)/2.0
    } else {
        if (which_pion == 1 || which_pion == 11 || which_pion == 3) {
            r = exclfit(t, thcm, phicm, q2, wsq, eps, kParam2021PiMinus, 1.0);
        } else {   // pi+
            r = exclfit(t, thcm, phicm, q2, wsq, eps, kParam2021PiPlus, 1.0);
        }
        sigexcl3 = r.sig;
    }
    if (excmn) *excmn = Excmn{r.sigl, r.sigt, r.sigtt, r.siglt};
    return sigexcl3;
}

// ------------------------------------------------------------------------------------------------ MAID
const std::array<const char*, 4> MaidTables::kFiles = {"maidpi0p.dat", "maidpi0n.dat", "maidpipn.dat", "maidpimp.dat"};

namespace {
/// One Fortran F edit-descriptor field read (Fw.d) from a record padded with blanks (PAD='YES', BLANK='NULL': the
/// blanks are ignored; an empty field is 0; without a decimal point d digits are implied).
double readFField(const std::string& rec, std::size_t pos, std::size_t w, int d) {
    std::string s;
    for (std::size_t i = pos; i < pos + w && i < rec.size(); ++i)
        if (rec[i] != ' ') s += rec[i];
    if (s.empty()) return 0.0;
    // std::from_chars (correctly rounded, as strtod: the same double) is several times faster than strtod on the
    // 1.9 million MAID numbers; anything it does not take whole (a '+' sign, ...) still goes through strtod.
    double v = 0.0;
    const char* b = s.data();
    const char* e = s.data() + s.size();
    const auto [ptr, ec] = std::from_chars(b, e, v);
    if (ec != std::errc() || ptr != e) {
        char* end = nullptr;
        v = std::strtod(s.c_str(), &end);
    }
    if (s.find('.') == std::string::npos && s.find_first_of("eEdD") == std::string::npos) v = v / std::pow(10.0, d);
    return v;
}
}  // namespace

// Fortran: physics_pion.f:sigmaid l.616-637 (read(9,'(f11.6,18f8.4)') per (iq, iw, ith))
bool MaidTables::load(const std::string& dir, std::string* err) {
    std::vector<double> tbl(static_cast<std::size_t>(kNipi) * kNq * kNw * kNth * kNcol, 0.0);
    for (int i = 1; i <= kNipi; ++i) {
        const std::string path = dir + "/" + kFiles[static_cast<std::size_t>(i - 1)];
        std::ifstream f(path);
        if (!f) {
            if (err) *err = "sigmaid: cannot open " + path;
            return false;
        }
        std::string rec;
        for (int iq = 1; iq <= kNq; ++iq)
            for (int iw = 1; iw <= kNw; ++iw)
                for (int ith = 1; ith <= kNth; ++ith) {
                    if (!std::getline(f, rec)) {
                        if (err) *err = "sigmaid: " + path + " ends before record (iq,iw,ith) = (" + std::to_string(iq) +
                                        "," + std::to_string(iw) + "," + std::to_string(ith) + ")";
                        return false;
                    }
                    if (!rec.empty() && rec.back() == '\r') rec.pop_back();
                    const std::size_t base =
                        static_cast<std::size_t>(((((i - 1) * kNq + (iq - 1)) * kNw + (iw - 1)) * kNth + (ith - 1)) * kNcol);
                    tbl[base] = readFField(rec, 0, 11, 6);
                    for (int j = 2; j <= kNcol; ++j)
                        tbl[base + static_cast<std::size_t>(j - 1)] = readFField(rec, 11 + 8 * static_cast<std::size_t>(j - 2), 8, 4);
                }
    }
    tbl_ = std::move(tbl);
    return true;
}

// Fortran: physics_pion.f:659-669
MaidBin sigmaidBin(double q2, double w, double costh) {
    static const double cthmin[6] = {-0.20, 0.20, 0.44, 0.63, 0.78, 0.90};
    static const double cthmax[6] = {0.20, 0.44, 0.63, 0.78, 0.90, 1.0};
    MaidBin b;
    b.iq = static_cast<int>((q2 + 0.1) / 0.2);   // Fortran int(): truncation
    b.iq = std::min(25, std::max(1, b.iq));
    b.iw = static_cast<int>((w - 1.090) / 0.020);
    b.iw = std::min(46, std::max(1, b.iw));
    b.ith = 0;
    for (int i = 1; i <= 6; ++i)
        if (costh >= cthmin[i - 1] && costh <= cthmax[i - 1]) b.ith = i;   // the last matching bin wins
    b.ith = std::min(6, std::max(1, b.ith));
    return b;
}

// Fortran: physics_pion.f:sigmaid (l.577-736)
SigmaidResult sigmaid(int ipi, double q2, double w, double e0, double costh, double phi, const MaidTables& maid) {
    const double am = 0.9383;
    SigmaidResult r;   // sig0 = sigz = sigez = sige = 0 (l.640-644)

    if (w < 1.08) return r;

    // get epsilon
    const double nu = (sq(w) - sq(am) + q2) / 2. / am;
    if (nu > e0) return r;
    const double ep = e0 - nu;
    const double sin2 = q2 / 4. / e0 / ep;
    if (sin2 <= 0.0 || sin2 > 1.) return r;
    const double eps = 1. / (1. + 2. * (1. + sq(nu) / q2) * sin2 / (1. - sin2));
    const double sthe = 2. * std::sqrt(sin2);                // sin(theta_electron)
    double pt = ep * sthe / std::sqrt(q2 + sq(nu));          // sin(theta_q)
    const double pl = std::sqrt(1. - sq(pt));                // cos(theta_q)
    const double pn = -pt * std::sin(phi);                   // normal component of the target field
    pt = pt * std::cos(phi);                                 // in-plane component

    const MaidBin b = sigmaidBin(q2, w, costh);
    auto T = [&](int j) { return maid.at(ipi, b.iq, b.iw, b.ith, j); };

    double wfact = 1.;
    if (w > 1.232) wfact = (w - 1.132) / 0.100;
    const double ST = T(1) / std::max(0.2, q2) / wfact;
    const double SL = T(2) * ST;
    const double STL = T(3) * ST;
    const double STT = T(4) * ST;
    const double STLZ = T(5) * ST;
    const double STTZ = T(6) * ST;
    const double STLPZ = T(7) * ST;
    const double STTPZ = T(8) * ST;
    const double STLP = T(9) * ST;
    const double STLX = T(10) * ST;
    const double STLY = T(11) * ST;
    const double STTX = T(12) * ST;
    const double STTY = T(13) * ST;
    const double STLPX = T(14) * ST;
    const double STLPY = T(15) * ST;
    const double STTPX = T(16) * ST;
    const double SLY = T(17) * ST;
    const double STY = T(18) * ST;

    const double CSF = std::cos(phi);
    const double SNF = std::sin(phi);
    const double CS2F = std::cos(2. * phi);
    const double SN2F = std::sin(2. * phi);
    // CST = costh, SNT = sqrt(1-costh**2) are set but not used (l.706-707)

    r.sig0 = ST + eps * SL + std::sqrt(2. * eps * (1. + eps)) * CSF * STL + eps * CS2F * STT;

    r.sigz = std::sqrt(2. * eps * (1. + eps)) * (pt * SNF * STLX + pn * CSF * STLY + pl * SNF * STLZ) +
             eps * (pt * SN2F * STTX + pn * CS2F * STTY + pl * SN2F * STTZ) + pn * (STY + eps * SLY);

    r.sigez = std::sqrt(2. * eps * (1. - eps)) * (pt * STLPX * CSF + pn * STLPY * SNF + pl * STLPZ * CSF) +
              std::sqrt(1. - sq(eps)) * (pt * STTPX + pl * STTPZ);

    r.sige = std::sqrt(2. * eps * (1. - eps)) * SNF * STLP;

    r.sigez = -1.0 * r.sigez;   // DIS sign convention (l.731)
    r.sigz = -1.0 * r.sigz;     // l.734
    return r;
}

// ------------------------------------------------------------------------------------------------ peepi
// Fortran: physics_pion.f:peepi (l.1-194)
double peepi(const Event& vertex, MainEvent& main, const PeepiConfig& cfg, const MaidTables& maid, std::ostream* log) {
    const FermiMomentum& f = main.pfermi;   // COMMON /pfermi/

    // l.47-60: photon-NUCLEON C.M. frame (jacobians.f)
    const reactions::CmTransform cm =
        reactions::transformToCm(vertex, main, f.pfer, f.pferx, f.pfery, f.pferz, f.efer, log);

    main.thetacm = cm.thetacm;
    main.phicm = cm.phicm;
    main.pcm = cm.phadcm;
    main.davejac = cm.jacobian;
    main.johnjac = cm.jac_old;   // approx. assuming collinear boost

    // l.74-79
    double tfcos = f.pferx * vertex.uq.x + f.pfery * vertex.uq.y + f.pferz * vertex.uq.z;
    if (tfcos - 1. > 0. && tfcos - 1. < 1.e-8) tfcos = 1.0;
    const double tfsin = std::sqrt(1. - sq(tfcos));

    const double s = sq(vertex.nu + f.efer) - sq(vertex.q + f.pfer * tfcos) - sq(f.pfer * tfsin);
    main.wcm = std::sqrt(s);

    // l.93 (JRA 2001): k_eq from wcm**2 (not s)
    const double k_eq = (sq(main.wcm) - sq(cfg.Mtar_struck)) / 2. / cfg.Mtar_struck;

    // l.126-127: P. Bosted's 2021 fit (the model choice of the Fortran; sig_blok/param04/param_3000 and the multipole
    // model are commented out, l.114-123 and 157-163)
    main.ntup.sigcm1 = sigParam2021(cm.thetacm, cm.phicm, main.t / 1.e6, vertex.Q2 / 1.e6, s / 1.e6, main.epsilon,
                                    cfg.which_pion, cfg.doing_pizero);
    double sigma_eepi = main.ntup.sigcm1;

    // l.133-154: MAID for W < 2.0 GeV, blended between 1.5 and 1.9 GeV
    if (main.wcm < 2000) {
        const double Q2gev = vertex.Q2 / 1.e6;
        const double Wgev = main.wcm / 1000.0;
        const double cthcm = std::cos(cm.thetacm);
        const double E0 = vertex.Ein / 1000.0;
        int ipi;
        if (cfg.doing_pizero) {
            ipi = 1;                                                 // pi0 from the proton
            if (cfg.which_pion == 1 || cfg.which_pion == 3) ipi = 2;   // pi0 from the neutron
        } else {
            ipi = 3;                                                                        // pi+
            if (cfg.which_pion == 1 || cfg.which_pion == 11 || cfg.which_pion == 3) ipi = 4;   // pi-
        }
        const double sig0 = sigmaid(ipi, Q2gev, Wgev, E0, cthcm, cm.phicm, maid).sig0;
        // mub/dOmega* -> mub/dt/dphi*
        main.ntup.sigcm2 = sig0 / cm.phadcm / cm.qstar / 2.;
        const double fac1 = std::min(1., std::max(0., (Wgev - 1.5) / 0.4));
        sigma_eepi = main.ntup.sigcm1 * fac1 + main.ntup.sigcm2 * (1 - fac1);
    } else if (!cfg.fortran_compat) {
        // FORTRAN-BUG(D-038 F10): original Fortran physics_pion.f:133-154: ntup%sigcm2 is assigned only inside
        // `if(main%wcm.lt.2000)`; corrected because for W >= 2 GeV it kept the value of an earlier event (or of the
        // central event): the MAID model is not evaluated there, so sigcm2 = 0.
        main.ntup.sigcm2 = 0.0;
    }
    main.ntup.sigcm = sigma_eepi;

    // l.186-190: flux factor with the non-collinear correction
    const double fac = 1. / (1. - f.pferz * f.pfer / f.efer) * cfg.Mtar_struck / f.efer;
    const double gtpr = constants::alpha / 2. / (pi * pi) * vertex.e.E / vertex.Ein * k_eq / vertex.Q2 / (1. - main.epsilon);

    return sigma_eepi * cm.jacobian * (gtpr * fac);   // ub/MeV^2/rad -> ub/sr -> ub/MeV/sr^2
}

}  // namespace pion
}  // namespace simc
