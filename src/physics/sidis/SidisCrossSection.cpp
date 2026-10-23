// src/physics/sidis/SidisCrossSection.cpp — semi_physics.f (peepiX, slacemcfit, Rhad_global), task P12.1.
// See include/simc/physics/sidis/SidisCrossSection.h. REAL(8) arithmetic in the source order (D-024).
#include "simc/physics/sidis/SidisCrossSection.h"

#include <cmath>
#include <cstdio>
#include <string>

#include "simc/core/Errors.h"
#include "simc/core/SimcConstants.h"
#include "simc/physics/sidis/F1F2In21.h"

namespace simc {
namespace sidis {

namespace {
inline double sq(double x) { return x * x; }
using constants::alpha;
using constants::hbarc;
using constants::pi;

// list-directed REAL(8) output, close to gfortran's (G0-like with 17 significant digits)
std::string fnum(double v) {
    char b[48];
    std::snprintf(b, sizeof b, "%25.17g", v);
    return b;
}

// semi_physics.f:107-112 PB fit of 9/20/2021 (DATA, REAL(8) literals)
constexpr double pf[13] = {0,       1.0424, -0.1714, 1.8960, -0.0307, 0.1636, -0.1272,
                           -4.2093, 5.0103, 2.7406,  -0.5778, 3.5292, 7.3910};
constexpr double pu[13] = {0,       0.7840, 0.2369, 1.4238, 0.1484, 0.1518, -1.2923,
                           -1.5710, 3.0305, 1.1995, 1.3553, 2.5868, 8.0666};
}  // namespace

// Fortran: semi_physics.f:slacemcfit (l.634-650)
double slacemcfit(double A, double x) {
    double r = 1.0;
    double Atmp = 1.0;
    if (A > 2.0) {
        Atmp = A;
        // x**2 ... x**8: gcc powi chain (x3 = x*x2, x4 = x2*x2, x5 = x2*x3, x6 = x3*x3, x7 = x3*x4, x8 = x4*x4)
        const double x2 = x * x, x3 = x * x2, x4 = x2 * x2, x5 = x2 * x3, x6 = x3 * x3, x7 = x3 * x4, x8 = x4 * x4;
        const double al = -0.070 + 2.189 * x - 24.667 * x2 + 145.291 * x3 - 497.237 * x4 + 1013.129 * x5 -
                          1208.393 * x6 + 775.767 * x7 - 205.872 * x8;
        const double lx = std::log(x);
        const double C = std::exp(0.017 + 0.018 * lx + 0.005 * (lx * lx));
        r = C * std::pow(Atmp, al);
    }
    return r;
}

// Fortran: semi_physics.f:Rhad_global (l.652-680)
double rhadGlobal(double A, double z) {
    const double Ahyd = 1.0;
    // FORTRAN-DIFF-CANDIDATE(D-051 F10): only A = 1 returns 1; deuterium gets the attenuation fit (0.955 at z = 0.5).
    if (A == Ahyd) return 1.0;
    const double Nzero = 0.98883 - 0.0038309 * A + 0.10841E-4 * (A * A);
    const double alphah = 0.31953E-01 - 0.18659E-02 * A + 0.51747E-05 * (A * A);
    const double Atmp = A < 83.8 ? A : 83.8;
    const double betah = 0.85475E-02 + 0.12763E-02 * Atmp - 0.24451E-05 * (Atmp * Atmp);
    return Nzero * std::pow(z, alphah) * std::pow(1 - z, betah);
}

// Fortran: semi_physics.f:peepiX (l.1-631)
double peepiX(const Event& vertex, const Event& vertex0, MainEvent& main, double& survivalprob, bool doing_cent,
              const PeepiXConfig& cfg, Cteq5Pdf& pdf, const FdssGrid* fdss, PeepiXState& st, std::ostream* log,
              PeepiXTrace* trace) {
    using constants::Mk;
    using constants::Mp;
    using constants::Mpi;
    using constants::Mpi0;
    const FermiMomentum& f = main.pfermi;   // COMMON /pfermi/
    constexpr double qu = 2. / 3., qd = -1. / 3., qs = -1. / 3.;   // l.121-124
    PeepiXTrace tr;

    // l.141: if(first) FINI=0 — the fDSS grid is (re)initialised on the first call that reaches fDSS; the grid object
    // is loaded once by the model, which is the same (one grid file per run).
    const double targA = cfg.targA;
    const double targZ = cfg.targZ;
    const double targN = targA - targZ;

    double Mpi_gev = Mpi / 1000.0;
    if (cfg.doing_pizero) Mpi_gev = Mpi0 / 1000.0;
    const double Mp_gev = Mp / 1000.0;

    double nu, qx, qy, qz, Q2, Eb, Eprime;
    double& pt2 = st.pt2;    // SAVEd (only matters for doing_cent with sigc_flag not 0/1)
    double& zhad = st.zhad;
    if (doing_cent) {
        // l.157-210 (dead code in SIMC: event.f:1523 is commented out)
        nu = vertex0.nu;
        qx = vertex0.uq.x * vertex0.q;
        qy = vertex0.uq.y * vertex0.q;
        qz = vertex0.uq.z * vertex0.q;
        Q2 = vertex0.Q2;
        Eb = vertex0.Ein;
        Eprime = vertex0.e.E;
        if (cfg.sigc_flag == 0) {   // binning in z
            pt2 = cfg.sigc_kin_ind * 1.e6;
            zhad = 0.0;
            for (int i = 1; i <= cfg.sigc_nbin; ++i) {
                const double klo = cfg.sigc_kin_min + (i - 1) * (cfg.sigc_kin_max - cfg.sigc_kin_min) / cfg.sigc_nbin;
                const double khi = cfg.sigc_kin_min + i * (cfg.sigc_kin_max - cfg.sigc_kin_min) / cfg.sigc_nbin;
                if (vertex.zhad > klo && vertex.zhad <= khi) zhad = klo + (khi - klo) / 2.;
            }
        } else if (cfg.sigc_flag == 1) {   // binning in pt2
            zhad = cfg.sigc_kin_ind;
            pt2 = 0.0;
            for (int i = 1; i <= cfg.sigc_nbin; ++i) {
                double klo = cfg.sigc_kin_min + (i - 1) * (cfg.sigc_kin_max - cfg.sigc_kin_min) / cfg.sigc_nbin;
                klo = klo * 1.e6;
                double khi = cfg.sigc_kin_min + i * (cfg.sigc_kin_max - cfg.sigc_kin_min) / cfg.sigc_nbin;
                khi = khi * 1.e6;
                if (vertex.pt2 > klo && vertex.pt2 <= khi) pt2 = klo + (khi - klo) / 2.;
            }
        }
        if (st.first_cent) {
            if (log) {
                *log << " Central Kinematics:\n";
                *log << " Ebeam (GeV):" << fnum(Eb / 1000) << "\n";
                *log << " nu (GeV)   :" << fnum(nu / 1000) << "\n";
                *log << " Q2 (GeV2)  :" << fnum(Q2 / 1e6) << "\n";
                *log << " theta_pq   :" << fnum(main.theta_pq) << "\n";
                if (cfg.sigc_flag == 0) {
                    *log << " Pt2 (GeV2) :" << fnum(pt2 / 1e6) << "\n";
                    *log << " Binning in z from" << fnum(cfg.sigc_kin_min) << "to" << fnum(cfg.sigc_kin_max) << "\n";
                } else if (cfg.sigc_flag == 1) {
                    *log << " z       :" << fnum(zhad) << "\n";
                    *log << " Binning in pt2 from" << fnum(cfg.sigc_kin_min) << "to" << fnum(cfg.sigc_kin_max)
                         << "GeV2\n";
                }
            }
            st.first_cent = false;
        }
    } else {
        nu = vertex.nu;
        qx = vertex.uq.x * vertex.q;
        qy = vertex.uq.y * vertex.q;
        qz = vertex.uq.z * vertex.q;
        Q2 = vertex.Q2;
        Eb = vertex.Ein;
        Eprime = vertex.e.E;
        pt2 = vertex.pt2;
        zhad = vertex.zhad;
    }
    // Epx, Epy, Epz (l.165-167, 219-221) are computed and never used

    double mhad;
    if (cfg.doing_semipi) {
        mhad = Mpi;
        if (cfg.doing_pizero && cfg.pizero_fixed)
            // FORTRAN-BUG(D-051 F4): original Fortran semi_physics.f:227: `mhad = Mpi` also for doing_pizero; corrected
            // because Mpi_gev (l.149) and dbase.f:194 (Mh) use Mpi0 for the pi0.
            mhad = Mpi0;
    } else if (cfg.doing_semika) {
        mhad = Mk;
    } else {
        if (log) *log << " semi_physics error: unknown hadron type\n";
        throw FortranStop("semi_physics.f:231-232 STOP: semi_physics error: unknown hadron type");
    }

    const double mtar = cfg.Mtar_struck;
    const double Ehad = zhad * nu;
    const double phad = std::sqrt(sq(Ehad) - sq(mhad));
    double cthpq;
    if (doing_cent)
        cthpq = std::sqrt(1. - pt2 / sq(phad));
    else
        cthpq = std::cos(vertex.theta_pq);

    // sx = (2.*Eb*mtar + mtar**2)/1.e6 (l.246) is never used
    double xbj;
    if (cfg.do_fermi) {   // xbj = Q2/(2 P.q)
        xbj = Q2 / 2. / (f.efer * nu - std::abs(f.pfer) * (f.pferx * qx + f.pfery * qy + f.pferz * qz));
        if (!doing_cent) main.ntup.xfermi = xbj;
        if (cfg.nonpositive_xbj_zero && !(xbj > 0.)) {
            // FORTRAN-BUG(D-051 F2): original Fortran semi_physics.f:247-248 passes the negative x_Bj of an unphysical
            // P.q <= 0 (large pfer along q) to Ctq5Pdf, which STOPs the run (Ctq5Pdf.f:80-83: no D run with
            // do_fermi = 1 can finish); corrected: the event gets cross section 0 (early return, as l.280-284).
            tr.xbj_raw = xbj;
            tr.xbj = xbj;
            tr.threshold_return = true;
            if (trace) *trace = tr;
            return 0.0;
        }
    } else {
        xbj = Q2 / 2. / mtar / nu;
    }
    tr.xbj_raw = xbj;
    if (xbj > 1.0) {
        if (log) *log << " XBj is too large!" << fnum(xbj) << "\n";
        xbj = 1.0;
    }

    // convert some stuff to GeV
    const double Q2gev = Q2 / 1.e6;
    double Qgev = std::sqrt(Q2gev);
    const double pt2gev = pt2 / 1.e6;

    double wsq = sq(Mp_gev) + Q2gev * (1. / xbj - 1.);
    const double w = std::sqrt(wsq);
    // zero below the 2 pi / pi K threshold (PB)
    const double mtargev = mtar / 1000.;
    const double nugev = nu / 1000.;
    // FORTRAN-DIFF-CANDIDATE(D-051 F7): a mass^2 multiplied by (1 - pt2 [GeV^2]); the intended missing-mass formula
    // is unknown, kept.
    const double mmpi2 = sq(mtargev) + 2. * mtargev * nugev * (1 - zhad) * (1 - pt2gev);
    tr.Mpi_gev = Mpi_gev;
    tr.mhad = mhad;
    tr.mtar = mtar;
    tr.xbj = xbj;
    tr.Q2gev = Q2gev;
    tr.Qgev = Qgev;
    tr.pt2gev = pt2gev;
    tr.wsq = wsq;
    tr.w = w;
    tr.mmpi2 = mmpi2;
    tr.zhad = zhad;
    tr.pt2 = pt2;
    tr.Ehad = Ehad;
    tr.phad = phad;
    tr.cthpq = cthpq;
    if (mmpi2 < sq(mtargev + mhad / 1000.)) {
        // FORTRAN-QUIRK(P12.1 Q1): early return: davejac, ntup%sigcm and the survival probability (decdist) keep the
        // values of the previous event (the weight is 0).
        tr.threshold_return = true;
        if (trace) *trace = tr;
        return 0.0;
    }

    // l.288-296: the first call of sqesub only reads the deuteron smearing tables (f1f2tables/*cdbonn*) into SAVEd
    // arrays that F1F2IN21 uses for A = 2 only; peepiX calls F1F2IN21 with A = 1: no effect, not ported.
    if (st.firstqe) st.firstqe = false;

    // Get the PDFs (l.299-320): SetCtq5(iset = 1) is the model's table (loaded at init)
    if (st.first) st.first = false;
    if (!pdf.loaded()) throw std::logic_error("peepiX: CTEQ5 table not loaded (SetCtq5)");
    const double u = pdf.pdf(1, xbj, Qgev, log);
    const double ubar = pdf.pdf(-1, xbj, Qgev, log);
    const double d = pdf.pdf(2, xbj, Qgev, log);
    const double dbar = pdf.pdf(-2, xbj, Qgev, log);
    const double s = pdf.pdf(3, xbj, Qgev, log);
    const double sbar = pdf.pdf(-3, xbj, Qgev, log);

    const double uA = targZ * u + targN * d;
    const double ubarA = targZ * ubar + targN * dbar;
    const double dA = targZ * d + targN * u;
    const double dbarA = targZ * dbar + targN * ubar;
    const double sA = targZ * s + targN * s;
    const double sbarA = targZ * sbar + targN * sbar;

    const double sum_sq = sq(qu) * (uA + ubarA) + sq(qd) * (dA + dbarA) + sq(qs) * (sA + sbarA);
    tr.u = u;
    tr.ubar = ubar;
    tr.d = d;
    tr.dbar = dbar;
    tr.s = s;
    tr.sbar = sbar;
    tr.sum_sq = sum_sq;
    tr.uA = uA;
    tr.ubarA = ubarA;
    tr.dA = dA;
    tr.dbarA = dbarA;
    tr.sA = sA;
    tr.sbarA = sbarA;
    tr.Qgev_after = Qgev;

    // Peter Bosted's parameterization (l.460-495) / fDSS (l.496-506)
    double u1 = 0, ub = 0, d1 = 0, db = 0, s1 = 0, sb = 0;
    if (cfg.doing_semipi) {
        const double xp = 2. * xbj / (1. + std::sqrt(1. + 4. * sq(xbj) * sq(Mp_gev) / Q2gev));
        const double zp = (zhad / 2.) * (xp / xbj) *
                          (1. + std::sqrt(1 - 4 * sq(xbj) * sq(Mp_gev) * (sq(Mpi_gev) + pt2gev) / sq(zhad) / sq(Q2gev)));
        const double sv = std::log(Q2gev / 2.);
        const double zp2 = zp * zp, zp3 = zp * zp2;
        double yf = pf[1] * std::pow(zp, pf[2] + pf[4] * sv + pf[9] / w) * std::pow(1. - zp, pf[3] + pf[5] * sv + pf[10] / w);
        yf = yf * (1. + pf[6] * zp + pf[7] * zp2 + pf[8] * zp3) * (1. + pf[11] / w + pf[12] / (w * w));
        double yu = pu[1] * std::pow(zp, pu[2] + pu[4] * sv + pu[9] / w) * std::pow(1. - zp, pu[3] + pu[5] * sv + pu[10] / w);
        yu = yu * (1. + pu[6] * zp + pu[7] * zp2 + pu[8] * zp3) * (1. + pu[11] / w + pu[12] / (w * w));
        if (cfg.doing_hplus) {
            u1 = yf;
            d1 = yu;
        } else {
            u1 = yu;
            d1 = yf;
        }
        ub = d1;
        db = u1;
        s1 = yu;
        sb = s1;
        if (cfg.doing_pizero) {   // just get average of pi+ and pi-
            u1 = (yf + yu) / 2.;
            d1 = (yf + yu) / 2.;
            ub = (yf + yu) / 2.;
            if (cfg.pizero_fixed)
                // FORTRAN-BUG(D-051 F4): original Fortran semi_physics.f:494: `d1 = (yf + yu) / 2.` (d1 a second time);
                // corrected because the comment and the three other lines average every light (anti)quark FF, and db
                // keeps the pi+ (doing_hplus) or pi- value of l.487 otherwise.
                db = (yf + yu) / 2.;
            else
                d1 = (yf + yu) / 2.;
        }
        tr.xp = xp;
        tr.zp = zp;
        tr.sv = sv;
        tr.yf = yf;
        tr.yu = yu;
    } else if (cfg.doing_semika) {
        if (!fdss || !fdss->loaded()) throw std::logic_error("peepiX: fDSS kaon grid not loaded");
        const int ICdss = cfg.doing_hplus ? 1 : -1;   // IHdss = 2 (kaons), IOdss = 1 (NLO)
        const FdssValues v = fdss->eval(ICdss, zhad, Q2gev, log);
        u1 = v.U;
        ub = v.UB;
        d1 = v.D;
        db = v.DB;
        s1 = v.S;
        sb = v.SB;
    }
    tr.u1 = u1;
    tr.ub = ub;
    tr.d1 = d1;
    tr.db = db;
    tr.s1 = s1;
    tr.sb = sb;

    double dsigdz;
    if (cfg.sum_sq_guard && !(sum_sq > 0.))
        // FORTRAN-BUG(D-051 F1): original Fortran semi_physics.f:508-513 divides by sum_sq without the guard of the old
        // code (l.378-400, commented out): all PDFs 0 (x_Bj clamped to 1, Q < Alambda) gives 0/0 = NaN in the weight;
        // corrected as the old guard: dsigdz = 0.
        dsigdz = 0.0;
    else
        dsigdz = (sq(qu) * uA * u1 + sq(qu) * ubarA * ub + sq(qd) * dA * d1 + sq(qd) * dbarA * db + sq(qs) * sA * s1 +
                  sq(qs) * sbarA * sb) /
                 sum_sq / zhad;

    // FORTRAN-QUIRK(D-051 F5): the input pt_b_param is not used (l.142 commented out); b is this z-dependent fit.
    const double b = 1. / (0.120 * sq(zhad) + 0.200);
    const double rhad = rhadGlobal(targA, zhad);
    const double sighad = rhad * dsigdz * b * std::exp(-b * pt2gev) / 2. / pi;

    wsq = sq(Mp_gev) + Q2gev * (1. / xbj - 1.);
    double F1, F2, F1p, F2p, F1n, F2n;
    // Updated to F1F2IN21 (Eric Christy's 2021 fit)
    f1f2in21(1.0, 1.0, Q2gev, wsq, F1, F2);
    F1p = F1;
    F2p = F2;
    f1f2in21(0.0, 1.0, Q2gev, wsq, F1, F2);
    F1n = F1;
    F2n = F2;
    const double emc = slacemcfit(targA, xbj);
    F1 = (targZ * F1p + targN * F1n) * emc;
    F2 = (targZ * F2p + targN * F2n) * emc;

    const double W1 = F1 / (mtar / 1000.);
    const double W2 = F2 / (nu / 1000.);
    const double sin2th2 = Q2 / 4. / Eb / Eprime;
    const double cos2th2 = 1. - sin2th2;
    const double W2coeff = cos2th2;
    const double sige = 4. * sq(alpha) * sq(Eprime / 1000) / sq(Q2gev) * (W2 * W2coeff + 2. * W1 * sin2th2);

    // dsig/(dOmega_e dE_e dz dpt**2 dPhi_had) in microbarn/GeV**3/sr
    const double sigsemi = sige * sighad * sq(hbarc / 1000.) * 10000.0;
    // to dsig/(dOmega_e dE_e dE_h dCos(theta) dPhi_had): 1/omega * 2*p_h**2*cos(theta)
    const double jacobian = 1. / (nu / 1000.) * 2. * sq(phad / 1000.) * cthpq;
    // 1.e6: microbarn/GeV^2 -> microbarn/MeV^2
    double sigma_eepiX = sigsemi * jacobian / 1.e6;
    // flux factor of the Fermi-smeared cross section (l.567-584)
    double fac;
    if (cfg.do_fermi)
        fac = 1. / (1. - f.pferz * f.pfer / f.efer) * mtar / f.efer;
    else
        fac = 1.0;
    sigma_eepiX = sigma_eepiX * fac;

    if (doing_cent) {
        main.johnjac = jacobian * 1000.0;
    } else {
        main.davejac = jacobian * 1000.0;
        main.ntup.sigcm = sighad;
    }
    const double peepix = sigma_eepiX;
    tr.dsigdz = dsigdz;
    tr.b = b;
    tr.rhad = rhad;
    tr.sighad = sighad;
    tr.F1p = F1p;
    tr.F2p = F2p;
    tr.F1n = F1n;
    tr.F2n = F2n;
    tr.emc = emc;
    tr.F1 = F1;
    tr.F2 = F2;
    tr.W1 = W1;
    tr.W2 = W2;
    tr.sin2th2 = sin2th2;
    tr.sige = sige;
    tr.sigsemi = sigsemi;
    tr.jacobian = jacobian;
    tr.fac = fac;
    tr.sigma = sigma_eepiX;

    // l.600-627: survival probability for main.weight when the decay is not simulated (applied to the weight by
    // complete_main for doing_semika only, event.f:1560)
    if (!cfg.doing_decay) {
        if (cfg.hadron_arm == 1) {
            st.zaero = -331.491;   // aero at 40.199 cm, last project=371.69 cm
        } else if (cfg.hadron_arm == 2) {
            st.zaero = -82.8;   // From Rick: aero at 263.2cm,last project=346(cal).
        } else if (cfg.hadron_arm == 3) {
            st.zaero = -183.;   // aero at 130cm,last project=313(S2)
        } else if (cfg.hadron_arm == 4) {
            st.zaero = -183.;
        } else if (cfg.survival_path_fixed) {
            // FORTRAN-BUG(D-051 F6, as D-043 F6): original Fortran semi_physics.f:611-620 has no branch for hadron_arm 5
            // (SHMS): zaero keeps the SAVEd value of an earlier call (0 in a SIMC run); corrected because the path must
            // not depend on earlier calls: 0 (full length to the last plane).
            st.zaero = 0.;
        }
        double pathlen;
        if (cfg.survival_path_fixed)
            // FORTRAN-BUG(D-051 F6, as D-043 F6): original Fortran semi_physics.f:621: zaero*(1+dx**2+dy**2);
            // corrected because the path length along the track is zaero*sqrt(1+dx**2+dy**2).
            pathlen = main.FP.p.path + st.zaero * std::sqrt(1 + sq(main.FP.p.dx) + sq(main.FP.p.dy));
        else
            pathlen = main.FP.p.path + st.zaero * (1 + sq(main.FP.p.dx) + sq(main.FP.p.dy));
        const double p_kaon = cfg.spec_p_P * (1. + main.SP.p.delta / 100.);
        double betak;
        if (cfg.survival_p_hadron)
            // FORTRAN-BUG(D-051 F6, as D-043 F5): original Fortran semi_physics.f:623: betak =
            // spec%p%P/sqrt(spec%p%P**2+Mh2) (the central momentum; p_kaon is computed and unused, against the comment
            // of l.607-608); corrected because the decay probability depends on the momentum of this hadron.
            betak = p_kaon / std::sqrt(sq(p_kaon) + cfg.Mh2);
        else
            betak = cfg.spec_p_P / std::sqrt(sq(cfg.spec_p_P) + cfg.Mh2);
        const double gammak = 1. / std::sqrt(1. - sq(betak));
        survivalprob = 1. / std::exp(pathlen / (cfg.ctau * betak * gammak));
        main.ntup.decdist = survivalprob;   // decdist in ntuple
        tr.pathlen = pathlen;
        tr.betak = betak;
        tr.gammak = gammak;
        tr.survival = survivalprob;
    }
    if (trace) *trace = tr;
    return peepix;
}

}  // namespace sidis
}  // namespace simc
