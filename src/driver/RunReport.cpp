// src/driver/RunReport.cpp — see include/simc/driver/RunReport.h.
// Fortran: simc.f l.360-432 (normalisation), l.462-532 (.geni), report (l.644-1142).
#include "simc/driver/RunReport.h"

#include <cmath>
#include <ctime>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "simc/core/Errors.h"
#include "simc/core/SimcConstants.h"
#include "simc/core/SpectrometerId.h"
#include "simc/driver/FortranWriter.h"
#include "simc/core/Units.h"

namespace simc {
namespace driver {

using fortran::List;
using fortran::Record;
using constants::degrad;

namespace {
// gfortran MAX(a, b) for reals: b if b > a or a is NaN
inline double fmaxF(double a, double b) { return (b > a || std::isnan(a)) ? b : a; }
inline double fminF(double a, double b) { return (b < a || std::isnan(a)) ? b : a; }
}  // namespace

// Fortran: simc.f l.360-432
Normalisation Normalise(const RunSetup& s, const RunTotals& t, double genvol) {
    Normalisation n;
    n.targetfac = s.targetfac;
    n.luminosity = s.luminosity;
    n.charge = s.params.EXPER.charge;
    // l.366 normfac = luminosity/ntried*nevent
    n.normfac = n.luminosity / static_cast<double>(t.ntried) * static_cast<double>(t.nevent);
    // l.374-394 (IReaction::genVolume)
    n.genvol = genvol;
    // l.396-398
    n.normfac = n.normfac * n.genvol;
    if (s.flags.doing_phsp) n.normfac = 1.0;
    n.wtcontribute = t.wtcontribute * n.normfac;

    // l.405-432 resolutions
    if (t.npasscuts > 1) {
        const double tmpnum = static_cast<double>(t.npasscuts);
        auto one = [&](const TwoArmSums::Arm& s1, const TwoArmSums::Arm& s2, TwoArmSums::Arm& ave, TwoArmSums::Arm& res) {
            ave.delta = s1.delta / tmpnum;
            ave.xptar = s1.xptar / tmpnum;
            ave.yptar = s1.yptar / tmpnum;
            ave.ytar = s1.ytar / tmpnum;
            auto r = [&](double a, double b) {
                const double m = a / tmpnum;
                return std::sqrt(fmaxF(0., (b / tmpnum) - m * m));
            };
            res.delta = r(s1.delta, s2.delta);
            res.xptar = r(s1.xptar, s2.xptar);
            res.yptar = r(s1.yptar, s2.yptar);
            res.ytar = r(s1.ytar, s2.ytar);
        };
        one(t.sumerr.e, t.sumerr2.e, n.aveerr.e, n.resol.e);
        one(t.sumerr.p, t.sumerr2.p, n.aveerr.p, n.resol.p);
    }
    return n;
}

std::string FortranCtime(long epoch_seconds) {
    std::time_t tt = static_cast<std::time_t>(epoch_seconds);
    char buf[64];
    std::tm tmv;
    localtime_r(&tt, &tmv);
    std::strftime(buf, sizeof buf, "%a %b %e %H:%M:%S %Y", &tmv);   // = ctime() without the newline
    return buf;
}

namespace {

class Out {
public:
    void line(const Record& r) { os_ << r.str() << "\n"; }
    void line(const std::string& s) { os_ << s << "\n"; }
    void list(const std::string& s) { line(List(s)); }   // write(iun,*) 'text'
    std::string str() const { return os_.str(); }

private:
    std::ostringstream os_;
};

// simc.f:report l.667-788: the reaction lines
void reactionLines(Out& o, const RunSetup& s) {
    const ReactionFlags& f = s.flags;
    const SimcParameters& p = s.params;
    const double A = s.targ.A;
    auto stop = [](const std::string& m) { throw FortranStop(m); };
    if (f.doing_eep) {
        if (f.doing_hyd_elast) o.list("              ****--------  H(e,e'p)  --------****");
        else if (f.doing_deuterium) o.list("              ****--------  D(e,e'p)  --------****");
        else if (f.doing_heavy) o.list("              ****--------  A(e,e'p)  --------****");
        else stop("I don't have ANY idea what (e,e'p) we're doing!!!");
    } else if (f.doing_semi) {
        if (f.doing_semipi) {
            if (A == 1) o.list(f.doing_hplus ? " ****--------  H(e,e'pi+)X  --------****" : " ****--------  H(e,e'pi-)X  --------****");
            else if (A == 2) o.list(f.doing_hplus ? " ****--------  D(e,e'pi+)X  --------****" : " ****--------  D(e,e'pi-)X  --------****");
            else if (A >= 3) o.list(f.doing_hplus ? " ****--------  A(e,e'pi+)X  --------****" : " ****--------  A(e,e'pi-)X  --------****");
            else stop("I don't have ANY idea what A(e,e'pi)X we're doing!!!");
        } else if (f.doing_semika) {
            if (A == 1) o.list(f.doing_hplus ? " ****--------  H(e,e'k+)X  --------****" : " ****--------  H(e,e'k-)X  --------****");
            else if (A == 2) o.list(f.doing_hplus ? " ****--------  D(e,e'k+)X  --------****" : " ****--------  D(e,e'k-)X  --------****");
            else stop("I don't have ANY idea what A(e,e'k)X we're doing!!!");
        } else {
            stop("I don't have ANY idea what A(e,e'x)X we're doing!!!");
        }
    } else if (f.doing_rho) {
        if (A == 1) o.list("              ****--------  H(e,e'rho)  --------****");
        else o.list("I am not set up for anything else yet!");
    } else if (f.doing_delta) {
        // FORTRAN-QUIRK simc.f:726-734: the delta titles go to unit 6, not to the .hist
        if (f.doing_hyddelta) std::cout << "  ****--------  H(e,e'p)pi  --------****\n";
        else if (f.doing_deutdelta) std::cout << "  ****--------  D(e,e'p)pi  --------****\n";
        else if (f.doing_hedelta) std::cout << "  ****--------  A(e,e'p)pi  --------****\n";
        else stop("I don't have ANY idea what (e,e'p)pi we're doing!!!");
    } else if (f.doing_pion) {
        if (f.doing_hydpi) {
            if (A == 1) o.list("              ****--------  H(e,e'pi)  --------****");
            else if (A >= 3) o.list("              ****--------  A(e,e'pi)  --------****");
        } else if (f.doing_deutpi) {
            o.list("              ****--------  D(e,e'pi)  --------****");
        } else if (f.doing_hepi) {
            o.list("              ****--------  A(e,e'pi)  --------****");
        } else {
            stop("I don't have ANY idea what (e,e'pi) we're doing!!!");
        }
        if (p.which_pion == 0 || p.which_pion == 1) o.list("              ****----  Default Final State ----****");
        else if (p.which_pion == 10) o.list("              ****----  Final State is A + pi ----****");
        else if (p.which_pion == 2 || p.which_pion == 3) o.list("              ****----  Final State is pi + Delta ----****");
    } else if (f.doing_kaon) {
        if (f.doing_hydkaon) {
            if (A == 1) o.list("              ****--------  H(e,e'K)  --------****");
            else if (A >= 3) o.list("              ****--------  A(e,e'K)  --------****");
        } else if (f.doing_deutkaon) {
            o.list("              ****--------  D(e,e'K)  --------****");
        } else if (f.doing_hekaon) {
            o.list("              ****--------  A(e,e'K)  --------****");
        } else {
            stop("I don't have ANY idea what (e,e'K) we're doing!!!");
        }
        switch (p.which_kaon) {
            case 0: o.list("              ****---- producing a LAMBDA ----****"); break;
            case 1: o.list("              ****---- producing a SIGMA0 ----****"); break;
            case 2: o.list("              ****---- producing a SIGMA- ----****"); break;
            case 10: o.list("              ****---- WITH BOUND LAMBDA  ----****"); break;
            case 11: o.list("              ****---- WITH BOUND SIGMA0  ----****"); break;
            case 12: o.list("              ****---- WITH BOUND SIGMA-  ----****"); break;
            default: stop("I don't have ANY idea what (e,e'K) we're doing!!!");
        }
    } else if (f.doing_phsp) {
        o.list("              ****--- PHASE SPACE - NO physics, NO radiation (may not work)---****");
    } else {
        stop("I don't have ANY idea what we're doing!!!");
    }
}

// '(9x,a12,'' = '',f15.4,2x,a10)'
Record kin1(const std::string& label, double v, const std::string& unit) {
    Record r;
    r.x(9).a(label, 12).lit(" = ").f(v, 15, 4).x(2).a(unit, 10);
    return r;
}
// '(9x,a12,'' = '',2(f15.4,2x),2x,a5)'
Record kin2(const std::string& label, double e, double p, const std::string& unit) {
    Record r;
    r.x(9).a(label, 12).lit(" = ").f(e, 15, 4).x(2).f(p, 15, 4).x(2).x(2).a(unit, 5);
    return r;
}
// '(17x,a10,'' = '',f15.4,2x,a9)'
Record cen1(const std::string& label, double v, const std::string& unit) {
    Record r;
    r.x(17).a(label, 10).lit(" = ").f(v, 15, 4).x(2).a(unit, 9);
    return r;
}
// 9911 format(2x,2(5x,a10,' = ',e12.6,1x,a5))
Record targ2(const std::string& l1, double v1, const std::string& u1, const std::string& l2, double v2,
             const std::string& u2) {
    Record r;
    r.x(2).x(5).a(l1, 10).lit(" = ").e(v1, 12, 6).x(1).a(u1, 5);
    r.x(5).a(l2, 10).lit(" = ").e(v2, 12, 6).x(1).a(u2, 5);
    return r;
}
// 9912 format(1x,a15,3f15.5,2x,a6)
Record targ3(const std::string& l, const TargetFunction& tf, const std::string& u) {
    Record r;
    r.x(1).a(l, 15).f(tf.ave, 15, 5).f(tf.min, 15, 5).f(tf.max, 15, 5).x(2).a(u, 6);
    return r;
}
// 9913 format(1x,a15,t25,f15.5,2x,a6)
Record targ1(const std::string& l, double v, const std::string& u) {
    Record r;
    r.x(1).a(l, 15).t(25).f(v, 15, 5).x(2).a(u, 6);
    return r;
}
// '(5x,3(2x,a19,''='',l2))' with 2 or 3 items
Record flags(std::initializer_list<std::pair<const char*, bool>> items) {
    Record r;
    r.x(5);
    for (const auto& it : items) r.x(2).a(it.first, 19).lit("=").l(it.second, 2);
    return r;
}
// 9917 format(1x,a18,t21,2f12.3,t50,2f10.3,2x,a5)
Record lim(const std::string& l, double umin, double umax, double lo, double hi, const std::string& u) {
    Record r;
    r.x(1).a(l, 18).t(21).f(umin, 12, 3).f(umax, 12, 3).t(50).f(lo, 10, 3).f(hi, 10, 3).x(2).a(u, 5);
    return r;
}
Record limHeader1() {   // '(t25,2x,a16,2x,t50,2x,a16,2x)'
    Record r;
    r.t(25).x(2).a("______used______", 16).x(2).t(50).x(2).a("_____found______", 16).x(2);
    return r;
}
Record limHeader2() {   // '(t25,2a10,t50,2a10)'
    Record r;
    r.t(25).a("min", 10).a("max", 10).t(50).a("lo", 10).a("hi", 10);
    return r;
}
// 9918 format(1x,a10,a12,t25,3f10.3,2x,a5)
Record slopLine(const std::string& a, const std::string& b, double used, double lo, double hi, const std::string& u) {
    Record r;
    r.x(1).a(a, 10).a(b, 12).t(25).f(used, 10, 3).f(lo, 10, 3).f(hi, 10, 3).x(2).a(u, 5);
    return r;
}
// 9914 format(1x,a18,' = ',4f11.3)
Record rad4(const std::string& l, std::initializer_list<double> v) {
    Record r;
    r.x(1).a(l, 18).lit(" = ");
    for (double x : v) r.f(x, 11, 3);
    return r;
}
// 9915 format(12x,a14,' = ',e16.6,1x,a6)
Record misc(const std::string& l, double v, const std::string& u) {
    Record r;
    r.x(12).a(l, 14).lit(" = ").e(v, 16, 6).x(1).a(u, 6);
    return r;
}
// 99165 format(2x,a22,2f12.5,2x,a)
Record recon(const std::string& l, double a, double b, const std::string& u) {
    Record r;
    r.x(2).a(l, 22).f(a, 12, 5).f(b, 12, 5).x(2).lit(u);
    return r;
}
// '(9x,a25,'' = '',2(2x,f15.4),a5)'
Record spedge(const std::string& l, const Cut& c, const std::string& u) {
    Record r;
    r.x(9).a(l, 25).lit(" = ").x(2).f(c.min, 15, 4).x(2).f(c.max, 15, 4).a(u, 5);
    return r;
}

}  // namespace

// Fortran: simc.f:report (l.644-1142)
std::string FormatHist(const RunSetup& s, const RunReport& r) {
    const SimcParameters& p = s.params;
    const ReactionFlags& f = s.flags;
    const RunTotals& t = r.totals;
    const Normalisation& n = r.norm;
    const reactions::Contrib& c = t.contrib;
    const CentralEvent& central = r.central;
    Out o;

    // Run time (l.663-664)
    o.line("");
    // character*30 timestring (left-justified, blank-padded) written with a30
    auto t30 = [](std::string v) { v.resize(30, ' '); return v; };
    { Record x; x.x(1).lit("BEGIN Time: ").a(t30(r.begin_time), 30); o.line(x); }
    { Record x; x.x(1).lit("END Time:   ").a(t30(r.end_time), 30); o.line(x); }

    // Kinematics (l.667-815)
    o.list("KINEMATICS:");
    reactionLines(o, s);
    o.line(kin1("Ebeam", s.Ebeam, "MeV"));
    o.line(kin1("(dE/E)beam", s.dEbeam / s.Ebeam, "(full wid)"));
    o.line(kin1("x-width", s.gen.xwid, "cm"));
    o.line(kin1("y-width", s.gen.ywid, "cm"));
    { Record x; x.x(9).a("fr_pattern", 12).lit(" = ").i(s.targ.fr_pattern, 15).x(2).a("1=square,2=circ", 16); o.line(x); }
    o.line(kin1("fr1", s.targ.fr1, "cm"));
    o.line(kin1("fr2", s.targ.fr2, "cm"));
    o.list(" ");
    { Record x; x.x(9).x(18).a("____E arm____", 15).x(2).a("____P arm____", 15).x(2); o.line(x); }
    o.line(kin2("angle", s.spec.e.theta * degrad, s.spec.p.theta * degrad, "deg"));
    o.line(kin2("momentum", s.spec.e.P, s.spec.p.P, "MeV/c"));
    o.line(kin2("x offset", s.spec.e.offset.x, s.spec.p.offset.x, "cm"));
    o.line(kin2("y offset", s.spec.e.offset.y, s.spec.p.offset.y, "cm"));
    o.line(kin2("z offset", s.spec.e.offset.z, s.spec.p.offset.z, "cm"));
    // FORTRAN-QUIRK (D-018): printed in rad under the label 'mr'
    o.line(kin2("xptar offset", s.spec.e.offset.xptar, s.spec.p.offset.xptar, "mr"));
    o.line(kin2("yptar offset", s.spec.e.offset.yptar, s.spec.p.offset.yptar, "mr"));

    o.list("                      VALUES FOR \"CENTRAL\" EVENT:");
    o.line(kin2("delta", central.e.delta, central.p.delta, "%"));
    o.line(kin2("xptar", central.e.xptar, central.p.xptar, "mr"));
    o.line(kin2("yptar", central.e.yptar, central.p.yptar, "mr"));
    o.line(cen1("Q2", central.Q2 / 1.e6, "(GeV/c)^2"));
    o.line(cen1("q", central.q, "MeV/c"));
    o.line(cen1("nu", central.nu, "MeV"));
    o.line(cen1("recon Em", central.Em, "MeV"));
    o.line(cen1("recon Pm", central.Pm, "MeV/c"));
    o.line(cen1("recon  W", central.W, "MeV/c"));
    o.line(cen1("recon MM", central.MM, "MeV/c"));

    // Target (l.837-875)
    const TargetInfo& tg = s.targ;
    o.list("TARGET specs:");
    o.line(targ2("A", tg.A, " ", "Z", tg.Z, " "));
    o.line(targ2("mass", tg.mass_amu, "amu", "mass", tg.M, "MeV"));
    o.line(targ2("Mrec", tg.mrec_amu, "amu", "Mrec", tg.Mrec, "MeV"));
    o.line(targ2("Mtar_struc", tg.Mtar_struck, "MeV", "Mrec_struc", tg.Mrec_struck, "MeV"));
    o.line(targ2("rho", tg.rho, "g/cm3", "thick", tg.thick, "g/cm2"));
    o.line(targ2("angle", tg.angle * degrad, "deg", "abundancy", tg.abundancy, "%"));
    o.line(targ2("X0", tg.X0, "g/cm2", "X0_cm", tg.X0_cm, "cm"));
    o.line(targ2("length", tg.length, "cm", "zoffset", tg.zoffset, "cm"));
    o.line(targ2("xoffset", tg.xoffset, "cm", "yoffset", tg.yoffset, "cm"));
    { Record x; x.t(12).a("__ave__", 15).a("__lo__", 15).a("__hi__", 15); o.line(x); }
    o.line(targ3("Coulomb", tg.Coulomb, "MeV"));
    o.line(targ3("Eloss_beam", tg.Eloss[0], "MeV"));
    o.line(targ3("Eloss_e", tg.Eloss[1], "MeV"));
    o.line(targ3("Eloss_p", tg.Eloss[2], "MeV"));
    o.line(targ3("teff_beam", tg.teff[0], "radlen"));
    o.line(targ3("teff_e", tg.teff[1], "radlen"));
    o.line(targ3("teff_p", tg.teff[2], "radlen"));
    o.line(targ1("musc_nsig_max", tg.musc_nsig_max, " "));
    o.line(targ1("musc_max_beam", tg.musc_max[0] * units::kMradPerRad, "mr"));
    o.line(targ1("musc_max_e", tg.musc_max[1] * units::kMradPerRad, "mr"));
    o.line(targ1("musc_max_p", tg.musc_max[2] * units::kMradPerRad, "mr"));

    // Flags (l.878-913). Values as the Fortran COMMON holds them after dbase.f.
    o.list("FLAGS:");
    o.line(flags({{"doing_eep", f.doing_eep}, {"doing_kaon", f.doing_kaon}, {"doing_pion", f.doing_pion}}));
    o.line(flags({{"doing_semi", f.doing_semi}, {"doing_rho", f.doing_rho}, {"doing_hplus", f.doing_hplus}}));
    o.line(flags({{"doing_semipi", f.doing_semipi}, {"doing_semika", f.doing_semika}, {"doing_pizero", f.doing_pizero}}));
    o.line(flags({{"doing_delta", f.doing_delta}, {"doing_phsp", f.doing_phsp}}));
    {
        Record x;
        x.x(5).x(2).a("which_pion", 19).lit("=").i(p.which_pion, 2);
        x.x(2).a("which_kaon", 19).lit("=").i(p.which_kaon, 2);
        x.x(2).a("pizero_ngamma", 19).lit("=").i(p.pizero_ngamma, 2);
        o.line(x);
    }
    o.line(flags({{"doing_hyd_elast", f.doing_hyd_elast}, {"doing_deuterium", f.doing_deuterium}, {"doing_heavy", f.doing_heavy}}));
    o.line(flags({{"doing_hydpi", f.doing_hydpi}, {"doing_deutpi", f.doing_deutpi}, {"doing_hepi", f.doing_hepi}}));
    o.line(flags({{"doing_hydkaon", f.doing_hydkaon}, {"doing_deutkaon", f.doing_deutkaon}, {"doing_hekaon", f.doing_hekaon}}));
    o.line(flags({{"doing_hydsemi", f.doing_hydsemi}, {"doing_deutsemi", f.doing_deutsemi}, {"do_fermi", f.do_fermi}}));
    o.line(flags({{"doing_hydrho", f.doing_hydrho}, {"doing_deutrho", f.doing_deutrho}, {"doing_herho", f.doing_herho}}));
    {
        Record x;
        x.x(5).x(2).a("mc_smear", 19).lit("=").l(s.sw.mc_smear, 2);
        x.x(2).a("electron_arm", 19).lit("=").i(p.electron_arm, 2);
        x.x(2).a("hadron_arm", 19).lit("=").i(p.hadron_arm, 2);
        o.line(x);
    }
    // deForest_flag is an integer printed with L (gfortran: T if non-zero)
    o.line(flags({{"using_Eloss", s.sw.using_Eloss}, {"using_Coulomb", s.sw.using_Coulomb}, {"deForest_flag", p.deForest_flag != 0}}));
    o.line(flags({{"correct_Eloss", s.sw.correct_Eloss}, {"correct_raster", p.correct_raster}, {"doing_decay", p.doing_decay}}));
    o.line(flags({{"using_E_arm_montecarlo", s.sw.using_E_arm_montecarlo},
                  {"using_P_arm_montecarlo", s.sw.using_P_arm_montecarlo},
                  {"use_benhar_sf", p.use_benhar_sf}}));
    // simc.f:907-909 (the Fortran also lists the retired code 6; rejected at input, D-032)
    if (Is(p.electron_arm, SpectrometerId::SHMS) || Is(p.hadron_arm, SpectrometerId::SHMS)) {
        Record x;
        x.x(7).a("use_first_cer", 19).lit("=").l(p.use_first_cer != 0, 2);
        o.line(x);
    }
    { Record x; x.x(7).a("ctau", 11).lit("=").f(p.ctau, 10, 3).a("cm", 4); o.line(x); }
    if (p.use_benhar_sf) { Record x; x.x(7).a("transparency", 12).lit("=").f(p.transparency, 8, 4); o.line(x); }

    // Counters (l.915-921)
    o.list("COUNTERS:");
    { Record x; x.x(12).lit("Ngen (request) = ").i(r.ngen, 10); o.line(x); }
    { Record x; x.x(12).lit("Ntried         = ").i(t.ntried, 10); o.line(x); }
    { Record x; x.x(12).lit("Ncontribute    = ").i(t.ncontribute, 10); o.line(x); }
    { Record x; x.x(12).lit("Nco_no_rad_prot= ").i(t.ncontribute_no_rad_proton, 10); o.line(x); }
    {
        const double pct = (100. * static_cast<double>(t.ncontribute_no_rad_proton)) /
                           fmaxF(static_cast<double>(t.ncontribute), 0.1e0);
        Record x;
        x.x(12).lit("-> %no_rad_prot= ").f(pct, 10, 3);
        o.line(x);
    }
    o.line("");
    { Record x; x.x(1).lit("INTEGRATED WEIGHTS (number of counts in delta/Em cuts!):"); o.line(x); }
    {
        Record x;
        x.x(1).lit("              MeV: wtcontr= ").e(n.wtcontribute / static_cast<double>(t.nevent), 16, 8);
        o.line(x);
    }

    // Radiative corrections (l.924-951)
    o.list("RADIATIVE CORRECTIONS:");
    if (!s.sw.using_rad) {
        Record x;
        x.x(1).a("using_rad", 14).lit("=").l(s.sw.using_rad, 3);
        o.line(x);
    } else {
        {
            Record x;
            x.x(1).a("use_expon", 14).lit("=").l(p.use_expon != 0, 3);
            x.x(1).a("include_hard", 14).lit("=").l(s.sw.include_hard, 3);
            x.x(1).a("calc_spence", 14).lit("=").l(s.sw.calculate_spence, 3);
            o.line(x);
        }
        {
            Record x;
            x.x(1).a("using_rad", 14).lit("=").l(s.sw.using_rad, 3);
            x.x(1).a("use_offshell_rad", 14).lit("=").l(p.use_offshell_rad, 3);
            o.line(x);
        }
        {
            Record x;
            x.x(1).a("rad_flag", 14).lit("=").i(s.sw.rad_flag, 3);
            x.x(1).a("extrad_flag", 14).lit("=").i(s.sw.extrad_flag, 3);
            x.x(1).a("one_tail", 14).lit("=").i(p.one_tail, 3);
            x.x(1).a("intcor_mode", 14).lit("=").i(p.intcor_mode, 3);
            o.line(x);
        }
        { Record x; x.x(1).a("dE_edge_test", 14).lit("=").f(s.params.dE_edge_test, 11, 3); o.line(x); }
        { Record x; x.x(1).a("Egamma_max", 14).lit("=").f(s.rad.Egamma_tot_max, 11, 3); o.line(x); }
        const CentralRad& cr = central.rad;
        o.list("Central Values:");
        o.line(rad4("hardcorfac", {cr.hardcorfac}));
        o.line(rad4("etatzai", {cr.etatzai}));
        o.line(rad4("frac(1:3)", {cr.frac[0], cr.frac[1], cr.frac[2]}));
        o.line(rad4("lambda(1:3)", {cr.lambda[0], cr.lambda[1], cr.lambda[2]}));
        o.line(rad4("bt(1:2)", {cr.bt[0], cr.bt[1]}));
        o.line(rad4("c_int(0:3)", {cr.c_int[0], cr.c_int[1], cr.c_int[2], cr.c_int[3]}));
        o.line(rad4("c_ext(0:3)", {cr.c_ext[0], cr.c_ext[1], cr.c_ext[2], cr.c_ext[3]}));
        o.line(rad4("c(0:3)", {cr.c[0], cr.c[1], cr.c[2], cr.c[3]}));
        o.line(rad4("g_int", {cr.g_int}));
        o.line(rad4("g_ext", {cr.g_ext}));
        o.line(rad4("g(0:3)", {cr.g[0], cr.g[1], cr.g[2], cr.g[3]}));
    }

    // Miscellaneous (l.954-973)
    o.list("MISCELLANEOUS:");
    o.list("Note that central.sigcc is for central delta,theta,phi in both spectrometers");
    o.list(" and may give non-physical kinematics, esp. for Hydrogen");
    o.list("Note also that AVE.sigcc is really AVER.weight (the two arenot exactly equal)");
    o.line(misc("CENTRAL.sigcc", central.sigcc, " "));
    o.line(misc("AVERAGE.sigcc", t.sum_sigcc / static_cast<double>(t.nevent), " "));
    o.line(misc("charge", n.charge, "mC"));
    o.line(misc("targetfac", n.targetfac, " "));
    o.line(misc("luminosity", n.luminosity, "ub^-1"));
    {
        const double h = constants::hbarc / 100000.;
        o.line(misc("luminosity", n.luminosity * (h * h), "GeV^2"));
    }
    o.line(misc("genvol", n.genvol, " "));
    o.line(misc("normfac", n.normfac, " "));
    if (f.doing_heavy) {
        Record x;
        x.x(12).lit("Theory file:  ").lit(s.theory.file);
        o.line(x);
    }
    if (r.random_seed != 0) { Record x; x.x(15).lit("Random Seed = ").i(r.random_seed, 10); o.line(x); }
    if (!r.random_state_file.empty() && r.random_state_file.find_first_not_of(' ') != std::string::npos) {
        Record x;
        x.lit(" Random State Save File:  ").lit(r.random_state_file.substr(0, r.random_state_file.find(' ')));
        o.line(x);
    }

    // Resolution summary (l.977-986)
    o.list("RECON SUMMARY:            Ave.Error   Resolution");
    const TwoArmSums& ae = n.aveerr;
    const TwoArmSums& rs = n.resol;
    o.line(recon("Electron arm: delta =", units::kPermillePerPercent * ae.e.delta, units::kPermillePerPercent * rs.e.delta, "x10^-3"));
    o.line(recon("xptar =", units::kMradPerRad * ae.e.xptar, units::kMradPerRad * rs.e.xptar, "mr"));
    o.line(recon("yptar =", units::kMradPerRad * ae.e.yptar, units::kMradPerRad * rs.e.yptar, "mr"));
    o.line(recon("ytar  =", units::kMmPerCm * ae.e.ytar, units::kMmPerCm * rs.e.ytar, "mm"));
    o.line(recon("Hadron arm:   delta =", units::kPermillePerPercent * ae.p.delta, units::kPermillePerPercent * rs.p.delta, "x10-3"));
    o.line(recon("xptar =", units::kMradPerRad * ae.p.xptar, units::kMradPerRad * rs.p.xptar, "mr"));
    o.line(recon("yptar =", units::kMradPerRad * ae.p.yptar, units::kMradPerRad * rs.p.yptar, "mr"));
    o.line(recon("ytar  =", units::kMmPerCm * ae.p.ytar, units::kMmPerCm * rs.p.ytar, "mm"));

    // Input spectrometer limits (l.991-1003)
    o.list("Input Spectrometer Limits:");
    o.line(spedge("SPedge.e.delta.min/max", s.SPedge.e.delta, "%"));
    o.line(spedge("SPedge.e.yptar.min/max", s.SPedge.e.yptar, "rad"));
    o.line(spedge("SPedge.e.xptar.min/max", s.SPedge.e.xptar, "rad"));
    o.line(spedge("SPedge.p.delta.min/max", s.SPedge.p.delta, "%"));
    o.line(spedge("SPedge.p.yptar.min/max", s.SPedge.p.yptar, "rad"));
    o.line(spedge("SPedge.p.xptar.min/max", s.SPedge.p.xptar, "rad"));

    // Generated quantities (l.1010-1050)
    const GenLimits& g = s.gen;
    o.list("Limiting VERTEX values (vertex.e/p.*,Em,Pm,Trec)");
    o.list("   USED limits are gen.e/p.*, and VERTEXedge.Em,Pm,Trec");
    o.line(limHeader1());
    o.line(limHeader2());
    o.line(lim("E arm  delta", g.e.delta.min, g.e.delta.max, c.gen.e.delta.lo, c.gen.e.delta.hi, "%"));
    o.line(lim("E arm  yptar", g.e.yptar.min * 1000., g.e.yptar.max * 1000., c.gen.e.yptar.lo * 1000., c.gen.e.yptar.hi * 1000., "mr"));
    o.line(lim("E arm  xptar", g.e.xptar.min * 1000., g.e.xptar.max * 1000., c.gen.e.xptar.lo * 1000., c.gen.e.xptar.hi * 1000., "mr"));
    o.line(lim("P arm  delta", g.p.delta.min, g.p.delta.max, c.gen.p.delta.lo, c.gen.p.delta.hi, "%"));
    o.line(lim("P arm  yptar", g.p.yptar.min * 1000., g.p.yptar.max * 1000., c.gen.p.yptar.lo * 1000., c.gen.p.yptar.hi * 1000., "mr"));
    o.line(lim("P arm  xptar", g.p.xptar.min * 1000., g.p.xptar.max * 1000., c.gen.p.xptar.lo * 1000., c.gen.p.xptar.hi * 1000., "mr"));
    o.line(lim("sumEgen", g.sumEgen.min, g.sumEgen.max, c.gen.sumEgen.lo, c.gen.sumEgen.hi, "MeV"));
    const EdgeTrue& V = s.VERTEXedge;
    o.line(lim("Trec", V.Trec.min, V.Trec.max, c.vertex.Trec.lo, c.vertex.Trec.hi, "MeV"));
    o.line(lim("Em", V.Em.min, V.Em.max, c.vertex.Em.lo, c.vertex.Em.hi, "MeV"));
    o.line(lim("Pm", V.Pm.min, V.Pm.max, c.vertex.Pm.lo, c.vertex.Pm.hi, "MeV/c"));
    if ((f.doing_deuterium || f.doing_pion || f.doing_kaon || f.doing_delta) && s.sw.using_rad)
        o.list("      *** NOTE: sumEgen.min only used in GENERATE_RAD");

    // True quantities (l.1053-1075)
    const EdgeTrue& E = s.edge;
    o.list("Limiting ORIGINAL values: orig.e/p.*,Em,Pm,Trec (no edge.* limits for Pm,Trec)");
    o.line(limHeader1());
    o.line(limHeader2());
    o.line(lim("E arm   E", E.e.E.min, E.e.E.max, c.tru.e.E.lo, c.tru.e.E.hi, "MeV"));
    o.line(lim("E arm  yptar", E.e.yptar.min * 1000., E.e.yptar.max * 1000., c.tru.e.yptar.lo * 1000., c.tru.e.yptar.hi * 1000., "mr"));
    o.line(lim("E arm  xptar", E.e.xptar.min * 1000., E.e.xptar.max * 1000., c.tru.e.xptar.lo * 1000., c.tru.e.xptar.hi * 1000., "mr"));
    o.line(lim("P arm      E", E.p.E.min, E.p.E.max, c.tru.p.E.lo, c.tru.p.E.hi, "MeV"));
    o.line(lim("P arm  yptar", E.p.yptar.min * 1000., E.p.yptar.max * 1000., c.tru.p.yptar.lo * 1000., c.tru.p.yptar.hi * 1000., "mr"));
    o.line(lim("P arm  xptar", E.p.xptar.min * 1000., E.p.xptar.max * 1000., c.tru.p.xptar.lo * 1000., c.tru.p.xptar.hi * 1000., "mr"));
    o.line(lim("Em", fmaxF(-999999.999e0, E.Em.min), fminF(999999.999e0, E.Em.max), c.tru.Em.lo, c.tru.Em.hi, "MeV"));
    o.line(lim("Pm", 0., 0., c.tru.Pm.lo, c.tru.Pm.hi, "MeV"));
    o.line(lim("Trec", 0., 0., c.tru.Trec.lo, c.tru.Trec.hi, "MeV"));

    // Spectrometer quantities (l.1078-1097)
    o.list("Limiting SPECTROMETER values");
    o.line(limHeader1());
    o.line(limHeader2());
    const auto& SE = s.SPedge;
    o.line(lim("E arm delta", SE.e.delta.min, SE.e.delta.max, c.SP.e.delta.lo, c.SP.e.delta.hi, "%"));
    o.line(lim("E arm  yptar", SE.e.yptar.min * 1000., SE.e.yptar.max * 1000., c.SP.e.yptar.lo * 1000., c.SP.e.yptar.hi * 1000., "mr"));
    o.line(lim("E arm  xptar", SE.e.xptar.min * 1000., SE.e.xptar.max * 1000., c.SP.e.xptar.lo * 1000., c.SP.e.xptar.hi * 1000., "mr"));
    o.line(lim("P arm  delta", SE.p.delta.min, SE.p.delta.max, c.SP.p.delta.lo, c.SP.p.delta.hi, "%"));
    o.line(lim("P arm  yptar", SE.p.yptar.min * 1000., SE.p.yptar.max * 1000., c.SP.p.yptar.lo * 1000., c.SP.p.yptar.hi * 1000., "mr"));
    o.line(lim("P arm  xptar", SE.p.xptar.min * 1000., SE.p.xptar.max * 1000., c.SP.p.xptar.lo * 1000., c.SP.p.xptar.hi * 1000., "mr"));

    // Radiation (l.1100-1113)
    if (s.sw.using_rad) {
        o.list("Limiting RADIATION values CONTRIBUTING to the (Em,Pm) distributions:");
        o.line(limHeader1());
        o.line(limHeader2());
        o.line(lim("Egamma(1)", 0., s.rad.Egamma1_max, c.rad.Egamma[0].lo, c.rad.Egamma[0].hi, "MeV"));
        o.line(lim("Egamma(2)", 0., s.rad.Egamma2_max, c.rad.Egamma[1].lo, c.rad.Egamma[1].hi, "MeV"));
        o.line(lim("Egamma(3)", 0., s.rad.Egamma3_max, c.rad.Egamma[2].lo, c.rad.Egamma[2].hi, "MeV"));
        o.line(lim("Egamma_total", 0., s.rad.Egamma_tot_max, c.rad.Egamma_total.lo, c.rad.Egamma_total.hi, "MeV"));
    }

    // Slops (l.1116-1137)
    const Slops& sl = t.slop;
    o.list("ACTUAL and LIMITING SLOP values used/obtained:");
    { Record x; x.t(25).a("__used__", 10).a("__min__", 10).a("__max__", 10); o.line(x); }
    o.line(slopLine("slop.MC  ", "E arm delta", sl.MC.e.delta.used, sl.MC.e.delta.lo, sl.MC.e.delta.hi, "%"));
    o.line(slopLine(" ", "E arm yptar", sl.MC.e.yptar.used * 1000., sl.MC.e.yptar.lo * 1000., sl.MC.e.yptar.hi * 1000., "mr"));
    o.line(slopLine(" ", "E arm xptar", sl.MC.e.xptar.used * 1000., sl.MC.e.xptar.lo * 1000., sl.MC.e.xptar.hi * 1000., "mr"));
    o.line(slopLine(" ", "P arm delta", sl.MC.p.delta.used, sl.MC.p.delta.lo, sl.MC.p.delta.hi, "%"));
    o.line(slopLine(" ", "P arm yptar", sl.MC.p.yptar.used * 1000., sl.MC.p.yptar.lo * 1000., sl.MC.p.yptar.hi * 1000., "mr"));
    o.line(slopLine(" ", "P arm xptar", sl.MC.p.xptar.used * 1000., sl.MC.p.xptar.lo * 1000., sl.MC.p.xptar.hi * 1000., "mr"));
    o.line(slopLine("slop.total", "Em", sl.total.Em.used, sl.total.Em.lo, sl.total.Em.hi, "MeV"));
    o.line(slopLine(" ", "Pm", 0., sl.total.Pm.lo, sl.total.Pm.hi, "MeV/c"));

    // write(iun,'(/)')
    o.line("");
    o.line("");
    return o.str();
}

namespace {
long cnt(const RunReport& r, const std::string& spec, const std::string& key) {
    for (const auto& sc : r.stops)
        if (sc.spectrometer == spec) {
            auto it = sc.counts.find(key);
            if (it == sc.counts.end()) throw std::runtime_error("FormatGeni: no counter '" + key + "' for " + spec);
            return it->second;
        }
    throw std::runtime_error("FormatGeni: no stop counters for " + spec);
}
bool has(const RunReport& r, const std::string& spec) {
    for (const auto& sc : r.stops)
        if (sc.spectrometer == spec) return true;
    return false;
}
}  // namespace

// Fortran: simc.f l.462-532
std::string FormatGeni(const RunSetup& s, const RunReport& r) {
    Out o;
    const int ea = s.params.electron_arm, ha = s.params.hadron_arm;
    auto row = [&](const std::string& label, const std::string& spec, std::initializer_list<const char*> keys) {
        Record x = List(label);
        for (const char* k : keys) x.listInt(cnt(r, spec, k));
        o.line(x);
    };
    // The Fortran blocks (by arm code) are written for the Fortran spectrometers; an arm served by another ISpectrometer
    // plug-in (e.g. examples/ideal_box) gets a generic block at the end (P7.4).
    auto uses = [&](SpectrometerId id) { return Is(ea, id) || Is(ha, id); };   // D-032 spectrometer table
    if (uses(SpectrometerId::HMS) && has(r, "hms")) {
        row("HMS Trials:           ", "hms", {"trials"});
        row("Slit hor/vert/corners ", "hms", {"slit_hor", "slit_vert", "slit_oct"});
        row("Q1 entrance/mid/exit  ", "hms", {"Q1_in", "Q1_mid", "Q1_out"});
        row("Q2 entrance/mid/exit  ", "hms", {"Q2_in", "Q2_mid", "Q2_out"});
        row("Q3 entrance/mid/exit  ", "hms", {"Q3_in", "Q3_mid", "Q3_out"});
        row("Dipole entrance/exit  ", "hms", {"D1_in", "D1_out"});
        row("Events reaching hut   ", "hms", {"hut"});
        row("DC1, DC2, Scin, Cal   ", "hms", {"dc1", "dc2", "scin", "cal"});
        row("Successes             ", "hms", {"successes"});
        o.line("");
    }
    if (uses(SpectrometerId::SHMS) && has(r, "shms")) {   // simc.f:510-511 also lists the retired code 6
        row("SHMS Trials:          ", "shms", {"trials"});
        row("HB phys entrance/mag entr/mag exit/phys exit  ", "shms", {"hb_in", "hb_men", "hb_mex", "hb_out"});
        row("Slit hor/vert/corners ", "shms", {"slit_hor", "slit_vert", "slit_oct"});
        row("Q1 phys entrance/mag entr/mid/mag exit/phys exit  ", "shms", {"q1_in", "q1_men", "q1_mid", "q1_mex", "q1_out"});
        row("Q2 phys entrance/mag entr/mid/mag exit/phys exit  ", "shms", {"q2_in", "q2_men", "q2_mid", "q2_mex", "q2_out"});
        // FORTRAN-QUIRK simc.f:519-520: the Q3 "mid" column prints shmsSTOP_q2_mid
        row("Q3 phys entrance/mag entr/mid/mag exit/phys exit       ", "shms", {"q3_in", "q3_men", "q2_mid", "q3_mex", "q3_out"});
        row("D1 entrance/flare/mid 1-2   ", "shms", {"D1_in", "D1_flr", "D1_mid1", "D1_mid2"});
        row("D1 mid 3-5            ", "shms", {"D1_mid3", "D1_mid4", "D1_mid5"});
        row("D1 mid 6-7/mag exit/phys exit         ", "shms", {"D1_mid6", "D1_mid7", "D1_mex", "D1_out"});
        row("Events reaching hut   ", "shms", {"hut"});
        // FORTRAN-QUIRK simc.f:526: label with four names, two values
        row("DC1, DC2, Scin, Cal   ", "shms", {"dc1", "dc2"});
        row("S1, S2, S3, Cal       ", "shms", {"s1", "s2", "s3", "cal"});
        row("Successes             ", "shms", {"successes"});
        o.line("");
    }
    // SOS/HRS/calorimeter arms (codes 2-4, 7, 8) have no Fortran-format block yet (P14); without a plug-in they are
    // refused at initialise (MonteCarlo::DefaultArm). Every non-Fortran spectrometer: one "key value" line per counter.
    for (const auto& sc : r.stops) {
        if (sc.spectrometer == "hms" || sc.spectrometer == "shms") continue;
        o.list(sc.spectrometer + " stop counters (" + sc.role + " arm, SIMC++ plug-in):");
        for (const auto& kv : sc.counts) {
            Record x = List(" " + kv.first);
            x.listInt(kv.second);
            o.line(x);
        }
        o.line("");
    }
    return o.str();
}

// ------------------------------------------------------------------------------------------ summary sections
RunSummary MakeRunSummary(const RunSetup& s, const RunReport& r) {
    const RunTotals& t = r.totals;
    const Normalisation& n = r.norm;
    RunSummary sum;
    sum.run_name = r.run_name;
    sum.ngen = r.ngen;
    sum.ntried = t.ntried;
    sum.ncontribute = t.ncontribute;
    sum.ncontribute_no_rad_proton = t.ncontribute_no_rad_proton;
    sum.nevent = t.nevent;
    sum.wtcontribute = n.wtcontribute;
    sum.sum_sigcc = t.sum_sigcc;
    sum.central_sigcc = r.central.sigcc;
    sum.charge = n.charge;
    sum.targetfac = n.targetfac;
    sum.luminosity = n.luminosity;
    sum.genvol = n.genvol;
    sum.normfac = n.normfac;
    sum.random_seed = r.random_seed;
    sum.stops = r.stops;

    Json& x = sum.extra;
    const CentralEvent& c = r.central;
    x["central"] = {{"Q2_GeV2", c.Q2 / 1.e6}, {"q", c.q}, {"nu", c.nu}, {"Em", c.Em}, {"Pm", c.Pm}, {"W", c.W}, {"MM", c.MM},
                    {"e", {{"delta", c.e.delta}, {"xptar_mr", c.e.xptar}, {"yptar_mr", c.e.yptar}}},
                    {"p", {{"delta", c.p.delta}, {"xptar_mr", c.p.xptar}, {"yptar_mr", c.p.yptar}}}};
    const TargetInfo& tg = s.targ;
    auto tf = [](const TargetFunction& v) { return Json::array({v.ave, v.min, v.max}); };
    x["target"] = {{"A", tg.A}, {"Z", tg.Z}, {"mass_amu", tg.mass_amu}, {"mass_MeV", tg.M}, {"Mrec_amu", tg.mrec_amu},
                   {"Mrec_MeV", tg.Mrec}, {"Mtar_struck", tg.Mtar_struck}, {"Mrec_struck", tg.Mrec_struck},
                   {"rho", tg.rho}, {"thick_g_cm2", tg.thick}, {"angle_deg", tg.angle * degrad},
                   {"abundancy", tg.abundancy}, {"X0", tg.X0}, {"X0_cm", tg.X0_cm}, {"length", tg.length},
                   {"zoffset", tg.zoffset}, {"xoffset", tg.xoffset}, {"yoffset", tg.yoffset},
                   {"musc_nsig_max", tg.musc_nsig_max}, {"musc_max_beam_mr", tg.musc_max[0] * units::kMradPerRad},
                   {"musc_max_e_mr", tg.musc_max[1] * units::kMradPerRad}, {"musc_max_p_mr", tg.musc_max[2] * units::kMradPerRad},
                   {"losses", {{"Coulomb", tf(tg.Coulomb)}, {"Eloss_beam", tf(tg.Eloss[0])}, {"Eloss_e", tf(tg.Eloss[1])},
                               {"Eloss_p", tf(tg.Eloss[2])}, {"teff_beam", tf(tg.teff[0])}, {"teff_e", tf(tg.teff[1])},
                               {"teff_p", tf(tg.teff[2])}}}};
    if (s.sw.using_rad) {
        const CentralRad& cr = c.rad;
        x["radiative"] = {{"rad_flag", s.sw.rad_flag}, {"extrad_flag", s.sw.extrad_flag}, {"one_tail", s.params.one_tail},
                          {"intcor_mode", s.params.intcor_mode},
                          {"dE_edge_test", s.params.dE_edge_test},
                          {"Egamma_max", s.rad.Egamma_tot_max},
                          {"central", {{"hardcorfac", cr.hardcorfac}, {"etatzai", cr.etatzai},
                                       {"frac", {cr.frac[0], cr.frac[1], cr.frac[2]}},
                                       {"lambda", {cr.lambda[0], cr.lambda[1], cr.lambda[2]}},
                                       {"bt", {cr.bt[0], cr.bt[1]}},
                                       {"c_int", {cr.c_int[0], cr.c_int[1], cr.c_int[2], cr.c_int[3]}},
                                       {"c_ext", {cr.c_ext[0], cr.c_ext[1], cr.c_ext[2], cr.c_ext[3]}},
                                       {"c", {cr.c[0], cr.c[1], cr.c[2], cr.c[3]}},
                                       {"g_int", cr.g_int}, {"g_ext", cr.g_ext},
                                       {"g", {cr.g[0], cr.g[1], cr.g[2], cr.g[3]}}}}};
    }
    auto rec = [&](double a, double b, const char* u) { return Json{{"ave_error", a}, {"resolution", b}, {"units", u}}; };
    const TwoArmSums& ae = n.aveerr;
    const TwoArmSums& rs = n.resol;
    x["recon_summary"] = {
        {"e", {{"delta", rec(units::kPermillePerPercent * ae.e.delta, units::kPermillePerPercent * rs.e.delta, "x10^-3")}, {"xptar", rec(units::kMradPerRad * ae.e.xptar, units::kMradPerRad * rs.e.xptar, "mr")},
               {"yptar", rec(units::kMradPerRad * ae.e.yptar, units::kMradPerRad * rs.e.yptar, "mr")}, {"ytar", rec(units::kMmPerCm * ae.e.ytar, units::kMmPerCm * rs.e.ytar, "mm")}}},
        {"p", {{"delta", rec(units::kPermillePerPercent * ae.p.delta, units::kPermillePerPercent * rs.p.delta, "x10-3")}, {"xptar", rec(units::kMradPerRad * ae.p.xptar, units::kMradPerRad * rs.p.xptar, "mr")},
               {"yptar", rec(units::kMradPerRad * ae.p.yptar, units::kMradPerRad * rs.p.yptar, "mr")}, {"ytar", rec(units::kMmPerCm * ae.p.ytar, units::kMmPerCm * rs.p.ytar, "mm")}}}};

    // limits tables (labels as hist2json writes them: white space collapsed)
    const reactions::Contrib& cb = t.contrib;
    auto row = [](const char* label, double umin, double umax, double lo, double hi, const char* u) {
        return Json{{"label", label}, {"used", {umin, umax}}, {"found", {lo, hi}}, {"units", u}};
    };
    const GenLimits& g = s.gen;
    const EdgeTrue& V = s.VERTEXedge;
    const EdgeTrue& E = s.edge;
    const auto& SE = s.SPedge;
    x["limits"]["vertex"] = Json::array({
        row("E arm delta", g.e.delta.min, g.e.delta.max, cb.gen.e.delta.lo, cb.gen.e.delta.hi, "%"),
        row("E arm yptar", g.e.yptar.min * 1000., g.e.yptar.max * 1000., cb.gen.e.yptar.lo * 1000., cb.gen.e.yptar.hi * 1000., "mr"),
        row("E arm xptar", g.e.xptar.min * 1000., g.e.xptar.max * 1000., cb.gen.e.xptar.lo * 1000., cb.gen.e.xptar.hi * 1000., "mr"),
        row("P arm delta", g.p.delta.min, g.p.delta.max, cb.gen.p.delta.lo, cb.gen.p.delta.hi, "%"),
        row("P arm yptar", g.p.yptar.min * 1000., g.p.yptar.max * 1000., cb.gen.p.yptar.lo * 1000., cb.gen.p.yptar.hi * 1000., "mr"),
        row("P arm xptar", g.p.xptar.min * 1000., g.p.xptar.max * 1000., cb.gen.p.xptar.lo * 1000., cb.gen.p.xptar.hi * 1000., "mr"),
        row("sumEgen", g.sumEgen.min, g.sumEgen.max, cb.gen.sumEgen.lo, cb.gen.sumEgen.hi, "MeV"),
        row("Trec", V.Trec.min, V.Trec.max, cb.vertex.Trec.lo, cb.vertex.Trec.hi, "MeV"),
        row("Em", V.Em.min, V.Em.max, cb.vertex.Em.lo, cb.vertex.Em.hi, "MeV"),
        row("Pm", V.Pm.min, V.Pm.max, cb.vertex.Pm.lo, cb.vertex.Pm.hi, "MeV/c")});
    x["limits"]["original"] = Json::array({
        row("E arm E", E.e.E.min, E.e.E.max, cb.tru.e.E.lo, cb.tru.e.E.hi, "MeV"),
        row("E arm yptar", E.e.yptar.min * 1000., E.e.yptar.max * 1000., cb.tru.e.yptar.lo * 1000., cb.tru.e.yptar.hi * 1000., "mr"),
        row("E arm xptar", E.e.xptar.min * 1000., E.e.xptar.max * 1000., cb.tru.e.xptar.lo * 1000., cb.tru.e.xptar.hi * 1000., "mr"),
        row("P arm E", E.p.E.min, E.p.E.max, cb.tru.p.E.lo, cb.tru.p.E.hi, "MeV"),
        row("P arm yptar", E.p.yptar.min * 1000., E.p.yptar.max * 1000., cb.tru.p.yptar.lo * 1000., cb.tru.p.yptar.hi * 1000., "mr"),
        row("P arm xptar", E.p.xptar.min * 1000., E.p.xptar.max * 1000., cb.tru.p.xptar.lo * 1000., cb.tru.p.xptar.hi * 1000., "mr"),
        row("Em", fmaxF(-999999.999e0, E.Em.min), fminF(999999.999e0, E.Em.max), cb.tru.Em.lo, cb.tru.Em.hi, "MeV"),
        row("Pm", 0., 0., cb.tru.Pm.lo, cb.tru.Pm.hi, "MeV"),
        row("Trec", 0., 0., cb.tru.Trec.lo, cb.tru.Trec.hi, "MeV")});
    x["limits"]["spectrometer"] = Json::array({
        row("E arm delta", SE.e.delta.min, SE.e.delta.max, cb.SP.e.delta.lo, cb.SP.e.delta.hi, "%"),
        row("E arm yptar", SE.e.yptar.min * 1000., SE.e.yptar.max * 1000., cb.SP.e.yptar.lo * 1000., cb.SP.e.yptar.hi * 1000., "mr"),
        row("E arm xptar", SE.e.xptar.min * 1000., SE.e.xptar.max * 1000., cb.SP.e.xptar.lo * 1000., cb.SP.e.xptar.hi * 1000., "mr"),
        row("P arm delta", SE.p.delta.min, SE.p.delta.max, cb.SP.p.delta.lo, cb.SP.p.delta.hi, "%"),
        row("P arm yptar", SE.p.yptar.min * 1000., SE.p.yptar.max * 1000., cb.SP.p.yptar.lo * 1000., cb.SP.p.yptar.hi * 1000., "mr"),
        row("P arm xptar", SE.p.xptar.min * 1000., SE.p.xptar.max * 1000., cb.SP.p.xptar.lo * 1000., cb.SP.p.xptar.hi * 1000., "mr")});
    if (s.sw.using_rad) {
        x["limits"]["radiation"] = Json::array({
            row("Egamma(1)", 0., s.rad.Egamma1_max, cb.rad.Egamma[0].lo, cb.rad.Egamma[0].hi, "MeV"),
            row("Egamma(2)", 0., s.rad.Egamma2_max, cb.rad.Egamma[1].lo, cb.rad.Egamma[1].hi, "MeV"),
            row("Egamma(3)", 0., s.rad.Egamma3_max, cb.rad.Egamma[2].lo, cb.rad.Egamma[2].hi, "MeV"),
            row("Egamma_total", 0., s.rad.Egamma_tot_max, cb.rad.Egamma_total.lo, cb.rad.Egamma_total.hi, "MeV")});
    }
    const Slops& sl = t.slop;
    auto srow = [](const char* label, double used, double lo, double hi, const char* u) {
        return Json{{"label", label}, {"used", used}, {"min", lo}, {"max", hi}, {"units", u}};
    };
    x["slop"]["MC"] = Json::array({
        srow("E arm delta", sl.MC.e.delta.used, sl.MC.e.delta.lo, sl.MC.e.delta.hi, "%"),
        srow("E arm yptar", sl.MC.e.yptar.used * 1000., sl.MC.e.yptar.lo * 1000., sl.MC.e.yptar.hi * 1000., "mr"),
        srow("E arm xptar", sl.MC.e.xptar.used * 1000., sl.MC.e.xptar.lo * 1000., sl.MC.e.xptar.hi * 1000., "mr"),
        srow("P arm delta", sl.MC.p.delta.used, sl.MC.p.delta.lo, sl.MC.p.delta.hi, "%"),
        srow("P arm yptar", sl.MC.p.yptar.used * 1000., sl.MC.p.yptar.lo * 1000., sl.MC.p.yptar.hi * 1000., "mr"),
        srow("P arm xptar", sl.MC.p.xptar.used * 1000., sl.MC.p.xptar.lo * 1000., sl.MC.p.xptar.hi * 1000., "mr")});
    x["slop"]["total"] = Json::array({srow("Em", sl.total.Em.used, sl.total.Em.lo, sl.total.Em.hi, "MeV"),
                                      srow("Pm", 0., sl.total.Pm.lo, sl.total.Pm.hi, "MeV/c")});
    // kinematics in the internal (converted) units, like the .hist prints them
    auto arm = [](const SpecArm& a) {
        return Json{{"angle_deg", a.theta * degrad}, {"momentum", a.P}, {"offset_x", a.offset.x}, {"offset_y", a.offset.y},
                    {"offset_z", a.offset.z}, {"offset_xptar_mr", a.offset.xptar}, {"offset_yptar_mr", a.offset.yptar}};
    };
    x["kinematics"] = {{"Ebeam", s.Ebeam}, {"dE_over_E", s.dEbeam / s.Ebeam}, {"x_width", s.gen.xwid},
                       {"y_width", s.gen.ywid}, {"fr_pattern", s.targ.fr_pattern}, {"fr1", s.targ.fr1},
                       {"fr2", s.targ.fr2}, {"e", arm(s.spec.e)}, {"p", arm(s.spec.p)}};
    // (the extra sections replace BuildSummaryJson's sections of the same name: give them complete)
    x["counters"] = {{"ngen", r.ngen},
                     {"ntried", t.ntried},
                     {"ncontribute", t.ncontribute},
                     {"nco_no_rad_prot", t.ncontribute_no_rad_proton},
                     {"pct_no_rad_proton", (100. * static_cast<double>(t.ncontribute_no_rad_proton)) /
                                               fmaxF(static_cast<double>(t.ncontribute), 0.1e0)}};
    x["run"] = {{"tag", r.run_name},
                {"producer", "simc++"},
                {"random_seed", r.random_seed},
                {"reaction", ReactionTitle(s.flags, s.params)},
                {"begin_time", r.begin_time},
                {"end_time", r.end_time}};
    if (s.flags.doing_heavy) x["run"]["theory_file"] = s.theory.file;
    return sum;
}

}  // namespace driver
}  // namespace simc
