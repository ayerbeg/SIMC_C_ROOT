// src/driver/RunReportRecords.cpp — part of the report writer (RunReport.h): the reaction lines of
// simc.f:report (l.667-788) and the formatted records (FORMAT statements) used by the .hist writer.
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
#include "RunReportInternal.h"

namespace simc {
namespace driver {
namespace report_detail {

using constants::degrad;

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
        // simc.f:677-720
        const char* tg = A == 1 ? "H" : (A == 2 ? "D" : "A");
        const bool known = A == 1 || A == 2 || A >= 3;
        auto title = [&](const char* h) {
            o.list(std::string(" ****--------  ") + tg + "(e,e'" + h + ")X  --------****");
        };
        if (f.doing_semipi) {
            if (f.doing_pizero && !s.fortran_compat && known)
                // FORTRAN-BUG(D-051 F12): original Fortran simc.f:679-700 has no pi0 line: a pi0 run is titled pi+ or pi-
                // by doing_hplus; corrected because dbase.f:701-753 announces the run as (e,e'pi0)X.
                title("pi0");
            else if (known)
                title(f.doing_hplus ? "pi+" : "pi-");
            else stop("I don't have ANY idea what A(e,e'pi)X we're doing!!!");
        } else if (f.doing_semika) {
            if (A == 1 || A == 2) title(f.doing_hplus ? "k+" : "k-");
            else if (A >= 3 && !s.fortran_compat)
                // FORTRAN-BUG(D-051 F14): original Fortran simc.f:701-717 has no A >= 3 branch: a kaon SIDIS run on a nucleus
                // STOPs ('I don't have ANY idea what A(e,e'k)X ...') after the whole event loop, when the .hist is written
                // (the .hist is truncated); corrected because the run itself is well defined (dbase.f:768-773 announces it,
                // peepiX handles A > 2): the A(e,e'k+-)X title is written and the report completes.
                title(f.doing_hplus ? "k+" : "k-");
            else stop("I don't have ANY idea what A(e,e'k)X we're doing!!!");
        } else {
            stop("I don't have ANY idea what A(e,e'x)X we're doing!!!");
        }
    } else if (f.doing_rho) {
        if (A == 1) o.list("              ****--------  H(e,e'rho)  --------****");
        else if (!s.fortran_compat && f.doing_deutrho)
            // FORTRAN-BUG(D-057 F15): original Fortran simc.f:721-725: `write(iun,*) 'I am not set up for anything else
            // yet!'` for every target but H; corrected because D and 3He rho runs are set up and run (dbase.f:867-871
            // announces them): the reaction title is written.
            o.list("              ****--------  D(e,e'rho)  --------****");
        else if (!s.fortran_compat && f.doing_herho)
            o.list("              ****--------  A(e,e'rho)  --------****");
        else o.list("I am not set up for anything else yet!");
    } else if (f.doing_delta) {
        // FORTRAN-BUG(D-018): original Fortran simc.f:727-732: `write(6,*) ' ****--------  H(e,e''p)pi  --------****'`
        // (and D, A); corrected because every other reaction title of report() goes to the report unit iun (the .hist),
        // so the (e,e'p)pi title was missing from the .hist. fortran_compat keeps unit 6.
        auto title = [&](const char* t) {
            if (s.fortran_compat) std::cout << " " << t << "\n";
            else o.list(t);
        };
        if (f.doing_hyddelta) title(" ****--------  H(e,e'p)pi  --------****");
        else if (f.doing_deutdelta) title(" ****--------  D(e,e'p)pi  --------****");
        else if (f.doing_hedelta) title(" ****--------  A(e,e'p)pi  --------****");
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
        // FORTRAN-BUG(D-038 F15): original Fortran simc.f:752: `else if (which_pion .eq. 10) then` prints the final-state
        // line for coherent pi+ only: for which_pion = 11 (coherent pi-) the .hist has no final-state line; corrected
        // because 11 is the same A + pi final state. fortran_compat prints nothing for 11.
        else if (p.which_pion == 10 || (p.which_pion == 11 && !s.fortran_compat))
            o.list("              ****----  Final State is A + pi ----****");
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


}  // namespace report_detail
}  // namespace driver
}  // namespace simc
