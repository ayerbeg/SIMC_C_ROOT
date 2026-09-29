// src/input/EffectiveValues.cpp
// Fortran rules reproduced here (all from dbase.f:dbase_read unless stated):
//   convert_to_logical (l.1129-1166), random_seed -1 -> time() (l.119),
//   doing_phsp => rad_flag=0 and reaction flags off (l.227-234),
//   dEbeam -> MeV (l.238), spec theta abs()/degrad (l.239-244),
//   targ%thick /1000 (l.451), targ%angle /degrad (l.453), Bangle/Bphi (l.454-455),
//   cryotarget angle forced to 0 (l.456-460), offsets xptar/yptar /1000 (l.472-475),
//   SPedge delta%min clamp and angles /1000 (l.479-488),
//   cuts%Em opened to +-1e6 (l.534-538), correct_Eloss forced off (l.542-545),
//   using_Coulomb off for Z=1 (l.547), extrad_flag default (init.f:radc_init l.621-626).
#include "simc/input/EffectiveValues.h"

#include <cmath>
#include <cstdio>
#include <sstream>

namespace simc {

namespace {

constexpr double kPi = 3.141592653589793;     // constants.inc
constexpr double kDegrad = 180.0 / kPi;       // constants.inc degrad

std::string num(double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.10g", v);
    return buf;
}
std::string tf(bool b) { return b ? "T" : "F"; }

}  // namespace

std::vector<EffectiveRow> EffectiveValues(const SimcParameters& in) {
    // Work on a copy that we transform exactly like dbase.f does.
    SimcParameters p = in;
    std::vector<EffectiveRow> rows;
    auto& table = ParameterTable();

    // Cross-variable rules first (they change what later rows show).
    std::string phsp_rule;
    if (p.doing_phsp) {
        p.rad_flag = 0;
        p.doing_pion = p.doing_kaon = p.doing_delta = p.doing_rho = false;
        phsp_rule = "doing_phsp: dbase.f:227-234";
    }
    const bool cryo = p.targ.Z < 2.4;
    const bool em_open = p.doing_pion || p.doing_kaon || p.doing_delta || (p.cuts.Em.min == p.cuts.Em.max);
    const bool fix_correct_eloss = p.correct_Eloss && !p.using_Eloss;
    const bool coulomb_off = std::lround(p.targ.Z) == 1 && p.using_Coulomb;
    int extrad = p.extrad_flag;
    std::string extrad_rule;
    if (extrad == 0) {
        if (p.rad_flag == 0) extrad = 3;
        else if (p.rad_flag >= 1 && p.rad_flag <= 3) extrad = 1;
        extrad_rule = "init.f:621-626 (0 => 3 if rad_flag=0, 1 if rad_flag=1..3)";
    }

    for (const auto& d : table) {
        EffectiveRow r;
        r.name = d.name;
        r.raw = FormatValue(d, in);
        r.effective = FormatValue(d, p);
        r.units = d.units;
        const std::string n = d.name;
        if (d.kind == ParamKind::Flag || d.kind == ParamKind::FlagArray) {
            r.rule = "convert_to_logical dbase.f:1129 (int>0 => .true.)";
            r.units = "logical";
        }
        if (n == "random_seed" && p.random_seed == -1) {
            r.effective = "time()";
            r.rule = "dbase.f:119";
        } else if (n == "dEbeam") {
            r.effective = num(p.Ebeam * p.dEbeam / 100.0);
            r.units = "MeV (full width)";
            r.rule = "dbase.f:238 dEbeam=Ebeam*dEbeam/100";
        } else if (n == "spec%e%theta" || n == "spec%p%theta") {
            const double th = (n == "spec%e%theta") ? p.spec.e.theta : p.spec.p.theta;
            r.effective = num(std::abs(th) / kDegrad);
            r.units = "rad";
            r.rule = "dbase.f:239-244 abs()/degrad";
        } else if (n == "targ%thick") {
            r.effective = num(p.targ.thick / 1000.0);
            r.units = "g/cm^2";
            r.rule = "dbase.f:451 /1000";
        } else if (n == "targ%angle") {
            double a = p.targ.angle / kDegrad;
            r.rule = "dbase.f:453 /degrad";
            if (cryo && std::abs(a) > 0.0001) {
                a = 0.0;
                r.rule += "; forced 0 for cryotarget (Z<2.4) dbase.f:456-460";
            }
            r.effective = num(a);
            r.units = "rad";
        } else if (n == "targ_Bangle" || n == "targ_Bphi") {
            r.effective = num((n == "targ_Bangle" ? p.targ_Bangle : p.targ_Bphi) / kDegrad);
            r.units = "rad";
            r.rule = "dbase.f:454-455 /degrad";
        } else if (n.find("offset%xptar") != std::string::npos || n.find("offset%yptar") != std::string::npos) {
            const SimcParameters::Offset& o = (n[5] == 'e') ? p.spec.e.offset : p.spec.p.offset;
            r.effective = num((n.find("xptar") != std::string::npos ? o.xptar : o.yptar) / 1000.0);
            r.units = "rad";
            r.rule = "dbase.f:472-475 /1000";
        } else if (n.rfind("SPedge%", 0) == 0) {
            const SimcParameters::ArmEdge& a = (n[7] == 'e') ? p.SPedge.e : p.SPedge.p;
            const bool is_min = n.compare(n.size() - 3, 3, "min") == 0;
            if (n.find("delta") != std::string::npos) {
                double v = is_min ? a.delta.min : a.delta.max;
                if (is_min && v <= -100.0) {
                    v = -99.99;
                    r.rule = "dbase.f:479-480 clamp to -99.99";
                }
                r.effective = num(v);
            } else {
                const ParamRange& rg = (n.find("yptar") != std::string::npos) ? a.yptar : a.xptar;
                r.effective = num((is_min ? rg.min : rg.max) / 1000.0);
                r.units = "rad";
                r.rule = "dbase.f:481-488 /1000 (later widened by slop in init.f:limits_init)";
            }
        } else if ((n == "cuts%Em%min" || n == "cuts%Em%max") && em_open) {
            r.effective = num(n == "cuts%Em%min" ? -1.0e6 : 1.0e6);
            r.rule = "dbase.f:534-538 opened to +-1e6 (pion/kaon/delta or min=max)";
        } else if (n == "correct_Eloss" && fix_correct_eloss) {
            r.effective = tf(false);
            r.rule = "dbase.f:542-545 forced .false. because using_Eloss is .false.";
        } else if (n == "using_Coulomb" && coulomb_off) {
            r.effective = tf(false);
            r.rule = "dbase.f:547 no Coulomb correction for Z=1";
        } else if (n == "extrad_flag" && !extrad_rule.empty()) {
            r.effective = std::to_string(extrad);
            r.rule = extrad_rule;
        } else if (!phsp_rule.empty() &&
                   (n == "rad_flag" || n == "doing_pion" || n == "doing_kaon" || n == "doing_delta" ||
                    n == "doing_rho")) {
            r.rule = phsp_rule;
        }
        rows.push_back(r);
    }
    return rows;
}

std::string FormatEffectiveTable(const std::vector<EffectiveRow>& rows) {
    std::ostringstream os;
    char buf[512];
    std::snprintf(buf, sizeof buf, "%-22s %-24s %-24s %-18s %s\n", "fortran_name", "input", "effective",
                  "units", "rule");
    os << buf;
    for (const auto& r : rows) {
        std::snprintf(buf, sizeof buf, "%-22s %-24s %-24s %-18s %s\n", r.name.c_str(), r.raw.c_str(),
                      r.effective.c_str(), r.units.c_str(), r.rule.c_str());
        os << buf;
    }
    return os.str();
}

}  // namespace simc
