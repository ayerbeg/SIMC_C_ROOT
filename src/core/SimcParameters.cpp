// src/core/SimcParameters.cpp
// Parameter table and JSON (de)serialisation for SimcParameters.
// Fortran: dbase.f:regallvars (l.967-1128) — the table below lists the 106 registered
// variables in the same order and with the same CTP names.
#include "simc/core/SimcParameters.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <sstream>

namespace simc {

namespace {

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

// Accessor helpers. Each macro expands to a captureless lambda converted to a function pointer.
#define P_INT(nm, member, units, fort)                                                          \
    ParamDescriptor{nm, ParamKind::Int, 1, units, fort,                                          \
                    [](SimcParameters& p) -> int* { return &p.member; }, nullptr, nullptr, nullptr}
#define P_FLAG(nm, member, units, fort)                                                         \
    ParamDescriptor{nm, ParamKind::Flag, 1, units, fort, nullptr,                                \
                    [](SimcParameters& p) -> bool* { return &p.member; }, nullptr, nullptr}
#define P_DBL(nm, member, units, fort)                                                          \
    ParamDescriptor{nm, ParamKind::Double, 1, units, fort, nullptr, nullptr,                     \
                    [](SimcParameters& p) -> double* { return &p.member; }, nullptr}
#define P_STR(nm, member, len, units, fort)                                                     \
    ParamDescriptor{nm, ParamKind::String, len, units, fort, nullptr, nullptr, nullptr,          \
                    [](SimcParameters& p) -> std::string* { return &p.member; }}
#define P_FLAGARR(nm, member, n, units, fort)                                                   \
    ParamDescriptor{nm, ParamKind::FlagArray, n, units, fort, nullptr,                           \
                    [](SimcParameters& p) -> bool* { return p.member.data(); }, nullptr, nullptr}
#define P_DBLARR(nm, member, n, units, fort)                                                    \
    ParamDescriptor{nm, ParamKind::DoubleArray, n, units, fort, nullptr, nullptr,                \
                    [](SimcParameters& p) -> double* { return p.member.data(); }, nullptr}

std::vector<ParamDescriptor> buildTable() {
    return {
        // --- dbase.f:984-992 ---
        P_FLAG("mc_smear", mc_smear, "flag", "mc_smear_int -> mc_smear (dbase.f:984)"),
        P_INT("electron_arm", electron_arm, "code", "electron_arm (dbase.f:985)"),
        P_INT("hadron_arm", hadron_arm, "code", "hadron_arm (dbase.f:986)"),
        P_INT("use_first_cer", use_first_cer, "flag", "use_first_cer (dbase.f:987)"),
        P_DBL("transparency", transparency, "", "transparency (dbase.f:988)"),
        P_FLAG("use_benhar_sf", use_benhar_sf, "flag", "use_benhar_sf_int -> use_benhar_sf (dbase.f:989)"),
        P_STR("extra_dbase_file", extra_dbase_file, 60, "file stem", "extra_dbase_file character*60 (dbase.f:990)"),
        P_STR("random_state_file", random_state_file, 80, "path", "random_state_file character*80 (dbase.f:991)"),
        P_INT("random_seed", random_seed, "", "random_seed (dbase.f:992)"),
        // --- EXPERIMENT dbase.f:996-1010 ---
        P_DBL("Ebeam", Ebeam, "MeV", "Ebeam (dbase.f:996)"),
        P_DBL("dEbeam", dEbeam, "% (full width)", "dEbeam (dbase.f:997)"),
        P_DBL("EXPER%charge", EXPER.charge, "mC", "EXPER%charge (dbase.f:998)"),
        P_FLAG("doing_kaon", doing_kaon, "flag", "doing_kaon_int (dbase.f:999)"),
        P_INT("which_kaon", which_kaon, "code", "which_kaon (dbase.f:1000)"),
        P_FLAG("doing_pion", doing_pion, "flag", "doing_pion_int (dbase.f:1001)"),
        P_INT("which_pion", which_pion, "code", "which_pion (dbase.f:1002)"),
        P_FLAG("doing_delta", doing_delta, "flag", "doing_delta_int (dbase.f:1003)"),
        P_FLAG("doing_semi", doing_semi, "flag", "doing_semi_int (dbase.f:1004)"),
        P_FLAG("doing_hplus", doing_hplus, "flag", "doing_hplus_int (dbase.f:1005)"),
        P_FLAG("doing_pizero", doing_pizero, "flag", "doing_pizero_int (dbase.f:1006)"),
        P_INT("pizero_ngamma", pizero_ngamma, "count", "pizero_ngamma (dbase.f:1007)"),
        P_FLAG("doing_rho", doing_rho, "flag", "doing_rho_int (dbase.f:1008)"),
        P_FLAG("doing_decay", doing_decay, "flag", "doing_decay_int (dbase.f:1009)"),
        P_DBL("ctau", ctau, "cm", "ctau (dbase.f:1010)"),
        // --- DEBUG dbase.f:1014 ---
        P_FLAGARR("debug", debug, 6, "flags", "debug_int(6) -> debug(6) (dbase.f:1014)"),
        // --- TARGET dbase.f:1018-1035 ---
        P_DBL("targ%A", targ.A, "", "targ%A (dbase.f:1018)"),
        P_DBL("targ%Z", targ.Z, "", "targ%Z (dbase.f:1019)"),
        P_DBL("targ%mass_amu", targ.mass_amu, "amu", "targ%mass_amu (dbase.f:1020)"),
        P_DBL("targ%mrec_amu", targ.mrec_amu, "amu", "targ%mrec_amu (dbase.f:1021)"),
        P_DBL("targ%rho", targ.rho, "g/cm^3", "targ%rho (dbase.f:1022)"),
        P_DBL("targ%thick", targ.thick, "mg/cm^2", "targ%thick (dbase.f:1023)"),
        P_DBL("targ%xoffset", targ.xoffset, "cm", "targ%xoffset (dbase.f:1024)"),
        P_DBL("targ%yoffset", targ.yoffset, "cm", "targ%yoffset (dbase.f:1025)"),
        P_INT("targ%fr_pattern", targ.fr_pattern, "code", "targ%fr_pattern (dbase.f:1026)"),
        P_DBL("targ%fr1", targ.fr1, "cm", "targ%fr1 (dbase.f:1027)"),
        P_DBL("targ%fr2", targ.fr2, "cm", "targ%fr2 (dbase.f:1028)"),
        P_DBL("targ%zoffset", targ.zoffset, "cm", "targ%zoffset (dbase.f:1029)"),
        P_DBL("targ%angle", targ.angle, "deg", "targ%angle (dbase.f:1030)"),
        P_DBL("targ%abundancy", targ.abundancy, "%", "targ%abundancy (dbase.f:1031)"),
        P_INT("targ%can", targ.can, "code", "targ%can (dbase.f:1032)"),
        P_DBL("targ_Bangle", targ_Bangle, "deg", "targ_Bangle (dbase.f:1033)"),
        P_DBL("targ_Bphi", targ_Bphi, "deg", "targ_Bphi (dbase.f:1034)"),
        P_DBL("targ_pol", targ_pol, "", "targ_pol (dbase.f:1035)"),
        // --- E_ARM_MAIN / P_ARM_MAIN dbase.f:1037-1041 ---
        P_DBL("spec%e%P", spec.e.P, "MeV/c", "spec%e%P (dbase.f:1039)"),
        P_DBL("spec%e%theta", spec.e.theta, "deg", "spec%e%theta (dbase.f:1040)"),
        P_DBL("spec%p%P", spec.p.P, "MeV/c", "spec%p%P (dbase.f:1044)"),
        P_DBL("spec%p%theta", spec.p.theta, "deg", "spec%p%theta (dbase.f:1045)"),
        // --- E_ARM_OFFSET / P_ARM_OFFSET dbase.f:1043-1057 ---
        P_DBL("gen%ywid", gen.ywid, "cm", "gen%ywid (dbase.f:1049)"),
        P_DBL("gen%xwid", gen.xwid, "cm", "gen%xwid (dbase.f:1050)"),
        P_DBL("spec%e%offset%x", spec.e.offset.x, "cm", "spec%e%offset%x (dbase.f:1051)"),
        P_DBL("spec%e%offset%y", spec.e.offset.y, "cm", "spec%e%offset%y (dbase.f:1052)"),
        P_DBL("spec%e%offset%z", spec.e.offset.z, "cm", "spec%e%offset%z (dbase.f:1053)"),
        P_DBL("spec%e%offset%xptar", spec.e.offset.xptar, "mrad", "spec%e%offset%xptar (dbase.f:1054)"),
        P_DBL("spec%e%offset%yptar", spec.e.offset.yptar, "mrad", "spec%e%offset%yptar (dbase.f:1055)"),
        P_DBL("spec%p%offset%x", spec.p.offset.x, "cm", "spec%p%offset%x (dbase.f:1059)"),
        P_DBL("spec%p%offset%y", spec.p.offset.y, "cm", "spec%p%offset%y (dbase.f:1060)"),
        P_DBL("spec%p%offset%z", spec.p.offset.z, "cm", "spec%p%offset%z (dbase.f:1061)"),
        P_DBL("spec%p%offset%xptar", spec.p.offset.xptar, "mrad", "spec%p%offset%xptar (dbase.f:1062)"),
        P_DBL("spec%p%offset%yptar", spec.p.offset.yptar, "mrad", "spec%p%offset%yptar (dbase.f:1063)"),
        // --- MISC2INT dbase.f:1059-1061 ---
        P_INT("use_expon", use_expon, "flag", "use_expon (radc.inc, dbase.f:1067)"),
        P_INT("one_tail", one_tail, "code", "one_tail (radc.inc, dbase.f:1068)"),
        P_INT("intcor_mode", intcor_mode, "code", "intcor_mode (radc.inc, dbase.f:1069)"),
        // --- SIMULATE dbase.f:1063-1097 ---
        P_INT("ngen", ngen, "events", "ngen (dbase.f:1073)"),
        P_FLAG("hard_cuts", hard_cuts, "flag", "hard_cuts_int (dbase.f:1074)"),
        P_FLAG("using_rad", using_rad, "flag", "using_rad_int (dbase.f:1075)"),
        P_INT("spect_mode", spect_mode, "code", "spect_mode (dbase.f:1076)"),
        P_FLAG("doing_phsp", doing_phsp, "flag", "doing_phsp_int (dbase.f:1077)"),
        P_DBL("cuts%Em%min", cuts.Em.min, "MeV", "cuts%Em%min (dbase.f:1078)"),
        P_DBL("cuts%Em%max", cuts.Em.max, "MeV", "cuts%Em%max (dbase.f:1079)"),
        P_FLAG("using_Eloss", using_Eloss, "flag", "using_Eloss_int (dbase.f:1080)"),
        P_FLAG("correct_Eloss", correct_Eloss, "flag", "correct_eloss_int (dbase.f:1081)"),
        P_FLAG("correct_raster", correct_raster, "flag", "correct_raster_int (dbase.f:1082)"),
        P_FLAG("using_HMScoll", using_HMScoll, "flag", "using_HMScoll_int (dbase.f:1083)"),
        P_FLAG("using_SHMScoll", using_SHMScoll, "flag", "using_SHMScoll_int (dbase.f:1084)"),
        P_INT("deForest_flag", deForest_flag, "code", "deForest_flag (dbase.f:1085)"),
        P_INT("rad_flag", rad_flag, "code", "rad_flag (radc.inc, dbase.f:1086)"),
        P_INT("extrad_flag", extrad_flag, "code", "extrad_flag (radc.inc, dbase.f:1087)"),
        P_DBLARR("lambda", lambda, 3, "", "lambda(3) (radc.inc, dbase.f:1088)"),
        P_FLAG("using_cit_generation", using_cit_generation, "flag", "using_cit_generation_int (dbase.f:1089)"),
        P_INT("Nntu", Nntu, "flag", "Nntu (dbase.f:1090)"),
        P_FLAG("using_Coulomb", using_Coulomb, "flag", "using_Coulomb_int (dbase.f:1091)"),
        P_DBL("dE_edge_test", dE_edge_test, "MeV", "dE_edge_test (dbase.f:1092)"),
        P_FLAG("use_offshell_rad", use_offshell_rad, "flag", "use_offshell_rad_int (dbase.f:1093)"),
        P_DBL("Egamma_gen_max", Egamma_gen_max, "MeV", "Egamma_gen_max (dbase.f:1094)"),
        P_FLAG("do_fermi", do_fermi, "flag", "do_fermi_int (dbase.f:1095)"),
        P_DBL("pt_b_param", pt_b_param, "GeV^-2", "pt_b_param (dbase.f:1096)"),
        P_INT("sigc_flag", sigc_flag, "code", "sigc_flag (dbase.f:1097)"),
        P_INT("sigc_nbin", sigc_nbin, "count", "sigc_nbin (dbase.f:1098)"),
        P_DBL("sigc_kin_min", sigc_kin_min, "z or GeV^2", "sigc_kin_min (dbase.f:1099)"),
        P_DBL("sigc_kin_max", sigc_kin_max, "z or GeV^2", "sigc_kin_max (dbase.f:1100)"),
        P_DBL("sigc_kin_ind", sigc_kin_ind, "GeV^2 or z", "sigc_kin_ind (dbase.f:1101)"),
        P_FLAG("using_tgt_field", using_tgt_field, "flag", "using_tgt_field_int (dbase.f:1102)"),
        P_STR("tgt_field_file", tgt_field_file, 60, "path", "tgt_field_file character*60 (dbase.f:1103)"),
        P_DBL("drift_to_cal", drift_to_cal, "cm", "drift_to_cal (dbase.f:1104)"),
        // --- E_ARM_ACCEPT / P_ARM_ACCEPT dbase.f:1099-1113 ---
        P_DBL("SPedge%e%delta%min", SPedge.e.delta.min, "%", "SPedge%e%delta%min (dbase.f:1108)"),
        P_DBL("SPedge%e%delta%max", SPedge.e.delta.max, "%", "SPedge%e%delta%max (dbase.f:1109)"),
        P_DBL("SPedge%e%yptar%min", SPedge.e.yptar.min, "mrad", "SPedge%e%yptar%min (dbase.f:1110)"),
        P_DBL("SPedge%e%yptar%max", SPedge.e.yptar.max, "mrad", "SPedge%e%yptar%max (dbase.f:1111)"),
        P_DBL("SPedge%e%xptar%min", SPedge.e.xptar.min, "mrad", "SPedge%e%xptar%min (dbase.f:1112)"),
        P_DBL("SPedge%e%xptar%max", SPedge.e.xptar.max, "mrad", "SPedge%e%xptar%max (dbase.f:1113)"),
        P_DBL("SPedge%p%delta%min", SPedge.p.delta.min, "%", "SPedge%p%delta%min (dbase.f:1117)"),
        P_DBL("SPedge%p%delta%max", SPedge.p.delta.max, "%", "SPedge%p%delta%max (dbase.f:1118)"),
        P_DBL("SPedge%p%yptar%min", SPedge.p.yptar.min, "mrad", "SPedge%p%yptar%min (dbase.f:1119)"),
        P_DBL("SPedge%p%yptar%max", SPedge.p.yptar.max, "mrad", "SPedge%p%yptar%max (dbase.f:1120)"),
        P_DBL("SPedge%p%xptar%min", SPedge.p.xptar.min, "mrad", "SPedge%p%xptar%min (dbase.f:1121)"),
        P_DBL("SPedge%p%xptar%max", SPedge.p.xptar.max, "mrad", "SPedge%p%xptar%max (dbase.f:1122)"),
    };
}

std::string fmtDouble(double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.10g", v);
    return buf;
}

}  // namespace

std::vector<std::string> ParamDescriptor::jsonPath() const {
    std::vector<std::string> out;
    std::string cur;
    for (const char* c = name; *c; ++c) {
        if (*c == '%') {
            out.push_back(cur);
            cur.clear();
        } else {
            cur += *c;
        }
    }
    out.push_back(cur);
    return out;
}

const std::vector<ParamDescriptor>& ParameterTable() {
    static const std::vector<ParamDescriptor> table = buildTable();
    return table;
}

const ParamDescriptor* FindParameter(const std::string& fortran_name) {
    const std::string key = lower(fortran_name);
    for (const auto& d : ParameterTable())
        if (lower(d.name) == key) return &d;
    return nullptr;
}

std::string FormatValue(const ParamDescriptor& d, const SimcParameters& cp) {
    auto& p = const_cast<SimcParameters&>(cp);
    std::ostringstream os;
    switch (d.kind) {
        case ParamKind::Int: os << *d.asInt(p); break;
        case ParamKind::Flag: os << (*d.asBool(p) ? "T" : "F"); break;
        case ParamKind::Double: os << fmtDouble(*d.asDouble(p)); break;
        case ParamKind::String: os << "'" << *d.asString(p) << "'"; break;
        case ParamKind::FlagArray:
            for (int i = 0; i < d.size; ++i) os << (i ? "," : "") << (d.asBool(p)[i] ? "T" : "F");
            break;
        case ParamKind::DoubleArray:
            for (int i = 0; i < d.size; ++i) os << (i ? "," : "") << fmtDouble(d.asDouble(p)[i]);
            break;
    }
    return os.str();
}

Json ToJson(const SimcParameters& cp) {
    auto& p = const_cast<SimcParameters&>(cp);
    Json j = Json::object();
    for (const auto& d : ParameterTable()) {
        Json* node = &j;
        auto path = d.jsonPath();
        for (size_t i = 0; i + 1 < path.size(); ++i) node = &(*node)[path[i]];
        Json& leaf = (*node)[path.back()];
        switch (d.kind) {
            case ParamKind::Int: leaf = *d.asInt(p); break;
            case ParamKind::Flag: leaf = *d.asBool(p); break;
            case ParamKind::Double: leaf = *d.asDouble(p); break;
            case ParamKind::String: leaf = *d.asString(p); break;
            case ParamKind::FlagArray:
                leaf = Json::array();
                for (int i = 0; i < d.size; ++i) leaf.push_back(d.asBool(p)[i]);
                break;
            case ParamKind::DoubleArray:
                leaf = Json::array();
                for (int i = 0; i < d.size; ++i) leaf.push_back(d.asDouble(p)[i]);
                break;
        }
    }
    return j;
}

namespace {

void warn(std::vector<std::string>* w, const std::string& msg) {
    if (w) w->push_back(msg);
}

// CTP rounding of a real value stored into an integer variable (CTP/daVar.h:131 floatToLong).
int ctpRound(double x) { return static_cast<int>(x > 0.0 ? x + 0.5 : x - 0.5); }

bool toNumber(const Json& v, double& out) {
    if (v.is_boolean()) { out = v.get<bool>() ? 1.0 : 0.0; return true; }
    if (v.is_number()) { out = v.get<double>(); return true; }
    return false;
}

void assignLeaf(const ParamDescriptor& d, const Json& v, SimcParameters& p, const std::string& key,
                std::vector<std::string>* w) {
    double x = 0.0;
    switch (d.kind) {
        case ParamKind::Int:
            if (!toNumber(v, x)) { warn(w, "parameter '" + key + "': expected a number, ignored"); return; }
            if (x != std::floor(x)) warn(w, "parameter '" + key + "': non-integer value rounded (CTP floatToLong)");
            *d.asInt(p) = ctpRound(x);
            return;
        case ParamKind::Flag:
            if (!toNumber(v, x)) { warn(w, "parameter '" + key + "': expected bool or number, ignored"); return; }
            // dbase.f:convert_to_logical — integer > 0 means .true.
            *d.asBool(p) = (v.is_boolean() ? v.get<bool>() : ctpRound(x) > 0);
            return;
        case ParamKind::Double:
            if (!toNumber(v, x)) { warn(w, "parameter '" + key + "': expected a number, ignored"); return; }
            *d.asDouble(p) = x;
            return;
        case ParamKind::String:
            if (!v.is_string()) { warn(w, "parameter '" + key + "': expected a string, ignored"); return; }
            *d.asString(p) = v.get<std::string>();
            if (static_cast<int>(d.asString(p)->size()) > d.size) {
                warn(w, "parameter '" + key + "': string longer than Fortran character*" +
                            std::to_string(d.size) + ", truncated");
                d.asString(p)->resize(d.size);
            }
            return;
        case ParamKind::FlagArray:
        case ParamKind::DoubleArray: {
            if (!v.is_array()) { warn(w, "parameter '" + key + "': expected an array, ignored"); return; }
            if (static_cast<int>(v.size()) > d.size)
                warn(w, "parameter '" + key + "': more than " + std::to_string(d.size) +
                            " elements, extra ignored");
            for (int i = 0; i < d.size && i < static_cast<int>(v.size()); ++i) {
                if (!toNumber(v[i], x)) { warn(w, "parameter '" + key + "': non-numeric element ignored"); continue; }
                if (d.kind == ParamKind::FlagArray)
                    d.asBool(p)[i] = v[i].is_boolean() ? v[i].get<bool>() : ctpRound(x) > 0;
                else
                    d.asDouble(p)[i] = x;
            }
            return;
        }
    }
}

void walk(const Json& node, const std::string& prefix, SimcParameters& p, std::vector<std::string>* w) {
    for (auto it = node.begin(); it != node.end(); ++it) {
        const std::string key = prefix.empty() ? it.key() : prefix + "%" + it.key();
        if (!it.key().empty() && it.key()[0] == '_') continue;  // "_comment" style keys are ignored
        if (const ParamDescriptor* d = FindParameter(key)) {
            assignLeaf(*d, it.value(), p, key, w);
        } else if (it.value().is_object()) {
            walk(it.value(), key, p, w);
        } else {
            warn(w, "unknown parameter '" + key + "' (not registered in dbase.f:regallvars), ignored");
        }
    }
}

}  // namespace

void FromJson(const Json& j, SimcParameters& p, std::vector<std::string>* warnings) {
    if (!j.is_object()) {
        warn(warnings, "parameters: expected a JSON object");
        return;
    }
    walk(j, "", p, warnings);
}

}  // namespace simc
